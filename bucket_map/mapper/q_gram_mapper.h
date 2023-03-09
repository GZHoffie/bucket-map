#ifndef BUCKET_MAP_Q_GRAM_MAP_H
#define BUCKET_MAP_Q_GRAM_MAP_H


#include "../indexer/bucket_fm_indexer.h"
#include "../utils.h"
#include "./mapper.h"
#include "./quality_filter.h"

#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <chrono>
#include <fstream>
#include <ranges>
#include <iostream>
#include <algorithm>

#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

using seqan3::operator""_shape;


template<unsigned int NUM_BUCKETS>
class fault_tolerate_filter {
    /**
     * @brief Novel data structure for fast bit-parallel pop-count of multiple BitSet objects.
     *        Eliminating those BitSets with more than `num_fault_tolerance` zeros.
     *        This can help filter out the correct buckets fast.
     */
private:
    unsigned int num_fault_tolerance;
    unsigned int allowed_max_candidate_buckets;
    std::vector<std::bitset<NUM_BUCKETS>> filters;

    std::vector<int> _set_bits(int index) {
        /**
         * @brief Find all set bits in filters[index]
         * @param index indicates which filter we want to check.
         */
        std::vector<int> res;
        if (index >= num_fault_tolerance || index < 0) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The input index "
                                 << index << " exceeds the fault tolerance level." << '\n';
            return res;
        } else if (filters[index].none()) {
            return res;
        }
        for (int i = filters[index]._Find_first(); i < NUM_BUCKETS; i = filters[index]._Find_next(i)) {
            res.push_back(i);
        }
        return res;
    }

public:
    fault_tolerate_filter(unsigned int fault, unsigned int num_buckets) {
        num_fault_tolerance = fault;
        for (int i = 0; i < num_fault_tolerance; i++) {
            std::bitset<NUM_BUCKETS> filter;
            filter.set();
            filters.push_back(filter);
        }
        allowed_max_candidate_buckets = num_buckets;
    }

    void reset() {
        for (int i = 0; i < num_fault_tolerance; i++) {
            filters[i].set();
        }
    }

    void read(std::bitset<NUM_BUCKETS> input) {
        /**
         * @brief Read a set of bits and filter out the indices at which they are 0.
         * @param input a bitset with size NUM_BUCKETS indicating whether a certain k-mer appear
         *              in each bucket. (0 at position i indicates that the k-mer doesn't exist 
         *              in bucket i)
         */
        // Record how many time we see a 0 for each bucket
        for(unsigned int i = 0; i < num_fault_tolerance - 1; i++) {
            filters[i] &= (filters[i+1] | input);
        }
        // Modify filters[size-1]
        filters[num_fault_tolerance - 1] &= input;
    }

    std::vector<int> best_results() {
        /**
         * @brief Return the buckets that contains the most number of k-mers.
         */
        std::vector<int> res;
        for (int i = num_fault_tolerance-1; i >= 0; i--) {
            res = _set_bits(i);
            if (!res.empty()) {
                return res;
            }
        }
        return res;
    }

    std::vector<int> ok_results() {
        /**
         * @brief Return the best `allowed_max_candidate_buckets` buckets.
         */
        std::vector<int> res;
        for (int i = num_fault_tolerance-1; i >= 0; i--) {
            if (i == 0 || (filters[i-1].count() > allowed_max_candidate_buckets && filters[i].count() != 0)) {
                res = _set_bits(i);
                return res;
            }
        }
        return res;
    }

    std::vector<int> all_results() {
        /**
         * @brief Return all the acceptable buckets.
         */
        return _set_bits(0);
    }

    int _check_bucket(unsigned int bucket_id) {
        /**
         * @brief Check how many errors actually present in the true bucket.
         */
        for (int i = num_fault_tolerance-1; i >= 0; i--) {
            if (filters[i][bucket_id]) {
                return i;
            }
        }
        return 0;
    }
};

template<unsigned int NUM_BUCKETS>
class distinguishability_filter {
    /**
     * @brief A filter that filter out Q-grams that appear in most of the buckets.
     */
private:
    unsigned int threshold;

public:
    std::vector<unsigned int> zeros;

    distinguishability_filter(float distinguishability) {
        /**
         * @brief Initializer of the filter.
         * @param distinguishability the percentage of zeros in the bitset for each Q-gram.
         */
        threshold = (unsigned int) (distinguishability * NUM_BUCKETS);
    }

    void read(const std::vector<std::bitset<NUM_BUCKETS>>& q_grams_index) {
        int valid_q_grams = 0;
        for (auto &i: q_grams_index) {
            unsigned int ones = i.count();
            if (NUM_BUCKETS - ones > threshold) {
                valid_q_grams++;
            }
            if (ones == 0) {
                // the q-gram doesn't appear at all
                zeros.push_back(NUM_BUCKETS);
            } else {
                zeros.push_back(NUM_BUCKETS - ones);
            }
        }
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of Q-grams with distinguishability >= " << ((float) threshold) / NUM_BUCKETS << ": " 
                             << valid_q_grams << " (" << ((float) valid_q_grams) / q_grams_index.size() * 100 << "%).\n";
    }

    bool is_highly_distinguishable(unsigned int kmer_hash) {
        return zeros[kmer_hash] >= threshold;
    }
};






template<unsigned int NUM_BUCKETS>
class q_gram_mapper : public mapper {
private:
    // q-gram index related information
    seqan3::shape q_gram_shape;
    std::vector<std::bitset<NUM_BUCKETS>> q_grams_index;
    unsigned int q;
    unsigned int size;

    // bucket index related information
    unsigned int bucket_length;
    unsigned int read_length;

    // mapper related information
    unsigned int num_samples;
    unsigned int num_fault_tolerance;
    unsigned int allowed_max_candidate_buckets;

    // filter that filter out the most possible bucket
    fault_tolerate_filter<NUM_BUCKETS>* filter;

    // Q-gram filters for map efficiency
    distinguishability_filter<NUM_BUCKETS>* dist_filter;
    // quality filter for avoiding sequencing errors
    unsigned int average_quality;


    // sampling k-mers
    Sampler* sampler;

    std::bitset<NUM_BUCKETS> _bitset_from_bytes(const std::vector<char>& buf) {
        /**
         * @brief Convert 8-byte chars to bitset.
         * Adopted from https://stackoverflow.com/a/7463972.
         */
        assert(buf.size() == ((NUM_BUCKETS + 7) >> 3));
        std::bitset<NUM_BUCKETS> result;
        for (unsigned int j = 0; j < NUM_BUCKETS; j++)
            result[j] = ((buf.at(j>>3) >> (j & 7)) & 1);
        return result;
    }

    struct _dna4_traits : seqan3::sequence_file_input_default_traits_dna {
        /**
         * @brief Syntax for reading the query file.
         * 
         */
        using sequence_alphabet = seqan3::dna4; // instead of dna5
 
        template <typename alph>
        using sequence_container = std::vector<alph>; // must be defined as a template!
    };


public:
    q_gram_mapper(unsigned int bucket_len, unsigned int read_len, seqan3::shape shape, 
                  unsigned int samples, unsigned int fault, float distinguishability,
                  unsigned int quality_threshold = 35, 
                  unsigned int num_candidate_buckets = 30) : mapper() {
        // initialize private variables
        bucket_length = bucket_len;
        read_length = read_len;
        
        q_gram_shape = shape;
        size = std::ranges::size(shape);
        q = shape.count();
        
        num_samples = samples;
        num_fault_tolerance = fault;
        allowed_max_candidate_buckets = num_candidate_buckets;

        // initialize filter
        filter = new fault_tolerate_filter<NUM_BUCKETS>(num_fault_tolerance, num_candidate_buckets);
        dist_filter = new distinguishability_filter<NUM_BUCKETS>(distinguishability);
        average_quality = quality_threshold * q;

        // Initialize sampler
        sampler = new Sampler(num_samples);
    }

    ~q_gram_mapper() {
        delete filter;
        delete dist_filter;
        delete sampler;
    }


    void load(std::filesystem::path const & index_directory, const std::string& indicator) {
        /**
         * @brief Look for index and pattern file inside the index_directory,
         *        read the files and store the values in class attribute.
         * @param index_directory the directory that stores the q_gram_count_file.
         */
        // q_grams_index should be empty
        if (!q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is not empty. Terminating load.\n";
            return;
        }
        // initialize q_gram index
        int total_q_grams = (int) pow(4, q);
        int num_chars_per_q_gram = (NUM_BUCKETS + 7) >> 3;

        // Read the index file
        std::ifstream file(index_directory / (indicator + ".qgram"));
        if (file) {
            Timer clock;
            clock.tick();

            // read several bytes from file at a time
            std::vector<char> buffer(num_chars_per_q_gram);
            for (unsigned int i = 0; i < total_q_grams; i++) {
                file.read(&buffer[0], sizeof(unsigned char) * num_chars_per_q_gram);
                q_grams_index.push_back(_bitset_from_bytes(buffer));
            }

            dist_filter->read(q_grams_index);

            // Complete the read
            clock.tock();
            seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for loading index files: " 
                                 << clock.elapsed_seconds() << " s." << '\n';
            seqan3::debug_stream << "[INFO]\t\t" << "Successfully loaded " 
                                 << index_directory / (indicator + ".qgram") << "." << '\n';
        }
    }


    std::vector<int> query(const std::vector<int>& q_gram_hash) {
        /**
         * @brief From `q_grams_index`, determine where the sequence may be coming from.
         * @param q_gram_hash the vector containing all hash values of q-grams in the
         *                    query sequence.
         * @returns a vector of integers indicating the possible regions that the sequence
         *          may belong to.
         */
        // q_grams_index should not be empty
        if (q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is empty. Cannot accept query.\n";
            std::vector<int> res;
            return res;
        }
        // Reset the filter
        filter->reset();

        // insert those samples into the filter
        for (int h : q_gram_hash) {
            filter->read(q_grams_index[h]);
        }
        return filter->best_results();
        //return filter->all_results();
        //return filter->ok_results();
    }

    std::vector<int> query_sequence(const std::vector<seqan3::dna4>& sequence,
                                    const std::vector<seqan3::phred42>& quality) {
        /**
         * @brief From `q_grams_index`, determine where the sequence may be coming from.
         * @param sequence A dna4 vector containing the query sequence.
         * @param quality the read quality of the sequence.
         * @returns a vector of integers indicating the possible regions that the sequence
         *          may belong to.
         */
        
        // get satisfactory k-mers that passes through quality filter and distinguishability filter
        auto kmers = sequence | seqan3::views::kmer_hash(q_gram_shape);
        auto kmer_qualities = quality | seqan3::views::kmer_quality(q_gram_shape);
        int num_kmers = kmers.size();
        //seqan3::debug_stream << std::views::iota(0, num_kmers) << "\n";

        auto good_kmers = std::ranges::iota_view{0, num_kmers} | std::views::filter([&](unsigned int i) {
                              return dist_filter->is_highly_distinguishable(kmers[i]) && kmer_qualities[i] >= average_quality;
                          }) | std::views::transform([&](unsigned int i) {
                              return kmers[i];
                          });
        std::vector<unsigned int> hash_values(good_kmers.begin(), good_kmers.end());

        // if not enough q-grams ramained to determine the exact location, simply ignore this query sequence.
        if (hash_values.size() < 0.2 * num_samples){
            std::vector<int> res;
            return res;
        }

        std::vector<int> samples;
        /*
        // Randomly sample `num_samples` q-grams for query.
        std::sample(hash_values.begin(), hash_values.end(), 
                    std::back_inserter(samples), num_samples,
                    std::mt19937{std::random_device{}()});
        */
        // Deterministically sample from the hash values.
        sampler->sample_deterministically(hash_values.size()-1);
        for (auto sample : sampler->samples) {
            samples.push_back(hash_values[sample]);
        }

        auto res = query(samples);
        if (res.size() > allowed_max_candidate_buckets) {
            res.clear();
        }
        
        return res;
        
    }



    std::pair<std::vector<std::vector<unsigned int>>, std::vector<std::vector<unsigned int>>>
    map(std::filesystem::path const & sequence_file) {
        /**
         * @brief Read a query fastq file and output the ids of the sequence that are mapped 
         *        to each file.
         */

        std::vector<std::vector<unsigned int>> res_orig;
        std::vector<std::vector<unsigned int>> res_rev_comp;
        seqan3::sequence_file_input<_dna4_traits> fin{sequence_file};
        // initialize returning result
        for (int i = 0; i < NUM_BUCKETS; i++) {
            std::vector<unsigned int> sequence_ids_orig;
            res_orig.push_back(sequence_ids_orig);
            std::vector<unsigned int> sequence_ids_rev_comp;
            res_rev_comp.push_back(sequence_ids_rev_comp);
        }

        Timer clock;
        clock.tick();
        for (auto & rec : fin) {
            // find if read is present in the buckets
            for (auto & bucket : query_sequence(rec.sequence(), rec.base_qualities())) {
                res_orig[bucket].push_back(num_records);
            }
            // find the reverse complement of the read in the buckets
            auto rec_rev_comp = rec.sequence() | std::views::reverse | seqan3::views::complement;
            auto quality_rev = rec.base_qualities() | std::views::reverse;
            std::vector<seqan3::phred42> qual_rev(quality_rev.begin(), quality_rev.end());
            seqan3::dna4_vector sequence_rev_comp(rec_rev_comp.begin(), rec_rev_comp.end());
            for (auto & bucket : query_sequence(sequence_rev_comp, qual_rev)) {
                res_rev_comp[bucket].push_back(num_records);
            }
            ++num_records;
        }
        clock.tock();
        float time = clock.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for bucket mapping: " 
                             << time << " s (" << time * 1000 * 1000 / num_records << " μs/seq).\n";
        return std::make_pair(res_orig, res_rev_comp);
    }

    std::vector<std::vector<int>> _query_file(std::filesystem::path sequence_file) {
        /**
         * * This function is just for benchmarking.
         * @brief Read a query fastq file and output the bucket ids each query belongs to.
         * TODO: include the quality information for fastq.
         */
        std::vector<std::vector<int>> res;
        seqan3::sequence_file_input<_dna4_traits> fin{sequence_file};
        Timer clock;
        clock.tick();
 
        for (auto & rec : fin) {
            std::vector<int> buckets = query_sequence(rec.sequence());
            res.push_back(buckets);
        }
        clock.tock();
        float time = clock.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for bucket query: " 
                             << time << " s (" << time * 1000 * 1000 / res.size() << " μs/seq).\n";
        return res;
    }

    void _check_ground_truth(const std::vector<std::vector<int>>& query_results, std::filesystem::path ground_truth_file) {
        /**
         * * This function is just for benchmarking.
         * @brief Check the performance of query results against the ground truth.
         *        Output necessary information.
         * @note The ground truth file must be the one generated together with sequence
         *       fastq file.
         * TODO: also check the exact location.
         */
        std::ifstream is(ground_truth_file);
        int bucket, exact_location;
        std::string cigar;

        // Test statistics
        int correct_map = 0;
        int total_bucket_numbers = 0;
        std::map<int, int> bucket_number_map;

        for (int i = 0; i < query_results.size(); i++) {
            is >> bucket >> exact_location >> cigar;
            std::vector<int> buckets = query_results[i];

            if (std::find(buckets.begin(), buckets.end(), bucket) != buckets.end()) {
                correct_map++;
            }
            total_bucket_numbers += buckets.size();
            ++bucket_number_map[buckets.size()];
        }

        // output information
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of sequences: " 
                             << query_results.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Correct bucket predictions: " 
                             << correct_map << " (" << ((float) correct_map) / query_results.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Average number of buckets returned: " 
                             << ((float) total_bucket_numbers) / query_results.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of sequences that have no candidate bucket: " 
                             << bucket_number_map[0] << " (" << ((float) bucket_number_map[0]) / query_results.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of uniquely mapped sequences: " 
                             << bucket_number_map[1] << " (" << ((float) bucket_number_map[1]) / query_results.size() * 100 << "%).\n";
        int small_bucket_numbers = 0;
        for (int i = 1; i <= 5; i++) {
            small_bucket_numbers += bucket_number_map[i];
        }
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of sequences mapped to <= 5 buckets: " 
                             << small_bucket_numbers << " (" << ((float) small_bucket_numbers) / query_results.size() * 100 << "%).\n";
        for (int i = 6; i <= 10; i++) {
            small_bucket_numbers += bucket_number_map[i];
        }
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of sequences mapped to <= 10 buckets: " 
                             << small_bucket_numbers << " (" << ((float) small_bucket_numbers) / query_results.size() * 100 << "%).\n";
    }

    void reset() {
        /**
         * @brief Release the memory, mainly used by the q_grams_index.
         * TODO: release other variables.
         */
        q_grams_index.clear();
        q_grams_index.shrink_to_fit();
    }
};


#endif