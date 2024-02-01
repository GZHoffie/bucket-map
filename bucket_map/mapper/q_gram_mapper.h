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
#include <seqan3/search/views/minimiser.hpp>

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

    std::vector<unsigned int> _set_bits(int index) {
        /**
         * @brief Find all set bits in filters[index]
         * @param index indicates which filter we want to check.
         */
        std::vector<unsigned int> res;
        if (index >= num_fault_tolerance || index < 0) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The input index "
                                 << index << " exceeds the fault tolerance level." << '\n';
            return res;
        } else if (filters[index].none()) {
            return res;
        }
        for (unsigned int i = filters[index]._Find_first(); i < NUM_BUCKETS; i = filters[index]._Find_next(i)) {
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
        for (unsigned int i = 0; i < num_fault_tolerance; i++) {
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

    std::vector<unsigned int> best_results() {
        /**
         * @brief Return the buckets that contains the most number of k-mers.
         */
        std::vector<unsigned int> res;
        for (int i = num_fault_tolerance-1; i >= 0; i--) {
            res = _set_bits(i);
            if (!res.empty()) {
                return res;
            }
        }
        return res;
    }

    std::vector<unsigned int> ok_results() {
        /**
         * @brief Return the best `allowed_max_candidate_buckets` buckets.
         */
        std::vector<unsigned int> res;
        for (unsigned int i = num_fault_tolerance-1; i >= 0; i--) {
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
    unsigned int k, q;
    unsigned int Q_BITS;
    std::vector<int>* kmer_to_index;

    int index_of_kmer(unsigned int k_mer_hash) {
        if (k_mer_hash < kmer_to_index->size()) return kmer_to_index->at(k_mer_hash);
        return -1;
    }

public:
    std::vector<unsigned int> zeros;

    distinguishability_filter(float distinguishability, unsigned int index_seed, unsigned int query_seed,
                              std::vector<int>* kmer_index) {
        /**
         * @brief Initializer of the filter.
         * @param distinguishability the percentage of zeros in the bitset for each Q-gram.
         */
        threshold = (unsigned int) (distinguishability * NUM_BUCKETS);
        q = index_seed;
        k = query_seed;
        Q_BITS = pow(4, q) - 1;

        kmer_to_index = kmer_index;
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
        for (unsigned int i = 0; i <= k - q; i++) {
            unsigned int h = ((kmer_hash >> (2 * i)) & Q_BITS);
                int index = index_of_kmer(h);
                if (index >= 0 && zeros[index] >= threshold) return true;
            }
        return false;
    }
};






template<unsigned int NUM_BUCKETS>
class q_gram_mapper : public mapper {
private:
    // q-gram index related information
    std::vector<std::bitset<NUM_BUCKETS>> q_grams_index;
    uint8_t q; // index seed length
    uint8_t k; // query seed length
    unsigned int Q_BITS; // bits representing the length of q-gram, used to find q-grams contained in a k-mer
    unsigned int size;

    // bucket index related information
    unsigned int bucket_length;
    unsigned int read_length;
    std::vector<int> kmer_to_index;

    // mapper related information
    unsigned int num_samples;
    unsigned int num_fault_tolerance;
    unsigned int allowed_max_candidate_buckets;

    // filter that filter out the most possible bucket
    fault_tolerate_filter<NUM_BUCKETS>* filter;

    // Q-gram filters for map efficiency
    distinguishability_filter<NUM_BUCKETS>* dist_filter;
    // quality filter for avoiding sequencing errors
    unsigned int min_base_quality;
    std::vector<bool> high_quality_kmers;


    // sampling k-mers
    Sampler* sampler;
    Sampler* segment_sampler;

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

    /**
     * @brief Find all high quality k-mers.
     * 
     * @param quality the vector of quality alphabets.
     * fills in `high_quality_kmers`: A vector with the i-th element representing whether
     *                                the i-th k-mer is high-quality.
     * TODO: support spaced k-mer.
     */
    void _high_quality_kmers(const std::vector<seqan3::phred94>& quality) {
        unsigned int consecutive_high_quality_base = 0, index = 0;
        high_quality_kmers.clear();
        for (auto & qual: quality) {
            if (qual.to_rank() >= min_base_quality) {
                consecutive_high_quality_base++;
                if (consecutive_high_quality_base >= k) {
                    high_quality_kmers.push_back(true);
                    continue;
                }
            } else {
                consecutive_high_quality_base = 0;
            }
            index++;
            if (index >= k) {
                high_quality_kmers.push_back(false);
            }
        }
        assert(high_quality_kmers.size() == quality.size() - k + 1);
    }


public:
    q_gram_mapper(unsigned int bucket_len, unsigned int read_len, 
                  uint8_t query_seed_length, uint8_t index_seed_length,
                  unsigned int samples, unsigned int fault, float distinguishability,
                  unsigned int quality_threshold = 35, 
                  unsigned int num_candidate_buckets = 30,
                  unsigned int num_segment_samples = 5) : mapper() {
        // initialize private variables
        bucket_length = bucket_len;
        read_length = read_len;
        
        k = query_seed_length;
        q = index_seed_length;
        Q_BITS = pow(4, q) - 1;
        seqan3::debug_stream << "[INFO]\t\tSet query seed length to be " << k << ", and index seed length " << q << ".\n";
        
        num_samples = samples;
        num_fault_tolerance = fault;
        allowed_max_candidate_buckets = num_candidate_buckets;

        // initialize filter
        filter = new fault_tolerate_filter<NUM_BUCKETS>(num_fault_tolerance, num_candidate_buckets);
        dist_filter = new distinguishability_filter<NUM_BUCKETS>(distinguishability, q, k, &kmer_to_index);
        min_base_quality = quality_threshold * k;

        // Initialize sampler
        sampler = new Sampler(num_samples);
        segment_sampler = new Sampler(num_segment_samples);
    }

    ~q_gram_mapper() {
        delete filter;
        delete dist_filter;
        delete sampler;
        delete segment_sampler;
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

        // Read the kmer index with FracMinHash file
        int sampled_q_grams = 0;
        std::ifstream frac_min_hash_file(index_directory / (indicator + ".kmers_index"));
        if (frac_min_hash_file) {
            int index;
            for (unsigned int i = 0; i < pow(4, q); i++) {
                frac_min_hash_file >> index;
                if (index >= 0) sampled_q_grams += 1;
                kmer_to_index.push_back(index);
            }
            seqan3::debug_stream << "[INFO]\t\t" << "Successfully loaded " 
                                 << index_directory / (indicator + ".kmers_index") << "." << '\n';
        }

        // initialize q_gram index
        int num_chars_per_q_gram = (NUM_BUCKETS + 7) >> 3;

        // Read the index file
        std::ifstream file(index_directory / (indicator + ".qgram"));
        if (file) {
            Timer clock;
            clock.tick();

            // read several bytes from file at a time
            std::vector<char> buffer(num_chars_per_q_gram);
            for (unsigned int i = 0; i < sampled_q_grams; i++) {
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

    int index_of_kmer(unsigned int k_mer_hash) {
        if (k_mer_hash < kmer_to_index.size()) return kmer_to_index[k_mer_hash];
        return -1;
    }


    std::vector<unsigned int> query(const std::vector<unsigned int>& kmer_hash) {
        /**
         * @brief From `q_grams_index`, determine where the sequence may be coming from.
         * @param kmer_hash the vector containing all hash values of k-mers in the
         *                  query sequence.
         * @returns a vector of integers indicating the possible regions that the sequence
         *          may belong to.
         */
        // q_grams_index should not be empty
        if (q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is empty. Cannot accept query.\n";
            std::vector<unsigned int> res;
            return res;
        }
        // Reset the filter
        filter->reset();

        // insert sampled k-mers into the filter
        for (unsigned int h : kmer_hash) {
            // find the buckets that contains all the smaller q-grams in the k-mer
            std::bitset<NUM_BUCKETS> bf_res;
            bf_res.set();
            for (unsigned int i = 0; i <= k - q; i++) {
                unsigned int q_gram_hash = ((h >> (2 * i)) & Q_BITS);
                int index = index_of_kmer(q_gram_hash);
                if (index >= 0) bf_res &= q_grams_index[index];
            }
            filter->read(bf_res);
        }
        return filter->best_results();
        //return filter->all_results();
        //return filter->ok_results();
    }

    std::pair<std::vector<unsigned int>, std::vector<unsigned int>> 
    query_sequence(const std::vector<seqan3::dna4>& sequence,
                   const std::vector<seqan3::phred94>& quality) {
        /**
         * @brief From `q_grams_index`, determine where the sequence may be coming from.
         * @param sequence A dna4 vector containing the query sequence.
         * @param quality the read quality of the sequence.
         * @returns two vectors of integers indicating the possible regions that the sequence
         *          may belong to. The first represent the buckets that the read can be mapped 
         *          directly, and the second is the ones where reads can be mapped after
         *          reverse complementing.
         */
        // initialize results
        std::vector<unsigned int> candidates_orig;
        std::vector<unsigned int> candidates_rev_comp;
        
        // get satisfactory k-mers that passes through quality filter and distinguishability filter
        auto kmers = sequence | seqan3::views::kmer_hash(seqan3::ungapped{k});
        auto kmer_qualities = quality | seqan3::views::kmer_quality(seqan3::ungapped{k});
        //_high_quality_kmers(quality);
        std::vector<unsigned int> qualities(kmer_qualities.begin(), kmer_qualities.end());
        int num_kmers = kmers.size();

        auto good_kmers = std::ranges::iota_view{0, num_kmers} | std::views::filter([&](int i) {
                              return dist_filter->is_highly_distinguishable(kmers[i]) && qualities[i] >= min_base_quality;
                          }) | std::views::transform([&](int i) {
                              return (unsigned int)kmers[i];
                          });
        std::vector<unsigned int> hash_values(good_kmers.begin(), good_kmers.end());

        // if not enough q-grams ramained to determine the exact location, simply ignore this query sequence.
        if (hash_values.size() < 0.2 * num_samples){
            return std::make_pair(candidates_orig, candidates_rev_comp);
        }

        std::vector<unsigned int> samples_orig;
        /*
        // Randomly sample `num_samples` q-grams for query.
        std::sample(hash_values.begin(), hash_values.end(), 
                    std::back_inserter(samples_orig), num_samples,
                    std::mt19937{std::random_device{}()});
        */
        // Deterministically sample from the hash values.
        sampler->sample_deterministically(hash_values.size()-1);
        for (auto sample : sampler->samples) {
            samples_orig.push_back(hash_values[sample]);
        }
        
        candidates_orig = query(samples_orig);

        // query the reverse complements of the sampled k-mers
        auto samples_rev_comp = samples_orig | std::views::transform([&](unsigned int hash) {
            return hash_reverse_complement(hash, k);
        });
        std::vector<unsigned int> samples_rev_comp_vec(samples_rev_comp.begin(), samples_rev_comp.end());
        candidates_rev_comp = query(samples_rev_comp_vec);
        if (candidates_orig.size() > allowed_max_candidate_buckets) {
            candidates_orig.clear();
        }
        if (candidates_rev_comp.size() > allowed_max_candidate_buckets) {
            candidates_rev_comp.clear();
        }
        
        return std::make_pair(candidates_orig, candidates_rev_comp);
        
    }


    std::pair<segments_t, segments_t>
    map(std::filesystem::path const & sequence_file) {
        /**
         * @brief Read a query fastq file and output the ids of the sequence that are mapped 
         *        to each file.
         */
        // benchmarking percentage of mapped reads.
        unsigned int mapped_reads = 0;
        unsigned int num_buckets_orig = 0, num_buckets_rev_comp = 0;
        
        seqan3::sequence_file_input<_phred94_traits> fin{sequence_file};
        // initialize returning result
        segments_t res_orig;
        segments_t res_rev_comp;
        for (int i = 0; i < NUM_BUCKETS; i++) {
            std::vector<segment_info_t> sequence_ids_orig;
            res_orig.push_back(sequence_ids_orig);
            std::vector<segment_info_t> sequence_ids_rev_comp;
            res_rev_comp.push_back(sequence_ids_rev_comp);
        }

        Timer clock;
        clock.tick();
        for (auto & rec : fin) {
            bool mapped = false;

            // the positions of the read where we want to query.
            std::vector<unsigned int> starting_positions{0};

            // if long read, break it down into several shorter reads.
            if (rec.sequence().size() > 2 * read_length) {
                segment_sampler->sample_deterministically(rec.sequence().size() - read_length - 1);
                starting_positions = segment_sampler->samples;
            }
            for (auto i : starting_positions) {
                // copy the segment of length read_length
                //TODO: optimize this by not copying the segment, but simply indicate the range.
                auto begin = i;
                auto end = std::min(i + read_length, (unsigned int)rec.sequence().size());
                seqan3::dna4_vector segment_sequence(rec.sequence().begin() + begin, rec.sequence().begin() + end);
                std::vector<seqan3::phred94> segment_quality(rec.base_qualities().begin() + begin, rec.base_qualities().begin() + end);

                // find if the segment is present in the buckets
                auto [buckets_orig, buckets_rev_comp] = query_sequence(segment_sequence, segment_quality);
                segment_info_t seg{num_records, (int)i};
                for (auto & b : buckets_orig) {
                    res_orig[b].push_back(seg);
                }
                for (auto & b_ : buckets_rev_comp) {
                    res_rev_comp[b_].push_back(seg);
                }
                if (!buckets_orig.empty() || !buckets_rev_comp.empty()) {
                    mapped = true;
                    num_buckets_orig += buckets_orig.size();
                    num_buckets_rev_comp += buckets_rev_comp.size();
                }
            }

            if (mapped) {
                ++mapped_reads;
            }
            ++num_records;
        }
        clock.tock();
        float time = clock.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for bucket mapping: " 
                             << time << " s (" << time * 1000 * 1000 / num_records << " μs/seq).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of reads that have at least one candidate bucket: " 
                             << mapped_reads << "  (" << ((float)mapped_reads) / num_records * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Average number of buckets an original read is mapped to: " 
                             << ((float) num_buckets_orig) / mapped_reads << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Average number of buckets a reverse complement of the read is mapped to: " 
                             << ((float) num_buckets_rev_comp) / mapped_reads << ".\n";
        return std::make_pair(res_orig, res_rev_comp);
    }


    std::vector<std::vector<unsigned int>> _query_file(std::filesystem::path sequence_file) {
        /**
         * * This function is just for benchmarking.
         * @brief Read a query fastq file and output the bucket ids each query belongs to.
         * TODO: include the quality information for fastq.
         */
        std::vector<std::vector<unsigned int>> res;
        seqan3::sequence_file_input<_phred94_traits> fin{sequence_file};
        Timer clock;
        clock.tick();
 
        for (auto & rec : fin) {
            auto [buckets_orig, buckets_rev_comp] = query_sequence(rec.sequence(), rec.base_qualities());
            res.push_back(buckets_orig);
        }
        clock.tock();
        float time = clock.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for bucket query: " 
                             << time << " s (" << time * 1000 * 1000 / res.size() << " μs/seq).\n";
        return res;
    }

    void _check_ground_truth(const std::vector<std::vector<unsigned int>>& query_results, std::filesystem::path ground_truth_file) {
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