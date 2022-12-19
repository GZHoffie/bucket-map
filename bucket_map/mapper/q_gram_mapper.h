#ifndef BUCKET_MAP_Q_GRAM_MAP_H
#define BUCKET_MAP_Q_GRAM_MAP_H


#include "../indexer/bucket_fm_indexer.h"
#include "../utils.h"
#include "./mapper.h"

#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <chrono>
#include <fstream>

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
    std::vector<std::bitset<NUM_BUCKETS>> filters;

    std::vector<int> _set_bits(int index) {
        /**
         * @brief Find all set bits in filters[index]
         * @param index indicates which filter we want to check.
         */
        std::vector<int> res;
        if (index >= num_fault_tolerance) {
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
    fault_tolerate_filter(unsigned int fault) {
        num_fault_tolerance = fault;
        for (int i = 0; i < num_fault_tolerance; i++) {
            std::bitset<NUM_BUCKETS> filter;
            filters.push_back(filter.set());
        }
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
};

template<unsigned int NUM_BUCKETS>
class distinguishability_filter {
    /**
     * @brief A filter that filter out Q-grams that appear in most of the buckets.
     */
private:
    unsigned int threshold;
    std::vector<unsigned int> zeros;

public:
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
            zeros.push_back(NUM_BUCKETS - ones);
        }
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of Q-grams with distinguishability >= " << ((float) threshold) / NUM_BUCKETS << ": " 
                             << valid_q_grams << " (" << ((float) valid_q_grams) / q_grams_index.size() * 100 << "%).\n";
    }

    std::function<bool(unsigned int)> get_filter() {
        return [&](unsigned int kmer_hash) {
            return zeros[kmer_hash] >= threshold;
        };
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

    // filter that filter out the most possible bucket
    fault_tolerate_filter<NUM_BUCKETS>* filter;

    // Q-gram filters for map efficiency
    distinguishability_filter<NUM_BUCKETS>* dist_filter;
    std::function<bool(unsigned int)> dist_view;
    //TODO: also consider quality


    std::bitset<NUM_BUCKETS> _bitset_from_bytes(const std::vector<unsigned char>& buf) {
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
                  unsigned int samples, unsigned int fault, float distinguishability) : mapper() {
        // initialize private variables
        bucket_length = bucket_len;
        read_length = read_len;
        
        q_gram_shape = shape;
        size = std::ranges::size(shape);
        q = shape.count();
        
        num_samples = samples;
        num_fault_tolerance = fault;

        // initialize filter
        filter = new fault_tolerate_filter<NUM_BUCKETS>(num_fault_tolerance);
        dist_filter = new distinguishability_filter<NUM_BUCKETS>(distinguishability);
    }

    ~q_gram_mapper() {
        delete filter;
        delete dist_filter;
    }


    void load(std::filesystem::path const & index_directory) {
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
        std::ifstream file(index_directory / "index.qgram");
        if (file) {
            Timer clock;
            clock.tick();
            // get length fo file
            file.seekg(0, file.end);
            int length = file.tellg();
            file.seekg(0, file.beg);
            
            // read the file into buffer
            char * buffer = new char[length];
            file.read(buffer, length);
            
            // load buffer into q_grams_index
            unsigned int start_index, end_index;
            for (unsigned int i = 0; i < total_q_grams; i++) {
                start_index = i * num_chars_per_q_gram;
                end_index = (i+1) * num_chars_per_q_gram;
                auto data = new std::vector<unsigned char>(buffer + start_index, buffer + end_index);
                q_grams_index.push_back(_bitset_from_bytes(*data));
                delete data;
            }
            delete[] buffer;
            dist_filter->read(q_grams_index);
            dist_view = dist_filter->get_filter();

            // Complete the read
            clock.tock();
            seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for loading index files: " 
                                 << clock.elapsed_seconds() << " s." << '\n';
            seqan3::debug_stream << "[INFO]\t\t" << "Successfully loaded " 
                                 << index_directory / "index.qgram" << "." << '\n';
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
    }

    std::vector<int> query_sequence(const std::vector<seqan3::dna4>& sequence) {
        /**
         * @brief From `q_grams_index`, determine where the sequence may be coming from.
         * @param sequence A dna4 vector containing the query sequence.
         * @returns a vector of integers indicating the possible regions that the sequence
         *          may belong to.
         * TODO: implement query with argument being a file containing short reads, together with quality.
         */
        auto q_gram = seqan3::views::kmer_hash(q_gram_shape);
        auto values = sequence | q_gram | std::views::filter(dist_view);
        std::vector<unsigned int> hash_values(values.begin(), values.end());

        // if not enough q-grams ramained to determine the exact location, simply ignore this query sequence.
        if (hash_values.size() < 0.5 * num_samples){
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
        float delta;
        if (num_samples == 1) {
            delta = 0;
        } else {
            delta = static_cast<float>(hash_values.size()-1) / (num_samples-1);
        }
        //std::cout << hash_values.size() << " " << delta << std::endl;
        for (int i = 0; i < num_samples; i++) {
            samples.push_back(hash_values[floor(i*delta)]);
        }
        return query(samples);
    }



    std::vector<std::vector<unsigned int>> map(std::filesystem::path const & sequence_file) {
        /**
         * @brief Read a query fastq file and output the ids of the sequence that are mapped 
         *        to each file.
         */

        std::vector<std::vector<unsigned int>> res;
        seqan3::sequence_file_input<_dna4_traits> fin{sequence_file};
        // initialize returning result
        for (int i = 0; i < NUM_BUCKETS; i++) {
            std::vector<unsigned int> sequence_ids;
            res.push_back(sequence_ids);
        }

        Timer clock;
        clock.tick();

        unsigned int index = 0;
        for (auto & rec : fin) {
            seqan3::debug_stream << index << ": " << rec << "\n";
            records.push_back(rec.sequence());
            for (auto & bucket : query_sequence(rec.sequence())) {
                res[bucket].push_back(index);
            }
            ++index;
        }
        clock.tock();
        float time = clock.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for bucket mapping: " 
                             << time << " s (" << time * 1000 * 1000 / res.size() << " μs/seq).\n";
        return res;
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

        // Test statistics
        int correct_map = 0;
        int total_bucket_numbers = 0;
        std::map<int, int> bucket_number_map;

        for (int i = 0; i < query_results.size(); i++) {
            is >> bucket >> exact_location;
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
};


#endif