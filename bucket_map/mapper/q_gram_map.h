#include "../index/bucket_index.h"
#include "../utils.h"

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
class q_gram_mapper {
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


    void _insert_into_bucket(std::vector<seqan3::dna4> sequence, unsigned int bucket_num) {
        /**
         * @brief Read the sequence, extract all the q-grams and store in `q_grams_index`.
         * @param sequence the corresponding sequence for the bucket.
         * @param bucket_num an integer indicating the position of the bucket. Should be
         *                   between 0 and NUM_BUCKETS - 1.
         */
        auto q_gram = seqan3::views::kmer_hash(q_gram_shape);

        // Extract all q_grams from sequence
        for (auto && value : sequence | q_gram) {
            q_grams_index[value].set(bucket_num);
        }
    }

    std::vector<unsigned char> _bitset_to_bytes(const std::bitset<NUM_BUCKETS>& bs){
        /**
         * @brief Convert bitset into 8-byte chars.
         * Adopted from https://stackoverflow.com/a/7463972.
         */
        std::vector<unsigned char> result((NUM_BUCKETS + 7) >> 3);
        for (unsigned int j = 0; j < NUM_BUCKETS; j++)
            result[j>>3] |= (bs[j] << (j & 7));
        return result;
    }

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
                  unsigned int samples, unsigned int fault) {
        // initialize private variables
        bucket_length = bucket_len;
        read_length = read_len;
        
        q_gram_shape = shape;
        size = std::ranges::size(shape);
        q = shape.count();
        seqan3::debug_stream << "[INFO]\t\t" << "Set q-gram shape to be: " 
                             << shape << " with number of effective characters: " << q << '\n';
        
        num_samples = samples;
        num_fault_tolerance = fault;

        // initialize filter
        filter = new fault_tolerate_filter<NUM_BUCKETS>(num_fault_tolerance);
    }

    ~q_gram_mapper() {
        delete filter;
    }

    void read(std::filesystem::path const & fasta_file_name) {
        /**
         * @brief Read the fasta file, store the q_gram of each bucket in `q_grams_index`.
         *        Store the information in the `index_directory`.
         * @param fasta_file_name the name of the file containing reference genome.
         */
        // q_grams_index should be empty
        if (!q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is not empty. Terminating read.\n";
            return;
        }
        Timer clock;
        clock.tick();
        // initialize q_gram index
        int total_q_grams = (int) pow(4, q);
        for (int i = 0; i < total_q_grams; i++) {
            std::bitset<NUM_BUCKETS> q_gram_bucket;
            q_grams_index.push_back(q_gram_bucket);
        }
        // Read the genome
        unsigned int bucket_num = 0;
        auto operation = [&](std::vector<seqan3::dna4> seq) {
            _insert_into_bucket(seq, bucket_num);
            bucket_num++;
        };
        iterate_through_buckets(fasta_file_name, bucket_length, read_length, operation);
        
        clock.tock();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for reading: " 
                             << clock.elapsed_seconds() << " s." << '\n';
        seqan3::debug_stream << "[INFO]\t\t" << "Total number of buckets: " 
                             << bucket_num << "." << '\n';
    }

    void store(std::filesystem::path const & index_directory) {
        /**
         * @brief Store the created index inside the index_directory.
         * @param index_directory the directory to store the q_gram_count_file.
         */
        // q_grams_index should not be empty
        if (q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is empty. No file is created.\n";
            return;
        }
        // Create directory if directory is not created yet.
        // Return if the index files already exist.
        if (!std::filesystem::create_directories(index_directory)) {
            seqan3::debug_stream << "[WARNING]\t" << "The specified index directory "
                                 << index_directory << " is already created." << '\n';
            for (const auto& entry : std::filesystem::directory_iterator(index_directory)) {
                if (entry.path().extension() == ".qgram" || entry.path().extension() == ".pattern") {
                    seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram file " << entry.path() << " already exists" 
                                         << " in the specified directory. Terminating store." << '\n';
                    return;
                }
            }
        }
        // Store the q-gram index in the directory
        Timer clock;
        clock.tick();
        std::ofstream index_file(index_directory / "index.qgram");
        for (const auto &i : q_grams_index) {
            std::vector<unsigned char> bytes =  _bitset_to_bytes(i);
            index_file.write((char *)&bytes[0], bytes.size());
        }
        clock.tock();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for storing index files: " 
                             << clock.elapsed_seconds() << " s." << '\n';
        seqan3::debug_stream << "[INFO]\t\t" << "The bucket q-gram index is stored in: " 
                             << index_directory / "index.qgram" << ".\n";
        
        // Store the q-gram pattern in the directory
        std::ofstream pattern_file(index_directory / "index.pattern");
        pattern_file << q_gram_shape;
        seqan3::debug_stream << "[INFO]\t\t" << "The q-gram shape is stored in: " 
                             << index_directory / "index.pattern" << "." << '\n';
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
            for (unsigned int i = 0; i < total_q_grams; i++) {
                unsigned int start_index = i * num_chars_per_q_gram;
                unsigned int end_index = (i+1) * num_chars_per_q_gram;
                auto data = new std::vector<unsigned char>(buffer + start_index, buffer + end_index);
                q_grams_index.push_back(_bitset_from_bytes(*data));
                delete data;
            }
            delete[] buffer;
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
        auto hash_values = sequence | q_gram;
        std::vector<int> samples;
        // Randomly sample `num_samples` q-grams for query.
        std::sample(hash_values.begin(), hash_values.end(), 
                    std::back_inserter(samples), num_samples,
                    std::mt19937{std::random_device{}()});
        return query(samples);
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