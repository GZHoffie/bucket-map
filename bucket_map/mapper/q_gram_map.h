#include "../index/bucket_index.h"

#include <string>
#include <vector>
#include <bitset>
#include <cmath>

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

        // initialize q_gram index
        int total_q_grams = (int) pow(4, q);
        for (int i = 0; i < total_q_grams; i++) {
            std::bitset<NUM_BUCKETS> q_gram_bucket;
            q_grams_index.push_back(q_gram_bucket);
        }

        // initialize filter
        filter = new fault_tolerate_filter<NUM_BUCKETS>(num_fault_tolerance);
    }

    ~q_gram_mapper() {
        delete filter;
    }

    bool read(std::filesystem::path const & fasta_file_name, 
              std::filesystem::path const & index_directory) {
        /**
         * @brief Read the fasta file, store the q_gram of each bucket in `q_grams_index`.
         *        Store the information in the `index_directory`.
         * @param fasta_file_name the name of the file containing reference genome.
         * @param index_directory the directory to store the q_gram_count_file.
         * @returns whether the read is successful.
         */
        // Create directory if directory is not created yet.
        if (!std::filesystem::create_directories(index_directory)) {
            seqan3::debug_stream << "[WARNING]\t" << "The specified index directory "
                                 << index_directory << " is already created." << '\n';
            for (const auto& entry : std::filesystem::directory_iterator(index_directory)) {
                if (entry.path().extension() == ".qgram") {
                    seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram file already exists in the " 
                                         << "specified directory. Terminating read." << '\n';
                    return false;
                }
            }
        }
        // Read the genome
        seqan3::sequence_file_input reference_genome{fasta_file_name};
        unsigned int bucket_num = 0;
        for (auto && record : reference_genome) {
            // Divide the record into buckets
            float total_length = (float) record.sequence().size();
            int num_buckets = (int) ceil(total_length / bucket_length);
            seqan3::debug_stream << "[INFO]\t\t" << record.id() << " with length " << (int) total_length
                                 << " divided into " << num_buckets << " buckets.\n";
            
            // read each bucket
            for (int i = 0; i < num_buckets; i++) {
                int start = i * bucket_length;
                int end = start + bucket_length + read_length;
                if (end > record.sequence().size()) {
                    end = record.sequence().size();
                }
                std::vector<seqan3::dna4> bucket_sequence(&record.sequence()[start], &record.sequence()[end]);
                _insert_into_bucket(bucket_sequence, bucket_num);
                bucket_num++;
            }
        }
        seqan3::debug_stream << "[INFO]\t\t" << "Total number of buckets: " 
                             << bucket_num << "." << '\n';
        
        // Store the q-gram index in the directory
        // TODO: Implement a more space-efficient way to store and load the q-gram index
        std::ofstream file(index_directory / "index.qgram");
        for (const auto &i : q_grams_index) {
            file << i << "\n";
        }
        seqan3::debug_stream << "[INFO]\t\t" << "The bucket q-gram index is stored in: " 
                             << index_directory / "index.qgram" << "." << '\n';
        return true;
    }


    std::vector<int> query(std::vector<seqan3::dna5> sequence) {
        /**
         * @brief From `q_grams_index`, determine where the sequence may be coming from.
         * @param sequence the DNA sequence to be mapped.
         * @returns a vector of integers indicating the possible regions that the sequence
         *          may belong to.
         */
        // Reset the filter
        filter->reset();

        // Extract the k-mers from sequence
        auto q_gram = seqan3::views::kmer_hash(q_gram_shape);

        // Extract all q_grams from sequence
        auto hash_values = sequence | q_gram;
        std::vector<int> samples;
        std::sample(hash_values.begin(), hash_values.end(), 
                    std::back_inserter(samples),
                    5, std::mt19937{std::random_device{}()});

        // insert those samples into the filter
        for (int h : samples) {
            filter->read(q_grams_index[h]);
        }
        return filter->best_results();

    }

};