#include "bucket_index.h"

#include <string>
#include <vector>
#include <bitset>
#include <cmath>

#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

using seqan3::operator""_shape;

template<unsigned int NUM_BUCKETS>
class q_gram_mapper {
private:
    seqan3::shape q_gram_shape;
    std::vector<std::bitset<NUM_BUCKETS>> q_grams_index;
    unsigned int q;
    unsigned int size;

    unsigned int bucket_length;
    unsigned int read_length;


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
    q_gram_mapper(unsigned int bucket_len, unsigned int read_len, seqan3::shape shape) {
        bucket_length = bucket_len;
        read_length = read_len;
        
        q_gram_shape = shape;
        size = std::ranges::size(shape);
        q = shape.count();
        seqan3::debug_stream << "[INFO]\t\t" << "Set q-gram shape to be: " 
                             << shape << " with number of effective characters: " << q << '\n';

        // initialize q_gram index
        int total_q_grams = (int) pow(4, q);
        for (int i = 0; i < total_q_grams; i++) {
            std::bitset<NUM_BUCKETS> q_gram_bucket;
            q_grams_index.push_back(q_gram_bucket);
        }
    }

    bool read(std::filesystem::path const & fasta_file_name, 
              std::filesystem::path const & q_gram_count_directory) {
        /**
         * @brief Read the fasta file, store the q_gram of each bucket in `q_grams_index`.
         *        Store the information in the `q_gram_count_directory`.
         * @param fasta_file_name the name of the file containing reference genome.
         * @param q_gram_count_directory the directory to store the q_gram_count_file.
         * @returns whether the read is successful.
         */
        // Create directory if directory is not created yet.
        if (!std::filesystem::create_directories(q_gram_count_directory)) {
            seqan3::debug_stream << "[WARNING]\t" << "The specified index directory "
                                 << q_gram_count_directory << " is already created." << '\n';
            for (const auto& entry : std::filesystem::directory_iterator(q_gram_count_directory)) {
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
        std::ofstream file(q_gram_count_directory / "index.qgram");
        for (const auto &i : q_grams_index) {
            file << i << "\n";
        }
        seqan3::debug_stream << "[INFO]\t\t" << "The bucket q-gram index is stored in: " 
                             << q_gram_count_directory / "index.qgram" << "." << '\n';
        return true;
    }

};