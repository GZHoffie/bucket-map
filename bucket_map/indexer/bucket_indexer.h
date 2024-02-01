#ifndef BUCKET_MAP_BUCKET_INDEXER_H
#define BUCKET_MAP_BUCKET_INDEXER_H

#include <string>
#include <vector>
#include <cmath>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/all.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/core/debug_stream.hpp>
 
#include <cereal/archives/binary.hpp>

#include "./indexer.h"

template<unsigned int NUM_BUCKETS>
class bucket_indexer : public indexer {
    /**
     * @brief A virtual class of indexer for bucket-map.
     *        General idea: create a separate index file for each bucket.
     */
protected:
    std::vector<std::string> bucket_id;                  // vector storing all bucket info
    std::vector<std::vector<seqan3::dna4>> bucket_seq;   // sequence for each bucket
    unsigned int bucket_length;                          // maximum length of each bucket
    unsigned int read_length;                            // maximum length of each short read
    std::string EXTENSION;                               // file name extension for the perticular indexer

    // parameters for FracMinHash
    std::function<std::size_t(std::size_t)> hash_function;
    unsigned int hash_threshold;

    // q-gram index related information
    std::vector<std::bitset<NUM_BUCKETS>> q_grams_index;
    uint8_t q;

    // fracMinHash related info, each element represent a k-mer.
    // if a k-mer is not sampled, kmer_to_index[k_mer] = -1.
    // otherwise, kmer_to_index[k_mer] = its index in q_grams_index.
    std::vector<int> kmer_to_index;
    unsigned int num_valid_kmers;


    void _insert_into_bucket(const seqan3::bitpacked_sequence<seqan3::dna4>& sequence, unsigned int bucket_num) {
        /**
         * @brief Read the sequence, extract all the q-grams and store in `q_grams_index`.
         * @param sequence the corresponding sequence for the bucket.
         * @param bucket_num an integer indicating the position of the bucket. Should be
         *                   between 0 and NUM_BUCKETS - 1.
         */
        // Extract all q_grams from sequence
        for (auto && value : sequence | seqan3::views::kmer_hash(seqan3::ungapped{q})) {
            int index = kmer_to_index[value];
            if (index >= 0) q_grams_index[index].set(bucket_num);
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


    void _store_q_gram_index(std::filesystem::path const & index_directory, const std::string& indicator) {
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
        if (!check_filename_in(index_directory, indicator + ".qgram")) {
            return;
        }
        // Store the q-gram index in the directory
        std::ofstream index_file(index_directory / (indicator + ".qgram"));
        for (const auto &i : q_grams_index) {
            std::vector<unsigned char> bytes =  _bitset_to_bytes(i);
            index_file.write((char *)&bytes[0], bytes.size());
        }
        index_file.close();
        seqan3::debug_stream << "[INFO]\t\t" << "The bucket q-gram index is stored in: " 
                             << index_directory / (indicator + ".qgram") << ".\n";
    }


    void _store_bucket_ids(std::filesystem::path const & index_directory, const std::string& indicator) {
        if (!check_filename_in(index_directory, indicator + ".bucket_id")) {
            return;
        }
        std::ofstream index_file(index_directory / (indicator + ".bucket_id"));
        for (const auto &i : bucket_id) {
            index_file << i << "\n";
        }
        index_file.close();
        seqan3::debug_stream << "[INFO]\t\t" << "The bucket ids are stored in: " 
                             << index_directory / (indicator + ".bucket_id") << ".\n";
    }

    void _store_sampled_kmers(std::filesystem::path const & index_directory, const std::string& indicator) {
        if (!check_filename_in(index_directory, indicator + ".kmers_index")) {
            return;
        }
        std::ofstream index_file(index_directory / (indicator + ".kmers_index"));
        for (const auto &i : kmer_to_index) {
            index_file << i << "\n";
        }
        index_file.close();
        seqan3::debug_stream << "[INFO]\t\t" << "The kmer indexes are stored in: " 
                             << index_directory / (indicator + ".kmers_index") << ".\n";
    }

    void _init_qgram_index() {
        // initialize q_gram index
        for (int i = 0; i < num_valid_kmers; i++) {
            std::bitset<NUM_BUCKETS> q_gram_bucket;
            q_grams_index.push_back(q_gram_bucket);
        }
    }

public:
    bucket_indexer(unsigned int bucket_len, unsigned int read_len, uint8_t index_seed_length, 
                   std::function<std::size_t(std::size_t)> min_hash_function, unsigned int min_hash_threshold) : indexer() {
        bucket_length = bucket_len;
        read_length = read_len;

        // q-gram index related information
        q = index_seed_length;
        seqan3::debug_stream << "[INFO]\t\t" << "Set index seed length to be: " << q << ".\n";

        // initialize FracMinHash
        int index = 0;
        for (unsigned int i = 0; i < pow(4, q); i++) {
            if (min_hash_function(i) <= min_hash_threshold) {
                kmer_to_index.push_back(index);
                index++;
            } else {
                kmer_to_index.push_back(-1);
            }
        }
        num_valid_kmers = index;
        seqan3::debug_stream << "[INFO]\t\t" << "Number of remaining k-mers after FracMinHash is " << num_valid_kmers << " out of " 
                             << pow(4, q) << " (" << ((float) num_valid_kmers) / pow(4, q) * 100 << "%).\n";
    }

    virtual ~bucket_indexer() = default;

    virtual void create_index(std::filesystem::path const & index_directory) {
        /**
         * @brief Create index files for the buckets and output to the index_directory. 
         * @remark needs to be run after we fill `bucket_id` and `bucket_seq`.
         */
    }

    unsigned int index(std::filesystem::path const & fasta_file_name, 
                       std::filesystem::path const & index_directory,
                       std::string const & indicator) {
        /**
         * @brief Read the fasta file, index each bucket.
         * @param fasta_file_name the name of the file containing reference genome.
         * @param index_directory the name of the directory to store the index. It has to be
         *                        either empty directory or not created yet.
         * @returns the total number of buckets.
         */
        // Create directory if directory is not created yet. Otherwise report the error and return.
        if (!check_extension_in(index_directory, EXTENSION) || !check_filename_in(index_directory, indicator + ".qgram")) {
            return 0;
        }
        Timer clock;
        clock.tick();

        // initialize index
        _init_qgram_index();

        unsigned int bucket_num = 0;
        auto operation = [&](const seqan3::bitpacked_sequence<seqan3::dna4>& seq, const std::string& id) {
            //bucket_seq.push_back(seq);
            bucket_id.push_back(id);
            _insert_into_bucket(seq, bucket_num);
            bucket_num++;
        };
        iterate_through_buckets(fasta_file_name, bucket_length, read_length, operation);
        //create_index(index_directory);

        // create q-gram index files
        _store_q_gram_index(index_directory, indicator);

        seqan3::debug_stream << "[INFO]\t\t" << "The number of buckets: " 
                             << bucket_id.size() << "." << '\n';
        // Store the bucket_id in the directory
        _store_bucket_ids(index_directory, indicator);
        _store_sampled_kmers(index_directory, indicator);

        clock.tock();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for creating and storing index files: " 
                             << clock.elapsed_seconds() << " s.\n";

        // Create and store 
        return bucket_id.size();
    }

    void reset() {
        /**
         * @brief Release the memory storing the sequences and index.
         * 
         */
        bucket_id.clear();
        bucket_id.shrink_to_fit();
        q_grams_index.clear();
        q_grams_index.shrink_to_fit();
        kmer_to_index.clear();
        kmer_to_index.shrink_to_fit();
    }
};


#endif