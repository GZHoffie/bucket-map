#ifndef BUCKET_MAP_H
#define BUCKET_MAP_H

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

#include "../utils.h"


class bucket_indexer {
private:
    std::vector<std::string> bucket_id;                  // vector storing all bucket info
    std::vector<std::vector<seqan3::dna4>> bucket_seq;   // sequence for each bucket
    unsigned int bucket_length;                          // maximum length of each bucket
    unsigned int read_length;                            // maximum length of each short read

public:
    bucket_indexer(unsigned int bucket_len, unsigned int read_len) {
        bucket_length = bucket_len;
        read_length = read_len;
    }

    virtual ~bucket_indexer() {}

    virtual void _create_index(std::filesystem::path const & index_directory) {
        /**
         * @brief Create index files for the buckets and output to the index_directory. 
         * @remark needs to be run after we fill `bucket_id` and `bucket_seq`.
         */
    }

    unsigned int index(std::filesystem::path const & fasta_file_name, 
                       std::filesystem::path const & index_directory) {
        /**
         * @brief Read the fasta file, index each bucket.
         * @param fasta_file_name the name of the file containing reference genome.
         * @param index_directory the name of the directory to store the index. It has to be
         *                        either empty directory or not created yet.
         * @returns the total number of buckets.
         */
        // Create directory if directory is not created yet.
        if (!std::filesystem::create_directories(index_directory)) {
            seqan3::debug_stream << "[WARNING]\t" << "The specified index directory "
                                 << index_directory << " is already created." << '\n';
            // count the number 
            unsigned int count = 0;
            for (const auto& entry : std::filesystem::directory_iterator(index_directory)) {
                if (entry.path().extension() == ".index")
                    ++count;
            }
            if (count > 0) {
                seqan3::debug_stream << "[INFO]\t\t" << "The number of index files in directory: " 
                                     << count << "." << '\n';
                return count;
            }
        }
        auto operation = [&](std::vector<seqan3::dna4> seq) {
            bucket_seq.push_back(seq);
        };
        iterate_through_buckets(fasta_file_name, bucket_length, read_length, operation);
        _create_index(index_directory);
        seqan3::debug_stream << "[INFO]\t\t" << "The number of index files created: " 
                             << bucket_id.size() << "." << '\n';
        
        // Store the bucket_id in the directory
        {
            std::ofstream os{index_directory / "bucket_id", std::ios::binary};
            cereal::BinaryOutputArchive oarchive{os};
            oarchive(bucket_id);
        }
        seqan3::debug_stream << "[INFO]\t\t" << "The bucket id are stored in: " 
                             << index_directory / "bucket_id" << "." << '\n';
        return bucket_id.size();
    }
};


#endif