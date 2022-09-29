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

#include "./bucket_index.h"



class bucket_fm_indexer : public bucket_indexer {
private:
    std::vector<std::string> bucket_id;                  // vector storing all bucket info
    std::vector<std::vector<seqan3::dna4>> bucket_seq;   // sequence for each bucket
    unsigned int bucket_length;                          // maximum length of each bucket
    unsigned int read_length;                            // maximum length of each short read



public:
    bucket_fm_indexer(unsigned int bucket_len, unsigned int read_len) : bucket_indexer(bucket_len, read_len) {}

    void _create_index(std::filesystem::path const & index_directory) {
        /**
         * @brief Create index files for the buckets. 
         * @remark needs to be run after we fill `bucket_id` and `bucket_seq`.
         */
        for (int i = 0; i < bucket_id.size(); i++) {
            seqan3::bi_fm_index index{bucket_seq[i]};
            // store the index in a file
            {
                std::ofstream os{index_directory / (std::to_string(i) + ".bfmi"), std::ios::binary};
                cereal::BinaryOutputArchive oarchive{os};
                oarchive(index);
            }
        }
    }
};


#endif