#ifndef BUCKET_MAP_BUCKET_FM_INDEX_H
#define BUCKET_MAP_BUCKET_FM_INDEX_H

#include "./bucket_indexer.h"

class bucket_fm_indexer : public bucket_indexer {
public:
    bucket_fm_indexer(unsigned int bucket_len, unsigned int read_len) : bucket_indexer(bucket_len, read_len) {
        EXTENSION = ".bfmi";
    }

    void create_index(std::filesystem::path const & index_directory) {
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