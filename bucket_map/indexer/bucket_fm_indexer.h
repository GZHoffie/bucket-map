#ifndef BUCKET_MAP_BUCKET_FM_INDEX_H
#define BUCKET_MAP_BUCKET_FM_INDEX_H

#include "./bucket_indexer.h"

template<unsigned int NUM_BUCKETS>
class bucket_fm_indexer : public bucket_indexer<NUM_BUCKETS> {
public:
    bucket_fm_indexer(unsigned int bucket_len, unsigned int read_len, seqan3::shape shape) : bucket_indexer<NUM_BUCKETS>(bucket_len, read_len, shape) {
        this->EXTENSION = ".bfmi";
    }

    void create_index(std::filesystem::path const & index_directory) {
        /**
         * @brief Create index files for the buckets. 
         * @remark needs to be run after we fill `bucket_id` and `bucket_seq`.
         */
        std::ofstream os{index_directory / ("index" + this->EXTENSION), std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        for (int i = 0; i < this->bucket_id.size(); i++) {
            seqan3::bi_fm_index index{this->bucket_seq[i]};
            oarchive(index);
        }
    }
};


#endif