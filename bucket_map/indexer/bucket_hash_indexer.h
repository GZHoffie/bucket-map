#ifndef BUCKET_MAP_BUCKET_HASH_INDEXER_H
#define BUCKET_MAP_BUCKET_HASH_INDEXER_H

#include "./bucket_indexer.h"

//#include <dlib/serialize.h>

template<unsigned int NUM_BUCKETS>
class bucket_hash_indexer : public bucket_indexer<NUM_BUCKETS> {
private:
    seqan3::shape index_shape;

public:
    bucket_hash_indexer(unsigned int bucket_len, unsigned int read_len, 
                        seqan3::shape bucket_shape, seqan3::shape locate_shape) : bucket_indexer<NUM_BUCKETS>(bucket_len, read_len, bucket_shape) {
        /**
         * @param bucket_shape the shape of q-gram for the mapper: to map reads to their buckets
         * @param locate_shape the shape of q-gram for the locator: to map reads to exact locations.
         */
        this->EXTENSION = ".bhi";
        index_shape = locate_shape;
    }

    void create_index(std::filesystem::path const & index_directory) {
        /**
         * @brief Create index files for the buckets. 
         * @remark needs to be run after we fill `bucket_id` and `bucket_seq`.
         */
        std::unordered_multimap<unsigned int, int> index;
        for (int i = 0; i < this->bucket_id.size(); i++) {
            index.clear();
            // insert the k-mers into the hash table
            auto values = this->bucket_seq[i] | seqan3::views::kmer_hash(index_shape);
            int offset = 0;
            for (auto hash : values) {
                // Only record the last appearance of the k-mer
                index.emplace(hash, offset);
                offset++;
            }
            std::ofstream os{index_directory / (std::to_string(i) + this->EXTENSION)};
            //dlib::serialize(index, os);
            os.close();
        }
    }
};


#endif