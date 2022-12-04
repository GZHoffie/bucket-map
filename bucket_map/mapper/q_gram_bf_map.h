#ifndef Q_GRAM_BF_MAP_H
#define Q_GRAM_BF_MAP_H

#include "./q_gram_map.h"

class q_gram_bloom_filter_mapper : public q_gram_mapper {
private:
    unsigned int K;
    unsigned int N;


public:
    q_gram_bloom_filter_mapper(unsigned int bucket_len, 
                               unsigned int read_len, 
                               seqan3::shape shape, 
                               unsigned int samples, 
                               unsigned int fault,
                               unsigned int k,
                               unsigned int n) : 
        q_gram_mapper(bucket_len, read_len, shape, samples, fault) {
        /**
         * @brief Variant of Q-gram mapper implemented using bloom filter.
         * @remark The `q_gram_mapper` class can be seen as a special case of this class
         *         where K=1 and the hash function being a direct addressing function.
         * @param k the number of hash functions used in the bloom filter.
         * @param n the number of bits used in each bloom filter.
         * @note In this bloom filter, m = `bucket_len` + `read_len`, k = `K`, n = `N`.
         *       The false positive rate is likely higher than that in the DAT.
         *       The hash functions are created using prime field method.
         */
        K = k;
        N = n;
        // initialize hash functions
        
    
    
    }

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