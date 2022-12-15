#ifndef BUCKET_MAP_BUCKET_LOCATOR_H
#define BUCKET_MAP_BUCKET_LOCATOR_H

#include "./locator.h"

#include <seqan3/search/fm_index/bi_fm_index.hpp>

class bucket_locator : public locator {
private:
    mapper* _m;
    std::vector<std::unordered_multimap<int, int>> index;
    unsigned int bucket_length;
    unsigned int read_length;
    seqan3::shape q_gram_shape;

    void _initialize_kmer_index(std::filesystem::path const & fasta_file_name) {
        auto operation = [&](std::vector<seqan3::dna4> seq, std::string id) {
            unsigned int offset = 0;
            auto values = seq | seqan3::views::kmer_hash(q_gram_shape);
            std::unordered_multimap<int, int> bucket_index;
            for (auto && hash : values) {
                bucket_index.insert(std::make_pair(hash, offset));
                offset++;
            }
            index.push_back(bucket_index);
        };
        iterate_through_buckets(fasta_file_name, bucket_length, read_length, operation);
    }


public:
    bucket_locator(indexer* ind, mapper* map, unsigned int bucket_len, unsigned int read_len, seqan3::shape shape) : locator(ind) {
        _m = map;
        bucket_length = bucket_len;
        read_length = read_len;
        q_gram_shape = shape;
    }

    unsigned int initialize(std::filesystem::path const & fasta_file_name, 
                            std::filesystem::path const & index_directory,
                            std::string const & indicator) {

        // load q-gram index to the mapper
        _m->load(index_directory);

        // load all index files
        _initialize_kmer_index(fasta_file_name);

        // create index files
        return locator::initialize(fasta_file_name, index_directory, indicator);
    }

    void locate(std::filesystem::path const & sequence_file, 
                std::filesystem::path const & index_file,
                std::filesystem::path const & sam_file) {
        
    }
};

#endif