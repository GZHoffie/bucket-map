#ifndef BUCKET_MAP_BUCKET_LOCATOR_H
#define BUCKET_MAP_BUCKET_LOCATOR_H

#include "./locator.h"

class bucket_locator : public locator {
private:
    mapper* _m;

    void _initialize_mapper(std::filesystem::path const & fasta_file_name);

public:
    bucket_locator(indexer* ind, mapper* map) : locator(ind) {
        _m = map;
    }



    void locate(std::filesystem::path const & sequence_file, 
                std::filesystem::path const & index_file,
                std::filesystem::path const & sam_file) {
    }
};

#endif