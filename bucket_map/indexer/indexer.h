#ifndef BUCKET_MAP_INDEXER_H
#define BUCKET_MAP_INDEXER_H

#include <seqan3/core/debug_stream.hpp>

#include "../utils.h"

class indexer {
    /**
     * @brief Virtual class for all indexer, which produce the indexing files that
     *        finds the exact locations of reads.
     */
public:
    indexer() {}

    /**
     * @brief Read the fasta file, and create index file in the index directory.
     * @param fasta_file_name the name of the file containing reference genome.
     * @param index_directory the name of the directory to store the index. It has to be
     *                        either empty directory or not created yet.
     * @param indicator a string indicating the name of output index file.
     * @returns how many index files we are generating.
     */
    virtual unsigned int index(std::filesystem::path const & fasta_file_name, 
                               std::filesystem::path const & index_directory,
                               std::string const & indicator) = 0;
        
    
};

#endif