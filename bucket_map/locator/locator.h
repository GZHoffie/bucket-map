#ifndef BUCKET_MAP_LOCATOR_H
#define BUCKET_MAP_LOCATOR_H

#include "../indexer/indexer.h"
#include "../mapper/mapper.h"
#include "../utils.h"


class locator {
    /**
     * @brief Virtual class for all locator that locates the *exact locations* of each
     *        short read.
     */
private:
    indexer* _i;

public:
    locator(indexer* ind) {
        _i = ind;
    }

    unsigned int initialize(std::filesystem::path const & fasta_file_name, 
                            std::filesystem::path const & index_directory,
                            std::string const & indicator) {
        /**
         * @brief Read the fasta file, and create index file in the index directory.
         * @param fasta_file_name the name of the file containing reference genome.
         * @param index_directory the name of the directory to store the index. It has to be
         *                        either empty directory or not created yet.
         * @param indicator a string indicating the name of output index file.
         * @returns how many index files we are generating.
         */
        return _i->index(fasta_file_name, index_directory, indicator);
    }
    
    virtual void locate(std::filesystem::path const & sequence_file, 
                        std::filesystem::path const & index_path,
                        std::filesystem::path const & sam_file) {
        /**
         * @brief Find the exact location of reads in a fastq file, and output
         *        the mapping result to a sam file.
         * @param sequence_file the path to the fastq file containing query reads.
         * @param index_file the path to the index file.
         * @param sam_file the path to output the sam file.
         */
    }
};


#endif