#ifndef BUCKET_MAP_FM_INDEX_H
#define BUCKET_MAP_FM_INDEX_H

#include <fstream>
 
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
 
#include <cereal/archives/binary.hpp>

#include "./indexer.h"

class fm_indexer : public indexer {
    /**
     * @brief The default indexer using FM-index (Burrows-Wheeler Transform).
     */

private:
    typedef struct reference_storage
    {
        std::vector<std::string> ids;
        std::vector<std::vector<seqan3::dna5>> seqs;
    } reference_storage_t;

    reference_storage_t* storage;

    void _read_reference(std::filesystem::path const & reference_path)
    {
        /**
         * @brief Read the reference genome (fasta file) and store the sequences in
         *        `storage` variable.
         */
        seqan3::sequence_file_input reference_in{reference_path};
        for (auto && record : reference_in)
        {
            storage->ids.push_back(record.id());
            storage->seqs.push_back(record.sequence());
        }
    }

    void _create_index(std::filesystem::path const & index_path)
    {
        /**
         * @brief Create the index file `index_path` from the `storage` variable.
         * @note must be called after `_read_reference` function.
         */
        seqan3::bi_fm_index index{storage->seqs};
        {
            std::ofstream os{index_path, std::ios::binary};
            cereal::BinaryOutputArchive oarchive{os};
            oarchive(index);
        }
        seqan3::debug_stream << "[INFO]\t\t" << "The FM-index file is stored in: " 
                             << index_path << "." << '\n';
    }
    

public:
    fm_indexer() : indexer() {
        storage = new reference_storage;
    }

    ~fm_indexer() {
        delete storage;
    }
    
    unsigned int index(std::filesystem::path const & fasta_file_name, 
                       std::filesystem::path const & index_directory,
                       std::string const & indicator) {
        // Create directory if directory is not created yet.
        if (check_filename_in(index_directory, indicator + ".fm_index")) {
            _read_reference(fasta_file_name);
            _create_index(index_directory / (indicator + ".fm_index"));
            return 1;
        }
        return 0;
    }
};

#endif