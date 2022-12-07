#ifndef FM_INDEX_H
#define FM_INDEX_H

#include <fstream>
 
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
 
#include <cereal/archives/binary.hpp>

class fm_index {
private:
    typedef struct reference_storage
    {
        std::vector<std::string> ids;
        std::vector<std::vector<seqan3::dna5>> seqs;
    } reference_storage_t;

    reference_storage_t* storage;

    void _read_reference(std::filesystem::path const & reference_path)
    {
        seqan3::sequence_file_input reference_in{reference_path};
        for (auto && record : reference_in)
        {
            storage->ids.push_back(record.id());
            storage->seqs.push_back(record.sequence());
        }
    }
    
    void _create_index(std::filesystem::path const & index_path)
    {
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
    fm_index() {
        storage = new reference_storage;
    }

    ~fm_index() {
        delete storage;
    }
    

    bool index(std::filesystem::path const & fasta_file_name, 
               std::filesystem::path const & index_directory,
               std::string const & indicator) {
        /**
         * @brief Read the fasta file, index each bucket.
         * @param fasta_file_name the name of the file containing reference genome.
         * @param index_directory the name of the directory to store the index. It has to be
         *                        either empty directory or not created yet.
         * @param indicator a string indicating the name of output index file.
         * @returns whether we successfully do the indexing.
         */
        // Create directory if directory is not created yet.
        if (!std::filesystem::create_directories(index_directory)) {
            seqan3::debug_stream << "[WARNING]\t" << "The specified index directory "
                                 << index_directory << " is already created." << '\n';
            // count the number 
            unsigned int count = 0;
            for (const auto& entry : std::filesystem::directory_iterator(index_directory)) {
                if (entry.path() == index_directory / indicator + ".fm_index") {
                    seqan3::debug_stream << "[ERROR]\t\t" << "The index file already exists in directory: " 
                                         << index_directory / indicator + ".fm_index" << "." << '\n';
                    return false;
                }
            }
        }
        _read_reference(fasta_file_name);
        _create_index(index_directory / indicator + ".fm_index");
        return true;
    }
};

#endif