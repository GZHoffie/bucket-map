#ifndef BUCKET_MAP_H
#define BUCKET_MAP_H

#include <string>
#include <vector>

#include <seqan3/search/fm_index/bi_fm_index.hpp>

/**
 * @brief A class storing the information of a bucket, namely the 
 * sequence that the bucket belongs to, as well as the index of bucket.
 */
class bucket_info {
    std::string sequence_id;
    unsigned int bucket_id;
};


class bucket_mapper {
private:
    std::string index_id;                          // a string indicating the handle to index files
    std::vector<bucket_info> buckets;              // vector storing all bucket info
    std::vector<seqan3::bi_fm_index> bucket_index; // index files for each buckets
    unsigned int bucket_length;                    // maximum length of each bucket

public:
    bucket_mapper(int length) {
        bucket_length = length;
    }

    int index(std::string fasta_file_name, bool output_index_file) {
        /**
         * @brief Read the fasta file, index each bucket and store in
         * `bucket_index`.
         * @param fasta_file_name the name of the file containing reference genome.
         * @param output_index_file whether we output the index file using cereal.
         * @returns the total number of buckets.
         */
        seqan3::sequence_file_input reference_genome{fasta_file_name};
        std::string previous_id = "";
        int index = 0;
        for (auto && record : reference_genome) {
            if (previous_id != record.id()) {
                previous_id = record.id();
                index = 0;
            }
            buckets.push_back(new bucket_info(record.id(), index));
            index++;
        }
    }

};


#endif