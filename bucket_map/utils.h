#ifndef BUCKET_MAP_UTILS_H
#define BUCKET_MAP_UTILS_H

#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/alphabet/all.hpp>

typedef struct config {
    unsigned int bucket_length;
    unsigned int read_length;
} config_t;


void iterate_through_buckets(std::filesystem::path const & fasta_file_name, int bucket_length, int read_length, 
                             std::function<void(std::vector<seqan3::dna4>)> op) {
    // Read the genome
    seqan3::sequence_file_input reference_genome{fasta_file_name};
    unsigned int bucket_num = 0;
    for (auto && record : reference_genome) {
        // Divide the record into buckets
        float total_length = (float) record.sequence().size();
        int num_buckets = (int) ceil(total_length / bucket_length);
        seqan3::debug_stream << "[INFO]\t\t" << record.id() << " with length " << (int) total_length
                             << " divided into " << num_buckets << " buckets.\n";
            
        // read each bucket
        for (int i = 0; i < num_buckets; i++) {
            int start = i * bucket_length;
            int end = start + bucket_length + read_length;
            if (end > record.sequence().size()) {
                end = record.sequence().size();
            }
            if (end - start <= read_length) {
                continue;
            }

            std::vector<seqan3::dna4> bucket_sequence(&record.sequence()[start], &record.sequence()[end]);
            op(bucket_sequence);
        }
    }
    seqan3::debug_stream << "[INFO]\t\t" << "Total number of buckets: " 
                         << bucket_num << "." << '\n';
}

#endif