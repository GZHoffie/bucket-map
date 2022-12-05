#ifndef BUCKET_MAP_UTILS_H
#define BUCKET_MAP_UTILS_H

#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include <chrono>
#include <thread>
#include <assert.h>

typedef struct config {
    unsigned int bucket_length;
    unsigned int read_length;
} config_t;


template <class DT = std::chrono::milliseconds,
          class ClockT = std::chrono::steady_clock>
class Timer
{
    /**
     * @brief Class for benchmark of runtime.
     *        Adopted from https://coliru.stacked-crooked.com/a/508ec779dea5a28d.
     */
    using timep_t = decltype(ClockT::now());
    
    timep_t _start = ClockT::now();
    timep_t _end = {};

public:
    void tick() { 
        _end = timep_t{};
        _start = ClockT::now(); 
    }
    
    void tock() {
        _end = ClockT::now(); 
    }
    
    template <class duration_t = DT>
    auto duration() const { 
        // Use gsl_Expects if your project supports it.
        assert(_end != timep_t{} && "Timer must toc before reading the time"); 
        return std::chrono::duration_cast<duration_t>(_end - _start); 
    }

    float elapsed_seconds() const {
        return ((float) this->duration().count()) / 1000;
    }
};

void iterate_through_buckets(std::filesystem::path const & fasta_file_name, int bucket_length, int read_length, 
                             std::function<void(std::vector<seqan3::dna4>)> op, bool print_info = false) {
    /**
     * @brief Util function that is used to iterate through all buckets in the reference genome, end execute
     *        some operation on each of the bucket.
     * @param fasta_file_name the path to the reference genome file.
     * @param op an std::function object that takes in a bucket (std::vector<seqan3::dna4>) and returns nothing.
     * @param print_info whether we verbosely print out the bucket information.
     */
    // Read the genome
    seqan3::sequence_file_input reference_genome{fasta_file_name};
    unsigned int _bucket_num = 0;
    for (auto && record : reference_genome) {
        // Divide the record into buckets
        float total_length = (float) record.sequence().size();
        int num_buckets = (int) ceil(total_length / bucket_length);
        if (print_info) {
            seqan3::debug_stream << "[INFO]\t\t" << record.id() << " with length " << (int) total_length
                                 << " divided into " << num_buckets << " buckets.\n";
        }
            
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
            ++_bucket_num;
        }
    }
    if (print_info) {
        seqan3::debug_stream << "[INFO]\t\t" << "Total number of buckets: " 
                             << _bucket_num << "." << '\n';
    }
}




#endif