#ifndef BUCKET_MAP_UTILS_H
#define BUCKET_MAP_UTILS_H

#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>

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
    
    // variable storing start and end time.
    timep_t _start = ClockT::now();
    timep_t _end = {};

    // total duration between each (_start, _end) pair
    DT _elapsed = DT::zero();

public:
    auto duration() const { 
        // Use gsl_Expects if your project supports it.
        assert(_end != timep_t{} && "Timer must toc before reading the time"); 
        return _elapsed; 
    }

    void tick() { 
        _end = timep_t{};
        _start = ClockT::now(); 
    }
    
    void tock() {
        _end = ClockT::now();
        _elapsed += std::chrono::duration_cast<DT>(_end - _start);
    }
    
    float elapsed_seconds() const {
        return ((float) _elapsed.count()) / 1000;
    }
};

void iterate_through_buckets(std::filesystem::path const & fasta_file_name, int bucket_length, int read_length, 
                             std::function<void(std::vector<seqan3::dna4>, std::string)> op, bool print_info = false) {
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
            op(bucket_sequence, record.id());
            ++_bucket_num;
        }
    }
    if (print_info) {
        seqan3::debug_stream << "[INFO]\t\t" << "Total number of buckets: " 
                             << _bucket_num << "." << '\n';
    }
}

bool check_extension_in(std::filesystem::path const & index_directory,
                            std::string ext) {
    /**
     * @brief Check if the specified directory exists. If it does, check if files with certain
     *        file extension name exist. If not, create the directory.
     * @param index_directory a path to a directory to be checked.
     * @param ext a string starting with a ".", indicating file extension name.
     * @returns true if the directory doesn't exist or no such file extension found. False otherwise.
     */
    if (!std::filesystem::create_directories(index_directory)) {
        for (const auto& entry : std::filesystem::directory_iterator(index_directory)) {
            if (entry.path().extension() == ext) {
                seqan3::debug_stream << "[ERROR]\t\t" << "The file with extension " << ext <<" already exists in directory: " 
                                     << entry.path() << "." << '\n';
                return false;
            }
        }    
    }
    return true;    
}

bool check_filename_in(std::filesystem::path const & index_directory,
                       std::string filename) {
    /**
     * @brief Check if the specified directory exists. If it does, check if files with certain
     *        name exist. If not, create the directory.
     * @param index_directory a path to a directory to be checked.
     * @param filename a string, including the extension name.
     * @returns true if the directory doesn't exist or no such file found. False otherwise.
     */
    if (!std::filesystem::create_directories(index_directory)) {
        for (const auto& entry : std::filesystem::directory_iterator(index_directory)) {
            if (entry.path() == index_directory / filename) {
                seqan3::debug_stream << "[ERROR]\t\t" << "The specified file already exists in directory: " 
                                     << index_directory / filename << "." << '\n';
                return false;
            }
        }    
    }
    return true;    
}

class Sampler {
private:
    // stored values for sample efficiency
    unsigned int last_upper_bound;
    unsigned int n;

public:
    std::vector<unsigned int> samples;

    Sampler(unsigned int num_samples) {
        last_upper_bound = 0;
        n = num_samples;
    }

    void sample_deterministically(unsigned int upper_bound) {
        /**
         * @brief Sample n numbers deterministically and uniformly between [0, upper_bound].
         * @param upper_bound the maximum number that can be sampled.
         */
        if (last_upper_bound == upper_bound) return;
        // do a resampling
        samples.clear();
        double delta; // difference between two samples
        if (n == 1) {
            delta = 0;
        } else {
            delta = static_cast<double>(upper_bound+1) / (n-1);
        }
        for (int i = 0; i < n-1; i++) {
            samples.push_back(floor(i*delta));
        }
        samples.push_back(upper_bound);
    }
};

struct _dna4_traits : seqan3::sequence_file_input_default_traits_dna {
    /**
     * @brief Syntax for reading the query file.
     */
    using sequence_alphabet = seqan3::dna4; // instead of dna5
 
    template <typename alph>
    using sequence_container = std::vector<alph>; // must be defined as a template!
};


#endif