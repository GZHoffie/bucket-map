#ifndef BUCKET_MAP_BUCKET_LOCATOR_H
#define BUCKET_MAP_BUCKET_LOCATOR_H

#include "./locator.h"
#include <cstdlib>

#include <seqan3/search/fm_index/bi_fm_index.hpp>

class bucket_locator : public locator {
private:
    mapper* _m;

    // sequence of each bucket
    std::vector<std::vector<seqan3::dna4>> bucket_seq;

    // information for buckets
    unsigned int bucket_length;
    unsigned int read_length;
    seqan3::shape q_gram_shape;

    // Parameters for verification
    int allowed_mismatch;
    int allowed_indel;

    // number of samples drawn
    int num_samples;

    // sampler of k-mers
    Sampler* sampler;

    // counter to calculate the possible starting positions
    std::unordered_map<unsigned int, unsigned int> counter;


    void _initialize_kmer_index(std::filesystem::path const & fasta_file_name) {
        /**
         * @brief Insert the sequence into `bucket_seq` variable.
         * @param fasta_file_name the path to the fasta file containing reference genome.
         */
        auto operation = [&](std::vector<seqan3::dna4> seq, std::string id) {
            bucket_seq.push_back(seq);
        };
        iterate_through_buckets(fasta_file_name, bucket_length, read_length, operation);
    }

    void _create_kmer_index(std::unordered_multimap<unsigned int, int>& index, const std::vector<seqan3::dna4>& seq) {
        /**
         * @brief Create the k-mer index given the sequence in bucket.
         * @param index the k-mer index storing the positions of each k-mer.
         * @param seq the sequence of the bucket.
         */
        index.clear();
        // insert the k-mers.
        auto values = seq | seqan3::views::kmer_hash(q_gram_shape);
        int offset = 0;
        for (auto hash : values) {
            // Only record the last appearance of the k-mer
            index.emplace(hash, offset);
            offset++;
        }
    }


    int _find_offset(const std::unordered_multimap<unsigned int, int>& bucket_kmer_index, int sequence_id) {
        /**
         * @brief Return the most probable offset of the given sequence in the bucket.
         *        return -1 if the sequence is considered not in the bucket.
         * @param bucket_kmer_index a pointer to the k-mer index for the target bucket.
         * @param sequence_id the index of sequence in the query file.
         */
        // load the sequence and pick up the k-mers
        auto record = _m->records[sequence_id];
        auto kmer_hash = record | seqan3::views::kmer_hash(q_gram_shape);

        // sample k-mers
        sampler->sample_deterministically(kmer_hash.size()-1);

        // reset counter
        counter.clear();

        // record the possible starting positions
        for (auto i : sampler->samples) {
            // find the positions of the k-mer
            auto range = bucket_kmer_index.equal_range(kmer_hash[i]);
            for (auto it = range.first; it != range.second; ++it) {
                auto position = it->second - i;
                for (int indel = -allowed_indel; indel <= allowed_indel; indel++) {
                    counter[position + indel]++;
                }
            }
        }

        // find potential good offset
        // We choose the smallest offset that contains a specific number of k-mers
        if (!counter.empty()) {
            auto res = std::max_element(counter.begin(), counter.end(), [](const auto &x, const auto &y) {
                                            return (x.second < y.second) || (x.second == y.second && x.first > y.first);
                                        });
            if (res->second >= num_samples - allowed_mismatch && res->first >= 0) {
                return res->first;
            }
        }
        
        return -1;
    }


public:
    bucket_locator(indexer* ind, mapper* map, unsigned int bucket_len, 
                   unsigned int read_len, seqan3::shape shape, 
                   float mismatch_rate, float indel_rate, 
                   unsigned int sample_size) : locator(ind) {
        _m = map;
        bucket_length = bucket_len;
        read_length = read_len;
        q_gram_shape = shape;

        allowed_mismatch = ceil(mismatch_rate * sample_size);
        allowed_indel = ceil(indel_rate * sample_size);

        // random number generator
        num_samples = sample_size;
        srand(time(NULL));

        // initialize counter
        counter.reserve(sample_size * allowed_indel);

        // initialize sampler
        sampler = new Sampler(num_samples);
    }

    ~bucket_locator() {
        delete sampler;
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


    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> _locate(std::filesystem::path const & sequence_file) {
        /**
         * * This function is just for benchmarking
         * @brief Locate the exact locations for the reads in the fastq file.
         */
        // initialize benchmark timers
        Timer index_timer;
        Timer query_timer;

        // map reads to buckets
        std::vector<std::vector<unsigned int>> sequence_ids = _m->map(sequence_file);
        unsigned int bucket_index = 0;

        // initialize result
        std::vector<std::vector<std::pair<unsigned int, unsigned int>>> res;
        for (int i = 0; i < _m->records.size(); i++) {
            std::vector<std::pair<unsigned int, unsigned int>> bucket{};
            res.push_back(bucket);
        }

        // create k-mer index
        std::unordered_multimap<unsigned int, int> bucket_kmer_index;

        for (int i = 0; i < sequence_ids.size(); i++) {
            auto & bucket_sequence = sequence_ids[i];

            // skip if no read is mapped to this bucket
            if (bucket_sequence.empty()) {
                continue;
            }
            
            // for each bucket, create the corresponding kmer index
            index_timer.tick();
            _create_kmer_index(bucket_kmer_index, bucket_seq[i]);
            index_timer.tock();
            
            // go through all sequences mapped to this bucket and check the offset.
            query_timer.tick();
            for (auto & id : bucket_sequence) {
                int offset = _find_offset(bucket_kmer_index, id);
                // the sequence is mapped to an exact location in the bucket
                if (offset > 0) {
                    res[id].push_back(std::make_pair(i, offset));
                }
            }
            query_timer.tock();
        }

        // print benchmark information
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total time used for building k-mer index for each bucket: " 
                             << index_timer.elapsed_seconds() << " s.\n";
        float time = query_timer.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total time used for finding exact location of the sequences: " 
                             << time << " s (" << time * 1000 * 1000 / _m->records.size() << " μs/seq).\n";

        return res;
    }

     void _check_ground_truth(const std::vector<std::vector<std::pair<unsigned int, unsigned int>>>& query_results, 
                              std::filesystem::path ground_truth_file) {
        /**
         * * This function is just for benchmarking.
         * @brief Check the performance of query results against the ground truth.
         *        Output necessary information.
         * @note The ground truth file must be the one generated together with sequence
         *       fastq file.
         */
        std::ifstream is(ground_truth_file);
        int bucket, exact_location;

        // Test statistics
        int correct_map = 0;
        int total_bucket_numbers = 0;
        std::map<int, int> bucket_number_map;

        double total_offset_error = 0;
        int almost_correct_offset = 0;

        for (int i = 0; i < query_results.size(); i++) {
            is >> bucket >> exact_location;
            std::vector<std::pair<unsigned int, unsigned int>> buckets = query_results[i];

            for (auto & pair : buckets) {
                // check the bucket id
                if (bucket == std::get<0>(pair)) {
                    correct_map++;
                    auto error = std::abs((double) std::get<1>(pair) - exact_location);
                    total_offset_error += error;
                    if (error <= allowed_indel) {
                        almost_correct_offset++;
                    }
                }
            }
            total_bucket_numbers += buckets.size();
            ++bucket_number_map[buckets.size()];
        }

        // output information
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of sequences: " 
                             << query_results.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Correct bucket predictions: " 
                             << correct_map << " (" << ((float) correct_map) / query_results.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "MAE of offset: " 
                             << total_offset_error / correct_map << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "(Almost) correct offset calculations: " 
                             << almost_correct_offset << " (" << ((float) almost_correct_offset) / query_results.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Average number of buckets returned: " 
                             << ((float) total_bucket_numbers) / query_results.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of sequences that have no candidate bucket: " 
                             << bucket_number_map[0] << " (" << ((float) bucket_number_map[0]) / query_results.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of uniquely mapped sequences: " 
                             << bucket_number_map[1] << " (" << ((float) bucket_number_map[1]) / query_results.size() * 100 << "%).\n";
        int small_bucket_numbers = 0;
        for (int i = 1; i <= 5; i++) {
            small_bucket_numbers += bucket_number_map[i];
        }
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of sequences mapped to <= 5 buckets: " 
                             << small_bucket_numbers << " (" << ((float) small_bucket_numbers) / query_results.size() * 100 << "%).\n";
        for (int i = 6; i <= 10; i++) {
            small_bucket_numbers += bucket_number_map[i];
        }
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of sequences mapped to <= 10 buckets: " 
                             << small_bucket_numbers << " (" << ((float) small_bucket_numbers) / query_results.size() * 100 << "%).\n";
    }
};

#endif