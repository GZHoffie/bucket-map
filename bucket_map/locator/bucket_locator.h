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

    // number of buckets we process each time
    unsigned int batch_size;


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

    void _create_kmer_index(std::vector<std::unordered_multimap<unsigned short, unsigned int>>& index, unsigned int batch_id) {
        /**
         * @brief Create the k-mer index given the sequence in bucket.
         * @param index the k-mer index storing the positions of each k-mer, initially empty.
         * @param batch_id the batch we are processing
         */
        for (unsigned int batch = batch_id * batch_size; batch < (batch_id + 1) * batch_size && batch < bucket_seq.size(); batch++) {
            std::unordered_multimap<unsigned short, unsigned int> kmer_index;
            auto kmer_hash = bucket_seq[batch] | seqan3::views::kmer_hash(q_gram_shape);
            unsigned int offset = 0;
            for (auto hash : kmer_hash) {
                kmer_index.emplace(hash, offset);
                offset++;
            }
            index.push_back(kmer_index);
        }
    }


    int _find_offset(const std::vector<std::unordered_multimap<unsigned short, unsigned int>>& bucket_kmer_index, 
                     unsigned short bucket_id, const std::vector<seqan3::dna4>& sequence) {
        /**
         * @brief Return the most probable offset of the given sequence in the bucket.
         *        return -1 if the sequence is considered not in the bucket.
         * @param bucket_kmer_index a pointer to the k-mer index for the target bucket.
         * @param sequence the sequence to be found in the query file.
         */
        // load the sequence and pick up the k-mers
        auto & index = bucket_kmer_index[bucket_id];
        auto kmer_hash = sequence | seqan3::views::kmer_hash(q_gram_shape);

        // sample k-mers
        sampler->sample_deterministically(kmer_hash.size()-1);

        // reset counter
        counter.clear();

        // record the possible starting positions
        for (auto i : sampler->samples) {
            // find the positions of the k-mer
            auto range = index.equal_range(kmer_hash[i]);
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
                   unsigned int sample_size, unsigned int buckets_per_batch = 1000) : locator(ind) {
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

        batch_size = buckets_per_batch;
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


    std::vector<std::vector<std::pair<unsigned short, unsigned int>>> _locate(std::filesystem::path const & sequence_file) {
        /**
         * * This function is just for benchmarking
         * @brief Locate the exact locations for the reads in the fastq file.
         */
        // initialize benchmark timers
        Timer index_timer;
        Timer query_timer;

        // map reads to buckets
        auto map_buckets = _m->map(sequence_file);
        unsigned int bucket_index = 0;

        // initialize result
        std::vector<std::vector<std::pair<unsigned short, unsigned int>>> res;
        for (int i = 0; i < _m->num_reads; i++) {
            std::vector<std::pair<unsigned short, unsigned int>> bucket{};
            res.push_back(bucket);
        }


        for (unsigned int batch_id = 0; batch_id < map_buckets.size(); batch_id++) {

            // get the mapped sequences for each batch
            auto & batch_sequence = map_buckets[batch_id];
            // skip if no read is mapped to this batch of buckets
            if (batch_sequence.empty()) continue;
            // create k-mer index
            std::vector<std::unordered_multimap<unsigned short, unsigned int>> batch_kmer_index;
            
            // for each bucket, create the corresponding kmer index
            index_timer.tick();
            _create_kmer_index(batch_kmer_index, batch_id);
            index_timer.tock();
            
            // go through all sequences mapped to this bucket and check the offset.
            query_timer.tick();
            // read the query file
            seqan3::sequence_file_input<_dna4_traits> fin{sequence_file};
            unsigned int read_index = 0;
            for (auto & rec : fin) {
                if (batch_sequence[read_index].empty()) {
                    read_index++;
                    continue;
                }
                // if the read is mapped to at least one bucket
                for (auto bucket_id : batch_sequence[read_index]) {
                    int offset = _find_offset(batch_kmer_index, bucket_id, rec.sequence());
                    if (offset > 0) {
                        res[read_index].push_back(std::make_pair(batch_id * batch_size + bucket_id, offset));
                    }
                }
                read_index++;
            }
            query_timer.tock();
        }

        // print benchmark information
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total time used for building k-mer index for each bucket: " 
                             << index_timer.elapsed_seconds() << " s.\n";
        float time = query_timer.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total time used for finding exact location of the sequences: " 
                             << time << " s (" << time * 1000 * 1000 / _m->num_reads << " μs/seq).\n";

        return res;
    }

     void _check_ground_truth(const std::vector<std::vector<std::pair<unsigned short, unsigned int>>>& query_results, 
                              std::filesystem::path ground_truth_file) {
        /**
         * * This function is just for benchmarking.
         * @brief Check the performance of query results against the ground truth.
         *        Output necessary information.
         * @note The ground truth file must be the one generated together with sequence
         *       fastq file.
         */
        std::ifstream is(ground_truth_file);
        unsigned short bucket; 
        unsigned int exact_location;

        // Test statistics
        int correct_map = 0;
        int total_bucket_numbers = 0;
        std::map<int, int> bucket_number_map;

        double total_offset_error = 0;
        int almost_correct_offset = 0;

        for (int i = 0; i < query_results.size(); i++) {
            is >> bucket >> exact_location;
            auto buckets = query_results[i];

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