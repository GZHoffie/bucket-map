#ifndef BUCKET_MAP_BUCKET_LOCATOR_H
#define BUCKET_MAP_BUCKET_LOCATOR_H

#include "./locator.h"
#include <cstdlib>

#include <seqan3/search/fm_index/bi_fm_index.hpp>

class bucket_locator : public locator {
private:
    mapper* _m;

    // k-mer index for finding offset
    std::vector<std::unordered_map<unsigned int, int>> index;

    // information for buckets
    unsigned int bucket_length;
    unsigned int read_length;
    seqan3::shape q_gram_shape;

    // Parameters for verification
    float allowed_substitution_rate;
    float allowed_indel_rate;

    void _initialize_kmer_index(std::filesystem::path const & fasta_file_name) {
        /**
         * @brief Generate k-mer index as maps (kmer_hash -> its offset in the bucket)
         *        for each bucket in the reference genome, and store in the `index` variable
         * @param fasta_file_name the path to the fasta file containing reference genome.
         */
        auto operation = [&](std::vector<seqan3::dna4> seq, std::string id) {
            seqan3::debug_stream << "initializing k-mer index for bucket " << id << ".\n";
            unsigned int offset = 0;
            auto values = seq | seqan3::views::kmer_hash(q_gram_shape);
            std::unordered_map<unsigned int, int> bucket_index;
            for (auto hash : values) {
                // Only record the last appearance of the k-mer
                bucket_index[hash] = offset;
                offset++;
            }
            index.push_back(bucket_index);
        };
        Timer clock;
        clock.tick();
        iterate_through_buckets(fasta_file_name, bucket_length, read_length, operation);
        clock.tock();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for loading k-mer index: " 
                             << clock.elapsed_seconds() << " s." << '\n';
    }


    int _find_offset(const std::unordered_map<unsigned int, int>* bucket_kmer_index, int sequence_id) {
        /**
         * @brief Return the most probable offset of the given sequence in the bucket.
         *        return -1 if the sequence is considered not in the bucket.
         * @param bucket_kmer_index a pointer to the k-mer index for the target bucket.
         * @param sequence_id the index of sequence in the query file.
         */
        auto record = _m->records[sequence_id];
        auto kmers_positions = record | seqan3::views::kmer_hash(q_gram_shape) | std::views::transform([&](int i) {
            if (bucket_kmer_index->contains(i)) {
                return bucket_kmer_index->at(i);
            }
            return -1;
        });

        // max errors allowed
        int allowed_mismatch = ceil(allowed_substitution_rate * record.size());
        int allowed_indel = ceil(allowed_indel_rate * record.size());

        // record the possible starting positions
        // TODO: try with sampling k-mers instead of using all k-mers.
        std::unordered_map<int, unsigned int> count;
        std::vector<int> starting_pos(kmers_positions.begin(), kmers_positions.end());
        for (int i = 0; i < starting_pos.size(); i++) {
            int pos = starting_pos[i] - i;
            if (pos < 0) continue;
            for (int indel = -allowed_indel; indel <= allowed_indel; indel++) {
                count[pos + indel]++;
            }
        }

        // find potential good offset
        // We choose the smallest offset that contains a specific number of k-mers
        int offset = INT_MAX;
        for (auto & candidate : count) {
            if (candidate.second >= record.size() - allowed_mismatch && candidate.first < offset && candidate.first >= 0) {
                offset = candidate.first;
            }
        }
        return offset == INT_MAX? -1 : offset;
    }


public:
    bucket_locator(indexer* ind, mapper* map, unsigned int bucket_len, 
                   unsigned int read_len, seqan3::shape shape, float mismatch_rate, float indel_rate) : locator(ind) {
        _m = map;
        bucket_length = bucket_len;
        read_length = read_len;
        q_gram_shape = shape;

        allowed_substitution_rate = mismatch_rate;
        allowed_indel_rate = indel_rate;
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
        // initialize result
        std::vector<std::vector<std::pair<unsigned int, unsigned int>>> res;
        for (int i = 0; i < _m->records.size(); i++) {
            std::vector<std::pair<unsigned int, unsigned int>> bucket;
            res.push_back(bucket);
        }
        // map reads to buckets
        std::vector<std::vector<unsigned int>> sequence_ids = _m->map(sequence_file);
        unsigned int bucket_index = 0;
        for (auto bucket_sequence : sequence_ids) {
            // for each bucket, get the corresponding kmer index
            std::unordered_map<unsigned int, int>* bucket_kmer_index = &index[bucket_index];
            // go through all sequences mapped to this bucket and check the offset.
            for (auto & id : bucket_sequence) {
                int offset = _find_offset(bucket_kmer_index, id);
                if (offset > 0) {
                    // the sequence is mapped to an exact location in the bucket
                    res[id].push_back(std::make_pair(bucket_index, offset));
                }
            }
            bucket_index++;
        }
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

        int total_offset_error = 0;

        for (int i = 0; i < query_results.size(); i++) {
            is >> bucket >> exact_location;
            std::vector<std::pair<unsigned int, unsigned int>> buckets = query_results[i];

            for (auto & pair : buckets) {
                // check the bucket id
                if (bucket == std::get<0>(pair)) {
                    correct_map++;
                    int error = std::get<1>(pair) - exact_location;
                    if (error > 0) {
                        total_offset_error += error;
                    } else {
                        total_offset_error -= error;
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
                             << ((float) total_offset_error) / query_results.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Average number of buckets returned: " 
                             << ((float) total_bucket_numbers) / query_results.size() << ".\n";
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