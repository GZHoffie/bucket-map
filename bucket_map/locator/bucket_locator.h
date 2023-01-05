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

    // sampled kmers from the query sequences
    std::vector<std::pair<unsigned int, std::vector<unsigned int>>> query_sequences;

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
    std::unordered_map<unsigned int, unsigned int> fuzzy_counter;
    std::unordered_map<unsigned int, unsigned int> exact_counter;


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


    std::vector<std::pair<unsigned int, unsigned int>> 
    _find_offset(const std::unordered_multimap<unsigned int, int>& bucket_kmer_index, int sequence_id) {
        /**
         * @brief Return the most probable offset of the given sequence in the bucket.
         * @param bucket_kmer_index a pointer to the k-mer index for the target bucket.
         * @param sequence_id the index of sequence in the query file.
         * @returns a vector containing the possible offsets and the obtained number of votes.
         */
        // load the sequence and pick up the k-mers
        auto record = query_sequences[sequence_id];
        auto kmer_hash = std::get<1>(record);

        // sample k-mers
        sampler->sample_deterministically(std::get<0>(record)-1);

        // initialize result
        std::vector<std::pair<unsigned int, unsigned int>> res;

        // reset counter
        fuzzy_counter.clear();
        exact_counter.clear();

        // record the possible starting positions
        for (auto i : sampler->samples) {
            // find the positions of the k-mer
            auto range = bucket_kmer_index.equal_range(kmer_hash[i]);
            for (auto it = range.first; it != range.second; ++it) {
                auto position = it->second - i;
                exact_counter[position]++;
                for (int indel = -allowed_indel; indel <= allowed_indel; indel++) {
                    fuzzy_counter[position + indel]++;
                }
            }
        }

        // find potential good offsets that are higher than the threshold
        std::erase_if(fuzzy_counter, [&](const auto& item) {
            auto const& [key, value] = item;
            return (value < num_samples - allowed_mismatch) || (key < 0);
        });

        // find potential good offset
        // We choose the smallest offset that contains a specific number of k-mers
        if (!fuzzy_counter.empty()) {
            //TODO: return all the valid offsets within a bucket. Currently only the best offset is returned.
            auto max = std::max_element(fuzzy_counter.begin(), fuzzy_counter.end(), [&](const auto &x, const auto &y) {
                                            return (x.second < y.second) || (x.second == y.second && exact_counter[x.first] < exact_counter[y.first]);
                                        });
            if (max->second >= num_samples - allowed_mismatch && max->first >= 0) {
                res.push_back(std::pair(max->first, max->second));
            }
        }
        return res;
    }

    void _prepare_read_query(std::filesystem::path const & sequence_file) {
        /**
         * @brief Read the query file and store the sampled k-mers from the reads in `query_sequences`.
         * @note Should be run before locating reads in the query sequence file.
         */
        seqan3::sequence_file_input<_dna4_traits> fin{sequence_file};

        for (auto & rec : fin) {
            auto kmers = rec.sequence() | seqan3::views::kmer_hash(q_gram_shape);

            // sample k-mers from the read
            sampler->sample_deterministically(kmers.size()-1);
            auto sampled_kmers = sampler->samples | std::views::transform([&](unsigned int i) {
                return kmers[i];
            });
            std::vector<unsigned int> sampled_kmer_vec(sampled_kmers.begin(), sampled_kmers.end());

            // store the sampled k-mers in `query_sequences`.
            query_sequences.push_back(std::make_pair(kmers.size(), sampled_kmer_vec));
        }
        

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
        fuzzy_counter.reserve(sample_size * allowed_indel);
        exact_counter.reserve(sample_size);

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

    typedef std::tuple<unsigned int, unsigned int, unsigned int> locate_t;

    std::vector<std::vector<locate_t>> 
    _locate(std::filesystem::path const & sequence_file) {
        /**
         * * This function is just for benchmarking
         * @brief Locate the exact locations for the reads in the fastq file.
         */
        // initialize benchmark timers
        Timer index_timer;
        Timer query_timer;

        // map reads to buckets
        auto sequence_ids = _m->map(sequence_file);
        unsigned int bucket_index = 0;
        //TODO: can delete the mapper at this point

        // prepare for the offset mapping
        _prepare_read_query(sequence_file);

        // initialize result
        std::vector<std::vector<locate_t>> res;
        for (int i = 0; i < query_sequences.size(); i++) {
            std::vector<locate_t> bucket{};
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
                for (auto &[offset, vote] : _find_offset(bucket_kmer_index, id)) {
                    res[id].push_back(std::make_tuple(i, offset, vote));
                }
            }
            query_timer.tock();
        }

        // print benchmark information
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total time used for building k-mer index for each bucket: " 
                             << index_timer.elapsed_seconds() << " s.\n";
        float time = query_timer.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total time used for finding exact location of the sequences: " 
                             << time << " s (" << time * 1000 * 1000 / query_sequences.size() << " μs/seq).\n";

        return res;
    }

     void _check_ground_truth(const std::vector<std::vector<locate_t>>& query_results, 
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
        std::string cigar;
        
        // Test statistics
        int correct_map = 0;
        int total_bucket_numbers = 0;
        std::map<int, int> bucket_number_map;

        double total_offset_error = 0;
        int almost_correct_offset = 0;

        for (int i = 0; i < query_results.size(); i++) {
            is >> bucket >> exact_location >> cigar;
            std::vector<locate_t> buckets = query_results[i];

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