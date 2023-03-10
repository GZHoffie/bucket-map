#ifndef BUCKET_MAP_BUCKET_LOCATOR_H
#define BUCKET_MAP_BUCKET_LOCATOR_H

#include "./locator.h"
#include "../mapper/quality_filter.h"
#include <cstdlib>

 
#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/all.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>


template <typename kmer_hash_t = unsigned int, typename index_t = uint16_t>
class query_sequences_storage {
private:
    typedef struct kmer_samples {
        index_t read_length;
        std::vector<kmer_hash_t> kmers;
        std::vector<index_t> indices;
    } kmer_samples_t;

    std::vector<kmer_samples_t> kmer_samples;
    unsigned int k;

public:

    query_sequences_storage(unsigned int kmer_size) {
        k = kmer_size;
    }

    void reset() {
        kmer_samples.clear();
        kmer_samples.shrink_to_fit();
    }

    kmer_hash_t get_num_kmers(unsigned int sequence_id) {
        return kmer_samples[sequence_id].num_kmers;
    }

    /**
     * @brief Get the recorded k-mer samples of the sequence.
     * 
     * @param sequence_id the index of sequence in the storage type.
     * @param reverse_complement whether we are checking the original read (false) or the 
     *                           reverse complement of the read (true)
     * @return std::vector<kmer_hash_t>& return the vector storing the sampled k-mers.
     */
    std::vector<kmer_hash_t> get_samples(unsigned int sequence_id, bool reverse_complement) {
        if (reverse_complement) {
            // find the reverse complement of the sampled k-mers
            auto kmers_rev_comp = kmer_samples[sequence_id].kmers | std::views::transform([&](unsigned int hash) {
                return hash_reverse_complement(hash, k);
            });
            std::vector<kmer_hash_t> res(kmers_rev_comp.begin(), kmers_rev_comp.end());
            return res;
        } else {
            return kmer_samples[sequence_id].kmers;
        }
    }

    /**
     * @brief Return the starting position of the sampled k-mers.
     * 
     * @param sequence_id the index of sequence in the storage type.
     * @param reverse_complement whether we are checking the original read (false) or the 
     *                           reverse complement of the read (true) 
     * @return std::vector<index_t>& A vector storing the indices of the sampled k-mers.
     */
    std::vector<index_t> get_indices(unsigned int sequence_id, bool reverse_complement) {
        if (reverse_complement) {
            auto kmers_indices = kmer_samples[sequence_id].indices | std::views::transform([&](unsigned int index) {
                return kmer_samples[sequence_id].read_length - k - index;
            });
            std::vector<index_t> res(kmers_indices.begin(), kmers_indices.end());
            return res;
        } else {
            return kmer_samples[sequence_id].indices;
        }
    }



    void push_back(index_t read_length, std::vector<index_t> kmer_indices, std::vector<kmer_hash_t> kmers) {
        kmer_samples_t samples(read_length, kmers, kmer_indices);
        kmer_samples.push_back(samples);
    }
};


class bucket_locator : public locator {
private:
    mapper* _m;

    // sequence of each bucket
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> bucket_seq;

    // sampled kmers from the query sequences
    query_sequences_storage<unsigned int, uint16_t>* records; // original sequence

    // information for buckets
    unsigned int bucket_length;
    unsigned int read_length;
    unsigned int average_quality;
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
    //std::unordered_map<unsigned int, unsigned int> exact_counter;

    // the path to the genome file
    std::filesystem::path genome_file_name;


    void _initialize_kmer_index(std::filesystem::path const & fasta_file_name) {
        /**
         * @brief Insert the sequence into `bucket_seq` variable.
         * @param fasta_file_name the path to the fasta file containing reference genome.
         */
        auto operation = [&](const seqan3::bitpacked_sequence<seqan3::dna4>& seq, const std::string& id) {
            bucket_seq.push_back(seq);
        };
        iterate_through_buckets(fasta_file_name, bucket_length, read_length, operation);
    }

    void _create_kmer_index(std::unordered_multimap<unsigned int, int>& index, const seqan3::bitpacked_sequence<seqan3::dna4>& seq) {
        /**
         * @brief Create the k-mer index given the sequence in bucket.
         * @param index the k-mer index storing the positions of each k-mer.
         * @param seq the sequence of the bucket.
         */
        index.clear();
        // insert the k-mers.
        int offset = 0;
        for (auto hash : seq | seqan3::views::kmer_hash(q_gram_shape)) {
            // Only record the last appearance of the k-mer
            index.emplace(hash, offset);
            offset++;
        }
    }


    std::pair<int, unsigned int> _find_offset(const std::unordered_multimap<unsigned int, int>& bucket_kmer_index, int sequence_id,
                                              query_sequences_storage<unsigned int>* sequence_records, bool reverse_complement) {
        /**
         * @brief Return the most probable offset of the given sequence in the bucket.
         * @param bucket_kmer_index a pointer to the k-mer index for the target bucket.
         * @param sequence_id the index of sequence in the query file.
         * @param reverse_complement whether we try to find the map of reverse complement of the read.
         * @returns A possible offset for the sequence. -1 if not found.
         */
        // load the sequence and pick up the k-mers
        auto record = sequence_records->get_samples(sequence_id, reverse_complement);

        // sample k-mers
        auto samples = sequence_records->get_indices(sequence_id, reverse_complement);
  
        // initialize result
        std::vector<std::pair<unsigned int, unsigned int>> res;

        // reset counter
        fuzzy_counter.clear();
        //exact_counter.clear();

        //seqan3::debug_stream << "Sequence " << sequence_id << ":\n";
        // record the possible starting positions
        std::unordered_set<unsigned int> votes;
        for (unsigned int i = 0; i < samples.size(); i++) {
            // find the positions of the k-mer
            //seqan3::debug_stream << "K-mer " << *(record + i) << ": ";
            auto range = bucket_kmer_index.equal_range(record[i]);
            for (auto it = range.first; it != range.second; ++it) {
                votes.clear();
                auto position = it->second - samples[i];
                for (int indel = -allowed_indel; indel <= allowed_indel; indel++) {
                    votes.emplace(position + indel);
                }
                for (auto & vote : votes) {
                    fuzzy_counter[vote]++;
                }
            }
            //seqan3::debug_stream << "\n";
        }
        
        
        //seqan3::debug_stream << fuzzy_counter << "\n";
        //seqan3::debug_stream << exact_counter << "\n";

        // find potential good offset
        // We choose the smallest offset that contains a specific number of k-mers
        if (!fuzzy_counter.empty()) {
            //TODO: return all the valid offsets within a bucket. Currently only the best offset is returned.
            auto max = std::max_element(fuzzy_counter.begin(), fuzzy_counter.end(), [&](const auto &x, const auto &y) {
                                            //return (x.second < y.second) || (x.second == y.second && exact_counter[x.first] < exact_counter[y.first]);
                                            return (x.second < y.second) || (x.second == y.second && x.first > y.first);
                                        });
            if (max->second >= num_samples - allowed_mismatch && max->first >= 0) {
                //seqan3::debug_stream << max->first << " " << max->second << "\n";
                return std::make_pair(max->first, max->second);
            } else {
                //seqan3::debug_stream << "Sequence: " << sequence_id << ", "<<  fuzzy_counter << "\n";
            }
        }
        return std::make_pair(-1, 0);
    }

    void _prepare_read_query(std::filesystem::path const & sequence_file) {
        /**
         * @brief Read the query file and store the sampled k-mers from the reads in `records`.
         * @note Should be run before locating reads in the query sequence file.
         */
        seqan3::sequence_file_input<_dna4_traits> fastq_file{sequence_file};

        for (auto & rec : fastq_file) {
            auto kmers = rec.sequence() | seqan3::views::kmer_hash(q_gram_shape);
            auto kmer_qualities = rec.base_qualities() | seqan3::views::kmer_quality(q_gram_shape);

            int num_kmers = kmers.size();

            // filter out low-quality k-mers
            auto good_kmers = std::ranges::iota_view{0, num_kmers} | std::views::filter([&](unsigned int i) {
                                  return kmer_qualities[i] >= average_quality;
                              });
            // sample kmers
            std::vector<uint16_t> good_indices(good_kmers.begin(), good_kmers.end());
            sampler->sample_deterministically(good_indices.size() - 1);
            auto sampled_indices = sampler->samples | std::views::transform([&](unsigned int i) { return good_indices[i]; });
            auto selected_kmers = sampled_indices | std::views::transform([&](unsigned int i) { return kmers[i]; });

            // convert to vectors and store them
            std::vector<uint16_t> sampled_indices_vec(sampled_indices.begin(), sampled_indices.end());
            std::vector<unsigned int> sampled_kmers_vec(selected_kmers.begin(), selected_kmers.end());
            records->push_back(rec.sequence().size(), sampled_indices_vec, sampled_kmers_vec);
        }
        

    }


public:
    bucket_locator(indexer* ind, mapper* map, unsigned int bucket_len, 
                   unsigned int read_len, seqan3::shape shape, 
                   float mismatch_rate, float indel_rate, 
                   unsigned int sample_size, unsigned int quality_threshold) : locator(ind) {
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
        //exact_counter.reserve(sample_size);

        // initialize sampler
        sampler = new Sampler(num_samples);

        average_quality = quality_threshold * q_gram_shape.count();
    }

    ~bucket_locator() {
        delete sampler;
        delete records;
    }

    unsigned int initialize(std::filesystem::path const & fasta_file_name, 
                            std::filesystem::path const & index_directory,
                            std::string const & indicator) {
        genome_file_name = fasta_file_name;

        // do the indexing
        auto res = locator::initialize(fasta_file_name, index_directory, indicator);

        // load q-gram index to the mapper
        _m->load(index_directory, indicator);

        // create index files
        return res;
    }

    void locate(std::filesystem::path const & sequence_file, 
                std::filesystem::path const & index_file,
                std::filesystem::path const & sam_file,
                unsigned int quality_threshold = 30) {
        /**
         * @brief Find the exact locations of the reads and output the mapping results
         *        to a sam file.
         * @param sequence_file the fastq file containing all the query reads.
         * @param index_file the `.bucket_id` file created by a `bucket_indexer` object,
         *                   which contains the name of each bucket.
         * @param sam_file the path to the output sam file.
         * @param quality_threshold the minimum quality score for an alignment to be valid.
         */
        // find the mapped locations of all the reads.
        auto locate_res = _locate(sequence_file);
        records->reset();

        // store the bucket information
        std::ifstream bucket_info(index_file);
        std::vector<std::string> bucket_name; // the name for the bucket
        std::vector<unsigned int> bucket_offsets; // the starting index of bucket in the chromosome.

        // store the sequence information for sam file headers.
        std::vector<std::string> ref_ids;
        std::vector<size_t> ref_lengths;

        std::string name;
        std::string last_bucket_name;

        unsigned int bucket_index = 0;
        for (int i = 0; i < bucket_seq.size(); i++) {
            std::getline(bucket_info, name);
            name = name.substr(0, name.find(' '));
            if (name != last_bucket_name) {
                if (bucket_index != 0) {
                    ref_ids.push_back(last_bucket_name);
                    ref_lengths.push_back(bucket_index * bucket_length); //FIXME: an upper bound on bucket length
                }
                last_bucket_name = name;
                bucket_index = 0;
            }
            bucket_name.push_back(name);
            bucket_offsets.push_back(bucket_index * bucket_length);
            bucket_index++;
        }
        if (bucket_index != 0) {
            ref_ids.push_back(last_bucket_name);
            ref_lengths.push_back(bucket_index * bucket_length); //FIXME: an upper bound on bucket length
        }

        // output to the sam file
        seqan3::sequence_file_input query_file_in{sequence_file};
 
        seqan3::sam_file_output sam_out{sam_file, ref_ids, ref_lengths,
                                        seqan3::fields<seqan3::field::seq,
                                                       seqan3::field::id,
                                                       seqan3::field::ref_id,
                                                       seqan3::field::ref_offset,
#ifdef BM_ALIGN
                                                       seqan3::field::cigar,
#endif
                                                       seqan3::field::qual,
                                                       seqan3::field::flag,
                                                       seqan3::field::mapq>{}};
 
#ifdef BM_ALIGN
        seqan3::configuration const align_config =
            seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                             seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                             seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                             seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
            | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_alignment{}
            | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{};
#endif

        unsigned int read_id = 0;
        unsigned int mapped_locations = 0;

        Timer align_timer;
        align_timer.tick();
        for (auto && record : query_file_in) {
            auto & query = record.sequence();
            for (auto & loc : locate_res[read_id]) {
                // get the components of the locate_t
                const auto & [bucket_id, offset, votes, is_original] = loc;

                // get the part of text that the read is mapped to
                auto start = bucket_seq[bucket_id].begin() + offset;
                size_t width = std::min(query.size() + 1 + allowed_indel, bucket_seq[bucket_id].size() - offset);

                std::vector<seqan3::dna4> text(start, start + width);

                // define flag value
                seqan3::sam_flag flag{seqan3::sam_flag::none};
                if (!is_original) {
                    // read is mapped to the reverse strand
                    flag |= seqan3::sam_flag::on_reverse_strand;
                }


#ifdef BM_ALIGN
                // do pairwise string alignment
                for (auto && alignment : seqan3::align_pairwise(std::tie(text, query), align_config)) {
                    size_t map_qual = 60u + alignment.score(); // TODO: consider the number of votes in the mapping quality.
                    if (map_qual < quality_threshold) {
                        continue;
                    }

                    auto cigar = seqan3::cigar_from_alignment(alignment.alignment());
                    size_t ref_offset = alignment.sequence1_begin_position() + bucket_offsets[bucket_id] + offset;
                    
                    
                    // output to the sam file
                    sam_out.emplace_back(query,
                                         record.id(),
                                         bucket_name[bucket_id],
                                         ref_offset,
                                         cigar,
                                         record.base_qualities(),
                                         flag,
                                         map_qual);
                    mapped_locations++;
                }
#else
                size_t map_qual = 6 * votes;
                size_t ref_offset = bucket_offsets[bucket_id] + offset;
                sam_out.emplace_back(query,
                                     record.id(),
                                     bucket_name[bucket_id],
                                     ref_offset,
                                     record.base_qualities(),
                                     flag,
                                     map_qual);
                mapped_locations++;
#endif
            }
            read_id++;
        }
        align_timer.tock();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total mapped locations: " 
                             << mapped_locations << " (" << (float) mapped_locations / read_id << " per sequence).\n";
        auto elapsed_seconds = align_timer.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total time used for alignment verification and output: " 
                             << elapsed_seconds << " s (" << (float) elapsed_seconds / mapped_locations * 1000 * 1000 << " μs per pairwise alignment).\n";
    }

    typedef std::tuple<unsigned int, // bucket id
                       unsigned int, // offset within the bucket
                       unsigned int, // number of votes get
                       bool          // whether the original read (true) or its reverse complement (false) is mapped
                       > locate_t;

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
        auto [sequence_ids_orig, sequence_ids_rev_comp] = _m->map(sequence_file);

        // reset mapper to release memory
        _m->reset();

        // load all index files
        _initialize_kmer_index(genome_file_name);

        // initialze query sequence storage
        records = new query_sequences_storage(q_gram_shape.count());

        unsigned int bucket_index = 0;

        // prepare for the offset mapping
        _prepare_read_query(sequence_file);

        // initialize result
        std::vector<std::vector<locate_t>> res;
        for (int i = 0; i < _m->num_records; i++) {
            std::vector<locate_t> bucket{};
            res.push_back(bucket);
        }

        // create k-mer index
        std::unordered_multimap<unsigned int, int> bucket_kmer_index;

        for (int i = 0; i < sequence_ids_orig.size(); i++) {
            auto & bucket_sequence_orig = sequence_ids_orig[i];
            auto & bucket_sequence_rev_comp = sequence_ids_rev_comp[i];

            // skip if no read is mapped to this bucket
            if (bucket_sequence_orig.empty() && bucket_sequence_rev_comp.empty()) {
                continue;
            }
            
            // for each bucket, create the corresponding kmer index
            index_timer.tick();
            //seqan3::debug_stream << "===========>>> Sequence " << i << "\n"; 
            _create_kmer_index(bucket_kmer_index, bucket_seq[i]);
            index_timer.tock();
            
            // go through all sequences mapped to this bucket and check the offset.
            query_timer.tick();
            /*
            for (auto & id : bucket_sequence_orig) {
                // check original read
                auto offset_res = _find_offset(bucket_kmer_index, id, records, false);
                auto & [offset, vote] = offset_res;
                if (offset > 0) {
                    res[id].push_back(std::make_tuple(i, offset, vote, true));
                }
            }
            */

            // check reverse complement of the reads read
            for (auto & id_rev : bucket_sequence_rev_comp) {
                auto offset_res = _find_offset(bucket_kmer_index, id_rev, records, true);
                auto & [offset, vote] = offset_res;
                if (offset > 0) {
                    res[id_rev].push_back(std::make_tuple(i, offset, vote, false));
                }
            }
            query_timer.tock();
        }

        // print benchmark information
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total time used for building k-mer index for each bucket: " 
                             << index_timer.elapsed_seconds() << " s.\n";
        float time = query_timer.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total time used for finding exact location of the sequences: " 
                             << time << " s (" << time * 1000 * 1000 / _m->num_records << " μs/seq).\n";

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
        int unique_correct_offset = 0;

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
                        if (buckets.size() == 1) {
                            unique_correct_offset++;
                        }
                    } else {
                        //seqan3::debug_stream << std::get<1>(pair) << " answer: " << exact_location << "\n";
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
        seqan3::debug_stream << "[BENCHMARK]\t" << "Correct bucket predictions among uniquely mapped reads: " 
                             << unique_correct_offset << " (" << ((float) unique_correct_offset) / bucket_number_map[1] * 100 << "%).\n";
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