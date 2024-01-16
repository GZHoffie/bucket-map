#ifndef BUCKET_MAP_BUCKET_LOCATOR_H
#define BUCKET_MAP_BUCKET_LOCATOR_H

#include "./locator.h"
#include "../mapper/quality_filter.h"
#include "../utils.h"
#include <cstdlib>

 
#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/all.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>


template <typename kmer_hash_t = unsigned int, typename index_t = uint16_t>
class query_sequences_storage {
private:
    // a map from segment id to its index used to store `kmers`, `indices` and `segment_length`.
    std::map<segment_info_t, unsigned int> segment_to_index;

    // kmer samples and their locations
    std::vector<std::vector<kmer_hash_t>> kmers;
    std::vector<std::vector<index_t>> indices;

    // length of each segment and each read.
    std::vector<unsigned int> segment_length;
    std::vector<unsigned int> read_length;

    unsigned int s;
    unsigned int current_index;

public:

    query_sequences_storage(unsigned int sample_size) {
        s = sample_size;
        current_index = 0;
    }

    void reset() {
        indices.clear();
        indices.shrink_to_fit();
        kmers.clear();
        kmers.shrink_to_fit();
        read_length.clear();
        read_length.shrink_to_fit();
        current_index = 0;
    }

    index_t get_segment_length(segment_info_t segment_id) {
        unsigned int segment_index = segment_to_index[segment_id];
        return segment_length[segment_index];
    }

    index_t get_read_length(segment_info_t segment_id) {
        return read_length[std::get<0>(segment_id)];
    }

    /**
     * @brief Get the recorded k-mer samples of the sequence.
     * 
     * @param segment_id the index of sequence in the storage type.
     */
    std::vector<kmer_hash_t>& get_samples(const segment_info_t& segment_id) {
        unsigned int segment_index = segment_to_index[segment_id];
        return kmers[segment_index];
    }

    /**
     * @brief Return the starting position of the sampled k-mers.
     * 
     * @param segment_id the index of sequence in the storage type.
     */
    std::vector<index_t>& get_indices(const segment_info_t& segment_id) {
        unsigned int segment_index = segment_to_index[segment_id];
        return indices[segment_index];
    }

    /**
     * @brief Store the segment information, including the segment id, length of the segment, and the k-mer samples.
     * 
    */
    void push_back(const segment_info_t& segment_id, unsigned int segment_len,
                   const std::vector<index_t>& kmer_indices, const std::vector<kmer_hash_t>& kmer_samples) {
        // check the size of the two vectors
        assert(kmer_indices.size() == s && kmer_samples.size() == s);
    
        indices.push_back(kmer_indices);
        kmers.push_back(kmer_samples);
        segment_length.push_back(segment_len);
        segment_to_index[segment_id] = current_index;
        
        current_index++;
    }


    void record_read_length(unsigned int read_len) {
        read_length.push_back(read_len);
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
    unsigned int min_base_quality;
    std::vector<bool> high_quality_kmers;
    uint8_t k; //size of k-mer

    // Parameters for verification
    int allowed_mismatch;
    int allowed_indel;
    float allowed_indel_rate;

    // number of samples drawn
    int num_samples;

    // sampler of k-mers
    Sampler* sampler;
    Sampler* segment_sampler;

    // counter to calculate the possible starting positions
    std::map<int, unsigned int> vote_counter;
    //std::unordered_map<unsigned int, unsigned int> exact_counter;

    // the path to the genome file
    std::filesystem::path genome_file_name;

    // type used to store mapped locations
    typedef std::tuple<unsigned int, // bucket id
                       int,          // segment offset within the bucket
                       unsigned int, // segment offset in the read
                       unsigned int, // number of votes get
                       bool          // whether the original read (true) or its reverse complement (false) is mapped
                       > locate_t;


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
        index.reserve(BM_BUCKET_LEN);
        // insert the k-mers.
        int offset = 0;
        for (auto hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{k})) {
            // Only record the last appearance of the k-mer
            index.emplace(hash, offset);
            offset++;
        }
    }

    /**
     * @brief Find all high quality k-mers. (Same as in the mapper)
     * 
     * @param quality the vector of quality alphabets.
     * fills in `high_quality_kmers`: A vector with the i-th element representing whether
     *                                the i-th k-mer is high-quality.
     * TODO: support spaced k-mer.
     */
    void _high_quality_kmers(const std::vector<seqan3::phred94>& quality) {
        unsigned int consecutive_high_quality_base = 0, index = 0;
        high_quality_kmers.clear();
        for (auto & qual: quality) {
            if (qual.to_rank() >= min_base_quality) {
                consecutive_high_quality_base++;
                if (consecutive_high_quality_base >= k) {
                    high_quality_kmers.push_back(true);
                    continue;
                }
            } else {
                consecutive_high_quality_base = 0;
            }
            index++;
            if (index >= k) {
                high_quality_kmers.push_back(false);
            }
        }
        assert(high_quality_kmers.size() == quality.size() - k + 1);
    }


    std::pair<int, unsigned int> _find_offset(const std::unordered_multimap<unsigned int, int>& bucket_kmer_index, segment_info_t segment,
                                              query_sequences_storage<unsigned int>* sequence_records, bool reverse_complement) {
        /**
         * @brief Return the most probable offset of the given sequence in the bucket.
         * @param bucket_kmer_index a pointer to the k-mer index for the target bucket.
         * @param sequence_id the index of sequence in the query file.
         * @param reverse_complement whether we try to find the map of reverse complement of the read.
         * @returns A possible offset for the sequence. -1 if not found.
         */
        // load the sequence and pick up the k-mers
        auto record = sequence_records->get_samples(segment);
        auto index = sequence_records->get_indices(segment);
        unsigned int length = sequence_records->get_segment_length(segment);

        // initialize result
        std::vector<std::pair<unsigned int, unsigned int>> res;

        // reset counter
        vote_counter.clear();

        // record the possible starting positions
        std::unordered_set<unsigned int> votes;
        std::vector<unsigned int> mapped_positions;

        for (unsigned int i = 0; i < num_samples; i++) {
            // start from the first sample in original string, or start from last sample for reverse string.
            unsigned int sample_index = i;
            if (reverse_complement) sample_index = num_samples - 1 - i;

            // find the positions of the k-mer (adjust if in the reverse string)
            unsigned int current_kmer = record[sample_index], current_index = index[sample_index];
            if (reverse_complement) {
                current_kmer = hash_reverse_complement(current_kmer, k);
                current_index = length - k - current_index;  // where the k-mer starts
            }

            // Cast the vote if found in the k-mer index
            auto range = bucket_kmer_index.equal_range(current_kmer);
            if (vote_counter.empty()) {
                // if there's no proposal, make a proposal
                for (auto it = range.first; it != range.second; ++it) {
                    auto position = it->second - current_index;  // where the k-mer propose as the starting point of the read
                    vote_counter[position]++;
                }
            } else {
                for (auto it = range.first; it != range.second; ++it) {
                    // if there is a close enough proposal, vote for it
                    bool voted = false;
                    auto position = it->second - current_index;

                    auto lower_bound = vote_counter.lower_bound(position - allowed_indel);
                    auto upper_bound = vote_counter.upper_bound(position + allowed_indel);

                    for (auto it = lower_bound; it != upper_bound; ++it) {
                        vote_counter[it->first]++;
                        voted = true;
                    }
                    
                    if (!voted) {
                        // otherwise, propose the new position
                        vote_counter[position]++;
                    }
                }

            }
        }

        //seqan3::debug_stream << segment << vote_counter << "\n";

        // find potential good offset
        // We choose the smallest offset that contains a specific number of k-mers
        if (!vote_counter.empty()) {
            auto max = std::max_element(vote_counter.begin(), vote_counter.end(), [&](const auto &x, const auto &y) {
                                            return (x.second < y.second) || (x.second == y.second && x.first > y.first);
                                        });
            if (max->second >= num_samples - allowed_mismatch && max->first >= 0) {
                // This position is a good fit...
                return std::make_pair(max->first, max->second);
            }
        }
        return std::make_pair(-1, 0);
    }

    void _prepare_read_query(std::filesystem::path const & sequence_file) {
        /**
         * @brief Read the query file and store the sampled k-mers from the reads in `records`.
         * @note Should be run before locating reads in the query sequence file.
         */
        seqan3::sequence_file_input<_phred94_traits> fastq_file{sequence_file};
        unsigned int read_index = 0;

        for (auto & rec : fastq_file) {

            std::vector<unsigned int> starting_positions{0};
            // if long read, break it down into several shorter reads.
            if (rec.sequence().size() > 2 * read_length) {
                segment_sampler->sample_deterministically(rec.sequence().size() - read_length - 1);
                starting_positions = segment_sampler->samples;
            }

            for (auto i : starting_positions) {
                // copy the segment of length read_length
                //TODO: optimize this by not copying the segment, but simply indicate the range.
                auto begin = i;
                auto end = std::min(i + read_length, (unsigned int)rec.sequence().size());
                seqan3::dna4_vector segment_sequence(rec.sequence().begin() + begin, rec.sequence().begin() + end);
                std::vector<seqan3::phred94> segment_quality(rec.base_qualities().begin() + begin, rec.base_qualities().begin() + end);

                // break the read down into k-mers
                auto kmers = segment_sequence | seqan3::views::kmer_hash(seqan3::ungapped{k});
                auto kmer_qualities = segment_quality | seqan3::views::kmer_quality(seqan3::ungapped{k});

                // copy the segment of length read_length
                int num_kmers = kmers.size();

                // filter out low-quality k-mers
                auto good_kmers = std::ranges::iota_view{0, num_kmers} | std::views::filter([&](unsigned int j) {
                                      return kmer_qualities[j] >= min_base_quality;
                                  });
                // sample kmers
                std::vector<uint16_t> good_indices(good_kmers.begin(), good_kmers.end());
                if (good_indices.empty()) { // consider all k-mers if none of the kmers are high-quality.
                    for (uint16_t j = 0; j < num_kmers; j++) good_indices.push_back(j);
                }
                sampler->sample_deterministically(good_indices.size() - 1);
                auto sampled_indices = sampler->samples | std::views::transform([&](unsigned int j) { return good_indices[j]; });
                auto selected_kmers = sampled_indices | std::views::transform([&](unsigned int j) { return kmers[j]; });

                // convert to vectors and store them
                std::vector<uint16_t> sampled_indices_vec(sampled_indices.begin(), sampled_indices.end());
                std::vector<unsigned int> sampled_kmers_vec(selected_kmers.begin(), selected_kmers.end());
                segment_info_t segment{read_index, (int)i};
                records->push_back(segment, segment_sequence.size(), sampled_indices_vec, sampled_kmers_vec);
            }
            // record the length of the entire read as well
            records->record_read_length(rec.sequence().size());
            read_index++;
        }
    }


    std::vector<locate_t> _filter_best_locations(std::vector<locate_t>& mapped_locations, unsigned int read_length) {
        /**
         * @brief find the best mapped locations. Let them vote for the best final location.
         */
        //seqan3::debug_stream << mapped_locations << "\n";
        // All locations cast their votes
        std::map<std::tuple<unsigned int, int, bool>, unsigned int> loc_votes;
        for (auto &[bucket_id, bucket_offset, segment_offset, votes, rev_comp] : mapped_locations) {
            if (loc_votes.empty()) {
                // propose as a new location
                //seqan3::debug_stream << "Proposed new!" << {bucket_id, bucket_offset, rev_comp} << "\n";
                loc_votes[{bucket_id, bucket_offset, rev_comp}] = votes;
            } else {
                bool found_close_loc = false;
                // if a close location is proposed, vote for that position.
                int lower_bound = bucket_offset - read_length * allowed_indel_rate;
                int upper_bound = bucket_offset + read_length * allowed_indel_rate;
                //seqan3::debug_stream << "Will merge if " << {lower_bound, upper_bound} << "\n";
                for (auto it = loc_votes.cbegin(); it != loc_votes.cend(); ++it) {
                    int proposed_offset = std::get<1>(it->first);

                    if (bucket_id == std::get<0>(it->first) && proposed_offset <= upper_bound 
                                                            && proposed_offset >= lower_bound
                                                            && std::get<2>(it->first) == rev_comp) {
                        //seqan3::debug_stream << "Merged." << "\n";
                        loc_votes[it->first] += votes;
                        found_close_loc = true;
                    }
                }
                if (!found_close_loc) {
                    loc_votes[{bucket_id, bucket_offset, rev_comp}] = votes;
                }
            }
        }

        //seqan3::debug_stream << loc_votes << "\n";

        // find the locations that gets the most votes.
        std::vector<locate_t> res;
        unsigned int max_votes = 0;
        for(auto it = loc_votes.cbegin(); it != loc_votes.cend(); ++it ) {
            if (it->second > max_votes) {
                res.clear();
                max_votes = it->second;
            }
            if (it->second == max_votes) {
                res.push_back(std::make_tuple(std::get<0>(it->first),
                                              std::get<1>(it->first),
                                              0,
                                              it->second,
                                              std::get<2>(it->first)));
            }
        }
        //seqan3::debug_stream << "Res: " << res << "\n";
        return res;
    }


public:
    bucket_locator(indexer* ind, mapper* map, unsigned int bucket_len, 
                   unsigned int read_len, uint8_t seed_len, 
                   float mismatch_rate, float indel_rate, 
                   unsigned int sample_size, unsigned int quality_threshold,
                   unsigned int num_segment_samples = 5) : locator(ind) {
        _m = map;
        bucket_length = bucket_len;
        read_length = read_len;
        k = seed_len;

        allowed_mismatch = ceil(mismatch_rate * sample_size);
        allowed_indel = ceil(indel_rate * read_len);
        allowed_indel_rate = indel_rate;

        // random number generator
        num_samples = sample_size;
        srand(time(NULL));

        // initialize sampler
        sampler = new Sampler(num_samples);
        segment_sampler = new Sampler(num_segment_samples);

        min_base_quality = quality_threshold * k;
    }

    ~bucket_locator() {
        delete sampler;
        delete segment_sampler;
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
        seqan3::sequence_file_input<_phred94_traits> query_file_in{sequence_file};
 
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

#ifndef BM_ALIGN
            //seqan3::debug_stream << read_id << "\n";
            locate_res[read_id] = _filter_best_locations(locate_res[read_id], record.sequence().size());
#endif

            // output to sam file
            for (auto & loc : locate_res[read_id]) {
                // get the components of the locate_t
                const auto & [bucket_id, offset, segment_offset, votes, is_original] = loc;

                // get the part of text that the read is mapped to
                auto start = bucket_seq[bucket_id].begin() + offset;
                size_t width = std::min(query.size() + 1 + (std::size_t)(allowed_indel_rate * record.sequence().size()), bucket_seq[bucket_id].size() - offset);
                
                // define flag value
                seqan3::sam_flag flag{seqan3::sam_flag::none};
                if (!is_original) {
                    // read is mapped to the reverse strand
                    flag |= seqan3::sam_flag::on_reverse_strand;
                }


#ifdef BM_ALIGN
                // pick up the text
                std::vector<seqan3::dna4> text(start, start + width);
                if (!is_original) {
                    auto text_rev_comp = text | std::views::reverse | seqan3::views::complement;
                    std::vector<seqan3::dna4> text_copy(text_rev_comp.begin(), text_rev_comp.end());
                    text = text_copy;
                }
                // do pairwise string alignment
                for (auto && alignment : seqan3::align_pairwise(std::tie(text, query), align_config)) {
                    size_t map_qual = 60u + alignment.score();
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
                size_t map_qual = std::min((unsigned int)60, 6 * votes);
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
        records = new query_sequences_storage(num_samples);

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
        bucket_kmer_index.reserve(BM_BUCKET_LEN);

        for (int i = 0; i < sequence_ids_orig.size(); i++) {
            auto & bucket_sequence_orig = sequence_ids_orig[i];
            auto & bucket_sequence_rev_comp = sequence_ids_rev_comp[i];
            //seqan3::debug_stream << i << bucket_sequence_orig << bucket_sequence_rev_comp << "\n";

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

            // check original read
            for (auto & id : bucket_sequence_orig) {
                auto offset_res = _find_offset(bucket_kmer_index, id, records, false);
                auto & [offset, vote] = offset_res;
                if (offset > 0) {
                    res[std::get<0>(id)].push_back(std::make_tuple(i, 
                                                                   offset - std::get<1>(id), 
                                                                   std::get<1>(id), 
                                                                   vote, true));
                }
            }
            
            // check reverse complement of the reads read
            for (auto & id_rev : bucket_sequence_rev_comp | std::views::reverse) {
                auto offset_res = _find_offset(bucket_kmer_index, id_rev, records, true);
                auto & [offset, vote] = offset_res;
                if (offset > 0) {
                    int segment_offset_ = records->get_read_length(id_rev) - std::get<1>(id_rev) - records->get_segment_length(id_rev);
                    res[std::get<0>(id_rev)].push_back(std::make_tuple(i, 
                                                                       offset - segment_offset_, 
                                                                       std::get<1>(id_rev), 
                                                                       vote, false));
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