#include <fstream>
#include <span>
 
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/all.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
 
#include <cereal/archives/binary.hpp>


class fm_index_mapper {
private:
    unsigned int num_errors_toleration;
    seqan3::bi_fm_index<seqan3::dna5, seqan3::text_layout::collection> index;


public:
    fm_index_mapper(int errors, std::filesystem::path const & index_path) {
        num_errors_toleration = errors;
        {
            std::ifstream is{index_path, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(index);
        }
    }

    ~fm_index_mapper(){}
    
    std::pair<std::vector<std::vector<size_t>>, std::vector<std::vector<size_t>>>
    map_reads(std::filesystem::path const & query_path)
    {
        // lists storing the mapping results
        std::vector<std::vector<size_t>> positions;
        std::vector<std::vector<size_t>> ids;

        seqan3::sequence_file_input query_file_in{query_path};
    
        seqan3::configuration const search_config =
            seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{num_errors_toleration}}
            | seqan3::search_cfg::hit_all_best{};
    
        seqan3::configuration const align_config =
            seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                             seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                             seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                             seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
            | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_alignment{}
            | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{};
    
        for (auto && record : query_file_in)
        {
            auto & query = record.sequence();
            std::vector<size_t> query_positions;
            std::vector<size_t> query_ids;
            for (auto && result : search(query, index, search_config))
            {
                size_t start = result.reference_begin_position() ? result.reference_begin_position() - 1 : 0;
                query_positions.push_back(start);
                query_ids.push_back(result.reference_id());
            }
        }
        return std::make_pair(positions, ids);
    }

};
 
