#include <seqan3/io/sam_file/all.hpp>
#include "../utils.h"

#include <string>

class sam_analyzer {
private:
    // position of the read
    typedef struct {
        bool reverse_complement;
        unsigned int sequence_id;
        unsigned int offset;
    } map_position_t;

    std::unordered_map<std::string, unsigned int> read_id_to_index;
    std::vector<bool> mapped_reads;
    std::vector<bool> correctly_mapped_reads;
    std::vector<std::vector<map_position_t>> answer;

    unsigned int offset_error_tolerance;

    /**
     * @brief Replace the space in id with an underscore.
     * 
     * @param id the id of the read
     * @return std::string the id with the ' ' character replaced with '_'.
     */
    std::string _space_to_underscore(std::string id) {
        std::replace(id.begin(), id.end(), ' ', '_');
        return id;
    }


    /**
     * @brief The id of the read may contain a slash, followed by a number.
     *        Some mappers would delete that slash and number. Therefore, we also
     *        delete that when storing in `read_id_to_index`.
     * 
     * @param id the id of the read
     * @return std::string the id with the '/' character and everything after it deleted.
     */
    std::string _remove_substring_after_slash_or_blank(std::string id) {
        std::size_t slash = id.find('/');
        std::size_t blank = id.find(' ');
        if (std::min(blank, slash) != std::string::npos) {
            std::string res(id.begin(), id.begin() + std::min(blank, slash));
            return res;
        } else {
            return id;
        }
    }

public:
    sam_analyzer(unsigned int error_tolerance = 5) {
        //read_id_to_index.clear();
        offset_error_tolerance = error_tolerance;
    }

    /**
     * @brief Specify a good mapper, and store the mapping positions in 
     *        `answer`.
     * @note must be run after `read_sequence_file`.
     * @param sam_path path to the good alignment file.
     */
    void read_best_alignment_file(std::filesystem::path sam_path) {
        seqan3::sam_file_input fin{sam_path};
        unsigned int recorded_answers = 0;

        for (auto && record : fin) {

            //seqan3::debug_stream << ref << ", " << pos << " | ";
            try {
                auto renamed_id = _remove_substring_after_slash_or_blank(_space_to_underscore(record.id()));
                unsigned int sequence_id = read_id_to_index.at(renamed_id);

                if (static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped)) {
                    // read is unmapped
                    continue;
                }
                // otherwise, store the mapping position in answer
                map_position_t pos;
                pos.reverse_complement = static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
                pos.sequence_id = record.reference_id().value();
                pos.offset = record.reference_position().value();
                answer[sequence_id].push_back(pos);
                recorded_answers++;

            } catch (const std::out_of_range& oor) {
                //seqan3::debug_stream << "key not found\n";
            }
        }
        seqan3::debug_stream << "[INFO]\t\tMapped sequences in the best alignment file: " << recorded_answers << ".\n";
    }

    /**
     * @brief If there is a `.position_ground_truth` file, read that and set as the correct answer
     * @note must be run after `read_sequence_file`.
     * @param ground_truth_path path to the good alignment file.
     */
    void read_ground_truth_file(std::filesystem::path ground_truth_path) {
        std::ifstream is(ground_truth_path);

        if (ground_truth_path.extension() == ".ground_truth") {
            // for my own simulator
            unsigned int origin; 
            unsigned int position;
            unsigned int reverse_complement;
            std::string cigar;

            unsigned int index = 0;

            while (is >> origin >> position >> reverse_complement >> cigar) {
                map_position_t pos;
                pos.reverse_complement = static_cast<bool>(reverse_complement);
                pos.sequence_id = origin;
                pos.offset = position;

                answer[index].push_back(pos);
                index++;
            }
        } else if (ground_truth_path.extension() == ".maf") {
            // for pbsim3 simulation tool
            std::string _i; // string to be ignored.
            unsigned int offset;
            std::string read_name;
            char rev_comp;
            while (is >> _i >> _i >> _i >> offset >> _i >> _i >> _i >> _i 
                      >> _i >> read_name >> _i >> _i >> rev_comp >> _i >> _i) {
            
                try {
                    unsigned int sequence_id = read_id_to_index.at(read_name);

                    map_position_t pos;
                    pos.reverse_complement = false;
                    if (rev_comp == '-') {
                        pos.reverse_complement = true;
                    }
                    unsigned int start = read_name.find('S') + 1;
                    unsigned int end = read_name.find('_');
                    pos.sequence_id = std::stoi(read_name.substr(start, end - start)) - 1;
                    pos.offset = offset;
                    answer[sequence_id].push_back(pos);
                } catch (const std::out_of_range& oor) {
                    //seqan3::debug_stream << "key not found\n";
                }
            }
        }
        
        
    }

    

    /**
     * @brief Read the fastq file, and record the sequence ids for identification 
     *        in the SAM file.
     * 
     * @param sequence_file path to the fastq file
     */
    void read_sequence_file(std::filesystem::path sequence_file) {
        seqan3::sequence_file_input<_phred94_traits> fin{sequence_file};
        unsigned int index = 0;

        for (auto & rec : fin) {
            auto renamed_id = _remove_substring_after_slash_or_blank(_space_to_underscore(rec.id()));
            read_id_to_index.emplace(renamed_id, index);
            index++;
        }

        // initialize private variables.
        for (int i = 0; i < index; i++) {
            mapped_reads.push_back(false);
            correctly_mapped_reads.push_back(false);
            std::vector<map_position_t> positions;
            answer.push_back(positions);
        }
    }

    void benchmark(std::filesystem::path sam_path) {
        seqan3::debug_stream << "[BENCHMARK]\t" << "============ Benchmarking sam file " << sam_path << " ============\n";
        // Benchmark statistics
        int mapped_locations = 0; // number of mapped locations
        std::fill(mapped_reads.begin(), mapped_reads.end(), false);
        std::fill(correctly_mapped_reads.begin(), correctly_mapped_reads.end(), false);
        unsigned int acceptable_maps = 0;

        seqan3::sam_file_input fin{sam_path};

        for (auto && record : fin) {

            //seqan3::debug_stream << ref << ", " << pos << " | ";
            try {
                auto renamed_id = _remove_substring_after_slash_or_blank(_space_to_underscore(record.id()));
                unsigned int sequence_id = read_id_to_index.at(renamed_id);

                if (static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped)) {
                    // read is unmapped
                    continue;
                }
                // otherwise, the read is mapped
                mapped_reads[sequence_id] = true;
                mapped_locations++;

                // get the mapping information
                bool reverse_comp = static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
                unsigned int ref_id = record.reference_id().value();
                unsigned int offset = record.reference_position().value();

                // check correctness against answer
                bool acceptable = false;
                for (map_position_t pos : answer[sequence_id]) {
                    if (//reverse_comp == pos.reverse_complement &&
                        ref_id == pos.sequence_id &&
                        std::abs((int)offset - (int)pos.offset) <= offset_error_tolerance) {
                        correctly_mapped_reads[sequence_id] = true;
                        acceptable = true;
                    }
                }
                if (acceptable) acceptable_maps++;

            } catch (const std::out_of_range& oor) {
                //seqan3::debug_stream << "key not found\n";
            } catch (const seqan3::format_error& fe) {

            }
        }

        // print out benchmark results
        unsigned int num_mapped_reads = std::count(mapped_reads.begin(), mapped_reads.end(), true);
        unsigned int num_correct_mapped_reads = std::count(correctly_mapped_reads.begin(), correctly_mapped_reads.end(), true);
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of mapped reads: " 
                             << num_mapped_reads << " (" << ((float) num_mapped_reads) / read_id_to_index.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of correctly mapped reads: " 
                             << num_correct_mapped_reads << " (" << ((float) num_correct_mapped_reads) / read_id_to_index.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of sequences: " 
                             << read_id_to_index.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of mapped locations returned: " 
                             << mapped_locations << " (" << ((float) mapped_locations) / num_mapped_reads << " per mapped read).\n";
        
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of acceptable mapped locations: " 
                             << acceptable_maps << " (precision: " << ((float) acceptable_maps) / mapped_locations * 100 << "%).\n";
    }

    void benchmark_directory(std::filesystem::path sam_directory) {
        for (const auto& entry : std::filesystem::directory_iterator(sam_directory)) {
            if (entry.path().extension() == ".sam") {
                benchmark(entry.path());
            }
        }    
    }
};


int main()
{
    sam_analyzer analyzer(2000);
    analyzer.read_sequence_file("/home/guzh/data/mapping/lr_simulated/sd_0001.fastq");
    //analyzer.read_best_alignment_file("/home/zhenhao/data/mapping/ecoli_simulated_golden.sam");
    analyzer.read_ground_truth_file("/home/guzh/data/mapping/lr_simulated/sd_0001.maf");
    //analyzer.read_best_alignment_file("/home/zhenhao/bucket-map/bucket_map/benchmark/output/bowtie2_map.sam");
    //analyzer.read_best_alignment_file("/home/zhenhao/bucket-map/bucket_map/benchmark/output/subread_map.sam");
    //analyzer.read_ground_truth_file("/mnt/d/genome/test/EGU_1500_10K.position_ground_truth");
    analyzer.benchmark_directory("/home/guzh/bucket-map/bucket_map/benchmark/long_read/output");
    return 0;

}