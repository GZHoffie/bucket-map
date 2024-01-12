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
        bool is_random = false;
    } map_position_t;

    std::unordered_map<std::string, unsigned int> read_id_to_index;
    std::unordered_map<std::string, unsigned int> sequence_id_to_index;
    std::vector<bool> mapped_reads;
    std::vector<bool> correctly_mapped_reads;
    std::vector<bool> mapped_random_reads;
    std::vector<bool> is_random_read;
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
     * Optional, only required for dwgsim files
     */
    void read_fasta_file(std::filesystem::path fasta_path) {
        seqan3::sequence_file_input fin{fasta_path};
        unsigned int index = 0;
 
        for (auto & record : fin) {
            sequence_id_to_index[_remove_substring_after_slash_or_blank(record.id())] = index;
            index++;
        }
        seqan3::debug_stream << sequence_id_to_index << "\n";
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
     * @param is_dwgsim if the sequence file is a dwgsim simulated fastq file. If yes, then
     *                  record the ground truth
     */
    void read_sequence_file(std::filesystem::path sequence_file, bool is_dwgsim = false) {
        seqan3::sequence_file_input<_phred94_traits> fin{sequence_file};
        unsigned int index = 0;

        for (auto & rec : fin) {
            auto renamed_id = _remove_substring_after_slash_or_blank(_space_to_underscore(rec.id()));
            read_id_to_index.emplace(renamed_id, index);

            if (is_dwgsim) {
                std::size_t pos = 0, prev = 0;
                std::vector<std::string> substrings;
                // split the read id according to deliminators
                while ((pos = renamed_id.find_first_of("_:", prev)) != std::string::npos) {
                    if (pos > prev)
                        substrings.push_back(renamed_id.substr(prev, pos-prev));
                        prev = pos + 1;
                }   
                if (prev < renamed_id.length())
                    substrings.push_back(renamed_id.substr(prev, std::string::npos));

                //seqan3::debug_stream << substrings << "\n";
                // record the ground truth
                map_position_t gt;
                gt.reverse_complement = static_cast<bool>(stoi(substrings[3]));
                gt.sequence_id = sequence_id_to_index[substrings[0]];
                gt.offset = stoi(substrings[1]);
                gt.is_random = static_cast<bool>(stoi(substrings[5]));
                is_random_read.push_back(gt.is_random);

                //seqan3::debug_stream << gt << "\n";
                std::vector<map_position_t> positions{gt};
                answer.push_back(positions);
            }
            index++;
        }

        // initialize private variables.
        for (int i = 0; i < index; i++) {
            mapped_reads.push_back(false);
            correctly_mapped_reads.push_back(false);
            mapped_random_reads.push_back(false);
        }
    }

    void benchmark(std::filesystem::path sam_path) {
        seqan3::debug_stream << "[BENCHMARK]\t" << "============ Benchmarking sam file " << sam_path << " ============\n";
        // Benchmark statistics
        int mapped_locations = 0; // number of mapped locations
        std::fill(mapped_reads.begin(), mapped_reads.end(), false);
        std::fill(correctly_mapped_reads.begin(), correctly_mapped_reads.end(), false);
        std::fill(mapped_random_reads.begin(), mapped_random_reads.end(), false);
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

                // if this read is random, stop evaluating
                if (is_random_read[sequence_id]) {
                    mapped_random_reads[sequence_id] = true;
                    continue;
                }

                // get the mapping information
                bool reverse_comp = static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
                unsigned int ref_id = record.reference_id().value();
                unsigned int offset = record.reference_position().value();

                // check correctness against answer
                bool acceptable = false;
                for (map_position_t pos : answer[sequence_id]) {
                    if (reverse_comp == pos.reverse_complement &&
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
        unsigned int num_random_reads = std::count(is_random_read.begin(), is_random_read.end(), true);
        unsigned int num_mapped_random_reads = std::count(mapped_random_reads.begin(), mapped_random_reads.end(), true);
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of reads: " 
                             << read_id_to_index.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of random reads: " 
                             << num_random_reads << ".\n";

        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of mapped reads: " 
                             << num_mapped_reads << " (" << ((float) num_mapped_reads) / (read_id_to_index.size() - num_random_reads) * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of correctly mapped reads (sensitivity): " 
                             << num_correct_mapped_reads << " (" << ((float) num_correct_mapped_reads) / (read_id_to_index.size() - num_random_reads) * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of mapped random reads (false positives): " 
                             << num_mapped_random_reads << " (" << ((float) num_mapped_random_reads) / (num_random_reads) * 100 << "%).\n";
        
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of mapped locations returned: " 
                             << mapped_locations << " (" << ((float) mapped_locations) / num_mapped_reads << " per mapped read).\n";
        
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of acceptable mapped locations (precision): " 
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
    sam_analyzer analyzer(10);
    analyzer.read_fasta_file("/home/zhenhao/mapping_data/GCA_004358405.1_ASM435840v1_genomic.fna");
    analyzer.read_sequence_file("/home/zhenhao/mapping_data/EColi_sim.bwa.read1.fastq", true);
    //analyzer.read_best_alignment_file("/home/zhenhao/data/mapping/ecoli_simulated_golden.sam");
    //analyzer.read_ground_truth_file("/home/guzh/data/mapping/lr_simulated/sd_0001.maf");

    //analyzer.read_best_alignment_file("/home/zhenhao/bucket-map/bucket_map/benchmark/output/bowtie2_map.sam");
    //analyzer.read_best_alignment_file("/home/zhenhao/bucket-map/bucket_map/benchmark/output/subread_map.sam");
    //analyzer.read_ground_truth_file("/mnt/d/genome/test/EGU_1500_10K.position_ground_truth");
    analyzer.benchmark_directory("/home/zhenhao/bucket-map/bucket_map/benchmark/short_read/output");
    return 0;

}