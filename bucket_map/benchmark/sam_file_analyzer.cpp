#include <seqan3/io/sam_file/all.hpp>
#include "../utils.h"

class sam_analyzer {
private:
    std::map<std::string, unsigned int> read_id_to_index;

public:
    sam_analyzer() {
        read_id_to_index.clear();
    }

    /**
     * @brief Read the fastq file, and record the sequence ids for identification 
     *        in the SAM file.
     * 
     * @param sequence_file path to the fastq file
     */
    void read_sequence_file(std::filesystem::path sequence_file) {
        seqan3::sequence_file_input fin{sequence_file};
        unsigned int index = 0;

        for (auto & rec : fin) {
            std::string id(rec.id());
            read_id_to_index.insert({id, index});
            index++;
        }
    }

    void benchmark(std::filesystem::path sam_path) {
        // Benchmark statistics
        int mapped_locations = 0; // number of mapped locations
        std::vector<bool> mapped_reads(read_id_to_index.size(), false); // whether the read is mapped to a location

        seqan3::sam_file_input fin{sam_path};

        for (auto && record : fin) {
            unsigned int sequence_id = read_id_to_index[record.id()];

            //seqan3::debug_stream << ref << ", " << pos << " | ";
            try {
                std::string id(record.id());
                unsigned int sequence_id = read_id_to_index[id];

                if (static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped)) {
                    // read is unmapped
                    continue;
                }
                // otherwise, the read is mapped
                mapped_reads[sequence_id] = true;
                mapped_locations++;
            } catch (const std::out_of_range& oor) {
                //seqan3::debug_stream << "\n";
            }
        }

        // print out benchmark results
        unsigned int num_mapped_reads = std::count(mapped_reads.begin(), mapped_reads.end(), true);
        //seqan3::debug_stream << "[BENCHMARK]\t" << "============ Benchmarking sam file " << sam_path << " ============\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of mapped reads: " 
                             << num_mapped_reads << " (" << ((float) num_mapped_reads) / read_id_to_index.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of sequences: " 
                             << read_id_to_index.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of mapped locations returned: " 
                             << mapped_locations << " (" << ((float) mapped_locations) / num_mapped_reads << " per mapped read).\n";
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
    sam_analyzer analyzer;
    analyzer.read_sequence_file("/mnt/d/genome/TS1.81.90.001.fq");

    analyzer.benchmark("/home/zhenhao/bucket-map/bucket_map/benchmark/output/bowtie2_map.sam");
    return 0;

}