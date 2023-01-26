#include <seqan3/io/sam_file/all.hpp>
#include "../utils.h"

class sam_analyzer {
private:
    std::vector<std::pair<int, int>> ground_truth;
    int allowed_error;

public:
    sam_analyzer(int error) {
        allowed_error = error;
    }

    void read(std::filesystem::path ground_truth_path) {
        std::ifstream is(ground_truth_path);
        int origin; 
        int position;
        std::string cigar;

        while (is >> origin >> position >> cigar) {
            ground_truth.push_back(std::make_pair(origin, position));
        }
    }

    void benchmark(std::filesystem::path sam_path) {
        // Benchmark statistics
        int mapped_locations = 0; // number of mapped locations
        std::vector<bool> mapped_reads(ground_truth.size(), false); // whether the read is mapped to a location
        std::vector<bool> correct_mapped_references(ground_truth.size(), false); // number of correctly mapped sequences in reference
        std::vector<bool> correct_mapped_positions(ground_truth.size(), false); // number of correctly mapped locations
        double total_absolute_error = 0; // difference between correctly mapped sequences and true offsets

        seqan3::sam_file_input fin{sam_path};

        for (auto && record : fin) {
            mapped_locations++;

            // get gound_truth
            unsigned int sequence_id = std::stoi(record.id());
            auto & [ref, pos] = ground_truth[sequence_id];

            //seqan3::debug_stream << ref << ", " << pos << " | ";
            try {
                // get predicted values
                int predict_ref = record.reference_id().value();
                int predict_pos = record.reference_position().value();

                mapped_reads[sequence_id] = true;

                //seqan3::debug_stream << predict_ref << ", " << predict_pos << "\n";

                // compare the two
                if (ref == predict_ref) {
                    correct_mapped_references[sequence_id] = true;
                    int error = std::abs(predict_pos - pos);
                    
                    if (error <= allowed_error && !correct_mapped_positions[sequence_id]) {
                        total_absolute_error += error;
                        correct_mapped_positions[sequence_id] = true;
                    }
                }
            } catch (std::bad_optional_access e) {
                //seqan3::debug_stream << "\n";
            }
        }

        // print out benchmark results
        unsigned int num_mapped_reads = std::count(mapped_reads.begin(), mapped_reads.end(), true);
        unsigned int correct_buckets = std::count(correct_mapped_references.begin(), correct_mapped_references.end(), true);
        unsigned int correct_locations = std::count(correct_mapped_positions.begin(), correct_mapped_positions.end(), true);
        seqan3::debug_stream << "[BENCHMARK]\t" << "============ Benchmarking sam file " << sam_path << " ============\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of mapped reads: " 
                             << num_mapped_reads << " (" << ((float) num_mapped_reads) / ground_truth.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of sequences: " 
                             << ground_truth.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Correct reference sequence predictions: " 
                             << correct_buckets << " (" << ((float) correct_buckets) / ground_truth.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of mapped locations returned: " 
                             << mapped_locations << " (" << ((float) mapped_locations) / ground_truth.size() << "/read).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "MAE of offset: " 
                             << total_absolute_error / correct_locations << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "(Almost) correct offset calculations: " 
                             << correct_locations << " (" << ((float) correct_locations) / ground_truth.size() * 100 << "%).\n";
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
    sam_analyzer analyzer(5);
    analyzer.read("/mnt/d/genome/test/sim_150_1M.position_ground_truth");

    analyzer.benchmark_directory("/home/zhenhao/bucket-map/bucket_map/benchmark/output");
    return 0;

}