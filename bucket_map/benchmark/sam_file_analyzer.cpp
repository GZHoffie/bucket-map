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

        while (is >> origin >> position) {
            ground_truth.push_back(std::make_pair(origin, position));
        }
    }

    void benchmark(std::filesystem::path sam_path) {
        // Benchmark statistics
        int mapped_locations = 0; // number of mapped locations
        int correct_mapped_references = 0; // number of correctly mapped sequences in reference
        int correct_mapped_positions = 0; // number of correctly mapped locations
        double total_absolute_error = 0; // difference between correctly mapped sequences and true offsets

        seqan3::sam_file_input fin{sam_path};

        for (auto && record : fin) {
            mapped_locations++;

            // get gound_truth
            int sequence_id = std::stoi(record.id());
            int ref = std::get<0>(ground_truth[sequence_id]);
            int pos = std::get<1>(ground_truth[sequence_id]);

            //seqan3::debug_stream << ref << ", " << pos << " | ";
            try {
                // get predicted values
                int predict_ref = record.reference_id().value();
                int predict_pos = record.reference_position().value();

                //seqan3::debug_stream << predict_ref << ", " << predict_pos << "\n";

                // compare the two
                if (ref == predict_ref) {
                    correct_mapped_references++;
                    int error = std::abs(predict_pos - pos);
                    total_absolute_error += error;
                    if (error <= allowed_error) {
                        correct_mapped_positions++;
                    }
                }
            } catch (std::bad_optional_access e) {
                //seqan3::debug_stream << "\n";
            }
        }

        // print out benchmark results
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of sequences: " 
                             << ground_truth.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Correct reference sequence predictions: " 
                             << correct_mapped_references << " (" << ((float) correct_mapped_references) / ground_truth.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of mapped locations returned: " 
                             << mapped_locations << " (" << ((float) mapped_locations) / ground_truth.size() << "/read).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "MAE of offset: " 
                             << total_absolute_error / correct_mapped_references << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "(Almost) correct offset calculations: " 
                             << correct_mapped_positions << " (" << ((float) correct_mapped_positions) / ground_truth.size() * 100 << "%).\n";
    }
};


int main()
{
    sam_analyzer analyzer(2);
    analyzer.read("/mnt/d/genome/test/sim.position_ground_truth");

    analyzer.benchmark("/home/zhenhao/bucket-map/bucket_map/benchmark/output/bowtie2_map.sam");
    return 0;

}