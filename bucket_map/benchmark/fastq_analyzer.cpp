#include "../utils.h"

class fastq_analyzer {
private:
    unsigned int max_length;
    unsigned int num_reads;
    unsigned int num_bases;
    double min_average_base_quality_per_read;
    double max_possible_errors_per_read;
    double average_base_quality;
    double possible_errors;

    // quality threshold to be considered an error
    unsigned int error_threshold;


public:
    fastq_analyzer(unsigned int threshold = 20) {
        max_length = 0;
        num_reads = 0;
        num_bases = 0;
        min_average_base_quality_per_read = 100;
        average_base_quality = 0;

        error_threshold = threshold;
    }


    void reset() {
        max_length = 0;
        num_reads = 0;
        num_bases = 0;
        min_average_base_quality_per_read = 100;
        average_base_quality = 0;
    }

    /**
     * @brief Read and output the relevant information about the reads in the fasta file.
     * 
     * @param fastq_path path to the fastq sequence file.
     */
    void read(const std::filesystem::path& fastq_path) {
        reset();
        seqan3::debug_stream << "[INFO]\t\tReading " << fastq_path << ".\n";
        seqan3::sequence_file_input<_phred94_traits> fin{fastq_path};
 
        for (auto & rec : fin) {
            max_length = std::max(max_length, (unsigned int)rec.sequence().size());
            num_reads++;
            num_bases += rec.sequence().size();


            double base_quality_sum = 0;
            double num_possible_errors = 0;
            for (auto & q : rec.base_qualities()) {
                base_quality_sum += q.to_rank();
                if (q.to_rank() < error_threshold) {
                    num_possible_errors++;
                }
            }
            auto average_qual = base_quality_sum / rec.sequence().size();
            min_average_base_quality_per_read = std::min(min_average_base_quality_per_read, average_qual);
            max_possible_errors_per_read = std::max(max_possible_errors_per_read, num_possible_errors);
            average_base_quality += base_quality_sum;
            possible_errors += num_possible_errors;
        }
        average_base_quality /= num_bases;
        possible_errors /= num_reads;
        
        // display information about the fastq file

        seqan3::debug_stream << "[BENCHMARK]\tTotal number of reads in the sequence file: " << num_reads << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\tTotal number of bases in the sequence file: " << num_bases << " ("
                             << ((double) num_bases) / num_reads << " per read in average).\n";
        seqan3::debug_stream << "[BENCHMARK]\tMaximum read length in the sequence file: " << max_length << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\tMinimum average base quality in each read: " << min_average_base_quality_per_read << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\tAverage base quality in the entire sequence file: " << average_base_quality << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\tMaximum number of bases with quality < " << error_threshold << ": " 
                             << max_possible_errors_per_read << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\tAverage number of bases with quality per read < " << error_threshold << ": "
                             << possible_errors << " (" << possible_errors / (num_bases / num_reads) * 100 << "\% estimated error rate).\n";
    }

};


int main() {
    fastq_analyzer analyzer(25);
    analyzer.read("/mnt/d/genome/DRR035999.fastq");
}