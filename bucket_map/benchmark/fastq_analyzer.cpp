#include "../utils.h"

class fastq_analyzer {
private:
    unsigned int max_length;
    unsigned int num_reads;
    unsigned int num_bases;
    double min_average_base_quality_per_read;
    double average_base_quality;


public:
    fastq_analyzer() {
        max_length = 0;
        num_reads = 0;
        num_bases = 0;
        min_average_base_quality_per_read = 0;
        average_base_quality = 0;
    }


    void reset() {
        max_length = 0;
        num_reads = 0;
        num_bases = 0;
        min_average_base_quality_per_read = 0;
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
        seqan3::sequence_file_input fin{fastq_path};
 
        for (auto & rec : fin) {
            max_length = std::max(max_length, rec.sequence().size());
            num_reads++;
            num_bases += rec.sequence().size();

            min_average_base_quality_per_read = std::min(min_average_base_quality_per_read)

            double base_quality_sum = 0;
            for (auto & q : rec.base_qualities()) {
                base_quality_sum += q.to_rank();
            }
            auto average_qual = base_quality_sum / rec.sequence().size();
            min_average_base_quality_per_read = std::min(min_average_base_quality_per_read, average_qual);
            average_base_quality += base_quality_sum;
        }
        average_base_quality /= num_bases;
        
        // display information about the fastq file

        seqan3::debug_stream << "[BENCHMARK]\tTotal number of reads in the sequence file: " << num_reads << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\tTotal number of bases in the sequence file: " << num_bases << " ("
                             << ((double) num_bases) / num_reads << " per read in average.\n";
        seqan3::debug_stream << "[BENCHMARK]\tMaximum read length in the sequence file: " << max_length << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\tMinimum average base quality in each read: " << min_average_base_quality_per_read << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\tAverage base quality in the entire sequence file: " << average_base_quality << ".\n";
    }

}


int main() {
    fastq_analyzer analyzer;
    analyzer.read("/mnt/d/genome/TS1.81.90.001.fq");
}