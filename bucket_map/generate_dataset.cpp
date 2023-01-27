#include "tools/short_read_simulator.h"


int main() {
    std::filesystem::path data_path = "/mnt/d/genome";
    std::filesystem::path genome_file = data_path / "GRCh38_latest_genomic.fasta";

    int bucket_length = 65536;
    int read_length = 300;

    short_read_simulator sim(bucket_length, read_length, 0.002, 0.00025, 0.00025);
    sim.read(genome_file);
    sim.generate_fastq_file(data_path / "test", "GRCH38_300_5M", 5000000);
}