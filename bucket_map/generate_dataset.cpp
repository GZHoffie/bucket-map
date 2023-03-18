#include "tools/short_read_simulator.h"


int main() {
    std::filesystem::path data_path = "/mnt/d/genome";
    std::filesystem::path genome_file = data_path / "Egu.v3.genome_f.fasta";

    int bucket_length = 65536;
    int read_length = 1500;

    short_read_simulator sim(bucket_length, read_length, 0.025, 0.025, 0.025);
    sim.read(genome_file);
    sim.generate_fastq_file(data_path / "test", "EGU_1500_10K", 10000);
}