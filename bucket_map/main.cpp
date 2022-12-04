#include "mapper/q_gram_map.h"
#include "simulation/short_read_simulator.h"
#include "./utils.h"

#include <map>
#include <chrono>

int main() {
    std::filesystem::path genome_file = "/mnt/c/data/Egu.v3.genome_f.fasta";

    int bucket_length = 50000;
    int read_length = 150;

    q_gram_mapper<34570> map(bucket_length, read_length, 0b1110100101001101_shape, 20, 10);
    //map.read(genome_file);
    //map.store("/mnt/c/data/index");
    map.load("/mnt/c/data/index");

    //short_read_simulator sim(bucket_length, read_length, 0.002, 0.00025, 0.00025);
    //sim.read(genome_file);
    //sim.generate_fastq_file("/mnt/d/data/test", "sim", 1000000);

    auto res = map._query_file("/mnt/c/data/test/sim.fastq");
    map._check_ground_truth(res, "/mnt/c/data/test/sim.ground_truth");
    
    

}