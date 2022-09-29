#include "mapper/q_gram_map.h"
#include "simulation/short_read_simulator.h"

#include <map>

int main() {
    std::filesystem::path genome_file = "/mnt/c/data/Egu.v3.genome_f.fasta";

    int bucket_length = 50000;
    int read_length = 150;

    //q_gram_mapper<34570> map(bucket_length, read_length, 0b1110100101001101_shape, 20, 10);
    //map.read(genome_file);
    //map.store("/home/zhenhao/mcomp-dissertation/build/sequence_sample");

    short_read_simulator sim(bucket_length, read_length, 0.002, 0.00025, 0.00025);
    sim.read(genome_file);
    sim.generate_fastq_file("/mnt/c/data/test", "sim", 1000000);
    /*
    int correct = 0;
    int total_size = 0;

    // Recording the distribution of size
    //std::map<int, int> sizes;

    for (int i = 0; i < 100000; i++) {
        auto sample = sim.sample();
        int bucket = std::get<1>(sample);
        std::vector<seqan3::dna4> sequence = std::get<0>(sample);
        std::vector<int> buckets = map.query_sequence(sequence);

        if (std::find(buckets.begin(), buckets.end(), bucket) != buckets.end()) {
            correct++;
        }
        total_size += buckets.size();
        //++sizes[buckets.size()];

    }
    std::cout << correct << std::endl;
    std::cout << ((float) total_size) / 100000 << std::endl;
    //for (const auto &[size, count] : sizes) {
    //    std::cout << size << "\t" << count << "\n";
    //}
    */

}