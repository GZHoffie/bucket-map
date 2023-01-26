#include "indexer/bucket_hash_indexer.h"
#include "mapper/q_gram_mapper.h"
#include "tools/short_read_simulator.h"
#include "./utils.h"

#include <map>
#include <chrono>

int main() {
    std::filesystem::path data_path = "/mnt/d/genome";
    std::filesystem::path genome_file = data_path / "Egu.v3.genome_f.fasta";
    seqan3::shape bucket_shape(0b110101110111_shape);
    seqan3::shape locate_shape(0b111111111_shape);

    int bucket_length = 65536;
    int read_length = 150;

    bucket_hash_indexer<26507> ind(bucket_length, read_length, bucket_shape, locate_shape);
    q_gram_mapper<26507> map(bucket_length, read_length, bucket_shape, 30, 6, 0.5);

    ind.index(genome_file, data_path, "index");
    map.load(data_path, "index");

    //short_read_simulator sim(bucket_length, read_length, 0.002, 0.00025, 0.00025);
    //sim.read(genome_file);
    //sim.generate_fastq_file(data_path / "test", "sim_150", 1000000);

    auto res = map._query_file(data_path / "test" / "sim.fastq");
    map._check_ground_truth(res, data_path / "test" / "sim.bucket_ground_truth");
    
    

}