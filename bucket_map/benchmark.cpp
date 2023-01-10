#include "indexer/bucket_hash_indexer.h"
#include "mapper/q_gram_mapper.h"
#include "locator/bucket_locator.h"
#include "tools/short_read_simulator.h"
#include "./utils.h"

#include <map>
#include <chrono>

int main() {
    std::filesystem::path data_path = "/mnt/d/genome";
    std::filesystem::path genome_file = data_path / "Egu.v3.genome_f.fasta";
    seqan3::shape bucket_shape(0b11011010101101_shape);
    seqan3::shape locate_shape(0b111111111_shape);

    int bucket_length = 65536;
    int read_length = 300;

    bucket_hash_indexer<26507> ind(bucket_length, read_length, bucket_shape, locate_shape);
    q_gram_mapper<26507> map(bucket_length, read_length, bucket_shape, 40, 8, 0.5);
    bucket_locator loc(&ind, &map, bucket_length, read_length, locate_shape, 0.01 * locate_shape.count(), 0.001 * std::ranges::size(locate_shape), 20);

    //ind.index(genome_file, data_path / "index_fm");
    //map.read(genome_file);
    //map.store(data_path / "index");
    //map.load(data_path / "index_10-mer");
    bool create_new_sample = true;

    if (create_new_sample) {
        short_read_simulator sim(bucket_length, read_length, 0.002, 0.00025, 0.00025);
        sim.read(genome_file);
        sim.generate_fastq_file(data_path / "test", "sim", 1000000);
    }
    

    //auto res = map._query_file(data_path / "test" / "sim.fastq");
    //map._check_ground_truth(res, data_path / "test" / "sim.bucket_ground_truth");

    loc.initialize(genome_file, data_path / "index", "test");
    auto res = loc._locate(data_path / "test" / "sim.fastq");
    loc._check_ground_truth(res, data_path / "test" / "sim.bucket_ground_truth");
    
    

}