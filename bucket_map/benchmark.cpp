#include "indexer/bucket_fm_indexer.h"
#include "mapper/q_gram_mapper.h"
#include "locator/bucket_locator.h"
#include "tools/short_read_simulator.h"
#include "./utils.h"

#include <map>
#include <chrono>

int main() {
    std::filesystem::path data_path = "/mnt/d/genome";
    std::filesystem::path genome_file = data_path / "Egu.v3.genome_f.fasta";
    seqan3::shape shape(0b1110100101001101_shape);

    int bucket_length = 65536;
    int read_length = 150;

    bucket_fm_indexer<26507> ind(bucket_length, read_length, shape);
    q_gram_mapper<26507> map(bucket_length, read_length, shape, 20, 5, 0.7);
    bucket_locator loc(&ind, &map, bucket_length, read_length, shape, 0.002 * 9, 0.00025 * 16);

    //ind.index(genome_file, data_path / "index_fm");
    //map.read(genome_file);
    //map.store(data_path / "index");
    //map.load(data_path / "index_fm");

    short_read_simulator sim(bucket_length, read_length, 0.002, 0.00025, 0.00025);
    sim.read(genome_file);
    sim.generate_fastq_file(data_path / "test", "sample", 10);

    //auto res = map._query_file(data_path / "test" / "sim.fastq");
    //map._check_ground_truth(res, data_path / "test" / "sim.ground_truth");

    loc.initialize(genome_file, data_path / "index_fm", "test");
    auto res = loc._locate(data_path / "test" / "sample.fastq");
    loc._check_ground_truth(res, data_path / "test" / "sample.ground_truth");
    
    

}