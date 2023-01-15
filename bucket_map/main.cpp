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
    std::filesystem::path output_path = "/home/zhenhao/bucket-map/bucket_map/benchmark/output";
    seqan3::shape bucket_shape(0b110101110111_shape);
    seqan3::shape locate_shape(0b111011011010111_shape);

    int bucket_length = 65536;
    int read_length = 300;

    bucket_hash_indexer<26507> ind(bucket_length, read_length, bucket_shape, locate_shape);
    q_gram_mapper<26507> map(bucket_length, read_length, bucket_shape, 30, 6, 0.5);
    bucket_locator loc(&ind, &map, bucket_length, read_length, locate_shape, 0.01 * locate_shape.count(), 0.001 * std::ranges::size(locate_shape), 10);

    //auto res = map._query_file(data_path / "test" / "sim.fastq");
    //map._check_ground_truth(res, data_path / "test" / "sim.bucket_ground_truth");

    loc.initialize(genome_file, data_path / "index", "test");
    loc.locate(data_path / "test" / "sim_illumina_1M.fastq", data_path / "index" / "index.bucket_id", output_path / "bucket_map.sam");
    
    

}