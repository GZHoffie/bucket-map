#include "indexer/bucket_hash_indexer.h"
#include "mapper/q_gram_mapper.h"
#include "tools/short_read_simulator.h"
#include "tools/hash_function_generator.h"
#include "./utils.h"

#include <map>
#include <chrono>

int main() {
    std::filesystem::path data_path = "/home/zhenhao/mapping_data";
    std::filesystem::path genome_file = data_path / "GRCh38_adjusted.fna";

    
    uint8_t query_seed_length = 12;
    uint8_t index_seed_length = 9;

    unsigned int max_read_length = 300;

    // mapper related arguments
    unsigned int mapper_sample_size = 15;
    float mapper_distinguishability_threshold = 0.5;
    unsigned int average_base_quality = 25;

    // locator related arguments
    float allowed_seed_miss_rate = 0.4;
    float locator_allowed_indel_rate = 0.02;
    float locator_sample_size = 10;
    unsigned int locator_quality_threshold = 40;

    // parameters for fracMinHash
    float frac_min_hash = 1;


    // initialize hash function for fracMinHash
    std::size_t HASH_TABLE_SIZE = 10000;
    hash_function_generator gen;
    auto min_hash_function = gen.generate(HASH_TABLE_SIZE);

    int bucket_length = 65536;
    int read_length = 150;

    bucket_hash_indexer<47000> ind(bucket_length, max_read_length, index_seed_length,
                                   min_hash_function, (unsigned int)(HASH_TABLE_SIZE * frac_min_hash));

    q_gram_mapper<47000> map(bucket_length, max_read_length, 
                                     query_seed_length, 
                                     index_seed_length,
                                     mapper_sample_size, 
                                     ceil(mapper_sample_size * allowed_seed_miss_rate),
                                     mapper_distinguishability_threshold,
                                     average_base_quality);

    //ind.index(genome_file, data_path, "index");
    map.load(data_path, "index");

    //short_read_simulator sim(bucket_length, read_length, 0.002, 0.00025, 0.00025);
    //sim.read(genome_file);
    //sim.generate_fastq_file(data_path / "test", "sim_10", 10);

    auto res = map._query_file(data_path / "test" / "sim_150.fastq");
    //seqan3::debug_stream << res << "\n";
    map._check_ground_truth(res, data_path / "test" / "sim_150.bucket_ground_truth");
    


    
    

}