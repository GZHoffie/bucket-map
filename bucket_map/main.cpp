#include "indexer/bucket_hash_indexer.h"
#include "mapper/q_gram_mapper.h"
#include "locator/bucket_locator.h"
#include "tools/short_read_simulator.h"
#include "./utils.h"

#include <map>
#include <chrono>

int main() {

    // check the alignment option
#ifdef BM_ALIGN
    seqan3::debug_stream << "[INFO]\t\tAllowing Smith-Waterman for alignment verifications.\n";
#else
    seqan3::debug_stream << "[INFO]\t\tNot using Smith-Waterman for alignment verifications.\n";
#endif

    std::filesystem::path data_path = "/mnt/d/genome";
    std::filesystem::path output_path = "/home/zhenhao/bucket-map/bucket_map/benchmark/output";
    seqan3::shape bucket_shape(0b110101110111_shape);
    seqan3::shape locate_shape(0b111011011010111_shape);

    int bucket_length = 65536;
    int read_length = 300;

#if defined(BM_BUCKET_NUM) && defined(BM_BUCKET_LEN)
    // build the indexer and mapper
    seqan3::debug_stream << "[INFO]\t\tInitializing indexer and mapper with bucket length: " 
                         << BM_BUCKET_LEN << ", and number of buckets: " << BM_BUCKET_NUM << ".\n";
    
    bucket_hash_indexer<BM_BUCKET_NUM> ind(BM_BUCKET_LEN, read_length, bucket_shape, locate_shape);
    q_gram_mapper<BM_BUCKET_NUM> map(BM_BUCKET_LEN, read_length, bucket_shape, 30, 6, 0.5);
#else
    seqan3::debug_stream << "[ERROR]\t\tThe definition of BM_BUCKET_NUM or BM_BUCKET_LEN is not found."
                         << " Rerun the CMake build process.\n";
    return -1;
#endif

    // initiate the locator instance
    bucket_locator loc(&ind, &map, bucket_length, read_length, locate_shape, 0.01 * locate_shape.count(), 0.001 * std::ranges::size(locate_shape), 10);

#if defined(BM_GENOME_PATH)
    // initialize the locator
    loc.initialize(BM_GENOME_PATH, data_path / "index", "test");
#else
    seqan3::debug_stream << "[ERROR]\t\tThe definition of BM_GENOME_FILE is not found. "
                         << "Rerun the CMake build process.\n";
    return -1;
#endif

    // output the map results to the SAM file.
    loc.locate(data_path / "test" / "sim_illumina_1M.fastq", data_path / "index" / "index.bucket_id", output_path / "bucket_align_map.sam");
}