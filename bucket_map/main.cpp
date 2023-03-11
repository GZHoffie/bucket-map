#include "indexer/bucket_hash_indexer.h"
#include "mapper/q_gram_mapper.h"
#include "locator/bucket_locator.h"
#include "tools/short_read_simulator.h"
#include "./utils.h"

#include <sharg/all.hpp>
#include <cmath>


struct cmd_arguments
{
    // running options
    bool only_indexer = false;

    std::filesystem::path fastq_path{};
    std::string index_indicator;
    std::filesystem::path output_sam_path{};

    std::string mapper_seed_shape = "111111111";
    std::string locator_seed_shape = "11111111111";

    unsigned int max_read_length = 300;

    // mapper related arguments
    unsigned int mapper_sample_size = 15;
    float mapper_distinguishability_threshold = 0.5;
    unsigned int average_base_quality = 30;

    // locator related arguments
    float allowed_seed_miss_rate = 0.4;
    float locator_allowed_indel_rate = 0.015;
    float locator_sample_size = 10;
    unsigned int locator_quality_threshold = 40;
};
 
void initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Gu Zhenhao";
    parser.info.short_description = "BucketMap, a Novel Hierarchical DNA Read Mapping Tool.";
    parser.info.version = "0.1.0";

    parser.add_flag(args.only_indexer,
                    sharg::config{.short_id = 'x',
                                  .long_id = "run-index",
                                  .description = "Only run the indexer, and store the output index in the current directory."});

    parser.add_option(args.fastq_path,
                      sharg::config{.short_id = 'q',
                                    .long_id = "query-file",
                                    .description = "The path to the FASTQ query file.",
                                    .validator = sharg::input_file_validator{{"fq", "fastq"}}});

    parser.add_option(args.index_indicator,
                      sharg::config{.short_id = 'i',
                                    .long_id = "index-indicator",
                                    .description = "The name of the index file.",
                                    .required = true});
    
    parser.add_option(args.output_sam_path,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output-file",
                                    .description = "The path to the output sam file.",
                                    .validator = sharg::output_file_validator{"sam"}});
    
    parser.add_option(args.mapper_seed_shape,
                      sharg::config{.short_id = 'm',
                                    .long_id = "mapper-seed",
                                    .description = "The shape of k-mer seed used in mapper, which assigns a read to its candidate bucket. Required to be symmetric."});
    
    parser.add_option(args.average_base_quality,
                      sharg::config{.short_id = 'b',
                                    .long_id = "average-base-quality",
                                    .description = "The average base quality in the read."});

    parser.add_option(args.locator_seed_shape,
                      sharg::config{.short_id = 'l',
                                    .long_id = "locator-seed",
                                    .description = "The shape of k-mer seed used in locator, which finds a read's exact location within the bucket. Required to be symmetric."});

    parser.add_option(args.max_read_length,
                      sharg::config{.short_id = 'r',
                                    .long_id = "read-len",
                                    .description = "The maximum length of DNA read to be mapped."});
    
    parser.add_option(args.mapper_sample_size,
                      sharg::config{.short_id = 's',
                                    .long_id = "mapper-samples",
                                    .description = "The number of k-mer samples taken by the mapper."});

    parser.add_option(args.mapper_distinguishability_threshold,
                      sharg::config{.short_id = 'd',
                                    .long_id = "distinguishability",
                                    .description = "The maximum percentage of buckets a sampled k-mer is allowed to be present."});
    
    parser.add_option(args.allowed_seed_miss_rate,
                      sharg::config{.short_id = 'e',
                                    .long_id = "max-error-rate",
                                    .description = "The maximum percentage of k-mer samples that are allowed to miss (1-epsilon)."});
    
    parser.add_option(args.locator_allowed_indel_rate,
                      sharg::config{.short_id = 'n',
                                    .long_id = "max-indel-rate",
                                    .description = "The maximum indel rate allowed."});

    parser.add_option(args.locator_sample_size,
                      sharg::config{.short_id = 'p',
                                    .long_id = "locator-samples",
                                    .description = "The number of k-mer samples taken by the locator."});
    
    parser.add_option(args.locator_quality_threshold,
                      sharg::config{.short_id = 'u',
                                    .long_id = "quality",
                                    .description = "The minimum quality score needed for an alignment to be in the output."});
}

seqan3::shape string_to_shape(const std::string& seed_shape) {
    /**
     * @brief Convert a string with binary digits to a seqan3::shape object.
     */
    std::bitset<64> binary_number(seed_shape);
    seqan3::shape res{seqan3::bin_literal{binary_number.to_ullong()}};
    return res;
}

int main(int argc, char ** argv) {
    // parse command line arguments
    sharg::parser parser{
#ifdef BM_ALIGN
        "bucketmap_align",
#else
        "bucketmap",
#endif
        argc, argv};

    cmd_arguments args{};
    initialise_parser(parser, args);

    try {
        parser.parse(); // trigger command line parsing
    }
    catch (sharg::parser_error const & ext)
    {
        seqan3::debug_stream << "[ERROR]\t\t" << ext.what() << "\n";
        return -1;
    }
 

    // check the alignment option
#ifdef BM_ALIGN
    seqan3::debug_stream << "[INFO]\t\tAllowing Smith-Waterman for alignment verifications.\n";
#else
    seqan3::debug_stream << "[INFO]\t\tNot using Smith-Waterman for alignment verifications.\n";
#endif


    std::filesystem::path data_path = "/mnt/d/genome";
    std::filesystem::path output_path = "/home/zhenhao/bucket-map/bucket_map/benchmark/output";
    seqan3::shape bucket_shape = string_to_shape(args.mapper_seed_shape);
    seqan3::shape locate_shape = string_to_shape(args.locator_seed_shape);


    // build the indexer and mapper
#if defined(BM_BUCKET_NUM) && defined(BM_BUCKET_LEN) && defined(BM_GENOME_PATH)
    seqan3::debug_stream << "[INFO]\t\tInitializing indexer and mapper with bucket length: " 
                         << BM_BUCKET_LEN << ", and number of buckets: " << BM_BUCKET_NUM << ".\n";
    
    bucket_hash_indexer<BM_BUCKET_NUM> ind(BM_BUCKET_LEN, 
                                           args.max_read_length, 
                                           bucket_shape, 
                                           locate_shape);
    if (args.only_indexer) {
        ind.index(BM_GENOME_PATH, std::filesystem::current_path(), args.index_indicator);
        return 0;
    }
    if (args.output_sam_path.empty()) {
        seqan3::debug_stream << "The output sam file is not set. Please set the output path using '-o' option.\n";
        return 1;
    }
    q_gram_mapper<BM_BUCKET_NUM> map(BM_BUCKET_LEN, 
                                     args.max_read_length, 
                                     bucket_shape, 
                                     args.mapper_sample_size, 
                                     ceil(args.mapper_sample_size * args.allowed_seed_miss_rate),
                                     args.mapper_distinguishability_threshold,
                                     args.average_base_quality);
    // initiate the locator instance
    bucket_locator loc(&ind, &map, 
                       BM_BUCKET_LEN, 
                       args.max_read_length, 
                       locate_shape, 
                       args.allowed_seed_miss_rate, 
                       args.locator_allowed_indel_rate,
                       args.locator_sample_size,
                       args.average_base_quality);
    
    // initialize the locator
    loc.initialize(BM_GENOME_PATH, std::filesystem::current_path(), args.index_indicator);
    
    // output the map results to the SAM file.
    loc.locate(args.fastq_path, std::filesystem::current_path() / (args.index_indicator + ".bucket_id"), 
               args.output_sam_path, args.locator_quality_threshold);

#else
    seqan3::debug_stream << "[ERROR]\t\tThe definition of BM_BUCKET_NUM, BM_BUCKET_LEN or BM_GENOME_FILE is not found."
                         << " Rerun the CMake build process.\n";
    return -1;
#endif

    return 0;
}