#ifndef BUCKET_MAP_MAPPER_H
#define BUCKET_MAP_MAPPER_H

class mapper {
    /**
     * @brief Virtual class for all mapper (specific to bucketMap), which map short reads
     *        to their corresponding bucket.
     */
public:
    mapper() {}

    unsigned int num_records = 0;
    
    /**
     * @brief Load the q-gram index file to the mapper.
     */
    virtual void load(std::filesystem::path const & index_file, const std::string & indicator) = 0;

    /**
     * @brief Read a query fastq file and output the ids of the sequence that are mapped 
     *        to each bucket.
     * 
     * @returns A pair of vectors, with each element in each vector representing a bucket. Each element is storing the indices
     *          of the reads that are mapped to the bucket. In the first vector stores the read that can be mapped directly,
     *          while the second vector is the reads that need to be reverse complemented.
     */
    virtual std::pair<std::vector<std::vector<unsigned int>>, std::vector<std::vector<unsigned int>>>
    map(std::filesystem::path const & sequence_file) = 0;

    /**
     * @brief Release the memory storing the sequences and index.
     */
    virtual void reset() = 0;
};

#endif