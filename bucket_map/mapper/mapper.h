#ifndef BUCKET_MAP_MAPPER_H
#define BUCKET_MAP_MAPPER_H

class mapper {
    /**
     * @brief Virtual class for all mapper (specific to bucketMap), which map short reads
     *        to their corresponding bucket.
     */
public:
    mapper() {}

    virtual std::vector<std::vector<int>> map(std::filesystem::path sequence_file) {
        /**
         * @brief Read a query fastq file and output the ids of the sequence that are mapped 
         *        to each bucket.
         */
    }
};

#endif