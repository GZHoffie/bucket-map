cmake_minimum_required(VERSION 3.4)
project(bucket_mapper CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Set environment variables
set(BUCKETMAP_GENOME_FILE ${BM_FASTA_FILE})
set(BUCKETMAP_BUCKET_LEN ${BM_BUCKET_LEN})

message("[INFO]\t\tThe path to the reference genome file: ${BUCKETMAP_GENOME_FILE}.")
message("[INFO]\t\tBucket length is set to be ${BUCKETMAP_BUCKET_LEN}.")

execute_process(
    COMMAND bash -c "
        awk -v bucket_l=${BUCKETMAP_BUCKET_LEN} '
        function ceil(value)
        {
            return (value == int(value)) ? value : int(value)+1
        }

        /^>/ { # header pattern detected
            if (seqlen){
                # print previous seqlen if exists 
                buckets += ceil(seqlen/bucket_l)
            }

            # initialize sequence
            seqlen = 0

            # skip further processing
            next
        }

        # accumulate sequence length
        {
            seqlen += length($0)
        }

        # remnant seqlen if exists
        END{
            if(seqlen){
                buckets += ceil(seqlen/bucket_l)
            }
            print buckets
        }
    ' ${BUCKETMAP_GENOME_FILE}
    "
    OUTPUT_VARIABLE BUCKETMAP_BUCKET_NUM
)
string(STRIP ${BUCKETMAP_BUCKET_NUM} BUCKETMAP_BUCKET_NUM)

message("[INFO]\t\tThe number of buckets is set to be ${BUCKETMAP_BUCKET_NUM}.")

# Set the global variables
set(BM_DEFINITIONS
    BM_BUCKET_NUM=${BUCKETMAP_BUCKET_NUM} 
    BM_GENOME_PATH=\"${BUCKETMAP_GENOME_FILE}\"
    BM_BUCKET_LEN=${BUCKETMAP_BUCKET_LEN})

include(FetchContent)
# dlib for serialization of the index files
#FetchContent_Declare(dlib
#    GIT_REPOSITORY https://github.com/davisking/dlib.git
#    GIT_TAG        v19.24
#)
#FetchContent_MakeAvailable(dlib)

# seqan3 for the basic skeleton for the mapper
FetchContent_Declare(seqan3
    GIT_REPOSITORY https://github.com/seqan/seqan3.git
    GIT_TAG        cdb6abe5e36b958e86ff1ba94561dac4955d5f14
)
FetchContent_MakeAvailable(seqan3)

# Sharg for argument parsing
FetchContent_Declare(sharg
    GIT_REPOSITORY https://github.com/seqan/sharg-parser.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(sharg)


#FetchContent_Declare(parasail
#  GIT_REPOSITORY https://github.com/jeffdaily/parasail.git
#  GIT_TAG        master
#)
#FetchContent_MakeAvailable(parasail)

# set common header files
set(INDEXER_FILES
    indexer/indexer.h
    indexer/bucket_indexer.h
    indexer/bucket_fm_indexer.h
    indexer/fm_indexer.h)

set(LOCATOR_FILES
    locator/locator.h
    locator/bucket_locator.h
    locator/fm_index_locator.h)

set(MAPPER_FILES
    mapper/mapper.h
    mapper/q_gram_mapper.h
    mapper/fm_index_mapper.h
    mapper/quality_filter.h)

set(UTIL_FILES
    utils.h
    tools/short_read_simulator.h
    tools/hash_function_generator.h)

set(SHARED_FILES
    ${UTIL_FILES}
    ${INDEXER_FILES}
    ${MAPPER_FILES}
    ${LOCATOR_FILES})

# build dataset generating tools
add_executable(generate_dataset generate_dataset.cpp ${UTIL_FILES})
target_link_libraries(generate_dataset PUBLIC seqan3::seqan3)

# build app with seqan3
add_executable(mapper_test mapper_test.cpp ${SHARED_FILES})
target_link_libraries(mapper_test PUBLIC seqan3::seqan3 cuco)


add_executable(fastq_analyzer benchmark/fastq_analyzer.cpp ${SHARED_FILES})
target_link_libraries(fastq_analyzer PUBLIC seqan3::seqan3)

#add_executable(quality_filter_test mapper/quality_filter.cpp ${SHARED_FILES})
#target_link_libraries(quality_filter_test PUBLIC seqan3::seqan3)

add_executable(bucketmap main.cpp ${SHARED_FILES})
target_compile_definitions(bucketmap PRIVATE ${BM_DEFINITIONS})
target_link_libraries(bucketmap PUBLIC seqan3::seqan3 sharg::sharg cuco)

add_executable(bucketmap_align main.cpp ${SHARED_FILES})
target_compile_definitions(bucketmap_align PRIVATE BM_ALIGN ${BM_DEFINITIONS})
target_link_libraries(bucketmap_align PUBLIC seqan3::seqan3 sharg::sharg)

#get_property(defs TARGET bucketmap_align PROPERTY COMPILE_DEFINITIONS)
#message("Defs: ${defs}")

# benchmark executables
add_executable(sam_benchmark benchmark/sam_file_analyzer.cpp)
target_link_libraries(sam_benchmark PUBLIC seqan3::seqan3)
