cmake_minimum_required(VERSION 3.4)
project(bucket_mapper CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# add seqan3 to search path
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan3/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package(seqan3 3.0 REQUIRED)
 
# set common header files
set(INDEXER_FILES
    indexer/indexer.h
    indexer/bucket_fm_indexer.h 
    indexer/bucket_indexer.h
    indexer/fm_indexer.h)

set(LOCATOR_FILES
    locator/locator.h
    locator/bucket_locator.h
    locator/fm_index_locator.h)

set(MAPPER_FILES
    mapper/mapper.h
    mapper/q_gram_mapper.h 
    mapper/q_gram_bf_mapper.h 
    mapper/fm_index_mapper.h)

set(UTIL_FILES
    utils.h
    tools/short_read_simulator.h
    tools/hash_function_generator.h)

set(SHARED_FILES
    ${UTIL_FILES}
    ${INDEXER_FILES}
    ${MAPPER_FILES}
    ${LOCATOR_FILES})

# build app with seqan3
add_executable(bucket_map main.cpp ${SHARED_FILES})
target_link_libraries(bucket_map seqan3::seqan3)
add_executable(bucket_map_benchmark benchmark.cpp ${SHARED_FILES})
target_link_libraries(bucket_map_benchmark seqan3::seqan3)
