cmake_minimum_required(VERSION 3.4)
project(bucket_mapper CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# add seqan3 to search path
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan3/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package(seqan3 3.0 REQUIRED)
 
# set common header files
set(SHARED_FILES
    utils.h 
    indexer/bucket_fm_index.h 
    indexer/bucket_hash_index.h
    indexer/fm_index.h
    mapper/q_gram_map.h 
    mapper/q_gram_bf_map.h 
    mapper/fm_index_map.h
    simulation/short_read_simulator.h)

# build app with seqan3
add_executable(bucket_map main.cpp ${SHARED_FILES})
target_link_libraries(bucket_map seqan3::seqan3)
add_executable(bucket_map_benchmark benchmark.cpp ${SHARED_FILES})
target_link_libraries(bucket_map_benchmark seqan3::seqan3)
