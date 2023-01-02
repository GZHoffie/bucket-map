cmake_minimum_required(VERSION 3.4)
project(bucket_mapper CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include(FetchContent)
FetchContent_Declare(dlib
    GIT_REPOSITORY https://github.com/davisking/dlib.git
    GIT_TAG        v19.24
)
FetchContent_MakeAvailable(dlib)

FetchContent_Declare(seqan3
    GIT_REPOSITORY git@github.com:seqan/seqan3.git
    GIT_TAG        "3.2.0"
)
FetchContent_MakeAvailable(seqan3)
 
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
target_link_libraries(bucket_map PUBLIC seqan3::seqan3 dlib::dlib)
add_executable(bucket_map_benchmark benchmark.cpp ${SHARED_FILES})
target_link_libraries(bucket_map_benchmark PUBLIC seqan3::seqan3 dlib::dlib)

# benchmark executables
add_executable(sam_benchmark benchmark/sam_file_analyzer.cpp)
target_link_libraries(sam_benchmark PUBLIC seqan3::seqan3)
