cmake_minimum_required(VERSION 3.4)
project(bucket_mapper CXX)
 
# add seqan3 to search path
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan3/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package(seqan3 3.0 REQUIRED)
 
# build app with seqan3
add_executable(bucket_map main.cpp bucket_map.h utils.h)
target_link_libraries(bucket_map seqan3::seqan3)