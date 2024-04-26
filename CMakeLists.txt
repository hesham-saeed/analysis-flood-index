cmake_minimum_required(VERSION 3.1.0)

project(Flood VERSION 1.00 LANGUAGES CXX)

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_COMPILE_FLAGS             "-W -Wall -pedantic -DLEVEL1_DCACHE_LINESIZE=${LEVEL1_DCACHE_LINESIZE} -DPAGESIZE=${PAGESIZE} -march=native -Wno-variadic-macros -Wno-gnu-zero-variadic-macro-arguments -Wno-gnu-label-as-value -Wno-vla-extension")
SET(CMAKE_C_FLAGS                   "${CMAKE_C_FLAGS} ${CMAKE_COMPILE_FLAGS}")
SET(CMAKE_CXX_FLAGS                 "-std=c++17 ${CMAKE_CXX_FLAGS} ${CMAKE_COMPILE_FLAGS}")
SET(CMAKE_CXX_FLAGS_DEBUG           "-ggdb3 -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize=address,undefined -fsanitize-address-use-after-scope")
SET(CMAKE_CXX_FLAGS_RELEASE         "-O2")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O2 -ggdb3")


include_directories(include)
include_directories(examples)
include_directories(third_party/boost/boost)
include_directories(third_party/grid_files)

add_executable(idx examples/example.cpp)
add_subdirectory(experiments)
