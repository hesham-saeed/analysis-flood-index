cmake_minimum_required(VERSION 3.1.0)

project(Flood VERSION 1.00 LANGUAGES CXX)

include_directories(include)
include_directories(examples)

add_executable(idx examples/example.cpp)

