cmake_minimum_required(VERSION 3.10)

project(NumericalMethods)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(include)
file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(NumericalMethods ${SRC_FILES})

set(TESTS_DIR "${CMAKE_SOURCE_DIR}/tests/")
set(RES_DIR "${CMAKE_SOURCE_DIR}/results/")


target_compile_definitions(NumericalMethods PRIVATE TESTS_DIR="${TESTS_DIR}")
target_compile_definitions(NumericalMethods PRIVATE RES_DIR="${RES_DIR}")