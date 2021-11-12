cmake_minimum_required(VERSION 3.15)

project(3YDMoeller VERSION 2021.01.25 LANGUAGES CXX)

add_executable(basic-test)

target_sources(basic-test PRIVATE test.cpp)

target_include_directories(basic-test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_features(basic-test PRIVATE cxx_std_11)

enable_testing()
add_test(NAME "basic-test" COMMAND basic-test)