cmake_minimum_required(VERSION 3.10)
project(ExtractRight)

find_package(Boost 1.66 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(Range-V3-VS2015/include)

set(gtest_force_shared_crt true)
add_subdirectory(googletest)

set(BENCHMARK_ENABLE_GTEST_TESTS false)
set(BENCHMARK_ENABLE_TESTING false)
add_subdirectory(benchmark)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_SUPPRESS_REGENERATION true)
add_compile_options(/EHsc /await)

add_executable(ExtractRight 
  ExtractRight/ExtractRight.cpp 
  ExtractRight/ExtractNoCheck.hpp
  ExtractRight/Common.hpp)
target_link_libraries(ExtractRight gtest benchmark)
