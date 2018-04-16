cmake_minimum_required(VERSION 3.10)
project(ExtractRight)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -fcoroutines-ts -Wno-logical-op-parentheses")
endif()

find_package(Boost 1.65 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

if(MSVC)
  include_directories(Range-V3-VS2015/include)
else()
  include_directories(range-v3/include)
endif()

set(gtest_force_shared_crt ON CACHE BOOL "")
add_subdirectory(googletest)

set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "")
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "")
add_subdirectory(benchmark)

include_directories(cppcoro/include)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_SUPPRESS_REGENERATION true)
if(MSVC)
  add_compile_options(/EHsc /await)
endif()
add_definitions(-DBOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE -D_SILENCE_CXX17_RESULT_OF_DEPRECATION_WARNING)

add_executable(ExtractRight 
  ExtractRight/ExtractRight.cpp 
  ExtractRight/ExtractNoCheck.hpp
  ExtractRight/Common.hpp)
target_link_libraries(ExtractRight gtest benchmark)
