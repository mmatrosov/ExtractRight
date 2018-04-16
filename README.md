# ExtractRight
Demo repository for Versatile C++ talk

# Installation

## Ubuntu 18.04

Set up clang

    sudo apt install clang
    sudo apt install libc++-dev
    sudo apt install libc++abi-dev

Install boost

    sudo apt install libboost-dev

Build sources

    git clone --recurse-submodules https://github.com/mmatrosov/ExtractRight.git
    cd ExtractRight
    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..

## Run

To save benchmark results to a .csv file, run as follows:

    ./ExtractRight --benchmark_out_format=csv --benchmark_out=benchmark.csv
