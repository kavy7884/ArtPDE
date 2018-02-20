cmake_minimum_required(VERSION 3.9)
project(ArtCFD)

set(CMAKE_CXX_STANDARD 11)

add_executable(ArtCFD main.cpp DataArray.hpp DataTable.hpp)