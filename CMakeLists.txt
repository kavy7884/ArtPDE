cmake_minimum_required(VERSION 3.9)
project(ArtCFD)

set(CMAKE_CXX_STANDARD 17)
set(EIGEN_DIR "/Users/kavy/Documents/Eigen/")
include_directories(${EIGEN_DIR})
add_executable(ArtCFD main.cpp geometry.hpp geometry_data.hpp point.hpp dimension_utility.hpp coordinate_basis_utility.hpp vertex.hpp cell.hpp geometry_data_reader.hpp geometry_data_reader_artpde.hpp geometry_data_algorithm.hpp edge.hpp)

target_link_libraries(ArtCFD ${EIGEN_DIR})
#target_link_libraries(ArtCFD)