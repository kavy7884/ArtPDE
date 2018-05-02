cmake_minimum_required(VERSION 3.9)
project(ArtCFD)

set(CMAKE_CXX_STANDARD 17)
set(EIGEN_DIR "/Users/kavy/Documents/Eigen/")
include_directories(${EIGEN_DIR})
add_executable(ArtCFD main.cpp)

target_link_libraries(ArtCFD ${EIGEN_DIR})
#target_link_libraries(ArtCFD)