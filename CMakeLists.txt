cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 11)

set(EIGEN_DIR "/Users/kavy/Documents/Eigen/")
include_directories(${EIGEN_DIR})

# #Include header zone
include_directories(Utility)
include_directories(Geometry)
include_directories(FunctionSpace)

# Run project zone
project(ArtPDE_kavy)
add_executable(ArtPDE_kavy main_kavy.cpp)

project(ArtPDE_marvin)
add_executable(ArtPDE_marvin main_marvin.cpp)