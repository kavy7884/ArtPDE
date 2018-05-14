cmake_minimum_required(VERSION 3.9)
project(ArtPDE)

set(CMAKE_CXX_STANDARD 11)

#set(EIGEN_DIR "/Users/Jeting/Documents/eigen/")
#include_directories(${EIGEN_DIR})

include_directories(Utility)
include_directories(Geometry)
include_directories(FunctionSpace)

#add_executable(ArtPDE main_Jeting.cpp)
#target_link_libraries(ArtPDE ${EIGEN_DIR})

add_executable(ArtPDE_jeting main_Jeting.cpp)
add_executable(ArtPDE_kavy main_kavy.cpp)
add_executable(ArtPDE_marvin main_marvin.cpp)

target_link_libraries(ArtPDE_jeting)
target_link_libraries(ArtPDE_kavy)
target_link_libraries(ArtPDE_marvin)