cmake_minimum_required(VERSION 3.9)

project(ArtPDE_jeting)
project(ArtPDE_marvin)
project(ArtPDE_kavy)

set(CMAKE_CXX_STANDARD 11)

include_directories(${CMAKE_SOURCE_DIR})

#add_definitions(-D_MSC_VER)

file(GLOB GeometrySource "Project/*.hpp")
file(GLOB GeometrySource "BasicData/*.hpp")
file(GLOB GeometrySource "Geometry/*.hpp")
file(GLOB GeometrySource "Geometry/inc/GeometricTree/*.hpp")
file(GLOB GeometrySource "Geometry/inc/GeometricData/*.hpp")
file(GLOB GeometrySource "Geometry/inc/GeometricAlgorithm/*.hpp")
file(GLOB GeometrySource "Geometry/inc/GeometricReader/*.hpp")
file(GLOB FunctionSpaceSource "FunctionSpace_Refactor/*.hpp")

set(AllSource ${GeometrySource} ${FunctionSpaceSource} )

add_executable(ArtPDE_jeting main_Jeting.cpp ${AllSource})
add_executable(ArtPDE_marvin main_marvin.cpp ${AllSource})
add_executable(ArtPDE_kavy main_kavy.cxx ${AllSource})