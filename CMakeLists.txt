cmake_minimum_required(VERSION 3.9)

project(ArtPDE_jeting)
project(ArtPDE_kavy)
project(ArtPDE_marvin)

set(CMAKE_CXX_STANDARD 11)

include_directories(${CMAKE_SOURCE_DIR})

#file(GLOB GeometrySource "Geometry/*.hpp")
file(GLOB GeometrySource "Geometry_New/*.hpp")
file(GLOB FunctionSpaceSource "FunctionSpace/*.hpp")

set(AllSource ${GeometrySource} ${FunctionSpaceSource} Geometry_New/geo_data.hpp Geometry_New/geo_data_factory.hpp Geometry_New/GeoTest2D.hpp)

add_executable(ArtPDE_jeting main_Jeting.cpp ${AllSource})
add_executable(ArtPDE_marvin main_marvin.cpp ${AllSource})
add_executable(ArtPDE_kavy main_kavy.cpp ${AllSource})