cmake_minimum_required(VERSION 3.9)
project(ArtCFD)

set(CMAKE_CXX_STANDARD 17)

add_executable(ArtCFD main.cpp DataArray.hpp DataTable.hpp Geometry.hpp GeometryMeshFemData.hpp IO_Uility.hpp ProjectUility.hpp GeometryMeshData.hpp CartesianVector.hpp DimensionUtility.hpp)