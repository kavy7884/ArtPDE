cmake_minimum_required(VERSION 3.9)
project(ArtCFD)

set(CMAKE_CXX_STANDARD 17)
set(EIGEN_DIR "/Users/kavy/Documents/Eigen/")
include_directories(${EIGEN_DIR})
add_executable(ArtCFD main.cpp DataArray.hpp DataTable.hpp Geometry.hpp IO_Uility.hpp ProjectUility.hpp CartesianVector.hpp DimensionUtility.hpp Domain.hpp Integration.hpp Operation.hpp DOF.hpp TrialSpace.hpp Approximation.hpp TestSpace.hpp BoundaryCondition.hpp NumericalMethodUtility.hpp GeometryData.hpp ShapeFunction.hpp DeformationMapping.hpp System.hpp DOF_Mannger.hpp Solver.hpp Point.hpp GeometryLoader.hpp)

target_link_libraries(ArtCFD ${EIGEN_DIR})
#target_link_libraries(ArtCFD)