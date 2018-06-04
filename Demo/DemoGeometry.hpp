//
// Created by Chingkai Chou on 5/31/18.
//

#ifndef ARTPDE_KAVY_DEMOGEOMETRY_HPP
#define ARTPDE_KAVY_DEMOGEOMETRY_HPP

#include <iostream>
#include "Project/art_project.hpp"
#include "Geometry/geometry.hpp"

void DemoGeometry() {
    using namespace art_pde;
    using namespace PositionVector;
    using namespace project;
    using namespace geometry::mesh_type::Dim2;

    auto proj = ArtProject::create("TestProj").build();

    using MeshDataType = GeometricData<ComputePositionVector<2>>;

    using GeoType = Geometry<GeometricAlgorithm<MeshDataType>, GeometricReader<MeshDataType>>;

    auto geo = GeoType::create();

    assert(geo->read(proj)); // Check loading status.


    std::cout<< geo->getTotalNum_Vertex() << std::endl;
}

#endif //ARTPDE_KAVY_DEMOGEOMETRY_HPP
