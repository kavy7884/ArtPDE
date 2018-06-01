//
// Created by Chingkai Chou on 5/31/18.
//

#ifndef ARTPDE_KAVY_DEMOGEOMETRY_HPP
#define ARTPDE_KAVY_DEMOGEOMETRY_HPP

#include <iostream>
#include "Geometry/geometry.hpp"

void DemoGeometry() {
    using namespace art_pde;

    const size_t Dim = 3;
    using MeshDataType = geometry::mesh_type::Data<Dim, PositionVector::ComputePositionVector>;
    using GeoType = geometry::Geometry<geometry::Algorithm<Dim, MeshDataType>, geometry::Reader<MeshDataType>>;
    GeoType geo;

    std::cout<< geo.getTotalNum_Vertex() << std::endl;

    std::cout<< geo.getTotalNum_Cell() << std::endl;
}

#endif //ARTPDE_KAVY_DEMOGEOMETRY_HPP
