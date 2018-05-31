//
// Created by Chingkai Chou on 5/31/18.
//

#ifndef ARTPDE_KAVY_DEMOGEOMETRY_HPP
#define ARTPDE_KAVY_DEMOGEOMETRY_HPP

#include <iostream>
#include "Geometry/geometry.hpp"

void DemoGeometry() {
    using namespace art_pde;

    using MeshDataType = geometry::mesh_type::Data<3, PositionVector::ComputePositionVector>;
    using GeoType = geometry::Geometry<geometry::Algorithm<MeshDataType>, geometry::Reader<MeshDataType>>;
    GeoType geo;

    std::cout<< geo.vertex << std::endl;
}

#endif //ARTPDE_KAVY_DEMOGEOMETRY_HPP
