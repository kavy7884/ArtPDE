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

    auto proj = project::ArtProject::create("TestProj").build();

    const size_t Dim = 3;
    using MeshDataType = geometry::mesh_type::Data<Dim, PositionVector::ComputePositionVector<Dim>>;

    using GeoType = geometry::Geometry<geometry::Algorithm<Dim, MeshDataType>, geometry::Reader<Dim, MeshDataType>>;
    GeoType geo;

    std::cout<< geo.getTotalNum_Vertex() << std::endl;

    std::cout<< geo.getTotalNum_Cell() << std::endl;

    geo.read(proj);
}

#endif //ARTPDE_KAVY_DEMOGEOMETRY_HPP
