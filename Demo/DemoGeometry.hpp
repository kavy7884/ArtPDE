//
// Created by Chingkai Chou on 5/31/18.
//

#ifndef ARTPDE_KAVY_DEMOGEOMETRY_HPP
#define ARTPDE_KAVY_DEMOGEOMETRY_HPP

#include <iostream>
#include "Project/art_project.hpp"
#include "Geometry/geometry.hpp"

template <typename GeometricDataType>
class UserDefineReader: public virtual GeometricDataType {
public:
public:
    UserDefineReader() : GeometricDataType() {
        std::cout << "UserDefineReader" << std::endl;
    }

    void read(){
        using PointType = typename GeometricDataType::type::PointType;
        using PtrPointType = typename GeometricDataType::type::PtrPointType;
        using VertexType = typename GeometricDataType::type::VertexType_traits::VertexType;
        using PtrVertexType = typename GeometricDataType::type::VertexType_traits::PtrVertexType;

        PtrPointType ptr_pt;
        PtrVertexType ptr_vertex;

        // Pt 1
        ptr_vertex = std::make_shared<VertexType>();
        ptr_pt = std::make_shared<PointType>();
        ptr_pt->setX(1.0); ptr_pt->setY(2.0);
        ptr_vertex->setPtr_data(ptr_pt);
        this->getTotalVec_PtrVertex().push_back(ptr_vertex);

        // Pt 2
        ptr_vertex = std::make_shared<VertexType>();
        ptr_pt = std::make_shared<PointType>();
        ptr_pt->setX(2.0); ptr_pt->setY(4.0);
        ptr_vertex->setPtr_data(ptr_pt);
        this->getTotalVec_PtrVertex().push_back(ptr_vertex);
    }
};

void DemoGeometry() {
    using namespace art_pde;
    using namespace PositionVector;
    using namespace project;
    using namespace geometry::mesh_type::Dim2;

//    auto proj = ArtProject::create("TestProj_SpecialCase").build();
    auto proj = ArtProject::create("TestProj").build();

    using MeshDataType = GeometricData<ComputePositionVector<2>>;

//    using GeoType = Geometry<MeshDataType>;
    using GeoType = Geometry<GeometricAlgorithm<MeshDataType>, GeometricReader<MeshDataType>>;
//    using GeoType = Geometry<GeometricAlgorithm<MeshDataType>, UserDefineReader<MeshDataType>>;

    auto geo = GeoType::create();

    assert(geo->read(proj)); // Check loading status.
//        geo->read();

    geo->merge();

    std::cout << ">> # of Vertex : ";
    std::cout << geo->getTotalNum_Vertex() << "\n";

    std::cout << ">> # of Edge : ";
    std::cout << geo->getTotalNum_Edge() << "\n";

    std::cout << ">> # of Face : ";
    std::cout << geo->getTotalNum_Face() << "\n";

    std::cout << ">> Listing all vertex" << "\n";
    for(auto &ptr_vertex: geo->c_getTotalVec_PtrVertex()){
        std::cout << "GeoType is " << ptr_vertex->getGeometric_Type() << " ; ";
        std::cout << *ptr_vertex << "\n";
    }

    std::cout << ">> Listing all edge" << "\n";
    for(auto &ptr_edge: geo->c_getTotalVec_PtrEdge()){
        std::cout << "GeoType is " << ptr_edge->getGeometric_Type() << " ; ";
        std::cout << *ptr_edge << "\n";
    }

    std::cout << ">> Listing all face" << "\n";
    for(auto &ptr_face: geo->c_getTotalVec_PtrFace()){
        std::cout << "GeoType is " << ptr_face->getGeometric_Type() << " ; ";
        std::cout << *ptr_face << "\n";
    }
//
//    std::cout << geo->c_getTotalVec_PtrFace()[1]->getGeometric_Type() << "\n";
}

#endif //ARTPDE_KAVY_DEMOGEOMETRY_HPP
