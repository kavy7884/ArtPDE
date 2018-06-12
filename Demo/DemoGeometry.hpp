//
// Created by Chingkai Chou on 5/31/18.
//

#ifndef ARTPDE_KAVY_DEMOGEOMETRY_HPP
#define ARTPDE_KAVY_DEMOGEOMETRY_HPP

#include <iostream>
#include <assert.h>
#include "Project/art_project.hpp"
#include "Geometry/geometry.hpp"

void DemoGeometry2D() {
    using namespace art_pde;
    using namespace position_vector;
    using namespace project;
    using namespace geometry::mesh_type::Dim2;

    // Step 1: Open the ArtProject!
    //- It will generate project folder tree under your run path (the default is: ./) automatically.
    //- (Four sub-folder: Geometry, Initial, Setting, Results generated simultaneously.)
    //-- e.g.: ArtProject::create("TestProj2D").setRunPath(".").build();
    //-- The project will generate: ./TestProj2D
    //-- e.g.: ArtProject::create("TestProj2D").setRunPath("/User/Kavy/Desktop").build();
    //-- The project will generate: /User/Kavy/Desktop/TestProj2D

    auto proj = ArtProject::create("TestProj2D").setRunPath(".").build();

    // Step 2: Defining mesh data and which Algorithm and Reader!

    using MeshDataType = GeometricData<ComputePositionVector<2>>;
    using GeoType = Geometry<GeometricAlgorithm<MeshDataType>, GeometricReader<MeshDataType>>;
    auto geo = GeoType::create();

    // Step 3: Reading data from ArtProject's geometry folder.
    //- two files: position.txt and connectivity.txt must existing and defining properly.
    //- if the loading process failure, the program will exit and show something error in this line.
    //- P.S.: the reading method provided by "GeometricReader" API class.
	geo->read(proj); // Check loading status.

    // Step 4: Merge geometry
    //- 2D case will merge Edge.
    geo->merge();

    // Step 5: Check all the geometry 2D data
    std::cout << ">> # of Vertex : ";
    std::cout << geo->getTotalNum_Vertex() << "\n";

    std::cout << ">> # of Edge : ";
    std::cout << geo->getTotalNum_Edge() << "\n";

    std::cout << ">> # of Face : ";
    std::cout << geo->getTotalNum_Face() << "\n";

    std::cout << ">> Listing all relations of vertex" << "\n";
    for(auto &ptr_vertex: geo->c_getTotalVec_PtrVertex()){
        std::cout << "GeoType is " << ptr_vertex->getGeometric_Type() << " ; ";
        std::cout << *ptr_vertex << "\n";
    }

    std::cout << ">> Listing all relations of edge" << "\n";
    for(auto &ptr_edge: geo->c_getTotalVec_PtrEdge()){
        std::cout << "GeoType is " << ptr_edge->getGeometric_Type() << " ; ";
        std::cout << *ptr_edge << "\n";
    }

    std::cout << ">> Listing all relations of face" << "\n";
    for(auto &ptr_face: geo->c_getTotalVec_PtrFace()){
        std::cout << "GeoType is " << ptr_face->getGeometric_Type() << " ; ";
        std::cout << *ptr_face << "\n";

        auto vec_ptr_vertex = ptr_face->c_getConnected_Vertex();
        std::cout << "Final face's vertex is : \n";
        for(auto &ptr_vertex: vec_ptr_vertex){
            std::cout << *ptr_vertex << std::endl;
        }
    }
}

void DemoGeometry3D() {
    using namespace art_pde;
    using namespace position_vector;
    using namespace project;
    using namespace geometry::mesh_type::Dim3;

    // Step 1: Open the ArtProject!
    //- It will generate project folder tree under your run path (the default is: ./) automatically.
    //- (Four sub-folder: Geometry, Initial, Setting, Results generated simultaneously.)
    //-- e.g.: ArtProject::create("TestProj3D").setRunPath(".").build();
    //-- The project will generate: ./TestProj3D
    //-- e.g.: ArtProject::create("TestProj3D").setRunPath("/User/Kavy/Desktop").build();
    //-- The project will generate: /User/Kavy/Desktop/TestProj3D

    auto proj = ArtProject::create("TestProj3D").build();

    // Step 2: Defining mesh data and which Algorithm and Reader!

    using MeshDataType = GeometricData<ComputePositionVector<3>>;
    using GeoType = Geometry<GeometricAlgorithm<MeshDataType>, GeometricReader<MeshDataType>>;
    auto geo = GeoType::create();

    // Step 3: Reading data from ArtProject's geometry folder.
    //- two files: position.txt and connectivity.txt must existing and defining properly.
    //- if the loading process failure, the program will exit and show something error in this line.
    //- P.S.: the reading method provided by "GeometricReader" API class.

    geo->read(proj); // Check loading status.

    // Step 4: Merge geometry
    //- 3D case will merge Edge and Face.
    geo->merge();

    // Step 5: Check all the geometry 2D data
    std::cout << ">> # of Vertex : ";
    std::cout << geo->getTotalNum_Vertex() << "\n";

    std::cout << ">> # of Edge : ";
    std::cout << geo->getTotalNum_Edge() << "\n";

    std::cout << ">> # of Face : ";
    std::cout << geo->getTotalNum_Face() << "\n";

    std::cout << ">> # of Cell : ";
    std::cout << geo->getTotalNum_Cell() << "\n";

    std::cout << ">> Listing all relations of vertex" << "\n";
    for(auto &ptr_vertex: geo->c_getTotalVec_PtrVertex()){
        std::cout << "GeoType is " << ptr_vertex->getGeometric_Type() << " ; ";
        std::cout << *ptr_vertex << "\n";
    }

    std::cout << ">> Listing all relations of edge" << "\n";
    for(auto &ptr_edge: geo->c_getTotalVec_PtrEdge()){
        std::cout << "GeoType is " << ptr_edge->getGeometric_Type() << " ; ";
        std::cout << *ptr_edge << "\n";
    }

    std::cout << ">> Listing all relations of face" << "\n";
    for(auto &ptr_face: geo->c_getTotalVec_PtrFace()){
        std::cout << "GeoType is " << ptr_face->getGeometric_Type() << " ; ";
        std::cout << *ptr_face << "\n";
    }

    std::cout << ">> Listing all relations of cell" << "\n";
    for(auto &ptr_cell: geo->c_getTotalVec_PtrCell()) {
        std::cout << "GeoType is " << ptr_cell->getGeometric_Type() << " ; ";
        std::cout << *ptr_cell << "\n";

        auto vec_ptr_vertex = ptr_cell->c_getConnected_Vertex();
        std::cout << "Final cell's vertex is : \n";
        for(auto &ptr_vertex: vec_ptr_vertex){
            std::cout << *ptr_vertex << std::endl;
        }
    }
}


//template <typename GeometricDataType>
//class UserDefineReader: public virtual GeometricDataType {
//public:
//public:
//    UserDefineReader() : GeometricDataType() {
//        std::cout << "UserDefineReader" << std::endl;
//    }
//
//    void read(){
//        using PointType = typename GeometricDataType::type::PointType;
//        using PtrPointType = typename GeometricDataType::type::PtrPointType;
//        using VertexType = typename GeometricDataType::type::VertexType_traits::VertexType;
//        using PtrVertexType = typename GeometricDataType::type::VertexType_traits::PtrVertexType;
//
//        PtrPointType ptr_pt;
//        PtrVertexType ptr_vertex;
//
//        // Pt 1
//        ptr_vertex = std::make_shared<VertexType>();
//        ptr_pt = std::make_shared<PointType>();
//        ptr_pt->setX(1.0); ptr_pt->setY(2.0);
//        ptr_vertex->setPtr_data(ptr_pt);
//        this->getTotalVec_PtrVertex().push_back(ptr_vertex);
//
//        // Pt 2
//        ptr_vertex = std::make_shared<VertexType>();
//        ptr_pt = std::make_shared<PointType>();
//        ptr_pt->setX(2.0); ptr_pt->setY(4.0);
//        ptr_vertex->setPtr_data(ptr_pt);
//        this->getTotalVec_PtrVertex().push_back(ptr_vertex);
//    }
//};

#endif //ARTPDE_KAVY_DEMOGEOMETRY_HPP
