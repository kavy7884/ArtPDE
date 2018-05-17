//
// Created by Chingkai Chou on 5/2/18.
//
#include <iostream>
#include <set>
#include "Utility/project_uility.hpp"
#include "Geometry/geometry.hpp"


int main() {


    auto proj_1 = ArtProject::create("TestProj").setRunPath(".").setDivideSlash("/").setInitialFolderName("Init").build();

    using GeoDataType =  art_pde::GeometryData<art_pde::MeshTypeMethod, art_pde::Dim2D, art_pde::CartesianCoordinate>;

    art_pde::Geometry<GeoDataType, art_pde::GeometryDataReaderArtPDE> geo;

    std::cout << "Vertex Num (before): " << geo.getNum_TotalVertex() << std::endl;
    std::cout << "Cell Num (before): " << geo.getNum_TotalCell() << std::endl;

    std::cout << geo.read(proj_1) << std::endl;

    std::cout << "Vertex Num (after): " << geo.getNum_TotalVertex() << std::endl;
    std::cout << "Cell Num (after): " << geo.getNum_TotalCell() << std::endl;


    std::cout << "List all vertex points: " << std::endl;
    for (size_t i = 0; i < geo.getNum_TotalVertex(); ++i) {
        std::cout << *geo.getVertex_PtrPoint(i) << std::endl;
    }

    std::cout << "List all cell center points: " << std::endl;
    for (size_t i = 0; i < geo.getNum_TotalCell(); ++i) {
        std::cout << *geo.getCell_Center_PtrPoint(i) << std::endl;
    }

    std::cout << "List all vertex points in each cell: " << std::endl;
    for (size_t i = 0; i < geo.getNum_TotalCell(); ++i) {
        std::cout << "Cell Type: " << GeoDataType::Type::GeoCellType::convertCellTypeInString(
                geo.getCell_CellDefineType(i) ) << " -> \t";
        auto vec_ptr_vertex_on_cell = geo.getRelation_Cell_Neighbor_Vertex(i);
        for(auto &ptr_vertex : vec_ptr_vertex_on_cell){
            std::cout << ptr_vertex->getPoint() << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << "List all vertex cell neighbor's center and type: " << std::endl;
    for (size_t i = 0; i < geo.getNum_TotalVertex(); ++i) {
        std::cout << "Cell " << i << " is: " << std::endl;
        auto vec_neighbor_cell = geo.getRelation_Vertex_Neighbor_Cell(i);
        for(auto &ptr_cell : vec_neighbor_cell){
            std::cout << ">>>> Cell Type is: ";
            std::cout << GeoDataType::Type::GeoCellType::convertCellTypeInString(ptr_cell->getCell_define_Type());
            std::cout << ", Cell center is: " << *ptr_cell->getPtr_cell_center_point()<< std::endl;
        }
    }

    std::cout << "List all edge center: " << std::endl;
    for (size_t i = 0; i < geo.getNum_TotalEdge(); ++i) {
        std::cout << *geo.getEdge_Center_PtrPoint(i) << std::endl;

    }


//
//    std::cout << "List all vertex connected edge number: " << std::endl;
//    for (size_t i = 0; i < geo.getTotal_VertexNum(); ++i) {
//        auto &ptr_vertex = geo.getVertex_PtrVertex(i);
//        auto &list_ptr_neighbor_edge = ptr_vertex->c_getList_ptr_neighbor_edge();
//        std::cout << "Edge num:  "<< list_ptr_neighbor_edge.size() << std::endl;
//        std::cout << ">> Edge neighbor cell num: " << std::endl;
//        for (auto & v: list_ptr_neighbor_edge) {
//            std::cout << v->getVec_ptr_neighbor_cell().size() << std::endl;
//        }
//
//    }


//    using GeoDataType =  art_pde::GeometryData<art_pde::MeshTypeMethod, art_pde::Dim2D, art_pde::CartesianCoordinate>;
//    using PointType = GeoDataType::Type::GeoPointType;
//    GeoDataType::Type::PtrGeoVertexType vertex_1, vertex_2, vertex_3, vertex_4;
//    GeoDataType::Type::PtrGeoEdgeType edge_1, edge_2, edge_3;
//    vertex_1 = std::make_shared<GeoDataType::Type::GeoVertexType>(PointType(0.0, 0.0));
//    vertex_2 = std::make_shared<GeoDataType::Type::GeoVertexType>(PointType(1.0, 0.0));
//    vertex_3 = std::make_shared<GeoDataType::Type::GeoVertexType>(PointType(0.0, -1.0));
//    vertex_4 = std::make_shared<GeoDataType::Type::GeoVertexType>(PointType(1.0, -1.0));
//
//    edge_1 = std::make_shared<art_pde::LineEdge<PointType>>(vertex_1, vertex_2);
//    edge_2 = std::make_shared<art_pde::LineEdge<PointType>>(vertex_2, vertex_1);
//    edge_3 = std::make_shared<art_pde::LineEdge<PointType>>(vertex_3, vertex_4);
//
//    std::cout << (*edge_1 == *edge_2) << std::endl;
//    std::cout << (*edge_2 == *edge_3) << std::endl;

//    std::shared_ptr<int> p1, p2;
//    std::set<std::shared_ptr<int>> s1, s2;
//    p1 = std::make_shared<int>();
//    p2 = std::make_shared<int>();
//
//    s1.insert(p1);
//    s2.insert(p2);
//    std::equal(s1.begin(), s1.end(), s2.begin());





//    std::shared_ptr<int> ptr, ptr1, ptr2;
//    std::set<std::shared_ptr<int>> set_ptr;
//
//    ptr = std::make_shared<int>(0);
//    ptr1 = ptr;
//    ptr2 = std::make_shared<int>(1);
//
//    std::cout << &ptr << std::endl;
//    std::cout << ptr << std::endl;
//    std::cout << &(*ptr) << std::endl;
//    std::cout << *ptr << std::endl;
//
//    set_ptr.insert(ptr);
//    set_ptr.insert(ptr1);
//    set_ptr.insert(ptr2);
//
//
//    std::cout << set_ptr.size() << std::endl;
//
//    std::cout << (ptr == ptr1) << std::endl;
//    std::cout << (ptr == ptr2) << std::endl;





//    using PointType = art_pde::Point<art_pde::Dim2D, art_pde::CartesianCoordinate>;
//    PointType tmp_pt;
//    tmp_pt += (PointType(4.0,8.0));
//    tmp_pt /= 2.0;
//    std::cout << tmp_pt << std::endl;



//    using PointType = art_pde::Point<art_pde::Dim2D, art_pde::CartesianCoordinate>;
//    using VertexType = art_pde::Vertex<PointType>;
//    using PrtVertexType = std::shared_ptr<VertexType>;
//
//    std::vector<PrtVertexType> vec_vertex;
//    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0,0)));
//    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0.5,0)));
//    vec_vertex.push_back(std::make_shared<VertexType>(PointType(1.0,0)));
//    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0,0.5)));
//    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0.5,0.5)));
//    vec_vertex.push_back(std::make_shared<VertexType>(PointType(1,0.5)));
//    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0,1)));
//    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0.5,1)));
//    vec_vertex.push_back(std::make_shared<VertexType>(PointType(1,1)));
//
//    std::cout << "Vertex Num: " << " " <<  vec_vertex.size() << std::endl;
//
//    using CellDefineType = art_pde::Cell<PointType>;
//    using PrtCellType = std::shared_ptr<CellDefineType>;
//    std::vector<PrtCellType> vec_cell;
//
//    art_pde::CellFactory<PointType , art_pde::Dim2D> cell_builder;
//    // Cell 1
//    cell_builder.clearVertex();
//    cell_builder.addVertex(vec_vertex[0]);
//    cell_builder.addVertex(vec_vertex[1]);
//    cell_builder.addVertex(vec_vertex[4]);
//    cell_builder.addVertex(vec_vertex[3]);
//    vec_cell.push_back(cell_builder.create());
//
//    // Cell 2
//    cell_builder.clearVertex();
//    cell_builder.addVertex(vec_vertex[1]);
//    cell_builder.addVertex(vec_vertex[2]);
//    cell_builder.addVertex(vec_vertex[5]);
//    cell_builder.addVertex(vec_vertex[4]);
//    vec_cell.push_back(cell_builder.create());
//
//    // Cell 3
//    cell_builder.clearVertex();
//    cell_builder.addVertex(vec_vertex[3]);
//    cell_builder.addVertex(vec_vertex[4]);
//    cell_builder.addVertex(vec_vertex[7]);
//    cell_builder.addVertex(vec_vertex[6]);
//    vec_cell.push_back(cell_builder.create());
//
//    // Cell 4
//    cell_builder.clearVertex();
//    cell_builder.addVertex(vec_vertex[4]);
//    cell_builder.addVertex(vec_vertex[5]);
//    cell_builder.addVertex(vec_vertex[8]);
//    cell_builder.addVertex(vec_vertex[7]);
//    vec_cell.push_back(cell_builder.create());
//
//    for (int i = 0; i < vec_cell.size(); ++i) {
//        std::cout << "Cell : " << i << ", Type = " << vec_cell[i]->getCellTypeInString() << std::endl;
//        std::cout << *vec_cell[i] << std::endl;
//    }










    return 0;
}