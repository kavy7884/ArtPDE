//
// Created by Chingkai Chou on 5/2/18.
//
#include <iostream>
#include "project_uility.hpp"
#include "geometry.hpp"

int main() {


    auto proj_1 = ArtProject::create("TestProj").setRunPath(".").setDivideSlash("/").setInitialFolderName("Init").build();

    using GeoDataType =  art_pde::GeometryData<art_pde::MeshTypeMethod, art_pde::Dim2D, art_pde::CartesianCoordinate>;

    art_pde::Geometry<GeoDataType, art_pde::GeometryDataReaderArtPDE> geo;

    std::cout << "Vertex Num (before): " << geo.getTotal_VertexNum() << std::endl;
    std::cout << "Cell Num (before): " << geo.getTotal_CellNum() << std::endl;
    std::cout << geo.read(proj_1) << std::endl;
    std::cout << "Vertex Num (after): " << geo.getTotal_VertexNum() << std::endl;
    std::cout << "Cell Num (after): " << geo.getTotal_CellNum() << std::endl;

    auto all_ptr_point_on_vertex = geo.getTotal_PtrPointOnVertex();

    std::cout << "List all vertex points: " << std::endl;
    for (size_t i = 0; i < geo.getTotal_VertexNum(); ++i) {
        std::cout << *all_ptr_point_on_vertex[i] << std::endl;
    }

    std::cout << "List all vertex points in each cell: " << std::endl;
    for (size_t i = 0; i < geo.getTotal_CellNum(); ++i) {
        auto vec_ptr_point_on_vertex_in_cell_id = geo.getCell_PtrPointOnVertex(i);
        std::cout << "Cell Type: " << GeoDataType::Type::GeoCellType::convertCellTypeInString(
                geo.getCell_CellType(i) ) << "\t";
        for(auto &ptr_pt : vec_ptr_point_on_vertex_in_cell_id){
            std::cout << *ptr_pt << "\t";
        }
        std::cout << std::endl;
    }



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
//    using CellType = art_pde::Cell<PointType>;
//    using PrtCellType = std::shared_ptr<CellType>;
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