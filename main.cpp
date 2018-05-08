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

    std::cout << geo.getTotalVertexNum() << std::endl;
    std::cout << geo.read(proj_1) << std::endl;
    std::cout << geo.getTotalVertexNum() << std::endl;



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
//    art_pde::CellBuilder<PointType , art_pde::Dim2D> cell_builder;
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