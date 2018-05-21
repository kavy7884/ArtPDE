//
// Created by Chingkai Chou on 5/20/18.
//

#ifndef ARTPDE_MARVIN_GEOTEST3D_HPP
#define ARTPDE_MARVIN_GEOTEST3D_HPP

#include "point3.hpp"
#include "geo_data.hpp"
#include "geo_data_factory.hpp"

class GeoTest2D{
public:
    using VertexType = Vertex<Point3>;
    using VecPtrVertex = std::vector<std::shared_ptr<VertexType>>;
    using EdgeType = Edge<Point3>;
    using VecPtrEdge = std::vector<std::shared_ptr<EdgeType>>;
    using FaceType = Face<Point3>;
    using VecPtrFace = std::vector<std::shared_ptr<FaceType>>;
    using CellType = Cell<Point3>;
    using VecPtrCell = std::vector<std::shared_ptr<CellType>>;

    GeoTest2D() {
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(2.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(2.0, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 2.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 2.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 1.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 1.0));
    }




private:
    VecPtrVertex vec_ptr_vertex;
    VecPtrVertex vec_ptr_cell;

    VecPtrEdge vec_ptr_merged_edge;
    VecPtrFace vec_ptr_merged_face;
};

#endif //ARTPDE_MARVIN_GEOTEST3D_HPP
