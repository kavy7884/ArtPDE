//
// Created by Chingkai Chou on 5/20/18.
//

#ifndef ARTPDE_MARVIN_GEOTEST3D_HPP
#define ARTPDE_MARVIN_GEOTEST3D_HPP

#include "point3.hpp"
#include "geo_data.hpp"
#include "geo_data_factory.hpp"

class GeoTest3D{
public:
    using VertexType = Vertex<Point3>;
    using VecPtrVertex = std::vector<std::shared_ptr<VertexType>>;
    using EdgeType = Edge<Point3>;
    using VecPtrEdge = std::vector<std::shared_ptr<EdgeType>>;
    using FaceType = Face<Point3>;
    using VecPtrFace = std::vector<std::shared_ptr<FaceType>>;
    using CellType = Cell<Point3>;
    using VecPtrCell = std::vector<std::shared_ptr<CellType>>;

    GeoTest3D() {
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 0.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(2.0, 0.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 1.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 1.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(2.0, 1.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 2.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 2.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 1.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 1.0));

        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 1.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 1.0, 1.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 1.0, 1.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 0.0, 1.0));

        vec_ptr_cell.push_back(std::make_shared<HexaCell<Point3>>(
                vec_ptr_vertex[0],
                vec_ptr_vertex[1],
                vec_ptr_vertex[2],
                vec_ptr_vertex[3],
                vec_ptr_vertex[4],
                vec_ptr_vertex[5],
                vec_ptr_vertex[6],
                vec_ptr_vertex[7]
        ));


    }




private:
    VecPtrVertex vec_ptr_vertex;
    VecPtrCell vec_ptr_cell;

    VecPtrEdge vec_ptr_merged_edge;
    VecPtrFace vec_ptr_merged_face;
};

#endif //ARTPDE_MARVIN_GEOTEST3D_HPP
