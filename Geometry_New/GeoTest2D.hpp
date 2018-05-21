//
// Created by Chingkai Chou on 5/18/18.
//

#ifndef ARTPDE_GEOTEST_HPP
#define ARTPDE_GEOTEST_HPP

#include <iostream>
#include <memory>
#include <vector>

#include "point3.hpp"
#include "geo_data.hpp"
#include "geo_data_factory.hpp"

class GeoTest2D{
    using VertexType = Vertex<Point3>;
    using VecPtrVertex = std::vector<std::shared_ptr<VertexType>>;
    using EdgeType = Edge<Point3>;
    using VecPtrEdge = std::vector<std::shared_ptr<EdgeType>>;
    using FaceType = Face<Point3>;
    using VecPtrFace = std::vector<std::shared_ptr<FaceType>>;

public:
    GeoTest2D() {
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.5, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.5, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.5, 0.5, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 0.5, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.5, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 1.0, 0.0));

        vec_ptr_face.push_back(std::make_shared<QuadFace<Point3>>(
                vec_ptr_vertex[0], vec_ptr_vertex[1], vec_ptr_vertex[4], vec_ptr_vertex[3])
        );
        vec_ptr_face.push_back(std::make_shared<QuadFace<Point3>>(
                vec_ptr_vertex[1], vec_ptr_vertex[2], vec_ptr_vertex[5], vec_ptr_vertex[4])
        );
        vec_ptr_face.push_back(std::make_shared<QuadFace<Point3>>(
                vec_ptr_vertex[3], vec_ptr_vertex[4], vec_ptr_vertex[7], vec_ptr_vertex[6])
        );
        vec_ptr_face.push_back(std::make_shared<QuadFace<Point3>>(
                vec_ptr_vertex[4], vec_ptr_vertex[5], vec_ptr_vertex[8], vec_ptr_vertex[7])
        );

    }

    void mergeEdge(){
        GeoTree_LayerMerge<VertexType, EdgeType> merge_algo(this->vec_ptr_vertex);
        vec_ptr_merged_edge = merge_algo.merge();
    }

    void checkOut(){

        std::cout << "List Vertex: " << std::endl;
        for (auto &ptr_vertex: vec_ptr_vertex) {
            std::cout << *ptr_vertex << std::endl;
        }

        auto vec_ptr_vertex = vec_ptr_face[0]->getVertex();
        std::cout << "List Face's Vertex: " << std::endl;
        for (size_t i = 0;  i < vec_ptr_face.size(); ++i) {
            std::cout << "Face " << i << " is : ";
            auto vec_ptr_vertex = vec_ptr_face[i]->getVertex();

            for (auto &ptr_vertex: vec_ptr_vertex){
                std::cout << *ptr_vertex ;
            }
            std::cout << std::endl;
        }

        std::cout << "Merged Edge Num: " << std::endl;
        std::cout << vec_ptr_merged_edge.size() << std::endl;

    }


private:
    VecPtrVertex vec_ptr_vertex;
    VecPtrFace vec_ptr_face;

    VecPtrEdge vec_ptr_merged_edge;
};

#endif //ARTPDE_GEOTEST_HPP
