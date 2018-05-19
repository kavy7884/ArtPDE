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
public:
    using VecPtrVertex = std::vector<std::shared_ptr<Vertex<Point3>>>;
    using VecPtrFace = std::vector<std::shared_ptr<Face<Point3>>>;
    using VecPtrEdge = std::vector<std::shared_ptr<Edge<Point3>>>;

    GeoTest2D() {

        vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.5, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(1.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.0, 0.5, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.5, 0.5, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(1.0, 0.5, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.0, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.5, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(1.0, 1.0, 0.0));

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

        std::cout << "List Vertex: " << std::endl;
        for (auto &ptr_vertex: vec_ptr_vertex) {
            std::cout << *ptr_vertex << std::endl;
        }


        std::cout << "List Face's Vertex: " << std::endl;

        for (size_t i = 0;  i < vec_ptr_face.size(); ++i) {
            std::cout << "Face " << i << " is : ";
            for (auto &ptr_edge: vec_ptr_face[i]->getChild()){
                std::cout << *(*ptr_edge->getChild().begin()) ;
            }
            std::cout << std::endl;
        }
    }


    void MergeEdge(){
        for (size_t i = 0; i < vec_ptr_vertex.size(); ++i) {
            auto & list_ptr_parent_edge = this->vec_ptr_vertex[i]->getParent();

            auto it_master = list_ptr_parent_edge.begin();

            while( it_master != list_ptr_parent_edge.end()){
                if((*it_master)->isMerged()){
                    ++it_master;
                    continue;
                }
                auto it_slave = it_master;
                ++it_slave;
                while( it_slave != list_ptr_parent_edge.end()) {
                    if((*(*it_master)) == (*(*it_slave))){
                        std::cout << "Merge" << std::endl;

                        (*it_master)->merge(*(*it_slave));

                        for(auto & replace_parent : (*it_slave)->getParent()){
                            replace_parent->replaceChild((*it_slave), (*it_master));
                        }

                        for(auto & replace_child : (*it_slave)->getChild()){
                            replace_child->eraseParent((*it_slave));
                        }
                        it_slave = it_master;
                    }
                    ++it_slave;
                }
                (*it_master)->setMerged(true);
                ++it_master;
            }
        }
    };


    VecPtrVertex vec_ptr_vertex;
    VecPtrFace vec_ptr_face;

    VecPtrEdge vec_ptr_merged_edge;

};

#endif //ARTPDE_GEOTEST_HPP
