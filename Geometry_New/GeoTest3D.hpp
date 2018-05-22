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

        // Single cube
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 1.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 0.0, 0.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 1.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 1.0, 1.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 1.0, 1.0));
        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 0.0, 1.0));

        auto cell = std::make_shared<HexaCell<Point3>>();
        cell->create(vec_ptr_vertex[0],
                     vec_ptr_vertex[1],
                     vec_ptr_vertex[2],
                     vec_ptr_vertex[3],
                     vec_ptr_vertex[4],
                     vec_ptr_vertex[5],
                     vec_ptr_vertex[6],
                     vec_ptr_vertex[7]);
        vec_ptr_cell.push_back(cell);

//        // Two cube
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 0.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(2.0, 0.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 1.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 1.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(2.0, 1.0, 0.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 0.0, 1.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 0.0, 1.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(2.0, 0.0, 1.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(0.0, 1.0, 1.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(1.0, 1.0, 1.0));
//        vec_ptr_vertex.push_back(std::make_shared<VertexType>(2.0, 1.0, 1.0));
//
//        auto cell = std::make_shared<HexaCell<Point3>>();
//        cell->create(vec_ptr_vertex[0],
//                     vec_ptr_vertex[1],
//                     vec_ptr_vertex[4],
//                     vec_ptr_vertex[3],
//                     vec_ptr_vertex[6],
//                     vec_ptr_vertex[7],
//                     vec_ptr_vertex[10],
//                     vec_ptr_vertex[9]);
//        vec_ptr_cell.push_back(cell);
//
//        cell = std::make_shared<HexaCell<Point3>>();
//        cell->create(vec_ptr_vertex[1],
//                     vec_ptr_vertex[2],
//                     vec_ptr_vertex[5],
//                     vec_ptr_vertex[4],
//                     vec_ptr_vertex[7],
//                     vec_ptr_vertex[8],
//                     vec_ptr_vertex[11],
//                     vec_ptr_vertex[10]);
//        vec_ptr_cell.push_back(cell);

    }

    void checkOut(){

        std::cout << vec_ptr_merged_edge.size() << std::endl;

        std::cout << "List Vertex: " << std::endl;
        for (auto &ptr_vertex: vec_ptr_vertex) {
            std::cout << *ptr_vertex << " , connected Edge num: " << ptr_vertex->getNum_ConnectedEdge() << std::endl;

        }

        for (auto &ptr_vertex: vec_ptr_vertex) {
            for(auto &ptr_edge : ptr_vertex->getVec_ptr_parents()){
                if(!ptr_edge->isLinked()){
                    vec_ptr_edge.push_back(ptr_edge);
                }
            }
        }

        std::cout << vec_ptr_edge.size() << std::endl;

        for(auto &ptr_edge: vec_ptr_edge){
            std::cout << "Edge connected Face num: " << ptr_edge->getNum_ConnectedFace() << std::endl;
            //std::cout <<ptr_edge->a << std::endl;
        }




//        std::cout << "List Vertex: " << std::endl;
//        for (auto &ptr_vertex: vec_ptr_vertex) {
//            std::cout << *ptr_vertex << std::endl;
//        }
//
//         //Check net memory
//        std::set<std::shared_ptr<EdgeType>> net_edge;
//        std::set<std::shared_ptr<FaceType>> net_face;
//
//        for (size_t i = 0; i < vec_ptr_cell.size(); ++i) {
//            auto it_face = vec_ptr_cell[i]->c_getPtr_list_ptr_childs()->begin();
//            while(it_face != vec_ptr_cell[i]->c_getPtr_list_ptr_childs()->end()){
//                net_face.insert(*it_face);
//                auto it_edge = (*it_face)->c_getPtr_list_ptr_childs()->begin();
//                while(it_edge != (*it_face)->c_getPtr_list_ptr_childs()->end()){
//                    net_edge.insert(*it_edge);
//
//                    ++it_edge;
//                }
//
//                ++it_face;
//            }
//
//        }
//        std::cout << net_face.size() << std::endl;
//        std::cout << net_edge.size() << std::endl;
//
////        std::cout << "List Vertex parent num: " << std::endl;
////        for (auto &ptr_vertex: vec_ptr_vertex) {
////            std::cout << *ptr_vertex << std::endl;
////            std::cout << ptr_vertex->getPtr_list_ptr_parents()->size()<< std::endl;
////        }
//
//        std::cout << "List Edge childs num: " << std::endl;
//        for (auto &ptr_edge: net_edge) {
////            std::cout << ptr_edge->getPtr_list_ptr_childs()->size()<< std::endl;
//            std::cout << ptr_edge->getPtr_list_ptr_parents()->size()<< std::endl;
//        }
//



    }

    void merge(){
        std::cout << "<Start> Merge Edge " << std::endl;
        mergeEdge();
        std::cout << "< End > Merge Edge " << std::endl;

        std::cout << "<Start> Merge Face " << std::endl;
        mergeFace();
        std::cout << "< End > Merge Face " << std::endl;
    }

protected:

    void mergeEdge(){
        GeoTree_LayerMerge<VertexType, EdgeType> merge_algo(this->vec_ptr_vertex);
        this->vec_ptr_merged_edge = merge_algo.merge();
    }

    void mergeFace(){
        GeoTree_LayerMerge<EdgeType, FaceType> merge_algo(this->vec_ptr_merged_edge);
        this->vec_ptr_merged_face = merge_algo.merge();
    }


private:
    VecPtrVertex vec_ptr_vertex;
    VecPtrEdge vec_ptr_edge;
    VecPtrCell vec_ptr_cell;

    VecPtrEdge vec_ptr_merged_edge;
    VecPtrFace vec_ptr_merged_face;
};

#endif //ARTPDE_MARVIN_GEOTEST3D_HPP
