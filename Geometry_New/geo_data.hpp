//
// Created by Chingkai Chou on 5/17/18.
//

#ifndef ARTPDE_GEO_DATA_HPP
#define ARTPDE_GEO_DATA_HPP

#include <ostream>
#include <set>
#include "geo_tree.hpp"
#include "point3.hpp"

template <typename Data> class Cell;
template <typename Data> class Face;
template <typename Data> class Edge;
template <typename Data> class Vertex;


template <typename Data>
class Cell:
        public GeoTree<Cell<Data>>,
        public GeoTree_Child<Face<Data>>,
        public std::enable_shared_from_this<Cell<Data>>{
public:
    Cell() :
            GeoTree<Cell<Data>>(TreeType::TreeHead),
            GeoTree_Child<Face<Data>>() {}

    void link_self() override {
        this->linked_to = this->shared_from_this();
    }
};

template <typename Data>
class Face:
        public GeoTree<Face<Data>>,
        public GeoTree_Parent<Cell<Data>>,
        public GeoTree_Child<Edge<Data>>,
        public std::enable_shared_from_this<Face<Data>>{
public:
    using PtrVertexType = std::shared_ptr<Vertex<Data>>;
    using VecPtrVertexType = std::vector<PtrVertexType>;
    using PtrEdgeType = std::shared_ptr<Edge<Data>>;
    using VecPtrEdgeType = std::vector<PtrEdgeType>;

    Face() :
            GeoTree<Face<Data>>(TreeType::TreeConnect),
            GeoTree_Parent<Cell<Data>>(),
            GeoTree_Child<Edge<Data>>() {}

    void link_self() override {
        this->linked_to = this->shared_from_this();
    }

    const VecPtrEdgeType& getConnected_Edge() const{
        return this->c_getVec_ptr_childs();
    }

    bool operator==(Face &rhs){
        if(this->getConnected_Edge().size() != rhs.getConnected_Edge().size())  return false;
        else{
            std::set<PtrEdgeType> set_self, set_rhs;
            for (auto &ptr_edge: this->getConnected_Edge()) {
                set_self.insert(ptr_edge->getLinked_to());
            }
            for (auto &ptr_edge: rhs.getConnected_Edge()) {
                set_rhs.insert(ptr_edge->getLinked_to());
            }
            if ( std::equal(set_self.cbegin(), set_self.cend(), set_rhs.cbegin() ))return true;
            else return false;
        }
        return true;
    }

    bool operator!=(const Face &rhs) {
        return !(rhs == *this);
    }

};

template <typename Data>
class Edge:
        public GeoTree<Edge<Data>>,
        public GeoTree_Parent<Face<Data>>,
        public GeoTree_Child<Vertex<Data>>,
        public std::enable_shared_from_this<Edge<Data>>{
public:
    using PtrVertexType = std::shared_ptr<Vertex<Data>>;
    using VecPtrVertexType = std::vector<PtrVertexType>;
    using PtrFaceType = std::shared_ptr<Face<Data>>;
    using ListPtrFaceType = std::list<PtrFaceType>;
    Edge() :
            GeoTree<Edge<Data>>(TreeType::TreeConnect),
            GeoTree_Parent<Face<Data>>(),
            GeoTree_Child<Vertex<Data>>() {
    }

    void link_self() override {
        this->linked_to = this->shared_from_this();
    }

    const VecPtrVertexType& getConnected_Vertex() const{
        return this->c_getVec_ptr_childs();
    }

    const ListPtrFaceType& getConnected_Face() {
        return this->getLinked_to()->c_getList_ptr_parents();
    }

    bool operator==(Edge &rhs){
        if(this->getConnected_Vertex().size() != rhs.getConnected_Vertex().size())  return false;
        else{
            std::set<PtrVertexType> set_self(this->getConnected_Vertex().cbegin(), this->getConnected_Vertex().cend());
            std::set<PtrVertexType> set_rhs(rhs.getConnected_Vertex().cbegin(), rhs.getConnected_Vertex().cend());
            if ( std::equal(set_self.cbegin(), set_self.cend(), set_rhs.cbegin() ))return true;
            else return false;
        }
    }

    bool operator!=(const Edge &rhs) {
        return !(rhs == *this);
    }
//
//    size_t getNum_ConnectedFace(){
//        size_t re_num = 0;
//        for (auto &ptr_face : this->getVec_ptr_parents()) {
//            if(!ptr_face->isLinked()){
//                ++re_num;
//            }
//        }
//        return re_num;
//    }

};


template <typename Data>
class Vertex :
        public GeoTree<Vertex<Data>>,
        public GeoTree_Parent<Edge<Data>>,
        public std::enable_shared_from_this<Vertex<Data>>{
public:
    using PtrDataType = std::shared_ptr<Data>;
    using PtrEdgeType = std::shared_ptr<Edge<Data>>;
    using ListPtrEdgeType = std::list<PtrEdgeType>;
    Vertex() :
            GeoTree<Vertex<Data>>(TreeType::TreeEnd),
            GeoTree_Parent<Edge<Data>>(){
    }

    Vertex(double x, double y, double z) :
            GeoTree<Vertex<Data>>(TreeType::TreeEnd),
            GeoTree_Parent<Edge<Data>>(){
        ptr_data = std::make_shared<Point3>(x, y, z);
    }

    void link_self() override {
        this->linked_to = this->shared_from_this();
    }

    const ListPtrEdgeType& getConnected_Edge(){
        return this->getLinked_to()->c_getList_ptr_parents();
    }

    const PtrDataType &getPtr_data() const {
        return ptr_data;
    }

    void setPtr_data(const PtrDataType &ptr_data) {
        Vertex::ptr_data = ptr_data;
    }

    friend std::ostream &operator<<(std::ostream &os, const Vertex<Data> &vertex) {
        os << *vertex.getPtr_data() ;
        return os;
    }

//    size_t getNum_ConnectedEdge(){
//        size_t re_num = 0;
//        for (auto &ptr_edge : this->getVec_ptr_parents()) {
//            if(!ptr_edge->isLinked()){
//                ++re_num;
//            }
//        }
//        return re_num;
//    }

private:
    PtrDataType ptr_data{nullptr};
};


#endif //ARTPDE_GEO_DATA_HPP
