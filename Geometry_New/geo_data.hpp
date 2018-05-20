//
// Created by Chingkai Chou on 5/17/18.
//

#ifndef ARTPDE_GEO_DATA_HPP
#define ARTPDE_GEO_DATA_HPP

#include <ostream>
#include <set>
#include "geo_tree.hpp"


template <typename Data> class Face;
template <typename Data> class Edge;
template <typename Data> class Vertex;


template <typename Data>
class Face:
        public GeoTree,
        public GeoTree_Child<Edge<Data>>{
public:
    Face() :
            GeoTree(TreeType::TreeHead),
            GeoTree_Child<Edge<Data>>() {}

};


template <typename Data>
class Edge:
        public GeoTree,
        public GeoTree_Parent<Face<Data>>,
        public GeoTree_Child<Vertex<Data>>{
public:
    using PtrVertexType = std::shared_ptr<Vertex<Data>>;
    Edge() :
            GeoTree(TreeType::TreeConnect),
            GeoTree_Parent<Face<Data>>(),
            GeoTree_Child<Vertex<Data>>() {}

    bool operator==(Edge &rhs){
        if(this->c_getPtr_list_ptr_childs()->size() != rhs.c_getPtr_list_ptr_childs()->size()) return false;
        else{
            std::set<PtrVertexType> set_self_vertex(this->c_getPtr_list_ptr_childs()->cbegin(), this->c_getPtr_list_ptr_childs()->cend());
            std::set<PtrVertexType> set_rhs_vertex(rhs.c_getPtr_list_ptr_childs()->cbegin(), rhs.c_getPtr_list_ptr_childs()->cend());
            if ( std::equal(set_self_vertex.cbegin(), set_self_vertex.cend(), set_rhs_vertex.cbegin() ))return true;
            else return false;
        }
    }

    bool operator!=(const Edge &rhs) {
        return !(rhs == *this);
    }

};


template <typename Data>
class Vertex :
        public GeoTree,
        public GeoTree_Parent<Edge<Data>>{
public:
    using PtrData = std::shared_ptr<Data>;
    Vertex() :
            GeoTree(TreeType::TreeEnd),
            GeoTree_Parent<Edge<Data>>(){}

    Vertex(double x, double y, double z) :
            GeoTree(TreeType::TreeEnd),
            GeoTree_Parent<Edge<Data>>(){
        ptr_data = std::make_shared<Point3>(x, y, z);
    }

    const PtrData &getPtr_data() const {
        return ptr_data;
    }

    void setPtr_data(const PtrData &ptr_data) {
        Vertex::ptr_data = ptr_data;
    }

    friend std::ostream &operator<<(std::ostream &os, const Vertex<Data> &vertex) {
        os << *vertex.getPtr_data() ;
        return os;
    }

private:
    PtrData ptr_data{nullptr};
};


#endif //ARTPDE_GEO_DATA_HPP
