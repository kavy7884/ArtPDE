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
class Cell:
        public GeoTree,
        public GeoTree_Child<Face<Data>>{
public:
    Cell() :
            GeoTree(TreeType::TreeHead),
            GeoTree_Child<Face<Data>>() {}

};

template <typename Data>
class Face:
        public GeoTree,
        public GeoTree_Parent<Cell<Data>>,
        public GeoTree_Child<Edge<Data>>,
        public std::enable_shared_from_this<Face<Data>>{
public:
    using PtrVertexType = std::shared_ptr<Vertex<Data>>;
    using VecPtrVertexType = std::vector<PtrVertexType>;
    using PtrEdgeType = std::shared_ptr<Edge<Data>>;

    Face() :
            GeoTree(TreeType::TreeConnect),
            GeoTree_Parent<Cell<Data>>(),
            GeoTree_Child<Edge<Data>>() {}

    bool operator==(Face &rhs){
        if(this->c_getPtr_list_ptr_childs()->size() != rhs.c_getPtr_list_ptr_childs()->size()) return false;
        else{
            std::set<PtrEdgeType> set_self(this->c_getPtr_list_ptr_childs()->cbegin(), this->c_getPtr_list_ptr_childs()->cend());
            std::set<PtrEdgeType> set_rhs(rhs.c_getPtr_list_ptr_childs()->cbegin(), rhs.c_getPtr_list_ptr_childs()->cend());
            if ( std::equal(set_self.cbegin(), set_self.cend(), set_rhs.cbegin() ))return true;
            else return false;
        }
    }

    bool operator!=(const Face &rhs) {
        return !(rhs == *this);
    }

    VecPtrVertexType getVertex(){
        VecPtrVertexType re_ptr_vertex;
        std::shared_ptr<Face<Data>> ptr_this_face = this->shared_from_this();
        auto it_edge = this->getPtr_list_ptr_childs()->begin();
        while(it_edge != this->getPtr_list_ptr_childs()->end()){
            if((*it_edge)->isMerged())
            {
                auto iter_vertex = (*it_edge)->getMergedGroup(ptr_this_face);
                re_ptr_vertex.push_back(*iter_vertex);
            }else{
                re_ptr_vertex.push_back(*(*it_edge)->c_getPtr_list_ptr_childs()->cbegin());
            }

            ++it_edge;
        }
        return re_ptr_vertex;
    }

};

template <typename Data>
class Edge:
        public GeoTree,
        public GeoTree_Parent<Face<Data>>,
        public GeoTree_Child<Vertex<Data>>{
public:
    using PtrVertexType = std::shared_ptr<Vertex<Data>>;
    using VecPtrVertexType = std::vector<PtrVertexType>;
    using PtrFaceType = std::shared_ptr<Face<Data>>;
    Edge() :
            GeoTree(TreeType::TreeConnect),
            GeoTree_Parent<Face<Data>>(),
            GeoTree_Child<Vertex<Data>>() {
        this->setNum_childs_per_group(2);
    }

    bool operator==(Edge &rhs){
        if(this->c_getPtr_list_ptr_childs()->size() != rhs.c_getPtr_list_ptr_childs()->size()) return false;
        else{
            std::set<PtrVertexType> set_self(this->c_getPtr_list_ptr_childs()->cbegin(), this->c_getPtr_list_ptr_childs()->cend());
            std::set<PtrVertexType> set_rhs(rhs.c_getPtr_list_ptr_childs()->cbegin(), rhs.c_getPtr_list_ptr_childs()->cend());
            if ( std::equal(set_self.cbegin(), set_self.cend(), set_rhs.cbegin() ))return true;
            else return false;
        }
    }

    bool operator!=(const Edge &rhs) {
        return !(rhs == *this);
    }

    typename
    std::list<PtrVertexType>::const_iterator getMergedGroup(const PtrFaceType &ptr_face){
        typename std::list<PtrVertexType>::const_iterator iter_child = this->c_getPtr_list_ptr_childs()->cbegin();

        size_t counter = 0;
        auto it_face = this->c_getPtr_list_ptr_parents()->cbegin();
        while(it_face != this->c_getPtr_list_ptr_parents()->cend()){
            if(*it_face == ptr_face){
                //std::cout << "Get" << std::endl;
                break;
            }
            ++it_face;
            ++counter;
        }

        counter *= this->getNum_childs_per_group();

        for (size_t i = 0; i < counter; ++i) {
            ++iter_child;
        }

        return iter_child;
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
