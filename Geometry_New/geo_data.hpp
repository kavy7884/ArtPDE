//
// Created by Chingkai Chou on 5/17/18.
//

#ifndef ARTPDE_GEO_DATA_HPP
#define ARTPDE_GEO_DATA_HPP

#include <ostream>
#include <set>
#include "geo_tree.hpp"

template <typename Data> class Cell;
template <typename Data> class Edge;
template <typename Data> class Vertex;

template <typename Data>
class Face: public GeoHead<Edge<Data>>{
public:
    Face() : GeoHead<Edge<Data>>(){}

};

template <typename Data>
class Edge: public GeoConnect<Face<Data>, Vertex<Data>>{
public:
    using PtrVertexType = std::shared_ptr<Vertex<Data>>;
    Edge() : GeoConnect<Face<Data>, Vertex<Data>>(){}

    bool operator==(Edge &rhs){
        if(this->getChild().size() != rhs.getChild().size()) return false;
        else{
            std::set<PtrVertexType> set_self_vertex(this->getChild().cbegin(), this->getChild().cend());
            std::set<PtrVertexType> set_rhs_vertex(rhs.getChild().cbegin(), rhs.getChild().cend());
            if ( std::equal(set_self_vertex.cbegin(), set_self_vertex.cend(), set_rhs_vertex.cbegin() ))return true;
            else return false;
        }
    }

    bool operator!=(const Edge &rhs) {
        return !(rhs == *this);
    }

};

template <typename Data>
class Vertex: public GeoEnd<Edge<Data>>{
public:
    Vertex(double x, double y, double z) : GeoEnd<Edge<Data>>(), point(x, y, z){}

    friend std::ostream &operator<<(std::ostream &os, const Vertex<Data> &vertex) {
        os << vertex.getPoint() ;
        return os;
    }

    Data getPoint() const {
        return point;
    }


private:
    Data point;
};



#endif //ARTPDE_GEO_DATA_HPP
