//
// Created by Chingkai Chou on 5/17/18.
//

#ifndef ARTPDE_MARVIN_GEO_DATA_HPP
#define ARTPDE_MARVIN_GEO_DATA_HPP

#include <ostream>
#include "geo_tree.hpp"

template <typename Data> class Cell;
template <typename Data> class Edge;
template <typename Data> class Vertex;

template <typename Data>
class Cell: public GeoHead<Edge<Data>>{
public:
    Cell() : GeoHead<Edge<Data>>(){}

};

template <typename Data>
class Edge: public GeoConnect<Cell<Data>, Vertex<Data>>{
public:
    Edge() : GeoConnect<Cell<Data>, Vertex<Data>>(){}

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



#endif //ARTPDE_MARVIN_GEO_DATA_HPP
