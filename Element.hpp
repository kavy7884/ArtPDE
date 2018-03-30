//
// Created by Chingkai Chou on 3/29/18.
//

#ifndef ARTCFD_ELEMENT_HPP
#define ARTCFD_ELEMENT_HPP

#include <memory>
#include "Point.hpp"
#include "Connectivity.hpp"
#include <Eigen/Dense>

enum class ElementGeoType{None, Line, Triangle, Quadrilateral, Tetrahedron, Pyramid, Prism, Hexahedron};

class Element{
public:
    Element() {}
};

template <class Dimension>
class GeoElement : virtual public Element{
public:
    using VertexConnectType = typename ElementList<Dimension>::PtrElementConnectivityType;
    GeoElement(const VertexConnectType &vertexConnect) : vertexConnect(vertexConnect) {
        calGeoElementType();
    }

    const VertexConnectType &getVertexConnect() const {
        return vertexConnect;
    }

    size_t &getElementId(){ return vertexConnect->getId();}

    ElementGeoType getGeoElementType() const;

    void calGeoElementType();
    void calGeoElementVolume();

protected:
    ElementGeoType geoElementType{ElementGeoType::None};
    VertexConnectType vertexConnect;
    double volume{0.0};
};

template <class Dimension>
void GeoElement<Dimension>::calGeoElementType(){ geoElementType = ElementGeoType::None;}

template<class Dimension>
ElementGeoType GeoElement<Dimension>::getGeoElementType() const {
    return geoElementType;
}

template <>
void GeoElement<Dim2D>::calGeoElementType(){
    if(vertexConnect->getVertexSize() == 3){ geoElementType = ElementGeoType::Triangle;}
    else if(vertexConnect->getVertexSize() == 4){ geoElementType = ElementGeoType::Quadrilateral;}
    else{geoElementType = ElementGeoType::None;}
}

template <>
void GeoElement<Dim3D>::calGeoElementType(){
    if(vertexConnect->getVertexSize() == 4){ geoElementType = ElementGeoType::Tetrahedron;}
    else if(vertexConnect->getVertexSize() == 5){ geoElementType = ElementGeoType::Pyramid;}
    else if(vertexConnect->getVertexSize() == 6){ geoElementType = ElementGeoType::Prism;}
    else if(vertexConnect->getVertexSize() == 8){ geoElementType = ElementGeoType::Hexahedron;}
    else{geoElementType = ElementGeoType::None;}
}

template <>
void GeoElement<Dim2D>::calGeoElementVolume(){
    const size_t& vertexSize = vertexConnect->getVertexSize();
    auto& vertex = vertexConnect->getVertex();
    
    volume = 0;
    
    for (size_t i = 0; i < vertexSize-1; ++i) {
        volume += (vertex[i]->x()*vertex[i+1]->y() - vertex[i]->y()*vertex[i+1]->x());
    }

    volume += (vertex[vertexSize-1]->x()*vertex[0]->y() - vertex[vertexSize-1]->y()*vertex[0]->x());
    
    volume *= 0.5;
    
    std::cout << "volume = " << volume << std::endl;
}

template <>
void GeoElement<Dim3D>::calGeoElementVolume(){
    
    
}


#endif //ARTCFD_ELEMENT_HPP
