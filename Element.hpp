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
        geoElementType = std::make_shared<ElementGeoType>(ElementGeoType::None);
        volume = std::make_shared<double>(0.0);
    }

    const VertexConnectType &getVertexConnect() const {
        return vertexConnect;
    }

    size_t &getElementId(){ return vertexConnect->getId();}
    size_t getElementId(){ return vertexConnect->getId();}

    ElementGeoType getGeoElementType() const;

    void calGeoElementType();
    void calGeoElementVolume();

    const std::shared_ptr<ElementGeoType> &getGeoElementType() const;

    const std::shared_ptr<double> &getVolume() const;

protected:
    ElementGeoType geoElementType{ElementGeoType::None};
    VertexConnectType vertexConnect;
    double volume{0.0};
    GeoElement(){}
    VertexConnectType vertexConnect{nullptr};
    std::shared_ptr<ElementGeoType> geoElementType{nullptr};
    std::shared_ptr<double> volume{nullptr};
};

template <class Dimension>
void GeoElement<Dimension>::calGeoElementType(){ geoElementType = ElementGeoType::None;}
class FemElement : public GeoElement<Dimension>{
public:
    using VertexConnectType = typename ElementList<Dimension>::PtrElementConnectivityType;
    FemElement(const VertexConnectType &vertexConnect) : GeoElement<Dimension>(vertexConnect) {}
    FemElement(const GeoElement<Dimension> &geoElement){copyDataFromGeoElement(geoElement);}
    FemElement<Dimension>& operator=(const GeoElement<Dimension>& geoElement){
        copyDataFromGeoElement(geoElement); return *this;
    }

private:
    void copyDataFromGeoElement(const GeoElement<Dimension> &geoElement);
};

template<class Dimension>
void FemElement<Dimension>::copyDataFromGeoElement(const GeoElement<Dimension> &geoElement) {
    this->vertexConnect = geoElement.getVertexConnect();
    this->geoElementType = geoElement.getGeoElementType();
    this->volume = geoElement.getVolume();
}

template<class Dimension>
ElementGeoType GeoElement<Dimension>::getGeoElementType() const {
const std::shared_ptr<ElementGeoType> &GeoElement<Dimension>::getGeoElementType() const {
    return geoElementType;
}

template<class Dimension>
const std::shared_ptr<double> &GeoElement<Dimension>::getVolume() const {
    return volume;
}

template <class Dimension>
void GeoElement<Dimension>::calGeoElementType(){ *this->geoElementType = ElementGeoType::None;}

template <>
void GeoElement<Dim2D>::calGeoElementType(){
    if(vertexConnect->getVertexSize() == 3){ geoElementType = ElementGeoType::Triangle;}
    else if(vertexConnect->getVertexSize() == 4){ geoElementType = ElementGeoType::Quadrilateral;}
    else{geoElementType = ElementGeoType::None;}
    if(vertexConnect->getVertexSize() == 3){ *this->geoElementType = ElementGeoType::Triangle;}
    else if(vertexConnect->getVertexSize() == 4){ *this->geoElementType = ElementGeoType::Quadrilateral;}
    else{*this->geoElementType = ElementGeoType::None;}
}

template <>
void GeoElement<Dim3D>::calGeoElementType(){
    if(vertexConnect->getVertexSize() == 4){ *this->geoElementType = ElementGeoType::Tetrahedron;}
    else if(vertexConnect->getVertexSize() == 5){ *this->geoElementType = ElementGeoType::Pyramid;}
    else if(vertexConnect->getVertexSize() == 6){ *this->geoElementType = ElementGeoType::Prism;}
    else if(vertexConnect->getVertexSize() == 8){ *this->geoElementType = ElementGeoType::Hexahedron;}
    else{*this->geoElementType = ElementGeoType::None;}
}

template <>
void GeoElement<Dim2D>::calGeoElementVolume(){
    const size_t& vertexSize = vertexConnect->getVertexSize();
    auto& vertex = vertexConnect->getVertex();
    
    volume = 0;

    *this->volume = 0;
    
    for (size_t i = 0; i < vertexSize-1; ++i) {
        volume += (vertex[i]->x()*vertex[i+1]->y() - vertex[i]->y()*vertex[i+1]->x());
        *this->volume += (vertex[i]->x()*vertex[i+1]->y() - vertex[i]->y()*vertex[i+1]->x());
    }

    volume += (vertex[vertexSize-1]->x()*vertex[0]->y() - vertex[vertexSize-1]->y()*vertex[0]->x());
    
    volume *= 0.5;
    *this->volume += (vertex[vertexSize-1]->x()*vertex[0]->y() - vertex[vertexSize-1]->y()*vertex[0]->x());

    *this->volume *= 0.5;
    
    std::cout << "volume = " << volume << std::endl;
    //std::cout << "volume = " << volume << std::endl;
}

<<<<<<< HEAD
//template <>
//void GeoElement<Dim3D>::calGeoElementVolume(){
//    if(vertexConnect->getVertexSize() == 4){
////        (axb.c)/6
//
//    }
//    else if(vertexConnect->getVertexSize() == 5){
////        2 tetra
//
//    }
//    else if(vertexConnect->getVertexSize() == 6){
////        3 tetra
//
//    }
//    else if(vertexConnect->getVertexSize() == 8){
////        6 tetra
//
//    }
//    else{
//        geoElementType = ElementGeoType::None;
//
//    }
//
//}
=======

template <>
void GeoElement<Dim3D>::calGeoElementVolume(){
    
    
}
>>>>>>> master


#endif //ARTCFD_ELEMENT_HPP
