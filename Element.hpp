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
        geoElementType = std::make_shared<ElementGeoType>(ElementGeoType::None);
        volume = std::make_shared<double>(0.0);
    }

    const VertexConnectType &getVertexConnect() const {
        return vertexConnect;
    }

    size_t getElementId(){ return vertexConnect->getId();}


    void calGeoElementType();
    void calGeoElementVolume();

    const std::shared_ptr<ElementGeoType> &getGeoElementType() const;

    const std::shared_ptr<double> &getVolume() const;

protected:
    GeoElement(){}
    VertexConnectType vertexConnect{nullptr};
    std::shared_ptr<ElementGeoType> geoElementType{nullptr};
    std::shared_ptr<double> volume{nullptr};
private:
    double calTetraVolume(size_t a, size_t b, size_t c, size_t d);
};

template <class Dimension>
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

    *this->volume = 0;
    
    for (size_t i = 0; i < vertexSize-1; ++i) {
        *this->volume += (vertex[i]->x()*vertex[i+1]->y() - vertex[i]->y()*vertex[i+1]->x());
    }

    *this->volume += (vertex[vertexSize-1]->x()*vertex[0]->y() - vertex[vertexSize-1]->y()*vertex[0]->x());

    *this->volume *= 0.5;
    
    //std::cout << "volume = " << volume << std::endl;
}

template <>
void GeoElement<Dim3D>::calGeoElementVolume(){
    const size_t& vertexSize = vertexConnect->getVertexSize();
    
    *this->volume = 0;
    
    if(vertexSize == 4){
        //Tetrahedron;
        *this->volume += calTetraVolume(0, 1, 2, 3);
    }
    else if(vertexSize == 5){
        //Pyramid;
        *this->volume += calTetraVolume(0, 1, 3, 4);
        *this->volume += calTetraVolume(2, 3, 1, 4);
    }
    else if(vertexSize == 6){
        //Prism;
        *this->volume += calTetraVolume(3, 4, 5, 0);
        *this->volume += calTetraVolume(0, 2, 5, 4);
        *this->volume += calTetraVolume(0, 1, 2, 4);
    }
    else if(vertexSize == 8){
        //Hexahedron;
        *this->volume += calTetraVolume(0, 1, 5, 3);
        *this->volume += calTetraVolume(0, 5, 4, 3);
        *this->volume += calTetraVolume(4, 5, 7, 3);
        *this->volume += calTetraVolume(1, 2, 6, 3);
        *this->volume += calTetraVolume(6, 5, 1, 3);
        *this->volume += calTetraVolume(3, 5, 6, 7);
    }
    else{
        *this->geoElementType = ElementGeoType::None;
        
    }
    
}

template<class Dimension>
double GeoElement<Dimension>::calTetraVolume(size_t a, size_t b, size_t c, size_t d){
    auto& vertex = vertexConnect->getVertex();
    
    Eigen::Vector3d v_ab(vertex[b]->x()-vertex[a]->x(), vertex[b]->y()-vertex[a]->y(), vertex[b]->z()-vertex[a]->z());
    Eigen::Vector3d v_ac(vertex[c]->x()-vertex[a]->x(), vertex[c]->y()-vertex[a]->y(), vertex[c]->z()-vertex[a]->z());
    Eigen::Vector3d v_ad(vertex[d]->x()-vertex[a]->x(), vertex[d]->y()-vertex[a]->y(), vertex[d]->z()-vertex[a]->z());

    return std::abs(v_ab.cross(v_ac).dot(v_ad)/ 6.0);
}

#endif //ARTCFD_ELEMENT_HPP
