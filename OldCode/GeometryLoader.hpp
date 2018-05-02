//
// Created by Chingkai Chou on 3/21/18.
//

#ifndef ARTCFD_GEOMETRYLOADER_H
#define ARTCFD_GEOMETRYLOADER_H

#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include "DimensionUtility.hpp"
#include "ProjectUility.hpp"
#include "Geometry.hpp"

enum class GeometrySourceFormat{ File, HDF5};

template <class Dimension, class NumericalMethodUtility>
class Geometry;

template <class Dimension, class NumericalMethodUtility>
class GeometryLoader{
public:

    GeometryLoader(const std::string &geoName);

    virtual void setProj(const std::shared_ptr<ArtProject> &proj) {
        GeometryLoader::proj = proj;
    }

    virtual void setGeoData(const std::shared_ptr<typename Geometry<Dimension, NumericalMethodUtility>::GeoType> & geoData);

    virtual bool load() = 0;

protected:
    std::string geoName;
    std::shared_ptr<ArtProject> proj{nullptr};
    std::shared_ptr<typename Geometry<Dimension, NumericalMethodUtility>::GeoType> geoData{nullptr};

};

template <class Dimension, class NumericalMethodUtility>
class GeometryFileLoaderBase : public GeometryLoader<Dimension, NumericalMethodUtility>{
public:
    GeometryFileLoaderBase(const std::string &geoName) : GeometryLoader<Dimension, NumericalMethodUtility>(geoName) {}

    bool load() override {
        return true;
    }

protected:
    bool load_BasicMeshData(){
        load_xNode();
        load_cElement30();
        return true;
    }

private:
    bool load_xNode();
    bool load_cElement30();
};


template <class Dimension, class NumericalMethodUtility>
class GeometryFileLoader;


template <class Dimension>
class GeometryFileLoader<Dimension, MeshTypeMethod> : public GeometryFileLoaderBase<Dimension, MeshTypeMethod>{
public:
    GeometryFileLoader(const std::string &geoName) : GeometryFileLoaderBase<Dimension, MeshTypeMethod>(geoName) {}

    bool load() override {
        this->load_BasicMeshData();
        return GeometryFileLoaderBase<Dimension, MeshTypeMethod>::load();
    }

};

template <class Dimension>
class GeometryFileLoader<Dimension, FEM> : public GeometryFileLoaderBase<Dimension, FEM>{
public:
    GeometryFileLoader(const std::string &geoName) : GeometryFileLoaderBase<Dimension, FEM>(geoName) {}

    bool load() override {
        this->load_BasicMeshData();
        return GeometryFileLoaderBase<Dimension, FEM>::load();
    }

};

template<class Dimension, class NumericalMethodUtility>
void GeometryLoader<Dimension, NumericalMethodUtility>::setGeoData(
        const std::shared_ptr<typename Geometry<Dimension, NumericalMethodUtility>::GeoType>& geoData){
    GeometryLoader::geoData = geoData;
}

template<class Dimension, class NumericalMethodUtility>
GeometryLoader<Dimension, NumericalMethodUtility>::GeometryLoader(const std::string &geoName):geoName(geoName) {};

template <class Dimension, class NumericalMethodUtility>
bool GeometryFileLoaderBase<Dimension, NumericalMethodUtility>::load_xNode(){
    std::string path_xNode = this->proj->getProjectGeometryPath() + this->geoName + this->proj->getSlash()
                            + "xNode.txt";

    std::ifstream fs;
    std::string bufferLine;
    fs.open(path_xNode);
    if(fs.is_open()){
        while(getline( fs, bufferLine )) {
            std::stringstream w(bufferLine);
            auto p_point = std::make_shared<Point<Dimension>>(this->geoData->xNode);
            for (size_t i = 0; i < Dimension::Dim; ++i) {
                w >> p_point->getDataByDim(i);
            }
            this->geoData->xNode->addPoint(p_point);
        }
    }else
        return false;

    std::cout << "Geometry loaded >> " << path_xNode << std::endl;
    return true;
}

template<class Dimension, class NumericalMethodUtility>
bool GeometryFileLoaderBase<Dimension, NumericalMethodUtility>::load_cElement30() {
    std::string path_cElement30 = this->proj->getProjectGeometryPath() + this->geoName + this->proj->getSlash()
                             + "cElement30.txt";

    std::ifstream fs;
    std::string bufferLine;
    size_t tempId;

    if(this->geoData->xNode == nullptr)
        return false;

    fs.open(path_cElement30);
    if(fs.is_open()) {
        while(getline( fs, bufferLine )) {
            std::stringstream w(bufferLine);
            auto p_elementConnect = std::make_shared<ElementConnectivity<Dimension>>();
            while(w >> tempId){
                auto p_point = this->geoData->xNode->getPoint(tempId);
                p_elementConnect->addVertex(p_point);
            }
            this->geoData->cElement30->addElementConnectivity(p_elementConnect);
        }
    }else
        return false;

    std::cout << "Geometry loaded >> " << path_cElement30 << std::endl;
    return true;
}


#endif //ARTCFD_GEOMETRYLOADER_H
