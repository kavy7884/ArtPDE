//
// Created by Chingkai Chou on 2/23/18.
//

#ifndef ARTCFD_GEOMETRY_HPP
#define ARTCFD_GEOMETRY_HPP

#include <iostream>
#include "DimensionUtility.hpp"
#include "NumericalMethodUtility.hpp"
#include "GeometryData.hpp"
#include "GeometryLoader.hpp"

template <class Dimension, class NumericalMethodUtility> class GeometryBuilder;

template <class Dimension, class NumericalMethodUtility>
class GeometryBase{
public:
    GeometryBase(const std::string & geoName): geoName(geoName){}

    virtual bool Init();

    virtual bool Load(std::shared_ptr<ArtProject> &project){ return true;};

    GeometrySourceFormat getFormat() const {
        return format;
    }

    void setFormat(GeometrySourceFormat format) {
        GeometryBase::format = format;
    };
protected:
    std::string geoName;
    GeometrySourceFormat format{GeometrySourceFormat::File};
    std::shared_ptr<GeometryLoader<Dimension, NumericalMethodUtility>> geoLoader{nullptr};
};

template<class Dimension, class NumericalMethodUtility>
bool GeometryBase<Dimension, NumericalMethodUtility>::Init() {
    if(format == GeometrySourceFormat::File){
        geoLoader = std::make_shared<GeometryFileLoader<Dimension, NumericalMethodUtility>>(geoName);
        return true;
    }else
        return false;

}

template <class Dimension, class NumericalMethodUtility>
class Geometry : public GeometryBase<Dimension, NumericalMethodUtility>{
public:
    Geometry(const std::string &geoName) : GeometryBase<Dimension, NumericalMethodUtility>(geoName) {}

public:
    using GeoType = GeometryData<Dimension>;
};

template <class Dimension>
class Geometry<Dimension, MeshTypeMethod> : public GeometryBase<Dimension, MeshTypeMethod> {
public:
    using GeoType = GeometryMeshData<Dimension>;

    Geometry(const std::string &geoName) : GeometryBase<Dimension, MeshTypeMethod>(geoName) {
        data = std::make_shared<GeoType>();
    }

    static GeometryBuilder<Dimension, MeshTypeMethod> create(const std::string & geoName);

    bool Load(std::shared_ptr<ArtProject> &project) override;

    std::shared_ptr<GeoType> &getData() {
        return data;
    }

private:
    std::shared_ptr<GeoType> data{nullptr};
};

template<class Dimension>
GeometryBuilder<Dimension, MeshTypeMethod> Geometry<Dimension, MeshTypeMethod>::create(const std::string &geoName) {
    return GeometryBuilder<Dimension, MeshTypeMethod>(geoName);
}

template<class Dimension>
bool Geometry<Dimension, MeshTypeMethod>::Load(std::shared_ptr<ArtProject> &project) {
    this->geoLoader->setProj(project);
    this->geoLoader->setGeoData(data);
    this->geoLoader->load();
    return GeometryBase<Dimension, MeshTypeMethod>::Load(project);
}


template <class Dimension>
class Geometry<Dimension, FEM> : public GeometryBase<Dimension, FEM> {
public:
    using GeoType = GeometryMeshFemData<Dimension>;

    Geometry(const std::string &geoName) : GeometryBase<Dimension, MeshTypeMethod>(geoName) {
        data = std::make_shared<GeoType>();
    }

    static GeometryBuilder<Dimension, FEM> create(const std::string & geoName);

    bool Load(std::shared_ptr<ArtProject> &project) override;

    std::shared_ptr<GeoType> &getData() {
        return data;
    }

private:
    std::shared_ptr<GeoType> data{nullptr};
};

template<class Dimension>
GeometryBuilder<Dimension, FEM> Geometry<Dimension, FEM>::create(const std::string &geoName) {
    return GeometryBuilder<Dimension, FEM>(geoName);
}

template<class Dimension>
bool Geometry<Dimension, FEM>::Load(std::shared_ptr<ArtProject> &project) {
    this->geoLoader->setProj(project);
    this->geoLoader->setGeoData(data);
    this->geoLoader->load();
    return GeometryBase<Dimension, FEM>::Load(project);
}


template <class Dimension, class NumericalMethodUtility>
class GeometryBuilder{
public:
    GeometryBuilder(const std::string & geoName) {
        p_geometry = std::make_shared<Geometry<Dimension, NumericalMethodUtility>>(geoName);
        p_geometry->Init();
    }

    GeometryBuilder & load(std::shared_ptr<ArtProject> &project){
        p_geometry->Load(project);
        return *this;
    }

    std::shared_ptr<typename Geometry<Dimension, NumericalMethodUtility>::GeoType> build(){
        return p_geometry->getData();
    };

    GeometryBuilder & setSourceFormat(GeometrySourceFormat format){
        p_geometry->setFormat(format);
        return *this;
    }

private:
    std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> p_geometry{nullptr};
};




#endif //ARTCFD_GEOMETRY_HPP


