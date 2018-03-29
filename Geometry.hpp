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
template <class Dimension, class NumericalMethodUtility> class GeometryType;

template <class Dimension, class NumericalMethodUtility>
class GeometryBase{
public:
    static GeometryBuilder<Dimension, NumericalMethodUtility> create(const std::string & geoName){
        return GeometryBuilder<Dimension, NumericalMethodUtility>(geoName);
    }

protected:
    GeometryBase(){};

    template <class Dimension_, class NumericalMethodUtility_> friend class GeometryBuilder;

    virtual bool Init();

    virtual bool Load(std::shared_ptr<ArtProject> &project);

    GeometrySourceFormat getFormat() const { return format; }

    void setFormat(GeometrySourceFormat format) { GeometryBase::format = format; };

    const std::string &getGeoName() const{ return geoName; }

    void setGeoName(const std::string &geoName){GeometryBase::geoName = geoName;}

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

template<class Dimension, class NumericalMethodUtility>
bool GeometryBase<Dimension, NumericalMethodUtility>::Load(std::shared_ptr<ArtProject> &project) {
    this->geoLoader->setProj(project);
    this->geoLoader->load();
    return true;
}

template <class Dimension, class NumericalMethodUtility>
class Geometry : public GeometryBase<Dimension, NumericalMethodUtility>, public GeometryType<Dimension, NumericalMethodUtility>{
    struct GeometryDummy {};
public:
    explicit Geometry(GeometryDummy):GeometryBase<Dimension, NumericalMethodUtility>() {
        data = std::make_shared<typename GeometryType<Dimension, NumericalMethodUtility>::GeoType>();
    }

private:
    template <class Dimension_, class NumericalMethodUtility_> friend class GeometryBuilder;

    static std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> make_shared_Geometry() {
        return std::make_shared<Geometry<Dimension, NumericalMethodUtility>>(GeometryDummy());
    }

    Geometry() = delete;

    bool Load(std::shared_ptr<ArtProject> &project) override;

    std::shared_ptr<typename GeometryType<Dimension, NumericalMethodUtility>::GeoType> &getData(){
        return data;
    };

    std::shared_ptr<typename GeometryType<Dimension, NumericalMethodUtility>::GeoType> data{nullptr};
};

template<class Dimension, class NumericalMethodUtility>
bool Geometry<Dimension, NumericalMethodUtility>::Load(std::shared_ptr<ArtProject> &project) {
    this->geoLoader->setGeoData(data);
    return GeometryBase<Dimension, NumericalMethodUtility>::Load(project);
}

template <class Dimension, class NumericalMethodUtility>
class GeometryBuilder{
public:
    GeometryBuilder(const std::string & geoName) {
        p_geometry = Geometry<Dimension, NumericalMethodUtility>::make_shared_Geometry();
        p_geometry->setGeoName(geoName);
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

template <class Dimension, class NumericalMethodUtility>
class GeometryType{ public: using GeoType = GeometryData<Dimension>; };

template <class Dimension>
class GeometryType<Dimension, MeshTypeMethod>{ public: using GeoType = GeometryMeshData<Dimension>; };

template <class Dimension>
class GeometryType<Dimension, FEM>{ public: using GeoType = GeometryMeshFemData<Dimension>; };


#endif //ARTCFD_GEOMETRY_HPP


