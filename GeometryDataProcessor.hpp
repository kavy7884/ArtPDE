//
// Created by Chingkai Chou on 3/29/18.
//

#ifndef ARTCFD_GEOMETRYDATAPROCESSOR_HPP
#define ARTCFD_GEOMETRYDATAPROCESSOR_HPP

#include <memory>
#include "GeometryData.hpp"
#include "Geometry.hpp"

template <class Dimension, class DataType>
class GeometryDataProcessor: public DataType{
public:
    GeometryDataProcessor(const std::shared_ptr<DataType> &pData) : pData(pData) {}

private:
    std::shared_ptr<DataType> pData{nullptr};
};

template <class Dimension>
class GeometryDataProcessor<Dimension, GeometryMeshData<Dimension>>: public GeometryDataProcessor<Dimension, GeometryData<Dimension>>{
public:
    using DataType = GeometryMeshData<Dimension>;
    GeometryDataProcessor(const std::shared_ptr<GeometryMeshData<Dimension>> &pData) : pData(pData),
            GeometryDataProcessor<Dimension, GeometryData<Dimension>>(pData) {};

    bool process(){
        this->calGeoMeshData();
        return true;
    }

protected:
    bool calGeoMeshData(){
        newGeoElement();
        calGeoElementType();
        calGeoElementVolume();
        return true;
    }

private:
    void newGeoElement();
    void calGeoElementType();
    void calGeoElementVolume();

    std::shared_ptr<DataType> pData{nullptr};
};

template <class Dimension>
class GeometryDataProcessor<Dimension, GeometryMeshFemData<Dimension>>: public GeometryDataProcessor<Dimension, GeometryMeshData<Dimension>> {
public:
    using DataType = GeometryMeshFemData<Dimension>;

    GeometryDataProcessor(const std::shared_ptr<GeometryMeshFemData<Dimension>> &pData) : pData(pData),
                                                                                       GeometryDataProcessor<Dimension, GeometryMeshData<Dimension>>(
                                                                                               pData) {};
    bool process() {
        this->calGeoMeshData();
        return true;
    }

private:
    std::shared_ptr<DataType> pData{nullptr};
};

template <class Dimension>
void GeometryDataProcessor<Dimension, GeometryMeshData<Dimension>>::newGeoElement(){
    pData->geoElement = std::vector<std::shared_ptr<GeoElement<Dimension>>>(pData->cElement30->size(), nullptr);
    for (size_t i = 0; i < pData->cElement30->size(); ++i) {
        pData->geoElement[i] = std::make_shared<GeoElement<Dimension>>(pData->cElement30->getElementConnectivity(i));
    }
}

template<class Dimension>
void GeometryDataProcessor<Dimension, GeometryMeshData<Dimension>>::calGeoElementType() {
    for (auto &element: pData->geoElement) element->calGeoElementType();
}

template<class Dimension>
void GeometryDataProcessor<Dimension, GeometryMeshData<Dimension>>::calGeoElementVolume() {
    for (auto &element: pData->geoElement) element->calGeoElementVolume();
}

#endif //ARTCFD_GEOMETRYDATAPROCESSOR_HPP
