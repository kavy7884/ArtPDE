//
// Created by Chingkai Chou on 3/8/18.
//

#ifndef ARTCFD_GEOMETRYDATA_HPP
#define ARTCFD_GEOMETRYDATA_HPP

#include <vector>
#include "Point.hpp"
#include "Connectivity.hpp"
#include "Element.hpp"

template <class Dimension>
class GeometryData{
public:
    GeometryData() {}
    size_t GeoDim{Dimension::Dim};

};

template <class Dimension>
class GeometryMeshData : public GeometryData<Dimension>{
public:
    GeometryMeshData(): GeometryData<Dimension>(){
        xNode = std::make_shared<PointList<Dimension>>();
        cElement30 = std::make_shared<ElementList<Dimension>>();
    };

    std::shared_ptr<PointList<Dimension>> xNode{nullptr};
    std::shared_ptr<ElementList<Dimension>> cElement30{nullptr};

    std::vector<std::shared_ptr<GeoElement<Dimension>>> geoElement;

private:

};

template <class Dimension>
class GeometryMeshFemData : public GeometryMeshData<Dimension>{
public:
    GeometryMeshFemData():GeometryMeshData<Dimension>(){};



};


#endif //ARTCFD_GEOMETRYDATA_HPP
