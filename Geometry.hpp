//
// Created by Chingkai Chou on 2/23/18.
//

#ifndef ARTCFD_GEOMETRY_HPP
#define ARTCFD_GEOMETRY_HPP

#include <iostream>
#include "DimensionUtility.hpp"
#include "NumericalMethodUtility.hpp"
#include "GeometryData.hpp"

template <class Dimension, class NumericalMethodUtility>
class Geometry {};

template <class Dimension>
class Geometry<Dimension, FEM> {
    using GeoType = GeometryMeshFemData<Dimension>;
public:

    Geometry() {
        GeoData = std::make_shared<GeoType>();
    }
    const std::shared_ptr<GeoType> &getGeoData() const {
        return GeoData;
    }
private:
    std::shared_ptr<GeoType> GeoData{nullptr};
};
#endif //ARTCFD_GEOMETRY_HPP
