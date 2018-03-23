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
class Geometry {
public:
    using GeoType = GeometryData<Dimension>;
};

template <class Dimension>
class Geometry<Dimension, MeshTypeMethod> {
public:
    using GeoType = GeometryMeshData<Dimension>;
    Geometry() {
        data = std::make_shared<GeoType>();
    }
    const std::shared_ptr<GeoType> &getData() const {
        return data;
    }

private:
    std::shared_ptr<GeoType> data{nullptr};
};

template <class Dimension>
class Geometry<Dimension, FEM> {
public:
    using GeoType = GeometryMeshFemData<Dimension>;
    Geometry() {
        data = std::make_shared<GeoType>();
    }
    const std::shared_ptr<GeoType> &getData() const {
        return data;
    }

private:
    std::shared_ptr<GeoType> data{nullptr};
};
#endif //ARTCFD_GEOMETRY_HPP
