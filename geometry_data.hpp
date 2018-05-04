//
// Created by Chingkai Chou on 5/2/18.
//

#ifndef ARTCFD_GEODATA_HPP
#define ARTCFD_GEODATA_HPP

#include "point.hpp"
#include "numerical_method_utility.hpp"

namespace art_pde {

    // Geometry data abstract define
    template <typename GeometryDataType, typename Dimension, typename CoordinateBasis>
    class GeometryData;

    // Geometry data MeshTypeMethod define
    template <typename Dimension, typename CoordinateBasis>
    class GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>{
    public:

        struct GeometryData_Traits{
            using PointType = art_pde::Point<Dimension, CoordinateBasis>;
            using GeoDataType = MeshTypeMethod;
        };

        GeometryData() {}


        typename GeometryData_Traits::PointType PT;




    };












}

#endif //ARTCFD_GEODATA_HPP
