//
// Created by Chingkai Chou on 6/2/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_READER_MESH_TYPE_HPP
#define ARTPDE_KAVY_GEOMETRY_READER_MESH_TYPE_HPP

#include "Geometry/inc/GeometricData/geometric_mesh_type_data.hpp"

namespace art_pde{ namespace geometry {

        template <typename GeoDataType> class ReadingMethod{};

        template <size_t Dimension, typename CalculatedPointType>
        class ReadingMethod<typename mesh_type::Data<Dimension, CalculatedPointType>>
                :public virtual mesh_type::Data<Dimension, CalculatedPointType>{
        public:
            ReadingMethod<typename mesh_type::Data<Dimension, CalculatedPointType>>()
                    : mesh_type::Data<Dimension, CalculatedPointType>(){}

        };

}}


#endif //ARTPDE_KAVY_GEOMETRY_READER_MESH_TYPE_HPP
