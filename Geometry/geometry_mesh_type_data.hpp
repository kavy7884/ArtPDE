//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_DATA_HPP
#define ARTPDE_KAVY_GEOMETRY_DATA_HPP

#include "inc/GeometricTree/geometric_tree_units.hpp"

namespace art_pde{ namespace geometry {
    namespace mesh_type{

        template <size_t Dimension, template<size_t> class CalculatedPointType>
        class Data{
        public:
            Data() {
                std::cout << "Constructor" << std::endl;
            }

            geometric_tree::Vertex<CalculatedPointType<Dimension>> vertex{1,2,3};
        };
    }
}}

#endif //ARTPDE_KAVY_GEOMETRY_DATA_HPP
