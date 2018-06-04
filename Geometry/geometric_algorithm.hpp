//
// Created by Chingkai Chou on 6/3/18.
//

#ifndef ARTPDE_KAVY_GEOMETRIC_ALGORITHM_HPP
#define ARTPDE_KAVY_GEOMETRIC_ALGORITHM_HPP

#include "geometric_data.hpp"

namespace art_pde{ namespace geometry {
        namespace mesh_type {
            namespace Dim2 {

                template <typename GeometricDataType>
                class GeometricAlgorithm: public virtual GeometricDataType{
                public:
                    GeometricAlgorithm():GeometricDataType(){
                        std::cout << "Dim2::GeometricAlgorithm" << std::endl;
                    }
                };

            }

            namespace Dim3 {

                template <typename GeometricDataType>
                class GeometricAlgorithm: public virtual GeometricDataType{
                public:
                    GeometricAlgorithm():GeometricDataType(){
                        std::cout << "Dim3::GeometricAlgorithm" << std::endl;
                    }
                };

            }
        }
}}

#endif //ARTPDE_KAVY_GEOMETRIC_ALGORITHM_HPP
