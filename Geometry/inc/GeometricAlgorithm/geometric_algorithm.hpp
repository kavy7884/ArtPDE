//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_ALGORITHM_HPP
#define ARTPDE_KAVY_GEOMETRY_ALGORITHM_HPP

#include "Geometry/inc/GeometricData/geometric_mesh_type_data.hpp"

namespace art_pde{ namespace geometry {

        template< size_t Dimension, typename ...DataType> class Algorithm;

        template<typename ...DataType>
        class Algorithm<1, DataType...> :public virtual DataType...{
        public:
            Algorithm<1, DataType...>(): DataType()...{
                std::cout << "Algorithm 1" << std::endl;
            }

        };

        template<typename ...DataType>
        class Algorithm<2, DataType...> :public virtual DataType...{
        public:
            Algorithm<2, DataType...>(): DataType()...{
                std::cout << "Algorithm 2" << std::endl;
            }

        };

        template<typename ...DataType>
        class Algorithm<3, DataType...> :public virtual DataType...{
        public:
            Algorithm<3, DataType...>(): DataType()...{
                std::cout << "Algorithm 3" << std::endl;
            }

        };

}}

#endif //ARTPDE_KAVY_GEOMETRY_ALGORITHM_HPP
