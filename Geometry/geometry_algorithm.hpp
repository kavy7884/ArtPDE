//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_ALGORITHM_HPP
#define ARTPDE_KAVY_GEOMETRY_ALGORITHM_HPP

#include "geometry_mesh_type_data.hpp"

namespace art_pde{ namespace geometry {


    template <typename DataType>
    class Algorithm: public virtual DataType{
    public:
        Algorithm(): DataType(){
            std::cout << "Algorithm" << std::endl;
        }


    };
}}

#endif //ARTPDE_KAVY_GEOMETRY_ALGORITHM_HPP
