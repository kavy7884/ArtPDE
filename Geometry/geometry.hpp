//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_HPP
#define ARTPDE_KAVY_GEOMETRY_HPP

#include "geometry_algorithm.hpp"
#include "geometry_reader.hpp"

namespace art_pde{ namespace geometry {

    template<typename AlgorithmType, typename ReaderType>
    class Geometry: public AlgorithmType, public ReaderType{
    public:
        Geometry() : AlgorithmType(), ReaderType(){}


    };
}}


#endif //ARTPDE_KAVY_GEOMETRY_HPP
