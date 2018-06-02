//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_HPP
#define ARTPDE_KAVY_GEOMETRY_HPP

#include "Project/art_project.hpp"
#include "Geometry/inc/GeometricAlgorithm/geometric_algorithm.hpp"
#include "Geometry/inc/GeometricReader/geometric_reader.hpp"

namespace art_pde{ namespace geometry {

    template<typename AlgorithmType, typename ReaderType>
    class Geometry: public AlgorithmType, public ReaderType{
    public:
        Geometry() : AlgorithmType(), ReaderType(){}

//        bool read(const std::shared_ptr<project::ArtProject>& art_project){
//            ReaderType::read(art_project);
//            return true;
//        }
    };
}}


#endif //ARTPDE_KAVY_GEOMETRY_HPP
