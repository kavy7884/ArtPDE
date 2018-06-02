//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_DATA_MESH_TYPE_READER_HPP
#define ARTPDE_KAVY_GEOMETRY_DATA_MESH_TYPE_READER_HPP

#include "Project/art_project.hpp"
#include "geometric_reader_mesh_type.hpp"

namespace art_pde{ namespace geometry {

    template< size_t Dimension, typename ...DataType>
    class Reader : public ReadingMethod<DataType...>{
    public:
        Reader(): ReadingMethod<DataType...>() {}

        bool read(const std::shared_ptr<project::ArtProject>& art_project){
            std::cout << "read art_proj" << std::endl;
            return true;
        }
    };

}}

#endif //ARTPDE_KAVY_GEOMETRY_DATA_MESH_TYPE_READER_HPP
