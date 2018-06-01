//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_DATA_MESH_TYPE_READER_HPP
#define ARTPDE_KAVY_GEOMETRY_DATA_MESH_TYPE_READER_HPP

namespace art_pde{ namespace geometry {

    enum class SourceType{ArtPDE};

    template <typename DataType>
    class Reader : public virtual DataType{
    public:
        Reader(): DataType(){}
    };
}}

#endif //ARTPDE_KAVY_GEOMETRY_DATA_MESH_TYPE_READER_HPP
