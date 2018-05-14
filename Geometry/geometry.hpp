//
// Created by Chingkai Chou on 5/2/18.
//

#ifndef ARTCFD_GEOMESH_HPP
#define ARTCFD_GEOMESH_HPP

#include "../Utility/numerical_method_utility.hpp"
#include "geometry_data.hpp"
#include "geometry_data_reader.hpp"
#include "geometry_data_reader_artpde.hpp"
#include "geometry_data_algorithm.hpp"

namespace art_pde {

    template<class GeometryDataType, template<class> class GeometryDataReaderType>
    class Geometry:
            virtual public GeometryDataType,
            public GeometryDataReaderType<GeometryDataType>,
            public GeometryDataAlgorithm<GeometryDataType>{
    public:
        Geometry() :
                GeometryDataType(),
                GeometryDataReaderType<GeometryDataType>(),
                GeometryDataAlgorithm<GeometryDataType>(){}

        bool read(const std::shared_ptr<ArtProject>& art_project){
            GeometryDataReaderType<GeometryDataType>::read(art_project);
            return GeometryDataReaderType<GeometryDataType>::getStatus();
        }

    };

}

#endif //ARTCFD_GEOMESH_HPP
