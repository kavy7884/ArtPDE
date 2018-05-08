//
// Created by Chingkai Chou on 5/7/18.
//

#ifndef ARTCFD_GEOMETRY_DATA_READER_HPP
#define ARTCFD_GEOMETRY_DATA_READER_HPP

#include "numerical_method_utility.hpp"
#include "geometry_data.hpp"

namespace art_pde {

    class GeometryReadFormat{};
    class ReadArtPDE : public GeometryReadFormat{};

    template<class GeometryDataType, typename GeometryReadFormat>
    class GeometryDataReader : virtual public GeometryDataType{
    public:
        GeometryDataReader(): GeometryDataType(){}

        bool getStatus(){
            return is_read_success;
        }

    protected:
        bool is_read_success{false};
    };
}


#endif //ARTCFD_GEOMETRY_DATA_READER_HPP
