//
// Created by Chingkai Chou on 5/2/18.
//

#ifndef ARTCFD_GEOMESH_HPP
#define ARTCFD_GEOMESH_HPP

#include "numerical_method_utility.hpp"
#include "geometry_data.hpp"

namespace art_pde {

    template <class NumericalMethodUtility>
    class Geometry : public GeometryData<NumericalMethodUtility> {
    public:
        Geometry(const size_t num_geo_dim) : GeometryData(num_geo_dim){}

    public:
        



    protected:


    private:

    };



}

#endif //ARTCFD_GEOMESH_HPP
