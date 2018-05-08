//
// Created by Chingkai Chou on 5/2/18.
//

#ifndef ARTCFD_GEODATA_HPP
#define ARTCFD_GEODATA_HPP

#include <vector>
#include <memory>
#include "numerical_method_utility.hpp"
#include "point.hpp"
#include "vertex.hpp"
#include "cell.hpp"

namespace art_pde {

    // Geometry data abstract define
    template <typename GeometryDescription, typename Dimension, typename CoordinateBasis>
    class GeometryData;

    // Geometry data MeshTypeMethod define
    template <typename Dimension, typename CoordinateBasis>
    class GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>{
    public:
        struct Type{
            using GeoDescribeType = MeshTypeMethod;
            using GeoDimType = Dimension;
            using GeoCoordType = CoordinateBasis;
            using GeoPointType = Point<Dimension, CoordinateBasis>;
            using PtrGeoPointType = std::shared_ptr<GeoPointType>;
            using VecPtrGeoPointType = std::vector<PtrGeoPointType>;
            using GeoVertexType = Vertex<GeoPointType>;
            using PtrGeoVertexType = std::shared_ptr<GeoVertexType>;
            using VecPtrGeoVertexType = std::vector<PtrGeoVertexType>;
            using GeoCellType = Cell<GeoPointType>;
            using PtrGeoCellType = std::shared_ptr<GeoCellType>;
            using VecPtrGeoCellType = std::vector<PtrGeoCellType>;
        };

        GeometryData(){}

        const size_t getTotalVertexNum() const{ return total_vertex.size();}




        double test{0.0};


    protected:
        typename Type::VecPtrGeoVertexType total_vertex;
        typename Type::VecPtrGeoCellType total_cell;
    };

}

#endif //ARTCFD_GEODATA_HPP
