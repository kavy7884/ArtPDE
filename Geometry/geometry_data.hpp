//
// Created by Chingkai Chou on 5/2/18.
//

#ifndef ARTCFD_GEODATA_HPP
#define ARTCFD_GEODATA_HPP

#include <vector>
#include <memory>
#include "../Utility/numerical_method_utility.hpp"
#include "point.hpp"
#include "vertex.hpp"
#include "edge.hpp"
#include "cell.hpp"
#include "geometry_data_algorithm.hpp"

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
            using GeoEdgeType = Edge<GeoPointType>;
            using PtrGeoEdgeType = std::shared_ptr<GeoEdgeType>;
            using VecPtrGeoEdgeType = std::vector<PtrGeoEdgeType>;
            using ListPtrGeoEdgeType = std::list<PtrGeoEdgeType>;
            using GeoCellType = Cell<GeoPointType>;
            using PtrGeoCellType = std::shared_ptr<GeoCellType>;
            using VecPtrGeoCellType = std::vector<PtrGeoCellType>;
        };

        GeometryData(){}

        const size_t getNum_TotalVertex() const{ return total_vertex.size();}
        const size_t getNum_TotalCell() const{ return total_cell.size();}

        const typename Type::PtrGeoPointType& getVertex_PtrPoint(const size_t vertex_id) const;

        const typename Type::PtrGeoPointType& getCell_Center_PtrPoint(const size_t cell_id) const;
        const typename Type::GeoCellType::CellDefineType& getCell_CellDefineType(const size_t cell_id) const { return total_cell[cell_id]->getCell_define_Type();};

        double test{0.0};

    protected:
        typename Type::VecPtrGeoVertexType total_vertex;
        typename Type::VecPtrGeoCellType total_cell;

    };


    template<typename Dimension, typename CoordinateBasis>
    const typename GeometryData<art_pde::MeshTypeMethod, Dimension, CoordinateBasis>::Type::PtrGeoPointType&
    GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>::getVertex_PtrPoint(const size_t vertex_id) const {
        return total_vertex[vertex_id]->getPtr_point();
    };

    template<typename Dimension, typename CoordinateBasis>
    const typename GeometryData<art_pde::MeshTypeMethod, Dimension, CoordinateBasis>::Type::PtrGeoPointType&
    GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>::getCell_Center_PtrPoint(const size_t cell_id) const {
        total_cell[cell_id]->calPtr_cell_center_point();
        return total_cell[cell_id]->getPtr_cell_center_point();
    };


}

#endif //ARTCFD_GEODATA_HPP
