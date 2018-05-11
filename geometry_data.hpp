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
#include "edge.hpp"
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
            using GeoEdgeType = Edge<GeoPointType>;
            using PtrGeoEdgeType = std::shared_ptr<GeoEdgeType>;
            using VecPtrGeoEdgeType = std::vector<PtrGeoEdgeType>;
            using GeoCellType = Cell<GeoPointType>;
            using PtrGeoCellType = std::shared_ptr<GeoCellType>;
            using VecPtrGeoCellType = std::vector<PtrGeoCellType>;
        };

        GeometryData(){}

        const size_t getTotal_VertexNum() const{ return total_vertex.size();}
        const size_t getTotal_CellNum() const{ return total_cell.size();}

        typename Type::VecPtrGeoPointType getTotal_VecPtrPointOnVertex() const;
        typename Type::VecPtrGeoPointType getCell_VecPtrPointOnVertex(const size_t cell_id) const;
        typename Type::GeoCellType::CellDefineType getCell_CellDefineType(const size_t cell_id) const { return total_cell[cell_id]->getCell_define_Type();};

        double test{0.0};


    protected:
        typename Type::VecPtrGeoVertexType total_vertex;
        typename Type::VecPtrGeoEdgeType total_edge;
        typename Type::VecPtrGeoCellType total_cell;
    };

    template<typename Dimension, typename CoordinateBasis>
    typename GeometryData<art_pde::MeshTypeMethod, Dimension, CoordinateBasis>::Type::VecPtrGeoPointType
    GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>::getTotal_VecPtrPointOnVertex() const {
        typename GeometryData<art_pde::MeshTypeMethod, Dimension, CoordinateBasis>::Type::VecPtrGeoPointType
                reVecPtrGeoPointType;

        for (size_t i = 0; i < this->getTotal_VertexNum(); ++i) {
            reVecPtrGeoPointType.push_back(
                    this->total_vertex[i]->getPtr_point()
            );
        }

        return reVecPtrGeoPointType;
    }

    template<typename Dimension, typename CoordinateBasis>
    typename GeometryData<art_pde::MeshTypeMethod, Dimension, CoordinateBasis>::Type::VecPtrGeoPointType
    GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>::getCell_VecPtrPointOnVertex(const size_t cell_id) const {
        typename GeometryData<art_pde::MeshTypeMethod, Dimension, CoordinateBasis>::Type::VecPtrGeoPointType
                reVecPtrGeoPointType;

        //ToDo - Range check

        auto vec_ptr_vertex_in_cell_id = this->total_cell[cell_id]->getVec_ptr_vetex();

        for (size_t i = 0; i < this->total_cell[cell_id]->getNumVertex(); ++i) {
            reVecPtrGeoPointType.push_back(
                    vec_ptr_vertex_in_cell_id[i]->getPtr_point()
            );
        }

        return reVecPtrGeoPointType;
    }

}

#endif //ARTCFD_GEODATA_HPP
