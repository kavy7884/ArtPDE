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

        bool isLoad_fundamental_data() const {
            return load_fundamental_data;
        }

        bool isMerged_edge() const {
            return merged_edge;
        }

        void setLoad_fundamental_data(bool load_fundamental_data) {
            GeometryData::load_fundamental_data = load_fundamental_data;
        }

        const size_t getNum_TotalVertex() const{ return total_vertex.size();}
        const size_t getNum_TotalCell() const{ return total_cell.size();}

        bool MergeEdge();

        double test{0.0};

    protected:
        typename Type::VecPtrGeoVertexType total_vertex;
        typename Type::VecPtrGeoCellType total_cell;
        bool load_fundamental_data{false};
        bool merged_edge{false};
    };

    template<typename Dimension, typename CoordinateBasis>
    bool GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>::MergeEdge() {

        for (int i = 0; i < this->getNum_TotalVertex(); ++i) {
            auto & list_ptr_neighbor_edge = this->total_vertex[i]->getList_ptr_neighbor_edge();

            auto it_1 = list_ptr_neighbor_edge.begin();
            while( it_1 != list_ptr_neighbor_edge.end()){

                if((*it_1)->is_Merged_edge()){
                    ++it_1;
                    continue;
                }

                auto it_2 = it_1;
                ++it_2;
                while( it_2 != list_ptr_neighbor_edge.end()) {

                    if((*it_2)->is_Merged_edge()){
                        ++it_2;
                        continue;
                    }

                    if((*(*it_1))==(*(*it_2)))
                    {
                        std::cout << "Merge" << std::endl;

                        (*it_1)->mergeEdge(*(*it_2));

                        for(auto & replace_cell : (*it_2)->getVec_ptr_neighbor_cell()){
                            replace_cell->replaceEdge((*it_2), (*it_1));
                        }

                        for(auto & ptr_vetex: (*it_2)->getVec_ptr_vetex()){
                            ptr_vetex->eraseEdge((*it_2));
                        }
                        break;
                    }
                    ++it_2;
                }
                ++it_1;
            }
        }

        this->merged_edge = true;
        return this->merged_edge;
    }

    // Geometry data MeshTypeMethod define
    template <typename Dimension, typename CoordinateBasis>
    class GeometryData<FEM, Dimension, CoordinateBasis> : public GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>{
    public:
        GeometryData(): GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>() {
            if(!this->isMerged_edge()){
                this->MergeEdge();
            }
        }


    };

//    template<typename Dimension, typename CoordinateBasis>
//    const typename GeometryData<art_pde::MeshTypeMethod, Dimension, CoordinateBasis>::Type::PtrGeoPointType&
//    GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>::getVertex_PtrPoint(const size_t vertex_id) const {
//        return this->total_vertex[vertex_id]->getPtr_point();
//    };
//
//    template<typename Dimension, typename CoordinateBasis>
//    const typename GeometryData<art_pde::MeshTypeMethod, Dimension, CoordinateBasis>::Type::PtrGeoPointType&
//    GeometryData<MeshTypeMethod, Dimension, CoordinateBasis>::getCell_Center_PtrPoint(const size_t cell_id) const {
//        this->total_cell[cell_id]->calPtr_cell_center_point();
//        return this->total_cell[cell_id]->getPtr_cell_center_point();
//    };


}

#endif //ARTCFD_GEODATA_HPP
