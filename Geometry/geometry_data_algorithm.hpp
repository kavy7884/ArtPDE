//
// Created by Chingkai Chou on 5/8/18.
//

#ifndef ARTCFD_GEOMETRY_DATA_ALGORITHM_HPP
#define ARTCFD_GEOMETRY_DATA_ALGORITHM_HPP

#include "set"

namespace art_pde {

    template<class GeometryDataType>
    class GeometryDataAlgorithm : virtual public GeometryDataType{
    public:
        using PtrGeoPointType = typename GeometryDataType::Type::PtrGeoPointType;
        using VecPtrGeoVertexType = typename GeometryDataType::Type::VecPtrGeoVertexType;
        using VecPtrGeoCellType = typename GeometryDataType::Type::VecPtrGeoCellType;
        using PtrGeoEdgeType = typename GeometryDataType::Type::PtrGeoEdgeType;
        using VecPtrGeoEdgeType = typename GeometryDataType::Type::VecPtrGeoEdgeType;

        GeometryDataAlgorithm() : GeometryDataType() {}

        const size_t getNum_TotalEdge(){
            if(!this->is_Cal_edge()){
                this->calEdge();
            }
            return this->total_edge.size();
        }

        const PtrGeoPointType& getEdge_Center_PtrPoint(const size_t edge_id);
        const VecPtrGeoVertexType& getRelation_Cell_Neighbor_Vertex(const size_t cell_id) const; // Relation 3-0
        const VecPtrGeoCellType& getRelation_Vertex_Neighbor_Cell(const size_t vertex_id) const; // Relation 0-3

        void calEdge();

        bool is_Cal_edge() const { return cal_edge;};

    protected:
        bool cal_edge{ false };
        void calEdgeMerge();

        VecPtrGeoEdgeType total_edge;
    };

    template<class GeometryDataType>
    const typename GeometryDataType::Type::VecPtrGeoVertexType&
    GeometryDataAlgorithm<GeometryDataType>::getRelation_Cell_Neighbor_Vertex(const size_t cell_id) const {
        return this->total_cell[cell_id]->getVec_ptr_vetex();
    };

    template<class GeometryDataType>
    const typename GeometryDataType::Type::VecPtrGeoCellType&
    GeometryDataAlgorithm<GeometryDataType>::getRelation_Vertex_Neighbor_Cell(const size_t vertex_id) const {
        return this->total_vertex[vertex_id]->getVec_ptr_neighbor_cell();
    };

    template<class GeometryDataType>
    void GeometryDataAlgorithm<GeometryDataType>::calEdge(){
        calEdgeMerge();

        std::set<PtrGeoEdgeType> edge_set;

        for (auto &ptr_cell: this->total_cell) {
            for (auto &ptr_edge: ptr_cell->getList_ptr_neighbor_edge()) {
                edge_set.insert(ptr_edge);
            }
        }

        this->total_edge.insert(this->total_edge.begin(), edge_set.begin(), edge_set.end());
        this->cal_edge = true;
    }


    template<class GeometryDataType>
    const typename GeometryDataType::Type::PtrGeoPointType&
    GeometryDataAlgorithm<GeometryDataType>::getEdge_Center_PtrPoint(const size_t edge_id) {
        if(!this->is_Cal_edge()){
            this->calEdge();
        }
        this->total_edge[edge_id]->calPtr_cell_center_point();
        return this->total_edge[edge_id]->getPtr_edge_center_point();
    };

    template<class GeometryDataType>
    void GeometryDataAlgorithm<GeometryDataType>::calEdgeMerge() {

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
    }


}

#endif //ARTCFD_GEOMETRY_DATA_ALGORITHM_HPP
