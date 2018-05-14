//
// Created by Chingkai Chou on 5/8/18.
//

#ifndef ARTCFD_GEOMETRY_DATA_ALGORITHM_HPP
#define ARTCFD_GEOMETRY_DATA_ALGORITHM_HPP

namespace art_pde {

    template<class GeometryDataType>
    class GeometryDataAlgorithm : virtual public GeometryDataType{
    public:

        GeometryDataAlgorithm() : GeometryDataType() {}

        typename GeometryDataType::Type::VecPtrGeoPointType getTotal_VecPtrPointOnCellCenter();
        const typename GeometryDataType::Type::VecPtrGeoCellType getVertex_VecPtrNeighborCell(const size_t vertex_id);
        void calEdge(){
            calEdgeMerge();
        }


    private:
        bool is_cal_vertex_neighbor_cell{ false };
        bool is_gen_cell_center_point{ false };
        bool is_cal_edge{ false };


        void genCellCenterPoint();
        void calVertexNeighborCell();
        void calEdgeMerge();
        

    };

    template<class GeometryDataType>
    typename GeometryDataType::Type::VecPtrGeoPointType GeometryDataAlgorithm<GeometryDataType>::getTotal_VecPtrPointOnCellCenter() {
        typename GeometryDataType::Type::VecPtrGeoPointType reVecPtrGeoPointType;

        if(!is_gen_cell_center_point) { genCellCenterPoint();}

        for (size_t i = 0; i < this->getTotal_CellNum(); ++i) {
            reVecPtrGeoPointType.push_back(
                    this->total_cell[i]->getPtr_cell_center_point()
            );
        }
        return reVecPtrGeoPointType;
    }

    template<class GeometryDataType>
    void GeometryDataAlgorithm<GeometryDataType>::genCellCenterPoint() {
        for (size_t i = 0; i < this->getTotal_CellNum(); ++i) {
            this->total_cell[i]->calPtr_cell_center_point();
        }
        is_gen_cell_center_point = true;
    }

    template<class GeometryDataType>
    void GeometryDataAlgorithm<GeometryDataType>::calVertexNeighborCell() {
        for (size_t i = 0; i < this->getTotal_CellNum(); ++i) {
            auto this_ptr_cell = this->total_cell[i];
            auto vec_ptr_vetex_in_cell = this_ptr_cell->getVec_ptr_vetex();
            for (size_t j = 0; j < vec_ptr_vetex_in_cell.size() ; ++j) {
                vec_ptr_vetex_in_cell[j]->addPtrNeighborCell(this_ptr_cell);
            }
        }
        is_cal_vertex_neighbor_cell = true;
    }

    template<class GeometryDataType>
    const typename GeometryDataType::Type::VecPtrGeoCellType
    GeometryDataAlgorithm<GeometryDataType>::getVertex_VecPtrNeighborCell(const size_t vertex_id) {
        if(!is_cal_vertex_neighbor_cell) calVertexNeighborCell();
        if(!is_gen_cell_center_point) genCellCenterPoint();
        return this->total_vertex[vertex_id]->getVec_ptr_neighbor_cell();
    }

    template<class GeometryDataType>
    void GeometryDataAlgorithm<GeometryDataType>::calEdgeMerge() {

        for (int i = 0; i < this->getTotal_VertexNum(); ++i) {
            auto & list_ptr_neighbor_edge = this->total_vertex[i]->getList_ptr_neighbor_edge();

            auto it_1 = list_ptr_neighbor_edge.begin();
            while( it_1 != list_ptr_neighbor_edge.end()){

                if((*it_1)->Is_merged()){
                    ++it_1;
                    continue;
                }

                auto it_2 = it_1;
                ++it_2;
                while( it_2 != list_ptr_neighbor_edge.end()) {

                    if((*it_2)->Is_merged()){
                        ++it_2;
                        continue;
                    }

                    if((*(*it_1))==(*(*it_2)))
                    {
                        //std::cout << "Merge" << std::endl;

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
