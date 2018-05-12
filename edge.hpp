//
// Created by Chingkai Chou on 5/10/18.
//

#ifndef ARTCFD_EDGE_HPP
#define ARTCFD_EDGE_HPP

#include <set>
#include <algorithm>
#include "vertex.hpp"
#include "cell.hpp"

namespace art_pde {

    template <typename PointType> class Vertex;
    template <typename PointType> class Cell;


    template <typename PointType>
    class Edge{
    public:
        using PtrPointType = std::shared_ptr<PointType>;
        using VertexType = art_pde::Vertex<PointType>;
        using PtrVertexType = std::shared_ptr<VertexType>;
        using VecPtrVertexType = std::vector<PtrVertexType>;
        using CellType = art_pde::Cell<PointType>;
        using PtrCellType = std::shared_ptr<CellType>;
        using VecPtrCellType = std::vector<PtrCellType>;
        using CellDefineType = typename CellType::CellDefineType;

        Edge() {}

        const size_t getVertexNum() const { return vec_ptr_vetex.size(); }

        CellDefineType getCell_define_Type() const {
            return cell_define_Type;
        }

        void setCell_define_Type(CellDefineType cell_define_Type) {
            Edge::cell_define_Type = cell_define_Type;
        }

        const VecPtrCellType &getVec_ptr_neighbor_cell() const {
            return vec_ptr_neighbor_cell;
        }

        void setVec_ptr_neighbor_cell(const VecPtrCellType &vec_ptr_neighbor_cell) {
            Edge::vec_ptr_neighbor_cell = vec_ptr_neighbor_cell;
        }

        const VecPtrVertexType &getVec_ptr_vetex() const;

        void addPtrNeighborCell(const PtrCellType &ptr_neighbor_cell){ vec_ptr_neighbor_cell.push_back(ptr_neighbor_cell); }

        bool operator==(const Edge &rhs) const;

        bool operator!=(const Edge &rhs) const {
            return !(rhs == *this);
        }

    protected:
        VecPtrVertexType vec_ptr_vetex;
        CellDefineType cell_define_Type {CellDefineType::None};
        VecPtrCellType vec_ptr_neighbor_cell;

    };

    template<typename PointType>
    bool Edge<PointType>::operator==(const Edge &rhs) const {
        if(this->getVertexNum() != rhs.getVertexNum()) return false;
        else{
            std::set<PtrVertexType> set_self_vertex(vec_ptr_vetex.cbegin(), vec_ptr_vetex.cend());
            std::set<PtrVertexType> set_rhs_vertex(rhs.getVec_ptr_vetex().cbegin(), rhs.getVec_ptr_vetex().cend());
            if ( std::equal(set_self_vertex.begin(), set_self_vertex.end(), set_rhs_vertex.begin() ))return true;
            else return false;
        }
    }

    template<typename PointType>
    const typename Edge<PointType>::VecPtrVertexType &Edge<PointType>::getVec_ptr_vetex() const {
        return vec_ptr_vetex;
    }

    template <typename PointType>
    class LineEdge : public Edge<PointType>{
        using VertexType = art_pde::Vertex<PointType>;
        using PtrVertexType = std::shared_ptr<VertexType>;
        using CellType = art_pde::Cell<PointType>;
        using CellDefineType = typename CellType::CellDefineType;

    public:
        LineEdge(const PtrVertexType &v_1, const PtrVertexType &v_2)
                :Edge<PointType>()
        {
            this->vec_ptr_vetex.push_back(v_1); this->vec_ptr_vetex.push_back(v_2);
            this->setCell_define_Type(CellDefineType::Line);
        }
    };

}

#endif //ARTCFD_EDGE_HPP
