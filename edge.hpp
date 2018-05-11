//
// Created by Chingkai Chou on 5/10/18.
//

#ifndef ARTCFD_EDGE_HPP
#define ARTCFD_EDGE_HPP

#include <list>
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
        using CellDefineType = typename CellType::CellDefineType;

        Edge() {}

        CellDefineType getCell_define_Type() const {
            return cell_define_Type;
        }

        void setCell_define_Type(CellDefineType cell_define_Type) {
            Edge::cell_define_Type = cell_define_Type;
        }

    protected:
        VecPtrVertexType vec_ptr_vetex;
        CellDefineType cell_define_Type {CellDefineType::None};

    };

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
