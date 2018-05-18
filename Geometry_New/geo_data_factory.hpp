//
// Created by Chingkai Chou on 5/18/18.
//

#ifndef ARTPDE_MARVIN_GEO_DATA_FACTORY_HPP
#define ARTPDE_MARVIN_GEO_DATA_FACTORY_HPP

#include "geo_data.hpp"

template <typename Data>
class QuadCell: public Cell<Data>{
public:
    using VertexType = Vertex<Data>;
    using PtrVertexType = std::shared_ptr<VertexType>;
    using VecPtrVertexType = std::vector<PtrVertexType>;
    using EdgeType = Edge<Data>;
    using PtrEdgeType = std::shared_ptr<EdgeType>;

    QuadCell(const PtrVertexType& v1, const PtrVertexType& v2, const PtrVertexType& v3, const PtrVertexType& v4)
            : Cell<Data>(){

        std::shared_ptr<QuadCell<Data>> cell(this);
        PtrEdgeType edge;

        edge = genEdge(v1, v2);
        edge->addParent(cell);
        this->addChild(edge);

        edge = genEdge(v2, v3);
        edge->addParent(cell);
        this->addChild(edge);

        edge = genEdge(v3, v4);
        edge->addParent(cell);
        this->addChild(edge);

        edge = genEdge(v4, v1);
        edge->addParent(cell);
        this->addChild(edge);

    }

private:
    VecPtrVertexType vec_ptr_vertex;


    PtrEdgeType genEdge(const PtrVertexType& v1, const PtrVertexType& v2){
        PtrEdgeType reEdge = std::make_shared<EdgeType>();
        v1->addParent(reEdge);
        v2->addParent(reEdge);
        reEdge->addChild(v1);
        reEdge->addChild(v2);
        return reEdge;
    }
};

#endif //ARTPDE_MARVIN_GEO_DATA_FACTORY_HPP
