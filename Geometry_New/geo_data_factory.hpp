//
// Created by Chingkai Chou on 5/18/18.
//

#ifndef ARTPDE_GEO_DATA_FACTORY_HPP
#define ARTPDE_GEO_DATA_FACTORY_HPP

#include "geo_data.hpp"

template <typename Data>
class QuadFace: public Face<Data>{
public:
    using VertexType = Vertex<Data>;
    using PtrVertexType = std::shared_ptr<VertexType>;
    using VecPtrVertexType = std::vector<PtrVertexType>;
    using EdgeType = Edge<Data>;
    using PtrEdgeType = std::shared_ptr<EdgeType>;

    QuadFace(const PtrVertexType& v1, const PtrVertexType& v2, const PtrVertexType& v3, const PtrVertexType& v4)
            : Face<Data>(){

        std::shared_ptr<QuadFace<Data>> face(this);
        PtrEdgeType edge;

        edge = genEdge(v1, v2);
        edge->addParent(face);
        this->addChild(edge);

        edge = genEdge(v2, v3);
        edge->addParent(face);
        this->addChild(edge);

        edge = genEdge(v3, v4);
        edge->addParent(face);
        this->addChild(edge);

        edge = genEdge(v4, v1);
        edge->addParent(face);
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

#endif //ARTPDE_GEO_DATA_FACTORY_HPP
