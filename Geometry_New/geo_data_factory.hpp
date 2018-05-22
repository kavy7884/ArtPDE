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

    QuadFace() : Face<Data>(){
    }

    void create(const PtrVertexType& v0, const PtrVertexType& v1, const PtrVertexType& v2, const PtrVertexType& v3){
        auto face = this->shared_from_this();

        PtrEdgeType edge;

        edge = genEdge(v0, v1);
        edge->addParent(face);
        this->addChild(edge);

        edge = genEdge(v1, v2);
        edge->addParent(face);
        this->addChild(edge);

        edge = genEdge(v2, v3);
        edge->addParent(face);
        this->addChild(edge);

        edge = genEdge(v3, v0);
        edge->addParent(face);
        this->addChild(edge);

    }

private:

    PtrEdgeType genEdge(const PtrVertexType& v0, const PtrVertexType& v1){
        PtrEdgeType reEdge = std::make_shared<EdgeType>();
        v0->addParent(reEdge);
        v1->addParent(reEdge);
        reEdge->addChild(v0);
        reEdge->addChild(v1);
        return reEdge;
    }
};

template <typename Data>
class HexaCell: public Cell<Data>{
public:
    using VertexType = Vertex<Data>;
    using PtrVertexType = std::shared_ptr<VertexType>;
    using FaceType = Face<Data>;
    using PtrFaceType = std::shared_ptr<FaceType>;

    HexaCell() : Cell<Data>(){
    }

    void create(const PtrVertexType& v0, const PtrVertexType& v1, const PtrVertexType& v2, const PtrVertexType& v3,
             const PtrVertexType& v4, const PtrVertexType& v5, const PtrVertexType& v6, const PtrVertexType& v7) {

        auto cell = this->shared_from_this();

        auto face = std::make_shared<QuadFace<Data>>();
        face->create(v0, v3, v2, v1);
        face->addParent(cell);
        this->addChild(face);

        face = std::make_shared<QuadFace<Data>>();
        face->create(v4, v5, v6, v7);
        face->addParent(cell);
        this->addChild(face);

        face = std::make_shared<QuadFace<Data>>();
        face->create(v0, v1, v5, v4);
        face->addParent(cell);
        this->addChild(face);

        face = std::make_shared<QuadFace<Data>>();
        face->create(v1, v2, v6, v5);
        face->addParent(cell);
        this->addChild(face);

        face = std::make_shared<QuadFace<Data>>();
        face->create(v2, v3, v7, v6);
        face->addParent(cell);
        this->addChild(face);

        face = std::make_shared<QuadFace<Data>>();
        face->create(v0, v4, v7, v3);
        face->addParent(cell);
        this->addChild(face);

    }

};


#endif //ARTPDE_GEO_DATA_FACTORY_HPP
