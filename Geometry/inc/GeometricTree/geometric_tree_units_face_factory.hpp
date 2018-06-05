//
// Created by Chingkai Chou on 6/4/18.
//

#ifndef ARTPDE_KAVY_GEOMETRIC_TREE_UNITS_FACTORY_HPP
#define ARTPDE_KAVY_GEOMETRIC_TREE_UNITS_FACTORY_HPP

#include "geometric_tree_units.hpp"

namespace art_pde{ namespace geometry {
        namespace geometric_tree {

            template <typename DataType>
            class Quadrilateral_Face: public Face<DataType>{
            public:
                using PtrVertexType = typename Face<DataType>::type::PtrVertexType;
                using EdgeType = Edge<DataType>;
                using PtrEdgeType = typename Face<DataType>::type::PtrEdgeType;

                Quadrilateral_Face(): Face<DataType>(GeometricType::Quadrilateral){}
                void create(const PtrVertexType& v1, const PtrVertexType& v2, const PtrVertexType& v3, const PtrVertexType& v4){
                    auto face = this->shared_from_this();

                    PtrEdgeType edge;

                    edge = genEdge(v1, v2);
                    edge->addParent(face);
                    face->addChild(edge);

                    edge = genEdge(v2, v3);
                    edge->addParent(face);
                    face->addChild(edge);

                    edge = genEdge(v3, v4);
                    edge->addParent(face);
                    face->addChild(edge);

                    edge = genEdge(v4, v1);
                    edge->addParent(face);
                    face->addChild(edge);
                }

                void create(const PtrEdgeType& e1, const PtrEdgeType& e2, const PtrEdgeType& e3, const PtrEdgeType& e4){
                    auto face = this->shared_from_this();

                    e1->addParent(face);
                    face->addChild(e1);

                    e2->addParent(face);
                    face->addChild(e2);

                    e3->addParent(face);
                    face->addChild(e3);

                    e4->addParent(face);
                    face->addChild(e4);
                }

            private:
                PtrEdgeType genEdge(const PtrVertexType& v1, const PtrVertexType& v2){
                    PtrEdgeType reEdge = std::make_shared<EdgeType>();
                    v1->addParent(reEdge);
                    v2->addParent(reEdge);
                    reEdge->addChild(v1);
                    reEdge->addChild(v2);
                    return reEdge;
                }

            };

            template <typename DataType>
            class Face_Factory{
            public:
                using VertexType = Vertex<DataType>;
                using PtrVertexType = std::shared_ptr<VertexType>;
                using VecPtrVertexType = std::vector<PtrVertexType>;
                using PtrFaceType = std::shared_ptr<Face<DataType>>;

                Face_Factory(){}
                void addVertex(PtrVertexType & ptr_vertex){ vec_ptr_vertex.push_back(ptr_vertex); }
                void clearVertex(){ vec_ptr_vertex.clear(); }

                PtrFaceType create(){

                    if(this->vec_ptr_vertex.size() == 3) {
                        // TODO - Triangle
                        auto ptr_face = std::make_shared<Face<DataType>>();
                        return ptr_face;
                    }
                    else if (this->vec_ptr_vertex.size() == 4){
                        auto ptr_quad_face = std::make_shared<Quadrilateral_Face<DataType>>();
                        ptr_quad_face->create(this->vec_ptr_vertex[0], this->vec_ptr_vertex[1], this->vec_ptr_vertex[2], this->vec_ptr_vertex[3]);

                        return ptr_quad_face;
                    }
                    else{
                        auto ptr_face = std::make_shared<Face<DataType>>();
                        return ptr_face;
                    }
                }

            private:
                VecPtrVertexType vec_ptr_vertex;
            };

        }
}}


#endif //ARTPDE_KAVY_GEOMETRIC_TREE_UNITS_FACTORY_HPP
