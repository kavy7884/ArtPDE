//
// Created by Chingkai Chou on 6/5/18.
//

#ifndef ARTPDE_KAVY_GEOMETRIC_TREE_UNITS_CELL_FACTORY_HPP
#define ARTPDE_KAVY_GEOMETRIC_TREE_UNITS_CELL_FACTORY_HPP

#include "memory"
#include "geometric_tree_units.hpp"
#include "geometric_tree_units_face_factory.hpp"

namespace art_pde{ namespace geometry {
        namespace geometric_tree {

            template <typename DataType>
            class Hexahedron_Cell: public Cell<DataType> {
            public:
                using PtrVertexType = typename Face<DataType>::type::PtrVertexType;
                using PtrEdgeType = typename Face<DataType>::type::PtrEdgeType;
                using VecPtrEdgeType = typename Face<DataType>::type::VecPtrEdgeType;
                using PtrFaceType = typename Cell<DataType>::type::PtrFaceType;
                using VecPtrFaceType = typename Cell<DataType>::type::VecPtrFaceType;

                Hexahedron_Cell(): Cell<DataType>(GeometricType::Hexahedron){}
                void create(const PtrVertexType& v1, const PtrVertexType& v2, const PtrVertexType& v3, const PtrVertexType& v4,
                            const PtrVertexType& v5, const PtrVertexType& v6, const PtrVertexType& v7, const PtrVertexType& v8){

                    auto cell = this->shared_from_this();
                    VecPtrEdgeType vec_ptr_edge(12);
                    VecPtrFaceType vec_ptr_face(6);

                    vec_ptr_edge[0] = genEdge(v1, v2);
                    vec_ptr_edge[1] = genEdge(v2, v3);
                    vec_ptr_edge[2] = genEdge(v3, v4);
                    vec_ptr_edge[3] = genEdge(v4, v1);
                    vec_ptr_edge[4] = genEdge(v5, v6);
                    vec_ptr_edge[5] = genEdge(v6, v7);
                    vec_ptr_edge[6] = genEdge(v7, v8);
                    vec_ptr_edge[7] = genEdge(v8, v5);
                    vec_ptr_edge[8] = genEdge(v1, v5);
                    vec_ptr_edge[9] = genEdge(v2, v6);
                    vec_ptr_edge[10] = genEdge(v3, v7);
                    vec_ptr_edge[11] = genEdge(v4, v8);


                    auto ptr_quad_face = std::make_shared<Quadrilateral_Face<DataType>>();
                    ptr_quad_face->create(vec_ptr_edge[0], vec_ptr_edge[1], vec_ptr_edge[2], vec_ptr_edge[3]);
                    vec_ptr_face[0] = ptr_quad_face;

                    ptr_quad_face = std::make_shared<Quadrilateral_Face<DataType>>();
                    ptr_quad_face->create(vec_ptr_edge[4], vec_ptr_edge[5], vec_ptr_edge[6], vec_ptr_edge[7]);
                    vec_ptr_face[1] = ptr_quad_face;

                    ptr_quad_face = std::make_shared<Quadrilateral_Face<DataType>>();
                    ptr_quad_face->create(vec_ptr_edge[0], vec_ptr_edge[9], vec_ptr_edge[4], vec_ptr_edge[8]);
                    vec_ptr_face[2] = ptr_quad_face;

                    ptr_quad_face = std::make_shared<Quadrilateral_Face<DataType>>();
                    ptr_quad_face->create(vec_ptr_edge[1], vec_ptr_edge[10], vec_ptr_edge[5], vec_ptr_edge[9]);
                    vec_ptr_face[3] = ptr_quad_face;

                    ptr_quad_face = std::make_shared<Quadrilateral_Face<DataType>>();
                    ptr_quad_face->create(vec_ptr_edge[2], vec_ptr_edge[10], vec_ptr_edge[6], vec_ptr_edge[11]);
                    vec_ptr_face[4] = ptr_quad_face;

                    ptr_quad_face = std::make_shared<Quadrilateral_Face<DataType>>();
                    ptr_quad_face->create(vec_ptr_edge[3], vec_ptr_edge[11], vec_ptr_edge[7], vec_ptr_edge[8]);
                    vec_ptr_face[5] = ptr_quad_face;

                    cell->addChild(vec_ptr_face[0]);
                    cell->addChild(vec_ptr_face[1]);
                    cell->addChild(vec_ptr_face[2]);
                    cell->addChild(vec_ptr_face[3]);
                    cell->addChild(vec_ptr_face[4]);
                    cell->addChild(vec_ptr_face[5]);
                }

            private:
                PtrEdgeType genEdge(const PtrVertexType& v1, const PtrVertexType& v2){
                    PtrEdgeType reEdge = std::make_shared<Edge<DataType>>();
                    v1->addParent(reEdge);
                    v2->addParent(reEdge);
                    reEdge->addChild(v1);
                    reEdge->addChild(v2);
                    return reEdge;
                }
            };

            template <typename DataType>
            class Cell_Factory{
            public:
                using VertexType = Vertex<DataType>;
                using PtrVertexType = std::shared_ptr<VertexType>;
                using VecPtrVertexType = std::vector<PtrVertexType>;
                using PtrCellType = std::shared_ptr<Cell<DataType>>;

                Cell_Factory(){}
                void addVertex(PtrVertexType & ptr_vertex){ vec_ptr_vertex.push_back(ptr_vertex); }
                void clearVertex(){ vec_ptr_vertex.clear(); }

                PtrCellType create(){

                    if(this->vec_ptr_vertex.size() == 4) {
                        // TODO - Tetra
                        auto ptr_cell = std::make_shared<Cell<DataType>>();
                        return ptr_cell;
                    }
                    else if (this->vec_ptr_vertex.size() == 8){
                        auto ptr_hexa_cell = std::make_shared<Hexahedron_Cell<DataType>>();
                        ptr_hexa_cell->create(this->vec_ptr_vertex[0], this->vec_ptr_vertex[1], this->vec_ptr_vertex[2], this->vec_ptr_vertex[3],
                                              this->vec_ptr_vertex[4], this->vec_ptr_vertex[5], this->vec_ptr_vertex[6], this->vec_ptr_vertex[7]);

                        return ptr_hexa_cell;
                    }
                    else{
                        auto ptr_cell = std::make_shared<Cell<DataType>>();
                        return ptr_cell;
                    }
                }

            private:
                VecPtrVertexType vec_ptr_vertex;
            };



        }
}}

#endif //ARTPDE_KAVY_GEOMETRIC_TREE_UNITS_CELL_FACTORY_HPP
