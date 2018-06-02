//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_DATA_HPP
#define ARTPDE_KAVY_GEOMETRY_DATA_HPP

#include "../GeometricTree/geometric_tree_units.hpp"

namespace art_pde{ namespace geometry {
        namespace mesh_type{
            using namespace geometric_tree;

            template <size_t Dimension, typename CalculatedPointType>
            class Data{
            public:
                struct type{
                    using VertexType = Vertex<CalculatedPointType>;
                    using VecPtrVertex = std::vector<std::shared_ptr<VertexType>>;
                    using EdgeType = Edge<CalculatedPointType>;
                    using VecPtrEdge = std::vector<std::shared_ptr<EdgeType>>;
                    using FaceType = Face<CalculatedPointType>;
                    using VecPtrFace = std::vector<std::shared_ptr<FaceType>>;
                    using CellType = Cell<CalculatedPointType>;
                    using VecPtrCell = std::vector<std::shared_ptr<CellType>>;
                };

                Data() {
                    std::cout << "Constructor" << std::endl;
                }

                const size_t getTotalNum_Vertex() const{
                    return this->vec_ptr_vertex.size();
                }

                const size_t getTotalNum_Cell() const{
                    return this->vec_ptr_cell.size();
                }

                const typename type::VecPtrVertex &getTotal_PtrVertex() const{
                    return this->vec_ptr_vertex;
                }

                const typename type::VecPtrCell &getTotal_PtrCell() const{
                    return this->vec_ptr_cell;
                }

            private:
                typename type::VecPtrVertex vec_ptr_vertex;
                typename type::VecPtrCell vec_ptr_cell;
            };
        }
}}

#endif //ARTPDE_KAVY_GEOMETRY_DATA_HPP
