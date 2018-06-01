//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_DATA_HPP
#define ARTPDE_KAVY_GEOMETRY_DATA_HPP

#include "inc/GeometricTree/geometric_tree_units.hpp"

namespace art_pde{ namespace geometry {
    namespace mesh_type{
        using namespace geometric_tree;

        template <size_t Dimension, template<size_t> class CalculatedPointType>
        class Data{
        public:
            struct type{
                using VertexType = Vertex<CalculatedPointType<Dimension>>;
                using VecPtrVertex = std::vector<std::shared_ptr<VertexType>>;
                using EdgeType = Edge<CalculatedPointType<Dimension>>;
                using VecPtrEdge = std::vector<std::shared_ptr<EdgeType>>;
                using FaceType = Face<CalculatedPointType<Dimension>>;
                using VecPtrFace = std::vector<std::shared_ptr<FaceType>>;
                using CellType = Cell<CalculatedPointType<Dimension>>;
                using VecPtrCell = std::vector<std::shared_ptr<CellType>>;
            };
            Data() {
                std::cout << "Constructor" << std::endl;
            }

            //geometric_tree::Vertex<CalculatedPointType<Dimension>> vertex{1,2,3};

        protected:
            typename type::VecPtrVertex vec_ptr_vertex;
            typename type::VecPtrCell vec_ptr_cell;
        };
    }
}}

#endif //ARTPDE_KAVY_GEOMETRY_DATA_HPP
