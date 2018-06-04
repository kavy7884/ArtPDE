//
// Created by Chingkai Chou on 6/3/18.
//

#ifndef ARTPDE_KAVY_GEOMETRIC_DATA_MESH_TYPE_HPP
#define ARTPDE_KAVY_GEOMETRIC_DATA_MESH_TYPE_HPP

#include <memory>
#include "../GeometricTree/geometric_tree_units.hpp"
namespace art_pde { namespace geometry { namespace mesh_type {
            namespace geometric_data{

                using namespace geometric_tree;

                template <typename PointType>
                class MeshTypeData_Vertex{
                public:
                    struct type {
                        using VertexType = Vertex<PointType>;
                        using PtrVertexType = std::shared_ptr<VertexType>;
                        using VecPtrVertexType = std::vector<PtrVertexType>;
                    };
                    MeshTypeData_Vertex(){
                        std::cout << "MeshTypeData_Vertex" << std::endl;
                    }

                    const size_t getTotalNum_Vertex() const {
                        return this->vec_ptr_vertex.size();
                    }

                    const typename type::VecPtrVertexType &c_getTotalVec_PtrVertex() const {
                        return this->vec_ptr_vertex;
                    }

                protected:
                    typename type::VecPtrVertexType &getTotalVec_PtrVertex() {
                        return this->vec_ptr_vertex;
                    }

                    typename type::VecPtrVertexType vec_ptr_vertex;
                };

                template <typename PointType>
                class MeshTypeData_Edge{
                public:
                    struct type {
                        using EdgeType = Edge<PointType>;
                        using PtrEdgeType = std::shared_ptr<EdgeType>;
                        using VecPtrEdgeType = std::vector<PtrEdgeType>;
                    };
                    MeshTypeData_Edge(){
                        std::cout << "MeshTypeData_Edge" << std::endl;
                    }

                    const size_t getTotalNum_Edge() const {
                        return this->vec_ptr_edge.size();
                    }

                    const typename type::VecPtrEdgeType &c_getTotalVec_PtrEdge() const {
                        return this->vec_ptr_edge;
                    }

                protected:
                    typename type::VecPtrEdgeType &getTotalVec_PtrEdge(){
                        return this->vec_ptr_edge;
                    }

                    typename type::VecPtrEdgeType vec_ptr_edge;
                };

                template <typename PointType>
                class MeshTypeData_Face{
                public:
                    struct type {
                        using FaceType = Face<PointType>;
                        using PtrFaceType = std::shared_ptr<FaceType>;
                        using VecPtrFaceType = std::vector<PtrFaceType>;
                    };
                    MeshTypeData_Face(){
                        std::cout << "MeshTypeData_Face" << std::endl;
                    }

                    const size_t getTotalNum_Face() const {
                        return this->vec_ptr_face.size();
                    }

                    const typename type::VecPtrFaceType &c_getTotalVec_PtrFace() const {
                        return this->vec_ptr_face;
                    }

                protected:
                    typename type::VecPtrFaceType &getTotalVec_PtrFace(){
                        return this->vec_ptr_face;
                    }

                    typename type::VecPtrFaceType vec_ptr_face;
                };

                template <typename PointType>
                class MeshTypeData_Cell{
                public:
                    struct type {
                        using CellType = Cell<PointType>;
                        using PtrCellType = std::shared_ptr<CellType>;
                        using VecPtrCellType = std::vector<PtrCellType>;
                    };
                    MeshTypeData_Cell(){
                        std::cout << "MeshTypeData_Cell" << std::endl;
                    }

                    const size_t getTotalNum_Cell() const {
                        return this->vec_ptr_cell.size();
                    }

                    const typename type::VecPtrCellType &c_getTotalVec_PtrCell() const {
                        return this->vec_ptr_cell;
                    }

                protected:
                    typename type::VecPtrCellType &getTotalVec_PtrCell() {
                        return this->vec_ptr_cell;
                    }

                    typename type::VecPtrCellType vec_ptr_cell;
                };

            }
}}}
#endif //ARTPDE_KAVY_GEOMETRIC_DATA_MESH_TYPE_HPP
