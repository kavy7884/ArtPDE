//
// Created by Chingkai Chou on 5/31/18.
//

#ifndef ARTPDE_KAVY_GEOMETRIC_TREE_UNITS_HPP
#define ARTPDE_KAVY_GEOMETRIC_TREE_UNITS_HPP

#include <set>
#include <vector>
#include "geometric_tree_components.hpp"

namespace art_pde{ namespace geometry {
    namespace geometric_tree {

        // -------- Tree unit pre-define <Start> -----------
        template <typename DataType> class Cell;
        template <typename DataType> class Face;
        template <typename DataType> class Edge;
        template <typename DataType> class Vertex;
        // -------- Tree unit pre-define < End > -----------

        // -------- Vertex <Start> -----------
        template <typename DataType>
        class Vertex :
                public GeometricTree<Vertex<DataType>>,
                public GeometricTreeParent<Edge<DataType>>,
                public std::enable_shared_from_this<Vertex<DataType>>{
        public:
            struct type{
                using PtrDataType = std::shared_ptr<DataType>;
                using PtrEdgeType = std::shared_ptr<Edge<DataType>>;
                using ListPtrEdgeType = std::list<PtrEdgeType>;
            };

            Vertex() :
                    GeometricTree<Vertex<DataType>>(GeometricType::Point),
                    GeometricTreeParent<Edge<DataType>>(){}

            Vertex(const std::initializer_list<double> &input_list) :
                    GeometricTree<Vertex<DataType>>(GeometricType::Point),
                    GeometricTreeParent<Edge<DataType>>(){
                this->ptr_data = std::make_shared<DataType>(input_list);
            }

            Vertex(const typename type::PtrDataType &ptr_data) :
                    GeometricTree<Vertex<DataType>>(GeometricType::Point),
                    GeometricTreeParent<Edge<DataType>>(){
                this->setPtr_data(ptr_data);
            }

            void link_self() override {
                this->linked_to = this->shared_from_this();
            }

            const typename type::ListPtrEdgeType& c_getConnected_Edge(){
                return this->getLinked_to()->c_getList_ptr_parents();
            }

            const typename type::PtrDataType &getPtr_data() const {
                return ptr_data;
            }

            void setPtr_data(const typename type::PtrDataType &ptr_data) {
                Vertex::ptr_data = ptr_data;
            }

            friend std::ostream &operator<<(std::ostream &os, const Vertex<DataType> &vertex) {
                os << *vertex.getPtr_data() ;
                return os;
            }

        private:
            typename type::PtrDataType ptr_data{nullptr};
        };
        // -------- Vertex < End > -----------

        // -------- Edge <Start> -----------
        template <typename DataType>
        class Edge:
                public GeometricTree<Edge<DataType>>,
                public GeometricTreeParent<Face<DataType>>,
                public GeometricTreeChild<Vertex<DataType>>,
                public std::enable_shared_from_this<Edge<DataType>>{
        public:
            struct type{
                using PtrVertexType = std::shared_ptr<Vertex<DataType>>;
                using VecPtrVertexType = std::vector<PtrVertexType>;
                using PtrFaceType = std::shared_ptr<Face<DataType>>;
                using ListPtrFaceType = std::list<PtrFaceType>;
            };

            Edge() :
                    GeometricTree<Edge<DataType>>(GeometricType::Line),
                    GeometricTreeParent<Face<DataType>>(),
                    GeometricTreeChild<Vertex<DataType>>() {
            }

            void link_self() override {
                this->linked_to = this->shared_from_this();
            }

            const typename type::VecPtrVertexType& c_getConnected_Vertex() const{
                return this->c_getVec_ptr_childs();
            }

            const typename type::ListPtrFaceType& c_getConnected_Face() {
                return this->getLinked_to()->c_getList_ptr_parents();
            }

            friend std::ostream &operator<<(std::ostream &os, const Edge<DataType> &edge) {
                auto it = edge.c_getConnected_Vertex().cbegin();
                os << "(" <<**it << " -> " << **(++it) << ")" ;
                return os;
            }

            bool operator==(Edge &rhs){
                if(this->c_getConnected_Vertex().size() != rhs.c_getConnected_Vertex().size())  return false;
                else{
                    std::set<typename type::PtrVertexType> set_self(this->c_getConnected_Vertex().cbegin(), this->c_getConnected_Vertex().cend());
                    std::set<typename type::PtrVertexType> set_rhs(rhs.c_getConnected_Vertex().cbegin(), rhs.c_getConnected_Vertex().cend());
                    if ( std::equal(set_self.cbegin(), set_self.cend(), set_rhs.cbegin() ))return true;
                    else return false;
                }
            }

            bool operator!=(const Edge &rhs) {
                return !(rhs == *this);
            }

        };
        // -------- Edge < End > -----------

        // -------- Face <Start> -----------
        template <typename DataType>
        class Face:
                public GeometricTree<Face<DataType>>,
                public GeometricTreeParent<Cell<DataType>>,
                public GeometricTreeChild<Edge<DataType>>,
                public std::enable_shared_from_this<Face<DataType>>{
        public:
            struct type{
                using PtrVertexType = std::shared_ptr<Vertex<DataType>>;
                using VecPtrVertexType = std::vector<PtrVertexType>;
                using PtrEdgeType = std::shared_ptr<Edge<DataType>>;
                using VecPtrEdgeType = std::vector<PtrEdgeType>;
                using ListPtrCellType = std::list<std::shared_ptr<Cell<DataType>>>;
            };

            Face() :
                    GeometricTree<Face<DataType>>(GeometricType::None),
                    GeometricTreeParent<Cell<DataType>>(),
                    GeometricTreeChild<Edge<DataType>>() {}

            Face(GeometricType geometric_type) :
                    GeometricTree<Face<DataType>>(geometric_type),
                    GeometricTreeParent<Cell<DataType>>(),
                    GeometricTreeChild<Edge<DataType>>() {}

            void link_self() override {
                this->linked_to = this->shared_from_this();
            }

            const typename type::VecPtrEdgeType& c_getConnected_Edge() const{
                return this->c_getVec_ptr_childs();
            }

            const typename type::ListPtrCellType& c_getConnected_Cell() {
                return this->getLinked_to()->c_getList_ptr_parents();
            }

            const typename type::VecPtrVertexType c_getConnected_Vertex() const{
                typename type::VecPtrVertexType re_vec_ptr_vertex;
                for(auto &ptr_edge: this->c_getConnected_Edge()){
                    re_vec_ptr_vertex.push_back(*ptr_edge->c_getConnected_Vertex().cbegin());
                }
                return re_vec_ptr_vertex;
            }

            friend std::ostream &operator<<(std::ostream &os, const Face<DataType> &face) {
                os << "[";
                auto it_ptr_edge = face.c_getConnected_Edge().cbegin();
                while(it_ptr_edge !=  (face.c_getConnected_Edge().cend() - 1)){
                    os << **it_ptr_edge << " , ";
                    ++it_ptr_edge;
                }
                os << **it_ptr_edge << "]";
                return os;
            }

            bool operator==(Face &rhs){
                if(this->c_getConnected_Edge().size() != rhs.c_getConnected_Edge().size())  return false;
                else{
                    std::set<typename type::PtrEdgeType> set_self, set_rhs;
                    for (auto &ptr_edge: this->c_getConnected_Edge()) {
                        set_self.insert(ptr_edge->getLinked_to());
                    }
                    for (auto &ptr_edge: rhs.c_getConnected_Edge()) {
                        set_rhs.insert(ptr_edge->getLinked_to());
                    }
                    if ( std::equal(set_self.cbegin(), set_self.cend(), set_rhs.cbegin() ))return true;
                    else return false;
                }
                return true;
            }

            bool operator!=(const Face &rhs) {
                return !(rhs == *this);
            }
        };
        // -------- Face <End> -----------

        // -------- Cell <Start> -----------
        template <typename DataType>
        class Cell:
                public GeometricTree<Cell<DataType>>,
                public GeometricTreeChild<Face<DataType>>,
                public std::enable_shared_from_this<Cell<DataType>>{
        public:
            struct type{
                using PtrVertexType = std::shared_ptr<Vertex<DataType>>;
                using VecPtrVertexType = std::vector<PtrVertexType>;
                using PtrFaceType = std::shared_ptr<Face<DataType>>;
                using VecPtrFaceType = std::vector<PtrFaceType>;
            };

            Cell() :
                    GeometricTree<Cell<DataType>>(GeometricType::None),
                    GeometricTreeChild<Face<DataType>>() {}

            Cell(GeometricType geometric_type) :
                    GeometricTree<Cell<DataType>>(geometric_type),
                    GeometricTreeChild<Face<DataType>>() {}

            void link_self() override {
                this->linked_to = this->shared_from_this();
            }

            const typename type::VecPtrFaceType& c_getConnected_Face() const{
                return this->c_getVec_ptr_childs();
            }

            const typename type::VecPtrVertexType c_getConnected_Vertex() const{
                typename type::VecPtrVertexType re_vec_ptr_vertex;

                // Hexa
                if(this->getGeometric_Type() == GeometricType::Hexahedron){
                    auto &face_1_vec_ptr_vertex = this->c_getConnected_Face()[0]->c_getConnected_Vertex();
                    auto &face_2_vec_ptr_vertex = this->c_getConnected_Face()[1]->c_getConnected_Vertex();
                    re_vec_ptr_vertex.insert(re_vec_ptr_vertex.end(), face_1_vec_ptr_vertex.begin(), face_1_vec_ptr_vertex.end());
                    re_vec_ptr_vertex.insert(re_vec_ptr_vertex.end(), face_2_vec_ptr_vertex.begin(), face_2_vec_ptr_vertex.end());
                }
                return re_vec_ptr_vertex;
            }

            friend std::ostream &operator<<(std::ostream &os, const Cell<DataType> &cell) {
                os << "{";
                auto it_ptr_face = cell.c_getConnected_Face().cbegin();
                while(it_ptr_face !=  (cell.c_getConnected_Face().cend() - 1)){
                    os << **it_ptr_face << " , ";
                    ++it_ptr_face;
                }
                os << **it_ptr_face << "}";
                return os;
            }


        };
        // -------- Cell < End > -----------
    }
}}

#endif //ARTPDE_KAVY_GEOMETRIC_TREE_UNITS_HPP
