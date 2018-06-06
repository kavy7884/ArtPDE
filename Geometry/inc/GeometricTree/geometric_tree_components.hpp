//
// Created by Chingkai Chou on 5/31/18.
//

#ifndef ARTPDE_KAVY_GEOMETRIC_TREE_COMPONENTS_HPP
#define ARTPDE_KAVY_GEOMETRIC_TREE_COMPONENTS_HPP

#include <iostream>
#include <memory>
#include <vector>
#include <list>

namespace art_pde{ namespace geometry {
    namespace geometric_tree{

        // -------- GeometricType <Start> -----------
        //enum class TreeType{TreeHead, TreeConnect, TreeEnd};
        enum class GeometricType{ None, Point, Line, Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid, Polyhedron };

        std::ostream &operator<<(std::ostream &os, const GeometricType &geometric_type) {
            switch(geometric_type) {
                case GeometricType::None    : os << "None"; break;
                case GeometricType::Point   : os << "Point"; break;
                case GeometricType::Line    : os << "Line"; break;
                case GeometricType::Triangle    : os << "Triangle"; break;
                case GeometricType::Quadrilateral    : os << "Quadrilateral"; break;
                case GeometricType::Tetrahedron    : os << "Tetrahedron"; break;
                case GeometricType::Hexahedron    : os << "Hexahedron"; break;
                case GeometricType::Prism    : os << "Prism"; break;
                case GeometricType::Pyramid    : os << "Pyramid"; break;
                case GeometricType::Polyhedron    : os << "Polyhedron"; break;
                default : os.setstate(std::ios_base::failbit);
            }
            return os;
        }

        // -------- TreeType < End > -----------

        // -------- GeometricType <Start> -----------
        template <typename SelfType>
        class GeometricTree{
        public:
            GeometricTree(GeometricType geometric_type) : geometric_type(geometric_type) {}

            GeometricType getGeometric_Type() const {
                return this->geometric_type;
            }

            void setGeometric_Type(GeometricType geometric_type){
                this->geometric_type = geometric_type;
            }

            void setLinked_to(const std::shared_ptr<SelfType> &linked_to) {
                this->linked_to = linked_to;
            }

            virtual void link_self() = 0;

            const std::shared_ptr<SelfType> &getLinked_to() {
                if(linked_to == nullptr) link_self();
                return linked_to;
            }

            const bool &getIsMerged() const { return this->is_merged;}
            void setIsMerged(bool merge_condition){this->is_merged = merge_condition;}

        protected:
            GeometricType geometric_type;
            bool is_merged{false};
            std::shared_ptr<SelfType> linked_to{nullptr};
        };
        // -------- GeometricTree < End > -----------

        // -------- GeometricTreeParent <Start> -----------
        template <typename ParentType>
        class GeometricTreeParent{
        public:
            struct type{
                using PtrParentType = std::shared_ptr<ParentType>;
                using ListPtrParentType = std::list<PtrParentType>;
            };

            GeometricTreeParent() {}

            void addParent(const typename type::PtrParentType &ptr_parent){
                list_ptr_parents.push_back(ptr_parent);
            }

            void clearParent(){
                list_ptr_parents.clear();
            }

            void eraseParent(const typename type::PtrParentType &ptr_parent){
                auto it = list_ptr_parents.begin();
                while (it != list_ptr_parents.end()){
                    if(*it == ptr_parent){
                        it = list_ptr_parents.erase(it);
                    }else{
                        ++it;
                    }
                }
            }

            const typename type::ListPtrParentType &c_getList_ptr_parents() const {
                return list_ptr_parents;
            }

            typename type::ListPtrParentType &getList_ptr_parents() {
                return list_ptr_parents;
            }

            void mergeParents(GeometricTreeParent& rhs){
                for(auto &ptr_parent: rhs.c_getList_ptr_parents()){
                    list_ptr_parents.push_back(ptr_parent);
                }
                rhs.clearParent();
            }

        protected:
            typename
            type::ListPtrParentType list_ptr_parents;
        };
        // -------- GeometricTreeParent < End > -----------

        // -------- GeometricTreeChild <Start> -----------
        template <typename ChildType>
        class GeometricTreeChild{
        public:
            struct type{
                using PtrChildType = std::shared_ptr<ChildType>;
                using VecPtrChildType = std::vector<PtrChildType>;
            };

            GeometricTreeChild() {}

            const typename type::VecPtrChildType &c_getVec_ptr_childs() const {
                return vec_ptr_childs;
            }

            typename type::VecPtrChildType &getVec_ptr_childs() {
                return vec_ptr_childs;
            }

            void addChild(const typename type::PtrChildType &ptr_child){
                vec_ptr_childs.push_back(ptr_child);
            }

        protected:
            typename
            type::VecPtrChildType vec_ptr_childs;
        };
        // -------- GeometricTreeChild <End> -----------

        // -------- GeometricTreeLayerMerge <Start> -----------
        template <typename MergeLayerType, typename BaseLayerType>
        class GeometricTreeLayerMerge{
        public:
            struct type{
                using PtrBaseLayerType = std::shared_ptr<BaseLayerType>;
                using VecPtrBaseLayerType = std::vector<PtrBaseLayerType>;
                using PtrMergeLayerType = std::shared_ptr<MergeLayerType>;
                using VecPtrMergeLayerType = std::vector<PtrMergeLayerType>;
            };

            GeometricTreeLayerMerge(typename type::VecPtrBaseLayerType &vec_ptr_base_layer_seed)
                    : vec_ptr_base_layer_seed(vec_ptr_base_layer_seed) {}

            bool merge(typename type::VecPtrMergeLayerType &vec_ptr_merged_layer){
                vec_ptr_merged_layer.clear();

                for (size_t i = 0; i < vec_ptr_base_layer_seed.size(); ++i) {
                    auto & list_ptr_parents = this->vec_ptr_base_layer_seed[i]->getLinked_to()->getList_ptr_parents();

                    auto it_master = list_ptr_parents.begin();

                    while( it_master != list_ptr_parents.end()){
                        if((*it_master)->getIsMerged()){
                            ++it_master;
                            continue;
                        }

                        auto it_slave = it_master;
                        ++it_slave;
                        while( it_slave != list_ptr_parents.end()) {
                            if((*(*it_master)) == (*(*it_slave))) {
                                //std::cout << "Merge" << std::endl;
                                (*it_slave)->setLinked_to((*it_master));
                                (*it_master)->mergeParents(*(*it_slave));
                                for (auto &ptr_child: (*it_slave)->getVec_ptr_childs()) {
                                    ptr_child->getLinked_to()->eraseParent((*it_slave));
                                }
                                it_slave = it_master;
                            }
                            ++it_slave;
                        }
                        (*it_master)->setIsMerged(true);
                        vec_ptr_merged_layer.push_back((*it_master));
                        ++it_master;
                    }
                }
                return true;
            }

        private:
            typename type::VecPtrBaseLayerType &vec_ptr_base_layer_seed;
            typename type::VecPtrMergeLayerType vec_ptr_merged_layer;
        };
        // -------- GeometricTreeLayerMerge < End > -----------
    }
}}

#endif //ARTPDE_KAVY_GEOMETRIC_TREE_COMPONENTS_HPP
