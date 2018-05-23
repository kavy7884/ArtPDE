//
// Created by Chingkai Chou on 5/17/18.
//

#ifndef ARTPDE_GEO_TREE_HPP
#define ARTPDE_GEO_TREE_HPP

#include <memory>
#include <vector>
#include <list>

enum class TreeType{TreeHead, TreeConnect, TreeEnd};

template <typename SelfType>
class GeoTree{
public:
    GeoTree(TreeType tree_type) : tree_type(tree_type) {}

    TreeType getTree_type() const {
        return tree_type;
    }

    void setLinked_to(const std::shared_ptr<SelfType> &linked_to) {
        this->linked_to = linked_to;
    }

    virtual void link_self() = 0;

    const std::shared_ptr<SelfType> &getLinked_to() {
        if(linked_to == nullptr) link_self();
        return linked_to;
    }

protected:
    TreeType tree_type;
    std::shared_ptr<SelfType> linked_to{nullptr};
};

template <typename ParentType>
class GeoTree_Parent{
public:
    struct type{
        using PtrParentType = std::shared_ptr<ParentType>;
        using ListPtrParentType = std::list<PtrParentType>;
    };

    GeoTree_Parent() {}

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

    void mergeParents(GeoTree_Parent& rhs){
        for(auto &ptr_parent: rhs.c_getList_ptr_parents()){
            list_ptr_parents.push_back(ptr_parent);
        }
        rhs.clearParent();
    }

protected:
    typename
    type::ListPtrParentType list_ptr_parents;

};

template <typename ChildType>
class GeoTree_Child{
public:
    struct type{
        using PtrChildType = std::shared_ptr<ChildType>;
        using VecPtrChildType = std::vector<PtrChildType>;
    };

    GeoTree_Child() {}

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

template <typename BaseLayerType, typename MergeLayerType>
class GeoTree_LayerMerge{
public:
    using PtrBaseLayerType = std::shared_ptr<BaseLayerType>;
    using VecPtrBaseLayerType = std::vector<PtrBaseLayerType>;
    using PtrMergeLayerType = std::shared_ptr<MergeLayerType>;
    using VecPtrMergeLayerType = std::vector<PtrMergeLayerType>;

    GeoTree_LayerMerge(VecPtrBaseLayerType &vec_ptr_base_layer_seed) : vec_ptr_base_layer_seed(
            vec_ptr_base_layer_seed) {}

    VecPtrMergeLayerType merge(){
        bool merged = false;

        for (size_t i = 0; i < vec_ptr_base_layer_seed.size(); ++i) {
            auto & list_ptr_parents = this->vec_ptr_base_layer_seed[i]->getLinked_to()->getList_ptr_parents();

            //std::cout << list_ptr_parents.size() << std::endl;

            auto it_master = list_ptr_parents.begin();

            while( it_master != list_ptr_parents.end()){

                auto it_slave = it_master;
                ++it_slave;
                merged = false;
                while( it_slave != list_ptr_parents.end()) {
                    if((*(*it_master)) == (*(*it_slave))) {
                        std::cout << "Merge" << std::endl;
                        merged = true;
                        (*it_slave)->setLinked_to((*it_master));
                        (*it_master)->mergeParents(*(*it_slave));
                        for (auto &ptr_child: (*it_slave)->getVec_ptr_childs()) {
                            ptr_child->getLinked_to()->eraseParent((*it_slave));
                        }
                        it_slave = it_master;
                    }
                    ++it_slave;
                }
                if(merged){
                    vec_ptr_merged_layer.push_back((*it_master));
                }

                ++it_master;
            }
        }
        return vec_ptr_merged_layer;
    }

private:
    VecPtrBaseLayerType &vec_ptr_base_layer_seed;
    VecPtrMergeLayerType vec_ptr_merged_layer;
};


#endif //ARTPDE_GEO_TREE_HPP
