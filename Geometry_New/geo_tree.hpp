//
// Created by Chingkai Chou on 5/17/18.
//

#ifndef ARTPDE_GEO_TREE_HPP
#define ARTPDE_GEO_TREE_HPP

#include <memory>
#include <vector>
#include <list>

class GeoTree{
public:
    enum class TreeType{TreeHead, TreeConnect, TreeEnd};
    GeoTree(TreeType tree_type) : tree_type(tree_type) {}

    TreeType getTree_type() const {
        return tree_type;
    }

    bool isLinked() const {
        return linked;
    }

    void setLinked_to(const std::shared_ptr <GeoTree> &linked_to) {
        this->linked = true;
        GeoTree::linked_to = linked_to;
    }

protected:
    TreeType tree_type;
    bool linked{false};
    std::shared_ptr<GeoTree> linked_to{nullptr};
};

template <typename ParentType>
class GeoTree_Parent{
public:
    struct type{
        using PtrParentType = std::shared_ptr<ParentType>;
        using VecPtrParentType = std::vector<PtrParentType>;
    };

    GeoTree_Parent() {}

    const typename type::VecPtrParentType &c_getVec_ptr_parents() const {
        return vec_ptr_parents;
    }

    const typename type::VecPtrParentType &getVec_ptr_parents() {
        return vec_ptr_parents;
    }

    void addParent(const typename type::PtrParentType &ptr_parent){
        vec_ptr_parents.push_back(ptr_parent);
    }

    void mergeParents(const GeoTree_Parent<ParentType> &rhs){
        for (auto $ptr_parent : rhs.c_getVec_ptr_parents()) {
            this->addParent($ptr_parent);
        }
    }

protected:
    typename
    type::VecPtrParentType vec_ptr_parents;

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

    const typename type::VecPtrChildType &getVec_ptr_childs() {
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

        for (size_t i = 0; i < vec_ptr_base_layer_seed.size(); ++i) {
            auto & vec_ptr_parents = this->vec_ptr_base_layer_seed[i]->getVec_ptr_parents();

            auto it_master = vec_ptr_parents.begin();

            while( it_master != vec_ptr_parents.end()){
                if((*it_master)->isLinked()){
                    ++it_master;
                    continue;
                }
                auto it_slave = it_master;
                ++it_slave;
                while( it_slave != vec_ptr_parents.end()) {
                    if((*it_slave)->isLinked()){
                        ++it_slave;
                        continue;
                    }
                    if((*(*it_master)) == (*(*it_slave))) {
                        std::cout << "Merge" << std::endl;
                        (*it_slave)->setLinked_to((*it_master));
                        (*it_master)->mergeParents(*(*it_slave));
                        vec_ptr_merged_layer.push_back((*it_master));
                    }
                    ++it_slave;
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
