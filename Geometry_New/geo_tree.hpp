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

    bool isMerged() const {
        return merged;
    }

    void setMerged(bool merged) {
        GeoTree::merged = merged;
    }

protected:
    TreeType tree_type;
    bool merged{false};
};

template <typename ParentType>
class GeoTree_Parent{
public:
    struct type{
        using PtrParentType = std::shared_ptr<ParentType>;
        using ListPtrParentType = std::list<PtrParentType>;
        using PtrListPtrParentType = std::shared_ptr<ListPtrParentType>;
    };

    GeoTree_Parent() {
        ptr_list_ptr_parents = std::make_shared<typename type::ListPtrParentType>();
    }

    const typename type::PtrListPtrParentType &c_getPtr_list_ptr_parents() const {
        return ptr_list_ptr_parents;
    }

    const typename type::PtrListPtrParentType &getPtr_list_ptr_parents() {
        return ptr_list_ptr_parents;
    }

    void addParent(const typename type::PtrParentType &ptr_parent){
        this->ptr_list_ptr_parents->push_back(ptr_parent);
    }

    void eraseParent(const typename type::PtrParentType &ptr_parent){
        auto it = this->ptr_list_ptr_parents->begin();
        while (it != this->ptr_list_ptr_parents->end()){
            if(*it == ptr_parent){
                it = this->ptr_list_ptr_parents->erase(it);
            }else{
                ++it;
            }
        }
    }

    void replaceParent(const typename type::PtrParentType &ptr_parent_old, const typename type::PtrParentType &ptr_parent_new){
        auto it = this->ptr_list_ptr_parents->begin();
        while (it != this->ptr_list_ptr_parents->end()) {
            if(*it == ptr_parent_old){
                *it = ptr_parent_new;
            }
            ++it;
        }
    }

    void mergeParents(const GeoTree_Parent<ParentType> &rhs){
        auto rhs_ptr_list_ptr_parents = rhs.c_getPtr_list_ptr_parents();
        auto it = rhs_ptr_list_ptr_parents->begin();
        while (it != rhs_ptr_list_ptr_parents->end()){
            this->ptr_list_ptr_parents->push_back(*it);
            ++it;
        }
    }

protected:
    typename
    type::PtrListPtrParentType ptr_list_ptr_parents{nullptr};
};

template <typename ChildType>
class GeoTree_Child{
public:
    struct type{
        using PtrChildType = std::shared_ptr<ChildType>;
        using ListPtrChildType = std::list<PtrChildType>;
        using PtrListPtrChildType = std::shared_ptr<ListPtrChildType>;
    };

    GeoTree_Child() {
        ptr_list_ptr_childs = std::make_shared<typename type::ListPtrChildType>();
    }

    const typename type::PtrListPtrChildType &c_getPtr_list_ptr_childs() const {
        return ptr_list_ptr_childs;
    }

    const typename type::PtrListPtrChildType &getPtr_list_ptr_childs() {
        return ptr_list_ptr_childs;
    }

    void addChild(const typename type::PtrChildType &ptr_child){
        this->ptr_list_ptr_childs->push_back(ptr_child);
    }

    void replaceChild(const typename type::PtrChildType &ptr_child_old, const typename type::PtrChildType &ptr_child_new){
        auto it = this->ptr_list_ptr_childs->begin();
        while (it != this->ptr_list_ptr_childs->end()) {
            if(*it == ptr_child_old){
                *it = ptr_child_new;
            }
            ++it;
        }
    }

    void mergeChilds(const GeoTree_Child<ChildType> &rhs){
        auto rhs_ptr_list_ptr_childs = rhs.c_getPtr_list_ptr_childs();
        auto it = rhs_ptr_list_ptr_childs->begin();
        while (it != rhs_ptr_list_ptr_childs->end()){
            this->ptr_list_ptr_childs->push_back(*it);
            ++it;
        }
    }

    size_t getNum_childs_per_group() const {
        return num_childs_per_group;
    }

    void setNum_childs_per_group(size_t num_childs_per_group) {
        GeoTree_Child::num_childs_per_group = num_childs_per_group;
    }

protected:
    typename
    type::PtrListPtrChildType ptr_list_ptr_childs{nullptr};

    size_t num_childs_per_group{0};
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
            auto & ptr_list_ptr_parents = this->vec_ptr_base_layer_seed[i]->getPtr_list_ptr_parents();

            auto it_master = ptr_list_ptr_parents->begin();

            while( it_master != ptr_list_ptr_parents->end()){
                if((*it_master)->isMerged()){
                    ++it_master;
                    continue;
                }
                auto it_slave = it_master;
                ++it_slave;
                while( it_slave != ptr_list_ptr_parents->end()) {
                    if((*(*it_master)) == (*(*it_slave))) {
                        std::cout << "Merge" << std::endl;

                        (*it_master)->mergeParents(*(*it_slave));
                        (*it_master)->mergeChilds(*(*it_slave));

                        for(auto & replace_parent : *(*it_slave)->getPtr_list_ptr_parents()){
                            replace_parent->replaceChild((*it_slave), (*it_master));
                        }

                        for(auto & replace_child : *(*it_slave)->getPtr_list_ptr_childs()){
                            replace_child->eraseParent((*it_slave));
                        }

                        if(!(*it_master)->isMerged()){
                            (*it_master)->setMerged(true);
                            vec_ptr_merged_layer.push_back(*it_master);
                        }
                        it_slave = it_master;
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
