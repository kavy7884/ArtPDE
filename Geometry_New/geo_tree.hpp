//
// Created by Chingkai Chou on 5/17/18.
//

#ifndef ARTPDE_GEO_TREE_HPP
#define ARTPDE_GEO_TREE_HPP

#include <memory>
#include <vector>
#include <list>


template <typename ChildType>
class GeoHead{
public:
    using PtrChildType = std::shared_ptr<ChildType>;
    using ListPtrChildType = std::list<PtrChildType>;
    GeoHead() {}

    void addChild(const PtrChildType & ptr_child){
        list_ptr_childs.push_back(ptr_child);
    }

    ListPtrChildType &getChild(){
        return list_ptr_childs;
    }

    void replaceChild(PtrChildType &old_child, PtrChildType &new_child){
        for(auto & ptr_child :list_ptr_childs){
            if(ptr_child == old_child){
                ptr_child = new_child;
                //break;
            }
        }
    }



protected:
    ListPtrChildType list_ptr_childs;

};

template <typename ParentType, typename ChildType>
class GeoConnect{
public:
    using PtrParentType = std::shared_ptr<ParentType>;
    using ListPtrParentType = std::list<PtrParentType>;
    using PtrChildType = std::shared_ptr<ChildType>;
    using ListPtrChildType = std::list<PtrChildType>;

    GeoConnect(){}

    void addParent(const PtrParentType & ptr_parent){
        this->list_ptr_parents.push_back(ptr_parent);
    }

    void addChild(const PtrChildType & ptr_child){
        this->list_ptr_childs.push_back(ptr_child);
    }

    ListPtrChildType &getChild(){
        return this->list_ptr_childs;
    }

    ListPtrParentType &getParent(){
        return this->list_ptr_parents;
    }

    bool isMerged() const {
        return this->merged;
    }

    void setMerged(bool merged) {
        GeoConnect::merged = merged;
    }

    void merge(GeoConnect &rhs){
        auto & rhs_list_ptr_parent = rhs.getParent();
        for (auto &rhs_ptr_parent: rhs_list_ptr_parent) {
            this->addParent(rhs_ptr_parent);
        }
    }

    void replaceChild(PtrChildType &old_child, PtrChildType &new_child){
        for(auto & ptr_child :list_ptr_childs){
            if(ptr_child == old_child){
                ptr_child = new_child;
                //break;
            }
        }
    }

    void replaceParent(PtrParentType &old_parent, PtrParentType &new_parent){
        for(auto & ptr_parent :list_ptr_parents){
            if(ptr_parent == old_parent){
                ptr_parent = new_parent;
                //break;
            }
        }
    }

    void eraseParent(PtrParentType & ptr_parent){
        auto it = list_ptr_parents.begin();
        while (it != list_ptr_parents.end()){
            if((*it) == ptr_parent){
                list_ptr_parents.erase(it);
                //break;
            }
            ++it;
        }
    }



protected:
    ListPtrParentType list_ptr_parents;
    ListPtrChildType list_ptr_childs;
    bool merged{false};

};

template <typename ParentType>
class GeoEnd {
public:
    using PtrParentType = std::shared_ptr<ParentType>;
    using ListPtrParentType = std::list<PtrParentType>;

    GeoEnd(){}

    void addParent(const PtrParentType & ptr_parent){
        list_ptr_parents.push_back(ptr_parent);
    }

    ListPtrParentType &getParent(){
        return list_ptr_parents;
    }

    void replaceParent(PtrParentType &old_parent, PtrParentType &new_parent){
        for(auto & ptr_parent :list_ptr_parents){
            if(ptr_parent == old_parent){
                ptr_parent = new_parent;
                //break;
            }
        }
    }

    void eraseParent(PtrParentType & ptr_parent){
        auto it = list_ptr_parents.begin();
        while (it != list_ptr_parents.end()){
            if((*it) == ptr_parent){
                list_ptr_parents.erase(it);
                //break;
            }
            ++it;
        }
    }

protected:
    ListPtrParentType list_ptr_parents;
};



#endif //ARTPDE_GEO_TREE_HPP
