//
// Created by Chingkai Chou on 5/17/18.
//

#ifndef ARTPDE_MARVIN_GEO_TREE_HPP
#define ARTPDE_MARVIN_GEO_TREE_HPP

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

    const ListPtrChildType &getChild(){
        return list_ptr_childs;
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
        list_ptr_parents.push_back(ptr_parent);
    }

    void addChild(const PtrChildType & ptr_child){
        list_ptr_childs.push_back(ptr_child);
    }

    const ListPtrChildType &getChild(){
        return list_ptr_childs;
    }

    bool isMerged() const {
        return merged;
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

protected:
    ListPtrParentType list_ptr_parents;
};



#endif //ARTPDE_MARVIN_GEO_TREE_HPP
