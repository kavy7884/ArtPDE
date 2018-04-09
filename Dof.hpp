//
// Created by Chingkai Chou on 4/4/18.
//

#ifndef ARTCFD_DOF_HPP
#define ARTCFD_DOF_HPP

#include <string>
#include <ostream>
#include "DofUnit.hpp"
#include "Point.hpp"

class DofBase{
public:
    DofBase(const std::string &dofName, size_t dofUnitSize) : dofName(dofName), dofUnitSize(dofUnitSize) {}
    DofBase(const std::string &dofName) : dofName(dofName) {}

    const std::string &getDofName() const {
        return dofName;
    }

    void setDofName(const std::string &dofName) {
        DofBase::dofName = dofName;
    }

    size_t getDofUnitSize() const {
        return dofUnitSize;
    }

    void setDofUnitSize(size_t dofDofUnitSize) {
        DofBase::dofUnitSize = dofDofUnitSize;
    }


protected:
    std::string dofName;
    size_t dofUnitSize{0};
};

template <class Dimension, class DofType>
class Dof : public DofBase{
public:
    Dof(const std::string &dofName, size_t dofUnitLength);
    Dof(const std::string &dofName, PointList<Dimension> &ptList);

    friend std::ostream &operator<<(std::ostream &os, const Dof &dof){
        for (size_t i = 0; i < dof.getDofUnitSize(); ++i) {
            os << *dof.c_getDofUnit(i) << std::endl;
        }
        return os;
    }

    std::shared_ptr<DofUnit<Dimension, DofType>> & getDofUnit(size_t &id){ return dofUnit[id];}

    const std::shared_ptr<DofUnit<Dimension, DofType>> & c_getDofUnit(size_t &id)const { return dofUnit[id];}

private:
    std::vector<std::shared_ptr<DofUnit<Dimension, DofType>>> dofUnit;
};

template<class Dimension, class DofType>
Dof<Dimension, DofType>::Dof(const std::string &dofName, size_t dofUnitLength):DofBase(dofName, dofUnitLength) {
    for (size_t i = 0; i < this->getDofUnitSize(); ++i) {
        dofUnit.emplace_back(std::make_shared<DofUnit<Dimension, DofType>>());
    }
}

template<class Dimension, class DofType>
Dof<Dimension, DofType>::Dof(const std::string &dofName, PointList<Dimension> &PointList) : DofBase(dofName){
    this->setDofUnitSize(PointList.size());
    for (size_t i = 0; i < this->getDofUnitSize(); ++i) {
        dofUnit.emplace_back(std::make_shared<DofUnit<Dimension, DofType>>(*PointList.getPoint(i)));
    }
}


#endif //ARTCFD_DOF_HPP
