//
// Created by Chingkai Chou on 4/4/18.
//

#ifndef ARTCFD_DOFUNIT_HPP
#define ARTCFD_DOFUNIT_HPP

#include <memory>
#include <vector>
#include "DimensionUtility.hpp"

class DofType{};
class DofScalar : public DofType{};
class DofVector : public DofType{};

class DofComponent{
public:
    DofComponent() {}

    double getData() const {
        return data;
    }

    void setData(double data) {
        DofComponent::data = data;
    }

private:
    double data{0.0};
};

class DofUnitBase{
public:
    DofUnitBase(size_t numOfComponent);

    size_t getComponentSize(){ return components.size(); };

    std::shared_ptr<DofComponent> &getComponent(size_t &direction);

protected:
    std::vector<std::shared_ptr<DofComponent>> components;
};

template <class Dimension, class DofType>
class DofUnit;

template <class Dimension>
class DofUnit<Dimension, DofScalar> : public DofUnitBase{
public:
    DofUnit() : DofUnitBase(1) {}

    std::shared_ptr<DofComponent> &getDof(){ return components[0];}
};

template <class Dimension> class DofUnit<Dimension, DofVector>;

template <>
class DofUnit<Dim1D, DofVector> : public DofUnitBase{
public:
    DofUnit() : DofUnitBase(1) {}
    std::shared_ptr<DofComponent> &getDof_x(){ return components[0];}
};

template <>
class DofUnit<Dim2D, DofVector> : public DofUnitBase{
public:
    DofUnit() : DofUnitBase(2) {}
    std::shared_ptr<DofComponent> &getDof_x(){ return components[0];}
    std::shared_ptr<DofComponent> &getDof_y(){ return components[1];}
};

template <>
class DofUnit<Dim3D, DofVector> : public DofUnitBase{
public:
    DofUnit() : DofUnitBase(3) {}
    std::shared_ptr<DofComponent> &getDof_x(){ return components[0];}
    std::shared_ptr<DofComponent> &getDof_y(){ return components[1];}
    std::shared_ptr<DofComponent> &getDof_z(){ return components[2];}
};

DofUnitBase::DofUnitBase(size_t numOfComponent) {
    for (size_t i = 0; i < numOfComponent; ++i) {
        components.emplace_back(std::make_shared<DofComponent>());
    }
    components.shrink_to_fit();
}

std::shared_ptr<DofComponent> &DofUnitBase::getComponent(size_t &direction) {
    if(direction >= getComponentSize()){
        //TODO - error exception
    }
    return components[direction];
}

#endif //ARTCFD_DOFUNIT_HPP
