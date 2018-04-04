//
// Created by Chingkai Chou on 4/4/18.
//

#ifndef ARTCFD_DOFUNIT_HPP
#define ARTCFD_DOFUNIT_HPP

#include <memory>
#include <vector>
#include <ostream>
#include "DimensionUtility.hpp"
#include "Point.hpp"

class DofType{};
class DofScalar : public DofType{};
class DofVector : public DofType{};

class DofComponent{
public:
    DofComponent() { data = std::make_shared<double>(0.0);}

    DofComponent(const std::shared_ptr<double> &data) : data(data) {}

    double getData() const {
        return *data;
    }

    void setData(double data) {
        *this->data = data;
    }

private:
    std::shared_ptr<double>data{nullptr};
};

class DofUnitBase{
public:
    DofUnitBase(size_t numOfComponent);

    DofUnitBase(size_t numOfComponent, PointBase &pt);

    size_t getComponentSize()const { return components.size(); };

    std::shared_ptr<DofComponent> &getComponent(size_t &direction);

    const std::shared_ptr<DofComponent> &c_getComponent(size_t &direction) const;

    friend std::ostream &operator<<(std::ostream &os, const DofUnitBase &obj);

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
    DofUnit(Point<Dim1D> & pt) : DofUnitBase(1, pt){}

    std::shared_ptr<DofComponent> &getDof_x(){ return components[0];}
};

template <>
class DofUnit<Dim2D, DofVector> : public DofUnitBase{
public:
    DofUnit() : DofUnitBase(2) {}
    DofUnit(Point<Dim2D> & pt) : DofUnitBase(2, pt){}

    std::shared_ptr<DofComponent> &getDof_x(){ return components[0];}
    std::shared_ptr<DofComponent> &getDof_y(){ return components[1];}
};

template <>
class DofUnit<Dim3D, DofVector> : public DofUnitBase{
public:
    DofUnit() : DofUnitBase(3) {}
    DofUnit(Point<Dim3D> & pt) : DofUnitBase(3, pt){}
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

const std::shared_ptr<DofComponent> &DofUnitBase::c_getComponent(size_t &direction) const{
    if(direction >= getComponentSize()){
        //TODO - error exception
    }
    return components[direction];
}

std::ostream &operator<<(std::ostream &os, const DofUnitBase &obj) {
    for (size_t i = 0; i < obj.getComponentSize(); ++i) {
        os << obj.c_getComponent(i)->getData() << "\t";
    }
    return os;
}

DofUnitBase::DofUnitBase(size_t numOfComponent, PointBase &pt) {
    for (size_t i = 0; i < numOfComponent; ++i) {
            components.emplace_back(std::make_shared<DofComponent>(pt.getDataPointerByDim(i)));
    }
    components.shrink_to_fit();
}

#endif //ARTCFD_DOFUNIT_HPP
