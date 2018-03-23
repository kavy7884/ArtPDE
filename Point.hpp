//
// Created by Chingkai Chou on 3/21/18.
//

#ifndef ARTCFD_POINT_HPP
#define ARTCFD_POINT_HPP

#include <vector>
#include <array>
#include <memory>
#include <ostream>
#include "DimensionUtility.hpp"

template <class Dimension>
class PointBase{
public:
    PointBase() {}

    size_t getId() const {
        return id;
    }

    void setId(size_t id) {
        PointBase::id = id;
    }

    double& getDataByDim(size_t &Dim){ return data[Dim];};

protected:
    size_t id{0};
    std::array<double, Dimension::Dim> data{0.0};
};

template <class Dimension>
class Point: public PointBase<Dimension>{};

template <class Dimension>
class PointList{
public:
    using PtrPointType = std::shared_ptr<Point<Dimension>>;
    PointList() {}
    size_t size() const { return data.size();}
    void addPoint(PtrPointType & p_point){
        p_point->setId(size());
        data.push_back(p_point);
    }
    PtrPointType & getPoint(size_t& id){ return data[id];}
    const PtrPointType & c_getPoint(size_t& id)const{ return data[id];}

    friend std::ostream &operator<<(std::ostream &os, const PointList<Dimension> &list){
        for (size_t i = 0; i < list.size(); ++i) {
            os << list.c_getPoint(i)->getId()<< "\t";
            for (size_t j = 0; j < Dimension::Dim; ++j) {
                os << list.c_getPoint(i)->getDataByDim(j) << "\t";
            }
            os << std::endl;
        }
        return os;
    }

private:
    std::vector<std::shared_ptr<Point<Dimension>>> data;
};


template <>
class Point<Dim1D> : public PointBase<Dim1D>{
    using OwnerType = std::shared_ptr<PointList<Dim1D>>;
public:
    Point(): PointBase<Dim1D>(){};
    Point(OwnerType &owner): PointBase<Dim1D>(), owner(owner){}
    double& x(){ return data[0];};
    const OwnerType &getOwner() const {
        return owner;
    }

private:
    OwnerType owner{nullptr};
};

template <>
class Point<Dim2D> : public PointBase<Dim2D>{
    using OwnerType = std::shared_ptr<PointList<Dim2D>>;
public:
    Point(): PointBase<Dim2D>(){};
    Point(OwnerType &owner): PointBase<Dim2D>(), owner(owner){}
    double& x(){ return data[0];};
    double& y(){ return data[1];};
    const OwnerType &getOwner() const {
        return owner;
    }

private:
    OwnerType owner{nullptr};
};


template <>
class Point<Dim3D> : public PointBase<Dim3D>{
    using OwnerType = std::shared_ptr<PointList<Dim3D>>;
public:
    Point(): PointBase<Dim3D>(){};
    Point(OwnerType &owner): PointBase<Dim3D>(), owner(owner){}
    double& x(){ return data[0];};
    double& y(){ return data[1];};
    double& z(){ return data[2];};
    const OwnerType &getOwner() const {
        return owner;
    }

private:
    OwnerType owner{nullptr};
};



#endif //ARTCFD_POINT_HPP
