//
// Created by Chingkai Chou on 3/21/18.
//

#ifndef ARTCFD_POINT_HPP
#define ARTCFD_POINT_HPP

#include <vector>
#include <memory>
#include <ostream>
#include "DimensionUtility.hpp"

class PointBase{
public:
    PointBase(const size_t dim): data(dim, nullptr){
        for (auto &d:data) { d = std::make_shared<double>(0.0);}
        data.shrink_to_fit();
    }

    size_t getId() const {
        return id;
    }

    void setId(size_t id) {
        PointBase::id = id;
    }

    double& getDataByDim(size_t &Dim){ return (*data[Dim]);}

    std::shared_ptr<double>& getDataPointerByDim(size_t &Dim){ return data[Dim];}

protected:
    size_t id{0};
    std::vector<std::shared_ptr<double>> data;
};

template <class Dimension>
class Point;

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
    PtrPointType & getPoint(size_t id){ return data[id];}
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
    std::vector<PtrPointType> data;
};


template <>
class Point<Dim1D> : public PointBase{
    using OwnerType = std::shared_ptr<PointList<Dim1D>>;
public:
    Point(): PointBase(Dim1D::Dim){};
    Point(OwnerType &owner): PointBase(Dim1D::Dim), owner(owner){}
    double& x(){ return *data[0];};
    const OwnerType &getOwner() const {
        return owner;
    }

private:
    OwnerType owner{nullptr};
};

template <>
class Point<Dim2D> : public PointBase{
    using OwnerType = std::shared_ptr<PointList<Dim2D>>;
public:
    Point(): PointBase(Dim2D::Dim){};
    Point(OwnerType &owner): PointBase(Dim2D::Dim), owner(owner){}
    double& x(){ return *data[0];};
    double& y(){ return *data[1];};
    const OwnerType &getOwner() const {
        return owner;
    }

private:
    OwnerType owner{nullptr};
};


template <>
class Point<Dim3D> : public PointBase{
    using OwnerType = std::shared_ptr<PointList<Dim3D>>;
public:
    Point(): PointBase(Dim3D::Dim){};
    Point(OwnerType &owner): PointBase(Dim3D::Dim), owner(owner){}
    double& x(){ return *data[0];};
    double& y(){ return *data[1];};
    double& z(){ return *data[2];};
    const OwnerType &getOwner() const {
        return owner;
    }

private:
    OwnerType owner{nullptr};
};



#endif //ARTCFD_POINT_HPP
