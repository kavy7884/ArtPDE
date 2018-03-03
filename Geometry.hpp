//
// Created by Chingkai Chou on 2/23/18.
//

#ifndef ARTCFD_GEOMETRY_HPP
#define ARTCFD_GEOMETRY_HPP

#include <iostream>
#include "DimensionUtility.hpp"
#include "IO_Uility.hpp"

template <class Dimension>
class GeometryData{
public:
    GeometryData() {}
    size_t GeoDim{Dimension::Dim};

    virtual bool readFile(IO_FileReader &IO_in) = 0;
    virtual bool writeFile(IO_FileWriter &IO_out) = 0;
};

template <class Dimension, template<class > class Data = GeometryData>
class Geometry {
public:
    Geometry(){
        GeoData = std::make_shared<Data<Dimension>>();
    };
    const std::shared_ptr<Data<Dimension>> &getGeoData() const {
        return GeoData;
    }
private:
    std::shared_ptr<Data<Dimension>> GeoData{nullptr};
};





//
//
//template <typename T>
//class CartesianVector{
//    size_t size{0};
//    std::string VectorName;
//    Dimension &Dim;
//public:
//    CartesianVector(Dimension &dim) : Dim(dim){}
//    std::shared_ptr<DataArray<T>> x{nullptr}, y{nullptr}, z{nullptr};
////    std::shared_ptr<DataArray<T>> getContentsByDimension(size_t Dim){
////        if(Dim == 0) return x;
////        else if(Dim == 1) return y;
////        else if(Dim == 2) return z;
////            // TODO - error exception
////        else return z;
////    }
//
//    void setComponentByDimension(size_t dim, std::shared_ptr<DataArray<T>> ptrArray){
//        if(dim < size_t(Dim) ){
//            if(dim == 0) x = ptrArray;
//            else if(dim == 1) y = ptrArray;
//            else if(dim == 2) z = ptrArray;
//        }
//        else{
//            // TODO - Exception
//        }
//    }
//
//    const std::shared_ptr<DataArray<T>>& getComponentByDimension(size_t dim) const {
//        if(dim < size_t(Dim) ){
//            if(dim == 0) return x;
//            else if(dim == 1) return y;
//            else if(dim == 2) return z;
//        }
//        else{
//            // TODO - Exception
//            return z;
//
//        }
//    }
//
//    size_t getSize() const {
//        return size;
//    }
//
//    void setSize(size_t size) {
//        CartesianVector::size = size;
//    }
//
//    const std::string &getName() const {
//        return VectorName;
//    }
//
//    void setName(const std::string &Name) {
//        CartesianVector::VectorName = Name;
//    }
//
//    friend std::ostream &operator<<(std::ostream &os, const CartesianVector &vector) {
//        os << "----------" << std::endl;
//        os << "Cartesian Vector Name: " << vector.getName() << std::endl;
//        os << "Cartesian Vector Size: " << vector.getSize() << std::endl;
//        os << "Cartesian Vector Contents: " << std::endl;
//        for (int i = 0; i < vector.getSize(); ++i) {
//            os << vector.x->at(i) << "\t";
//            if(vector.y != nullptr) os << vector.y->at(i) << "\t";
//            if(vector.z != nullptr) os << vector.z->at(i) << "\t";
//            os << std::endl;
//        }
//        os << "----------" << std::endl;
//        return os;
//    }
//};
//







#endif //ARTCFD_GEOMETRY_HPP
