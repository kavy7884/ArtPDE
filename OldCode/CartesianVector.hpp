//
// Created by Chingkai Chou on 3/2/18.
//

#ifndef ARTCFD_CARTESIANVECTOR_HPP
#define ARTCFD_CARTESIANVECTOR_HPP

#include <array>
#include "DimensionUtility.hpp"
#include "DataArray.hpp"

template <class Dimension, typename T>
class CartesianVectorBasic{
public:
    CartesianVectorBasic() {}

    size_t getSize() const {
        return size;
    }

    void setSize(size_t size) {
        CartesianVectorBasic::size = size;
    }

    const std::string &getVectorName() const {
        return VectorName;
    }

    void setVectorName(const std::string &VectorName) {
        CartesianVectorBasic::VectorName = VectorName;
    }

    std::shared_ptr<DataArray<T>> &getDataComponent(size_t Dim){
        return this->data[Dim];
    }

    template <class _Dimension, typename _T>
    friend std::ostream &operator<<(std::ostream &os, CartesianVectorBasic<_Dimension, _T> &vector) {
        os << "----------" << std::endl;
        os << "Cartesian Vector Name: " << vector.getVectorName() << std::endl;
        os << "Cartesian Vector Size: " << vector.getSize() << std::endl;
        os << "Cartesian Vector Contents: " << std::endl;
        for (size_t i = 0; i < vector.getSize(); ++i) {
            for (size_t d = 0; d < _Dimension::Dim; ++d) {
                os << vector.getDataComponent(d)->at(i) << "\t";
            }
            os << std::endl;
        }
        os << "----------" << std::endl;
        return os;
    }

protected:
    size_t size{0};
    std::string VectorName;
    std::array<std::shared_ptr<DataArray<T>>, Dimension::Dim> data{nullptr};

};

template <class Dimension, typename T>
class CartesianVector :public CartesianVectorBasic<Dimension, T>{
public:
    CartesianVector(): CartesianVectorBasic<Dimension, T>() {};

};

template <typename T>
class CartesianVector<Dim1D, T> :public CartesianVectorBasic<Dim1D, T>{
public:
    CartesianVector(): CartesianVectorBasic<Dim1D, T>() {};
    std::shared_ptr<DataArray<T>> &x(){
        return this->data[0];
    }

    T &x(size_t Id){
        return this->data[0]->at(Id);
    }

};


template <typename T>
class CartesianVector<Dim2D, T> :public CartesianVectorBasic<Dim2D, T>{
public:
    CartesianVector(): CartesianVectorBasic<Dim2D, T>() {};

    std::shared_ptr<DataArray<T>> &x(){
        return this->data[0];
    }
    std::shared_ptr<DataArray<T>> &y(){
        return this->data[1];
    }

    T &x(size_t Id){
        return this->data[0]->at(Id);
    }
    T &y(size_t Id){
        return this->data[1]->at(Id);
    }
};

template <typename T>
class CartesianVector<Dim3D, T> :public CartesianVectorBasic<Dim3D, T>{
public:
    CartesianVector(): CartesianVectorBasic<Dim3D, T>() {};

    std::shared_ptr<DataArray<T>> &x(){
        return this->data[0];
    }
    std::shared_ptr<DataArray<T>> &y(){
        return this->data[1];
    }
    std::shared_ptr<DataArray<T>> &z(){
        return this->data[2];
    }

    T &x(size_t Id){
        return this->data[0]->at(Id);
    }
    T &y(size_t Id){
        return this->data[1]->at(Id);
    }
    T &z(size_t Id){
        return this->data[2]->at(Id);
    }
};



#endif //ARTCFD_CARTESIANVECTOR_HPP
