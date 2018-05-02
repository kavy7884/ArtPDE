//
// Created by Chingkai Chou on 2/16/18.
//

#ifndef ARTCFD_DATAARRAY_HPP
#define ARTCFD_DATAARRAY_HPP


#include <vector>
#include <string>
#include <sstream>
#include <ostream>
#include <stdexcept>

template <typename T> class DataArrayBuilder;

template <typename T>
class DataArray {
    template <typename S> friend class DataArrayBuilder;
    size_t arraySize{0};
    T *array{nullptr};
    std::string arrayName;

    void resetArray(const size_t &length);
    void checkSize(const size_t &pos);

public:
    DataArray(){};
    DataArray(const DataArray& source);
    DataArray(DataArray&& source);
    virtual ~DataArray();

    const size_t& size() const { return arraySize;}

    const std::string &getName() const;
    void setName(const std::string &arrayName);

    static DataArrayBuilder<T> create(const std::string & arrayName){
        return DataArrayBuilder<T>(arrayName);
    }

    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, const DataArray<U> &array);

    DataArray& operator=(const DataArray& source);
    DataArray& operator=(DataArray&& source);

    T& at(const size_t &pos);
    T& operator()(const size_t & pos) { return at(pos); }

    class iterator: public std::iterator<std::forward_iterator_tag, T>{
        T* ptr_;
    public:
        iterator(T* ptr) : ptr_(ptr) {};
        iterator(const iterator& source) : ptr_(source.ptr_) {}
        iterator& operator++() {++ptr_;return *this;}
        iterator operator++(int) {iterator tmp(*this); operator++(); return tmp;}
        T& operator*() { return *ptr_; }
        T* operator->() { return ptr_; }
        bool operator==(const iterator& rhs) const {return ptr_==rhs.ptr_;}
        bool operator!=(const iterator& rhs) const {return ptr_!=rhs.ptr_;}
    };

    class const_iterator: public std::iterator<std::forward_iterator_tag, T>{
        T* ptr_;
    public:
        const_iterator(T* ptr) : ptr_(ptr) {};
        const_iterator(const iterator& source) : ptr_(source.ptr_) {}
        const_iterator& operator++() {++ptr_;return *this;}
        const_iterator operator++(int) {const_iterator tmp(*this); operator++(); return tmp;}
        const T& operator*() { return *ptr_; }
        const T* operator->() { return ptr_; }
        bool operator==(const const_iterator& rhs) const {return ptr_==rhs.ptr_;}
        bool operator!=(const const_iterator& rhs) const {return ptr_!=rhs.ptr_;}
    };

    iterator begin()
    {
        return iterator(array);
    }

    iterator end()
    {
        return iterator(array + arraySize);
    }

    const_iterator cbegin() const
    {
        return const_iterator(array);
    }

    const_iterator cend() const
    {
        return const_iterator(array + arraySize);
    }


};

template <typename T>
class DataArrayBuilder {
    std::shared_ptr<DataArray<T>> shPtr_dataArray;
public:
    DataArrayBuilder(const std::string & arrayName);

    DataArrayBuilder<T>& dataFrom(std::shared_ptr<std::vector<T>> source);
    DataArrayBuilder<T>& dataFrom(std::vector<T> source);

    std::shared_ptr<DataArray<T>> build(){
        return shPtr_dataArray;
    }
};

template<typename T>
DataArrayBuilder<T>::DataArrayBuilder(const std::string & arrayName) {
    shPtr_dataArray = std::make_shared<DataArray<T>>();
    shPtr_dataArray->setName(arrayName);
}

template<typename T>
DataArrayBuilder<T>& DataArrayBuilder<T>::dataFrom(std::shared_ptr<std::vector<T>> source) {
    (*shPtr_dataArray).arraySize = (*source).size();
    (*shPtr_dataArray).array = new T[(*source).size()]();
    for (size_t i = 0; i < (*source).size(); ++i) {
        (*shPtr_dataArray).array[i] = (*source)[i];
    }
    return *this;
}

template<typename T>
DataArrayBuilder<T> &DataArrayBuilder<T>::dataFrom(std::vector<T> source) {
    (*shPtr_dataArray).arraySize = source.size();
    (*shPtr_dataArray).array = new T[source.size()]();
    for (size_t i = 0; i < source.size(); ++i) {
        (*shPtr_dataArray).array[i] = source[i];
    }
    return *this;
}

template<typename T>
DataArray<T>::~DataArray() {
    if(array != nullptr) {
        delete[](array);
        arraySize = 0;
    };
}

template<typename T>
const std::string &DataArray<T>::getName() const {
    return arrayName;
}

template<typename T>
void DataArray<T>::setName(const std::string &arrayName) {
    DataArray::arrayName = arrayName;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const DataArray<T> &array) {
    os << "----------" << std::endl;
    os << "Array Name: " << array.arrayName << std::endl;
    os << "Array Size: " << array.arraySize << std::endl;
    os << "Array Contents: " << std::endl;
    for (int i = 0; i < array.arraySize; ++i) {
        os << array.array[i] << std::endl;
    }
    os << "----------" << std::endl;
    return os;
}

template<typename T>
DataArray<T> &DataArray<T>::operator=(const DataArray<T> &source) {

    if(&source != this){
        resetArray(source.size());
        std::copy(&source.array[0], &source.array[0] + source.size(), &array[0]);
        arraySize = source.size();
        arrayName = source.getName();
        return *this;
    }
    return *this;
}

template<typename T>
void DataArray<T>::resetArray(const size_t &length) {
    if(array != nullptr) {
        delete[](array);
        array = nullptr;
    }
    array = new T[length];
}

template<typename T>
DataArray<T>::DataArray(const DataArray<T> &source) {
    arrayName = source.getName();
    arraySize = source.size();
    array = new T[arraySize];
    std::copy(&source.array[0], &source.array[0] + source.size(), &array[0]);
}

template<typename T>
DataArray<T> &DataArray<T>::operator=(DataArray<T> &&source) {
    arrayName = std::move(source.arrayName);
    arraySize = source.size();
    source.arraySize = 0;
    array = std::move(source.array);
    source.array = nullptr;

    return *this;
}

template<typename T>
DataArray<T>::DataArray(DataArray<T> &&source) {
    if(&source != this) {
        arrayName = std::move(source.arrayName);
        arraySize = source.size();
        source.arraySize = 0;
        array = std::move(source.array);
        source.array = nullptr;
    }
}

template<typename T>
T &DataArray<T>::at(const size_t &pos) {
    checkSize(pos);
    return array[pos];
}

template<typename T>
void DataArray<T>::checkSize(const size_t &pos) {
    if (pos >= arraySize || pos < 0) {
        std::ostringstream errMsg;
        errMsg << ">> Error<DataArray>::access Pos = " << pos
               << " out of range! ( Size = " << size() << " )";
        throw std::out_of_range(errMsg.str());
    }
}


#endif //ARTCFD_DATAARRAY_HPP
