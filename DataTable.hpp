//
// Created by Chingkai Chou on 2/17/18.
//

#ifndef ARTCFD_DATATABLE_HPP
#define ARTCFD_DATATABLE_HPP

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

template <typename T>
using StdVectorTensor2D = std::vector<std::vector<T>>;

template <typename T> class DataTableBuilder;
template <typename T> class DataTable_colArray;

template <typename T>
class DataTable {
    template <typename S> friend class DataTableBuilder;
    template <typename V> friend class DataTable_colArray;
    std::string tableName;
    size_t tableRowSize{0};
    T *values{nullptr};
    size_t *columnPtr{nullptr};

    void resetMemory(const size_t &rowSize, const size_t &valueSize);
    void checkRowRange(const size_t &row_pos) const;
    void checkColRange(const size_t &row_pos, const size_t &col_pos) const;
public:
    DataTable(){};
    DataTable(const DataTable& source);
    DataTable(DataTable&& source);
    virtual ~DataTable();

    size_t getRowSize() const;
    size_t getColSize(const size_t &row_pos) const;
    size_t getElementNumber() const;

    const std::string &getName() const {
        return tableName;
    }

    void setName(const std::string &tableName) {
        DataTable::tableName = tableName;
    }

    static DataTableBuilder<T> create(const std::string & tableName){
        return DataTableBuilder<T>(tableName);
    }

    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, const DataTable<U> &table);

    T& at(const size_t &pos_row, const size_t &pos_col);
    T& operator()(const size_t &pos_row, const size_t &pos_col) { return at(pos_row, pos_col); }

    DataTable& operator=(const DataTable& source);
    DataTable& operator=(DataTable&& source);

    DataTable_colArray<T> row(const size_t &rowNum);

};

template <typename T>
class DataTable_colArray{
public:

private:
    size_t rowNum{0}, colSize{0};
    DataTable<T> &tableData;
public:
    DataTable_colArray(const size_t &rowNumber, DataTable<T> & refTable)
            :rowNum(rowNumber), tableData(refTable){colSize = tableData.getColSize(rowNumber);}

    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, const DataTable_colArray<U> &array);;

    size_t getRowNum() const {
        return rowNum;
    }

    T& col(const size_t &pos_col);

    const size_t& size() const;

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
        return iterator(tableData.values);
    }

    iterator end()
    {
        return iterator(tableData.values + size());
    }

    const_iterator cbegin() const
    {
        return const_iterator(tableData.values);
    }

    const_iterator cend() const
    {
        return const_iterator(tableData.values + size());
    }
};

template<typename T>
const size_t &DataTable_colArray<T>::size() const{
    return colSize;
}

template<typename T>
T &DataTable_colArray<T>::col(const size_t &pos_col) {
    return tableData.at(rowNum, pos_col);
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const DataTable_colArray<T> &array) {
    os << "----------" << std::endl;
    os << "Row :" << array.rowNum << " in TableName: " << array.tableData.getName()<< std::endl;
    os << "Col Size: " << array.size() << std::endl;
    os << "Col Contents: " << std::endl;
    for (auto it = array.cbegin(); it != array.cend(); ++it) {
        os << *it << "\t";
    }
    os << std::endl;
    os << "----------" << std::endl;
    return os;
}

template <typename T>
class DataTableBuilder {
    std::shared_ptr<DataTable<T>> shPtr_dataTable;

public:
    DataTableBuilder(const std::string & tableName);

    std::shared_ptr<DataTable<T>> build(){
        return shPtr_dataTable;
    }

    DataTableBuilder<T>& dataFrom(std::shared_ptr<StdVectorTensor2D<T>> source);
    DataTableBuilder<T>& dataFrom(StdVectorTensor2D<T> source);

};

template<typename T>
DataTableBuilder<T>::DataTableBuilder(const std::string &tableName) {
    shPtr_dataTable = std::make_shared<DataTable<T>>();
    shPtr_dataTable->setName(tableName);
}

template<typename T>
DataTableBuilder<T> &DataTableBuilder<T>::dataFrom(std::shared_ptr<StdVectorTensor2D<T>> source) {
    (*shPtr_dataTable).tableRowSize = (*source).size();
    (*shPtr_dataTable).columnPtr = new size_t[(*source).size()+1];
    size_t numElement = 0;
    (*shPtr_dataTable).columnPtr[0] = 0;
    for (size_t i = 0; i < (*source).size(); ++i) {
        numElement += (*source)[i].size();
        (*shPtr_dataTable).columnPtr[i+1] = numElement;
    }
    (*shPtr_dataTable).values = new T[numElement];
    for (size_t j = 0; j < (*source).size(); ++j) {
        std::copy((*source)[j].begin(), (*source)[j].end(), (*shPtr_dataTable).values + (*shPtr_dataTable).columnPtr[j]);
    }
    return *this;
}

template<typename T>
DataTableBuilder<T> &DataTableBuilder<T>::dataFrom(StdVectorTensor2D<T> source) {
    (*shPtr_dataTable).tableRowSize = source.size();
    (*shPtr_dataTable).columnPtr = new size_t[source.size()+1];
    size_t numElement = 0;
    (*shPtr_dataTable).columnPtr[0] = 0;
    for (size_t i = 0; i < source.size(); ++i) {
        numElement += source[i].size();
        (*shPtr_dataTable).columnPtr[i+1] = numElement;
    }
    (*shPtr_dataTable).values = new T[numElement];
    for (size_t j = 0; j < source.size(); ++j) {
        std::copy(source[j].begin(), source[j].end(), (*shPtr_dataTable).values + (*shPtr_dataTable).columnPtr[j]);
    }
    return *this;
}

template<typename T>
DataTable<T>::~DataTable() {
    if(values != nullptr) {
        delete[](values);
    };
    if(columnPtr != nullptr) {
        delete[](columnPtr);
    };
    tableRowSize = 0;
}

template<typename T>
size_t DataTable<T>::getRowSize() const {
    return tableRowSize;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const DataTable<T> &table) {
    os << "----------" << std::endl;
    os << "Table Name: " << table.tableName << std::endl;
    os << "Table RowSize: " << table.getRowSize() << std::endl;
    os << "Table Contents: " << std::endl;
    for (size_t i = 0; i < table.getRowSize(); ++i){
        for (size_t j = table.columnPtr[i]; j < table.columnPtr[i+1]; ++j) {
            std::cout << table.values[j];
            if(j < table.columnPtr[i+1]-1){ std::cout << "\t";}
        }
        std::cout << std::endl;
    }

    os << "----------" << std::endl;
    return os;
}

template<typename T>
T &DataTable<T>::at(const size_t &pos_row, const size_t &pos_col) {
    checkRowRange(pos_row);
    checkColRange(pos_row, pos_col);
    return values[columnPtr[pos_row] + pos_col];
}

template<typename T>
DataTable<T>::DataTable(const DataTable<T> &source) {
    if(&source != this) {
        tableName = source.getName();
        tableRowSize = source.getRowSize();
        columnPtr = new size_t[source.getRowSize()+1];
        values = new T[source.columnPtr[source.getRowSize()]];
        std::copy(&source.columnPtr[0], &source.columnPtr[0] + source.getRowSize()+1, &columnPtr[0]);
        std::copy(&source.values[0], &source.values[0] + source.columnPtr[source.getRowSize()], &values[0]);
    }
}

template<typename T>
DataTable<T>::DataTable(DataTable<T> &&source) {
    if(&source != this) {
        tableName = std::move(source.tableName);
        tableRowSize = source.getRowSize();
        source.tableRowSize = 0;
        columnPtr = std::move(source.columnPtr);
        values = std::move(source.values);
        source.columnPtr = nullptr;
        source.values = nullptr;
    }
}

template<typename T>
DataTable<T> &DataTable<T>::operator=(const DataTable<T> &source) {
    if(&source != this){
        resetMemory(source.getRowSize(), source.getElementNumber());
        std::copy(&source.columnPtr[0], &source.columnPtr[0] + source.getRowSize()+1, &columnPtr[0]);
        std::copy(&source.values[0], &source.values[0] + source.columnPtr[source.getRowSize()], &values[0]);
        tableRowSize = source.getRowSize();
        tableName = source.getName();
        return *this;
    }
    return *this;
}

template<typename T>
DataTable<T> &DataTable<T>::operator=(DataTable<T> &&source) {
    if(&source != this){
        tableName = std::move(source.tableName);
        tableRowSize = source.getRowSize();
        source.tableRowSize = 0;
        columnPtr = std::move(source.columnPtr);
        values = std::move(source.values);
        source.columnPtr = nullptr;
        source.values = nullptr;
    }
    return *this;
}

template<typename T>
void DataTable<T>::resetMemory(const size_t &rowSize, const size_t &valueSize) {
    if(columnPtr != nullptr) {
        delete[](columnPtr);
        columnPtr = nullptr;
    }
    if(values != nullptr) {
        delete[](values);
        values = nullptr;
    }
    columnPtr = new size_t[rowSize + 1];
    values = new T[valueSize];
}

template<typename T>
size_t DataTable<T>::getElementNumber() const {
    if(columnPtr != nullptr || getRowSize() >0){
        return columnPtr[getRowSize()];
    }
    return 0;
}

template<typename T>
DataTable_colArray<T> DataTable<T>::row(const size_t &rowNum) {
    checkRowRange(rowNum);
    return DataTable_colArray<T>(rowNum, *this);
}

template<typename T>
void DataTable<T>::checkRowRange(const size_t &row_pos) const{
    if (row_pos >= getRowSize() || row_pos < 0) {
        std::ostringstream errMsg;
        errMsg << ">> Error<DataTable>::access Row_Pos = " << row_pos
               << " out of range! ( Row Size = " << getRowSize() << " )";
        throw std::out_of_range(errMsg.str());
    }
}

template<typename T>
void DataTable<T>::checkColRange(const size_t &row_pos, const size_t &col_pos) const{
    if (col_pos >= getColSize(row_pos) || col_pos < 0) {
        std::ostringstream errMsg;
        errMsg << ">> Error<DataTable>::access Col_Pos = " << col_pos
               << " out of range! ( Col Size = " << getColSize(row_pos) << " )";
        throw std::out_of_range(errMsg.str());
    }
}

template<typename T>
size_t DataTable<T>::getColSize(const size_t &row_pos) const {
    checkRowRange(row_pos);
    return columnPtr[row_pos + 1] - columnPtr[row_pos];
}

#endif //ARTCFD_DATATABLE_HPP
