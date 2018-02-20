//
// Created by Chingkai Chou on 2/17/18.
//

#ifndef ARTCFD_DATATABLE_HPP
#define ARTCFD_DATATABLE_HPP

#include <vector>
#include <string>
#include <ostream>

template <typename T>
using StdVectorTensor2D = std::vector<std::vector<T>>;

template <typename T> class DataTableBuilder;

template <typename T>
class DataTable {
    template <typename S> friend class DataTableBuilder;
    std::string tableName;
    size_t tableRowSize;
    T *values;
    size_t *columnPtr;
public:
    DataTable();
    virtual ~DataTable();

    size_t getRowSize() const;

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


};

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
DataTable<T>::DataTable() {
    tableRowSize = 0;
    values = nullptr;
    columnPtr = nullptr;
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


#endif //ARTCFD_DATATABLE_HPP
