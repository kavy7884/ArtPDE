//
// Created by Chingkai Chou on 3/8/18.
//

#ifndef ARTCFD_GEOMETRYDATA_HPP
#define ARTCFD_GEOMETRYDATA_HPP

#include "IO_Uility.hpp"
#include "CartesianVector.hpp"
#include "DataArray.hpp"
#include "DataTable.hpp"

enum class ElementType{None, Line, Triangle, Quadrilateral, Tetrahedron, Pyramid, Prism, Hexahedron};

template <class Dimension>
class GeometryData{
public:
    GeometryData() {}
    size_t GeoDim{Dimension::Dim};

    virtual bool readFile(IO_FileReader &IO_in) = 0;
    virtual bool writeFile(IO_FileWriter &IO_out) = 0;

    virtual bool DataProcessor() = 0;
};

template <class Dimension>
class GeometryMeshData : public GeometryData<Dimension>{
public:
    GeometryMeshData(): GeometryData<Dimension>(){};

    std::shared_ptr<CartesianVector<Dimension, double>> xNode{nullptr};
    std::shared_ptr<DataTable<size_t>> cElement30{nullptr};
    std::shared_ptr<DataArray<ElementType>> typeElement30{nullptr};

    virtual bool readFile(IO_FileReader &IO_in) override;
    virtual bool writeFile(IO_FileWriter &IO_out) override  {return true;} ;

    virtual bool DataProcessor() override {
        calElementType();
        return true;
    };
private:
    bool calElementType();

};

template <class Dimension>
class GeometryMeshFemData : public GeometryMeshData<Dimension>{
public:
    GeometryMeshFemData():GeometryMeshData<Dimension>(){};
    virtual bool readFile(IO_FileReader &IO_in) override;
    virtual bool writeFile(IO_FileWriter &IO_out)override{return true;};
    virtual bool DataProcessor() override{
        return GeometryMeshData<Dimension>::DataProcessor();
    }

};

template<class Dimension>
bool GeometryMeshData<Dimension>::readFile(IO_FileReader &IO_in) {
    std::ifstream fs;
    std::string fileName, meshPath, filePostFix;
    meshPath = IO_in.getProjectMeshPath();

    std::string bufferLine;
    std::vector<double> outLineDouble;
    std::vector<std::vector<double >> loadPointsTemp{this->GeoDim, std::vector<double >{}};
    std::vector<size_t> outLineSizeT;
    std::vector<std::vector<size_t >> loadConnectionTemp;

    // Load xNode
    fileName = "xNode";
    filePostFix = ".txt";
    fs.open(meshPath+fileName+filePostFix);
    if(fs.is_open()){

        for (auto &v : loadPointsTemp) {
            v.clear();
        }
        while(getline( fs, bufferLine )){
            IO_Basic::splitLineString<double>(bufferLine, outLineDouble);
            for (size_t i = 0; i < this->GeoDim; ++i) {
                loadPointsTemp[i].push_back(outLineDouble[i]);
            }
        }
        fs.close();

        xNode = std::make_shared<CartesianVector<Dimension, double>>();
        xNode->setVectorName(fileName);
        xNode->setSize(loadPointsTemp[0].size());
        for (size_t j = 0; j < this->GeoDim; ++j) {
            std::ostringstream issName;
            issName <<fileName<< "_" << j;
            xNode->getDataComponent(j) = DataArray<double>::create(issName.str()).dataFrom(loadPointsTemp[j]).build();
        }
        std::cout << "Loaded : " << fileName << std::endl;
    }else{
        //Todo - Exception
        std::cout << ">> Error -> Loaded : " << fileName << " Fail!" << std::endl;
    }


    // Load cElement30
    fileName = "cElement30";
    fs.open(meshPath+fileName+filePostFix);

    if(fs.is_open()){
        for (auto &v : loadConnectionTemp) {
            v.clear();
        }
        while(getline( fs, bufferLine )){
            IO_Basic::splitLineString<size_t>(bufferLine, outLineSizeT);
            loadConnectionTemp.push_back(outLineSizeT);
        }
        fs.close();
        cElement30 = DataTable<size_t>::create(fileName).dataFrom(loadConnectionTemp).build();

        std::cout << "Loaded : " << fileName << std::endl;
    }else{
        //Todo - Exception
        std::cout << ">> Error -> Loaded : " << fileName << " Fail!" << std::endl;
    }

    return true;
}

template<class Dimension>
bool GeometryMeshData<Dimension>::calElementType() {
    std::vector<ElementType> type;
    for (size_t i = 0; i < cElement30->getRowSize(); ++i) {
        //std::cout << cElement30->row(i).size() << std::endl;
        if(Dimension::Dim == 2){
            if(cElement30->row(i).size() == 3) type.push_back(ElementType::Triangle);
            else if(cElement30->row(i).size() == 4) type.push_back(ElementType::Quadrilateral);
        }
        else if(Dimension::Dim == 3){

        }
        else{
            type.push_back(ElementType::None);
        }
    }

    typeElement30 = DataArray<ElementType>::create("typeElement30").dataFrom(type).build();

    return true;
}

template<class Dimension>
bool GeometryMeshFemData<Dimension>::readFile(IO_FileReader &IO_in) {
    return GeometryMeshData<Dimension>::readFile(IO_in);
}

#endif //ARTCFD_GEOMETRYDATA_HPP
