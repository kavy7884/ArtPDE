//
// Created by Chingkai Chou on 3/8/18.
//

#ifndef ARTCFD_GEOMETRYDATA_HPP
#define ARTCFD_GEOMETRYDATA_HPP

#include "Point.hpp"

enum class ElementType{None, Line, Triangle, Quadrilateral, Tetrahedron, Pyramid, Prism, Hexahedron};

template <class Dimension>
class GeometryData{
public:
    GeometryData() {}
    size_t GeoDim{Dimension::Dim};

//    virtual bool DataProcessor() = 0;
};

template <class Dimension>
class GeometryMeshData : public GeometryData<Dimension>{
public:
    GeometryMeshData(): GeometryData<Dimension>(){
        xNode = std::make_shared<PointList<Dimension>>();
    };

    std::shared_ptr<PointList<Dimension>> xNode{nullptr};

//    std::shared_ptr<DataTable<size_t>> cElement30{nullptr};
    //std::shared_ptr<DataArray<ElementType>> typeElement30{nullptr};

//
//    virtual bool DataProcessor() override {
//        calElementType();
//        return true;
//    };
private:
//    bool calElementType();

};

template <class Dimension>
class GeometryMeshFemData : public GeometryMeshData<Dimension>{
public:
    GeometryMeshFemData():GeometryMeshData<Dimension>(){};

//    std::shared_ptr<DataArray<size_t>> constraintId{nullptr};
//    std::shared_ptr<DataArray<double>> constraintValue{nullptr};

//    virtual bool DataProcessor() override{
//        return GeometryMeshData<Dimension>::DataProcessor();
//    }

};

//template<class Dimension>
//bool GeometryMeshData<Dimension>::readFile(IO_FileReader &IO_in) {
//    std::ifstream fs;
//    std::string fileName, meshPath, filePostFix;
//    meshPath = IO_in.getProjectMeshPath();
//
//    std::string bufferLine;
//    size_t bufferSizeT;
//    double bufferDouble;
//    std::vector<double> outLineDouble;
//    std::vector<std::vector<double >> loadPointsTemp{this->GeoDim, std::vector<double >{}};
//    std::vector<size_t> outLineSizeT;
//    std::vector<std::vector<size_t >> loadConnectionTemp;
//
//    // Load xNode
//    fileName = "xNode";
//    filePostFix = ".txt";
//    fs.open(meshPath+fileName+filePostFix);
//    if(fs.is_open()){
//
//        for (auto &v : loadPointsTemp) {
//            v.clear();
//        }
//        while(getline( fs, bufferLine )){
//            IO_Basic::splitLineString<double>(bufferLine, outLineDouble);
//            for (size_t i = 0; i < this->GeoDim; ++i) {
//                loadPointsTemp[i].push_back(outLineDouble[i]);
//            }
//        }
//        fs.close();
//
//        xNode = std::make_shared<CartesianVector<Dimension, double>>();
//        xNode->setVectorName(fileName);
//        xNode->setSize(loadPointsTemp[0].size());
//        for (size_t j = 0; j < this->GeoDim; ++j) {
//            std::ostringstream issName;
//            issName <<fileName<< "_" << j;
//            xNode->getDataComponent(j) = DataArray<double>::create(issName.str()).dataFrom(loadPointsTemp[j]).build();
//        }
//        std::cout << "Loaded : " << fileName << std::endl;
//    }else{
//        //Todo - Exception
//        std::cout << ">> Error -> Loaded : " << fileName << " Fail!" << std::endl;
//    }
//
//
//    // Load cElement30
//    fileName = "cElement30";
//    fs.open(meshPath+fileName+filePostFix);
//
//    if(fs.is_open()){
//        for (auto &v : loadConnectionTemp) {
//            v.clear();
//        }
//        while(getline( fs, bufferLine )){
//            IO_Basic::splitLineString<size_t>(bufferLine, outLineSizeT);
//            loadConnectionTemp.push_back(outLineSizeT);
//        }
//        fs.close();
//        cElement30 = DataTable<size_t>::create(fileName).dataFrom(loadConnectionTemp).build();
//
//        std::cout << "Loaded : " << fileName << std::endl;
//    }else{
//        //Todo - Exception
//        std::cout << ">> Error -> Loaded : " << fileName << " Fail!" << std::endl;
//    }
//
//
//
//    return true;
//}
//
//template<class Dimension>
//bool GeometryMeshData<Dimension>::calElementType() {
//    std::vector<ElementType> type;
//    for (size_t i = 0; i < cElement30->getRowSize(); ++i) {
//        //std::cout << cElement30->row(i).size() << std::endl;
//        if(Dimension::Dim == 2){
//            if(cElement30->row(i).size() == 3) type.push_back(ElementType::Triangle);
//            else if(cElement30->row(i).size() == 4) type.push_back(ElementType::Quadrilateral);
//        }
//        else if(Dimension::Dim == 3){
//
//        }
//        else{
//            type.push_back(ElementType::None);
//        }
//    }
//
//    typeElement30 = DataArray<ElementType>::create("typeElement30").dataFrom(type).build();
//
//    return true;
//}
//
//template<class Dimension>
//bool GeometryMeshFemData<Dimension>::readFile(IO_FileReader &IO_in) {
//    GeometryMeshData<Dimension>::readFile(IO_in);
//
//    std::ifstream fs;
//    std::string fileName, meshPath, filePostFix;
//    meshPath = IO_in.getProjectMeshPath();
//    size_t bufferSizeT;
//    double bufferDouble;
//    std::vector<double> outLineDouble;
//    std::vector<size_t> outLineSizeT;
//
//    // Load constraint
//    fileName = "constraint";
//    filePostFix = ".txt";
//    fs.open(meshPath+fileName+filePostFix);
//    if(fs.is_open()){
//        outLineSizeT.clear();
//        outLineDouble.clear();
//        while(fs >> bufferSizeT >> bufferDouble){
//            outLineSizeT.push_back(bufferSizeT);
//            outLineDouble.push_back(bufferDouble);
//        }
//
//        fs.close();
//
//        constraintId = DataArray<size_t >::create(fileName).dataFrom(outLineSizeT).build();
//        constraintValue = DataArray<double>::create(fileName).dataFrom(outLineDouble).build();
//
//        std::cout << "Loaded : " << fileName << std::endl;
//    }else{
//        //Todo - Exception
//        std::cout << ">> Error -> Loaded : " << fileName << " Fail!" << std::endl;
//    }
//
//    return true;
//}

#endif //ARTCFD_GEOMETRYDATA_HPP
