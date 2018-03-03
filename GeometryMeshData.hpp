//
// Created by Chingkai Chou on 3/2/18.
//

#ifndef ARTCFD_GEOMETRYMESHDATA_HPP
#define ARTCFD_GEOMETRYMESHDATA_HPP

#include "DataTable.hpp"
#include "CartesianVector.hpp"
#include "Geometry.hpp"

template <class Dimension>
class GeometryMeshData : public GeometryData<Dimension>{
public:
    GeometryMeshData(): GeometryData<Dimension>(){};

    std::shared_ptr<CartesianVector<Dimension, double>> xNode{nullptr};
    std::shared_ptr<DataTable<size_t>> cElement30{nullptr};

    virtual bool readFile(IO_FileReader &IO_in) override;
    virtual bool writeFile(IO_FileWriter &IO_out) override  {return true;} ;

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

#endif //ARTCFD_GEOMETRYMESHDATA_HPP
