//
// Created by Chingkai Chou on 2/24/18.
//

#ifndef ARTCFD_GEOMETRYFEM_IO_HPP
#define ARTCFD_GEOMETRYFEM_IO_HPP

#include "GeometryDataFEM.hpp"

template <class IO_FileReader>
class GeometryReader{
    Geometry<GeometryDataFEM> &Geo;
    IO_FileReader &IO_reader;
public:
    GeometryReader(Geometry<GeometryDataFEM> &Geo, IO_FileReader &MesgIO) : Geo(Geo), IO_reader(MesgIO){}
    bool read();

};

template<class IO_FileReader>
bool GeometryReader<IO_FileReader>::read() {

    std::ifstream fs;
    std::string fileName, meshPath, filePostFix;
    meshPath = IO_reader.getProjectMeshPath();

    // Load xNode
    fileName = "xNode";
    filePostFix = ".txt";
    fs.open(meshPath+fileName+filePostFix);
    std::string bufferLine;
    std::vector<double> outLineDouble;
    std::vector<std::vector<double >> loadPointsTemp{size_t(Geo.GeoDim), std::vector<double >{}};

    for (auto &v : loadPointsTemp) {
        v.clear();
    }
    while(getline( fs, bufferLine )){
        IO_Basic::splitLineString<double>(bufferLine, outLineDouble);
        for (size_t i = 0; i < size_t(Geo.GeoDim); ++i) {
            loadPointsTemp[i].push_back(outLineDouble[i]);
        }
    }
    fs.close();

    Geo.GeoData->xNode = std::make_shared<PositionVector<double>>();
    Geo.GeoData->xNode->setPositionName(fileName);
    Geo.GeoData->xNode->setSize(loadPointsTemp[0].size());
    for (size_t j = 0; j < size_t(Geo.GeoDim); ++j) {
        std::ostringstream issName;
        issName <<fileName<< "_" << j;
        Geo.GeoData->xNode->getContentsByDimension(j) = DataArray<double>::create(issName.str()).dataFrom(loadPointsTemp[j]).build();
    }
    std::cout << "Loaded : " << fileName << std::endl;


    fileName = "cElement30";
    fs.open(meshPath+fileName+filePostFix);
    std::vector<size_t> outLineSizeT;
    std::vector<std::vector<size_t >> loadConnectionTemp;
    for (auto &v : loadConnectionTemp) {
        v.clear();
    }
    while(getline( fs, bufferLine )){
        IO_Basic::splitLineString<size_t>(bufferLine, outLineSizeT);
        loadConnectionTemp.push_back(outLineSizeT);
    }
    fs.close();
    Geo.GeoData->cElement30 = DataTable<size_t>::create(fileName).dataFrom(loadConnectionTemp).build();

    std::cout << "Loaded : " << fileName << std::endl;





    return false;
}

template <class IO_FileWriter>
class GeometryWriter{
    Geometry<GeometryDataFEM> &Geo;
    IO_FileWriter &IO_writer;
public:
    GeometryWriter(Geometry<GeometryDataFEM> &Geo, IO_FileWriter &MesgIO) : Geo(Geo), IO_writer(MesgIO){}
    bool write() {
        std::fstream fs;
        std::string fileName;

        fileName += IO_writer.getProjectResultsPath()+ "test.txt";
        fs.open(fileName, std::fstream::out);
        fs << "Haha" << std::endl;
        fs.close();
        return true;}

};

#endif //ARTCFD_GEOMETRYFEM_IO_HPP
