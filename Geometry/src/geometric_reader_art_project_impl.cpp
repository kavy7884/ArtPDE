//
// Created by Chingkai Chou on 6/3/18.
//

#include "../inc/GeometricReader/geometric_reader_art_project.hpp"

template<typename GeometricReaderType, size_t Dimension>
bool ArtProjectReader<GeometricReaderType, Dimension>::read() {
    std::cout << ">> (Start) Reading ArtPDE project format files..." << std::endl;

    if(!readPosition()) return false;

    if(Dimension == 2){
        // Loading Face
    }

    if(Dimension == 3){
        // Loading Cell
    }

    std::cout << ">> ( End ) Reading ArtPDE project format files..." << std::endl;
    return true;
}

template<typename GeometricReaderType, size_t Dimension>
bool ArtProjectReader<GeometricReaderType, Dimension>::readPosition() {
    std::cout << ">>>> (Start) Loading Position..." << std::endl;
    std::string file_name;
    file_name += art_project->getProjectGeometryPath();
    file_name += "position.txt";

//                    std::cout << file_name << std::endl;

    using PointType = typename GeometricReaderType::type::PointType;
    using PtrPointType = typename GeometricReaderType::type::PtrPointType;
    using VertexType = typename GeometricReaderType::type::VertexType_traits::VertexType;
    using PtrVertexType = typename GeometricReaderType::type::VertexType_traits::PtrVertexType;

    std::ifstream fs;
    std::string bufferLine;
    double temp_pt[Dimension];
    PtrPointType ptr_point;
    PtrVertexType ptr_vertex;

    fs.open(file_name);
    if(fs.is_open()){
        while(getline( fs, bufferLine )) {
            std::stringstream w(bufferLine);

            for(size_t i = 0; i< Dimension; ++i)
                w >> temp_pt[i];

            ptr_point = std::make_shared<PointType>();
            ptr_vertex = std::make_shared<VertexType>();

            for(size_t i = 0; i< Dimension; ++i)
                ptr_point->setDataById(0, temp_pt[i]);

            ptr_data->getTotalVec_PtrVertex().push_back(ptr_vertex);
        }
    }else
        return false;

    std::cout << ">>>> ( End ) Loading Position..." << std::endl;
    return true;
}

