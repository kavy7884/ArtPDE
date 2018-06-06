//
// Created by Chingkai Chou on 6/3/18.
//

#include "../inc/GeometricReader/geometric_reader_art_project.hpp"

template<typename GeometricReaderType, size_t Dimension>
bool ArtProjectReader<GeometricReaderType, Dimension>::read() {
    std::cout << ">> (Start) Reading ArtPDE project format files..." << std::endl;

    if(!readPosition()) return false;

    if(Dimension == 2){
        if(!readFace()) return false;
    }

    if(Dimension == 3){
        if(!readCell()) return false;
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
    auto &vec_ptr_vertex = ptr_data->getTotalVec_PtrVertex();

    fs.open(file_name);
    if(fs.is_open()){
        while(getline( fs, bufferLine )) {
            std::stringstream w(bufferLine);

            for(size_t i = 0; i< Dimension; ++i)
                w >> temp_pt[i];

            ptr_vertex = std::make_shared<VertexType>();
            ptr_point = std::make_shared<PointType>();

            for(size_t i = 0; i< Dimension; ++i){
                ptr_point->setDataById(i, temp_pt[i]);
            }
            ptr_vertex->setPtr_data(ptr_point);

            vec_ptr_vertex.push_back(ptr_vertex);
        }
    }else
        return false;

    std::cout << ">>>> ( End ) Loading Position..." << std::endl;
    return true;
}

template<typename GeometricReaderType, size_t Dimension>
bool ArtProjectReader<GeometricReaderType, Dimension>::readFace() {

    using namespace geometric_tree;
    using PointType = typename GeometricReaderType::type::PointType;
    using PtrPointType = typename GeometricReaderType::type::PtrPointType;
    using VertexType = typename GeometricReaderType::type::VertexType_traits::VertexType;
    using PtrVertexType = typename GeometricReaderType::type::VertexType_traits::PtrVertexType;

    std::cout << ">>>> (Start) Loading Connectivity..." << std::endl;
    std::string file_name;
    file_name += art_project->getProjectGeometryPath();
    file_name += "connectivity.txt";

    //std::cout << file_name << std::endl;

    std::ifstream fs;
    std::string bufferLine;
    size_t bufferId;

    auto &vec_ptr_vertex = ptr_data->getTotalVec_PtrVertex();
    auto &vec_ptr_face = ptr_data->getTotalVec_PtrFace();
    Face_Factory<PointType> face_factory;

    fs.open(file_name);
    if(fs.is_open()){

        while(getline( fs, bufferLine )) {
            std::stringstream w(bufferLine);

            face_factory.clearVertex();

            while(w >> bufferId){
                face_factory.addVertex(vec_ptr_vertex[bufferId]);
            }

            vec_ptr_face.push_back(
                    face_factory.create()
            );
        }
    }else
        return false;


    std::cout << ">>>> ( End ) Loading Connectivity..." << std::endl;
    return true;
}

template<typename GeometricReaderType, size_t Dimension>
bool ArtProjectReader<GeometricReaderType, Dimension>::readCell() {
    using namespace geometric_tree;
    using PointType = typename GeometricReaderType::type::PointType;
    using PtrPointType = typename GeometricReaderType::type::PtrPointType;
    using VertexType = typename GeometricReaderType::type::VertexType_traits::VertexType;
    using PtrVertexType = typename GeometricReaderType::type::VertexType_traits::PtrVertexType;

    std::cout << ">>>> (Start) Loading Connectivity..." << std::endl;
    std::string file_name;
    file_name += art_project->getProjectGeometryPath();
    file_name += "connectivity.txt";

    //std::cout << file_name << std::endl;

    std::ifstream fs;
    std::string bufferLine;
    size_t bufferId;

    auto &vec_ptr_vertex = ptr_data->getTotalVec_PtrVertex();
    auto &vec_ptr_cell = ptr_data->getTotalVec_PtrCell();
    Cell_Factory<PointType> cell_factory;

    fs.open(file_name);
    if(fs.is_open()){

        while(getline( fs, bufferLine )) {
            std::stringstream w(bufferLine);

            cell_factory.clearVertex();

            while(w >> bufferId){
                cell_factory.addVertex(vec_ptr_vertex[bufferId]);
            }

            vec_ptr_cell.push_back(
                    cell_factory.create()
            );
        }
    }else
        return false;

    return true;
}

