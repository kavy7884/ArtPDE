//
// Created by Chingkai Chou on 5/7/18.
//

#ifndef ARTCFD_GEOMETRY_DATA_READER_ARTPDE_HPP
#define ARTCFD_GEOMETRY_DATA_READER_ARTPDE_HPP

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "geometry_data_reader.hpp"
#include "project_uility.hpp"

namespace art_pde {

    template <class GeometryDataType, typename GeometryDescription, typename Dimension, typename CoordinateBasis>
    class GeometryDataReaderArtPDE_Implementation;

    template<class GeometryDataType>
    class GeometryDataReaderArtPDE : public GeometryDataReader<GeometryDataType, ReadArtPDE> {
    public:

        using ImplementClass = GeometryDataReaderArtPDE_Implementation<GeometryDataType,
                typename GeometryDataType::Type::GeoDescribeType,
                typename GeometryDataType::Type::GeoDimType,
                typename GeometryDataType::Type::GeoCoordType>;

        GeometryDataReaderArtPDE() : GeometryDataReader<GeometryDataType, ReadArtPDE>(){}

        void read(const std::shared_ptr<ArtProject>& ptr_artpde_project){

            std::cout << ">> (Start) Reading ArtPDE project format files..." << std::endl;

            ImplementClass im(this);
            im.setPtr_artpde_project(ptr_artpde_project);
            if(im.load()){
                this->is_read_success = true;
            }else{
                this->is_read_success = false;
            }

            std::cout << ">> ( End ) Reading ArtPDE project format files..." << std::endl;
        }

    private:
        friend ImplementClass;

        typename GeometryDataType::Type::VecPtrGeoVertexType& getTotalVertex(){ return this->total_vertex; }
        typename GeometryDataType::Type::VecPtrGeoCellType& getTotalCell(){ return this->total_cell; }
    };

    template <class GeometryDataType, typename GeometryDescription, typename Dimension, typename CoordinateBasis>
    class GeometryDataReaderArtPDE_Implementation;

    template <class GeometryDataType>
    class GeometryDataReaderArtPDE_Implementation<GeometryDataType, MeshTypeMethod, Dim2D, CartesianCoordinate>{
    public:
        using ArtProjectType = ArtProject;
        using PtrArtProjectType = std::shared_ptr<ArtProjectType>;

        GeometryDataReaderArtPDE_Implementation(GeometryDataReaderArtPDE<GeometryDataType> *loader) : loader(loader) {}

        void setPtr_artpde_project(const PtrArtProjectType &ptr_artpde_project) {
            this->ptr_artpde_project = ptr_artpde_project;
        }

        bool load(){
            loader->test = 3.0;

            if(!load_xNode()) return false;
            if(!load_cElement30()) return false;

            return true;
        }

    private:
        bool load_xNode();
        bool load_cElement30();

        GeometryDataReaderArtPDE<GeometryDataType> *loader{nullptr};
        PtrArtProjectType ptr_artpde_project{nullptr};
    };

    template<class GeometryDataType>
    bool
    GeometryDataReaderArtPDE_Implementation<GeometryDataType, MeshTypeMethod, Dim2D, CartesianCoordinate>::load_xNode() {
        std::cout << ">>>> (Start) Loading xNode..." << std::endl;
        std::string file_name;
        file_name += ptr_artpde_project->getProjectGeometryPath();
        file_name += "xNode.txt";

//        std::cout << file_name << std::endl;

        std::ifstream fs;
        std::string bufferLine;
        double pt[Dim2D::k_NumDim];
        typename GeometryDataType::Type::VecPtrGeoVertexType& total_vertex = loader->getTotalVertex();

        fs.open(file_name);
        if(fs.is_open()){
            while(getline( fs, bufferLine )) {
                std::stringstream w(bufferLine);
                w >> pt[0];
                w >> pt[1];
                total_vertex.push_back(
                        std::make_shared<typename GeometryDataType::Type::GeoVertexType>(
                                typename GeometryDataType::Type::GeoPointType(pt[0], pt[1])
                        )
                );
            }
        }else
            return false;

        std::cout << ">>>> ( End ) Loading xNode..." << std::endl;
        return true;
    }

    template<class GeometryDataType>
    bool
    GeometryDataReaderArtPDE_Implementation<GeometryDataType, MeshTypeMethod, Dim2D, CartesianCoordinate>::load_cElement30() {
        std::cout << ">>>> (Start) Loading xElement30..." << std::endl;
        std::string file_name;
        file_name += ptr_artpde_project->getProjectGeometryPath();
        file_name += "cElement30.txt";

        //std::cout << file_name << std::endl;

        std::ifstream fs;
        std::string bufferLine;
        size_t bufferId;
        CellFactory<typename GeometryDataType::Type::GeoPointType, typename GeometryDataType::Type::GeoDimType> cell_factory;
        typename GeometryDataType::Type::VecPtrGeoVertexType& total_vertex = loader->getTotalVertex();
        typename GeometryDataType::Type::VecPtrGeoCellType& total_cell = loader->getTotalCell();

        fs.open(file_name);
        if(fs.is_open()){
            while(getline( fs, bufferLine )) {
                std::stringstream w(bufferLine);

                cell_factory.clearVertex();
                while(w >> bufferId){
                    cell_factory.addVertex(total_vertex[bufferId]);
                }

                total_cell.push_back(
                        cell_factory.create()
                );
            }
        }else
            return false;

        std::cout << ">>>> ( End ) Loading xElement30..." << std::endl;
        return true;
    }

    template <class GeometryDataType>
    class GeometryDataReaderArtPDE_Implementation<GeometryDataType, MeshTypeMethod, Dim3D, CartesianCoordinate>{
    public:
        using ArtProjectType = ArtProject;
        using PtrArtProjectType = std::shared_ptr<ArtProjectType>;

        GeometryDataReaderArtPDE_Implementation(GeometryDataReaderArtPDE<GeometryDataType> *loader) : loader(loader) {}

        void setPtr_artpde_project(const PtrArtProjectType &ptr_artpde_project) {
            this->ptr_artpde_project = ptr_artpde_project;
        }

        bool load(){

            loader->test = 3.0;

            if(!load_xNode()) return false;
            if(!load_cElement30()) return false;

            return true;
        }

    private:
        bool load_xNode();
        bool load_cElement30();


        GeometryDataReaderArtPDE<GeometryDataType> *loader{nullptr};
        PtrArtProjectType ptr_artpde_project{nullptr};
    };

    template<class GeometryDataType>
    bool
    GeometryDataReaderArtPDE_Implementation<GeometryDataType, MeshTypeMethod, Dim3D, CartesianCoordinate>::load_xNode() {
        std::cout << ">>>> (Start) Loading xNode..." << std::endl;
        std::string file_name;
        file_name += ptr_artpde_project->getProjectGeometryPath();
        file_name += "xNode.txt";

//        std::cout << file_name << std::endl;

        std::ifstream fs;
        std::string bufferLine;
        double pt[Dim3D::k_NumDim];
        typename GeometryDataType::Type::VecPtrGeoVertexType& total_vertex = loader->getTotalVertex();

        fs.open(file_name);
        if(fs.is_open()){
            while(getline( fs, bufferLine )) {
                std::stringstream w(bufferLine);
                w >> pt[0];
                w >> pt[1];
                w >> pt[2];
                total_vertex.push_back(
                        std::make_shared<typename GeometryDataType::Type::GeoVertexType>(
                                typename GeometryDataType::Type::GeoPointType(pt[0], pt[1], pt[3])
                        )
                );
            }
        }else
            return false;

        std::cout << ">>>> ( End ) Loading xNode..." << std::endl;
        return true;
    }

    template<class GeometryDataType>
    bool
    GeometryDataReaderArtPDE_Implementation<GeometryDataType, MeshTypeMethod, Dim3D, CartesianCoordinate>::load_cElement30() {
        std::cout << ">>>> (Start) Loading xElement30..." << std::endl;
        std::string file_name;
        file_name += ptr_artpde_project->getProjectGeometryPath();
        file_name += "cElement30.txt";

        //std::cout << file_name << std::endl;

        std::ifstream fs;
        std::string bufferLine;
        size_t bufferId;
        CellFactory<typename GeometryDataType::Type::GeoPointType, typename GeometryDataType::Type::GeoDimType> cell_factory;
        typename GeometryDataType::Type::VecPtrGeoVertexType& total_vertex = loader->getTotalVertex();
        typename GeometryDataType::Type::VecPtrGeoCellType& total_cell = loader->getTotalCell();

        fs.open(file_name);
        if(fs.is_open()){
            while(getline( fs, bufferLine )) {
                std::stringstream w(bufferLine);

                cell_factory.clearVertex();
                while(w >> bufferId){
                    cell_factory.addVertex(total_vertex[bufferId]);
                }

                total_cell.push_back(
                        cell_factory.create()
                );
            }
        }else
            return false;

        std::cout << ">>>> ( End ) Loading xElement30..." << std::endl;
        return true;
    }

}
#endif //ARTCFD_GEOMETRY_DATA_READER_ARTPDE_HPP
