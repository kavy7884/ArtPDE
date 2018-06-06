//
// Created by Chingkai Chou on 6/2/18.
//

#ifndef ARTPDE_KAVY_GEOMETRIC_READER_ART_PROJECT_HPP
#define ARTPDE_KAVY_GEOMETRIC_READER_ART_PROJECT_HPP

#include <memory>
#include <fstream>
#include <sstream>
#include "../../../Project/art_project.hpp"
#include "../GeometricData/geometric_data_mesh_type.hpp"

namespace art_pde{ namespace geometry {
        namespace mesh_type {
            namespace geometric_reader{

                template <typename GeometricReaderType, size_t Dimension>
                class ArtProjectReaderBase{
                public:
                    ArtProjectReaderBase(const std::shared_ptr<project::ArtProject> &input_ptr_proj)
                            : art_project(input_ptr_proj){};

                    virtual bool read() = 0;

                protected:
                    using VecPtrVertexType = typename GeometricReaderType::type::VertexType_traits::VecPtrVertexType;
                    bool readPosition(VecPtrVertexType &vec_ptr_vertex);

                    std::shared_ptr<project::ArtProject> art_project{nullptr};
                };

                template <typename GeometricReaderType, size_t Dimension> class ArtProjectReader;

                template <typename GeometricReaderType>
                class ArtProjectReader<GeometricReaderType, 2>
                        : public ArtProjectReaderBase<GeometricReaderType, 2>{
                public:

                    ArtProjectReader<GeometricReaderType, 2>(const std::shared_ptr<project::ArtProject> &input_ptr_proj,
                                                             GeometricReaderType *ptr_data) :
                            ArtProjectReaderBase<GeometricReaderType, 2>(input_ptr_proj), ptr_data(ptr_data) {}

                    bool read(){
                        if(!this->readPosition(ptr_data->getTotalVec_PtrVertex())) return false;
                        if(!this->readFace()) return false;
                        return true;
                    }

                private:
                    bool readFace();
//                    bool readCell();

                    GeometricReaderType *ptr_data;
                };


                template <typename GeometricReaderType>
                class ArtProjectReader<GeometricReaderType, 3>
                        : public ArtProjectReaderBase<GeometricReaderType, 3>{
                public:
                    ArtProjectReader<GeometricReaderType, 3>(const std::shared_ptr<project::ArtProject> &input_ptr_proj,
                                                             GeometricReaderType *ptr_data) :
                            ArtProjectReaderBase<GeometricReaderType, 3>(input_ptr_proj), ptr_data(ptr_data) {}

                    bool read(){
                        if(!this->readPosition(ptr_data->getTotalVec_PtrVertex())) return false;
                        if(!this->readCell()) return false;
                        return true;
                    }

                private:
                    bool readCell();

                    GeometricReaderType *ptr_data;
                };

                #include "../../src/geometric_reader_art_project_impl.cpp"
            }
        }
}}

#endif //ARTPDE_KAVY_GEOMETRIC_READER_ART_PROJECT_HPP
