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
                class ArtProjectReader{
                public:
                    ArtProjectReader(GeometricReaderType *input_ptr_data)
                            :ptr_data(input_ptr_data){}

                    void setArtProject(const std::shared_ptr<project::ArtProject> &input_ptr_proj){
                        this->art_project = input_ptr_proj;
                    }

                    bool read();

                private:
                    bool readPosition();
                    bool readFace();

                    GeometricReaderType *ptr_data;
                    std::shared_ptr<project::ArtProject> art_project{nullptr};
                };



#include "../../src/geometric_reader_art_project_impl.cpp"
            }
        }
}}

#endif //ARTPDE_KAVY_GEOMETRIC_READER_ART_PROJECT_HPP
