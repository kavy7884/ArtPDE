//
// Created by Chingkai Chou on 6/3/18.
//

#ifndef ARTPDE_KAVY_GEOMETRIC_READER_HPP
#define ARTPDE_KAVY_GEOMETRIC_READER_HPP

#include "inc/GeometricReader/geometric_reader_art_project.hpp"

namespace art_pde{ namespace geometry {
        namespace mesh_type {

            namespace Dim2 {

                template <typename GeometricDataType>
                class GeometricReader: public virtual GeometricDataType{
                public:
                    GeometricReader(): GeometricDataType(){
                        //std::cout << "Dim2::GeometricReader" << std::endl;
                    }

                    bool read(const std::shared_ptr<project::ArtProject> &input_ptr_proj){
                        //std::cout << "Dim2::GeometricReader->read()" << std::endl;
                        geometric_reader::ArtProjectReader<GeometricReader<GeometricDataType>, 2>
                                reader(input_ptr_proj, this);
                        return reader.read();
                    }

                protected:
                    template <typename T, size_t > friend class geometric_reader::ArtProjectReader;

                };
            }

            namespace Dim3 {

                template <typename GeometricDataType>
                class GeometricReader: public virtual GeometricDataType{
                public:
                    GeometricReader(): GeometricDataType(){
                        //std::cout << "Dim3::GeometricReader" << std::endl;
                    }

                    bool read(const std::shared_ptr<project::ArtProject> &input_ptr_proj){
                        //std::cout << "Dim3::GeometricReader->read()" << std::endl;
                        geometric_reader::ArtProjectReader<GeometricReader<GeometricDataType>, 3>
                                reader(input_ptr_proj, this);
                        return reader.read();
                    }

                protected:
                    template <typename T, size_t > friend class geometric_reader::ArtProjectReader;

                };
            }

        }
}}

#endif //ARTPDE_KAVY_GEOMETRIC_READER_HPP
