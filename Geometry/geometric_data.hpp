//
// Created by Chingkai Chou on 6/3/18.
//

#ifndef ARTPDE_KAVY_GEOMETRICDATA_HPP
#define ARTPDE_KAVY_GEOMETRICDATA_HPP

#include "inc/GeometricData/geometric_data_mesh_type.hpp"

namespace art_pde{ namespace geometry {
        namespace mesh_type {
            namespace Dim2 {

                using namespace geometric_data;

                template <typename BasicPointType>
                class GeometricData
                        :public MeshTypeData_Vertex<BasicPointType>,
                         public MeshTypeData_Edge<BasicPointType>,
                         public MeshTypeData_Face<BasicPointType>{
                public:
                    struct type{
                        using PointType = BasicPointType;
                        using PtrPointType = std::shared_ptr<BasicPointType>;
                        using VertexType_traits = typename MeshTypeData_Vertex<BasicPointType>::type;
                        using EdgeType_traits = typename MeshTypeData_Edge<BasicPointType>::type;
                        using FaceType_traits = typename MeshTypeData_Face<BasicPointType>::type;
                    };

                    GeometricData():
                            MeshTypeData_Vertex<BasicPointType>(),
                            MeshTypeData_Edge<BasicPointType>(),
                            MeshTypeData_Face<BasicPointType>(){
                        //std::cout << "Dim2::GeometricData" << std::endl;
                    }
                };

            }

            namespace Dim3 {

                using namespace geometric_data;

                template <typename BasicPointType>
                class GeometricData
                        :public MeshTypeData_Vertex<BasicPointType>,
                         public MeshTypeData_Edge<BasicPointType>,
                         public MeshTypeData_Face<BasicPointType>,
                         public MeshTypeData_Cell<BasicPointType>{
                public:
                    struct type{
                        using PointType = BasicPointType;
                        using PtrPointType = std::shared_ptr<BasicPointType>;
                        using VertexType_traits = typename MeshTypeData_Vertex<BasicPointType>::type;
                        using EdgeType_traits = typename MeshTypeData_Edge<BasicPointType>::type;
                        using FaceType_traits = typename MeshTypeData_Face<BasicPointType>::type;
                        using CellType_traits = typename MeshTypeData_Cell<BasicPointType>::type;
                    };
                    GeometricData():
                            MeshTypeData_Vertex<BasicPointType>(),
                            MeshTypeData_Edge<BasicPointType>(),
                            MeshTypeData_Face<BasicPointType>(),
                            MeshTypeData_Cell<BasicPointType>(){
                        //std::cout << "Dim3::GeometricData" << std::endl;
                    }
                };
            }
        }
}}

#endif //ARTPDE_KAVY_GEOMETRICDATA_HPP
