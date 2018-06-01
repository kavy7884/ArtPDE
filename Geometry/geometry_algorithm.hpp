//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_ALGORITHM_HPP
#define ARTPDE_KAVY_GEOMETRY_ALGORITHM_HPP

#include "geometry_mesh_type_data.hpp"

namespace art_pde{ namespace geometry {

    template<typename DataType> class Algorithm_Base;

    template <size_t Dimension, template<size_t> class CalculatedPointType>
    class Algorithm_Base<mesh_type::Data<Dimension, CalculatedPointType>>:
            public virtual mesh_type::Data<Dimension, CalculatedPointType>{
    public:
        struct type{
            using DataType = mesh_type::Data<Dimension, CalculatedPointType>;
            using VecPtrVertex = typename DataType::type::VecPtrVertex;
            using VecPtrCell = typename DataType::type::VecPtrCell;
        };

        Algorithm_Base(): mesh_type::Data<Dimension, CalculatedPointType>(){
            std::cout << "Algorithm" << std::endl;
        }

        const size_t getTotalNum_Vertex() const{
            return this->vec_ptr_vertex.size();
        }

        virtual const size_t getTotalNum_Cell() const{
            return this->vec_ptr_cell.size();
        }

        const typename type::VecPtrVertex &getTotal_PtrVertex() const{
            return this->vec_ptr_vertex;
        }

        virtual const typename type::VecPtrCell &getTotal_PtrCell() const{
            return this->vec_ptr_cell;
        }

    };

    template <size_t Dimension, typename DataType> class Algorithm;

    template <typename DataType>
    class Algorithm<1, DataType> : public Algorithm_Base<DataType>{
    public:
        Algorithm<1, DataType>(): Algorithm_Base<DataType>(){}

        const size_t getTotalNum_Cell() const = delete;
        const typename DataType::type::VecPtrCell &getTotal_PtrCell() const = delete;

    };

    template <typename DataType>
    class Algorithm<2, DataType> : public Algorithm_Base<DataType>{
    public:
        Algorithm<2, DataType>(): Algorithm_Base<DataType>(){}

        const size_t getTotalNum_Cell() const = delete;
        const typename DataType::type::VecPtrCell &getTotal_PtrCell() const = delete;

    };

    template <typename DataType>
    class Algorithm<3, DataType> : public Algorithm_Base<DataType>{
    public:
        Algorithm<3, DataType>(): Algorithm_Base<DataType>(){}

    };
}}

#endif //ARTPDE_KAVY_GEOMETRY_ALGORITHM_HPP
