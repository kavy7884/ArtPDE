//
// Created by Chingkai Chou on 5/25/18.
//

#ifndef ARTPDE_POSITION_VECTOR_HPP
#define ARTPDE_POSITION_VECTOR_HPP

#include <iostream>
#include <array>
#include <memory>
#include <algorithm>

namespace art_pde {
    namespace PositionVector{

        // -------- PointData <Start> -----------
        template <size_t Dimension>
        class PointData{
        public:
            typedef std::array<double,Dimension> ArrayType;
            typedef std::shared_ptr<ArrayType> PtrArrayType;

            PointData<Dimension>& operator=(const PointData<Dimension>& other);

            template <size_t Dimension_>
            friend std::ostream &operator<<(std::ostream &os, const PointData<Dimension_> &point_data);

        protected:
            PointData(){ this->newData(); }
            void newData();
            void addDataByList(std::initializer_list<double> &v);
            PtrArrayType data{nullptr};
        };
        #include "./src/position_vector_pointdata_impl.cpp"
        // -------- PointData <End> -----------

        // -------- CartesianAPI <Start> -----------
        template <size_t Dimension, bool Authority> class CartesianAPI;
        #include "./src/position_vector_cartesian_impl.cpp"
        // -------- CartesianAPI <End> -----------

        // -------- Real apply class <Start> -----------
        template <size_t Dimension>
        class ViewPositionVector:
                public CartesianAPI<Dimension, false>{
        public:
            ViewPositionVector(): PointData<Dimension>(){}
        };

        template <size_t Dimension>
        class ComputePositionVector:
                public CartesianAPI<Dimension, true>{
        public:
            ComputePositionVector(): PointData<Dimension>(){}

            ComputePositionVector<Dimension>& operator=(const ViewPositionVector<Dimension>& other)
            {
                PointData<Dimension>::operator=(other);
                return *this;
            }
        };
        // -------- Real apply class <End> -----------


    }
}




#endif //ARTPDE_POSITION_VECTOR_HPP
