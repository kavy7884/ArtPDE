//
// Created by Chingkai Chou on 5/25/18.
//

#ifndef ARTPDE_POSITION_VECTOR_HPP
#define ARTPDE_POSITION_VECTOR_HPP

#include <iostream>
#include <array>
#include <memory>
#include <algorithm>
#include <assert.h>

namespace art_pde{ namespace position_vector{

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
            PointData(const std::initializer_list<double> &v){ this->newData(); this->addDataByList(v);}
            void newData();
            void addDataByList(const std::initializer_list<double> &v);
            const double &getDataById(const size_t &id) const;
            void setDataById(const size_t &id, const double &value);
            void setDataById(const size_t &id, double &value);

            PtrArrayType data{nullptr}; //real data storage
        };
#include "./src/position_vector_pointdata_impl.cpp"
// -------- PointData <End> -----------

// -------- CartesianAPI <Start> -----------
#include "./src/position_vector_cartesian_impl.cpp"
    // -------- CartesianAPI <End> -----------

    // -------- Real apply class <Start> -----------
    template <size_t Dimension>
    class ViewPositionVector:
            public CartesianAPI<Dimension, false>{
    public:
        ViewPositionVector(): PointData<Dimension>(){}
        ViewPositionVector(const std::initializer_list<double> &v): PointData<Dimension>(v){}

        const double &getDataById(const size_t id) const{
            return PointData<Dimension>::getDataById(id);
        }

    };

    template <size_t Dimension>
    class ComputePositionVector:
            public CartesianAPI<Dimension, true>{
    public:
        ComputePositionVector(): PointData<Dimension>(){}
        ComputePositionVector(const std::initializer_list<double> &v): PointData<Dimension>(v){}

        ComputePositionVector<Dimension>& operator=(const ViewPositionVector<Dimension>& other)
        {
            PointData<Dimension>::operator=(other);
            return *this;
        }

        const double &getDataById(const size_t id) const{
            return PointData<Dimension>::getDataById(id);
        }

        void setDataById(const size_t id, const double &value){
            return PointData<Dimension>::setDataById(id, value);
        }

        void setDataById(const size_t id, double &value){
            return PointData<Dimension>::setDataById(id, value);
        }

    };
    // -------- Real apply class <End> -----------

    template <size_t Dimension>
    class PositionVector{
    public:
        static std::shared_ptr<ViewPositionVector<Dimension>> createViewPoint(){
            std::shared_ptr<ViewPositionVector<Dimension>> re_ptr;
            re_ptr = std::make_shared<ViewPositionVector<Dimension>>();
            return re_ptr;
        }
        static std::shared_ptr<ViewPositionVector<Dimension>> createViewPoint(const std::initializer_list<double> &input_list){
            std::shared_ptr<ViewPositionVector<Dimension>> re_ptr;
            re_ptr = std::make_shared<ViewPositionVector<Dimension>>(input_list);
            return re_ptr;
        }

        static std::shared_ptr<ComputePositionVector<Dimension>> createComputePoint(){
            std::shared_ptr<ComputePositionVector<Dimension>> re_ptr;
            re_ptr = std::make_shared<ComputePositionVector<Dimension>>();
            return re_ptr;
        }
        static std::shared_ptr<ComputePositionVector<Dimension>> createComputePoint(const std::initializer_list<double> &input_list){
            std::shared_ptr<ComputePositionVector<Dimension>> re_ptr;
            re_ptr = std::make_shared<ComputePositionVector<Dimension>>(input_list);
            return re_ptr;
        }

        PositionVector() = delete;
    };
}}




#endif //ARTPDE_POSITION_VECTOR_HPP
