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

        template <unsigned Dimension>
        class PositionVector{
            typedef std::array<double,Dimension> ArrayType;
            typedef std::shared_ptr<ArrayType> PtrArrayType;
        public:
            PositionVector<Dimension>() {
                newData();
            }

            PositionVector<Dimension>(std::initializer_list<double> &input_list) {
                newData();
                addDataByList(input_list);
            }

            friend std::ostream &operator<<(std::ostream &os, const PositionVector<Dimension> &vector) {
                os << "[ ";
                for (size_t i = 0; i < vector.data->size() - 1; ++i) {
                    os << vector.data->at(i) << " ";
                }
                os << vector.data->at(vector.data->size() - 1);
                os << " ] ";
                return os;
            }

        protected:
            PtrArrayType data{nullptr};

            void newData(){
                this->data = nullptr;
                this->data = std::make_shared<ArrayType>();
            }

            void addDataByList(std::initializer_list<double> &v){
                assert(v.size() <= Dimension);
                std::copy(v.begin(), v.end(), data->begin());
            }
        };

        template <unsigned Dimension, typename T> class CartesianReadable;
        template <unsigned Dimension, typename T> class CartesianWritable;

        #include "./position_vector_cartesian_detail.cpp"

    }
}




#endif //ARTPDE_POSITION_VECTOR_HPP
