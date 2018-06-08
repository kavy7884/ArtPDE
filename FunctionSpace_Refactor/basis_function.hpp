//
// Created by Chingkai Chou on 6/7/18.
//

#ifndef ARTPDE_KAVY_BASIS_FUNCTION_HPP
#define ARTPDE_KAVY_BASIS_FUNCTION_HPP

#include <array>
#include "../BasicData/position_vector.hpp"

namespace art_pde{ namespace function_space{
        namespace isoparametric{

            enum class BasisOrder{ Liner, Quadratic};
            std::ostream &operator<<(std::ostream &os, const BasisOrder &basis_order) {
                switch(basis_order) {
                    case BasisOrder::Liner    : os << "Liner"; break;
                    case BasisOrder::Quadratic   : os << "Quadratic"; break;
                    default : os.setstate(std::ios_base::failbit);
                }
                return os;
            }

            template <typename GeometryType>
            class BasisFunction{
            public:
                struct type{

                };

                BasisFunction(const std::shared_ptr<GeometryType> &geo): geo(geo){
                    std::cout << "isoparametric::BasisFunction" << std::endl;
                    std::cout << "Geometry address is: " << geo << std::endl;
                }

                BasisFunction(const std::shared_ptr<GeometryType> &geo, BasisOrder basis_order)
                        : geo(geo), basis_order(basis_order){
                    std::cout << "isoparametric::BasisFunction (with order) " << std::endl;
                    std::cout << "Geometry address is: " << geo << std::endl;
                }

                void testBasisFunc(){
                    std::cout << "Basis func test." << std::endl;
                    std::cout << "Geometry dimension is: " << geo->getDimension() << std::endl;
                    std::cout << "Basis function order is: " << basis_order << std::endl;

                    std::cout << "new a point" << std::endl;
                    using namespace position_vector;
                    auto ptr_point = PositionVector<GeometryType::Dimension>::createComputePoint();
                    std::cout << *ptr_point << std::endl;

                }

            protected:
                std::shared_ptr<GeometryType> geo{nullptr};
                BasisOrder basis_order{BasisOrder::Liner};
            };
        }
}}

#endif //ARTPDE_KAVY_BASIS_FUNCTION_HPP
