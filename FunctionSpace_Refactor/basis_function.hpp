//
// Created by Chingkai Chou on 6/7/18.
//

#ifndef ARTPDE_KAVY_BASIS_FUNCTION_HPP
#define ARTPDE_KAVY_BASIS_FUNCTION_HPP

namespace art_pde{ namespace function_space{
        namespace isoparametric{

            enum class BasisOrder{ Liner, Quadratic};

            template <typename GeometryType>
            class BasisFunction{
            public:
                BasisFunction(const std::shared_ptr<GeometryType> &geo): geo(geo){
                    std::cout << "isoparametric::BasisFunction" << std::endl;
                    std::cout << geo << std::endl;
                }

                BasisFunction(const std::shared_ptr<GeometryType> &geo, BasisOrder basis_order)
                        : geo(geo), basis_order(basis_order){
                    std::cout << "isoparametric::BasisFunction (with order) " << std::endl;
                    std::cout << geo << std::endl;
                }

                void testBasisFunc(){
                    std::cout << "Basis func test." << std::endl;
                }

            protected:
                std::shared_ptr<GeometryType> geo{nullptr};
                BasisOrder basis_order{BasisOrder::Liner};
            };
        }
}}

#endif //ARTPDE_KAVY_BASIS_FUNCTION_HPP
