//
// Created by Chingkai Chou on 6/7/18.
//

#ifndef ARTPDE_KAVY_FUNCTION_SPACE_HPP
#define ARTPDE_KAVY_FUNCTION_SPACE_HPP

#include <memory>
#include "basis_function.hpp"
namespace art_pde{

    template <typename BasisFucntionType>
    class FunctionSpace : public BasisFucntionType{
        struct FunctionSpaceDummy {};
    public:
        template <typename ...Args>
        static std::shared_ptr<FunctionSpace> create(Args&& ...args){
            return make_shared_FunctionSpace(std::forward<Args>(args)...);
        }

        template <typename ...Args>
        explicit FunctionSpace(FunctionSpaceDummy, Args&& ...args): BasisFucntionType(std::forward<Args>(args)...) {}


    protected:
        template <typename ...Args>
        static std::shared_ptr<FunctionSpace> make_shared_FunctionSpace(Args&& ...args) {
            return std::make_shared<FunctionSpace>(FunctionSpaceDummy(), std::forward<Args>(args)...);
        }


    };
}


#endif //ARTPDE_KAVY_FUNCTION_SPACE_HPP
