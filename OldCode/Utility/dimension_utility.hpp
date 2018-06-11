//
// Created by Chingkai Chou on 5/3/18.
//

#ifndef ARTCFD_DIMENSION_UTILITY_HPP
#define ARTCFD_DIMENSION_UTILITY_HPP

#include <cstdlib>
namespace art_pde {
    class Dimension{
        const static size_t k_NumDim = 0;
    };
    class Dim1D : public Dimension{
    public:
        const static size_t k_NumDim = 1;
    };
    class Dim2D : public Dimension{
    public:
        const static size_t k_NumDim = 2;
    };
    class Dim3D : public Dimension{
    public:
        const static size_t k_NumDim = 3;
    };
}

#endif //ARTCFD_DIMENSION_UTILITY_HPP
