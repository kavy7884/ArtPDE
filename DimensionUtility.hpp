//
// Created by Chingkai Chou on 3/2/18.
//

#ifndef ARTCFD_DIMENSIONUTILITY_HPP
#define ARTCFD_DIMENSIONUTILITY_HPP
class Dimension{
public:
    static const size_t Dim{0};
};
class Dim1D : public Dimension{
public:
    static const size_t Dim{1};
};
class Dim2D : public Dimension{
public:
    static const size_t Dim{2};
};
class Dim3D : public Dimension{
public:
    static const size_t Dim{3};
};
#endif //ARTCFD_DIMENSIONUTILITY_HPP
