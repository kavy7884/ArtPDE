//
// Created by Chingkai Chou on 3/8/18.
//

#ifndef ARTCFD_SHAPEFUNCTION_HPP
#define ARTCFD_SHAPEFUNCTION_HPP

#include "Eigen/Dense"

enum class ShapeFuncType{ Lagrange, Spline};
enum class ShapeFuncOrder{ Linear, Quadratic};

class ShapeFunction{
public:
    ShapeFunction() {}
    Eigen::RowVector4d ShapeFunctionQ4_N(double &xi, double &eta);
    Eigen::Matrix<double, 2, 4> ShapeFunctionQ4_dN_dxi(double &xi, double &eta);
};


class ShapeFunctionLagrange{

};

Eigen::RowVector4d ShapeFunction::ShapeFunctionQ4_N(double &xi, double &eta) {
    Eigen::RowVector4d shape;
    shape(0) = (1.0-xi)*(1.0-eta) / 4.0;
    shape(1) = (1.0+xi)*(1.0-eta) / 4.0;
    shape(2) = (1.0+xi)*(1.0+eta) / 4.0;
    shape(3) = (1.0-xi)*(1.0+eta) / 4.0;
    return shape;
}

Eigen::Matrix<double, 2, 4> ShapeFunction::ShapeFunctionQ4_dN_dxi(double &xi, double &eta) {
    Eigen::Matrix<double, 2, 4> naturalDerivatives;
    naturalDerivatives(0,0) = -(1.0-eta)/4.0;
    naturalDerivatives(0,1) = (1.0-eta)/4.0;
    naturalDerivatives(0,2) = (1.0+eta)/4.0;
    naturalDerivatives(0,3) = -(1.0+eta)/4.0;
    naturalDerivatives(1,0) = -(1.0-xi)/4.0;
    naturalDerivatives(1,1) = -(1.0+xi)/4.0;
    naturalDerivatives(1,2) = (1.0+xi)/4.0;
    naturalDerivatives(1,3) = (1.0-xi)/4.0;

    return naturalDerivatives;
}

#endif //ARTCFD_SHAPEFUNCTION_HPP
