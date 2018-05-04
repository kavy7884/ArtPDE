//
// Created by Chingkai Chou on 5/2/18.
//

#ifndef ARTCFD_POINT_HPP
#define ARTCFD_POINT_HPP

#include "dimension_utility.hpp"
#include "coordinate_basis_utility.hpp"

namespace art_pde {

    // Point abstract define
    template <typename Dimension, typename PointCoordinateBasis>
    class Point;

    // Point partial define in Cartesian coordinate case
    template <typename Dimension>
    class Point<Dimension, CartesianCoordinate>;

    // Point specialized in 1D Cartesian coordinate case
    template <>
    class Point<Dim1D, CartesianCoordinate>{
    public:
        Point() {}
        Point(double x) : x(x) {}

        double getX() const {
            return x;
        }

        void setX(double x) {
            Point::x = x;
        }

    private:
        double x{0.0};
    };

    // Point specialized in 2D Cartesian coordinate case
    template <>
    class Point<Dim2D, CartesianCoordinate>{
    public:
        Point() {}
        Point(double x, double y) : x(x), y(y) {}

        double getX() const {
            return x;
        }

        void setX(double x) {
            Point::x = x;
        }

        double getY() const {
            return y;
        }

        void setY(double y) {
            Point::y = y;
        }

    private:
        double x{0.0}, y{0.0};
    };

    // Point specialized in 3D Cartesian coordinate case
    template <>
    class Point<Dim3D, CartesianCoordinate>{
    public:
        Point() {}
        Point(double x, double y, double z) : x(x), y(y), z(z) {}

        double getX() const {
            return x;
        }

        void setX(double x) {
            Point::x = x;
        }

        double getY() const {
            return y;
        }

        void setY(double y) {
            Point::y = y;
        }

        double getZ() const {
            return z;
        }

        void setZ(double z) {
            Point::z = z;
        }

    private:
        double x{0.0}, y{0.0}, z{0.0};
    };

}

#endif //ARTCFD_POINT_HPP
