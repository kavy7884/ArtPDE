//
// Created by Chingkai Chou on 5/2/18.
//

#ifndef ARTCFD_POINT_HPP
#define ARTCFD_POINT_HPP

#include <ostream>
#include <string>
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
        using PointDim = Dim1D;
        Point() {}
        Point(const Point& point): x(point.getX()){}
        Point(double x) : x(x) {}

        double getX() const {
            return x;
        }

        void setX(double x) {
            Point::x = x;
        }

        inline Point& operator+=(const Point& rhs){
            x += rhs.getX();
            return *this;
        }

        inline Point& operator/=(const double & double_rhs){
            x /= double_rhs;
            return *this;
        }

        friend inline Point operator+(Point lhs, const Point& rhs){
            lhs += rhs;
            return lhs;
        }

        friend std::ostream &operator<<(std::ostream &os, const Point<Dim1D, CartesianCoordinate> &point) {
            os << " [ " << point.x << " ] ";
            return os;
        }

    private:
        double x{0.0};
    };

    // Point specialized in 2D Cartesian coordinate case
    template <>
    class Point<Dim2D, CartesianCoordinate>{
    public:
        using PointDim = Dim2D;
        Point() {}
        Point(const Point& point): x(point.getX()), y(point.getY()){}
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

        inline Point& operator+=(const Point& rhs){
            x += rhs.getX();
            y += rhs.getY();
            return *this;
        }

        inline Point& operator/=(const double & double_rhs){
            x /= double_rhs;
            y /= double_rhs;
            return *this;
        }

        friend inline Point operator+(Point lhs, const Point& rhs){
            lhs += rhs;
            return lhs;
        }

        friend std::ostream &operator<<(std::ostream &os, const Point<Dim2D, CartesianCoordinate> &point) {
            os << " [ " << point.x << ", "<< point.y << " ] ";
            return os;
        }

    private:
        double x{0.0}, y{0.0};
    };

    // Point specialized in 3D Cartesian coordinate case
    template <>
    class Point<Dim3D, CartesianCoordinate>{
    public:
        using PointDim = Dim3D;
        Point() {}
        Point(const Point& point): x(point.getX()), y(point.getY()), z(point.getZ()){}
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

        inline Point& operator+=(const Point& rhs){
            x += rhs.getX();
            y += rhs.getY();
            z += rhs.getZ();
            return *this;
        }

        inline Point& operator/=(const double & double_rhs){
            x /= double_rhs;
            y /= double_rhs;
            z /= double_rhs;
            return *this;
        }

        friend inline Point operator+(Point lhs, const Point& rhs){
            lhs += rhs;
            return lhs;
        }

        friend std::ostream &operator<<(std::ostream &os, const Point<Dim3D, CartesianCoordinate> &point) {
            os << " [ " << point.x << ", "<< point.y << ", "<< point.z << " ] ";
            return os;
        }

    private:
        double x{0.0}, y{0.0}, z{0.0};
    };

}

#endif //ARTCFD_POINT_HPP
