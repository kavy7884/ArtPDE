//
// Created by Chingkai Chou on 5/17/18.
//

#ifndef ARTPDE_MARVIN_POINT3_HPP
#define ARTPDE_MARVIN_POINT3_HPP

#include <ostream>

class Point3{
public:
    Point3(double x, double y, double z) : x(x), y(y), z(z) {}

    double x, y, z;

    friend std::ostream &operator<<(std::ostream &os, const Point3 &point3) {
        os << "[ " << point3.x << " , " << point3.y << " , " << point3.z << "] ";
        return os;
    }
};

#endif //ARTPDE_MARVIN_POINT3_HPP
