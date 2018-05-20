//
// Created by Chingkai Chou on 5/2/18.
//

#include "Geometry_New/GeoTest2D.hpp"

int main() {
    GeoTest2D test2D;

    std::cout << "<Start> Merge Edge " << std::endl;
    test2D.mergeEdge();
    std::cout << "< End > Merge Edge " << std::endl;
    test2D.checkOut();

    return 0;
}