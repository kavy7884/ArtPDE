//
// Created by Chingkai Chou on 5/2/18.
//

#include <iostream>
//#include "Geometry_New/GeoTest2D.hpp"
//#include "Geometry_New/GeoTest3D.hpp"
//
//#include "Geometry_New/geo_data.hpp"

#include "BasicData/position_vector.hpp"

int main() {

    using namespace art_pde::PositionVector;
    const unsigned Dim = 3;

    auto ptr_view_pt = PositionVector<Dim>::createViewPoint({1.0,2.0,3.0});
    decltype(PositionVector<Dim>::createViewPoint()) ptr_view_pt_1;
    ptr_view_pt_1 = ptr_view_pt; //shallow copy

    std::cout << *ptr_view_pt <<  std::endl;
    std::cout << *ptr_view_pt_1 <<  std::endl;

    auto ptr_compute_pt = PositionVector<Dim>::createComputePoint();
    auto ptr_compute_pt_1 = PositionVector<Dim>::createComputePoint();
    auto ptr_compute_pt_2 = PositionVector<Dim>::createComputePoint();

    ptr_compute_pt_1 = ptr_compute_pt;//shallow copy

    std::cout << *ptr_compute_pt <<  std::endl;
    std::cout << *ptr_compute_pt_1 <<  std::endl;

    ptr_compute_pt_1->setX(2.0); ptr_compute_pt_1->setY(3.0);

    std::cout << *ptr_compute_pt <<  std::endl;
    std::cout << *ptr_compute_pt_1 <<  std::endl;

    (*ptr_compute_pt_2) = (*ptr_compute_pt); //deep copy

    ptr_compute_pt_2->setX(1.0); ptr_compute_pt_2->setY(1.0); ptr_compute_pt_2->setZ(1.0);

    std::cout << *ptr_compute_pt <<  std::endl;
    std::cout << *ptr_compute_pt_2 <<  std::endl;

    auto ptr_view_pt_2 = PositionVector<Dim>::createViewPoint({9.0, 9.0, 9.0});
    (*ptr_compute_pt_2) = (*ptr_view_pt_2);
    std::cout << *ptr_compute_pt_2 <<  std::endl;

    return 0;
}