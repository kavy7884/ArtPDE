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
    const unsigned Dim1 = 1;
    const unsigned Dim2 = 2;
    const unsigned Dim3 = 3;
////
//    using ViewPoint3 = CartesianReadable<Dim3,PositionVector<Dim3>>;
////    using CompPoint3 = CartesianWritable<Dim3,PositionVector<Dim3>>;
//////
////    ViewPoint3 view_pt3({1.0, 2.0, 3.0});
////    ViewPoint3 view_pt3_1(view_pt3);
////    std::cout << view_pt3_1 << std::endl;
//
////    CartesianReadable<Dim3,PositionVector<Dim3>> a({1.0, 2.0, 3.0});
////    CartesianReadable<Dim3,PositionVector<Dim3>> b(a);
////    std::cout << b << std::endl;
//
//    Test<Dim3> tt;
//    std::cout << tt.a <<  std::endl;
//    tt.haha1();
//    std::cout << tt.a <<  std::endl;
//    tt.haha2();
//    std::cout << tt.a <<  std::endl;

    ViewPositionVector<Dim3> view_pt;
    std::cout << view_pt<<  std::endl;

    ComputePositionVector<Dim3> comp_pt, comp_pt_1;
    comp_pt.setX(1.0);
    std::cout << comp_pt <<  std::endl;
    comp_pt_1 = comp_pt;
    std::cout << comp_pt_1 <<  std::endl;
    comp_pt_1 = view_pt;
    std::cout << comp_pt_1 <<  std::endl;



//    CompPoint3 comp_pt3;
//
//    std::cout << comp_pt3 << std::endl;

//    PositionVector<Dim3> Pt{1.0, 2.0, 3.0};
//    std::cout << Pt << std::endl;


//
//    using ComputePoint2 = CartesianWritable<Dim2,PositionVector<Dim2>>;
//    using ViewPoint2 = CartesianReadable<Dim2,PositionVector<Dim2>>;
//    ComputePoint2 compute_pt2_1, compute_pt2_2(2.0, 3.0);
//    compute_pt2_1.setX(1.0); compute_pt2_1.setY(1.0);
//    std::cout<< compute_pt2_1 << std::endl;
//    std::cout<< compute_pt2_2 << std::endl;
//
//    using ComputePoint1 = CartesianWritable<Dim1,PositionVector<Dim1>>;
//    ComputePoint1 compute_pt1(100.0);
//    compute_pt1.setX(500.0);
//    std::cout<< compute_pt1 << std::endl;

//    ComputePoint2 a;
//    ViewPoint2 b;
//    a.setX(2.0);
//    a = b;
//    std::cout<< a << std::endl;
//    std::cout<< b << std::endl;

//    ViewPoint2 v1, v2;
//    v1 = v2;

//
//    PositionOperator<PositionVector<Dim3>> a, b;
//    a = b;








//    PositionVector<Dim> pt;

    //test<2> tt = {1,2};


//    GeoTest2D test2D;
//
//    test2D.merge();
//
//    test2D.checkOut();


//    GeoTest3D test3D;
//
//    test3D.merge();
//
//    test3D.checkOut();

//    std::list<int> tt{1, 2, 3};
//
//    auto it_1 = tt.begin();
//    auto it_2 = tt.end();
//
//    ++it_1;
//    ++it_1;
//    ++it_1;
//    std::cout << *it_1 << std::endl;
//    std::cout << *it_2 << std::endl;

//
//    std::cout << it_2 - it_1 << std::endl;



//    auto tt = std::make_shared<Test>();
//    auto gg = tt->sharedThis();
//
//    gg->aa = 3;
//
//    std::cout << tt->aa << std::endl;

//    std::cout << "<Start> Merge Edge " << std::endl;
//    test2D.mergeEdge();
//    std::cout << "< End > Merge Edge " << std::endl;
//    test2D.checkOut();

    return 0;
}