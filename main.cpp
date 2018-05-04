//
// Created by Chingkai Chou on 5/2/18.
//
#include <iostream>
#include "project_uility.hpp"
#include "dimension_utility.hpp"
#include "point.hpp"
#include "geometry_data.hpp"


//
//template <typename A, typename B>
//class test1{
//
//public:
//    struct test1_traits{
//        using type_A = A;
//        using type_B = B;
//    };
//
//    test1() {}
//
//    A a;
//    B b;
//};
//
//template <typename C, typename test1_type>
//class test2{
//public:
//    using test1_typeA = typename test1_type::test1_traits::type_A;
//    using test1_typeB = typename test1_type::test1_traits::type_B;
//
//    struct test2_traits{
//        using test1_trait = test1_type;
//        using type_C = C;
//    };
//
//    test2(const test1<test1_typeA, test1_typeB> &tt1) : tt1(tt1) {}
//
//    C c;
//    test1<test1_typeA, test1_typeB> tt1;
//};
//
//template <typename D, typename test2_type>
//class test3{
//
//public:
//    using test2_typeC = typename test2_type::test2_traits::type_C;
//    using test1_traits = typename test2_type::test2_traits::test1_trait;
//
//    test3(const test2<test2_typeC, test1_traits> &tt2) : tt2(tt2) {}
//
//    D d;
//    test2<test2_typeC, test1_traits> tt2;
//};

int main() {


    auto proj_1 = ArtProject::create("TestProj").setRunPath(".").setDivideSlash("/").setInitialFolderName("Init").build();

//    using T1Type = test1<int, double >;
//    T1Type t1;
//
//    using T2Type = test2<size_t , T1Type::test1_traits >;
//    T2Type t2(t1);
//
//    using T3Type = test3<double , T2Type::test2_traits >;
//    T2Type t3(t2);

    art_pde::GeometryData< MeshTypeMethod, art_pde::Dim2D, art_pde::CartesianCoordinate> data;
    std::cout << data.PT.getX() << " " <<  data.PT.getY() << std::endl;







//    art_pde::Geometry<FEM> geo(2);


//
//    test ttt(P);
//
//
//
//    const size_t i = 2;
//
//    Dimension<i>::DimClass aa;
//    std::cout << aa.kNumDim << std::endl;
//
//    std::array<double, aa.kNumDim> bb;




    return 0;
}