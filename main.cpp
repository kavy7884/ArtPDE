//
// Created by Chingkai Chou on 5/2/18.
//
#include <iostream>
#include "project_uility.hpp"
#include "dimension_utility.hpp"
#include "point.hpp"
#include "vertex.hpp"
#include "cell.hpp"
//#include "geometry_data.hpp"


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

//    using GeoData = art_pde::GeometryData< MeshTypeMethod, art_pde::Dim2D, art_pde::CartesianCoordinate>;
//    GeoData data;
//    std::cout << data.PT.getX() << " " <<  data.PT.getY() << std::endl;
//
//    using PtType = GeoData::GeometryData_Traits::PointType;
//    art_pde::Vertex<PtType> V1(PtType(1.0, 2.0));
//    art_pde::Vertex<PtType> V2(std::make_shared<PtType>(3.0, 4.0));
//
//    std::cout << V1.getPoint().getX() << " " <<  V1.getPoint().getY() << std::endl;
//    std::cout << V2.getPtr_point()->getX() << " " <<  V2.getPtr_point()->getY()<< std::endl;

    using PointType = art_pde::Point<art_pde::Dim2D, art_pde::CartesianCoordinate>;
    using VertexType = art_pde::Vertex<PointType>;
    using PrtVertexType = std::shared_ptr<VertexType>;

    std::vector<PrtVertexType> vec_vertex;
    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0,0)));
    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0.5,0)));
    vec_vertex.push_back(std::make_shared<VertexType>(PointType(1.0,0)));
    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0,0.5)));
    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0.5,0.5)));
    vec_vertex.push_back(std::make_shared<VertexType>(PointType(1,0.5)));
    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0,1)));
    vec_vertex.push_back(std::make_shared<VertexType>(PointType(0.5,1)));
    vec_vertex.push_back(std::make_shared<VertexType>(PointType(1,1)));

    std::cout << "Vertex Num: " << " " <<  vec_vertex.size() << std::endl;

    using CellType = art_pde::Cell<PointType>;
    using PrtCellType = std::shared_ptr<CellType>;
    std::vector<PrtCellType> vec_cell;

    art_pde::CellBuilder<PointType , art_pde::Dim2D> cell_builder;
    // Cell 1
    cell_builder.clearVertex();
    cell_builder.addVertex(vec_vertex[0]);
    cell_builder.addVertex(vec_vertex[1]);
    cell_builder.addVertex(vec_vertex[4]);
    cell_builder.addVertex(vec_vertex[3]);
    vec_cell.push_back(cell_builder.create());

    // Cell 2
    cell_builder.clearVertex();
    cell_builder.addVertex(vec_vertex[1]);
    cell_builder.addVertex(vec_vertex[2]);
    cell_builder.addVertex(vec_vertex[5]);
    cell_builder.addVertex(vec_vertex[4]);
    vec_cell.push_back(cell_builder.create());

    // Cell 3
    cell_builder.clearVertex();
    cell_builder.addVertex(vec_vertex[3]);
    cell_builder.addVertex(vec_vertex[4]);
    cell_builder.addVertex(vec_vertex[7]);
    cell_builder.addVertex(vec_vertex[6]);
    vec_cell.push_back(cell_builder.create());

    // Cell 4
    cell_builder.clearVertex();
    cell_builder.addVertex(vec_vertex[4]);
    cell_builder.addVertex(vec_vertex[5]);
    cell_builder.addVertex(vec_vertex[8]);
    cell_builder.addVertex(vec_vertex[7]);
    vec_cell.push_back(cell_builder.create());

    for (int i = 0; i < vec_cell.size(); ++i) {
        std::cout << "Cell : " << i << ", Type = " << vec_cell[i]->getCellTypeInString() << std::endl;
        std::cout << *vec_cell[i] << std::endl;
    }
    


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