//
// Created by Chingkai Chou on 6/7/18.
//

#ifndef ARTPDE_KAVY_DEMOFUNCTIONSPACE_HPP
#define ARTPDE_KAVY_DEMOFUNCTIONSPACE_HPP

#include "Project/art_project.hpp"
#include "Geometry/geometry.hpp"
#include "FunctionSpace_Refactor/function_space.hpp"

void DemoFunctionSpace(){
    using namespace art_pde;
    using namespace PositionVector;
    using namespace project;
    using namespace geometry::mesh_type::Dim2;
    using namespace function_space::isoparametric;

    auto proj = ArtProject::create("TestProj2D").setRunPath(".").build();

    // Step 2: Defining mesh data and which Algorithm and Reader!

    using MeshDataType = GeometricData<ComputePositionVector<2>>;
    using GeoType = Geometry<GeometricAlgorithm<MeshDataType>, GeometricReader<MeshDataType>>;
    auto geo = GeoType::create();

    // Step 3: Reading data from ArtProject's geometry folder.
    //- two files: position.txt and connectivity.txt must existing and defining properly.
    //- if the loading process failure, the program will exit and show something error in this line.
    //- P.S.: the reading method provided by "GeometricReader" API class.

    assert(geo->read(proj)); // Check loading status.

    // Step 4: Merge geometry
    //- 2D case will merge Edge.
    geo->merge();


    using BasisFuncType = BasisFunction<GeoType>;
    using FuncSpaceType = FunctionSpace<BasisFuncType>;

    auto func_space_1 = FuncSpaceType::create(geo);
    auto func_space_2 = FuncSpaceType::create(geo, BasisOrder::Quadratic);

    func_space_1->testBasisFunc();
    func_space_2->testBasisFunc();

}

#endif //ARTPDE_KAVY_DEMOFUNCTIONSPACE_HPP
