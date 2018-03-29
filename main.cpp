#include <iostream>
#include "ProjectUility.hpp"
#include "Geometry.hpp"
//#include "Domain.hpp"
//#include "Solver.hpp"

#include "Point.hpp"
#include "GeometryLoader.hpp"

//void DemoDataArray();
//void DemoDataTable();


int main() {
    auto proj_1 = ArtProject::create("TestProj").setRunPath(".").setDivideSlash("/").setInitialFolderName("Init").build();

    auto quadMeshForIntgration = Geometry<Dim2D, MeshTypeMethod>::create("Test2D_Quad").setSourceFormat(GeometrySourceFormat::File)
                                                   .load(proj_1).build();

//    auto blmMeshForApproximation = Geometry<Dim2D, MeshTypeMethod>::create("Test2D_BLM").load(proj_1).build();

    auto geoSetting = Geometry<Dim2D, FEM>::create("Test2D_BLM").setSourceFormat(GeometrySourceFormat::File);
    auto geoLoading = geoSetting.load(proj_1);
    auto geoPreProcess = geoLoading.preprocess();
    auto blmMeshForApproximation = geoPreProcess.build();


    std::cout << "---- Int Mesh ----" << std::endl;
    std::cout << *(quadMeshForIntgration->xNode) << std::endl;
    std::cout << *(quadMeshForIntgration->cElement30) << std::endl;

    std::cout << "---- Approxi Mesh ----" << std::endl;
    std::cout << *(blmMeshForApproximation->xNode) << std::endl;
    std::cout << *(blmMeshForApproximation->cElement30) << std::endl;


//    auto one2Dpoint = blmMeshForApproximation->xNode->getPoint(1);
//    one2Dpoint->x(); one2Dpoint->y();
//
//    auto eleConnect = blmMeshForApproximation->cElement30->getElementConnectivity(3);
//    eleConnect->getVertex()[1]->x();



//    // ArtCFD project
//    ProjectArt proj("TestProj");

//
//    // Geometry set up
//    auto pGeo = std::make_shared<Geometry<Dim2D,FEM>>();
//
//    // FEM Geometry Read
//    IO_FileReader IO_Reader(proj);
//    pGeo->getGeoData()->readFile(IO_Reader);
//    pGeo->getGeoData()->DataProcessor();
//
//    // DOF create
//    auto dof = std::make_shared<DOF<Dim2D, ScalarDOF>>("T", pGeo->getGeoData()->xNode->getSize());
//
//    // System create
//    auto sys = std::make_shared<System<Dim2D, FEM, ScalarDOF>>(pGeo, dof);
//
//    // Domain create (and add system)
//    auto domain = std::make_shared<Domain<Dim2D, FEM>>();
//    domain->addSystem(std::static_pointer_cast<SystemBase<Dim2D, FEM>>(sys));
//
//    // Problem solve
//    auto solver = std::make_shared<Solver<Dim2D, FEM>>(domain);
//    solver->solve();
//
//    // Output DOF
//    IO_FileWriter IO_writer(proj);
//    dof->writeFile(IO_writer, pGeo->getGeoData()->xNode);


//    DemoDataArray();
//    DemoDataTable();

    return 0;
}

//void DemoDataArray(){
//    // DataArray Demo
//    std::cout<< "======== Start to Demo DataArray ========" <<std::endl;
//    // Two kind of data source will be created
//    std::shared_ptr<std::vector<double>> data1 = std::make_shared<std::vector<double>>(5,5.5);
//    std::vector<int> data2{1,2,3,4,5,6,7,8,9};
//
//    // Influence data builder construct memory-safe pointer with 1-D array
//    auto dataArrayOut1 = DataArray<double>::create("Array1").dataFrom(data1).build();
//    auto dataArrayOut2 = DataArray<int>::create("Array2").dataFrom(data2).build();
//
//    // Show array results
//    std::cout<< (*dataArrayOut1) <<std::endl;
//    std::cout<< (*dataArrayOut2) <<std::endl;
//
//    // Copy (Constructor & Assignment)
//    DataArray<double> testCopyArray1{};
//    testCopyArray1 = (*dataArrayOut1);
//    testCopyArray1.setName("testCopyArray1");
//
//    DataArray<int> testCopyArray2(*dataArrayOut2);
//    testCopyArray2.setName("testCopyArray2");
//
//    std::cout<< testCopyArray1 <<std::endl;
//    std::cout<< testCopyArray2 <<std::endl;
//
//    // Move (Constructor & Assignment)
//    decltype(testCopyArray1) testMoveArray1;
//    testMoveArray1 = std::move(testCopyArray1);
//    testMoveArray1.setName("testCopyArray1(Move)");
//
//    std::cout << "After Move:" << std::endl;
//    std::cout<< testMoveArray1 <<std::endl;
//    std::cout << "Original Data:" << std::endl;
//    std::cout<< testCopyArray1 <<std::endl;
//
//    decltype(testCopyArray2) testMoveArray2(std::move(testCopyArray2));
//    testMoveArray2.setName("testCopyArray2(Move)");
//    std::cout << "After Move:" << std::endl;
//    std::cout<< testMoveArray2 <<std::endl;
//    std::cout << "Original Data:" << std::endl;
//    std::cout<< testCopyArray2 <<std::endl;
//
//    // Access element
//    size_t pos = 3;
//    std::cout<<testMoveArray1.at(pos)<<std::endl;
//    std::cout<<testMoveArray2(5)<<std::endl;
//
//    // Iterator
//    std::cout<< "Iterator Access testMoveArray1" << std::endl;
//    for (auto i = testMoveArray1.begin(); i != testMoveArray1.end(); ++i) {
//        std::cout<< *i << std::endl;
//    }
//
//    std::cout<< "Const Iterator Access testMoveArray2" << std::endl;
//    for (auto i = testMoveArray2.cbegin(); i != testMoveArray2.cend(); ++i) {
//        std::cout<< *i << std::endl;
//    }
//
//    // Error
//    try {
//        std::cout<<testMoveArray2(100)<<std::endl;
//    }
//    catch(std::exception &err){
//        std::cout << err.what() <<std::endl;
//    }
//}
//
//void DemoDataTable(){
//    // DataTable Demo
//    std::cout<< "======== Start to Demo DataTable ========" <<std::endl;
//    // Two kind of data source will be created
//    std::shared_ptr<StdVectorTensor2D<int>> data1 = std::make_shared<StdVectorTensor2D<int>>();
//    (*data1).emplace_back(std::vector<int>{1,2,3});
//    (*data1).emplace_back(std::vector<int>{4,5,6});
//    (*data1).emplace_back(std::vector<int>{7,8,9});
//
//    StdVectorTensor2D<double> data2{};
//    data2.emplace_back(std::vector<double>{1.1,1.2,1.3});
//    data2.emplace_back(std::vector<double>{2.1,2.2,2.3,2.4,2.5,2.6,2.7});
//    data2.emplace_back(std::vector<double>{3.1,3.2});
//    data2.emplace_back(std::vector<double>{4e-1,5e-1,6e-1});
//
//    // Influence data builder construct memory-safe pointer with 2-D table
//    auto dataTableOut1 = DataTable<int>::create("Table1").dataFrom(data1).build();
//    auto dataTableOut2 = DataTable<double>::create("Table2").dataFrom(data2).build();
//
//    // Show table results
//    std::cout<< (*dataTableOut1) <<std::endl;
//    std::cout<< (*dataTableOut2) <<std::endl;
//
//    // Access element
//    std::cout<< "Access dataTableOut1(2,0)" << std::endl;
//    std::cout<< (*dataTableOut1).at(2, 0) <<std::endl;
//    std::cout<< "Access dataTableOut2(3,1)" << std::endl;
//    std::cout<< (*dataTableOut2)(3,1) <<std::endl;
//    std::cout<< "Access dataTableOut2(3,1)" << std::endl;
//    std::cout<< (*dataTableOut2).row(3).col(2) << std::endl;
//
//    // Copy (Constructor & Assignment)
//    DataTable<int> testCopyTable1{};
//    testCopyTable1 = (*dataTableOut1);
//    testCopyTable1.setName("testCopyTable1");
//
//    DataTable<double> testCopyTable2(*dataTableOut2);
//    testCopyTable2.setName("testCopyTable2");
//
//    std::cout<< testCopyTable1 <<std::endl;
//    std::cout<< testCopyTable2 <<std::endl;
//
//    // Move (Constructor & Assignment)
//    decltype(testCopyTable1) testMoveTable1;
//    testMoveTable1 = std::move(testCopyTable1);
//    testMoveTable1.setName("testCopyTable1(Move)");
//
//    std::cout << "After Move:" << std::endl;
//    std::cout<< testMoveTable1 <<std::endl;
//    std::cout << "Original Data:" << std::endl;
//    std::cout<< testCopyTable1 <<std::endl;
//
//    decltype(testCopyTable2) testMoveTable2(std::move(testCopyTable2));
//    testMoveTable2.setName("testCopyTable2(Move)");
//    std::cout << "After Move:" << std::endl;
//    std::cout<< testMoveTable2 <<std::endl;
//    std::cout << "Original Data:" << std::endl;
//    std::cout<< testCopyTable2 <<std::endl;
//
//    // Iterator
//    std::cout<< "Iterator Access row 1 in testMoveTable1" << std::endl;
//    for (auto i = testMoveTable1.row(1).begin(); i != testMoveTable1.row(1).end(); ++i) {
//        std::cout<< *i << std::endl;
//    }
//
//    std::cout<< "Iterator Access row 1 in testMoveTable2" << std::endl;
//    for (auto i = testMoveTable2.row(1).cbegin(); i != testMoveTable2.row(1).cend(); ++i) {
//        std::cout<< *i << std::endl;
//    }
//
//    // Isolated column
//    auto TableOut2_C0 = (*dataTableOut2).row(0);
//    auto TableOut2_C1 = (*dataTableOut2).row(1);
//    auto TableOut2_C2 = (*dataTableOut2).row(2);
//
//    std::cout << TableOut2_C0 << std::endl;
//    std::cout << TableOut2_C1 << std::endl;
//    std::cout << TableOut2_C2 << std::endl;
//
//    // Error
//    try {
//        std::cout<< "Try access testMoveTable1(100,1)"<<std::endl;
//        std::cout<<testMoveTable1(100,1)<<std::endl;
//    }
//    catch(std::exception &err){
//        std::cout << err.what() <<std::endl;
//    }
//
//    try {
//        std::cout<< "Try access testMoveTable2(1,100)"<<std::endl;
//        std::cout<<testMoveTable2(1,100)<<std::endl;
//    }
//    catch(std::exception &err){
//        std::cout << err.what() <<std::endl;
//    }
//
//    try {
//        std::cout<< "Try access testMoveTable2(100,100)"<<std::endl;
//        std::cout<<testMoveTable2(100,100)<<std::endl;
//    }
//    catch(std::exception &err){
//        std::cout << err.what() <<std::endl;
//    }
//
//}