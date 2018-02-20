#include <iostream>
#include "DataArray.hpp"
#include "DataTable.hpp"

//void DemoDataArray();
int main() {
//    DemoDataArray();

    // DataTable Demo
    // Two kind of data source will be created
    std::shared_ptr<StdVectorTensor2D<int>> data1 = std::make_shared<StdVectorTensor2D<int>>();
    (*data1).emplace_back(std::vector<int>{1,2,3});
    (*data1).emplace_back(std::vector<int>{4,5,6});
    (*data1).emplace_back(std::vector<int>{7,8,9});

    StdVectorTensor2D<double> data2{};
    data2.emplace_back(std::vector<double>{1.1,1.2,1.3});
    data2.emplace_back(std::vector<double>{2.1,2.2,2.3,2.4,2.5,2.6,2.7});
    data2.emplace_back(std::vector<double>{3.1,3.2});
    data2.emplace_back(std::vector<double>{4e-1,5e-1,6e-1});

    // Influence data builder construct memory-safe pointer with 2-D table
    auto dataTableOut1 = DataTable<int>::create("Table1").dataFrom(data1).build();
    auto dataTableOut2 = DataTable<double>::create("Table2").dataFrom(data2).build();

    // Show table results
    std::cout<< (*dataTableOut1) <<std::endl;
    std::cout<< (*dataTableOut2) <<std::endl;

    return 0;
}

void DemoDataArray(){
    // DataArray Demo
    // Two kind of data source will be created
    std::shared_ptr<std::vector<double>> data1 = std::make_shared<std::vector<double>>(5,5.5);
    std::vector<int> data2{1,2,3,4,5,6,7,8,9};

    // Influence data builder construct memory-safe pointer with 1-D array
    auto dataArrayOut1 = DataArray<double>::create("Array1").dataFrom(data1).build();
    auto dataArrayOut2 = DataArray<int>::create("Array2").dataFrom(data2).build();

    // Show array results
    std::cout<< (*dataArrayOut1) <<std::endl;
    std::cout<< (*dataArrayOut2) <<std::endl;

    // Copy (Constructor & Assignment)
    DataArray<double> testCopyArray1{};
    testCopyArray1 = (*dataArrayOut1);
    testCopyArray1.setName("testCopyArray1");

    DataArray<int> testCopyArray2(*dataArrayOut2);
    testCopyArray2.setName("testCopyArray2");

    std::cout<< testCopyArray1 <<std::endl;
    std::cout<< testCopyArray2 <<std::endl;

    // Move (Constructor & Assignment)
    decltype(testCopyArray1) testMoveArray1;
    testMoveArray1 = std::move(testCopyArray1);
    testMoveArray1.setName("testCopyArray1(Move)");

    std::cout << "After Move:" << std::endl;
    std::cout<< testMoveArray1 <<std::endl;
    std::cout << "Original Data:" << std::endl;
    std::cout<< testCopyArray1 <<std::endl;

    decltype(testCopyArray2) testMoveArray2(std::move(testCopyArray2));
    testMoveArray2.setName("testCopyArray2(Move)");
    std::cout << "After Move:" << std::endl;
    std::cout<< testMoveArray2 <<std::endl;
    std::cout << "Original Data:" << std::endl;
    std::cout<< testCopyArray2 <<std::endl;

    // Access element
    size_t pos = 3;
    std::cout<<testMoveArray1.at(pos)<<std::endl;
    std::cout<<testMoveArray2(5)<<std::endl;

    // Iterator
    std::cout<< "Iterator Access testMoveArray1" << std::endl;
    for (auto i = testMoveArray1.begin(); i != testMoveArray1.end(); ++i) {
        std::cout<< *i << std::endl;
    }

    std::cout<< "Const Iterator Access testMoveArray2" << std::endl;
    for (auto i = testMoveArray2.cbegin(); i != testMoveArray2.cend(); ++i) {
        std::cout<< *i << std::endl;
    }

    // Error
    try {
        std::cout<<testMoveArray2(100)<<std::endl;
    }
    catch(std::exception &err){
        std::cout << err.what() <<std::endl;
    }
}