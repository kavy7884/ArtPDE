//
// Created by Chingkai Chou on 5/2/18.
//
#include <iostream>
#include <memory>
#include <vector>
#include "Geometry_New/point3.hpp"
#include "Geometry_New/geo_data.hpp"
#include "Geometry_New/geo_data_factory.hpp"


int main() {
    using VecPtrVertex = std::vector<std::shared_ptr<Vertex<Point3>>>;
    using VecPtrCell = std::vector<std::shared_ptr<Cell<Point3>>>;

    VecPtrVertex vec_ptr_vertex;
    VecPtrCell vec_ptr_cell;

    vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.0, 0.0, 0.0));
    vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.5, 0.0, 0.0));
    vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(1.0, 0.0, 0.0));
    vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.0, 0.5, 0.0));
    vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.5, 0.5, 0.0));
    vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(1.0, 0.5, 0.0));
    vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.0, 1.0, 0.0));
    vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(0.5, 1.0, 0.0));
    vec_ptr_vertex.push_back(std::make_shared<Vertex<Point3>>(1.0, 1.0, 0.0));

    vec_ptr_cell.push_back(std::make_shared<QuadCell<Point3>>(
            vec_ptr_vertex[0], vec_ptr_vertex[1], vec_ptr_vertex[4], vec_ptr_vertex[3])
    );
    vec_ptr_cell.push_back(std::make_shared<QuadCell<Point3>>(
            vec_ptr_vertex[1], vec_ptr_vertex[2], vec_ptr_vertex[5], vec_ptr_vertex[4])
    );
    vec_ptr_cell.push_back(std::make_shared<QuadCell<Point3>>(
            vec_ptr_vertex[3], vec_ptr_vertex[4], vec_ptr_vertex[6], vec_ptr_vertex[6])
    );
    vec_ptr_cell.push_back(std::make_shared<QuadCell<Point3>>(
            vec_ptr_vertex[4], vec_ptr_vertex[5], vec_ptr_vertex[8], vec_ptr_vertex[7])
    );

    std::cout << "List Vertex: " << std::endl;
    for (auto &ptr_vertex: vec_ptr_vertex) {
        std::cout << *ptr_vertex << std::endl;
    }


    std::cout << "List Cell's Vertex: " << std::endl;

    for (size_t i = 0;  i < vec_ptr_cell.size(); ++i) {
        std::cout << "Cell " << i << " is : ";
        for (auto &ptr_edge: vec_ptr_cell[i]->getChild()){
            std::cout << *(*ptr_edge->getChild().begin()) ;
        }
        std::cout << std::endl;
    }







    return 0;
}