//
// Created by Chingkai Chou on 5/4/18.
//

#ifndef ARTCFD_CELL_HPP
#define ARTCFD_CELL_HPP

#include <memory>
#include <vector>
#include <ostream>
#include <list>
#include "vertex.hpp"
#include "edge.hpp"
#include "dimension_utility.hpp"


namespace art_pde {

    template <typename PointType> class Vertex;

    class CellInterface{
    public:
        virtual void genEdgeData() = 0;
    };

    // Geometry cell define
    template <typename PointType>
    class Cell: public CellInterface{
    public:
        // Geometry cell type define
        enum class CellDefineType{ None, Line, Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid, Polyhedron };

    protected:
        using PtrPointType = std::shared_ptr<PointType>;
        using VertexType = art_pde::Vertex<PointType>;
        using PtrVertexType = std::shared_ptr<VertexType>;
        using VecPtrVertexType = std::vector<PtrVertexType>;
        using EdgeType = Edge<PointType>;
        using PtrEdgeType = std::shared_ptr<EdgeType>;
        using ListPtrEdgeType = std::list<PtrEdgeType>;

    public:

        Cell() {}

        const CellDefineType& getCell_define_Type() const;

        void setCell_define_Type(CellDefineType cell_define_Type);

        const VecPtrVertexType &getVec_ptr_vetex() const {
            return vec_ptr_vetex;
        }

        void setVec_ptr_vetex(const VecPtrVertexType &ptr_vetex) {
            Cell::vec_ptr_vetex = ptr_vetex;
        }

        const size_t getNumVertex() const{ return vec_ptr_vetex.size(); };

        std::string getCellTypeInString();

        void calPtr_cell_center_point();

        const PtrPointType &getPtr_cell_center_point() const{
            return ptr_cell_center_point;
        }

        template <typename PointType_>
        friend std::ostream &operator<<(std::ostream &os, const Cell<PointType_> &cell) {
            auto vec_ptr_vetex = cell.getVec_ptr_vetex();
            for (int i = 0; i < cell.getNumVertex(); ++i) {
                os << *vec_ptr_vetex[i] << "\t";
            }
            return os;
        }

        static std::string convertCellTypeInString(CellDefineType cellType_);

        void genEdgeData() override {};

    protected:
        CellDefineType cell_define_Type {CellDefineType::None};
        VecPtrVertexType vec_ptr_vetex;
        PtrPointType ptr_cell_center_point {nullptr};
        ListPtrEdgeType list_ptr_neighbor_edge;

    };

    template<typename PointType>
    std::string Cell<PointType>::getCellTypeInString() {
        return Cell<PointType>::convertCellTypeInString(cell_define_Type);
    }

    template<typename PointType>
    std::string Cell<PointType>::convertCellTypeInString(CellDefineType cellType_) {
        std::string re_string;
        switch (cellType_){
            case CellDefineType::Line : re_string = "Line"; break;
            case CellDefineType::Triangle : re_string = "Triangle"; break;
            case CellDefineType::Quadrilateral : re_string = "Quadrilateral"; break;
            case CellDefineType::Tetrahedron : re_string = "Tetrahedron"; break;
            case CellDefineType::Hexahedron : re_string = "Hexahedron"; break;
            default:
                re_string = "None"; break;
        }
        return re_string;
    }

    template<typename PointType>
    void Cell<PointType>::calPtr_cell_center_point() {
        ptr_cell_center_point = std::make_shared<PointType>();
        for (size_t i = 0; i < vec_ptr_vetex.size(); ++i) {
            *ptr_cell_center_point += vec_ptr_vetex[i]->getPoint();
        }
        *ptr_cell_center_point /= double(vec_ptr_vetex.size());
    }

    template<typename PointType>
    void Cell<PointType>::setCell_define_Type(Cell::CellDefineType cell_define_Type) {
        Cell::cell_define_Type = cell_define_Type;
    }

    template<typename PointType>
    const typename Cell<PointType>::CellDefineType &Cell<PointType>::getCell_define_Type() const {
        return cell_define_Type;
    }


    // Geometry triangleCell define
    template <typename PointType>
    class TriangleCell : public Cell<PointType>{
        using VertexType = art_pde::Vertex<PointType>;
        using PtrVertexType = std::shared_ptr<VertexType>;
        using EdgeType = Edge<PointType>;
        using PtrEdgeType = std::shared_ptr<EdgeType>;
    public:
        TriangleCell(const PtrVertexType &v_1, const PtrVertexType &v_2, const PtrVertexType &v_3)
                : Cell<PointType>()
        {
            this->vec_ptr_vetex.push_back(v_1); this->vec_ptr_vetex.push_back(v_2); this->vec_ptr_vetex.push_back(v_3);
            this->vec_ptr_vetex.shrink_to_fit();
            Cell<PointType>::setCell_define_Type(Cell<PointType>::CellDefineType::Triangle);
        }

        void genEdgeData() override {
            this->list_ptr_neighbor_edge.push_back(genPtrLineEdge(this->vec_ptr_vetex[0], this->vec_ptr_vetex[1]));
            this->list_ptr_neighbor_edge.push_back(genPtrLineEdge(this->vec_ptr_vetex[1], this->vec_ptr_vetex[2]));
            this->list_ptr_neighbor_edge.push_back(genPtrLineEdge(this->vec_ptr_vetex[2], this->vec_ptr_vetex[0]));
        };

    private:
        PtrEdgeType genPtrLineEdge(const PtrVertexType &v_1, const PtrVertexType &v_2){
            PtrEdgeType re_ptr_edge = std::make_shared<LineEdge<PointType>>(v_1, v_2);
            v_1->addPtrNeighborEdge(re_ptr_edge);
            v_2->addPtrNeighborEdge(re_ptr_edge);
            re_ptr_edge->addPtrNeighborCell(std::shared_ptr<TriangleCell<PointType>>(this));
            return re_ptr_edge;
        }
    };

    // Geometry quadrilateralCell define
    template <typename PointType>
    class QuadrilateralCell : public Cell<PointType>{
        using VertexType = art_pde::Vertex<PointType>;
        using PtrVertexType = std::shared_ptr<VertexType>;
        using EdgeType = Edge<PointType>;
        using PtrEdgeType = std::shared_ptr<EdgeType>;
    public:
        QuadrilateralCell(const PtrVertexType &v_1, const PtrVertexType &v_2,
                       const PtrVertexType &v_3, const PtrVertexType &v_4)
                : Cell<PointType>()
        {
            this->vec_ptr_vetex.push_back(v_1); this->vec_ptr_vetex.push_back(v_2);
            this->vec_ptr_vetex.push_back(v_3); this->vec_ptr_vetex.push_back(v_4);
            this->vec_ptr_vetex.shrink_to_fit();
            Cell<PointType>::setCell_define_Type(Cell<PointType>::CellDefineType::Quadrilateral);
        }

        void genEdgeData() override {
            this->list_ptr_neighbor_edge.push_back(genPtrLineEdge(this->vec_ptr_vetex[0], this->vec_ptr_vetex[1]));
            this->list_ptr_neighbor_edge.push_back(genPtrLineEdge(this->vec_ptr_vetex[1], this->vec_ptr_vetex[2]));
            this->list_ptr_neighbor_edge.push_back(genPtrLineEdge(this->vec_ptr_vetex[2], this->vec_ptr_vetex[3]));
            this->list_ptr_neighbor_edge.push_back(genPtrLineEdge(this->vec_ptr_vetex[3], this->vec_ptr_vetex[0]));
        };

    private:
        PtrEdgeType genPtrLineEdge(const PtrVertexType &v_1, const PtrVertexType &v_2){
            PtrEdgeType re_ptr_edge = std::make_shared<LineEdge<PointType>>(v_1, v_2);
            v_1->addPtrNeighborEdge(re_ptr_edge);
            v_2->addPtrNeighborEdge(re_ptr_edge);
            re_ptr_edge->addPtrNeighborCell(std::shared_ptr<QuadrilateralCell<PointType>>(this));
            return re_ptr_edge;
        }
    };

    // Geometry cell builder define (base)
    template <typename PointType, typename Dimension>
    class CellFactoryBase{
        using VertexType = art_pde::Vertex<PointType>;
        using PtrVertexType = std::shared_ptr<VertexType>;
        using VecPtrVertexType = std::vector<PtrVertexType>;
        using PtrCellType = typename std::shared_ptr<Cell<PointType>>;
    public:
        CellFactoryBase() {}

        void addVertex(PtrVertexType & ptr_vertex){ vec_ptr_vertex.push_back(ptr_vertex); }
        void clearVertex(){ vec_ptr_vertex.clear(); }
        const size_t getNumVertex() const { return vec_ptr_vertex.size(); };

        virtual PtrCellType create() = 0;

    protected:
        VecPtrVertexType vec_ptr_vertex;

    };

    // Geometry cell builder abstract define
    template <typename PointType, typename Dimension>
    class CellFactory;

    // Geometry cell builder partial define in 2D case
    template <typename PointType>
    class CellFactory<PointType, Dim2D> :public CellFactoryBase<PointType, Dim2D>{
    public:
        using PtrCellType = typename std::shared_ptr<Cell<PointType>>;
        CellFactory(): CellFactoryBase<PointType, Dim2D>(){}

        PtrCellType create() override {

            PtrCellType ptr_cell{nullptr};

            if(this->getNumVertex() == 3) {
                ptr_cell = std::make_shared<TriangleCell<PointType>>(
                        this->vec_ptr_vertex[0], this->vec_ptr_vertex[1], this->vec_ptr_vertex[2]);
            }
            else if (this->getNumVertex() == 4){
                ptr_cell = std::make_shared<QuadrilateralCell<PointType>>(
                        this->vec_ptr_vertex[0], this->vec_ptr_vertex[1], this->vec_ptr_vertex[2], this->vec_ptr_vertex[3]);
            }
            else{
                ptr_cell = std::make_shared<Cell<PointType>>();
            }

            ptr_cell->genEdgeData();

            return ptr_cell;
        }

    };

    // Geometry cell builder partial define in 3D case
    template <typename PointType>
    class CellFactory<PointType, Dim3D> :public CellFactoryBase<PointType, Dim3D>{
    public:
        using PtrCellType = typename std::shared_ptr<Cell<PointType>>;
        CellFactory(): CellFactoryBase<PointType, Dim3D>(){}

        PtrCellType create() override {

            PtrCellType ptr_cell{nullptr};

            ptr_cell = std::make_shared<Cell<PointType>>();
            // TODO - Add 3D cell implementation.
            return ptr_cell;
        }
    };

}


#endif //ARTCFD_CELL_HPP
