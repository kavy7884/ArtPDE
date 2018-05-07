//
// Created by Chingkai Chou on 5/4/18.
//

#ifndef ARTCFD_CELL_HPP
#define ARTCFD_CELL_HPP

#include <memory>
#include <vector>
#include <ostream>
#include "vertex.hpp"
#include "dimension_utility.hpp"

namespace art_pde {

    // Geometry cell type define
    enum class CellType{ None, Line, Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid, Polyhedron };

    // Geometry cell define
    template <typename PointType>
    class Cell{
    public:
        using VertexType = art_pde::Vertex<PointType>;
        using PtrVertexType = std::shared_ptr<VertexType>;
        using VecPtrVertexType = std::vector<PtrVertexType>;
        Cell() {}

        CellType getCellType() const {
            return cellType;
        }

        void setCellType(CellType cellType) {
            Cell::cellType = cellType;
        }

        const VecPtrVertexType &getVec_ptr_vetex() const {
            return vec_ptr_vetex;
        }

        void setVec_ptr_vetex(const VecPtrVertexType &ptr_vetex) {
            Cell::vec_ptr_vetex = ptr_vetex;
        }

        const size_t getNumVertex() const{ return vec_ptr_vetex.size(); };

        std::string getCellTypeInString();

        template <typename PointType_>
        friend std::ostream &operator<<(std::ostream &os, const Cell<PointType_> &cell) {
            auto vec_ptr_vetex = cell.getVec_ptr_vetex();
            for (int i = 0; i < cell.getNumVertex(); ++i) {
                os << *vec_ptr_vetex[i] << "\t";
            }
            return os;
        }

    protected:
        CellType cellType{CellType::None};
        VecPtrVertexType vec_ptr_vetex;
    };

    template<typename PointType>
    std::string Cell<PointType>::getCellTypeInString() {
        std::string re_string;
        switch (cellType){
            case CellType::Line : re_string = "Line"; break;
            case CellType::Triangle : re_string = "Triangle"; break;
            case CellType::Quadrilateral : re_string = "Quadrilateral"; break;
            case CellType::Tetrahedron : re_string = "Tetrahedron"; break;
            case CellType::Hexahedron : re_string = "Hexahedron"; break;
            default:
                re_string = "None"; break;
        }
        return re_string;
    }

    // Geometry triangleCell define
    template <typename PointType>
    class TriangleCell : public Cell<PointType>{
    public:
        using PtrVertexType = typename Cell<PointType>::PtrVertexType;
        TriangleCell(const PtrVertexType &v_1, const PtrVertexType &v_2, const PtrVertexType &v_3): Cell<PointType>() {
            this->vec_ptr_vetex.push_back(v_1); this->vec_ptr_vetex.push_back(v_2); this->vec_ptr_vetex.push_back(v_3);
            this->vec_ptr_vetex.shrink_to_fit();
            Cell<PointType>::setCellType(CellType::Triangle);
        }
    };

    // Geometry quadrilateralCell define
    template <typename PointType>
    class QuadrilateralCell : public Cell<PointType>{
    public:
        using PtrVertexType = typename Cell<PointType>::PtrVertexType;
        QuadrilateralCell(const PtrVertexType &v_1, const PtrVertexType &v_2,
                       const PtrVertexType &v_3, const PtrVertexType &v_4): Cell<PointType>() {
            this->vec_ptr_vetex.push_back(v_1); this->vec_ptr_vetex.push_back(v_2);
            this->vec_ptr_vetex.push_back(v_3); this->vec_ptr_vetex.push_back(v_4);
            this->vec_ptr_vetex.shrink_to_fit();
            Cell<PointType>::setCellType(CellType::Quadrilateral);
        }
    };

    // Geometry cell builder define (base)
    template <typename PointType, typename Dimension>
    class CellBuilderBase{

    public:
        using VertexType = Vertex<PointType>;
        using PtrVertexType = typename Cell<PointType>::PtrVertexType;
        using VecPtrVertexType = typename Cell<PointType>::VecPtrVertexType;
        using PtrCellType = typename std::shared_ptr<Cell<PointType>>;

        CellBuilderBase() {}

        void addVertex(PtrVertexType & ptr_vertex){ vec_ptr_vertex.push_back(ptr_vertex); }
        void clearVertex(){ vec_ptr_vertex.clear(); }
        const size_t getNumVertex() const { return vec_ptr_vertex.size(); };

        virtual PtrCellType create() = 0;

    protected:
        VecPtrVertexType vec_ptr_vertex;

    };

    // Geometry cell builder abstract define
    template <typename PointType, typename Dimension>
    class CellBuilder;

    // Geometry cell builder partial define in 2D case
    template <typename PointType>
    class CellBuilder<PointType, Dim2D> :public CellBuilderBase<PointType, Dim2D>{
    public:
        using PtrCellType = typename std::shared_ptr<Cell<PointType>>;
        CellBuilder(): CellBuilderBase<PointType, Dim2D>(){}

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

            return ptr_cell;
        }

    };

    // Geometry cell builder partial define in 3D case
    template <typename PointType>
    class CellBuilder<PointType, Dim3D> :public CellBuilderBase<PointType, Dim3D>{
    public:
        using PtrCellType = typename std::shared_ptr<Cell<PointType>>;
        CellBuilder(): CellBuilderBase<PointType, Dim3D>(){}

        PtrCellType create() override {

            PtrCellType ptr_cell{nullptr};

            ptr_cell = std::make_shared<Cell<PointType>>();
            // TODO - Add 3D cell implementation.
            return ptr_cell;
        }
    };

}


#endif //ARTCFD_CELL_HPP
