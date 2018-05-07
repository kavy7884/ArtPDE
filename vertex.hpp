//
// Created by Chingkai Chou on 5/3/18.
//

#ifndef ARTCFD_VERTEX_HPP
#define ARTCFD_VERTEX_HPP

#include <memory>
#include <ostream>

namespace art_pde {

    // Geometry vertex define
    template <typename PointType>
    class Vertex{
    public:
        using PtrPointType = std::shared_ptr<PointType>;
        Vertex();
        Vertex(const PointType& point);
        Vertex(const PtrPointType& ptr_point);

        const PointType &getPoint() const;
        const PtrPointType &getPtr_point() const;
        void setPtr_point(const PtrPointType &ptr_point);

        template <typename PointType_>
        friend std::ostream &operator<<(std::ostream &os, const Vertex<PointType_> &vertex){
            os << vertex.getPoint();
            return os;
        }

    private:
        PtrPointType ptr_point{nullptr};
    };

    template<typename PointType>
    Vertex<PointType>::Vertex() {
        ptr_point = std::make_shared<PointType>();
    }

    template<typename PointType>
    Vertex<PointType>::Vertex(const PointType &point) {
        ptr_point = std::make_shared<PointType>(point);
    }

    template<typename PointType>
    Vertex<PointType>::Vertex(const Vertex::PtrPointType &ptr_point) {
        Vertex<PointType>::ptr_point = ptr_point;
    }

    template<typename PointType>
    const typename Vertex<PointType>::PtrPointType &Vertex<PointType>::getPtr_point() const {
        return ptr_point;
    }

    template<typename PointType>
    void Vertex<PointType>::setPtr_point(const Vertex<PointType>::PtrPointType &ptr_point) {
        Vertex::ptr_point = ptr_point;
    }

    template<typename PointType>
    const PointType &Vertex<PointType>::getPoint() const {
        return *ptr_point;
    }

}

#endif //ARTCFD_VERTEX_HPP
