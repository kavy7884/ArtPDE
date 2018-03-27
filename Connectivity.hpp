//
// Created by Chingkai Chou on 3/26/18.
//

#ifndef ARTCFD_CONNECTIVITY_HPP
#define ARTCFD_CONNECTIVITY_HPP

#include <vector>
#include <memory>
#include <ostream>
#include "Point.hpp"

template <class Dimension>
class Connectivity{
public:
    Connectivity() {}

    size_t getId() const {
        return id;
    }

    void setId(size_t id) {
        Connectivity::id = id;
    }

    friend std::ostream &operator<<(std::ostream &os, const Connectivity &connectivity){
        for (auto &pt : connectivity.points) {
            os << pt->getId() << "\t";
        }
        return os;
    }

protected:
    size_t getPointsSize() const { return points.size();}

    void addPoint(std::shared_ptr<Point<Dimension>> &pointData);

    size_t id{0};
    std::vector<std::shared_ptr<Point<Dimension>>> points;
};

template<class Dimension>
void Connectivity<Dimension>::addPoint(std::shared_ptr<Point<Dimension>> &pointData) {
    return points.push_back(pointData);
}

template <class Dimension>
class ElementConnectivity :public Connectivity<Dimension>{
public:
    ElementConnectivity() : Connectivity<Dimension>(){}

    size_t getVertexSize() const { return this->getPointsSize();}

    void addVertex(std::shared_ptr<Point<Dimension>> &vertexData){ return this->addPoint(vertexData);}

    friend std::ostream &operator<<(std::ostream &os, const ElementConnectivity &connectivity) {
        os << static_cast<const Connectivity<Dimension> &>(connectivity);
        return os;
    }
};

template <class Dimension>
class ElementList{
public:
    using PtrElementConnectivityType = std::shared_ptr<ElementConnectivity<Dimension>>;

    ElementList() {}

    size_t size() const { return data.size();}
    void addElementConnectivity(PtrElementConnectivityType & p_elementConnectivity){
        p_elementConnectivity->setId(size());
        data.push_back(p_elementConnectivity);
    }
    PtrElementConnectivityType & getElementConnectivity(size_t& id){ return data[id];}
    const PtrElementConnectivityType & c_getElementConnectivity(size_t& id) const { return data[id];}

    friend std::ostream &operator<<(std::ostream &os, const ElementList &list) {
        for (auto &connectivity : list.data) {
            os << connectivity->getId() << "\t" << *connectivity << std::endl;
        }
        return os;
    }

private:
    std::vector<PtrElementConnectivityType> data;
};

#endif //ARTCFD_CONNECTIVITY_HPP
