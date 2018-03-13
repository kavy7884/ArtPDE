//
// Created by Chingkai Chou on 3/9/18.
//

#ifndef ARTCFD_DEFORMATIONMAPPING_HPP
#define ARTCFD_DEFORMATIONMAPPING_HPP

#include "Geometry.hpp"
#include "ShapeFunction.hpp"

template <class Dimension, class NumericalMethodUtility>
class DeformationMappingBase{
public:
    DeformationMappingBase(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &GeoData,
                           const std::shared_ptr<ShapeFunction> &shapeFunction) : GeoData(GeoData),
                                                                                  shapeFunction(shapeFunction) {}

protected:
    std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> GeoData;
    std::shared_ptr<ShapeFunction> shapeFunction;
};

template <class Dimension, class NumericalMethodUtility>
class DeformationMapping : public DeformationMappingBase<Dimension, NumericalMethodUtility>{
public:

    DeformationMapping(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &GeoData,
                       const std::shared_ptr<ShapeFunction> &shapeFunction);

    void setLocalNode(size_t &elementId);
    void cal_dx_dxi(double &xi, double &eta);
    void cal_dxi_dx(double &xi, double &eta);
    void cal_jacobin();

    const Eigen::MatrixXd &getLocalNode() const;

    const Eigen::MatrixXd &get_dx_dxi() const;

    const Eigen::MatrixXd &get_dxi_dx() const;

    double getJacobin() const;



private:
    Eigen::MatrixXd localNode, dx_dxi, dxi_dx;
    double jacobin{0.0};

};

template<class Dimension, class NumericalMethodUtility>
void DeformationMapping<Dimension, NumericalMethodUtility>::cal_dx_dxi(double &xi, double &eta) {
    auto naturalDerivatives = this->shapeFunction->ShapeFunctionQ4_dN_dxi(xi, eta);
    dx_dxi = naturalDerivatives * localNode;
}

template<class Dimension, class NumericalMethodUtility>
void DeformationMapping<Dimension, NumericalMethodUtility>::cal_dxi_dx(double &xi, double &eta) {
    dxi_dx = dx_dxi.inverse();
}

template<class Dimension, class NumericalMethodUtility>
void DeformationMapping<Dimension, NumericalMethodUtility>::cal_jacobin() {
    jacobin = dxi_dx.determinant();
}

template<class Dimension, class NumericalMethodUtility>
const Eigen::MatrixXd &DeformationMapping<Dimension, NumericalMethodUtility>::getLocalNode() const {
    return localNode;
}

template<class Dimension, class NumericalMethodUtility>
const Eigen::MatrixXd &DeformationMapping<Dimension, NumericalMethodUtility>::get_dx_dxi() const {
    return dx_dxi;
}

template<class Dimension, class NumericalMethodUtility>
const Eigen::MatrixXd &DeformationMapping<Dimension, NumericalMethodUtility>::get_dxi_dx() const {
    return dxi_dx;
}

template<class Dimension, class NumericalMethodUtility>
double DeformationMapping<Dimension, NumericalMethodUtility>::getJacobin() const {
    return jacobin;
}

template<class Dimension, class NumericalMethodUtility>
void DeformationMapping<Dimension, NumericalMethodUtility>::setLocalNode(size_t &elementId) {
    auto eleId = this->GeoData->getGeoData()->cElement30->row(elementId);
    auto N = eleId.size();
    localNode.resize(N, 2);
    for (int i = 0; i < N; ++i) {
        localNode(i, 0) = this->GeoData->getGeoData()->xNode->x(i);
        localNode(i, 1) = this->GeoData->getGeoData()->xNode->y(i);
    }
}

template<class Dimension, class NumericalMethodUtility>
DeformationMapping<Dimension, NumericalMethodUtility>::DeformationMapping(
        const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &GeoData,
        const std::shared_ptr<ShapeFunction> &shapeFunction):DeformationMappingBase<Dimension, NumericalMethodUtility>(GeoData, shapeFunction) {}

#endif //ARTCFD_DEFORMATIONMAPPING_HPP
