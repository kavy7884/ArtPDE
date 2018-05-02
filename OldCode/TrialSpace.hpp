//
// Created by Chingkai Chou on 3/7/18.
//

#ifndef ARTCFD_TRIALBASIS_HPP
#define ARTCFD_TRIALBASIS_HPP

#include <memory>
#include "Geometry.hpp"
#include "ShapeFunction.hpp"
#include "DeformationMapping.hpp"
#include "Eigen/Dense"

template <class Dimension, class NumericalMethodUtility>
class TrialSpaceBase{
public:
    TrialSpaceBase(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &GeoData,
                   const std::shared_ptr<ShapeFunction> &shapeFunction,
                   const std::shared_ptr<DeformationMapping<Dimension, NumericalMethodUtility>> &mappingFunction)
            : GeoData(GeoData), shapeFunction(shapeFunction), mappingFunction(mappingFunction) {}

    TrialSpaceBase(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &GeoData,
                   const std::shared_ptr<DeformationMapping<Dimension, NumericalMethodUtility>> &mappingFunction)
            : GeoData(GeoData), mappingFunction(mappingFunction) {}

protected:
    std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> GeoData;
    std::shared_ptr<ShapeFunction> shapeFunction;
    std::shared_ptr<DeformationMapping<Dimension, NumericalMethodUtility>> mappingFunction{nullptr};
};


template <class Dimension, class NumericalMethodUtility>
class TrialSpace :public TrialSpaceBase<Dimension, NumericalMethodUtility>{
public:

    TrialSpace(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &GeoData,
               const std::shared_ptr<ShapeFunction> &shapeFunction,
               const std::shared_ptr<DeformationMapping<Dimension, NumericalMethodUtility>> &mappingFunction);

    bool ElementInit(size_t &elementId);

    Eigen::MatrixXd cal_N(double &xi, double &eta);
    Eigen::MatrixXd cal_dN_dx(double &xi, double &eta);

    std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> getGeoData(){ return this->GeoData;};
    std::shared_ptr<ShapeFunction> getShapeFunction(){ return this->shapeFunction;};
    std::shared_ptr<DeformationMapping<Dimension, NumericalMethodUtility>> getMappingFunction(){ return this->mappingFunction;};

};

template<class Dimension, class NumericalMethodUtility>
bool TrialSpace<Dimension, NumericalMethodUtility>::ElementInit(size_t &elementId) {
    this->mappingFunction->setLocalNode(elementId);
    return true;
}

template<class Dimension, class NumericalMethodUtility>
Eigen::MatrixXd TrialSpace<Dimension, NumericalMethodUtility>::cal_N(double &xi, double &eta) {
    return this->shapeFunction->ShapeFunctionQ4_N(xi, eta);
}

template<class Dimension, class NumericalMethodUtility>
Eigen::MatrixXd TrialSpace<Dimension, NumericalMethodUtility>::cal_dN_dx(double &xi, double &eta) {
    auto naturalDerivatives = this->shapeFunction->ShapeFunctionQ4_dN_dxi(xi, eta);
    this->mappingFunction->cal_dx_dxi(xi, eta);
    this->mappingFunction->cal_dxi_dx(xi, eta);
    this->mappingFunction->cal_jacobin();

    return (this->mappingFunction->get_dxi_dx() * naturalDerivatives);
//    return naturalDerivatives;
}

template<class Dimension, class NumericalMethodUtility>
TrialSpace<Dimension, NumericalMethodUtility>::TrialSpace(
        const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &GeoData,
        const std::shared_ptr<ShapeFunction> &shapeFunction,
        const std::shared_ptr<DeformationMapping<Dimension, NumericalMethodUtility>> &mappingFunction)
        :TrialSpaceBase<Dimension, NumericalMethodUtility>(
        GeoData, shapeFunction, mappingFunction) {}

#endif //ARTCFD_TRIALBASIS_HPP
