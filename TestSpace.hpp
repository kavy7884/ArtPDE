//
// Created by Chingkai Chou on 3/7/18.
//

#ifndef ARTCFD_TESTBASIS_HPP
#define ARTCFD_TESTBASIS_HPP

#include "TrialSpace.hpp"

template <class Dimension, class NumericalMethodUtility>
class TestSpace: public TrialSpace<Dimension, NumericalMethodUtility>{
public:
    TestSpace(TrialSpace<Dimension, NumericalMethodUtility> &trialSpace)
    :TrialSpace<Dimension, NumericalMethodUtility>(trialSpace.getGeoData(), trialSpace.getShapeFunction(), trialSpace.getMappingFunction()){}


    TestSpace(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &GeoData,
                   const std::shared_ptr<ShapeFunction> &shapeFunction,
                   const std::shared_ptr<DeformationMapping<Dimension, NumericalMethodUtility>> &mappingFunction)
            : TrialSpace<Dimension, NumericalMethodUtility>(GeoData, shapeFunction, mappingFunction) {}

    Eigen::MatrixXd cal_W(double &xi, double &eta){
        this->cal_N(xi, eta);
    }
    Eigen::MatrixXd cal_dW_dx(double &xi, double &eta){
        this->cal_dN_dx(xi, eta);
    }
};

#endif //ARTCFD_TESTBASIS_HPP
