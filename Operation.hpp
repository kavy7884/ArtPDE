//
// Created by Chingkai Chou on 3/5/18.
//

#ifndef ARTCFD_OPERATOR_HPP
#define ARTCFD_OPERATOR_HPP

#include "Approximation.hpp"
#include "TestSpace.hpp"
#include "Eigen/Dense"

template <class Dimension, class NumericalMethodUtility, class DOF_Type>
class Operation{
public:
    Operation(const std::shared_ptr<TestSpace<Dimension, NumericalMethodUtility>> &testSpace,
              const std::shared_ptr<Approximation<Dimension, NumericalMethodUtility, DOF_Type>> &dofApproximation)
            : testSpace(testSpace), dofApproximation(dofApproximation) {}

    bool Init(){return true;};
    bool ElementLoopInit(size_t &elementId);
    bool CalIntegralPoint(double &xi, double &eta, double &w);
    bool ElementLoopEnd(size_t &elementId){return true;};
    bool End(){return true;};

    const Eigen::MatrixXd &getStiffness() const {
        return stiffness;
    }

    const Eigen::MatrixXd &getForce() const {
        return force;
    }

protected:
    std::shared_ptr<Approximation<Dimension, NumericalMethodUtility, DOF_Type>> dofApproximation;
    std::shared_ptr<TestSpace<Dimension, NumericalMethodUtility>> testSpace;

    Eigen::MatrixXd stiffness, force;

};


//template <class Dimension, class NumericalMethodUtility, class DOF_Type>
//class Operation_Laplace: public Operation<Dimension, NumericalMethodUtility, DOF_Type>{};


template <class Dimension, class NumericalMethodUtility, class DOF_Type>
bool Operation<Dimension, NumericalMethodUtility, DOF_Type>::ElementLoopInit(size_t &elementId){
    dofApproximation->getTrialSpace()->ElementInit(elementId);
    //std::cout << dofApproximation->getTrialSpace()->getMappingFunction()->getLocalNode() << std::endl;
    return true;
};

template <class Dimension, class NumericalMethodUtility, class DOF_Type>
bool Operation<Dimension, NumericalMethodUtility, DOF_Type>::CalIntegralPoint(double &xi, double &eta, double &w){

    auto du_dx = dofApproximation->getTrialSpace()->cal_dN_dx(xi, eta);
    auto dv_dx = dofApproximation->getTrialSpace()->cal_dN_dx(xi, eta);

    std::cout << du_dx << std::endl;
    std::cout << dv_dx << std::endl;
    return true;
};


#endif //ARTCFD_OPERATOR_HPP
