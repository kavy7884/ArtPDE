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
    bool CalIntegralPoint(size_t &elementId, double &xi, double &eta, double &w,
                          Eigen::MatrixXd &systemStiffness, Eigen::VectorXd &systemForce);
    bool ElementLoopEnd(size_t &elementId){return true;};
    bool End(){ return true; };

protected:
    std::shared_ptr<Approximation<Dimension, NumericalMethodUtility, DOF_Type>> dofApproximation;
    std::shared_ptr<TestSpace<Dimension, NumericalMethodUtility>> testSpace;

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
bool Operation<Dimension, NumericalMethodUtility, DOF_Type>::CalIntegralPoint(size_t &elementId, double &xi, double &eta, double &w, Eigen::MatrixXd &systemStiffness, Eigen::VectorXd &systemForce){

    auto du_dx = dofApproximation->getTrialSpace()->cal_dN_dx(xi, eta);
    auto dv_dx = testSpace->cal_dN_dx(xi, eta);
    auto jacobin = dofApproximation->getTrialSpace()->getMappingFunction()->getJacobin();
    auto localStiffness = (dv_dx.transpose() * du_dx) * jacobin * w;

    auto localDOF = dofApproximation->getTrialSpace()->getGeoData()->getGeoData()->cElement30->row(elementId);

    size_t rowId{0}, colId{0};
    for (size_t i = 0; i < localDOF.size(); ++i) {
        rowId = localDOF.col(i);
        for (size_t j = 0; j < localDOF.size(); ++j){
            colId = localDOF.col(j);
            systemStiffness(rowId, colId) += localStiffness(i, j);
        }
    }

    return true;
};


#endif //ARTCFD_OPERATOR_HPP
