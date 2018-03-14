//
// Created by Chingkai Chou on 3/11/18.
//

#ifndef ARTCFD_SYSTEM_HPP
#define ARTCFD_SYSTEM_HPP

#include <memory>
#include "Approximation.hpp"
#include "TestSpace.hpp"
#include "Integration.hpp"
#include "Eigen/Dense"

template <class Dimension, class NumericalMethodUtility>
class SystemBase{
public:
    SystemBase(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &geoData) : geoData(geoData) {}

    std::shared_ptr<DOF_Base> &getDofBase() {
            return dofBase;
    }

    virtual bool Assembly() = 0;

protected:
    std::shared_ptr<DOF_Base> dofBase;
    std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> geoData;
    std::shared_ptr<ShapeFunction> shapeFunction;
    std::shared_ptr<DeformationMapping<Dimension, NumericalMethodUtility>> mapping;
    std::shared_ptr<TrialSpace<Dimension, NumericalMethodUtility>> trialSpace;
    std::shared_ptr<TestSpace<Dimension, NumericalMethodUtility>> testSpace;
    Eigen::MatrixXd systemStiffness;
    Eigen::VectorXd systemForce;

};

template <class Dimension, class NumericalMethodUtility, class DOF_Type>
class System : public SystemBase<Dimension, NumericalMethodUtility>{
public:
    System(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &geoData,
           const std::shared_ptr<DOF<Dimension, DOF_Type>> &dof) : SystemBase<Dimension, NumericalMethodUtility>(geoData), dof(dof) {
        this->dofBase = dof;
        this->shapeFunction = std::make_shared<ShapeFunction>();
        this->mapping = std::make_shared<DeformationMapping<Dimension, NumericalMethodUtility>>(this->geoData, this->shapeFunction);
        this->trialSpace = std::make_shared<TrialSpace<Dimension, NumericalMethodUtility>>(this->geoData, this->shapeFunction, this->mapping);
        this->dofApproximation = std::make_shared<Approximation<Dimension, NumericalMethodUtility, DOF_Type>>(this->trialSpace, this->dof);
        this->testSpace = std::make_shared<TestSpace<Dimension, NumericalMethodUtility>>(*(this->trialSpace));
        this->problemOperation = std::make_shared<Operation<Dimension, NumericalMethodUtility, DOF_Type>>(this->testSpace, this->dofApproximation);
        this->integration = std::make_shared<Integration<Dimension, NumericalMethodUtility, DOF_Type>>(this->geoData, this->problemOperation);
        this->systemStiffness.resize(dof->size(), dof->size());
        this->systemForce.resize(dof->size());
    }

    virtual bool Assembly() override ;
private:
    std::shared_ptr<DOF<Dimension, DOF_Type>> dof;
    std::shared_ptr<Approximation<Dimension, NumericalMethodUtility, DOF_Type>> dofApproximation;
    std::shared_ptr<Integration<Dimension, NumericalMethodUtility, DOF_Type>> integration;
    std::shared_ptr<Operation<Dimension, NumericalMethodUtility, DOF_Type>> problemOperation;
};

template <class Dimension, class NumericalMethodUtility, class DOF_Type>
bool System<Dimension, NumericalMethodUtility, DOF_Type>::Assembly() {

    std::cout << "Integral Cell Data Init" << std::endl;
    auto intGeo = this->integration->getIntGeoData()->getGeoData();
    auto intPt = this->integration->getIntPtData();
    size_t integralCellNum = intGeo->cElement30->getRowSize();

    std::cout << "Operation Data Init" << std::endl;
    auto operate = this->problemOperation;

    std::cout << "Operate Init" << std::endl;
    operate->Init();

    std::cout << "Cell Loop Start" << std::endl;
    for (size_t cellId = 0; cellId < integralCellNum; ++cellId) {
        std::cout << "Cell Id: " << cellId << std::endl;
        operate->ElementLoopInit(cellId);
        for (size_t q = 0; q < intPt->QuadraturePointN; ++q) {
            operate->CalIntegralPoint(intPt->xi[q], intPt->eta[q], intPt->w[q]);
            //std::cout << "xi: " << intPt->xi[q] << ", eta: " << intPt->eta[q] << ", w: " <<intPt->w[q] << std::endl;
        }
        operate->ElementLoopEnd(cellId);
    }

    std::cout << "Operate End" << std::endl;
    operate->End();

    return true;
}

#endif //ARTCFD_SYSTEM_HPP
