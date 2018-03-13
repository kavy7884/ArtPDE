//
// Created by Chingkai Chou on 3/11/18.
//

#ifndef ARTCFD_SYSTEM_HPP
#define ARTCFD_SYSTEM_HPP

#include <memory>
#include "Approximation.hpp"
#include "TestSpace.hpp"
#include "Integration.hpp"

template <class Dimension, class NumericalMethodUtility>
class SystemBase{
public:
    SystemBase(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &geoData) : geoData(geoData) {}

protected:
    std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> geoData;
    std::shared_ptr<ShapeFunction> shapeFunction;
    std::shared_ptr<DeformationMapping<Dimension, NumericalMethodUtility>> mapping;
    std::shared_ptr<TrialSpace<Dimension, NumericalMethodUtility>> trialSpace;
    std::shared_ptr<TestSpace<Dimension, NumericalMethodUtility>> testSpace;

};

template <class Dimension, class NumericalMethodUtility, class DOF_Type>
class System : public SystemBase<Dimension, NumericalMethodUtility>{
public:
    System(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &geoData,
           const std::shared_ptr<DOF<Dimension, DOF_Type>> &dof) : SystemBase<Dimension, NumericalMethodUtility>(geoData), dof(dof) {
        this->shapeFunction = std::make_shared<ShapeFunction>();
        this->mapping = std::make_shared<DeformationMapping<Dimension, NumericalMethodUtility>>(this->geoData, this->shapeFunction);
        this->trialSpace = std::make_shared<TrialSpace<Dimension, NumericalMethodUtility>>(this->geoData, this->shapeFunction, this->mapping);
        this->dofApproximation = std::make_shared<Approximation<Dimension, NumericalMethodUtility, DOF_Type>>(this->trialSpace, this->dof);
        this->testSpace = std::make_shared<TestSpace<Dimension, NumericalMethodUtility>>(*(this->trialSpace));
        this->problemOperation = std::make_shared<Operation<Dimension, NumericalMethodUtility, DOF_Type>>(this->testSpace, this->dofApproximation);
        this->integration = std::make_shared<Integration<Dimension, NumericalMethodUtility, DOF_Type>>(this->geoData, this->problemOperation);
    }

private:

    std::shared_ptr<DOF<Dimension, DOF_Type>> dof;
    std::shared_ptr<Approximation<Dimension, NumericalMethodUtility, DOF_Type>> dofApproximation;
    std::shared_ptr<Integration<Dimension, NumericalMethodUtility, DOF_Type>> integration;
    std::shared_ptr<Operation<Dimension, NumericalMethodUtility, DOF_Type>> problemOperation;
};

#endif //ARTCFD_SYSTEM_HPP
