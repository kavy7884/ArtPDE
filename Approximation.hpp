//
// Created by Chingkai Chou on 3/7/18.
//

#ifndef ARTCFD_INTERPOLATION_HPP
#define ARTCFD_INTERPOLATION_HPP

#include <memory>
#include "TrialSpace.hpp"
#include "Dof_New.hpp"

template <class Dimension, class NumericalMethodUtility, class DOF_Type>
class Approximation{
public:
    Approximation(const std::shared_ptr<TrialSpace<Dimension, NumericalMethodUtility>> &trialSpace,
                  const std::shared_ptr<DOF<Dimension, DOF_Type>> &dof) : trialSpace(trialSpace), dof(dof) {}

    const std::shared_ptr<TrialSpace<Dimension, NumericalMethodUtility>> &getTrialSpace() const {
        return trialSpace;
    }

    const std::shared_ptr<DOF<Dimension, DOF_Type>> &getDof() const {
        return dof;
    }


private:
    std::shared_ptr<TrialSpace<Dimension, NumericalMethodUtility>> trialSpace;
    std::shared_ptr<DOF<Dimension, DOF_Type>> dof;

};

template <class Dimension, class NumericalMethodUtility, class DOF_Type>
class Interpolation : public Approximation<Dimension, NumericalMethodUtility, DOF_Type >{
public:
    Interpolation(const std::shared_ptr<TrialSpace<Dimension, NumericalMethodUtility>> &trialSpace,
                  const std::shared_ptr<DOF<Dimension, DOF_Type>> &dof) : Approximation<Dimension, NumericalMethodUtility, DOF_Type >(trialSpace, dof) {}

};



#endif //ARTCFD_INTERPOLATION_HPP
