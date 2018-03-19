//
// Created by Chingkai Chou on 3/5/18.
//

#ifndef ARTCFD_DOMAIN_HPP
#define ARTCFD_DOMAIN_HPP

#include <vector>
#include "System.hpp"
#include "DOF_Mannger.hpp"
#include "Eigen/Dense"

template <class Dimension, class NumericalMethodUtility>
class Domain{
public:
    Domain() {}

    bool addSystem(std::shared_ptr<SystemBase<Dimension, NumericalMethodUtility>> system);
    bool Assembly(){
        LHS.resize(dofMannger.getTotalDof(), dofMannger.getTotalDof());
        RHS.resize(dofMannger.getTotalDof());
        for (size_t i = 0; i < systemData.size(); ++i) {
            systemData[i]->Assembly();
        }

        // TODO: Add matrix by system

        LHS = systemData[0]->getSystemStiffness();
        RHS = systemData[0]->getSystemForce();

        return true;
    }

    Eigen::MatrixXd &getLHS();

    Eigen::VectorXd &getRHS();

    void DofSolutionWriteBack(Eigen::VectorXd &sol);

private:
    std::vector<std::shared_ptr<SystemBase<Dimension, NumericalMethodUtility>>> systemData;
    DOF_Mannger dofMannger;
    Eigen::MatrixXd LHS;
    Eigen::VectorXd RHS;

};

template<class Dimension, class NumericalMethodUtility>
bool Domain<Dimension, NumericalMethodUtility>::addSystem(
        std::shared_ptr<SystemBase<Dimension, NumericalMethodUtility>> system) {

    systemData.push_back(system);
    dofMannger.addDofData(system->getDofBase());

    //std::cout << dofMannger.getTotalDof() << std::endl;
    return true;
}

template<class Dimension, class NumericalMethodUtility>
Eigen::MatrixXd &Domain<Dimension, NumericalMethodUtility>::getLHS() {
    return LHS;
}

template<class Dimension, class NumericalMethodUtility>
Eigen::VectorXd &Domain<Dimension, NumericalMethodUtility>::getRHS() {
    return RHS;
}

template<class Dimension, class NumericalMethodUtility>
void Domain<Dimension, NumericalMethodUtility>::DofSolutionWriteBack(Eigen::VectorXd &sol) {
    size_t systemNum = dofMannger.getSystemId().size() - 1;
    for (size_t i = 0; i < systemNum; ++i) {
        dofMannger.getDofGroupData()[i]->DofDataAssignment(sol, dofMannger.getSystemId()[i], dofMannger.getSystemId()[i + 1]);
    }
}


#endif //ARTCFD_DOMAIN_HPP
