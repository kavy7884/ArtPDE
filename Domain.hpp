//
// Created by Chingkai Chou on 3/5/18.
//

#ifndef ARTCFD_DOMAIN_HPP
#define ARTCFD_DOMAIN_HPP

#include <vector>
#include "System.hpp"
#include "DOF_Mannger.hpp"

template <class Dimension, class NumericalMethodUtility>
class Domain{
public:
    Domain() {}

    bool addSystem(std::shared_ptr<SystemBase<Dimension, NumericalMethodUtility>> system);
    bool Assembly(){
        for (size_t i = 0; i < systemData.size(); ++i) {
            systemData[i]->Assembly();
        }
        return true;
    };

private:
    std::vector<std::shared_ptr<SystemBase<Dimension, NumericalMethodUtility>>> systemData;
    DOF_Mannger dofMannger;

};

template<class Dimension, class NumericalMethodUtility>
bool Domain<Dimension, NumericalMethodUtility>::addSystem(
        std::shared_ptr<SystemBase<Dimension, NumericalMethodUtility>> system) {

    systemData.push_back(system);
    dofMannger.addDofData(system->getDofBase());

    //std::cout << dofMannger.getTotalDof() << std::endl;
    return true;
}

#endif //ARTCFD_DOMAIN_HPP
