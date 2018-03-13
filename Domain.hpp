//
// Created by Chingkai Chou on 3/5/18.
//

#ifndef ARTCFD_DOMAIN_HPP
#define ARTCFD_DOMAIN_HPP

#include <vector>
#include "System.hpp"

template <class Dimension, class NumericalMethodUtility>
class Domain{
public:
    Domain() {}

    bool addSystem(std::shared_ptr<SystemBase<Dimension, NumericalMethodUtility>> system);

private:
    std::vector<std::shared_ptr<SystemBase<Dimension, NumericalMethodUtility>>> systemData;


};

template<class Dimension, class NumericalMethodUtility>
bool Domain<Dimension, NumericalMethodUtility>::addSystem(
        std::shared_ptr<SystemBase<Dimension, NumericalMethodUtility>> system) {

    systemData.push_back(system);
    return true;
}

#endif //ARTCFD_DOMAIN_HPP
