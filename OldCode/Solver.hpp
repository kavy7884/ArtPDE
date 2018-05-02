//
// Created by Chingkai Chou on 3/18/18.
//

#ifndef ARTCFD_SOLVER_HPP
#define ARTCFD_SOLVER_HPP

#include "Domain.hpp"

template <class Dimension, class NumericalMethodUtility>
class Solver{
public:
    Solver(const std::shared_ptr<Domain<Dimension, NumericalMethodUtility>> &pDomain) : pDomain(pDomain) {}
    bool solve();

private:
    std::shared_ptr<Domain<Dimension, NumericalMethodUtility>> pDomain{nullptr};
};

template<class Dimension, class NumericalMethodUtility>
bool Solver<Dimension, NumericalMethodUtility>::solve() {

    pDomain->Assembly();

    auto LHS = pDomain->getLHS();
    auto RHS = pDomain->getRHS();

    Eigen::VectorXd sol = LHS.fullPivLu().solve(RHS);


    pDomain->DofSolutionWriteBack(sol);

    //std::cout << RHS << std::endl;

    return true;
}

#endif //ARTCFD_SOLVER_HPP
