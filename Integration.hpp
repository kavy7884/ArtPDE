//
// Created by Chingkai Chou on 3/5/18.
//

#ifndef ARTCFD_INTEGRATE_HPP
#define ARTCFD_INTEGRATE_HPP

#include <memory>
#include "TestSpace.hpp"
#include "Operation.hpp"

template<class Dimension>
class GaussianQuadrature{
public:
    GaussianQuadrature() {
        xi = new double[QuadraturePointN]();
        eta = new double[QuadraturePointN]();
        w = new double[QuadraturePointN]();

        xi[0] = 0.7745966692414834;     eta[0] = 0.0000000000000000;    w[0] = -0.7745966692414834;
        xi[1] = 0.7745966692414834;     eta[1] = 0.0000000000000000;    w[1] = -0.7745966692414834;
        xi[2] = 0.7745966692414834;     eta[2] = 0.0000000000000000;    w[2] = -0.7745966692414834;
        xi[3] = 0.7745966692414834;     eta[3] = 0.7745966692414834;    w[3] = 0.7745966692414834;
        xi[4] = 0.0000000000000000;     eta[4] = 0.0000000000000000;    w[4] = 0.0000000000000000;
        xi[5] = -0.7745966692414834;    eta[5] = -0.7745966692414834;   w[5] = -0.7745966692414834;
        xi[6] = 0.3086419753086406;     eta[6] = 0.4938271604938261;    w[6] = 0.3086419753086406;
        xi[7] = 0.4938271604938261;     eta[7] = 0.7901234567901234;    w[7] = 0.4938271604938261;
        xi[8] = 0.3086419753086406;     eta[8] = 0.4938271604938261;    w[8] = 0.3086419753086406;
    }

    virtual ~GaussianQuadrature() {
        delete[](xi);
        delete[](eta);
        delete[](w);
    }

    size_t QuadraturePointN{9};
    double *xi, *eta, *w;
};

template <class Dimension, class NumericalMethodUtility, class DOF_Type>
class Integration{
public:
    Integration(const std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> &intGeoData,
                const std::shared_ptr<Operation<Dimension, NumericalMethodUtility, DOF_Type>> &problemOperation)
            : intGeoData(intGeoData), problemOperation(problemOperation) {
        intPtData = std::make_shared<GaussianQuadrature<Dimension>>();
    }

private:
    std::shared_ptr<Geometry<Dimension, NumericalMethodUtility>> intGeoData;
    std::shared_ptr<GaussianQuadrature<Dimension>> intPtData;
    std::shared_ptr<Operation<Dimension, NumericalMethodUtility, DOF_Type>> problemOperation;
};


#endif //ARTCFD_INTEGRATE_HPP
