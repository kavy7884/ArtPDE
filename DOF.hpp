//
// Created by Chingkai Chou on 3/6/18.
//

#ifndef ARTCFD_DOF_HPP
#define ARTCFD_DOF_HPP

#include <string>
#include "Eigen/Dense"
#include "DimensionUtility.hpp"

class DOF_Type{};
class ScalarDOF : public DOF_Type{};
class VectorDOF : public DOF_Type{};

class DOF_Base{
public:
    DOF_Base(const std::string &dofName) : dofName(dofName) {}
    size_t size(){ return groupNum*lengthDOF;}
    size_t getGroupNum(){return groupNum;};
    size_t getLengthDOF(){return lengthDOF;};
    std::string& getDOFName(){ return dofName;}

protected:
    std::string dofName;
    size_t groupNum{0}, lengthDOF{0};
    std::shared_ptr<Eigen::VectorXd> dofData{nullptr};
};

template <class Dimension, class DOF_Type>
class DOF{};

template <class Dimension>
class DOF<Dimension, ScalarDOF> : public DOF_Base{
public:
    DOF(const std::string &dofName, size_t length) : DOF_Base(dofName) {
        groupNum = 1;
        lengthDOF = length;
        dofData = std::make_shared<Eigen::VectorXd>();
        dofData->resize(lengthDOF);
    }
};

template <class Dimension>
class DOF<Dimension, VectorDOF> : public DOF_Base{
public:
    DOF(const std::string &dofName, size_t length) : DOF_Base(dofName) {
        groupNum = Dimension::Dim;
        lengthDOF = length;
        dofData = std::make_shared<Eigen::VectorXd>();
        dofData->resize(lengthDOF * groupNum);
    }
};


#endif //ARTCFD_DOF_HPP
