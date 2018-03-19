//
// Created by Chingkai Chou on 3/6/18.
//

#ifndef ARTCFD_DOF_HPP
#define ARTCFD_DOF_HPP

#include <string>
#include <ostream>
#include "Eigen/Dense"
#include "DimensionUtility.hpp"
#include "IO_Uility.hpp"

class DOF_Type{};
class ScalarDOF : public DOF_Type{};
class VectorDOF : public DOF_Type{};

class DOF_Base{
public:
    DOF_Base(const std::string &dofName) : dofName(dofName) {}
    size_t size() const{ return groupNum*lengthDOF;}
    size_t getGroupNum() const {return groupNum;};
    size_t getLengthDOF() const {return lengthDOF;};
    const std::string& getDOFName() const { return dofName;}

    void DofDataAssignment(Eigen::VectorXd &sol, size_t beginId, size_t endId);

    const std::shared_ptr<Eigen::VectorXd> &getDofData() const;

protected:
    std::string dofName;
    size_t groupNum{0}, lengthDOF{0};
    std::shared_ptr<Eigen::VectorXd> dofData{nullptr};
};

void DOF_Base::DofDataAssignment(Eigen::VectorXd &sol, size_t beginId, size_t endId) {
    if (dofData != nullptr){
        size_t length = (endId - beginId);
        for (size_t i = 0; i < length; ++i) {
            (*dofData)(i) = sol[beginId + i];
        }
    }
}

const std::shared_ptr<Eigen::VectorXd> &DOF_Base::getDofData() const {
    return dofData;
}

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

    friend std::ostream &operator<<(std::ostream &os, const DOF<Dimension, ScalarDOF> &dof) {
        os << "----------" << std::endl;
        os << "DOF name: " << dof.getDOFName() <<std::endl;
        for (size_t i = 0; i < dof.size(); ++i) {
            os << dof.getDofData()->data()[i] << std::endl;
        }
        os << "----------" << std::endl;
        return os;
    }

    void writeFile(IO_FileWriter &IO, std::shared_ptr<CartesianVector<Dimension, double>> &node){
        std::ofstream fs;
        std::string fileName, resultsPath, filePostFix;
        resultsPath = IO.getProjectResultsPath();

        if(this->size() == node->getSize()){
            fileName = this->getDOFName();
            filePostFix = ".txt";
            fs.open(resultsPath+fileName+filePostFix);
            for (size_t i = 0; i < node->getSize(); ++i) {
                for (size_t d = 0; d < Dimension::Dim; ++d) {
                    fs << node->getDataComponent(d)->at(i) << "\t";
                }
                fs << dofData->data()[i];
                fs << std::endl;
            }
            fs.close();
        }
        else return;
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
