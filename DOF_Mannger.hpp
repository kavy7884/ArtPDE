//
// Created by Chingkai Chou on 3/12/18.
//

#ifndef ARTCFD_DOF_MANNGER_HPP
#define ARTCFD_DOF_MANNGER_HPP
#include <vector>
#include <memory>
#include "Dof_New.hpp"

class DOF_Mannger{
public:
    DOF_Mannger() {
        groupId.push_back(0);
        systemId.push_back(0);
    }
    bool addDofData(std::shared_ptr<DOF_Base> dof);
    size_t &getTotalDof(){ return tailId;}

    const std::vector<size_t> &getGroupId() const;

    const std::vector<size_t> &getSystemId() const;

    const std::vector<std::shared_ptr<DOF_Base>> &getDofGroupData() const;


private:
    std::vector<std::shared_ptr<DOF_Base>> dofGroupData;
    std::vector<size_t> groupId, systemId;
    size_t tailId{0};

};

bool DOF_Mannger::addDofData(std::shared_ptr<DOF_Base> dof) {
    auto tailBeforeAdd = tailId;
    auto groupNum = dof->getGroupNum();
    for (size_t i = 0; i < groupNum; ++i) {
        tailId += dof->getLengthDOF();
        groupId.push_back(tailId);
    }
    systemId.push_back(tailId);
    dofGroupData.push_back(dof);
    return ((tailId - tailBeforeAdd) == 0);
}

const std::vector<size_t> &DOF_Mannger::getGroupId() const {
    return groupId;
}

const std::vector<size_t> &DOF_Mannger::getSystemId() const {
    return systemId;
}

const std::vector<std::shared_ptr<DOF_Base>> &DOF_Mannger::getDofGroupData() const {
    return dofGroupData;
}

#endif //ARTCFD_DOF_MANNGER_HPP
