//
// Created by Chingkai Chou on 3/12/18.
//

#ifndef ARTCFD_DOF_MANNGER_HPP
#define ARTCFD_DOF_MANNGER_HPP
#include <vector>
#include <memory>
#include "DOF.hpp"

class DOF_Mannger{
public:
    DOF_Mannger() {
        groupId.push_back(0);
        systemId.push_back(0);
    }
    bool addDofData(std::shared_ptr<DOF_Base> dof);
    size_t &getTotalDof(){ return tailId;}


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
    return ((tailId - tailBeforeAdd) == 0);
}

#endif //ARTCFD_DOF_MANNGER_HPP
