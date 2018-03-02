//
// Created by Chingkai Chou on 2/24/18.
//

#ifndef ARTCFD_PROJECTUILITY_HPP
#define ARTCFD_PROJECTUILITY_HPP

#include <string>
#include "IO_Uility.hpp"

class ProjectArt{
    std::string runPath{"./"};
    std::string projectName{};
    std::string projectMeshFolderName{"Mesh/"};
    std::string projectAnalysisFolderName{"Analysis/"};
    std::string projectResultsFolderName{"Results/"};
public:
    const std::string &getProjectMeshFolderName() const {
        return projectMeshFolderName;
    }

    void setProjectMeshFolderName(const std::string &projectMeshFolderName) {
        ProjectArt::projectMeshFolderName = projectMeshFolderName;
    }

    const std::string &getProjectAnalysisFolderName() const {
        return projectAnalysisFolderName;
    }

    void setProjectAnalysisFolderName(const std::string &projectAnalysisFolderName) {
        ProjectArt::projectAnalysisFolderName = projectAnalysisFolderName;
    }

    const std::string &getProjectResultsFolderName() const {
        return projectResultsFolderName;
    }

    void setProjectResultsFolderName(const std::string &projectResultsFolderName) {
        ProjectArt::projectResultsFolderName = projectResultsFolderName;
    }

public:
    const std::string &getProjectName() const {
        return projectName;
    }

    void setProjectName(const std::string &projectName) {
        ProjectArt::projectName = projectName;
    }

    const std::string &getRunPath() const {
        return runPath;
    }

    void setRunPath(const std::string &RunPath) {
        ProjectArt::runPath = RunPath;
    }

public:
    ProjectArt(const std::string &projectName) : projectName(projectName) {}

    bool init();
};

bool ProjectArt::init() {

    return false;
}


#endif //ARTCFD_PROJECTUILITY_HPP
