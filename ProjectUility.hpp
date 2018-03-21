//
// Created by Chingkai Chou on 2/24/18.
//

#ifndef ARTCFD_PROJECTUILITY_HPP
#define ARTCFD_PROJECTUILITY_HPP

#include <string>
#include <memory>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

class ArtProjectBuilder;

class ArtProject
{
    struct ArtProjectDummy {};
public:
    static ArtProjectBuilder create(const std::string & projectName);

    const std::string &getRunPath() const {
        return runPath;
    }

    const std::string &getProjectName() const {
        return projectName;
    }

    const std::string getProjectPath() const {
        std::string out = runPath + slash + projectName + slash;
        return out;
    }

    const std::string &getProjectGeometryFolderName() const {
        return projectGeometryFolderName;
    }

    const std::string getProjectGeometryPath() const {
        std::string out = getProjectPath() + projectGeometryFolderName + slash;
        return out;
    }

    const std::string &getProjectSettingFolderName() const {
        return projectSettingFolderName;
    }

    const std::string getProjectSettingPath() const {
        std::string out = getProjectPath() + projectSettingFolderName + slash;
        return out;
    }

    const std::string &getProjectInitialFolderName() const {
        return projectInitialFolderName;
    }

    const std::string getProjectInitialPath() const {
        std::string out = getProjectPath() + projectInitialFolderName + slash;
        return out;
    }

    const std::string &getProjectResultsFolderName() const {
        return projectResultsFolderName;
    }

    const std::string getProjectResultsPath() const {
        std::string out = getProjectPath() + projectResultsFolderName + slash;
        return out;
    }

    const std::string &getSlash() const {
        return slash;
    }

    bool Init();

    static std::shared_ptr<ArtProject> make_shared_ArtProject() {
        return std::make_shared<ArtProject>(ArtProjectDummy());
    }

    explicit ArtProject(ArtProjectDummy) {};

private:
    friend class ArtProjectBuilder;
    ArtProject() = delete;

    void setRunPath(const std::string &runPath_) {
        runPath = runPath_;
    }

    void setProjectName(const std::string &projectName) {
        ArtProject::projectName = projectName;
    }

    void setProjectGeometryFolderName(const std::string &projectGeometryFolderName) {
        ArtProject::projectGeometryFolderName = projectGeometryFolderName;
    }

    void setProjectSettingFolderName(const std::string &projectSettingFolderName) {
        ArtProject::projectSettingFolderName = projectSettingFolderName;
    }

    void setProjectInitialFolderName(const std::string &projectInitialFolderName) {
        ArtProject::projectInitialFolderName = projectInitialFolderName;
    }

    void setProjectResultsFolderName(const std::string &projectResultsFolderName) {
        ArtProject::projectResultsFolderName = projectResultsFolderName;
    }

    void setSlash(const std::string &slash) {
        ArtProject::slash = slash;
    }

    bool checkFolder(const std::string &folderPath);

    std::string runPath{"."};
    std::string projectName{"ArtProject"};
    std::string projectGeometryFolderName{"Geometry"};
    std::string projectSettingFolderName{"Setting"};
    std::string projectInitialFolderName{"Initial"};
    std::string projectResultsFolderName{"Results"};
    std::string slash{"/"};
};

class ArtProjectBuilder
{
public:
    ArtProjectBuilder(const std::string & projectName){
        p_project = ArtProject::make_shared_ArtProject();
        p_project->setProjectName(projectName);
    }

    std::shared_ptr<ArtProject> build(){
        p_project->Init();
        //TODO: Error Exception
        return p_project;
    }

    ArtProjectBuilder &setRunPath(const std::string &runPath){
        p_project->setRunPath(runPath);
        return *this;
    }

    ArtProjectBuilder & setGeometryFolderName(const std::string &geometryFolderName){
        p_project->setProjectGeometryFolderName(geometryFolderName);
        return *this;
    }

    ArtProjectBuilder & setSettingFolderName(const std::string &settingFolderName){
        p_project->setProjectSettingFolderName(settingFolderName);
        return *this;
    }

    ArtProjectBuilder & setInitialFolderName(const std::string &initialFolderName){
        p_project->setProjectInitialFolderName(initialFolderName);
        return *this;
    }

    ArtProjectBuilder & setResultsFolderName(const std::string &resultsFolderName){
        p_project->setProjectResultsFolderName(resultsFolderName);
        return *this;
    }

    ArtProjectBuilder & setDivideSlash(const std::string &slash){
        p_project->setSlash(slash);
        return *this;
    }

private:
    std::shared_ptr<ArtProject> p_project;
};

ArtProjectBuilder ArtProject::create(const std::string & projectName){
    return ArtProjectBuilder(projectName);
}

bool ArtProject::checkFolder(const std::string &folderPath) {

    //Check folder exist or not!
    if ((access(folderPath.c_str(), 0)) == -1)
    {
        if(mkdir(folderPath.c_str(), 0777) == 0){
            return true;
        }else{
            return false;
        }
    }
    return true;
}

bool ArtProject::Init() {
    if(!checkFolder(getProjectPath())) return false;
    if(!checkFolder(getProjectGeometryPath())) return false;
    if(!checkFolder(getProjectInitialPath())) return false;
    if(!checkFolder(getProjectSettingPath())) return false;
    if(!checkFolder(getProjectResultsPath())) return false;

    return true;
}


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
