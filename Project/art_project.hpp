//
// Created by Chingkai Chou on 6/1/18.
//

#ifndef ARTPDE_KAVY_PROJECT_HPP
#define ARTPDE_KAVY_PROJECT_HPP

#include <iostream>
#include <string>
#include <memory>

#if defined( _MSC_VER )
#include <direct.h>
#include <windows.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

namespace art_pde{ namespace project {

        class ArtProjectBuilder;

        class ArtProject
        {
            struct ArtProjectDummy {};
        public:
            explicit ArtProject(ArtProjectDummy) {};

            static ArtProjectBuilder create(const std::string & projectName);

            static std::shared_ptr<ArtProject> make_shared_ArtProject() {
                return std::make_shared<ArtProject>(ArtProjectDummy());
            }

            bool Init();

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

            int checkFolder(const std::string &folderPath);

            std::string runPath{"."};
            std::string projectName{"ArtPDE"};
            std::string projectGeometryFolderName{"Geometry"};
            std::string projectSettingFolderName{"Setting"};
            std::string projectInitialFolderName{"Initial"};
            std::string projectResultsFolderName{"Results"};
			#if defined( _MSC_VER )
			std::string slash{ '\\' };
			#else
			std::string slash{ '/' };
			#endif
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

        #include "src/art_project_impl.cpp"
}}

#endif //ARTPDE_KAVY_PROJECT_HPP
