//
// Created by Chingkai Chou on 2/23/18.
//

#ifndef ARTCFD_IOBASE_HPP
#define ARTCFD_IOBASE_HPP

#include <string>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "ProjectUility.hpp"

enum class IO_State{ReadIn, WriteOut};

class IO_Basic{
public:
    IO_Basic(ProjectArt &proj) : proj(proj) {
        //TODO - error exception
        init();
    }
    std::string getProjectRootPath();
    std::string getProjectMeshPath();
    std::string getProjectAnalysisPath();
    std::string getProjectResultsPath();

    bool CheckFolder(std::string folderPath);

    virtual bool init();

    template <typename T>
    static bool splitLineString(std::string bufferLine, std::vector<T> &outLine);

protected:
    ProjectArt & proj;
};

std::string IO_Basic::getProjectRootPath() {
    return std::string(proj.getRunPath() + proj.getProjectName() + "/");
}

std::string IO_Basic::getProjectMeshPath() {
    return std::string(proj.getRunPath() + proj.getProjectName() + "/" + proj.getProjectMeshFolderName() + "/");
}

std::string IO_Basic::getProjectAnalysisPath() {
    return std::string(proj.getRunPath() + proj.getProjectName() + "/" + proj.getProjectAnalysisFolderName() + "/");
}

std::string IO_Basic::getProjectResultsPath() {
    return std::string(proj.getRunPath() + proj.getProjectName() + "/" + proj.getProjectResultsFolderName() + "/");
}

bool IO_Basic::CheckFolder(std::string folderPath) {
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

bool IO_Basic::init() {
    if(CheckFolder(getProjectRootPath()))
        return true;
    else
        return false;
}

template<typename T>
bool IO_Basic::splitLineString(std::string bufferLine, std::vector<T> &outLine) {
    std::istringstream os(bufferLine);
    T temp;
    outLine.clear();
    while (os>>temp){
        outLine.push_back(temp);
    }
    return true;
}

class IO_FileMode : public IO_Basic{
public:
    IO_FileMode(ProjectArt &proj) : IO_Basic(proj) {
        //TODO - error exception
        init();
    }

    virtual bool init() override ;

};

bool IO_FileMode::init() {
    if(CheckFolder(getProjectMeshPath()) && CheckFolder(getProjectAnalysisPath()) && CheckFolder(getProjectResultsPath()))
        return true;
    else
        return false;
}

class IO_FileReader : public IO_FileMode{
public:
    IO_FileReader(ProjectArt &proj) : IO_FileMode(proj){}

protected:
    IO_State state{IO_State::ReadIn};


};

class IO_FileWriter : public IO_FileMode{
public:
    IO_FileWriter(ProjectArt &proj) : IO_FileMode(proj){}

protected:
    IO_State state{IO_State::WriteOut};
};




#endif //ARTCFD_IOBASE_HPP
