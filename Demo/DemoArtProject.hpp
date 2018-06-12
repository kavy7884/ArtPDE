//
// Created by Chingkai Chou on 6/1/18.
//

#ifndef ARTPDE_KAVY_DEMOARTPROJECT_HPP
#define ARTPDE_KAVY_DEMOARTPROJECT_HPP

#include "Project/art_project.hpp"

void DemoArtProject() {
    using namespace art_pde::project;
    auto proj_1 = ArtProject::create("TestProj").setRunPath(".").setInitialFolderName("Init").build();

    auto proj_2 = ArtProject::create("OutSideProj").setRunPath("./Proj").build();
}

#endif //ARTPDE_KAVY_DEMOARTPROJECT_HPP
