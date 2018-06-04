//
// Created by Chingkai Chou on 5/30/18.
//

#ifndef ARTPDE_KAVY_GEOMETRY_HPP
#define ARTPDE_KAVY_GEOMETRY_HPP

#include <memory>
#include "Project/art_project.hpp"
#include "geometric_data.hpp"
#include "geometric_algorithm.hpp"
#include "geometric_reader.hpp"

namespace art_pde{

    template<typename ...T>
    class Geometry: public virtual T...{
        struct GeometryDummy {};
    public:
        explicit Geometry(GeometryDummy): T()...{};

        static std::shared_ptr<Geometry> create(){
            return make_shared_ArtProject();
        }

    protected:
        static std::shared_ptr<Geometry> make_shared_ArtProject() {
            return std::make_shared<Geometry>(GeometryDummy());
        }

    private:
        Geometry() = delete;
    };
}


#endif //ARTPDE_KAVY_GEOMETRY_HPP
