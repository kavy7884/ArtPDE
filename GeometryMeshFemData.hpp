//
// Created by Chingkai Chou on 2/23/18.
//

#ifndef ARTCFD_GEOMETRYFEM_HPP
#define ARTCFD_GEOMETRYFEM_HPP

#include "GeometryMeshData.hpp"

template <class Dimension>
class GeometryMeshFemData : public GeometryMeshData<Dimension>{
public:
    GeometryMeshFemData():GeometryMeshData<Dimension>(){};
    virtual bool readFile(IO_FileReader &IO_in) override;
    virtual bool writeFile(IO_FileWriter &IO_out)override{return true;};

};

template<class Dimension>
bool GeometryMeshFemData<Dimension>::readFile(IO_FileReader &IO_in) {
    return GeometryMeshData<Dimension>::readFile(IO_in);
}


#endif //ARTCFD_GEOMETRYFEM_HPP
