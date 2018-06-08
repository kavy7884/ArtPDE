
#ifndef ARTCFD_FUNCTIONSPACE_HPP
#define ARTCFD_FUNCTIONSPACE_HPP

// Std Lib Include Zone
#include <unordered_map>

// ArtPDE Lib Include Zone

#include "ShapeFunctionFactory.h"

namespace art_pde{

	/*template<>
	class FunctionSpace{
	public:

	private:
		NonZeroBasisId int2basis;
		std::shared_ptr<GeometryData<zzz>> mesh;
		ShapeFunctionFactory<zzz> shape;
	};
    */

	class IntegrationMethod;

	template<class Unit>
	class IntUnit{
	private:
		
		std::size_t global_id_;
		std::shared_ptr<IntegrationMethod> method;
	};
}

#endif 


