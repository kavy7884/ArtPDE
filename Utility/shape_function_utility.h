
#ifndef ARTCFD_SHAPE_FUNCTION_NAME_UTILITY_HPP
#define ARTCFD_SHAPE_FUNCTION_NAME_UTILITY_HPP

// Std Lib Include Zone

// ArtPDE Lib Include Zone


namespace art_pde{

	struct FunctionType{};

	struct Lagrange :public FunctionType{};

	struct Serendipity :public FunctionType{};

	struct Multiquadric :public FunctionType{};
}


#endif //ARTCFD_SHAPE_FUNCTION_NAME_UTILITY_HPP

