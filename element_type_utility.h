
#ifndef ARTCFD_ELEMENT_TYPE_UTILITY_HPP
#define ARTCFD_ELEMENT_TYPE_UTILITY_HPP

// Std Lib Include Zone

// ArtPDE Lib Include Zone


namespace art_pde{

	struct ElementType{};

	struct LinearQuadrilateral :public ElementType{ static const int NUM = 4; };

	struct LinearTriangle :public ElementType{ static const int NUM = 3; };

	struct ScatterPoint{};

}




#endif //ARTCFD_ELEMENT_TYPE_UTILITY_HPP
