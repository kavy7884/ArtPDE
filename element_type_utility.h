
#ifndef ARTCFD_ELEMENT_TYPE_UTILITY_HPP
#define ARTCFD_ELEMENT_TYPE_UTILITY_HPP

// Std Lib Include Zone

// ArtPDE Lib Include Zone


namespace art_pde{

	struct ElementType{};

	// 2D

	struct Q4 :public ElementType{ static const int NUM = 4; };

	struct Q8 :public ElementType{ static const int NUM = 8; };

	struct Q9 :public ElementType{ static const int NUM = 9; };

	struct T3 :public ElementType{ static const int NUM = 3; };

	struct T6 :public ElementType{ static const int NUM = 6; };

	// 3D

	struct Hexa8 :public ElementType{ static const int NUM = 8; };

	struct Tetra4 :public ElementType{ static const int NUM = 4; };

	struct Prism5 :public ElementType{ static const int NUM = 5; };

	struct ScatterPoint{};

}




#endif //ARTCFD_ELEMENT_TYPE_UTILITY_HPP
