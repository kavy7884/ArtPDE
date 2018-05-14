
#ifndef ARTCFD_TYPELIST_HPP
#define ARTCFD_TYPELIST_HPP

// Std Lib Include Zone

// ArtPDE Lib Include Zone


namespace art_pde{

	struct Nulltype{};

	template<class T, class U>
	struct Typelist{
		using Start = T;
		using End = U;
	};

	# define TypeList1(T1) Typelist<T1,Nulltype>
	# define TypeList2(T1, T2) Typelist<T1, TypeList1(T2)>
	# define TypeList3(T1, T2, T3) Typelist<T1, TypeList2(T2, T3)>

	template<class T, template<class> class Unit>
	class Hierarchy :public Unit<T>{
		using Base = Unit<T>;
	};

}


#endif //ARTCFD_TYPELIST_HPP

