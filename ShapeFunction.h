
#ifndef ARTCFD_SHAPEFUNCTION_HPP
#define ARTCFD_SHAPEFUNCTION_HPP

// Std Lib Include Zone
#include <array>
#include <vector>
#include <memory>
#include <cmath>

// ArtPDE Lib Include Zone
#include "Typelist.h"
#include "dimension_utility.hpp"
#include "numerical_method_utility.hpp"
#include "shape_function_name_utility.h"
#include "element_type_utility.h"
#include "Point.hpp"
#include "Eigen/Dense"

namespace art_pde{

	template<class T>
	class SingletonBase{
	public:
		static T& instance(){
			static T instance_;
			return instance_;
		}
	protected:
		SingletonBase(){}
		SingletonBase(const T&){}
		SingletonBase& operator=(const T&){}
	};

	template<class Tlist>
	class ShapeFunction{};

	template<>
	class ShapeFunction< TypeList3(Dim2D, Multiquadric, ScatterPoint) >{

	public:

		using PoinT = Point<Dim2D, CartesianCoordinate>;
		using PtrPoinT = std::shared_ptr<PoinT>;

		Eigen::RowVectorXd evaluate_shape(const PtrPoinT center, const std::vector<PtrPoinT>& support){

			std::size_t length = support.size();
			std::vector<double> weight(length);
			Eigen::MatrixXd rbf(length, length);
			for (int i = 0; i < length; ++i){
				for (int j = 0; j < length; ++j){
					rbf(i, j) = mq(norm(support[i], support[j]));
				}
			}
			
			Eigen::RowVectorXd rbfi(length);
			for (int j = 0; j < length; ++j){
				rbfi(j) = mq(center, support[j]);
			}
			return rbfi * (rbf.inverse());
		}

		ShapeFunction< TypeList3(Dim2D, Multiquadric, ScatterPoint) >& setC(double c_){
			c = c_;
			return *this;
		}

	private:

		double norm(PtrPoinT current, PtrPoinT another){
			double dx = current->getX() - another->getX();
			double dy = current->getY() - another->getY();
			return sqrt(dx*dx + dy*dy);
		}
		double mq(double radius){
			return sqrt(radius*radius + c*c);
		}
		double mq(PtrPoinT current, PtrPoinT another){
			double dx = current->getX() - another->getX();
			double dy = current->getY() - another->getY();
			return sqrt(dx*dx + dy*dy + c*c);
		}
		double c = 1.0;

	};

	template<>
	class ShapeFunction< TypeList3(Dim2D, LagrangePoly, LinearQuadrilateral) >
	{
	public:

		std::array<double, 4>& evaluate_shape(const Point<Dim2D, CartesianCoordinate>& point){

			xi = point.getX();
			eta = point.getY();
			shape_[0] = 0.25*(1 - xi)*(1 - eta);
			shape_[1] = 0.25*(1 + xi)*(1 - eta);
			shape_[2] = 0.25*(1 + xi)*(1 + eta);
			shape_[3] = 0.25*(1 - xi)*(1 + eta);

			return shape_;
		}

		std::array<std::array<double, 2>, 4>& evaluate_grad(const Point<Dim2D, CartesianCoordinate>& point){

			xi = point.getX();
			eta = point.getY();
			grad_[0][0] = -0.25*(1 - eta);
			grad_[0][1] = 0.25*(1 - eta);
			grad_[0][2] = 0.25*(1 + eta);
			grad_[0][3] = -0.25*(1 + eta);

			grad_[1][0] = -0.25*(1 - xi);
			grad_[1][1] = -0.25*(1 + xi);
			grad_[1][2] = 0.25*(1 + xi);
			grad_[1][3] = 0.25*(1 - xi);

			return grad_;
		}

	private:

		double xi;
		double eta;
		std::array<double, 4> shape_;
		std::array<std::array<double, 2>, 4> grad_;

	};

	template<>
	class ShapeFunction< TypeList3(Dim2D, LagrangePoly, LinearTriangle) > {

	public:

		std::array<double, 3>& evaluate_shape(Point<Dim2D, CartesianCoordinate>& point){

			xi = point.getX();
			eta = point.getY();
			shape_[0] = 1 - xi - eta;
			shape_[1] = xi;
			shape_[2] = eta;

			return shape_;
		}

		std::array<std::array<double, 2>, 3>& evaluate_grad(Point<Dim2D, CartesianCoordinate>& point){

			grad_[0][0] = -1.0;
			grad_[0][1] = 1.0;
			grad_[0][2] = 0.0;

			grad_[1][0] = -1.0;
			grad_[1][1] = 0.0;
			grad_[1][2] = 1.0;

			return grad_;
		}

	private:

		double xi;
		double eta;
		std::array<double, 3> shape_;
		std::array<std::array<double, 2>, 3> grad_;

	};


	template<class Tlist>
	class ShapeFunctionFactory{};


	template<class DimensionT, class FunctionNameT>
	class ShapeFunctionFactory< TypeList3(DimensionT, FunctionNameT, ElementType) >{

		template <class ElementT>
		using ReturnType = ShapeFunction< TypeList3(DimensionT, FunctionNameT, ElementT) >;

	public:

		template <class ElementT>
		ReturnType<ElementT>& getInstance(ElementT&){
			return Hierarchy< ReturnType<ElementT>, SingletonBase> ::instance();
		}
	};

	template<class DimensionT, class FunctionNameT>
	class ShapeFunctionFactory< TypeList3(DimensionT, FunctionNameT, ScatterPoint) >{

		using ReturnType = ShapeFunction< TypeList3(DimensionT, FunctionNameT, ScatterPoint) >;

	public:

		ReturnType& getInstance(){
			return Hierarchy< ReturnType, SingletonBase> ::instance();
		}
	};


}


/*namespace art_pde
{

using Position = Point<Dim1D, CartesianCoordinate>;

template<typename T>
using MatrixDynamic = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

class ShapeFunctionType{
public:
std::string element_type;
std::string function_name;
int derivatives_order;

bool operator==(const ShapeFunctionType& type){
return (element_type == type.element_type
&& function_name == type.function_name
&& derivatives_order == type.derivatives_order);
}
};

template<typename T>
class ShapeFunction{
public:
MatrixDynamic<T> operator()(const Position& position){
return this->shape(position);
}
protected:
virtual MatrixDynamic<T> shape(const Position& position) = 0;
};

template<typename T>
class IsoLinear2DQuadrilateral :public ShapeFunction<T>{
protected:
MatrixDynamic<T> shape(const Position& position){
T xi = position[0];
T eta = position[1];
MatrixDynamic<T> temp(1, 4);
temp.setZero();
temp(0, 0) = 0.25*(1 - xi)*(1 - eta);
temp(0, 1) = 0.25*(1 + xi)*(1 - eta);
temp(0, 2) = 0.25*(1 + xi)*(1 + eta);
temp(0, 3) = 0.25*(1 - xi)*(1 + eta);
return temp;
}
};

template<typename T>
class IsoLinear2DQuadrilateral_dNdx :public ShapeFunction<T>{
protected:
MatrixDynamic<T> shape(const Position& position){
T xi = position[0];
T eta = position[1];
MatrixDynamic<T> temp(2, 4);
temp.setZero();
temp(0, 0) = -0.25*(1 - eta);
temp(0, 1) = 0.25*(1 - eta);
temp(0, 2) = 0.25*(1 + eta);
temp(0, 3) = -0.25*(1 + eta);

temp(1, 0) = -0.25*(1 - xi);
temp(1, 1) = -0.25*(1 + xi);
temp(1, 2) = 0.25*(1 + xi);
temp(1, 3) = 0.25*(1 - xi);
return temp;
}
};

template<typename T>
class IsoLinear3DHexa :public ShapeFunction<T>{
protected:
MatrixDynamic<T> shape(const Position& position){
T xi = coordinates[0];
T eta = coordinates[1];
T zeta = coordinates[2];
MatrixDynamic<T> temp(1, 8);
temp.setZero();
temp(0, 0) = (1 / 8.)*(1 - xi)*(1 - eta)*(1 - zeta);
temp(0, 1) = (1 / 8.)*(1 + xi)*(1 - eta)*(1 - zeta);
temp(0, 2) = (1 / 8.)*(1 + xi)*(1 + eta)*(1 - zeta);
temp(0, 3) = (1 / 8.)*(1 - xi)*(1 + eta)*(1 - zeta);
temp(0, 4) = (1 / 8.)*(1 - xi)*(1 - eta)*(1 + zeta);
temp(0, 5) = (1 / 8.)*(1 + xi)*(1 - eta)*(1 + zeta);
temp(0, 6) = (1 / 8.)*(1 + xi)*(1 + eta)*(1 + zeta);
temp(0, 7) = (1 / 8.)*(1 - xi)*(1 + eta)*(1 + zeta);
return temp;
}
};

}

namespace std{
template<>
struct hash<art_pde::ShapeFunctionType>{
std::size_t operator()(const art_pde::ShapeFunctionType& type) const{
return(std::hash<std::string>()(type.element_type)
^ (std::hash<std::string>()(type.function_name) << 1) >> 1)
^ (std::hash<int>()(type.derivatives_order) << 1);
}
};
}*/

#endif //ARTCFD_SHAPEFUNCTION_HPP

/*
std::unordered_map<ShapeFunctionType, ShapeFunction<double>> table = {
{ { "Quadrilateral", "Lagrange", 0 },  IsoLinear2DQuadrilateral<double>() },
{ { "Quadrilateral", "Lagrange", 1 }, IsoLinear2DQuadrilateral_dNdx<double>() },
{ { "Hexahedron", "Lagrange", 0 }, IsoLinear3DHexa<double>() }
};
*/