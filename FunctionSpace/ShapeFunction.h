
#ifndef ARTCFD_SHAPEFUNCTION_HPP
#define ARTCFD_SHAPEFUNCTION_HPP

// Std Lib Include Zone
#include <vector>
#include <memory>
#include <cmath>

// ArtPDE Lib Include Zone
#include "Typelist.h"
#include "../Utility/dimension_utility.hpp"
#include "../Utility/numerical_method_utility.hpp"
#include "../Utility/shape_function_name_utility.h"
#include "../Utility/element_type_utility.h"
#include "../Geometry/Point.hpp"

// External Lib Include Zone
#include "/Users/Jeting/Documents/eigen/Eigen/Dense"

namespace art_pde{

	template<class T>
	class SingletonHolder{
	public:
		static T& instance(){
			static T instance_;
			return instance_;
		}
	private:
		SingletonHolder(){}
		SingletonHolder(const SingletonHolder&){}
		SingletonHolder& operator=(const SingletonHolder&){}
	};

	template<class Tlist>
	class ShapeFunction;

	template<class DimensionT, class ElementTypeT>
	class ShapeFunction<TypeList3(DimensionT, ElementTypeT, LagrangePoly)>{
	public:
		virtual std::vector<double>&
			evaluate_shape(const Point<DimensionT, IsoparametricCoordinate>& point) = 0;

		virtual std::vector<std::vector<double>>&
			evaluate_grad(const Point<DimensionT, IsoparametricCoordinate>& point) = 0;

	protected:

		ShapeFunction(){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
		virtual ~ShapeFunction(){}

	};

	template<class DimensionT, class ElementTypeT>
	using LagrangePolyBase = ShapeFunction<TypeList3(DimensionT, ElementTypeT, LagrangePoly)>;

	template<>
	class ShapeFunction< TypeList3(Dim2D, Q4, LagrangePoly) > :
		public LagrangePolyBase<Dim2D, ElementType>
	{
	public:

		friend class SingletonHolder<ShapeFunction< TypeList3(Dim2D, Q4, LagrangePoly) > >;

		virtual std::vector<double>&
			evaluate_shape(const Point<Dim2D, IsoparametricCoordinate>& point) override
		{
			xi = point.getXi();
			eta = point.getEta();
			shape_[0] = 0.25*(1 - xi)*(1 - eta);
			shape_[1] = 0.25*(1 + xi)*(1 - eta);
			shape_[2] = 0.25*(1 + xi)*(1 + eta);
			shape_[3] = 0.25*(1 - xi)*(1 + eta);

			return shape_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_grad(const Point<Dim2D, IsoparametricCoordinate>& point) override
		{
			xi = point.getXi();
			eta = point.getEta();
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
		std::vector<double> shape_;
		std::vector<std::vector<double>> grad_;

		ShapeFunction() :shape_(4), grad_(2, std::vector<double>(4)){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<>
	class ShapeFunction< TypeList3(Dim2D, Q8, LagrangePoly) > :
		public LagrangePolyBase<Dim2D, ElementType>
	{
	public:

		friend class SingletonHolder<ShapeFunction< TypeList3(Dim2D, Q8, LagrangePoly) > >;

		virtual std::vector<double>&
			evaluate_shape(const Point<Dim2D, IsoparametricCoordinate>& point) override
		{
			xi = point.getXi();
			eta = point.getEta();
			shape_[0] = 0.25*(1. - xi)*(1. - eta)*(-xi - eta - 1.);
			shape_[1] = 0.25*(1. + xi)*(1. - eta)*(xi - eta - 1.);
			shape_[2] = 0.25*(1. + xi)*(1. + eta)*(xi + eta - 1.);
			shape_[3] = 0.25*(1. - xi)*(1. + eta)*(eta - xi - 1.);
			shape_[4] = 0.5 * (1. - xi * xi)*(1. - eta);
			shape_[5] = 0.5 * (1. + xi)*(1. - eta * eta);
			shape_[6] = 0.5 * (1. - xi * xi)*(1. + eta);
			shape_[7] = 0.5 * (1. - xi)*(1. - eta * eta);

			return shape_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_grad(const Point<Dim2D, IsoparametricCoordinate>& point) override
		{
			xi = point.getXi();
			eta = point.getEta();

			grad_[0][0] = 0.25*(-eta * eta + eta - 2. * xi*eta + 2. * xi);
			grad_[0][1] = 0.25*(-eta + 2. * xi - 2. * xi*eta + eta * eta);
			grad_[0][2] = 0.25*(eta + 2. * xi + 2. * xi*eta + eta * eta);
			grad_[0][3] = 0.25*(-eta + 2. * xi + 2. * xi*eta - eta * eta);
			grad_[0][4] = -xi*(1. - eta);
			grad_[0][5] = 0.5 * (1. - eta * eta);
			grad_[0][6] = -xi*(1. + eta);
			grad_[0][7] = -0.5 * (1. - eta * eta);

			grad_[1][0] = 0.25*(2. * eta - 2. * xi*eta + xi - xi * xi);
			grad_[1][1] = 0.25*(-xi + 2. * eta - xi * xi + 2. * xi*eta);
			grad_[1][2] = 0.25*(xi + 2. * eta + xi * xi + 2. * xi*eta);
			grad_[1][3] = 0.25*(-xi + 2. * eta + xi * xi - 2. * xi*eta);
			grad_[1][4] = -0.5 * (1. - xi * xi);
			grad_[1][5] = -eta*(1. + xi);
			grad_[1][6] = 0.5 * (1. - xi * xi);
			grad_[1][7] = -eta*(1. - xi);

			return grad_;
		}

	private:
		double xi;
		double eta;
		std::vector<double> shape_;
		std::vector<std::vector<double>> grad_;

		ShapeFunction() :shape_(8), grad_(2, std::vector<double>(8)){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<>
	class ShapeFunction< TypeList3(Dim2D, Q9, LagrangePoly) > :
		public LagrangePolyBase<Dim2D, ElementType>
	{
	public:

		friend class SingletonHolder<ShapeFunction< TypeList3(Dim2D, Q9, LagrangePoly) > >;

		virtual std::vector<double>&
			evaluate_shape(const Point<Dim2D, IsoparametricCoordinate>& point) override
		{
			xi = point.getXi();
			eta = point.getEta();
			shape_[0] = 0.25 * (1. - xi)*(1. - eta)*xi*eta;
			shape_[1] = 0.25*(1. + xi)*(1. - eta)*(-xi*eta);
			shape_[2] = 0.25*(1. + xi)*(1. + eta)*xi*eta;
			shape_[3] = 0.25*(1. - xi)*(1. + eta)*(-xi*eta);
			shape_[4] = 0.5* (1. - xi * xi)*(1. - eta)*(-eta);
			shape_[5] = 0.5* (1. + xi)*(1. - eta * eta)*(xi);
			shape_[6] = 0.5* (1. - xi * xi)*(1. + eta)*eta;
			shape_[7] = 0.5* (1. - xi)*(1. - eta * eta)*(-xi);
			shape_[8] = (1. - xi * xi)*(1. - eta * eta);

			return shape_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_grad(const Point<Dim2D, IsoparametricCoordinate>& point) override
		{
			xi = point.getXi();
			eta = point.getEta();

			grad_[0][0] = 0.25*(eta - eta * eta - 2. * xi*eta + 2. * xi*eta * eta);
			grad_[0][1] = 0.25*(-eta + eta * eta - 2. * xi*eta + 2. * xi*eta * eta);
			grad_[0][2] = 0.25*(eta + eta * eta + 2. * xi*eta + 2. * xi*eta * eta);
			grad_[0][3] = 0.25*(-eta - eta * eta + 2. * xi*eta + 2. * xi*eta * eta);
			grad_[0][4] = 0.5 * (2. * xi*eta - 2. * xi*eta * eta);
			grad_[0][5] = 0.5 * (1. - eta * eta + 2. * xi - 2. * xi*eta * eta);
			grad_[0][6] = 0.5 * (-2. * xi*eta - 2. * xi*eta * eta);
			grad_[0][7] = 0.5 * (-1. + 2. * xi + eta * eta - 2. * xi*eta * eta);
			grad_[0][8] = -2.0 * xi*(1. - eta * eta);

				
			grad_[1][0] = 0.25*(xi - xi * xi - 2. * xi*eta + 2. * eta*xi * xi);
			grad_[1][1] = 0.25*(-xi - xi * xi + 2. * xi*eta + 2. * eta*xi * xi);
			grad_[1][2] = 0.25*(xi + xi * xi + 2. * xi*eta + 2. * eta*xi * xi);
			grad_[1][3] = 0.25*(-xi + xi * xi - 2. * xi*eta + 2. * eta*xi * xi);
			grad_[1][4] = 0.5 * (-1. + xi * xi + 2. * eta - 2. * eta*xi * xi);
			grad_[1][5] = 0.5 * (-2. * xi*eta - 2. * eta*xi * xi);
			grad_[1][6] = 0.5 * (1. - xi * xi + 2. * eta - 2. * eta*xi * xi);
			grad_[1][7] = 0.5 * (2. * xi*eta - 2. * eta*xi * xi);
			grad_[1][8] = -2.0 * eta*(1. - xi * xi);

			return grad_;
		}

	private:
		double xi;
		double eta;
		std::vector<double> shape_;
		std::vector<std::vector<double>> grad_;

		ShapeFunction() :shape_(9), grad_(2, std::vector<double>(9)){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<>
	class ShapeFunction< TypeList3(Dim2D, T3, LagrangePoly) > :
		public LagrangePolyBase<Dim2D, ElementType>
	{
	public:

		friend class SingletonHolder < ShapeFunction< TypeList3(Dim2D, T3, LagrangePoly)> >;

		virtual std::vector<double>&
			evaluate_shape(const Point<Dim2D, IsoparametricCoordinate>& point) override
		{
			xi = point.getXi();
			eta = point.getEta();
			shape_[0] = 1 - xi - eta;
			shape_[1] = xi;
			shape_[2] = eta;

			return shape_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_grad(const Point<Dim2D, IsoparametricCoordinate>& point) override
		{
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
		std::vector<double> shape_;
		std::vector<std::vector<double>> grad_;

		ShapeFunction() :shape_(3), grad_(2, std::vector<double>(3)){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};


	template<>
	class ShapeFunction< TypeList3(Dim3D, Hexa8, LagrangePoly) > :
		public LagrangePolyBase<Dim3D, ElementType>
	{
	public:

		friend class SingletonHolder < ShapeFunction< TypeList3(Dim3D, Hexa8, LagrangePoly)> >;

		virtual std::vector<double>&
			evaluate_shape(const Point<Dim3D, IsoparametricCoordinate>& point) override
		{
			xi = point.getXi();
			eta = point.getEta();
			zeta = point.getZeta();

			shape_[0] = (1 / 8.)*(1 - xi)*(1 - eta)*(1 - zeta);
			shape_[1] = (1 / 8.)*(1 + xi)*(1 - eta)*(1 - zeta);
			shape_[2] = (1 / 8.)*(1 + xi)*(1 + eta)*(1 - zeta);
			shape_[3] = (1 / 8.)*(1 - xi)*(1 + eta)*(1 - zeta);
			shape_[4] = (1 / 8.)*(1 - xi)*(1 - eta)*(1 + zeta);
			shape_[5] = (1 / 8.)*(1 + xi)*(1 - eta)*(1 + zeta);
			shape_[6] = (1 / 8.)*(1 + xi)*(1 + eta)*(1 + zeta);
			shape_[7] = (1 / 8.)*(1 - xi)*(1 + eta)*(1 + zeta);

			return shape_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_grad(const Point<Dim3D, IsoparametricCoordinate>& point) override
		{
			grad_[0][0] = -(1 / 8.)*(1 - eta)*(1 - zeta);
			grad_[0][1] = (1 / 8.)*(1 - eta)*(1 - zeta);
			grad_[0][2] = (1 / 8.)*(1 + eta)*(1 - zeta);
			grad_[0][3] = -(1 / 8.)*(1 + eta)*(1 - zeta);
			grad_[0][4] = -(1 / 8.)*(1 - eta)*(1 + zeta);
			grad_[0][5] = (1 / 8.)*(1 - eta)*(1 + zeta);
			grad_[0][6] = (1 / 8.)*(1 + eta)*(1 + zeta);
			grad_[0][7] = -(1 / 8.)*(1 + eta)*(1 + zeta);

			grad_[1][0] = -(1 / 8.)*(1 - xi)*(1 - zeta);
			grad_[1][1] = -(1 / 8.)*(1 + xi)*(1 - zeta);
			grad_[1][2] = (1 / 8.)*(1 + xi)*(1 - zeta);
			grad_[1][3] = (1 / 8.)*(1 - xi)*(1 - zeta);
			grad_[1][4] = -(1 / 8.)*(1 - xi)*(1 + zeta);
			grad_[1][5] = -(1 / 8.)*(1 + xi)*(1 + zeta);
			grad_[1][6] = (1 / 8.)*(1 + xi)*(1 + zeta);
			grad_[1][7] = (1 / 8.)*(1 - xi)*(1 + zeta);

			grad_[2][0] = -(1 / 8.)*(1 - xi)*(1 - eta);
			grad_[2][1] = -(1 / 8.)*(1 + xi)*(1 - eta);
			grad_[2][2] = -(1 / 8.)*(1 + xi)*(1 + eta);
			grad_[2][3] = -(1 / 8.)*(1 - xi)*(1 + eta);
			grad_[2][4] = (1 / 8.)*(1 - xi)*(1 - eta);
			grad_[2][5] = (1 / 8.)*(1 + xi)*(1 - eta);
			grad_[2][6] = (1 / 8.)*(1 + xi)*(1 + eta);
			grad_[2][7] = (1 / 8.)*(1 - xi)*(1 + eta);

			return grad_;
		}

	private:

		double xi;
		double eta;
		double zeta;
		std::vector<double> shape_;
		std::vector<std::vector<double>> grad_;

		ShapeFunction() :shape_(8), grad_(3, std::vector<double>(8)){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};



	//////////////////////////////////////////////////////////Meshfreeeee////////////////////////////////////////////////////////////


	template<class DimensionT>
	class ShapeFunction<TypeList2(DimensionT, ScatterPoint)>{
	protected:
		// some reusable function for RBFs

		virtual ~ShapeFunction(){}
		ShapeFunction(){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<class DimensionT>
	using RBF = ShapeFunction<TypeList2(DimensionT, ScatterPoint)>;

	template<>
	class ShapeFunction< TypeList3(Dim2D, ScatterPoint, Multiquadric) > :
		public RBF<Dim2D>
	{
	public:

		using THIS = ShapeFunction< TypeList3(Dim2D, ScatterPoint, Multiquadric) >;
		using PoinT = Point<Dim2D, CartesianCoordinate>;
		using PtrPoinT = std::shared_ptr<PoinT>;

		friend class SingletonHolder<THIS>;

		Eigen::RowVectorXd  evaluate_shape(const PtrPoinT center, const std::vector<PtrPoinT>& support)
		{
			std::size_t length = support.size();
			std::vector<double> weight(length);
			Eigen::MatrixXd rbf(length, length);
			for (std::size_t i = 0; i < length; ++i){
				for (std::size_t j = 0; j < length; ++j){
					rbf(i, j) = mq(norm(support[i], support[j]));
				}
			}

			Eigen::RowVectorXd rbfi(length);
			for (std::size_t j = 0; j < length; ++j){
				rbfi(j) = mq(center, support[j]);
			}
			return rbfi * (rbf.inverse());
		}

		std::vector<std::vector<double>> evaluate_grad(const Point<Dim2D, CartesianCoordinate>& point)
		{
			return std::vector<std::vector<double>>();
		}

		std::vector<double> evaluate_laplace(const Point<Dim2D, CartesianCoordinate>& point)
		{
			return std::vector<double>();
		}

		THIS& setC(double c_){
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
		double c{ 1.0 };
		ShapeFunction(){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};


}


#endif //ARTCFD_SHAPEFUNCTION_HPP










