
#ifndef ARTCFD_SHAPEFUNCTION_HPP
#define ARTCFD_SHAPEFUNCTION_HPP

// Std Lib Include Zone
#include <vector>
#include <memory>
#include <cmath>

// ArtPDE Lib Include Zone
#include "dimension_utility.hpp"
#include "numerical_method_utility.hpp"
#include "shape_function_utility.h"
#include "element_type_utility.h"
#include "Point.hpp"
#include "Eigen/Dense"

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

	template<class ...>
	class ShapeFunction;

	template<class DimensionT>
	class ShapeFunction<DimensionT, ElementType, Lagrange>{
		using PointType = Point<DimensionT, CartesianCoordinate>;
	public:
		virtual std::vector<double>& 
			evaluate_shape(const PointType& iso_point) = 0;

		virtual std::vector<std::vector<double>>& 
			evaluate_dNdxi(const PointType& iso_point) = 0;

		virtual std::vector<std::vector<double>>& 
			evaluate_dNdx (const PointType& iso_point, const std::vector<PointType>& elem_nodes) = 0;

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) = 0;

		virtual double
			evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) = 0;

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) = 0;

	protected:
		ShapeFunction(){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
		virtual ~ShapeFunction(){}
	};

	template<class DimensionT>
	class ShapeFunction<DimensionT, ElementType, Serendipity>{
		using PointType = Point<DimensionT, CartesianCoordinate>;
	public:
		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) = 0;

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) = 0;

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) = 0;

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) = 0;

		virtual double
			evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) = 0;

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) = 0;

	protected:

		ShapeFunction(){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
		virtual ~ShapeFunction(){}

	};

	template<class DimensionT>
	using LagrangeType = ShapeFunction<DimensionT, ElementType, Lagrange>;

	template<class DimensionT>
	using SerendipityType = ShapeFunction<DimensionT, ElementType, Serendipity>;


	//////////////////////////////////////////////////////////Meshfreeeee////////////////////////////////////////////////////////////


	template<class DimensionT>
	class ShapeFunction<DimensionT, ScatterPoint>{
	protected:
		// some reusable function for RBFs

		virtual ~ShapeFunction(){}
		ShapeFunction(){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<class DimensionT>
	using RBF = ShapeFunction<DimensionT, ScatterPoint>;

	template<>
	class ShapeFunction< Dim2D, ScatterPoint, Multiquadric > :
		public RBF<Dim2D>
	{
	public:

		using THIS = ShapeFunction< Dim2D, ScatterPoint, Multiquadric >;
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










