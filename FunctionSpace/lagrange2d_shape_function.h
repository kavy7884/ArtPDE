
#ifndef ARTCFD_LAGRANGE2DSHAPEFUNCTION_H
#define ARTCFD_LAGRANGE2DSHAPEFUNCTION_H

// ArtPDE Lib Include Zone
#include "ShapeFunction.h"


namespace art_pde{


	template<>
	class ShapeFunction< Dim2D, Q4, Lagrange > :
		public LagrangeType<Dim2D>
	{

	public:

		using PointType = Point<Dim2D, CartesianCoordinate>;

		friend class SingletonHolder<ShapeFunction< Dim2D, Q4, Lagrange > >;

		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			N_[0] = 0.25*(1 - xi)*(1 - eta);
			N_[1] = 0.25*(1 + xi)*(1 - eta);
			N_[2] = 0.25*(1 + xi)*(1 + eta);
			N_[3] = 0.25*(1 - xi)*(1 + eta);

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			dNdxi_[0][0] = -0.25*(1 - eta);
			dNdxi_[0][1] = 0.25*(1 - eta);
			dNdxi_[0][2] = 0.25*(1 + eta);
			dNdxi_[0][3] = -0.25*(1 + eta);

			dNdxi_[1][0] = -0.25*(1 - xi);
			dNdxi_[1][1] = -0.25*(1 + xi);
			dNdxi_[1][2] = 0.25*(1 + xi);
			dNdxi_[1][3] = 0.25*(1 - xi);

			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 4; ++j){
				for (int i = 0; i < 2; ++i){
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();

			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();

			Jacobian_[0][0] = 0.25*((1 - eta)*(x2 - x1) + (1 + eta)*(x3 - x4));
			Jacobian_[1][0] = 0.25*((1 - xi)*(x4 - x1) + (1 + xi)*(x3 - x2));
			Jacobian_[0][1] = 0.25*((1 - eta)*(y2 - y1) + (1 + eta)*(y3 - y4));
			Jacobian_[1][1] = 0.25*((1 - xi)*(y4 - y1) + (1 + xi)*(y3 - y2));

			return Jacobian_;
		}

		virtual double evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_Jacobian(iso_point, elem_nodes);
			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] - 
				            Jacobian_[0][1] * Jacobian_[1][0];
			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_detJacobian(iso_point, elem_nodes);
			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*Jacobian_[1][1];
			inv_Jacobian_[0][1] = -(1. / det_Jacobian_)*Jacobian_[0][1];
			inv_Jacobian_[1][0] = -(1. / det_Jacobian_)*Jacobian_[1][0];
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*Jacobian_[0][0];

			return inv_Jacobian_;
		}

	private:
		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 2x4
		std::vector<std::vector<double>> dNdx_;  // 2x4
		std::vector<std::vector<double>> Jacobian_; // 2x2
		std::vector<std::vector<double>> inv_Jacobian_; // 2x2
		double det_Jacobian_;

		ShapeFunction() :
			N_(4),
			dNdxi_(2, std::vector<double>(4)),
			dNdx_(2, std::vector<double>(4)),
			Jacobian_(2, std::vector<double>(2)),
			inv_Jacobian_(2, std::vector<double>(2)){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};


	template<>
	class ShapeFunction< Dim2D, Q9, Lagrange > :
		public LagrangeType<Dim2D>
	{
	public:

		using PointType = Point<Dim2D, CartesianCoordinate>;

		friend class SingletonHolder<ShapeFunction< Dim2D, Q9, Lagrange > >;

		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			N_[0] = 0.25 * (1. - xi)*(1. - eta)*xi*eta;
			N_[1] = 0.25*(1. + xi)*(1. - eta)*(-xi*eta);
			N_[2] = 0.25*(1. + xi)*(1. + eta)*xi*eta;
			N_[3] = 0.25*(1. - xi)*(1. + eta)*(-xi*eta);
			N_[4] = 0.5* (1. - xi * xi)*(1. - eta)*(-eta);
			N_[5] = 0.5* (1. + xi)*(1. - eta * eta)*(xi);
			N_[6] = 0.5* (1. - xi * xi)*(1. + eta)*eta;
			N_[7] = 0.5* (1. - xi)*(1. - eta * eta)*(-xi);
			N_[8] = (1. - xi * xi)*(1. - eta * eta);

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			dNdxi_[0][0] = 0.25*(eta - eta * eta - 2. * xi*eta + 2. * xi*eta * eta);
			dNdxi_[0][1] = 0.25*(-eta + eta * eta - 2. * xi*eta + 2. * xi*eta * eta);
			dNdxi_[0][2] = 0.25*(eta + eta * eta + 2. * xi*eta + 2. * xi*eta * eta);
			dNdxi_[0][3] = 0.25*(-eta - eta * eta + 2. * xi*eta + 2. * xi*eta * eta);
			dNdxi_[0][4] = 0.5 * (2. * xi*eta - 2. * xi*eta * eta);
			dNdxi_[0][5] = 0.5 * (1. - eta * eta + 2. * xi - 2. * xi*eta * eta);
			dNdxi_[0][6] = 0.5 * (-2. * xi*eta - 2. * xi*eta * eta);
			dNdxi_[0][7] = 0.5 * (-1. + 2. * xi + eta * eta - 2. * xi*eta * eta);
			dNdxi_[0][8] = -2.0 * xi*(1. - eta * eta);


			dNdxi_[1][0] = 0.25*(xi - xi * xi - 2. * xi*eta + 2. * eta*xi * xi);
			dNdxi_[1][1] = 0.25*(-xi - xi * xi + 2. * xi*eta + 2. * eta*xi * xi);
			dNdxi_[1][2] = 0.25*(xi + xi * xi + 2. * xi*eta + 2. * eta*xi * xi);
			dNdxi_[1][3] = 0.25*(-xi + xi * xi - 2. * xi*eta + 2. * eta*xi * xi);
			dNdxi_[1][4] = 0.5 * (-1. + xi * xi + 2. * eta - 2. * eta*xi * xi);
			dNdxi_[1][5] = 0.5 * (-2. * xi*eta - 2. * eta*xi * xi);
			dNdxi_[1][6] = 0.5 * (1. - xi * xi + 2. * eta - 2. * eta*xi * xi);
			dNdxi_[1][7] = 0.5 * (2. * xi*eta - 2. * eta*xi * xi);
			dNdxi_[1][8] = -2.0 * eta*(1. - xi * xi);

			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 9; ++j){
				for (int i = 0; i < 2; ++i){
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_dNdxi(iso_point);
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();

			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();
			const double& x5 = elem_nodes[4].getX();
			const double& x6 = elem_nodes[5].getX();
			const double& x7 = elem_nodes[6].getX();
			const double& x8 = elem_nodes[7].getX();
			const double& x9 = elem_nodes[8].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();
			const double& y5 = elem_nodes[4].getY();
			const double& y6 = elem_nodes[5].getY();
			const double& y7 = elem_nodes[6].getY();
			const double& y8 = elem_nodes[7].getY();
			const double& y9 = elem_nodes[8].getY();

			Jacobian_[0][0] = 
				2. * x9*xi*(eta * eta - 1.) - 
				(x8*(eta * eta - 1.)*(4. * xi - 2.)) / 4. - 
				(x6*(eta * eta - 1.)*(4. * xi + 2.)) / 4. + 
				(eta*x1*(2. * xi - 1.)*(eta - 1.)) / 4. + 
				(eta*x2*(2. * xi + 1.)*(eta - 1.)) / 4. + 
				(eta*x3*(2. * xi + 1.)*(eta + 1.)) / 4. + 
				(eta*x4*(2. * xi - 1.)*(eta + 1.)) / 4. - 
				eta*x5*xi*(eta - 1.) - 
				eta*x7*xi*(eta + 1.);

			Jacobian_[1][0] = 
				2. * xi*y9*(eta * eta - 1.) - 
				(y8*(eta * eta - 1.)*(4. * xi - 2.)) / 4. - 
				(y6*(eta * eta - 1.)*(4. * xi + 2.)) / 4. + 
				(eta*y1*(2. * xi - 1.)*(eta - 1.)) / 4. + 
				(eta*y2*(2. * xi + 1.)*(eta - 1.)) / 4. + 
				(eta*y3*(2. * xi + 1.)*(eta + 1.)) / 4. +
				(eta*y4*(2. * xi - 1.)*(eta + 1.)) / 4. - 
				eta*xi*y5*(eta - 1.) -
				eta*xi*y7*(eta + 1.);

			Jacobian_[0][1] = 
				2. * eta*x9*(xi * xi - 1.) - 
				(x5*(2. * eta - 1.)*(2. * xi * xi - 2.)) / 4. - 
				(x7*(2. * eta + 1.)*(2. * xi * xi - 2.)) / 4. + 
				(x1*xi*(2. * eta - 1.)*(xi - 1.)) / 4. + 
				(x2*xi*(2. * eta - 1.)*(xi + 1.)) / 4. + 
				(x3*xi*(2. * eta + 1.)*(xi + 1.)) / 4. + 
				(x4*xi*(2. * eta + 1.)*(xi - 1.)) / 4. - 
				eta*x6*xi*(xi + 1.) - eta*x8*xi*(xi - 1.);

			Jacobian_[1][1] = 2. * eta*y9*(xi * xi - 1.) - 
				(y5*(2. * eta - 1.)*(2. * xi * xi - 2.)) / 4. - 
				(y7*(2. * eta + 1.)*(2. * xi * xi - 2.)) / 4. + 
				(xi*y1*(2. * eta - 1.)*(xi - 1.)) / 4. + 
				(xi*y2*(2. * eta - 1.)*(xi + 1.)) / 4. + 
				(xi*y3*(2. * eta + 1.)*(xi + 1.)) / 4. + 
				(xi*y4*(2. * eta + 1.)*(xi - 1.)) / 4. - 
				eta*xi*y6*(xi + 1.) - 
				eta*xi*y8*(xi - 1.);
			return Jacobian_;
		}

		virtual double evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_Jacobian(iso_point, elem_nodes);
			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] -
				Jacobian_[0][1] * Jacobian_[1][0];
			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_detJacobian(iso_point, elem_nodes);
			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*Jacobian_[1][1];
			inv_Jacobian_[0][1] = -(1. / det_Jacobian_)*Jacobian_[0][1];
			inv_Jacobian_[1][0] = -(1. / det_Jacobian_)*Jacobian_[1][0];
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*Jacobian_[0][0];

			return inv_Jacobian_;
		}

	private:

		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 2x9
		std::vector<std::vector<double>> dNdx_;  // 2x9
		std::vector<std::vector<double>> Jacobian_; // 2x2
		std::vector<std::vector<double>> inv_Jacobian_; // 2x2
		double det_Jacobian_;

		ShapeFunction() :
			N_(9),
			dNdxi_(2, std::vector<double>(9)),
			dNdx_(2, std::vector<double>(9)),
			Jacobian_(2, std::vector<double>(2)),
			inv_Jacobian_(2, std::vector<double>(2)){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<>
	class ShapeFunction< Dim2D, T3, Lagrange > :
		public LagrangeType<Dim2D>
	{
	public:

		using PointType = Point<Dim2D, CartesianCoordinate>;

		friend class SingletonHolder < ShapeFunction< Dim2D, T3, Lagrange > >;

		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			N_[0] = 1 - xi - eta;
			N_[1] = xi;
			N_[2] = eta;

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			dNdxi_[0][0] = -1.0;
			dNdxi_[0][1] = 1.0;
			dNdxi_[0][2] = 0.0;

			dNdxi_[1][0] = -1.0;
			dNdxi_[1][1] = 0.0;
			dNdxi_[1][2] = 1.0;

			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 3; ++j){
				for (int i = 0; i < 2; ++i){
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1];
				}
			}
			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_dNdxi(iso_point);

			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();

			Jacobian_[0][0] = x2 - x1;
			Jacobian_[1][0] = y2 - y1;
			Jacobian_[0][1] = x3 - x1;
			Jacobian_[1][1] = y3 - y1;

			return Jacobian_;
		}

		virtual double evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_Jacobian(iso_point, elem_nodes);
			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] -
				Jacobian_[0][1] * Jacobian_[1][0];
			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_detJacobian(iso_point, elem_nodes);
			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*Jacobian_[1][1];
			inv_Jacobian_[0][1] = -(1. / det_Jacobian_)*Jacobian_[0][1];
			inv_Jacobian_[1][0] = -(1. / det_Jacobian_)*Jacobian_[1][0];
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*Jacobian_[0][0];

			return inv_Jacobian_;
		}

	private:

		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 2x3
		std::vector<std::vector<double>> dNdx_;  // 2x3
		std::vector<std::vector<double>> Jacobian_; // 2x2
		std::vector<std::vector<double>> inv_Jacobian_; // 2x2
		double det_Jacobian_;

		ShapeFunction() :
			N_(3),
			dNdxi_(2, std::vector<double>(3)),
			dNdx_(2, std::vector<double>(3)),
			Jacobian_(2, std::vector<double>(2)),
			inv_Jacobian_(2, std::vector<double>(2)){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};


	
}


#endif //ARTCFD_LAGRANGE2DSHAPEFUNCTION_H



