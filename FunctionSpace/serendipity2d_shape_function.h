
#ifndef ARTCFD_SERENDIPITY2DSHAPEFUNCTION_H
#define ARTCFD_SERENDIPITY2DSHAPEFUNCTION_H

// ArtPDE Lib Include Zone
#include "ShapeFunction.h"


namespace art_pde{


	template<>
	class ShapeFunction< Dim2D, Q8, Lagrange > :
		public LagrangeType<Dim2D>
	{
	public:

		using PointType = Point<Dim2D, CartesianCoordinate>;

		friend class SingletonHolder<ShapeFunction< Dim2D, Q9, Lagrange > >;

		virtual std::vector<double>&
			N(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			N_[0] = 0.25*(1. - xi)*(1. - eta)*(-xi - eta - 1.);
			N_[1] = 0.25*(1. + xi)*(1. - eta)*(xi - eta - 1.);
			N_[2] = 0.25*(1. + xi)*(1. + eta)*(xi + eta - 1.);
			N_[3] = 0.25*(1. - xi)*(1. + eta)*(eta - xi - 1.);
			N_[4] = 0.5 * (1. - xi * xi)*(1. - eta);
			N_[5] = 0.5 * (1. + xi)*(1. - eta * eta);
			N_[6] = 0.5 * (1. - xi * xi)*(1. + eta);
			N_[7] = 0.5 * (1. - xi)*(1. - eta * eta);

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			dNdxi(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			dNdxi_[0][0] = 0.25*(-eta * eta + eta - 2. * xi*eta + 2. * xi);
			dNdxi_[0][1] = 0.25*(-eta + 2. * xi - 2. * xi*eta + eta * eta);
			dNdxi_[0][2] = 0.25*(eta + 2. * xi + 2. * xi*eta + eta * eta);
			dNdxi_[0][3] = 0.25*(-eta + 2. * xi + 2. * xi*eta - eta * eta);
			dNdxi_[0][4] = -xi*(1. - eta);
			dNdxi_[0][5] = 0.5 * (1. - eta * eta);
			dNdxi_[0][6] = -xi*(1. + eta);
			dNdxi_[0][7] = -0.5 * (1. - eta * eta);

			dNdxi_[1][0] = 0.25*(2. * eta - 2. * xi*eta + xi - xi * xi);
			dNdxi_[1][1] = 0.25*(-xi + 2. * eta - xi * xi + 2. * xi*eta);
			dNdxi_[1][2] = 0.25*(xi + 2. * eta + xi * xi + 2. * xi*eta);
			dNdxi_[1][3] = 0.25*(-xi + 2. * eta + xi * xi - 2. * xi*eta);
			dNdxi_[1][4] = -0.5 * (1. - xi * xi);
			dNdxi_[1][5] = -eta*(1. + xi);
			dNdxi_[1][6] = 0.5 * (1. - xi * xi);
			dNdxi_[1][7] = -eta*(1. - xi);

			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			dNdxi(iso_point);
			invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 8; ++j){
				for (int i = 0; i < 2; ++i){
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			dNdxi(iso_point);
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

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();
			const double& y5 = elem_nodes[4].getY();
			const double& y6 = elem_nodes[5].getY();
			const double& y7 = elem_nodes[6].getY();
			const double& y8 = elem_nodes[7].getY();

			Jacobian_[0][0] = 
				(x8-x6)*(eta * eta / 2. - 0.5) + 
				x5*xi*(eta - 1.) - 
				x7*xi*(eta + 1.) - 
				(x1*(eta + 2. * xi)*(eta - 1.)) / 4. + 
				(x2*(eta - 2. * xi)*(eta - 1.)) / 4. + 
				(x3*(eta + 2. * xi)*(eta + 1.)) / 4. -
				(x4*(eta - 2. * xi)*(eta + 1.)) / 4.;

			Jacobian_[1][0] = 
				(y8-y6)*(eta * eta / 2. - 0.5) + 
				xi*y5*(eta - 1.) - 
				xi*y7*(eta + 1.) - 
				(y1*(eta + 2. * xi)*(eta - 1.)) / 4. + 
				(y2*(eta - 2. * xi)*(eta - 1.)) / 4. + 
				(y3*(eta + 2. * xi)*(eta + 1.)) / 4. - 
				(y4*(eta - 2. * xi)*(eta + 1.)) / 4.;

			Jacobian_[0][1] = 
				(x5-x7)*(xi * xi / 2. - 0.5) -
				eta*x6*(xi + 1.) + 
				eta*x8*(xi - 1.) - 
				(x1*(2. * eta + xi)*(xi - 1.)) / 4. + 
				(x3*(2. * eta + xi)*(xi + 1.)) / 4. + 
				(x2*(xi + 1.)*(2. * eta - xi)) / 4. - 
				(x4*(xi - 1.)*(2. * eta - xi)) / 4.;

			Jacobian_[1][1] = 
				(y5-y7)*(xi * xi / 2. - 0.5) - 
				eta*y6*(xi + 1.) + 
				eta*y8*(xi - 1.) - 
				(y1*(2. * eta + xi)*(xi - 1.)) / 4. + 
				(y3*(2. * eta + xi)*(xi + 1.)) / 4. + 
				(y2*(xi + 1.)*(2. * eta - xi)) / 4. - 
				(y4*(xi - 1.)*(2. * eta - xi)) / 4.;
		}

		virtual double detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			Jacobian(iso_point, elem_nodes);
			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] -
				Jacobian_[0][1] * Jacobian_[1][0];
			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			detJacobian(iso_point, elem_nodes);
			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*Jacobian_[1][1];
			inv_Jacobian_[0][1] = -(1. / det_Jacobian_)*Jacobian_[0][1];
			inv_Jacobian_[1][0] = -(1. / det_Jacobian_)*Jacobian_[1][0];
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*Jacobian_[0][0];
		}

	private:

		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 2x8
		std::vector<std::vector<double>> dNdx_;  // 2x8
		std::vector<std::vector<double>> Jacobian_; // 2x2
		std::vector<std::vector<double>> inv_Jacobian_; // 2x2
		double det_Jacobian_;

		ShapeFunction() :
			N_(8),
			dNdxi_(2, std::vector<double>(8)),
			dNdx_(2, std::vector<double>(8)),
			Jacobian_(2, std::vector<double>(2)),
			inv_Jacobian_(2, std::vector<double>(2)){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};


}


#endif //ARTCFD_SERENDIPITY2DSHAPEFUNCTION_H



