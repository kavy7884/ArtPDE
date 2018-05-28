
#ifndef ARTCFD_LAGRANGE3DSHAPEFUNCTION_H
#define ARTCFD_LAGRANGE3DSHAPEFUNCTION_H


// ArtPDE Lib Include Zone
#include "ShapeFunction.h"


namespace art_pde{

	template<>
	class ShapeFunction< Dim3D, Tetra4, Lagrange > :
		public LagrangeType<Dim3D>
	{
	public:

		using PointType = Point<Dim3D, CartesianCoordinate>;

		friend class SingletonHolder < ShapeFunction< Dim3D, Tetra4, Lagrange > >;

		virtual std::vector<double>& 
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();
			N_[0] = xi;
			N_[1] = eta;
			N_[2] = zeta;
			N_[3] = 1.0 - xi - eta - zeta;

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			dNdxi_[0][0] = 1.0;
			dNdxi_[0][1] = 0.0;
			dNdxi_[0][2] = 0.0;
			dNdxi_[0][3] = -1.0;

			dNdxi_[1][0] = 0.0;
			dNdxi_[1][1] = 1.0;
			dNdxi_[1][2] = 0.0;
			dNdxi_[1][3] = -1.0;

			dNdxi_[2][0] = 0.0;
			dNdxi_[2][1] = 0.0;
			dNdxi_[2][2] = 1.0;
			dNdxi_[2][3] = -1.0;

			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 4; ++j){
				for (int i = 0; i < 3; ++i){
					dNdx_[i][j] = 
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1] +
						dNdxi_[2][j] * inv_Jacobian_[i][2];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>& 
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();

			const double& z1 = elem_nodes[0].getZ();
			const double& z2 = elem_nodes[1].getZ();
			const double& z3 = elem_nodes[2].getZ();
			const double& z4 = elem_nodes[3].getZ();

			Jacobian_[0][0] = x1 - x4;
			Jacobian_[1][0] = y1 - y4;
			Jacobian_[2][0] = z1 - z4;
			Jacobian_[0][1] = x2 - x4;
			Jacobian_[1][1] = y2 - y4;
			Jacobian_[2][1] = z2 - z4;
			Jacobian_[0][2] = x3 - x4;
			Jacobian_[1][2] = y3 - y4;
			Jacobian_[2][2] = z3 - z4;
			return Jacobian_;
		}

		virtual double 
			evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();

			const double& z1 = elem_nodes[0].getZ();
			const double& z2 = elem_nodes[1].getZ();
			const double& z3 = elem_nodes[2].getZ();
			const double& z4 = elem_nodes[3].getZ();
			det_Jacobian_ =
				x1*y2*z3 - x1*y3*z2 - x2*y1*z3 +
				x2*y3*z1 + x3*y1*z2 - x3*y2*z1 -
				x1*y2*z4 + x1*y4*z2 + x2*y1*z4 -
				x2*y4*z1 - x4*y1*z2 + x4*y2*z1 +
				x1*y3*z4 - x1*y4*z3 - x3*y1*z4 +
				x3*y4*z1 + x4*y1*z3 - x4*y3*z1 -
				x2*y3*z4 + x2*y4*z3 + x3*y2*z4 -
				x3*y4*z2 - x4*y2*z3 + x4*y3*z2;

			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>& 
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();

			const double& z1 = elem_nodes[0].getZ();
			const double& z2 = elem_nodes[1].getZ();
			const double& z3 = elem_nodes[2].getZ();
			const double& z4 = elem_nodes[3].getZ();
			evaluate_detJacobian(iso_point, elem_nodes);
			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*((y2 - y4)*(z3 - z4) - (y3 - y4)*(z2 - z4));
			inv_Jacobian_[1][0] = (1. / det_Jacobian_)*((y3 - y4)*(z1 - z4) - (y1 - y4)*(z3 - z4));
			inv_Jacobian_[2][0] = (1. / det_Jacobian_)*((y1 - y4)*(z2 - z4) - (y2 - y4)*(z1 - z4));
			inv_Jacobian_[0][1] = (1. / det_Jacobian_)*((x3 - x4)*(z2 - z4) - (x2 - x4)*(z3 - z4));
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*((x1 - x4)*(z3 - z4) - (x3 - x4)*(z1 - z4));
			inv_Jacobian_[2][1] = (1. / det_Jacobian_)*((x2 - x4)*(z1 - z4) - (x1 - x4)*(z2 - z4));
			inv_Jacobian_[0][2] = (1. / det_Jacobian_)*((x2 - x4)*(y3 - y4) - (x3 - x4)*(y2 - y4));
			inv_Jacobian_[1][2] = (1. / det_Jacobian_)*((x3 - x4)*(y1 - y4) - (x1 - x4)*(y3 - y4));
			inv_Jacobian_[2][2] = (1. / det_Jacobian_)*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4));

			return inv_Jacobian_;
		}

	private:

		

		
		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 3x4
		std::vector<std::vector<double>> dNdx_;  // 3x4
		std::vector<std::vector<double>> Jacobian_; // 3x3
		std::vector<std::vector<double>> inv_Jacobian_; // 3x3
		double det_Jacobian_;

		ShapeFunction() :
			N_(4),
			dNdxi_(3, std::vector<double>(4)),
			dNdx_ (3, std::vector<double>(4)),
			Jacobian_(3, std::vector<double>(3)),
			inv_Jacobian_(3, std::vector<double>(3)){}
		ShapeFunction(const ShapeFunction&){}
		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<>
	class ShapeFunction< Dim3D, Hexa8, Lagrange > :
		public LagrangeType<Dim3D>
	{
	public:

		using PointType = Point<Dim3D, CartesianCoordinate>;

		friend class SingletonHolder < ShapeFunction< Dim3D, Hexa8, Lagrange > >;
		
		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();
			N_[0] = (1 / 8.)*(1 - xi)*(1 - eta)*(1 - zeta);
			N_[1] = (1 / 8.)*(1 + xi)*(1 - eta)*(1 - zeta);
			N_[2] = (1 / 8.)*(1 + xi)*(1 + eta)*(1 - zeta);
			N_[3] = (1 / 8.)*(1 - xi)*(1 + eta)*(1 - zeta);
			N_[4] = (1 / 8.)*(1 - xi)*(1 - eta)*(1 + zeta);
			N_[5] = (1 / 8.)*(1 + xi)*(1 - eta)*(1 + zeta);
			N_[6] = (1 / 8.)*(1 + xi)*(1 + eta)*(1 + zeta);
			N_[7] = (1 / 8.)*(1 - xi)*(1 + eta)*(1 + zeta);

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();
			dNdxi_[0][0] = -(1 / 8.)*(1 - eta)*(1 - zeta);
			dNdxi_[0][1] = (1 / 8.)*(1 - eta)*(1 - zeta);
			dNdxi_[0][2] = (1 / 8.)*(1 + eta)*(1 - zeta);
			dNdxi_[0][3] = -(1 / 8.)*(1 + eta)*(1 - zeta);
			dNdxi_[0][4] = -(1 / 8.)*(1 - eta)*(1 + zeta);
			dNdxi_[0][5] = (1 / 8.)*(1 - eta)*(1 + zeta);
			dNdxi_[0][6] = (1 / 8.)*(1 + eta)*(1 + zeta);
			dNdxi_[0][7] = -(1 / 8.)*(1 + eta)*(1 + zeta);

			dNdxi_[1][0] = -(1 / 8.)*(1 - xi)*(1 - zeta);
			dNdxi_[1][1] = -(1 / 8.)*(1 + xi)*(1 - zeta);
			dNdxi_[1][2] = (1 / 8.)*(1 + xi)*(1 - zeta);
			dNdxi_[1][3] = (1 / 8.)*(1 - xi)*(1 - zeta);
			dNdxi_[1][4] = -(1 / 8.)*(1 - xi)*(1 + zeta);
			dNdxi_[1][5] = -(1 / 8.)*(1 + xi)*(1 + zeta);
			dNdxi_[1][6] = (1 / 8.)*(1 + xi)*(1 + zeta);
			dNdxi_[1][7] = (1 / 8.)*(1 - xi)*(1 + zeta);

			dNdxi_[2][0] = -(1 / 8.)*(1 - xi)*(1 - eta);
			dNdxi_[2][1] = -(1 / 8.)*(1 + xi)*(1 - eta);
			dNdxi_[2][2] = -(1 / 8.)*(1 + xi)*(1 + eta);
			dNdxi_[2][3] = -(1 / 8.)*(1 - xi)*(1 + eta);
			dNdxi_[2][4] = (1 / 8.)*(1 - xi)*(1 - eta);
			dNdxi_[2][5] = (1 / 8.)*(1 + xi)*(1 - eta);
			dNdxi_[2][6] = (1 / 8.)*(1 + xi)*(1 + eta);
			dNdxi_[2][7] = (1 / 8.)*(1 - xi)*(1 + eta);

			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 8; ++j){
				for (int i = 0; i < 3; ++i){
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1] +
						dNdxi_[2][j] * inv_Jacobian_[i][2];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

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

			const double& z1 = elem_nodes[0].getZ();
			const double& z2 = elem_nodes[1].getZ();
			const double& z3 = elem_nodes[2].getZ();
			const double& z4 = elem_nodes[3].getZ();
			const double& z5 = elem_nodes[4].getZ();
			const double& z6 = elem_nodes[5].getZ();
			const double& z7 = elem_nodes[6].getZ();
			const double& z8 = elem_nodes[7].getZ();

			Jacobian_[0][0] = (1 / 8.)*
				((x2 - x1)*(eta - 1)*(zeta - 1) + 
				(x4 - x3)*(eta + 1)*(zeta - 1) +
				(x5 - x6)*(eta - 1)*(zeta + 1) + 
				(x7 - x8)*(eta + 1)*(zeta + 1));
			Jacobian_[1][0] = (1 / 8.)*
				((y2 - y1)*(eta - 1)*(zeta - 1) +
				(y4 - y3)*(eta + 1)*(zeta - 1) +
				(y5 - y6)*(eta - 1)*(zeta + 1) +
				(y7 - y8)*(eta + 1)*(zeta + 1));
			Jacobian_[2][0] = (1 / 8.)*(
				(z2 - z1)*(eta - 1)*(zeta - 1) + 
				(z4 - z3)*(eta + 1)*(zeta - 1) + 
				(z5 - z6)*(eta - 1)*(zeta + 1) + 
				(z7 - z8)*(eta + 1)*(zeta + 1));
			Jacobian_[0][1] = (1 / 8.)*(
				(x2 - x3)*(xi + 1)*(zeta - 1) + 
				(x4 - x1)*(xi - 1)*(zeta - 1) + 
				(x5 - x8)*(xi - 1)*(zeta + 1) + 
				(x7 - x6)*(xi + 1)*(zeta + 1));
			Jacobian_[1][1] = (1 / 8.)*(
				(y2 - y3)*(xi + 1)*(zeta - 1) + 
				(y4 - y1)*(xi - 1)*(zeta - 1) + 
				(y5 - y8)*(xi - 1)*(zeta + 1) + 
				(y7 - y6)*(xi + 1)*(zeta + 1));
			Jacobian_[2][1] = (1 / 8.)*(
				(z2 - z3)*(xi + 1)*(zeta - 1) + 
				(z4 - z1)*(xi - 1)*(zeta - 1) + 
				(z5 - z8)*(xi - 1)*(zeta + 1) + 
				(z7 - z6)*(xi + 1)*(zeta + 1));
			Jacobian_[0][2] = (1 / 8.)*(
				(x2 - x6)*(eta - 1)*(xi + 1) + 
				(x5 - x1)*(eta - 1)*(xi - 1) + 
				(x7 - x3)*(eta + 1)*(xi + 1) + 
				(x4 - x8)*(eta + 1)*(xi - 1));
			Jacobian_[1][2] = (1 / 8.)*(
				(y2 - y6)*(eta - 1)*(xi + 1) + 
				(y5 - y1)*(eta - 1)*(xi - 1) + 
				(y7 - y3)*(eta + 1)*(xi + 1) + 
				(y4 - y8)*(eta + 1)*(xi - 1));
			Jacobian_[2][2] = (1 / 8.)*(
				(z2 - z6)*(eta - 1)*(xi + 1) + 
				(z5 - z1)*(eta - 1)*(xi - 1) + 
				(z7 - z3)*(eta + 1)*(xi + 1) + 
				(z4 - z8)*(eta + 1)*(xi - 1));
			return Jacobian_;
		}

		virtual double
			evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_Jacobian(iso_point, elem_nodes);
			
			det_Jacobian_ = Jacobian_[0][0]*Jacobian_[1][1]*Jacobian_[2][2] + 
				     Jacobian_[1][0]*Jacobian_[2][1]*Jacobian_[0][2] + 
				     Jacobian_[0][1]*Jacobian_[1][2]*Jacobian_[2][0] - 
				     Jacobian_[0][2]*Jacobian_[1][1]*Jacobian_[2][0] - 
				     Jacobian_[0][1]*Jacobian_[1][0]*Jacobian_[2][2] - 
				     Jacobian_[1][2]*Jacobian_[2][1]*Jacobian_[0][0];
				

			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			
			evaluate_detJacobian(iso_point, elem_nodes);
			double& J11 = Jacobian_[0][0];
			double& J21 = Jacobian_[1][0];
			double& J31 = Jacobian_[2][0];
			double& J12 = Jacobian_[0][1];
			double& J22 = Jacobian_[1][1];
			double& J32 = Jacobian_[2][1];
			double& J13 = Jacobian_[0][2];
			double& J23 = Jacobian_[1][2];
			double& J33 = Jacobian_[2][2];

			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*(J22*J33 - J32*J23);
			inv_Jacobian_[1][0] = (1. / det_Jacobian_)*(J31*J23 - J21*J33);
			inv_Jacobian_[2][0] = (1. / det_Jacobian_)*(J21*J32 - J31*J22);
			inv_Jacobian_[0][1] = (1. / det_Jacobian_)*(J32*J13 - J12*J33);
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*(J11*J33 - J31*J13);
			inv_Jacobian_[2][1] = (1. / det_Jacobian_)*(J31*J12 - J11*J32);
			inv_Jacobian_[0][2] = (1. / det_Jacobian_)*(J12*J23 - J22*J13);
			inv_Jacobian_[1][2] = (1. / det_Jacobian_)*(J21*J13 - J11*J23);
			inv_Jacobian_[2][2] = (1. / det_Jacobian_)*(J11*J22 - J21*J12);

			return inv_Jacobian_;
		}

	private:

		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 3x8
		std::vector<std::vector<double>> dNdx_;  // 3x8
		std::vector<std::vector<double>> Jacobian_; // 3x3
		std::vector<std::vector<double>> inv_Jacobian_; // 3x3
		double det_Jacobian_;

		ShapeFunction() :
			N_(8),
			dNdxi_(3, std::vector<double>(8)),
			dNdx_(3, std::vector<double>(8)),
			Jacobian_(3, std::vector<double>(3)),
			inv_Jacobian_(3, std::vector<double>(3)){}

		ShapeFunction(const ShapeFunction&){}

		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<>
	class ShapeFunction< Dim3D, Pyramid5, Lagrange > :
		public LagrangeType<Dim3D>
	{
	public:

		using PointType = Point<Dim3D, CartesianCoordinate>;

		friend class SingletonHolder < ShapeFunction< Dim3D, Pyramid5, Lagrange > >;

		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();
			N_[0] = (1. / 8.)*(1 - xi)*(1 - eta)*(1 - zeta);
			N_[1] = (1. / 8.)*(1 + xi)*(1 - eta)*(1 - zeta);
			N_[2] = (1. / 8.)*(1 + xi)*(1 + eta)*(1 - zeta);
			N_[3] = (1. / 8.)*(1 - xi)*(1 + eta)*(1 - zeta);
			N_[4] = (1. / 2.)*(1 + zeta);

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();
			dNdxi_[0][0] = -(1 / 8.)*(1 - eta)*(1 - zeta);
			dNdxi_[0][1] = (1 / 8.)*(1 - eta)*(1 - zeta);
			dNdxi_[0][2] = (1 / 8.)*(1 + eta)*(1 - zeta);
			dNdxi_[0][3] = -(1 / 8.)*(1 + eta)*(1 - zeta);
			dNdxi_[0][4] = 0.0;
			
			dNdxi_[1][0] = -(1 / 8.)*(1 - xi)*(1 - zeta);
			dNdxi_[1][1] = -(1 / 8.)*(1 + xi)*(1 - zeta);
			dNdxi_[1][2] = (1 / 8.)*(1 + xi)*(1 - zeta);
			dNdxi_[1][3] = (1 / 8.)*(1 - xi)*(1 - zeta);
			dNdxi_[1][4] = 0.0;
			
			dNdxi_[2][0] = -(1 / 8.)*(1 - xi)*(1 - eta);
			dNdxi_[2][1] = -(1 / 8.)*(1 + xi)*(1 - eta);
			dNdxi_[2][2] = -(1 / 8.)*(1 + xi)*(1 + eta);
			dNdxi_[2][3] = -(1 / 8.)*(1 - xi)*(1 + eta);
			dNdxi_[2][4] = 0.5;
		
			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 5; ++j){
				for (int i = 0; i < 3; ++i){
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1] +
						dNdxi_[2][j] * inv_Jacobian_[i][2];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();
			const double& x5 = elem_nodes[4].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();
			const double& y5 = elem_nodes[4].getY();

			const double& z1 = elem_nodes[0].getZ();
			const double& z2 = elem_nodes[1].getZ();
			const double& z3 = elem_nodes[2].getZ();
			const double& z4 = elem_nodes[3].getZ();
			const double& z5 = elem_nodes[4].getZ();

			Jacobian_[0][0] = (1 / 8.)*((x2 - x1)*(eta - 1)*(zeta - 1) + (x4 - x3)*(eta + 1)*(zeta - 1));
			Jacobian_[1][0] = (1 / 8.)*((y2 - y1)*(eta - 1)*(zeta - 1) + (y4 - y3)*(eta + 1)*(zeta - 1));
			Jacobian_[2][0] = (1 / 8.)*((z2 - z1)*(eta - 1)*(zeta - 1) + (z4 - z3)*(eta + 1)*(zeta - 1));
			Jacobian_[0][1] = (1 / 8.)*((x2 - x3)*(xi + 1)*(zeta - 1) + (x4 - x1)*(xi - 1)*(zeta - 1));
			Jacobian_[1][1] = (1 / 8.)*((y2 - y3)*(xi + 1)*(zeta - 1) + (y4 - y1)*(xi - 1)*(zeta - 1));
			Jacobian_[2][1] = (1 / 8.)*((z2 - z3)*(xi + 1)*(zeta - 1) + (z4 - z1)*(xi - 1)*(zeta - 1));
			Jacobian_[0][2] = (1 / 2.)*x5 + (1 / 8.) * (x2*(eta - 1)*(xi + 1) - x1*(eta - 1)*(xi - 1) - x3*(eta + 1)*(xi + 1) + x4*(eta + 1)*(xi - 1));
			Jacobian_[1][2] = (1 / 2.)*y5 + (1 / 8.) * (y2*(eta - 1)*(xi + 1) - y1*(eta - 1)*(xi - 1) - y3*(eta + 1)*(xi + 1) + y4*(eta + 1)*(xi - 1));
			Jacobian_[2][2] = (1 / 2.)*z5 + (1 / 8.) * (z2*(eta - 1)*(xi + 1) - z1*(eta - 1)*(xi - 1) - z3*(eta + 1)*(xi + 1) + z4*(eta + 1)*(xi - 1));
			return Jacobian_;
		}

		virtual double
			evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_Jacobian(iso_point, elem_nodes);


			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] * Jacobian_[2][2] +
				Jacobian_[1][0] * Jacobian_[2][1] * Jacobian_[0][2] +
				Jacobian_[0][1] * Jacobian_[1][2] * Jacobian_[2][0] -
				Jacobian_[0][2] * Jacobian_[1][1] * Jacobian_[2][0] -
				Jacobian_[0][1] * Jacobian_[1][0] * Jacobian_[2][2] -
				Jacobian_[1][2] * Jacobian_[2][1] * Jacobian_[0][0];


			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_detJacobian(iso_point, elem_nodes);
			double& J11 = Jacobian_[0][0];
			double& J21 = Jacobian_[1][0];
			double& J31 = Jacobian_[2][0];
			double& J12 = Jacobian_[0][1];
			double& J22 = Jacobian_[1][1];
			double& J32 = Jacobian_[2][1];
			double& J13 = Jacobian_[0][2];
			double& J23 = Jacobian_[1][2];
			double& J33 = Jacobian_[2][2];

			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*(J22*J33 - J32*J23);
			inv_Jacobian_[1][0] = (1. / det_Jacobian_)*(J31*J23 - J21*J33);
			inv_Jacobian_[2][0] = (1. / det_Jacobian_)*(J21*J32 - J31*J22);
			inv_Jacobian_[0][1] = (1. / det_Jacobian_)*(J32*J13 - J12*J33);
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*(J11*J33 - J31*J13);
			inv_Jacobian_[2][1] = (1. / det_Jacobian_)*(J31*J12 - J11*J32);
			inv_Jacobian_[0][2] = (1. / det_Jacobian_)*(J12*J23 - J22*J13);
			inv_Jacobian_[1][2] = (1. / det_Jacobian_)*(J21*J13 - J11*J23);
			inv_Jacobian_[2][2] = (1. / det_Jacobian_)*(J11*J22 - J21*J12);

			return inv_Jacobian_;
		}

	private:

		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 3x5
		std::vector<std::vector<double>> dNdx_;  // 3x5
		std::vector<std::vector<double>> Jacobian_; // 3x3
		std::vector<std::vector<double>> inv_Jacobian_; // 3x3
		double det_Jacobian_;

		ShapeFunction() :
			N_(5),
			dNdxi_(3, std::vector<double>(5)),
			dNdx_(3, std::vector<double>(5)),
			Jacobian_(3, std::vector<double>(3)),
			inv_Jacobian_(3, std::vector<double>(3)){}

		ShapeFunction(const ShapeFunction&){}

		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<>
	class ShapeFunction< Dim3D, Prism6, Lagrange > :
		public LagrangeType<Dim3D>
	{
	public:

		using PointType = Point<Dim3D, CartesianCoordinate>;

		friend class SingletonHolder < ShapeFunction< Dim3D, Prism6, Lagrange > >;

		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();
			N_[0] = (1 - xi)*(1 - eta - zeta);
			N_[1] = (1 - xi)*zeta;
			N_[2] = (1 - xi)*eta;
			N_[3] = xi*(1 - eta - zeta);
			N_[4] = xi*zeta;
			N_[5] = xi*eta;

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();
			dNdxi_[0][0] = eta + zeta - 1;
			dNdxi_[0][1] = -zeta;
			dNdxi_[0][2] = -eta;
			dNdxi_[0][3] = 1 - eta - zeta;
			dNdxi_[0][4] = zeta;
			dNdxi_[0][5] = eta;

			dNdxi_[1][0] = xi - 1;
			dNdxi_[1][1] = 0.0;
			dNdxi_[1][2] = 1 - xi;
			dNdxi_[1][3] = -xi;
			dNdxi_[1][4] = 0.0;
			dNdxi_[1][5] = xi;

			dNdxi_[2][0] = xi - 1;
			dNdxi_[2][1] = 1 - xi;
			dNdxi_[2][2] = 0.0;
			dNdxi_[2][3] = -xi;
			dNdxi_[2][4] = xi;
			dNdxi_[2][5] = 0.0;

			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 5; ++j){
				for (int i = 0; i < 3; ++i){
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1] +
						dNdxi_[2][j] * inv_Jacobian_[i][2];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();
			const double& x5 = elem_nodes[4].getX();
			const double& x6 = elem_nodes[5].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();
			const double& y5 = elem_nodes[4].getY();
			const double& y6 = elem_nodes[5].getY();

			const double& z1 = elem_nodes[0].getZ();
			const double& z2 = elem_nodes[1].getZ();
			const double& z3 = elem_nodes[2].getZ();
			const double& z4 = elem_nodes[3].getZ();
			const double& z5 = elem_nodes[4].getZ();
			const double& z6 = elem_nodes[5].getZ();

			Jacobian_[0][0] = (x1 - x4)*(eta + zeta - 1) + (x6 - x3)*eta + (x5 - x2)*zeta;
			Jacobian_[1][0] = (y1 - y4)*(eta + zeta - 1) + (y6 - y3)*eta + (y5 - y2)*zeta;
			Jacobian_[2][0] = (z1 - z4)*(eta + zeta - 1) + (z6 - z3)*eta + (z5 - z2)*zeta;
			Jacobian_[0][1] = (x6 - x4)*xi + (x1 - x3)*(xi - 1);
			Jacobian_[1][1] = (y6 - y4)*xi + (y1 - y3)*(xi - 1);
			Jacobian_[2][1] = (z6 - z4)*xi + (z1 - z3)*(xi - 1);
			Jacobian_[0][2] = (x5 - x4)*xi + (x1 - x2)*(xi - 1);
			Jacobian_[1][2] = (y5 - y4)*xi + (y1 - y2)*(xi - 1);
			Jacobian_[2][2] = (z5 - z4)*xi + (z1 - z2)*(xi - 1);
			return Jacobian_;
		}

		virtual double
			evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_Jacobian(iso_point, elem_nodes);


			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] * Jacobian_[2][2] +
				Jacobian_[1][0] * Jacobian_[2][1] * Jacobian_[0][2] +
				Jacobian_[0][1] * Jacobian_[1][2] * Jacobian_[2][0] -
				Jacobian_[0][2] * Jacobian_[1][1] * Jacobian_[2][0] -
				Jacobian_[0][1] * Jacobian_[1][0] * Jacobian_[2][2] -
				Jacobian_[1][2] * Jacobian_[2][1] * Jacobian_[0][0];


			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_detJacobian(iso_point, elem_nodes);
			double& J11 = Jacobian_[0][0];
			double& J21 = Jacobian_[1][0];
			double& J31 = Jacobian_[2][0];
			double& J12 = Jacobian_[0][1];
			double& J22 = Jacobian_[1][1];
			double& J32 = Jacobian_[2][1];
			double& J13 = Jacobian_[0][2];
			double& J23 = Jacobian_[1][2];
			double& J33 = Jacobian_[2][2];

			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*(J22*J33 - J32*J23);
			inv_Jacobian_[1][0] = (1. / det_Jacobian_)*(J31*J23 - J21*J33);
			inv_Jacobian_[2][0] = (1. / det_Jacobian_)*(J21*J32 - J31*J22);
			inv_Jacobian_[0][1] = (1. / det_Jacobian_)*(J32*J13 - J12*J33);
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*(J11*J33 - J31*J13);
			inv_Jacobian_[2][1] = (1. / det_Jacobian_)*(J31*J12 - J11*J32);
			inv_Jacobian_[0][2] = (1. / det_Jacobian_)*(J12*J23 - J22*J13);
			inv_Jacobian_[1][2] = (1. / det_Jacobian_)*(J21*J13 - J11*J23);
			inv_Jacobian_[2][2] = (1. / det_Jacobian_)*(J11*J22 - J21*J12);

			return inv_Jacobian_;
		}

	private:

		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 3x5
		std::vector<std::vector<double>> dNdx_;  // 3x5
		std::vector<std::vector<double>> Jacobian_; // 3x3
		std::vector<std::vector<double>> inv_Jacobian_; // 3x3
		double det_Jacobian_;

		ShapeFunction() :
			N_(6),
			dNdxi_(3, std::vector<double>(6)),
			dNdx_(3, std::vector<double>(6)),
			Jacobian_(3, std::vector<double>(3)),
			inv_Jacobian_(3, std::vector<double>(3)){}

		ShapeFunction(const ShapeFunction&){}

		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<>
	class ShapeFunction< Dim3D, Tetra10, Lagrange > :
		public LagrangeType<Dim3D>
	{
	public:

		using PointType = Point<Dim3D, CartesianCoordinate>;

		friend class SingletonHolder < ShapeFunction< Dim3D, Tetra10, Lagrange > >;

		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			N_[0] = xi*(2.0 * xi - 1.0);
			N_[1] = eta*(2.0 * eta - 1.0);
			N_[2] = zeta*(2.0 * zeta - 1.0);
			N_[3] = (1.0 - xi - eta - zeta)*(2.0 * (1.0 - xi - eta - zeta) - 1.0);
			N_[4] = 4.0 * xi*eta;
			N_[5] = 4.0 * eta*zeta;
			N_[6] = 4.0 * zeta*xi;
			N_[7] = 4.0 * xi*(1.0 - xi - eta - zeta);
			N_[8] = 4.0 * eta*(1.0 - xi - eta - zeta);
			N_[9] = 4.0 * zeta*(1.0 - xi - eta - zeta);

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			dNdxi_[0][0] = 4. * xi - 1.;
			dNdxi_[0][1] = 0.;
			dNdxi_[0][2] = 0.;
			dNdxi_[0][3] = 4. * (xi + eta + zeta) - 3.;
			dNdxi_[0][4] = 4. * eta;
			dNdxi_[0][5] = 0.;
			dNdxi_[0][6] = 4. * zeta;
			dNdxi_[0][7] = -4. * (2. * xi + eta + zeta - 1.);
			dNdxi_[0][8] = -4. * eta;
			dNdxi_[0][9] = -4. * zeta;
			

			dNdxi_[1][0] = 0.;
			dNdxi_[1][1] = 4. * eta - 1.;
			dNdxi_[1][2] = 0.;
			dNdxi_[1][3] = 4. * (xi + eta + zeta) - 3.;
			dNdxi_[1][4] = 4. * xi;
			dNdxi_[1][5] = 4. * zeta;
			dNdxi_[1][6] = 0.;
			dNdxi_[1][7] = -4. * xi;
			dNdxi_[1][8] = -4. * (xi + 2. * eta + zeta - 1.);
			dNdxi_[1][9] = -4. * zeta;

			dNdxi_[2][0] = 0.;
			dNdxi_[2][1] = 0.;
			dNdxi_[2][2] = 4. * zeta - 1.;
			dNdxi_[2][3] = 4. * (xi + eta + zeta) - 3.;
			dNdxi_[2][4] = 0.;
			dNdxi_[2][5] = 4. * eta;
			dNdxi_[2][6] = 4. * xi;
			dNdxi_[2][7] = -4. * xi;
			dNdxi_[2][8] = -4. * eta;
			dNdxi_[2][9] = -4. * (xi + eta + 2. * zeta - 1.);

			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 10; ++j){
				for (int i = 0; i < 3; ++i){
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1] +
						dNdxi_[2][j] * inv_Jacobian_[i][2];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();
			const double& x5 = elem_nodes[4].getX();
			const double& x6 = elem_nodes[5].getX();
			const double& x7 = elem_nodes[6].getX();
			const double& x8 = elem_nodes[7].getX();
			const double& x9 = elem_nodes[8].getX();
			const double& x10 = elem_nodes[9].getX();
			

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();
			const double& y5 = elem_nodes[4].getY();
			const double& y6 = elem_nodes[5].getY();
			const double& y7 = elem_nodes[6].getY();
			const double& y8 = elem_nodes[7].getY();
			const double& y9 = elem_nodes[8].getY();
			const double& y10 = elem_nodes[9].getY();

			const double& z1 = elem_nodes[0].getZ();
			const double& z2 = elem_nodes[1].getZ();
			const double& z3 = elem_nodes[2].getZ();
			const double& z4 = elem_nodes[3].getZ();
			const double& z5 = elem_nodes[4].getZ();
			const double& z6 = elem_nodes[5].getZ();
			const double& z7 = elem_nodes[6].getZ();
			const double& z8 = elem_nodes[7].getZ();
			const double& z9 = elem_nodes[8].getZ();
			const double& z10 = elem_nodes[9].getZ();

			Jacobian_[0][0] = x1*(4. * xi - 1.) + 4. * eta*x5 - 4. * eta*x9 + 4. * x7*zeta - 4. * x10*zeta + x4*(4. * eta + 4. * xi + 4. * zeta - 3.) - x8*(4. * eta + 8. * xi + 4. * zeta - 4.);
			Jacobian_[1][0] = y1*(4. * xi - 1.) + 4. * eta*y5 - 4. * eta*y9 + 4. * y7*zeta - 4. * y10*zeta + y4*(4. * eta + 4. * xi + 4. * zeta - 3.) - y8*(4. * eta + 8. * xi + 4. * zeta - 4.);
			Jacobian_[2][0] = z1*(4. * xi - 1.) + 4. * eta*z5 - 4. * eta*z9 + 4. * z7*zeta - 4. * z10*zeta + z4*(4. * eta + 4. * xi + 4. * zeta - 3.) - z8*(4. * eta + 8. * xi + 4. * zeta - 4.);
			Jacobian_[0][1] = 4. * x5*xi - 4. * x8*xi + 4. * x6*zeta - 4. * x10*zeta + x4*(4. * eta + 4. * xi + 4. * zeta - 3.) - x9*(8. * eta + 4. * xi + 4. * zeta - 4.) + x2*(4. * eta - 1.);
			Jacobian_[1][1] = 4. * xi*y5 - 4. * xi*y8 + 4. * y6*zeta - 4. * y10*zeta + y4*(4. * eta + 4. * xi + 4. * zeta - 3.) - y9*(8. * eta + 4. * xi + 4. * zeta - 4.) + y2*(4. * eta - 1.);
			Jacobian_[2][1] = 4. * xi*z5 - 4. * xi*z8 + 4. * z6*zeta - 4. * z10*zeta + z4*(4. * eta + 4. * xi + 4. * zeta - 3.) - z9*(8. * eta + 4. * xi + 4. * zeta - 4.) + z2*(4. * eta - 1.);
			Jacobian_[0][2] = x3*(4. * zeta - 1.) + 4. * eta*x6 - 4. * eta*x9 + 4. * x7*xi - 4. * x8*xi + x4*(4. * eta + 4. * xi + 4. * zeta - 3.) - x10*(4. * eta + 4. * xi + 8. * zeta - 4.);
			Jacobian_[1][2] = y3*(4. * zeta - 1.) + 4. * eta*y6 - 4. * eta*y9 + 4. * xi*y7 - 4. * xi*y8 + y4*(4. * eta + 4. * xi + 4. * zeta - 3.) - y10*(4. * eta + 4. * xi + 8. * zeta - 4.);
			Jacobian_[2][2] = z3*(4. * zeta - 1.) + 4. * eta*z6 - 4. * eta*z9 + 4. * xi*z7 - 4. * xi*z8 + z4*(4. * eta + 4. * xi + 4. * zeta - 3.) - z10*(4. * eta + 4. * xi + 8. * zeta - 4.);
			return Jacobian_;
		}

		virtual double
			evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_Jacobian(iso_point, elem_nodes);


			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] * Jacobian_[2][2] +
				Jacobian_[1][0] * Jacobian_[2][1] * Jacobian_[0][2] +
				Jacobian_[0][1] * Jacobian_[1][2] * Jacobian_[2][0] -
				Jacobian_[0][2] * Jacobian_[1][1] * Jacobian_[2][0] -
				Jacobian_[0][1] * Jacobian_[1][0] * Jacobian_[2][2] -
				Jacobian_[1][2] * Jacobian_[2][1] * Jacobian_[0][0];


			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_detJacobian(iso_point, elem_nodes);
			double& J11 = Jacobian_[0][0];
			double& J21 = Jacobian_[1][0];
			double& J31 = Jacobian_[2][0];
			double& J12 = Jacobian_[0][1];
			double& J22 = Jacobian_[1][1];
			double& J32 = Jacobian_[2][1];
			double& J13 = Jacobian_[0][2];
			double& J23 = Jacobian_[1][2];
			double& J33 = Jacobian_[2][2];

			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*(J22*J33 - J32*J23);
			inv_Jacobian_[1][0] = (1. / det_Jacobian_)*(J31*J23 - J21*J33);
			inv_Jacobian_[2][0] = (1. / det_Jacobian_)*(J21*J32 - J31*J22);
			inv_Jacobian_[0][1] = (1. / det_Jacobian_)*(J32*J13 - J12*J33);
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*(J11*J33 - J31*J13);
			inv_Jacobian_[2][1] = (1. / det_Jacobian_)*(J31*J12 - J11*J32);
			inv_Jacobian_[0][2] = (1. / det_Jacobian_)*(J12*J23 - J22*J13);
			inv_Jacobian_[1][2] = (1. / det_Jacobian_)*(J21*J13 - J11*J23);
			inv_Jacobian_[2][2] = (1. / det_Jacobian_)*(J11*J22 - J21*J12);

			return inv_Jacobian_;
		}

	private:

		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 3x10
		std::vector<std::vector<double>> dNdx_;  // 3x10
		std::vector<std::vector<double>> Jacobian_; // 3x3
		std::vector<std::vector<double>> inv_Jacobian_; // 3x3
		double det_Jacobian_;

		ShapeFunction() :
			N_(10),
			dNdxi_(3, std::vector<double>(10)),
			dNdx_(3, std::vector<double>(10)),
			Jacobian_(3, std::vector<double>(3)),
			inv_Jacobian_(3, std::vector<double>(3)){}

		ShapeFunction(const ShapeFunction&){}

		ShapeFunction& operator=(const ShapeFunction&){}
	};
	
	template<>
	class ShapeFunction< Dim3D, Hexa20, Lagrange > :
		public LagrangeType<Dim3D>
	{
	public:

		using PointType = Point<Dim3D, CartesianCoordinate>;

		friend class SingletonHolder < ShapeFunction< Dim3D, Hexa20, Lagrange > >;

		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			N_[0] = -0.125*(1. - xi)*(1. - eta)*(1. - zeta)*(2. + xi + eta + zeta);
			N_[1] = -0.125*(1. + xi)*(1. - eta)*(1. - zeta)*(2. - xi + eta + zeta);
			N_[2] = -0.125*(1. + xi)*(1. + eta)*(1. - zeta)*(2. - xi - eta + zeta);
			N_[3] = -0.125*(1. - xi)*(1. + eta)*(1. - zeta)*(2. + xi - eta + zeta);
			N_[4] = -0.125*(1. - xi)*(1. - eta)*(1. + zeta)*(2. + xi + eta - zeta);
			N_[5] = -0.125*(1. + xi)*(1. - eta)*(1. + zeta)*(2. - xi + eta - zeta);
			N_[6] = -0.125*(1. + xi)*(1. + eta)*(1. + zeta)*(2. - xi - eta - zeta);
			N_[7] = -0.125*(1. - xi)*(1. + eta)*(1. + zeta)*(2. + xi - eta - zeta);
			N_[8] = 0.25*(1. - xi * xi)*(1. - eta)*(1. - zeta);
			N_[9] = 0.25*(1. + xi)*(1. - eta * eta)*(1. - zeta);
			N_[10] = 0.25*(1. - xi * xi)*(1. + eta)*(1. - zeta);
			N_[11] = 0.25*(1. - xi)*(1. - eta * eta)*(1. - zeta);
			N_[12] = 0.25*(1. - xi * xi)*(1. - eta)*(1. + zeta);
			N_[13] = 0.25*(1. + xi)*(1. - eta * eta)*(1. + zeta);
			N_[14] = 0.25*(1. - xi * xi)*(1. + eta)*(1. + zeta);
			N_[15] = 0.25*(1. - xi)*(1. - eta * eta)*(1. + zeta);
			N_[16] = 0.25*(1. - xi)*(1. - eta)*(1. - zeta * zeta);
			N_[17] = 0.25*(1. + xi)*(1. - eta)*(1. - zeta * zeta);
			N_[18] = 0.25*(1. + xi)*(1. + eta)*(1. - zeta * zeta);
			N_[19] = 0.25*(1. - xi)*(1. + eta)*(1. - zeta * zeta);

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();
			
			dNdxi_[0][0] = 0.125*(eta - 1)*(zeta - 1)*(2 * xi + eta + zeta + 1);
			dNdxi_[0][1] = -0.125*(eta - 1)*(zeta - 1)*(-2 * xi + eta + zeta + 1);
			dNdxi_[0][2] = -0.125*(eta + 1)*(zeta - 1)*(2 * xi + eta - zeta - 1);
			dNdxi_[0][3] = 0.125*(eta + 1)*(zeta - 1)*(-2 * xi + eta - zeta - 1);
			dNdxi_[0][4] = -0.125*(eta - 1)*(zeta + 1)*(2 * xi + eta - zeta + 1);
			dNdxi_[0][5] = 0.125*(eta - 1)*(zeta + 1)*(-2 * xi + eta - zeta + 1);
			dNdxi_[0][6] = 0.125*(eta + 1)*(zeta + 1)*(2 * xi + eta + zeta - 1);
			dNdxi_[0][7] = -0.125*(eta + 1)*(zeta + 1)*(-2 * xi + eta + zeta - 1);
			dNdxi_[0][8] = -0.5*xi*(eta - 1)*(zeta - 1);
			dNdxi_[0][9] = 0.25*(eta*eta - 1)*(zeta - 1);
			dNdxi_[0][10] = 0.5*xi*(eta + 1)*(zeta - 1);
			dNdxi_[0][11] = -0.25*(eta*eta - 1)*(zeta - 1);
			dNdxi_[0][12] = 0.5*xi*(eta - 1)*(zeta + 1);
			dNdxi_[0][13] = -0.25*(eta*eta - 1)*(zeta + 1);
			dNdxi_[0][14] = -0.5*xi*(eta + 1)*(zeta + 1);
			dNdxi_[0][15] = 0.25*(eta*eta - 1)*(zeta + 1);
			dNdxi_[0][16] = -0.25*(eta - 1)*(zeta*zeta - 1);
			dNdxi_[0][17] = 0.25*(eta - 1)*(zeta*zeta - 1);
			dNdxi_[0][18] = -0.25*(eta + 1)*(zeta*zeta - 1);
			dNdxi_[0][19] = 0.25*(eta + 1)*(zeta*zeta - 1);

			dNdxi_[1][0] = 0.125*(xi - 1)*(zeta - 1)*(xi + 2 * eta + zeta + 1);
			dNdxi_[1][1] = 0.125*(xi + 1)*(zeta - 1)*(xi - 2 * eta - zeta - 1);
			dNdxi_[1][2] = -0.125*(xi + 1)*(zeta - 1)*(xi + 2 * eta - zeta - 1);
			dNdxi_[1][3] = -0.125*(xi - 1)*(zeta - 1)*(xi - 2 * eta + zeta - 1);
			dNdxi_[1][4] = -0.125*(xi - 1)*(zeta + 1)*(xi + 2 * eta - zeta + 1);
			dNdxi_[1][5] = -0.125*(xi + 1)*(zeta + 1)*(xi - 2 * eta + zeta - 1);
			dNdxi_[1][6] = 0.125*(xi + 1)*(zeta + 1)*(xi + 2 * eta + zeta - 1);
			dNdxi_[1][7] = 0.125*(xi - 1)*(zeta + 1)*(xi - 2 * eta - zeta + 1);
			dNdxi_[1][8] = -0.25*(xi*xi - 1)*(zeta - 1);
			dNdxi_[1][9] = 0.5*(xi + 1)*eta*(zeta - 1);
			dNdxi_[1][10] = 0.25*(xi*xi - 1)*(zeta - 1);
			dNdxi_[1][11] = -0.5*(xi - 1)*eta*(zeta - 1);
			dNdxi_[1][12] = 0.25*(xi*xi - 1)*(zeta + 1);
			dNdxi_[1][13] = -0.5*(xi + 1)*eta*(zeta + 1);
			dNdxi_[1][14] = -0.25*(xi*xi - 1)*(zeta + 1);
			dNdxi_[1][15] = 0.5*(xi - 1)*eta*(zeta + 1);
			dNdxi_[1][16] = -0.25*(xi - 1)*(zeta*zeta - 1);
			dNdxi_[1][17] = 0.25*(xi + 1)*(zeta*zeta - 1);
			dNdxi_[1][18] = -0.25*(xi + 1)*(zeta*zeta - 1);
			dNdxi_[1][19] = 0.25*(xi - 1)*(zeta*zeta - 1);
			
			dNdxi_[2][0] = 0.125*(xi - 1)*(eta - 1)*(xi + eta + 2 * zeta + 1);
			dNdxi_[2][1] = 0.125*(xi + 1)*(eta - 1)*(xi - eta - 2 * zeta - 1);
			dNdxi_[2][2] = -0.125*(xi + 1)*(eta + 1)*(xi + eta - 2 * zeta - 1);
			dNdxi_[2][3] = -0.125*(xi - 1)*(eta + 1)*(xi - eta + 2 * zeta + 1);
			dNdxi_[2][4] = -0.125*(xi - 1)*(eta - 1)*(xi + eta - 2 * zeta + 1);
			dNdxi_[2][5] = -0.125*(xi + 1)*(eta - 1)*(xi - eta + 2 * zeta - 1);
			dNdxi_[2][6] = 0.125*(xi + 1)*(eta + 1)*(xi + eta + 2 * zeta - 1);
			dNdxi_[2][7] = 0.125*(xi - 1)*(eta + 1)*(xi - eta - 2 * zeta + 1);
			dNdxi_[2][8] = -0.25*(xi*xi - 1)*(eta - 1);
			dNdxi_[2][9] = 0.25*(xi + 1)*(eta*eta - 1);
			dNdxi_[2][10] = 0.25*(xi*xi - 1)*(eta + 1);
			dNdxi_[2][11] = -0.25*(xi - 1)*(eta*eta - 1);
			dNdxi_[2][12] = 0.25*(xi*xi - 1)*(eta - 1);
			dNdxi_[2][13] = -0.25*(xi + 1)*(eta*eta - 1);
			dNdxi_[2][14] = -0.25*(xi*xi - 1)*(eta + 1);
			dNdxi_[2][15] = 0.25*(xi - 1)*(eta*eta - 1);
			dNdxi_[2][16] = -0.5*(xi - 1)*(eta - 1)*zeta;
			dNdxi_[2][17] = 0.5*(xi + 1)*(eta - 1)*zeta;
			dNdxi_[2][18] = -0.5*(xi + 1)*(eta + 1)*zeta;
			dNdxi_[2][19] = 0.5*(xi - 1)*(eta + 1)*zeta;

			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 10; ++j){
				for (int i = 0; i < 3; ++i){
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1] +
						dNdxi_[2][j] * inv_Jacobian_[i][2];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();
			const double& x5 = elem_nodes[4].getX();
			const double& x6 = elem_nodes[5].getX();
			const double& x7 = elem_nodes[6].getX();
			const double& x8 = elem_nodes[7].getX();
			const double& x9 = elem_nodes[8].getX();
			const double& x10 = elem_nodes[9].getX();
			const double& x11 = elem_nodes[10].getX();
			const double& x12 = elem_nodes[11].getX();
			const double& x13 = elem_nodes[12].getX();
			const double& x14 = elem_nodes[13].getX();
			const double& x15 = elem_nodes[14].getX();
			const double& x16 = elem_nodes[15].getX();
			const double& x17 = elem_nodes[16].getX();
			const double& x18 = elem_nodes[17].getX();
			const double& x19 = elem_nodes[18].getX();
			const double& x20 = elem_nodes[19].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();
			const double& y5 = elem_nodes[4].getY();
			const double& y6 = elem_nodes[5].getY();
			const double& y7 = elem_nodes[6].getY();
			const double& y8 = elem_nodes[7].getY();
			const double& y9 = elem_nodes[8].getY();
			const double& y10 = elem_nodes[9].getY();
			const double& y11 = elem_nodes[10].getY();
			const double& y12 = elem_nodes[11].getY();
			const double& y13 = elem_nodes[12].getY();
			const double& y14 = elem_nodes[13].getY();
			const double& y15 = elem_nodes[14].getY();
			const double& y16 = elem_nodes[15].getY();
			const double& y17 = elem_nodes[16].getY();
			const double& y18 = elem_nodes[17].getY();
			const double& y19 = elem_nodes[18].getY();
			const double& y20 = elem_nodes[19].getY();

			const double& z1 = elem_nodes[0].getZ();
			const double& z2 = elem_nodes[1].getZ();
			const double& z3 = elem_nodes[2].getZ();
			const double& z4 = elem_nodes[3].getZ();
			const double& z5 = elem_nodes[4].getZ();
			const double& z6 = elem_nodes[5].getZ();
			const double& z7 = elem_nodes[6].getZ();
			const double& z8 = elem_nodes[7].getZ();
			const double& z9 = elem_nodes[8].getZ();
			const double& z10 = elem_nodes[9].getZ();
			const double& z11 = elem_nodes[10].getZ();
			const double& z12 = elem_nodes[11].getZ();
			const double& z13 = elem_nodes[12].getZ();
			const double& z14 = elem_nodes[13].getZ();
			const double& z15 = elem_nodes[14].getZ();
			const double& z16 = elem_nodes[15].getZ();
			const double& z17 = elem_nodes[16].getZ();
			const double& z18 = elem_nodes[17].getZ();
			const double& z19 = elem_nodes[18].getZ();
			const double& z20 = elem_nodes[19].getZ();

			Jacobian_[0][0] =
				  (x10 - x12)*(eta * eta / 4. - 1. / 4.)*(zeta - 1.)
				+ (x16 - x14)*(eta * eta / 4. - 1. / 4.)*(zeta + 1.)
				+ (x18 - x17)*(eta / 4. - 1. / 4.)*(zeta * zeta - 1.)
				+ (x20 - x19)*(eta / 4. + 1. / 4.)*(zeta * zeta - 1.)
				- x3*(eta / 8. + 1. / 8.)*(zeta - 1.)*(eta + 2. * xi - zeta - 1.)
				+ x4*(eta / 8. + 1. / 8.)*(zeta - 1.)*(eta - 2. * xi - zeta - 1.)
				- x5*(eta / 8. - 1. / 8.)*(zeta + 1.)*(eta + 2. * xi - zeta + 1.)
				+ x6*(eta / 8. - 1. / 8.)*(zeta + 1.)*(eta - 2. * xi - zeta + 1.)
				+ x1*(eta / 8. - 1. / 8.)*(zeta - 1.)*(eta + 2. * xi + zeta + 1.)
				- x2*(eta / 8. - 1. / 8.)*(zeta - 1.)*(eta - 2. * xi + zeta + 1.)
				+ x7*(eta / 8. + 1. / 8.)*(zeta + 1.)*(eta + 2. * xi + zeta - 1.)
				- x8*(eta / 8. + 1. / 8.)*(zeta + 1.)*(eta - 2. * xi + zeta - 1.)
				- (x9*xi*(eta - 1.)*(zeta - 1.)) / 2.
				+ (x11*xi*(eta + 1.)*(zeta - 1.)) / 2.
				+ (x13*xi*(eta - 1.)*(zeta + 1.)) / 2.
				- (x15*xi*(eta + 1.)*(zeta + 1.)) / 2.;

			Jacobian_[1][0] = 
				  (y10 - y12)*(eta * eta / 4. - 1. / 4.)*(zeta - 1.)
				+ (y16 - y14)*(eta * eta / 4. - 1. / 4.)*(zeta + 1.)
				+ (y18 - y17)*(eta / 4. - 1. / 4.)*(zeta * zeta - 1.)
				+ (y20 - y19)*(eta / 4. + 1. / 4.)*(zeta * zeta - 1.)
				- y3*(eta / 8. + 1. / 8.)*(zeta - 1.)*(eta + 2. * xi - zeta - 1.)
				+ y4*(eta / 8. + 1. / 8.)*(zeta - 1.)*(eta - 2. * xi - zeta - 1.)
				- y5*(eta / 8. - 1. / 8.)*(zeta + 1.)*(eta + 2. * xi - zeta + 1.)
				+ y6*(eta / 8. - 1. / 8.)*(zeta + 1.)*(eta - 2. * xi - zeta + 1.)
				+ y1*(eta / 8. - 1. / 8.)*(zeta - 1.)*(eta + 2. * xi + zeta + 1.)
				- y2*(eta / 8. - 1. / 8.)*(zeta - 1.)*(eta - 2. * xi + zeta + 1.)
				+ y7*(eta / 8. + 1. / 8.)*(zeta + 1.)*(eta + 2. * xi + zeta - 1.)
				- y8*(eta / 8. + 1. / 8.)*(zeta + 1.)*(eta - 2. * xi + zeta - 1.)
				- (xi*y9*(eta - 1.)*(zeta - 1.)) / 2.
				+ (xi*y11*(eta + 1.)*(zeta - 1.)) / 2.
				+ (xi*y13*(eta - 1.)*(zeta + 1.)) / 2.
				- (xi*y15*(eta + 1.)*(zeta + 1.)) / 2.;

			Jacobian_[2][0] = 
			   	  (z10 - z12)*(eta * eta / 4. - 1. / 4.)*(zeta - 1.)
				+ (z16 - z14)*(eta * eta / 4. - 1. / 4.)*(zeta + 1.)
				+ (z18 - z17)*(eta / 4. - 1. / 4.)*(zeta * zeta - 1.)
				+ (z20 - z19)*(eta / 4. + 1. / 4.)*(zeta * zeta - 1.)
				- z3*(eta / 8. + 1. / 8.)*(zeta - 1.)*(eta + 2. * xi - zeta - 1.) 
				- z4*(eta / 8. + 1. / 8.)*(zeta - 1.)*(2. * xi - eta + zeta + 1.)
				- z5*(eta / 8. - 1. / 8.)*(zeta + 1.)*(eta + 2. * xi - zeta + 1.)
				+ z6*(eta / 8. - 1. / 8.)*(zeta + 1.)*(eta - 2. * xi - zeta + 1.)
				+ z1*(eta / 8. - 1. / 8.)*(zeta - 1.)*(eta + 2. * xi + zeta + 1.)
				- z2*(eta / 8. - 1. / 8.)*(zeta - 1.)*(eta - 2. * xi + zeta + 1.)
				+ z7*(eta / 8. + 1. / 8.)*(zeta + 1.)*(eta + 2. * xi + zeta - 1.)
				- z8*(eta / 8. + 1. / 8.)*(zeta + 1.)*(eta - 2. * xi + zeta - 1.)
				- (xi*z9*(eta - 1.)*(zeta - 1.)) / 2.
				+ (xi*z11*(eta + 1.)*(zeta - 1.)) / 2.
				+ (xi*z13*(eta - 1.)*(zeta + 1.)) / 2.
				- (xi*z15*(eta + 1.)*(zeta + 1.)) / 2.;

			Jacobian_[0][1] = 
				  (x11 - x9)*(xi * xi / 4. - 1. / 4.)*(zeta - 1.)
				+ (x13 - x15)*(xi * xi / 4. - 1. / 4.)*(zeta + 1.)
				+ (x18 - x19)*(xi / 4. + 1. / 4.)*(zeta * zeta - 1.)
				+ (x20 - x17)*(xi / 4. - 1. / 4.)*(zeta * zeta - 1.)
				+ eta*x10*(xi / 2. + 1. / 2.)*(zeta - 1.)
				- eta*x12*(xi / 2. - 1. / 2.)*(zeta - 1.)
				- eta*x14*(xi / 2. + 1. / 2.)*(zeta + 1.)
				+ eta*x16*(xi / 2. - 1. / 2.)*(zeta + 1.)
				- x2*(xi / 8. + 1. / 8.)*(zeta - 1.)*(2. * eta - xi + zeta + 1.)
				- x3*(xi / 8. + 1. / 8.)*(zeta - 1.)*(2. * eta + xi - zeta - 1.)
				- x5*(xi / 8. - 1. / 8.)*(zeta + 1.)*(2. * eta + xi - zeta + 1.)
				- x8*(xi / 8. - 1. / 8.)*(zeta + 1.)*(2. * eta - xi + zeta - 1.)
				+ x4*(xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta - xi - zeta + 1.)
				+ x6*(xi / 8. + 1. / 8.)*(zeta + 1.)*(2. * eta - xi - zeta + 1.)
				+ x1*(xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta + xi + zeta + 1.)
				+ x7*(xi / 8. + 1. / 8.)*(zeta + 1.)*(2. * eta + xi + zeta - 1.);

			Jacobian_[1][1] = 
				  (y11 - y9)*(xi * xi / 4. - 1. / 4.)*(zeta - 1.)
				+ (y13 - y15)*(xi * xi / 4. - 1. / 4.)*(zeta + 1.)
				+ (y18 - y19)*(xi / 4. + 1. / 4.)*(zeta * zeta - 1.)
				+ (y20 - y17)*(xi / 4. - 1. / 4.)*(zeta * zeta - 1.)
				+ y1*(xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta + xi + zeta + 1.)
				+ y7*(xi / 8. + 1. / 8.)*(zeta + 1.)*(2. * eta + xi + zeta - 1.)
				+ eta*y10*(xi / 2. + 1. / 2.)*(zeta - 1.)
				- eta*y12*(xi / 2. - 1. / 2.)*(zeta - 1.)
				- eta*y14*(xi / 2. + 1. / 2.)*(zeta + 1.)
				+ eta*y16*(xi / 2. - 1. / 2.)*(zeta + 1.)
				- y2*(xi / 8. + 1. / 8.)*(zeta - 1.)*(2. * eta - xi + zeta + 1.)
				- y3*(xi / 8. + 1. / 8.)*(zeta - 1.)*(2. * eta + xi - zeta - 1.)
				- y5*(xi / 8. - 1. / 8.)*(zeta + 1.)*(2. * eta + xi - zeta + 1.)
				- y8*(xi / 8. - 1. / 8.)*(zeta + 1.)*(2. * eta - xi + zeta - 1.)
				+ y4*(xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta - xi - zeta + 1.)
				+ y6*(xi / 8. + 1. / 8.)*(zeta + 1.)*(2. * eta - xi - zeta + 1.);

			Jacobian_[2][1] = 
				  (z11 - z9)*(xi * xi / 4. - 1. / 4.)*(zeta - 1.)
				+ (z13 - z15)*(xi * xi / 4. - 1. / 4.)*(zeta + 1.)
				+ (z18 - z19)*(xi / 4. + 1. / 4.)*(zeta * zeta - 1.)
				+ (z20 - z17)*(xi / 4. - 1. / 4.)*(zeta * zeta - 1.)
				+ z1*(xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta + xi + zeta + 1.)
				+ z7*(xi / 8. + 1. / 8.)*(zeta + 1.)*(2. * eta + xi + zeta - 1.)
				+ eta*z10*(xi / 2. + 1. / 2.)*(zeta - 1.)
				- eta*z12*(xi / 2. - 1. / 2.)*(zeta - 1.)
				- eta*z14*(xi / 2. + 1. / 2.)*(zeta + 1.)
				+ eta*z16*(xi / 2. - 1. / 2.)*(zeta + 1.)
				- z2*(xi / 8. + 1. / 8.)*(zeta - 1.)*(2. * eta - xi + zeta + 1.)
				- z3*(xi / 8. + 1. / 8.)*(zeta - 1.)*(2. * eta + xi - zeta - 1.)
				- z5*(xi / 8. - 1. / 8.)*(zeta + 1.)*(2. * eta + xi - zeta + 1.)
				- z8*(xi / 8. - 1. / 8.)*(zeta + 1.)*(2. * eta - xi + zeta - 1.)
				+ z4*(xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta - xi - zeta + 1.)
				+ z6*(xi / 8. + 1. / 8.)*(zeta + 1.)*(2. * eta - xi - zeta + 1.);

			Jacobian_[0][2] = 
				  (x10-x14)*(eta * eta - 1.)*(xi / 4. + 1. / 4.)
				+ (x13-x9)*(xi * xi / 4. - 1. / 4.)*(eta - 1.)
				+ (x11-x15)*(xi * xi / 4. - 1. / 4.)*(eta + 1.)
				+ (x16-x12)*(eta * eta - 1.)*(xi / 4. - 1. / 4.)
				- x2*(xi / 8. + 1. / 8.)*(eta - 1.)*(eta - xi + 2. * zeta + 1.)
				- x4*(xi / 8. - 1. / 8.)*(eta + 1.)*(xi - eta + 2. * zeta + 1.)
				+ x6*(xi / 8. + 1. / 8.)*(eta - 1.)*(eta - xi - 2. * zeta + 1.)
				- x8*(xi / 8. - 1. / 8.)*(eta + 1.)*(eta - xi + 2. * zeta - 1.)
				- x17*zeta*(xi / 2. - 1. / 2.)*(eta - 1.) + x18*zeta*(xi / 2. + 1. / 2.)*(eta - 1.)
				- x19*zeta*(xi / 2. + 1. / 2.)*(eta + 1.) + x20*zeta*(xi / 2. - 1. / 2.)*(eta + 1.)
				+ x1*(xi / 8. - 1. / 8.)*(eta - 1.)*(eta + xi + 2. * zeta + 1.)
				- x3*(xi / 8. + 1. / 8.)*(eta + 1.)*(eta + xi - 2. * zeta - 1.)
				- x5*(xi / 8. - 1. / 8.)*(eta - 1.)*(eta + xi - 2. * zeta + 1.)
				+ x7*(xi / 8. + 1. / 8.)*(eta + 1.)*(eta + xi + 2. * zeta - 1.);

			Jacobian_[1][2] = 
				(y10 - y14)*(eta * eta - 1.)*(xi / 4. + 1. / 4.)
				+ (y13 - y9)*(xi * xi / 4. - 1. / 4.)*(eta - 1.)
				+ (y11-y15)*(xi * xi / 4. - 1. / 4.)*(eta + 1.)
				+ (y16-y12)*(eta * eta - 1.)*(xi / 4. - 1. / 4.)
				- y2*(xi / 8. + 1. / 8.)*(eta - 1.)*(eta - xi + 2. * zeta + 1.)
				- y4*(xi / 8. - 1. / 8.)*(eta + 1.)*(xi - eta + 2. * zeta + 1.)
				+ y6*(xi / 8. + 1. / 8.)*(eta - 1.)*(eta - xi - 2. * zeta + 1.)
				- y8*(xi / 8. - 1. / 8.)*(eta + 1.)*(eta - xi + 2. * zeta - 1.)
				- y17*zeta*(xi / 2. - 1. / 2.)*(eta - 1.)
				+ y18*zeta*(xi / 2. + 1. / 2.)*(eta - 1.)
				- y19*zeta*(xi / 2. + 1. / 2.)*(eta + 1.)
				+ y20*zeta*(xi / 2. - 1. / 2.)*(eta + 1.)
				+ y1*(xi / 8. - 1. / 8.)*(eta - 1.)*(eta + xi + 2. * zeta + 1.)
				- y3*(xi / 8. + 1. / 8.)*(eta + 1.)*(eta + xi - 2. * zeta - 1.)
				- y5*(xi / 8. - 1. / 8.)*(eta - 1.)*(eta + xi - 2. * zeta + 1.)
				+ y7*(xi / 8. + 1. / 8.)*(eta + 1.)*(eta + xi + 2. * zeta - 1.);

			Jacobian_[2][2] = 
				  (z10-z14)*(eta * eta - 1.)*(xi / 4. + 1. / 4.)
				+ (z13-z9)*(xi * xi / 4. - 1. / 4.)*(eta - 1.)
				+ (z11-z15)*(xi * xi / 4. - 1. / 4.)*(eta + 1.)			
				+ (z16-z12)*(eta * eta - 1.)*(xi / 4. - 1. / 4.)
				- z2*(xi / 8. + 1. / 8.)*(eta - 1.)*(eta - xi + 2. * zeta + 1.)
				- z4*(xi / 8. - 1. / 8.)*(eta + 1.)*(xi - eta + 2. * zeta + 1.)
				+ z6*(xi / 8. + 1. / 8.)*(eta - 1.)*(eta - xi - 2. * zeta + 1.)
				- z8*(xi / 8. - 1. / 8.)*(eta + 1.)*(eta - xi + 2. * zeta - 1.)
				- z17*zeta*(xi / 2. - 1. / 2.)*(eta - 1.)
				+ z18*zeta*(xi / 2. + 1. / 2.)*(eta - 1.)
				- z19*zeta*(xi / 2. + 1. / 2.)*(eta + 1.)
				+ z20*zeta*(xi / 2. - 1. / 2.)*(eta + 1.)
				+ z1*(xi / 8. - 1. / 8.)*(eta - 1.)*(eta + xi + 2. * zeta + 1.)
				- z3*(xi / 8. + 1. / 8.)*(eta + 1.)*(eta + xi - 2. * zeta - 1.)
				- z5*(xi / 8. - 1. / 8.)*(eta - 1.)*(eta + xi - 2. * zeta + 1.)
				+ z7*(xi / 8. + 1. / 8.)*(eta + 1.)*(eta + xi + 2. * zeta - 1.);
			return Jacobian_;
		}

		virtual double
			evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_Jacobian(iso_point, elem_nodes);


			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] * Jacobian_[2][2] +
				Jacobian_[1][0] * Jacobian_[2][1] * Jacobian_[0][2] +
				Jacobian_[0][1] * Jacobian_[1][2] * Jacobian_[2][0] -
				Jacobian_[0][2] * Jacobian_[1][1] * Jacobian_[2][0] -
				Jacobian_[0][1] * Jacobian_[1][0] * Jacobian_[2][2] -
				Jacobian_[1][2] * Jacobian_[2][1] * Jacobian_[0][0];


			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_detJacobian(iso_point, elem_nodes);
			double& J11 = Jacobian_[0][0];
			double& J21 = Jacobian_[1][0];
			double& J31 = Jacobian_[2][0];
			double& J12 = Jacobian_[0][1];
			double& J22 = Jacobian_[1][1];
			double& J32 = Jacobian_[2][1];
			double& J13 = Jacobian_[0][2];
			double& J23 = Jacobian_[1][2];
			double& J33 = Jacobian_[2][2];

			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*(J22*J33 - J32*J23);
			inv_Jacobian_[1][0] = (1. / det_Jacobian_)*(J31*J23 - J21*J33);
			inv_Jacobian_[2][0] = (1. / det_Jacobian_)*(J21*J32 - J31*J22);
			inv_Jacobian_[0][1] = (1. / det_Jacobian_)*(J32*J13 - J12*J33);
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*(J11*J33 - J31*J13);
			inv_Jacobian_[2][1] = (1. / det_Jacobian_)*(J31*J12 - J11*J32);
			inv_Jacobian_[0][2] = (1. / det_Jacobian_)*(J12*J23 - J22*J13);
			inv_Jacobian_[1][2] = (1. / det_Jacobian_)*(J21*J13 - J11*J23);
			inv_Jacobian_[2][2] = (1. / det_Jacobian_)*(J11*J22 - J21*J12);

			return inv_Jacobian_;
		}

	private:

		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 3x20
		std::vector<std::vector<double>> dNdx_;  // 3x20
		std::vector<std::vector<double>> Jacobian_; // 3x3
		std::vector<std::vector<double>> inv_Jacobian_; // 3x3
		double det_Jacobian_;

		ShapeFunction() :
			N_(20),
			dNdxi_(3, std::vector<double>(20)),
			dNdx_(3, std::vector<double>(20)),
			Jacobian_(3, std::vector<double>(3)),
			inv_Jacobian_(3, std::vector<double>(3)){}

		ShapeFunction(const ShapeFunction&){}

		ShapeFunction& operator=(const ShapeFunction&){}
	};

	template<>
	class ShapeFunction< Dim3D, Pyramid13, Lagrange > :
		public LagrangeType<Dim3D>
	{
	public:

		using PointType = Point<Dim3D, CartesianCoordinate>;

		friend class SingletonHolder < ShapeFunction< Dim3D, Pyramid13, Lagrange > >;

		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			N_[0] = -(1. / 16.)*(1. - xi)*(1. - eta)*(1. - zeta)*(4. + 3. * xi + 3. * eta + 2. * xi*eta + 2. * zeta + xi * zeta + eta * zeta + 2. * xi*eta*zeta);
			N_[1] = - (1. / 16.)*(1. + xi)*(1. - eta)*(1. - zeta)*(4. - 3. * xi + 3. * eta - 2. * xi*eta + 2. * zeta - xi * zeta + eta * zeta - 2. * xi*eta*zeta);
			N_[2] = - (1. / 16.)*(1. + xi)*(1. + eta)*(1. - zeta)*(4. - 3. * xi - 3. * eta + 2. * xi*eta + 2. * zeta - xi * zeta - eta * zeta + 2. * xi*eta*zeta);
			N_[3] = - (1. / 16.)*(1. - xi)*(1. + eta)*(1. - zeta)*(4. + 3. * xi - 3. * eta - 2. * xi*eta + 2. * zeta + xi * zeta - eta * zeta - 2. * xi*eta*zeta);
			N_[4] = (1. / 2.)*zeta*(1. + zeta);
			N_[5] = (1. / 8.)*(1. - xi * xi)*(1. - eta)*(1. - zeta)*(2. + eta + eta * zeta);
			N_[6] = (1. / 8.)*(1. + xi)*(1. - eta * eta)*(1. - zeta)*(2. - xi - xi * zeta);
			N_[7] = (1. / 8.)*(1. - xi * xi)*(1. + eta)*(1. - zeta)*(2. - eta - eta * zeta);
			N_[8] = (1. / 8.)*(1. - xi)*(1. - eta * eta)*(1. - zeta)*(2. + xi + xi * zeta);
			N_[9] = (1. / 4.)*(1. - xi)*(1. - eta)*(1. - zeta * zeta);
			N_[10] = (1. / 4.)*(1. + xi)*(1. - eta)*(1. - zeta * zeta);
			N_[11] = (1. / 4.)*(1. + xi)*(1. + eta)*(1. - zeta * zeta);
			N_[12] = (1. / 4.)*(1. - xi)*(1. + eta)*(1. - zeta * zeta);

			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			dNdxi_[0][0] = (1. / 16)*(1. - eta)*(1. - zeta)*(1. + 6. * xi + eta + 4. * xi*eta + zeta + 2. * xi*zeta - eta * zeta + 4. * xi*eta*zeta);
			dNdxi_[0][1] = -(1. / 16)*(1. - eta)*(1. - zeta)*(1. - 6. * xi + eta - 4. * xi*eta + zeta - 2. * xi*zeta - eta * zeta - 4. * xi*eta*zeta);
			dNdxi_[0][2] = -(1. / 16)*(1. + eta)*(1. - zeta)*(1. - 6. * xi - eta + 4. * xi*eta + zeta - 2. * xi*zeta + eta * zeta + 4. * xi*eta*zeta);
			dNdxi_[0][3] = (1. / 16)*(1. + eta)*(1. - zeta)*(1. + 6. * xi - eta - 4. * xi*eta + zeta + 2. * xi*zeta + eta * zeta - 4. * xi*eta*zeta);
			dNdxi_[0][4] = 0;
			dNdxi_[0][5] = -(1. / 4.)*xi*(1. - eta)*(1. - zeta)*(2. + eta + eta * zeta);
			dNdxi_[0][6] = (1. / 8.)*(1. - eta * eta)*(1. - zeta)*(1. - 2. * xi - zeta - 2. * xi*zeta);
			dNdxi_[0][7] = -(1. / 4.)*xi*(1. + eta)*(1. - zeta)*(2. - eta - eta * zeta);
			dNdxi_[0][8] = -(1. / 8.)*(1. - eta * eta)*(1. - zeta)*(1. + 2. * xi - zeta + 2. * xi*zeta);
			dNdxi_[0][9] = -(1. / 4.)*(1. - eta)*(1. - zeta * zeta);
			dNdxi_[0][10] = (1. / 4.)*(1. - eta)*(1. - zeta * zeta);
			dNdxi_[0][11] = (1. / 4.)*(1. + eta)*(1. - zeta * zeta);
			dNdxi_[0][12] = -(1. / 4.)*(1. + eta)*(1. - zeta * zeta);

			dNdxi_[1][0] = (1. / 16)*(1. - xi)*(1. - zeta)*(1. + xi + 6. * eta + 4. * xi*eta + zeta - xi * zeta + 2. * eta*zeta + 4. * xi*eta*zeta);
			dNdxi_[1][1] = (1. / 16)*(1. + xi)*(1. - zeta)*(1. - xi + 6. * eta - 4. * xi*eta + zeta + xi * zeta + 2. * eta*zeta - 4. * xi*eta*zeta);
			dNdxi_[1][2] = -(1. / 16)*(1. + xi)*(1. - zeta)*(1. - xi - 6. * eta + 4. * xi*eta + zeta + xi * zeta - 2. * eta*zeta + 4. * xi*eta*zeta);
			dNdxi_[1][3] = -(1. / 16)*(1. - xi)*(1. - zeta)*(1. + xi - 6. * eta - 4. * xi*eta + zeta - xi * zeta - 2. * eta*zeta - 4. * xi*eta*zeta);
			dNdxi_[1][4] = 0;
			dNdxi_[1][5] = -(1. / 8.)* (1. - xi * xi)*(1. - zeta)*(1. + 2. * eta - zeta + 2. * eta*zeta);
			dNdxi_[1][6] = -(1. / 4.)*(1. + xi)*eta*(1. - zeta)*(2. - xi - xi * zeta);
			dNdxi_[1][7] = (1. / 8.)*(1. - xi * xi)*(1. - zeta)*(1. - 2. * eta - zeta - 2. * eta*zeta);
			dNdxi_[1][8] = -(1. / 4.)*(1. - xi)*eta*(1. - zeta)*(2. + xi + xi * zeta);
			dNdxi_[1][9] = -(1. / 4.)*(1. - xi)*(1. - zeta * zeta);
			dNdxi_[1][10] = -(1. / 4.)*(1. + xi)*(1. - zeta * zeta);
			dNdxi_[1][11] = (1. / 4.)*(1. + xi)*(1. - zeta * zeta);
			dNdxi_[1][12] = (1. / 4.)*(1. - xi)*(1. - zeta * zeta);

			dNdxi_[2][0] = (1. / 8.)*(1. - xi)*(1. - eta)*(1. + xi + eta + 2. * zeta + xi * zeta + eta * zeta + 2. * xi*eta*zeta);
			dNdxi_[2][1] = (1. / 8.)*(1. + xi)*(1. - eta)*(1. - xi + eta + 2. * zeta - xi * zeta + eta * zeta - 2. * xi*eta*zeta);
			dNdxi_[2][2] = (1. / 8.)*(1. + xi)*(1. + eta)*(1. - xi - eta + 2. * zeta - xi * zeta - eta * zeta + 2. * xi*eta*zeta);
			dNdxi_[2][3] = (1. / 8.)*(1. - xi)*(1. + eta)*(1. + xi - eta + 2. * zeta + xi * zeta - eta * zeta - 2. * xi*eta*zeta);
			dNdxi_[2][4] = 0.5 + zeta;
			dNdxi_[2][5] = -(1. / 4.)*(1. - xi * xi)*(1. - eta)*(1. + eta * zeta);
			dNdxi_[2][6] = -(1. / 4.)*(1. + xi)*(1. - eta * eta)*(1. - xi * zeta);
			dNdxi_[2][7] = -(1. / 4.)*(1. - xi * xi)*(1. + eta)*(1. - eta * zeta);
			dNdxi_[2][8] = -(1. / 4.)*(1. - xi)*(1. - eta * eta)*(1. + xi * zeta);
			dNdxi_[2][9] = -(1. / 2.)*(1. - xi)*(1. - eta)*zeta;
			dNdxi_[2][10] = -(1. / 2.)*(1. + xi)*(1. - eta)*zeta;
			dNdxi_[2][11] = -(1. / 2.)*(1. + xi)*(1. + eta)*zeta;
			dNdxi_[2][12] = -(1. / 2.)*(1. - xi)*(1. + eta)*zeta;
	

			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 10; ++j) {
				for (int i = 0; i < 3; ++i) {
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1] +
						dNdxi_[2][j] * inv_Jacobian_[i][2];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();
			const double& x5 = elem_nodes[4].getX();
			const double& x6 = elem_nodes[5].getX();
			const double& x7 = elem_nodes[6].getX();
			const double& x8 = elem_nodes[7].getX();
			const double& x9 = elem_nodes[8].getX();
			const double& x10 = elem_nodes[9].getX();
			const double& x11 = elem_nodes[10].getX();
			const double& x12 = elem_nodes[11].getX();
			const double& x13 = elem_nodes[12].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();
			const double& y5 = elem_nodes[4].getY();
			const double& y6 = elem_nodes[5].getY();
			const double& y7 = elem_nodes[6].getY();
			const double& y8 = elem_nodes[7].getY();
			const double& y9 = elem_nodes[8].getY();
			const double& y10 = elem_nodes[9].getY();
			const double& y11 = elem_nodes[10].getY();
			const double& y12 = elem_nodes[11].getY();
			const double& y13 = elem_nodes[12].getY();

			const double& z1 = elem_nodes[0].getZ();
			const double& z2 = elem_nodes[1].getZ();
			const double& z3 = elem_nodes[2].getZ();
			const double& z4 = elem_nodes[3].getZ();
			const double& z5 = elem_nodes[4].getZ();
			const double& z6 = elem_nodes[5].getZ();
			const double& z7 = elem_nodes[6].getZ();
			const double& z8 = elem_nodes[7].getZ();
			const double& z9 = elem_nodes[8].getZ();
			const double& z10 = elem_nodes[9].getZ();
			const double& z11 = elem_nodes[10].getZ();
			const double& z12 = elem_nodes[11].getZ();
			const double& z13 = elem_nodes[12].getZ();

			Jacobian_[0][0] = 
				 (x11 - x10) * (eta / 4. - 1. / 4.)*(zeta * zeta - 1.)
				+(x13 - x12) * (eta / 4. + 1. / 4.)*(zeta * zeta - 1.)
				+ x2 * (eta / 16. - 1. / 16.)*(zeta - 1.)*(6 * xi - eta - zeta + 4. * eta*xi + eta * zeta + 2. * xi*zeta + 4. * eta*xi*zeta - 1.)
				- x7 * (eta * eta / 8. - 1. / 8.)*(zeta - 1.)*(2. * xi + zeta + 2. * xi*zeta - 1.)
				+ x1 * (eta / 16. - 1. / 16.)*(zeta - 1.)*(eta + 6 * xi + zeta + 4. * eta*xi - eta * zeta + 2. * xi*zeta + 4. * eta*xi*zeta + 1.)
				- x9 * (eta * eta / 8. - 1. / 8.)*(zeta - 1.)*(2. * xi - zeta + 2. * xi*zeta + 1.)
				+ x3 * (eta / 16. + 1. / 16.)*(zeta - 1.)*(zeta - 6 * xi - eta + 4. * eta*xi + eta * zeta - 2. * xi*zeta + 4. * eta*xi*zeta + 1.)
				- x4 * (eta / 16. + 1. / 16.)*(zeta - 1.)*(6 * xi - eta + zeta - 4. * eta*xi + eta * zeta + 2. * xi*zeta - 4. * eta*xi*zeta + 1.)
				- (x6*xi*(eta - 1.)*(zeta - 1.)*(eta + eta * zeta + 2.)) / 4.
				- (x8*xi*(eta + 1.)*(zeta - 1.)*(eta + eta * zeta - 2.)) / 4.;

			Jacobian_[1][0] = 
				  (y11- y10) * (eta / 4. - 1. / 4.)*(zeta * zeta - 1.)
				+ (y13- y12) * (eta / 4. + 1. / 4.)*(zeta * zeta - 1.)
				+ y2 * (eta / 16. - 1. / 16.)*(zeta - 1.)*(6 * xi - eta - zeta + 4. * eta*xi + eta * zeta + 2. * xi*zeta + 4. * eta*xi*zeta - 1.)
				- y7 * (eta * eta / 8. - 1. / 8.)*(zeta - 1.)*(2. * xi + zeta + 2. * xi*zeta - 1.) 
				+ y1 * (eta / 16. - 1. / 16.)*(zeta - 1.)*(eta + 6 * xi + zeta + 4. * eta*xi - eta * zeta + 2. * xi*zeta + 4. * eta*xi*zeta + 1.)
				- y9 * (eta * eta / 8. - 1. / 8.)*(zeta - 1.)*(2. * xi - zeta + 2. * xi*zeta + 1.)
				+ y3 * (eta / 16. + 1. / 16.)*(zeta - 1.)*(zeta - 6 * xi - eta + 4. * eta*xi + eta * zeta - 2. * xi*zeta + 4. * eta*xi*zeta + 1.)
				- y4 * (eta / 16. + 1. / 16.)*(zeta - 1.)*(6 * xi - eta + zeta - 4. * eta*xi + eta * zeta + 2. * xi*zeta - 4. * eta*xi*zeta + 1.)
				- (xi*y6*(eta - 1.)*(zeta - 1.)*(eta + eta * zeta + 2.)) / 4.
				- (xi*y8*(eta + 1.)*(zeta - 1.)*(eta + eta * zeta - 2.)) / 4.;

			Jacobian_[2][0] = 
				  (z11- z10) * (eta / 4. - 1. / 4.)*(zeta * zeta - 1.)
				+ (z13- z12) * (eta / 4. + 1. / 4.)*(zeta * zeta - 1.) 
				+ z2 * (eta / 16. - 1. / 16.)*(zeta - 1.)*(6 * xi - eta - zeta + 4. * eta*xi + eta * zeta + 2. * xi*zeta + 4. * eta*xi*zeta - 1.)
				- z7 * (eta * eta / 8. - 1. / 8.)*(zeta - 1.)*(2. * xi + zeta + 2. * xi*zeta - 1.)
				+ z1 * (eta / 16. - 1. / 16.)*(zeta - 1.)*(eta + 6 * xi + zeta + 4. * eta*xi - eta * zeta + 2. * xi*zeta + 4. * eta*xi*zeta + 1.)
				- z9 * (eta * eta / 8. - 1. / 8.)*(zeta - 1.)*(2. * xi - zeta + 2. * xi*zeta + 1.)
				+ z3 * (eta / 16. + 1. / 16.)*(zeta - 1.)*(zeta - 6 * xi - eta + 4. * eta*xi + eta * zeta - 2. * xi*zeta + 4. * eta*xi*zeta + 1.)
				- z4 * (eta / 16. + 1. / 16.)*(zeta - 1.)*(6 * xi - eta + zeta - 4. * eta*xi + eta * zeta + 2. * xi*zeta - 4. * eta*xi*zeta + 1.)
				- (xi*z6*(eta - 1.)*(zeta - 1.)*(eta + eta * zeta + 2.)) / 4.
				- (xi*z8*(eta + 1.)*(zeta - 1.)*(eta + eta * zeta - 2.)) / 4.;

			Jacobian_[0][1] = 
				  (x11- x12) * (xi / 4. + 1. / 4.)*(zeta * zeta - 1.)
				+ (x13- x10) * (xi / 4. - 1. / 4.)*(zeta * zeta - 1.)
				- x8 * (xi * xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta + zeta + 2. * eta*zeta - 1.)
				+ x4 * (xi / 16. - 1. / 16.)*(zeta - 1.)*(6 * eta - xi - zeta + 4. * eta*xi + 2. * eta*zeta + xi * zeta + 4. * eta*xi*zeta - 1.)
				- x6 * (xi * xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta - zeta + 2. * eta*zeta + 1.)
				+ x1 * (xi / 16. - 1. / 16.)*(zeta - 1.)*(6 * eta + xi + zeta + 4. * eta*xi + 2. * eta*zeta - xi * zeta + 4. * eta*xi*zeta + 1.)
				- x2 * (xi / 16. + 1. / 16.)*(zeta - 1.)*(6 * eta - xi + zeta - 4. * eta*xi + 2. * eta*zeta + xi * zeta - 4. * eta*xi*zeta + 1.)
				+ x3 * (xi / 16. + 1. / 16.)*(zeta - 1.)*(zeta - xi - 6 * eta + 4. * eta*xi - 2. * eta*zeta + xi * zeta + 4. * eta*xi*zeta + 1.)
				- eta * x7*(xi / 4. + 1. / 4.)*(zeta - 1.)*(xi + xi * zeta - 2.)
				- eta * x9*(xi / 4. - 1. / 4.)*(zeta - 1.)*(xi + xi * zeta + 2.);

			Jacobian_[1][1] = 
				  (y11- y12) * (xi / 4. + 1. / 4.)*(zeta * zeta - 1.)
				+ (y13- y10) * (xi / 4. - 1. / 4.)*(zeta * zeta - 1.)
				- y8 * (xi * xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta + zeta + 2. * eta*zeta - 1.)
				+ y4 * (xi / 16. - 1. / 16.)*(zeta - 1.)*(6 * eta - xi - zeta + 4. * eta*xi + 2. * eta*zeta + xi * zeta + 4. * eta*xi*zeta - 1.)
				- y6 * (xi * xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta - zeta + 2. * eta*zeta + 1.)
				+ y1 * (xi / 16. - 1. / 16.)*(zeta - 1.)*(6 * eta + xi + zeta + 4. * eta*xi + 2. * eta*zeta - xi * zeta + 4. * eta*xi*zeta + 1.)
				- y2 * (xi / 16. + 1. / 16.)*(zeta - 1.)*(6 * eta - xi + zeta - 4. * eta*xi + 2. * eta*zeta + xi * zeta - 4. * eta*xi*zeta + 1.)
				+ y3 * (xi / 16. + 1. / 16.)*(zeta - 1.)*(zeta - xi - 6 * eta + 4. * eta*xi - 2. * eta*zeta + xi * zeta + 4. * eta*xi*zeta + 1.)
				- eta * y7*(xi / 4. + 1. / 4.)*(zeta - 1.)*(xi + xi * zeta - 2.)
				- eta * y9*(xi / 4. - 1. / 4.)*(zeta - 1.)*(xi + xi * zeta + 2.);

			Jacobian_[2][1] = 
				  (z11- z12) * (xi / 4. + 1. / 4.)*(zeta * zeta - 1.)
				+ (z13- z10) * (xi / 4. - 1. / 4.)*(zeta * zeta - 1.)
				- z8 * (xi * xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta + zeta + 2. * eta*zeta - 1.)
				+ z4 * (xi / 16. - 1. / 16.)*(zeta - 1.)*(6 * eta - xi - zeta + 4. * eta*xi + 2. * eta*zeta + xi * zeta + 4. * eta*xi*zeta - 1.)
				- z6 * (xi * xi / 8. - 1. / 8.)*(zeta - 1.)*(2. * eta - zeta + 2. * eta*zeta + 1.)
				+ z1 * (xi / 16. - 1. / 16.)*(zeta - 1.)*(6 * eta + xi + zeta + 4. * eta*xi + 2. * eta*zeta - xi * zeta + 4. * eta*xi*zeta + 1.)
				- z2 * (xi / 16. + 1. / 16.)*(zeta - 1.)*(6 * eta - xi + zeta - 4. * eta*xi + 2. * eta*zeta + xi * zeta - 4. * eta*xi*zeta + 1.)
				+ z3 * (xi / 16. + 1. / 16.)*(zeta - 1.)*(zeta - xi - 6 * eta + 4. * eta*xi - 2. * eta*zeta + xi * zeta + 4. * eta*xi*zeta + 1.)
				- eta * z7*(xi / 4. + 1. / 4.)*(zeta - 1.)*(xi + xi * zeta - 2.)
				- eta * z9*(xi / 4. - 1. / 4.)*(zeta - 1.)*(xi + xi * zeta + 2.);

			Jacobian_[0][2] = 
				x5 * (zeta + 1. / 2.)
				+ x1 * (xi / 8. - 1. / 8.)*(eta - 1.)*(eta + xi + 2. * zeta + eta * zeta + xi * zeta + 2. * eta*xi*zeta + 1.)
				- x3 * (xi / 8. + 1. / 8.)*(eta + 1.)*(eta + xi - 2. * zeta + eta * zeta + xi * zeta - 2. * eta*xi*zeta - 1.)
				- x10 * zeta*(xi / 2. - 1. / 2.)*(eta - 1.)
				+ x11 * zeta*(xi / 2. + 1. / 2.)*(eta - 1.)
				- x12 * zeta*(xi / 2. + 1. / 2.)*(eta + 1.)
				+ x13 * zeta*(xi / 2. - 1. / 2.)*(eta + 1.)
				- x2 * (xi / 8. + 1. / 8.)*(eta - 1.)*(eta - xi + 2. * zeta + eta * zeta - xi * zeta - 2. * eta*xi*zeta + 1.)
				- x4 * (xi / 8. - 1. / 8.)*(eta + 1.)*(xi - eta + 2. * zeta - eta * zeta + xi * zeta - 2. * eta*xi*zeta + 1.)
				- x6 * (xi * xi / 4. - 1. / 4.)*(eta*zeta + 1.)*(eta - 1.)
				- x8 * (xi * xi / 4. - 1. / 4.)*(eta*zeta - 1.)*(eta + 1.)
				- x7 * (eta * eta - 1.)*(xi / 4. + 1. / 4.)*(xi*zeta - 1.) 
				- x9 * (eta * eta - 1.)*(xi / 4. - 1. / 4.)*(xi*zeta + 1.);

			Jacobian_[1][2] = 
				y5 * (zeta + 1. / 2.)
				+ y1 * (xi / 8. - 1. / 8.)*(eta - 1.)*(eta + xi + 2. * zeta + eta * zeta + xi * zeta + 2. * eta*xi*zeta + 1.)
				- y3 * (xi / 8. + 1. / 8.)*(eta + 1.)*(eta + xi - 2. * zeta + eta * zeta + xi * zeta - 2. * eta*xi*zeta - 1.) 
				- y10 * zeta*(xi / 2. - 1. / 2.)*(eta - 1.)
				+ y11 * zeta*(xi / 2. + 1. / 2.)*(eta - 1.)
				- y12 * zeta*(xi / 2. + 1. / 2.)*(eta + 1.)
				+ y13 * zeta*(xi / 2. - 1. / 2.)*(eta + 1.)
				- y2 * (xi / 8. + 1. / 8.)*(eta - 1.)*(eta - xi + 2. * zeta + eta * zeta - xi * zeta - 2. * eta*xi*zeta + 1.)
				- y4 * (xi / 8. - 1. / 8.)*(eta + 1.)*(xi - eta + 2. * zeta - eta * zeta + xi * zeta - 2. * eta*xi*zeta + 1.)
				- y6 * (xi * xi / 4. - 1. / 4.)*(eta*zeta + 1.)*(eta - 1.)
				- y8 * (xi * xi / 4. - 1. / 4.)*(eta*zeta - 1.)*(eta + 1.)
				- y7 * (eta * eta - 1.)*(xi / 4. + 1. / 4.)*(xi*zeta - 1.)
				- y9 * (eta * eta - 1.)*(xi / 4. - 1. / 4.)*(xi*zeta + 1.);

			Jacobian_[2][2] = 
				z5 * (zeta + 1. / 2.)
				+ z1 * (xi / 8. - 1. / 8.)*(eta - 1.)*(eta + xi + 2. * zeta + eta * zeta + xi * zeta + 2. * eta*xi*zeta + 1.)
				- z3 * (xi / 8. + 1. / 8.)*(eta + 1.)*(eta + xi - 2. * zeta + eta * zeta + xi * zeta - 2. * eta*xi*zeta - 1.)
				- z10 * zeta*(xi / 2. - 1. / 2.)*(eta - 1.)
				+ z11 * zeta*(xi / 2. + 1. / 2.)*(eta - 1.)
				- z12 * zeta*(xi / 2. + 1. / 2.)*(eta + 1.)
				+ z13 * zeta*(xi / 2. - 1. / 2.)*(eta + 1.)
				- z2 * (xi / 8. + 1. / 8.)*(eta - 1.)*(eta - xi + 2. * zeta + eta * zeta - xi * zeta - 2. * eta*xi*zeta + 1.)
				- z4 * (xi / 8. - 1. / 8.)*(eta + 1.)*(xi - eta + 2. * zeta - eta * zeta + xi * zeta - 2. * eta*xi*zeta + 1.)
				- z6 * (xi * xi / 4. - 1. / 4.)*(eta*zeta + 1.)*(eta - 1.)
				- z8 * (xi * xi / 4. - 1. / 4.)*(eta*zeta - 1.)*(eta + 1.)
				- z7 * (eta * eta - 1.)*(xi / 4. + 1. / 4.)*(xi*zeta - 1.)
				- z9 * (eta * eta - 1.)*(xi / 4. - 1. / 4.)*(xi*zeta + 1.);

			return Jacobian_;
		}

		virtual double
			evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_Jacobian(iso_point, elem_nodes);


			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] * Jacobian_[2][2] +
				Jacobian_[1][0] * Jacobian_[2][1] * Jacobian_[0][2] +
				Jacobian_[0][1] * Jacobian_[1][2] * Jacobian_[2][0] -
				Jacobian_[0][2] * Jacobian_[1][1] * Jacobian_[2][0] -
				Jacobian_[0][1] * Jacobian_[1][0] * Jacobian_[2][2] -
				Jacobian_[1][2] * Jacobian_[2][1] * Jacobian_[0][0];


			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_detJacobian(iso_point, elem_nodes);
			double& J11 = Jacobian_[0][0];
			double& J21 = Jacobian_[1][0];
			double& J31 = Jacobian_[2][0];
			double& J12 = Jacobian_[0][1];
			double& J22 = Jacobian_[1][1];
			double& J32 = Jacobian_[2][1];
			double& J13 = Jacobian_[0][2];
			double& J23 = Jacobian_[1][2];
			double& J33 = Jacobian_[2][2];

			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*(J22*J33 - J32 * J23);
			inv_Jacobian_[1][0] = (1. / det_Jacobian_)*(J31*J23 - J21 * J33);
			inv_Jacobian_[2][0] = (1. / det_Jacobian_)*(J21*J32 - J31 * J22);
			inv_Jacobian_[0][1] = (1. / det_Jacobian_)*(J32*J13 - J12 * J33);
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*(J11*J33 - J31 * J13);
			inv_Jacobian_[2][1] = (1. / det_Jacobian_)*(J31*J12 - J11 * J32);
			inv_Jacobian_[0][2] = (1. / det_Jacobian_)*(J12*J23 - J22 * J13);
			inv_Jacobian_[1][2] = (1. / det_Jacobian_)*(J21*J13 - J11 * J23);
			inv_Jacobian_[2][2] = (1. / det_Jacobian_)*(J11*J22 - J21 * J12);

			return inv_Jacobian_;
		}

	private:

		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 3x13
		std::vector<std::vector<double>> dNdx_;  // 3x13
		std::vector<std::vector<double>> Jacobian_; // 3x3
		std::vector<std::vector<double>> inv_Jacobian_; // 3x3
		double det_Jacobian_;

		ShapeFunction() :
			N_(13),
			dNdxi_(3, std::vector<double>(13)),
			dNdx_(3, std::vector<double>(13)),
			Jacobian_(3, std::vector<double>(3)),
			inv_Jacobian_(3, std::vector<double>(3)) {}

		ShapeFunction(const ShapeFunction&) {}

		ShapeFunction& operator=(const ShapeFunction&) {}
	};

	template<>
	class ShapeFunction< Dim3D, Prism15, Lagrange > :
		public LagrangeType<Dim3D>
	{
	public:

		using PointType = Point<Dim3D, CartesianCoordinate>;

		friend class SingletonHolder < ShapeFunction< Dim3D, Prism15, Lagrange > >;

		virtual std::vector<double>&
			evaluate_shape(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			N_[0] = 0.5*(xi + 1.)*(1. - eta - zeta)*(xi - 2. * eta - 2. * zeta);
			N_[1] = 0.5*eta*(1. + xi)*(2. * eta + xi - 2.);
			N_[2] = 0.5*zeta*(1. + xi)*(xi + 2. * zeta - 2.);
			N_[3] = 0.5*(xi - 1.)*(1. - eta - zeta)*(xi + 2. * eta + 2. * zeta);
			N_[4] = 0.5*eta*(1. - xi)*(2. * eta - xi - 2.);
			N_[5] = 0.5*zeta*(1. - xi)*(2. * zeta - xi - 2.);
			N_[6] = 2. * eta*(1. - eta - zeta)*(1. + xi);
			N_[7] = 2. * eta*zeta*(1. + xi);
			N_[8] = 2. * zeta*(1. - eta - zeta)*(1. + xi);
			N_[9] = 2. * eta*(1. - eta - zeta)*(1. - xi);
			N_[10] = 2. * eta*zeta*(1. - xi);
			N_[11] = 2. * zeta*(1. - eta - zeta)*(1. - xi);
			N_[12] = (1. - eta - zeta)*(1. - xi * xi);
			N_[13] = eta*(1. - xi * xi);
			N_[14] = zeta * (1. - xi * xi);
			
			return N_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdxi(const PointType& iso_point) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			dNdxi_[0][0] = (eta + zeta - 1.)*(eta + zeta - xi - 0.5);
			dNdxi_[0][1] = eta * (xi + eta - 0.5);
			dNdxi_[0][2] = zeta * (xi + zeta - 0.5);
			dNdxi_[0][3] = -(eta + zeta - 1.)*(xi + eta + zeta - 0.5);
			dNdxi_[0][4] = eta * (xi - eta + 0.5);
			dNdxi_[0][5] = zeta * (xi - zeta + 0.5);
			dNdxi_[0][6] = 2. * eta*(1. - eta - zeta);
			dNdxi_[0][7] = 2. * eta*zeta;
			dNdxi_[0][8] = 2. * zeta*(1. - eta - zeta);
			dNdxi_[0][9] = 2. * eta*(eta + zeta - 1.);
			dNdxi_[0][10] = -2. * eta*zeta;
			dNdxi_[0][11] = 2. * zeta*(eta + zeta - 1.);
			dNdxi_[0][12] = 2. * xi*(eta + zeta - 1.);
			dNdxi_[0][13] = -2. * xi*eta;
			dNdxi_[0][14] = -2. * xi*zeta;

			dNdxi_[1][0] = xi * (-0.5*xi + 2. * eta + 2. * zeta - 1.5) + 2. * eta + 2. * zeta - 1.;
			dNdxi_[1][1] = 0.5*(xi + 1.)*(xi + 4 * eta - 2.);
			dNdxi_[1][2] = 0;
			dNdxi_[1][3] = xi * (1.5 - 0.5*xi) + (2. - 2. * xi)*eta + (2. - 2. * xi)*zeta - 1.;
			dNdxi_[1][4] = 2. * (1. - xi)*eta + 0.5*(xi + 1.)*xi - 1.;
			dNdxi_[1][5] = 0.0;
			dNdxi_[1][6] = -2. * (xi + 1.)*(2. * eta + zeta - 1.);
			dNdxi_[1][7] = 2. * (xi + 1.)*zeta;
			dNdxi_[1][8] = -2. * (xi + 1.)*zeta;
			dNdxi_[1][9] = 2. * (xi - 1.)*(2. * eta + zeta - 1.);
			dNdxi_[1][10] = 2. * (1. - xi)*zeta;
			dNdxi_[1][11] = 2. * (xi - 1.)*zeta;
			dNdxi_[1][12] = xi * xi - 1.;
			dNdxi_[1][13] = 1. - xi * xi;
			dNdxi_[1][14] = 0.0;

			dNdxi_[2][0] = xi*(-0.5*xi + 2. * eta + 2. * zeta - 1.5) + 2. * eta + 2. * zeta - 1.;
			dNdxi_[2][1] = 0.0;
			dNdxi_[2][2] = 0.5*(xi + 1.)*(xi + 4 * zeta - 2.);
			dNdxi_[2][3] = xi*(1.5 - 0.5*xi) + (2. - 2. * xi)*eta + (2. - 2. * xi)*zeta - 1.;
			dNdxi_[2][4] = 0.0;
			dNdxi_[2][5] = 2. * (1. - xi)*zeta + (0.5*xi + 0.5)*xi - 1.;
			dNdxi_[2][6] = -2. * (xi + 1.)*eta;
			dNdxi_[2][7] = 2. * (xi + 1.)*eta;
			dNdxi_[2][8] = -2. * (xi + 1.)*(eta + 2. * zeta - 1.);
			dNdxi_[2][9] = 2. * (xi - 1.)*eta;
			dNdxi_[2][10] = 2. * (1. - xi)*eta;
			dNdxi_[2][11] = 2. * (xi - 1.)*(eta + 2. * zeta - 1.);
			dNdxi_[2][12] = xi*xi - 1.;
			dNdxi_[2][13] = 0.0;
			dNdxi_[2][14] = 1. - xi*xi;


			return dNdxi_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			evaluate_dNdxi(iso_point);
			evaluate_invJacobian(iso_point, elem_nodes);
			for (int j = 0; j < 10; ++j) {
				for (int i = 0; i < 3; ++i) {
					dNdx_[i][j] =
						dNdxi_[0][j] * inv_Jacobian_[i][0] +
						dNdxi_[1][j] * inv_Jacobian_[i][1] +
						dNdxi_[2][j] * inv_Jacobian_[i][2];
				}
			}

			return dNdx_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			const double& xi = iso_point.getX();
			const double& eta = iso_point.getY();
			const double& zeta = iso_point.getZ();

			const double& x1 = elem_nodes[0].getX();
			const double& x2 = elem_nodes[1].getX();
			const double& x3 = elem_nodes[2].getX();
			const double& x4 = elem_nodes[3].getX();
			const double& x5 = elem_nodes[4].getX();
			const double& x6 = elem_nodes[5].getX();
			const double& x7 = elem_nodes[6].getX();
			const double& x8 = elem_nodes[7].getX();
			const double& x9 = elem_nodes[8].getX();
			const double& x10 = elem_nodes[9].getX();
			const double& x11 = elem_nodes[10].getX();
			const double& x12 = elem_nodes[11].getX();
			const double& x13 = elem_nodes[12].getX();
			const double& x14 = elem_nodes[13].getX();
			const double& x15 = elem_nodes[14].getX();

			const double& y1 = elem_nodes[0].getY();
			const double& y2 = elem_nodes[1].getY();
			const double& y3 = elem_nodes[2].getY();
			const double& y4 = elem_nodes[3].getY();
			const double& y5 = elem_nodes[4].getY();
			const double& y6 = elem_nodes[5].getY();
			const double& y7 = elem_nodes[6].getY();
			const double& y8 = elem_nodes[7].getY();
			const double& y9 = elem_nodes[8].getY();
			const double& y10 = elem_nodes[9].getY();
			const double& y11 = elem_nodes[10].getY();
			const double& y12 = elem_nodes[11].getY();
			const double& y13 = elem_nodes[12].getY();
			const double& y14 = elem_nodes[13].getY();
			const double& y15 = elem_nodes[14].getY();

			const double& z1 = elem_nodes[0].getZ();
			const double& z2 = elem_nodes[1].getZ();
			const double& z3 = elem_nodes[2].getZ();
			const double& z4 = elem_nodes[3].getZ();
			const double& z5 = elem_nodes[4].getZ();
			const double& z6 = elem_nodes[5].getZ();
			const double& z7 = elem_nodes[6].getZ();
			const double& z8 = elem_nodes[7].getZ();
			const double& z9 = elem_nodes[8].getZ();
			const double& z10 = elem_nodes[9].getZ();
			const double& z11 = elem_nodes[10].getZ();
			const double& z12 = elem_nodes[11].getZ();
			const double& z13 = elem_nodes[12].getZ();
			const double& z14 = elem_nodes[13].getZ();
			const double& z15 = elem_nodes[14].getZ();
			      
			Jacobian_[0][0] = 
				x6 * zeta*(xi - zeta + 1 / 2)
				+ x1 * (eta + zeta - 1)*(eta - xi + zeta - 1 / 2)
				+ eta * x2*(eta + xi - 1 / 2)
				- 2 * eta*x7*(eta + zeta - 1)
				+ 2 * eta*x10*(eta + zeta - 1)
				+ 2 * x13*xi*(eta + zeta - 1)
				- 2 * x9*zeta*(eta + zeta - 1)
				+ 2 * x12*zeta*(eta + zeta - 1)
				+ x3 * zeta*(xi + zeta - 1 / 2)
				- 2 * eta*x14*xi
				+ 2 * eta*x8*zeta
				- 2 * eta*x11*zeta
				- 2 * x15*xi*zeta
				- x4 * (eta + zeta - 1)*(eta + xi + zeta - 1 / 2)
				+ eta * x5*(xi - eta + 1 / 2);

			Jacobian_[1][0] = 
				y6 * zeta*(xi - zeta + 1. / 2.)
				+ y1 * (eta + zeta - 1.)*(eta - xi + zeta - 1. / 2.)
				+ eta * y2*(eta + xi - 1. / 2.)
				- 2. * eta*y7*(eta + zeta - 1.)
				+ 2. * eta*y10*(eta + zeta - 1.)
				+ 2. * xi*y13*(eta + zeta - 1.)
				- 2. * y9*zeta*(eta + zeta - 1.)
				+ 2. * y12*zeta*(eta + zeta - 1.)
				+ y3 * zeta*(xi + zeta - 1. / 2.)
				- 2. * eta*xi*y14
				+ 2. * eta*y8*zeta 
				- 2. * eta*y11*zeta 
				- 2. * xi*y15*zeta 
				- y4 * (eta + zeta - 1.)*(eta + xi + zeta - 1. / 2.)
				+ eta * y5*(xi - eta + 1. / 2.)
				;

			Jacobian_[2][0] = 
				z6 * zeta*(xi - zeta + 1. / 2.)
				+ z1 * (eta + zeta - 1.)*(eta - xi + zeta - 1. / 2.)
				+ eta * z2*(eta + xi - 1. / 2.)
				- 2. * eta*z7*(eta + zeta - 1.)
				+ 2. * eta*z10*(eta + zeta - 1.)
				+ 2. * xi*z13*(eta + zeta - 1.)
				- 2. * z9*zeta*(eta + zeta - 1.)
				+ 2. * z12*zeta*(eta + zeta - 1.)
				+ z3 * zeta*(xi + zeta - 1. / 2.)
				- 2. * eta*xi*y14
				+ 2. * eta*z8*zeta
				- 2. * eta*z11*zeta
				- 2. * xi*y15*zeta
				- z4 * (eta + zeta - 1.)*(eta + xi + zeta - 1. / 2.)
				+ eta * z5*(xi - eta + 1. / 2.);

			Jacobian_[0][1] = 
				x13 * (xi * xi - 1.)
				- x14 * (xi * xi - 1.)
				- x5 * (eta*(2. * xi - 2.) - xi * (xi / 2. + 1. / 2.) + 1.)
				- x4 * (xi*(xi / 2. - 1.5) + zeta * (2. * xi - 2.) + eta * (2. * xi - 2.) + 1.)
				+ x1 * (2. * eta + 2. * zeta + xi * (2. * eta - xi / 2. + 2. * zeta - 1.5) - 1.)
				+ x8 * zeta*(2. * xi + 2.)
				- x9 * zeta*(2. * xi + 2.)
				- x11 * zeta*(2. * xi - 2.)
				+ x12 * zeta*(2. * xi - 2.)
				+ x2 * (xi / 2. + 1. / 2.)*(4. * eta + xi - 2.)
				- x7 * (2. * xi + 2.)*(2. * eta + zeta - 1.)
				+ x10 * (2. * xi - 2.)*(2. * eta + zeta - 1.);

			Jacobian_[1][1] = 
				y13 * (xi * xi - 1.)
				- y14 * (xi * xi - 1.)
				- y5 * (eta*(2. * xi - 2.) - xi * (xi / 2. + 1. / 2.) + 1.)
				- y4 * (xi*(xi / 2. - 1.5) + zeta * (2. * xi - 2.) + eta * (2. * xi - 2.) + 1.)
				+ y1 * (2. * eta + 2. * zeta + xi * (2. * eta - xi / 2. + 2. * zeta - 1.5) - 1.)
				+ y8 * zeta*(2. * xi + 2.)
				- y9 * zeta*(2. * xi + 2.)
				- y11 * zeta*(2. * xi - 2.)
				+ y12 * zeta*(2. * xi - 2.)
				+ y2 * (xi / 2. + 1. / 2.)*(4. * eta + xi - 2.)
				- y7 * (2. * xi + 2.)*(2. * eta + zeta - 1.)
				+ y10 * (2. * xi - 2.)*(2. * eta + zeta - 1.);

			Jacobian_[2][1] =
				z13 * (xi * xi - 1.)
				- y14 * (xi * xi - 1.)
				- z5 * (eta*(2. * xi - 2.) - xi * (xi / 2. + 1. / 2.) + 1.)
				- z4 * (xi*(xi / 2. - 1.5) + zeta * (2. * xi - 2.) + eta * (2. * xi - 2.) + 1.)
				+ z1 * (2. * eta + 2. * zeta + xi * (2. * eta - xi / 2. + 2. * zeta - 1.5) - 1.)
				+ z8 * zeta*(2. * xi + 2.)
				- z9 * zeta*(2. * xi + 2.)
				- z11 * zeta*(2. * xi - 2.)
				+ z12 * zeta*(2. * xi - 2.)
				+ z2 * (xi / 2. + 1. / 2.)*(4. * eta + xi - 2.)
				- z7 * (2. * xi + 2.)*(2. * eta + zeta - 1.)
				+ z10 * (2. * xi - 2.)*(2. * eta + zeta - 1.);
				
			Jacobian_[0][2] =
				x13 * (xi * xi - 1.)
				- x15 * (xi * xi - 1.)
				- x4 * (xi*(xi / 2. - 1.5) + zeta * (2. * xi - 2.) + eta * (2. * xi - 2.) + 1.)
				- x6 * (zeta*(2. * xi - 2.) - xi * (xi / 2. + 1. / 2.) + 1.)
				+ x1 * (2. * eta + 2. * zeta + xi * (2. * eta - xi / 2. + 2. * zeta - 1.5) - 1.)
				- eta * x7*(2. * xi + 2.)
				+ eta * x8*(2. * xi + 2.)
				+ eta * x10*(2. * xi - 2.)
				- eta * x11*(2. * xi - 2.)
				- x9 * (2. * xi + 2.)*(eta + 2. * zeta - 1.)
				+ x12 * (2. * xi - 2.)*(eta + 2. * zeta - 1.)
				+ x3 * (xi / 2. + 1. / 2.)*(xi + 4. * zeta - 2.);

			Jacobian_[1][2] =
				y13 * (xi * xi - 1.)
				- y15 * (xi * xi - 1.)
				- y4 * (xi*(xi / 2. - 1.5) + zeta * (2. * xi - 2.) + eta * (2. * xi - 2.) + 1.)
				- y6 * (zeta*(2. * xi - 2.) - xi * (xi / 2. + 1. / 2.) + 1.)
				+ y1 * (2. * eta + 2. * zeta + xi * (2. * eta - xi / 2. + 2. * zeta - 1.5) - 1.)
				- eta * y7*(2. * xi + 2.)
				+ eta * y8*(2. * xi + 2.)
				+ eta * y10*(2. * xi - 2.)
				- eta * y11*(2. * xi - 2.)
				- y9 * (2. * xi + 2.)*(eta + 2. * zeta - 1.)
				+ y12 * (2. * xi - 2.)*(eta + 2. * zeta - 1.)
				+ y3 * (xi / 2. + 1. / 2.)*(xi + 4. * zeta - 2.);

			Jacobian_[2][2] =
				z13 * (xi * xi - 1.)
				- y15 * (xi * xi - 1.)
				- z4 * (xi*(xi / 2. - 1.5)+ zeta * (2. * xi - 2.)+ eta * (2. * xi - 2.) + 1.)
				- z6 * (zeta*(2. * xi - 2.) - xi * (xi / 2. + 1. / 2.) + 1.)
				+ z1 * (2. * eta + 2. * zeta + xi * (2. * eta - xi / 2. + 2. * zeta - 1.5) - 1.)
				- eta * z7*(2. * xi + 2.)
				+ eta * z8*(2. * xi + 2.)
				+ eta * z10*(2. * xi - 2.)
				- eta * z11*(2. * xi - 2.)
				- z9 * (2. * xi + 2.)*(eta + 2. * zeta - 1.)
				+ z12 * (2. * xi - 2.)*(eta + 2. * zeta - 1.)
				+ z3 * (xi / 2. + 1. / 2.)*(xi + 4. * zeta - 2.)
				;
			return Jacobian_;
		}

		virtual double
			evaluate_detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_Jacobian(iso_point, elem_nodes);


			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] * Jacobian_[2][2] +
				Jacobian_[1][0] * Jacobian_[2][1] * Jacobian_[0][2] +
				Jacobian_[0][1] * Jacobian_[1][2] * Jacobian_[2][0] -
				Jacobian_[0][2] * Jacobian_[1][1] * Jacobian_[2][0] -
				Jacobian_[0][1] * Jacobian_[1][0] * Jacobian_[2][2] -
				Jacobian_[1][2] * Jacobian_[2][1] * Jacobian_[0][0];


			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			evaluate_invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			evaluate_detJacobian(iso_point, elem_nodes);
			double& J11 = Jacobian_[0][0];
			double& J21 = Jacobian_[1][0];
			double& J31 = Jacobian_[2][0];
			double& J12 = Jacobian_[0][1];
			double& J22 = Jacobian_[1][1];
			double& J32 = Jacobian_[2][1];
			double& J13 = Jacobian_[0][2];
			double& J23 = Jacobian_[1][2];
			double& J33 = Jacobian_[2][2];

			inv_Jacobian_[0][0] = (1. / det_Jacobian_)*(J22*J33 - J32 * J23);
			inv_Jacobian_[1][0] = (1. / det_Jacobian_)*(J31*J23 - J21 * J33);
			inv_Jacobian_[2][0] = (1. / det_Jacobian_)*(J21*J32 - J31 * J22);
			inv_Jacobian_[0][1] = (1. / det_Jacobian_)*(J32*J13 - J12 * J33);
			inv_Jacobian_[1][1] = (1. / det_Jacobian_)*(J11*J33 - J31 * J13);
			inv_Jacobian_[2][1] = (1. / det_Jacobian_)*(J31*J12 - J11 * J32);
			inv_Jacobian_[0][2] = (1. / det_Jacobian_)*(J12*J23 - J22 * J13);
			inv_Jacobian_[1][2] = (1. / det_Jacobian_)*(J21*J13 - J11 * J23);
			inv_Jacobian_[2][2] = (1. / det_Jacobian_)*(J11*J22 - J21 * J12);

			return inv_Jacobian_;
		}

	private:

		std::vector<double> N_;
		std::vector<std::vector<double>> dNdxi_; // 3x15
		std::vector<std::vector<double>> dNdx_;  // 3x15
		std::vector<std::vector<double>> Jacobian_; // 3x3
		std::vector<std::vector<double>> inv_Jacobian_; // 3x3
		double det_Jacobian_;

		ShapeFunction() :
			N_(15),
			dNdxi_(3, std::vector<double>(15)),
			dNdx_(3, std::vector<double>(15)),
			Jacobian_(3, std::vector<double>(3)),
			inv_Jacobian_(3, std::vector<double>(3)) {}

		ShapeFunction(const ShapeFunction&) {}

		ShapeFunction& operator=(const ShapeFunction&) {}
	};



}




#endif //ARTCFD_LAGRANGE3DSHAPEFUNCTION_H










