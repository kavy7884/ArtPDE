
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
			N(const PointType& iso_point) override
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
			dNdxi(const PointType& iso_point) override
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
			dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			dNdxi(iso_point);
			invJacobian(iso_point, elem_nodes);
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
			Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
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
			detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
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
			invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
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
			detJacobian(iso_point, elem_nodes);
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
			N(const PointType& iso_point) override
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
			dNdxi(const PointType& iso_point) override
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
			dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			dNdxi(iso_point);
			invJacobian(iso_point, elem_nodes);
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
			Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
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
			detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			Jacobian(iso_point, elem_nodes);
			
			det_Jacobian_ = Jacobian_[0][0]*Jacobian_[1][1]*Jacobian_[2][2] + 
				     Jacobian_[1][0]*Jacobian_[2][1]*Jacobian_[0][2] + 
				     Jacobian_[0][1]*Jacobian_[1][2]*Jacobian_[2][0] - 
				     Jacobian_[0][2]*Jacobian_[1][1]*Jacobian_[2][0] - 
				     Jacobian_[0][1]*Jacobian_[1][0]*Jacobian_[2][2] - 
				     Jacobian_[1][2]*Jacobian_[2][1]*Jacobian_[0][0];
				

			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			
			detJacobian(iso_point, elem_nodes);
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
			N(const PointType& iso_point) override
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
			dNdxi(const PointType& iso_point) override
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
			dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			dNdxi(iso_point);
			invJacobian(iso_point, elem_nodes);
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
			Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
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
			detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			Jacobian(iso_point, elem_nodes);


			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] * Jacobian_[2][2] +
				Jacobian_[1][0] * Jacobian_[2][1] * Jacobian_[0][2] +
				Jacobian_[0][1] * Jacobian_[1][2] * Jacobian_[2][0] -
				Jacobian_[0][2] * Jacobian_[1][1] * Jacobian_[2][0] -
				Jacobian_[0][1] * Jacobian_[1][0] * Jacobian_[2][2] -
				Jacobian_[1][2] * Jacobian_[2][1] * Jacobian_[0][0];


			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			detJacobian(iso_point, elem_nodes);
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
			N(const PointType& iso_point) override
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
			dNdxi(const PointType& iso_point) override
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
			dNdx(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{
			dNdxi(iso_point);
			invJacobian(iso_point, elem_nodes);
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
			Jacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
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
			detJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			Jacobian(iso_point, elem_nodes);


			det_Jacobian_ = Jacobian_[0][0] * Jacobian_[1][1] * Jacobian_[2][2] +
				Jacobian_[1][0] * Jacobian_[2][1] * Jacobian_[0][2] +
				Jacobian_[0][1] * Jacobian_[1][2] * Jacobian_[2][0] -
				Jacobian_[0][2] * Jacobian_[1][1] * Jacobian_[2][0] -
				Jacobian_[0][1] * Jacobian_[1][0] * Jacobian_[2][2] -
				Jacobian_[1][2] * Jacobian_[2][1] * Jacobian_[0][0];


			return det_Jacobian_;
		}

		virtual std::vector<std::vector<double>>&
			invJacobian(const PointType& iso_point, const std::vector<PointType>& elem_nodes) override
		{

			detJacobian(iso_point, elem_nodes);
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

}


#endif //ARTCFD_LAGRANGE3DSHAPEFUNCTION_H










