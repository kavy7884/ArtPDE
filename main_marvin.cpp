

#include <iostream>
# include "FunctionSpace\lagrange3d_shape_function.h"
# include "FunctionSpace\ShapeFunctionFactory.h"

# include <random>

using namespace art_pde;

auto& Linearfactory2 =
SingletonHolder<ShapeFunctionFactory<LagrangeType<Dim2D>, ElementType2D>>::instance();

auto& Linearfactory =
SingletonHolder<ShapeFunctionFactory<LagrangeType<Dim3D>, ElementType3D>>::instance();

double rand1(void){
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	return dist(mt);
}

double rand2(void){
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(-1.0, 1.0);
	return dist(mt);
}

class PolyQ4{
public:
	using PointType = Point<Dim2D, CartesianCoordinate>;

	PolyQ4() :
		basis(Linearfactory2.getInstance(ElementType2D::Q4))
	{

		NUM = 4;
		elem_iso_coor = {
			PointType(-1, -1),
			PointType(1, -1),
			PointType(1, 1),
			PointType(-1, 1)
		};
		elem_phy_coor = {
			PointType(0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.5 + 0.2*rand1()),
			PointType(0.2*rand1(), 0.5 + 0.2*rand1()),
		};
		for (int i = 0; i < NUM; ++i){
			coeff.push_back(100.0*rand2());
		}


	}

	void testDelta(){
		for (int i = 0; i < NUM; ++i){
			auto& N = basis.evaluate_shape(elem_iso_coor[i]);
			for (int j = 0; j < NUM; ++j){
				std::cout << N[j] << " ,";
			}
			std::cout << "\n";
		}
		int x;
	}

	void testUnity(void){

		PointType trial_iso(0.5*rand1(), 0.5*rand1());

		auto& N = basis.evaluate_shape(trial_iso);
		double val = 0.0;
		for (int i = 0; i < NUM; ++i){
			val += N[i];
		}

		std::cout << val << "\n";

	}

	void testIsoInterpolation(void){

		PointType trial_iso(0.5*rand1(), 0.5*rand1());

		auto& N = basis.evaluate_shape(trial_iso);

		double approx = 0.0;
		for (int i = 0; i < NUM; ++i){
			approx += N[i] * polyValue(elem_iso_coor[i]);
		}

		double exact = polyValue(trial_iso);

		std::cout << exact*exact - approx*approx << "\n";
	}

	void testIsoGrad(void){

		PointType trial_iso(rand2(), rand2());

		auto& dNdxi = basis.evaluate_dNdxi(trial_iso);

		std::vector<double> approx(2, 0.0);
		for (int i = 0; i < NUM; ++i){
			approx[0] += dNdxi[0][i] * polyValue(elem_iso_coor[i]);
			approx[1] += dNdxi[1][i] * polyValue(elem_iso_coor[i]);

		}

		auto exact = polyGrad(trial_iso);

		std::cout <<
			(exact[0] * exact[0] - approx[0] * approx[0]) << ", " <<
			(exact[1] * exact[1] - approx[1] * approx[1]) << "\n";
	}

	void testPhyGrad(void){

		PointType trial_iso(rand2(), rand2());

		auto& N = basis.evaluate_shape(trial_iso);

		double trial_phy_x{ 0.0 };
		double trial_phy_y{ 0.0 };

		for (int i = 0; i < NUM; ++i){
			trial_phy_x += N[i] * elem_phy_coor[i].getX();
			trial_phy_y += N[i] * elem_phy_coor[i].getY();
		}

		PointType trial_phy(trial_phy_x, trial_phy_y);

		auto& dNdx = basis.evaluate_dNdx(trial_iso, elem_phy_coor);

		double approx = 0.0;

		for (int i = 0; i < NUM; ++i){
			approx += dNdx[i] * polyValue
		}

		auto exact = polyGrad(trial_phy);

		std::cout <<
			(exact[0] * exact[0] - approx[0] * approx[0]) << ", " <<
			(exact[1] * exact[1] - approx[1] * approx[1]) << "\n";
	}

private:

	double polyValue(PointType& coor){
		double x = coor.getX();
		double y = coor.getY();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * x*y;
	}

	std::vector<double> polyGrad(PointType& coor){
		double x = coor.getX();
		double y = coor.getY();
		std::vector<double> grad(2);
		grad[0] = coeff[1] + coeff[3] * y;
		grad[1] = coeff[2] + coeff[3] * x;
		return grad;
	}
	PointType centroid(std::vector<PointType>& pts){
		PointType avg(0.0, 0.0);
		for (auto i : pts){
			avg += i;
		}
		avg /= static_cast<double>(NUM);
		return avg;
	}
	std::vector<double> coeff;
	std::vector<PointType> elem_iso_coor;
	std::vector<PointType> elem_phy_coor;
	LagrangeType<Dim2D>& basis;
	int NUM;
};

class PolyTET4{
public:
	using PointType = Point<Dim3D, CartesianCoordinate>;

	PolyTET4() :
		basis(Linearfactory.getInstance(ElementType3D::Tetra4))
	{
		elem_iso_coor = {
			PointType(1, 0, 0),
			PointType(0, 1, 0),
			PointType(0, 0, 1),
			PointType(0, 0, 0)
		};
		elem_phy_coor = {
			PointType(0.2*rand1(), 0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.2*rand1(), 0.2*rand1()),
			PointType(0.2*rand1(), 0.5 + 0.2*rand1(), 0.2*rand1()),
			PointType(0.2*rand1(), 0.2*rand1(), 0.5 + 0.2*rand1()),
		};
		for (int i = 0; i < elem_iso_coor.size(); ++i){
			coeff.push_back(100.0*rand2());
		}

		NUM = 4;
	}

	void testDelta(){
		for (int i = 0; i < NUM; ++i){
			auto& N = basis.evaluate_shape(elem_iso_coor[i]);
			for (int j = 0; j < NUM; ++j){
				std::cout << N[j] << " ,";
			}
			std::cout << "\n";
		}
	}

	void testUnity(void){

		PointType trial_iso(0.5*rand1(), 0.5*rand1(), 0.5*rand1());

		auto& N = basis.evaluate_shape(trial_iso);
		double val = 0.0;
		for (int i = 0; i < NUM; ++i){
			val += N[i];
		}

		std::cout << val << "\n";

	}

	void testIsoInterpolation(void){

		PointType trial_iso(0.5*rand1(), 0.5*rand1(), 0.5*rand1());

		auto& N = basis.evaluate_shape(trial_iso);

		double approx = 0.0;
		for (int i = 0; i < NUM; ++i){
			approx += N[i] * polyValue(elem_iso_coor[i]);
		}

		double exact = polyValue(trial_iso);

		std::cout << exact*exact - approx*approx << "\n";
	}

	void testIsoGrad(void){

		PointType trial_iso(0.5*rand2(), 0.5*rand2(), 0.5*rand2());

		auto& dNdxi = basis.evaluate_dNdxi(trial_iso);

		std::vector<double> approx(3, 0.0);
		for (int i = 0; i < NUM; ++i){
			approx[0] += dNdxi[0][i] * polyValue(elem_iso_coor[i]);
			approx[1] += dNdxi[1][i] * polyValue(elem_iso_coor[i]);
			approx[2] += dNdxi[2][i] * polyValue(elem_iso_coor[i]);
		}

		auto exact = polyGrad(trial_iso);

		std::cout <<
			(exact[0] * exact[0] - approx[0] * approx[0]) << ", " <<
			(exact[1] * exact[1] - approx[1] * approx[1]) << ", " <<
			(exact[2] * exact[2] - approx[2] * approx[2]) << "\n";
	}

	void testPhyGrad(void){

		PointType trial_phy = centroid(elem_iso_coor);

		auto& dNdx = basis.evaluate_dNdx(trial_phy, elem_phy_coor);

		std::vector<double> approx(3, 0.0);

		for (int i = 0; i < NUM; ++i){
			approx[0] += dNdx[0][i] * polyValue(elem_phy_coor[i]);
			approx[1] += dNdx[1][i] * polyValue(elem_phy_coor[i]);
			approx[2] += dNdx[2][i] * polyValue(elem_phy_coor[i]);
		}

		auto exact = polyGrad(trial_phy);

		std::cout <<
			(exact[0] * exact[0] - approx[0] * approx[0]) << ", " <<
			(exact[1] * exact[1] - approx[1] * approx[1]) << ", " <<
			(exact[2] * exact[2] - approx[2] * approx[2]) << "\n";
	}

private:
	double polyValue(PointType& coor){
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * z;
	}

	std::vector<double> polyGrad(PointType& coor){
		std::vector<double> grad(3);
		grad[0] = coeff[1];
		grad[1] = coeff[2];
		grad[2] = coeff[3];
		return grad;
	}
	PointType centroid(std::vector<PointType>& pts){
		PointType avg(0.0, 0.0, 0.0);
		for (auto i : pts){
			avg += i;
		}
		avg /= static_cast<double>(NUM);
		return avg;
	}
	std::vector<double> coeff;
	std::vector<PointType> elem_iso_coor;
	std::vector<PointType> elem_phy_coor;
	LagrangeType<Dim3D>& basis;
	int NUM;
};

class PolyT3{
public:
	using PointType = Point<Dim2D, CartesianCoordinate>;
	PolyT3() :
		basis(Linearfactory2.getInstance(ElementType2D::T3))
	{

		NUM = 3;
		elem_iso_coor = {
			PointType(0, 0),
			PointType(1, 0),
			PointType(0, 1)
		};
		elem_phy_coor = {
			PointType(0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.2*rand1()),
			PointType(0.2*rand1(), 0.5 + 0.2*rand1())
		};
		for (int i = 0; i < NUM; ++i){
			coeff.push_back(100.0*rand2());
		}


	}

	void testDelta(){
		for (int i = 0; i < NUM; ++i){
			auto& N = basis.evaluate_shape(elem_iso_coor[i]);
			for (int j = 0; j < NUM; ++j){
				std::cout << N[j] << " ,";
			}
			std::cout << "\n";
		}
		int x;
	}

	void testUnity(void){

		PointType trial_iso(0.5*rand1(), 0.5*rand1());

		auto& N = basis.evaluate_shape(trial_iso);
		double val = 0.0;
		for (int i = 0; i < NUM; ++i){
			val += N[i];
		}

		std::cout << val << "\n";

	}

	void testIsoInterpolation(void){

		PointType trial_iso(0.5*rand1(), 0.5*rand1());

		auto& N = basis.evaluate_shape(trial_iso);

		double approx = 0.0;
		for (int i = 0; i < NUM; ++i){
			approx += N[i] * polyValue(elem_iso_coor[i]);
		}

		double exact = polyValue(trial_iso);

		std::cout << exact*exact - approx*approx << "\n";
	}

	void testIsoGrad(void){

		PointType trial_iso(rand1(), rand1());

		auto& dNdxi = basis.evaluate_dNdxi(trial_iso);

		std::vector<double> approx(2, 0.0);
		for (int i = 0; i < NUM; ++i){
			approx[0] += dNdxi[0][i] * polyValue(elem_iso_coor[i]);
			approx[1] += dNdxi[1][i] * polyValue(elem_iso_coor[i]);

		}

		auto exact = polyGrad(trial_iso);

		std::cout <<
			(exact[0] * exact[0] - approx[0] * approx[0]) << ", " <<
			(exact[1] * exact[1] - approx[1] * approx[1]) << "\n";
	}

	void testPhyGrad(void){

		PointType trial_phy = centroid(elem_iso_coor);

		auto& dNdx = basis.evaluate_dNdx(trial_phy, elem_phy_coor);

		std::vector<double> approx(2, 0.0);

		for (int i = 0; i < NUM; ++i){
			approx[0] += dNdx[0][i] * polyValue(elem_phy_coor[i]);
			approx[1] += dNdx[1][i] * polyValue(elem_phy_coor[i]);

		}

		auto exact = polyGrad(trial_phy);

		std::cout <<
			(exact[0] * exact[0] - approx[0] * approx[0]) << ", " <<
			(exact[1] * exact[1] - approx[1] * approx[1]) << "\n";
	}

private:
	double polyValue(PointType& coor){
		double x = coor.getX();
		double y = coor.getY();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y;
	}

	std::vector<double> polyGrad(PointType& coor){
		double x = coor.getX();
		double y = coor.getY();
		std::vector<double> grad(2);
		grad[0] = coeff[1];
		grad[1] = coeff[2];
		return grad;
	}
	PointType centroid(std::vector<PointType>& pts){
		PointType avg(0.0, 0.0);
		for (auto i : pts){
			avg += i;
		}
		avg /= static_cast<double>(NUM);
		return avg;
	}
	std::vector<double> coeff;
	std::vector<PointType> elem_iso_coor;
	std::vector<PointType> elem_phy_coor;
	LagrangeType<Dim2D>& basis;
	int NUM;
};

int main() {

	PolyQ4 TET4;
	TET4.testDelta();
	TET4.testUnity();
	TET4.testIsoInterpolation();
	TET4.testIsoGrad();
	TET4.testPhyGrad();



	system("pause");
	return 0;
}
