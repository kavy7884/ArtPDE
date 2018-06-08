
#include <iostream>
#include <Eigen/Dense>
#include <random>

#include "FunctionSpace\basis_function_lagrange.h"


namespace  dim2 = art_pde::function_space::isoparametric::Dim2D;
namespace  dim3 = art_pde::function_space::isoparametric::Dim3D;

class Dim2D{ public: const static int num = 2; };
class Dim3D{ public: const static int num = 3; };

class Coor2{
public:
	using dimension_tag = ::Dim2D;
	Coor2() = default;
	Coor2(double x_, double y_) :x(x_), y(y_){
	}
	double const getX() const { return x; }
	double const getY() const { return y; }
private:
	double x{ 0.0 }, y{ 0.0 };
};

class Coor3{
public:
	using dimension_tag = ::Dim3D;
	Coor3() = default;
	Coor3(double x_, double y_, double z_) :
		x(x_), y(y_), z(z_){
	}
	double const getX()const{ return x; }
	double const getY()const{ return y; }
	double const getZ()const{ return z; }
private:
	double x{ 0.0 }, y{ 0.0 }, z{ 0.0 };
};


auto& basis_func2d =
art_pde::function_space::SingletonHolder<dim2::BasisFunctionFactory<Coor2>>::instance();

auto& basis_func3d =
art_pde::function_space::SingletonHolder<dim3::BasisFunctionFactory<Coor3>>::instance();


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

template<class PointT, template<class> class BasisT >
class TestBase{
public:

	using PointType = PointT;

	TestBase(BasisT<PointType>& basis_) :basis(basis_){
		dim = PointT::dimension_tag::num;
	}

	void testDelta(){
		std::cout << "testDelta=====================================\n\n";
		for (size_t i = 0; i < NUM; ++i){
			auto& N = basis.evaluate_shape(elem_iso_coor[i]);
			for (size_t j = 0; j < NUM; ++j){
				std::cout << N[j] << " ,";
			}
			std::cout << "\n\n";
		}
	}

	void testUnity(void){

		PointType trial_iso = getIsoPoint();

		auto& N = basis.evaluate_shape(trial_iso);
		double val = 0.0;
		for (size_t i = 0; i < NUM; ++i){
			val += N[i];
		}

		std::cout << "testUnity=====================================\n\n";
		std::cout << val << "\n\n";

	}

	void testIsoInterpolation(void){

		PointType trial_iso = getIsoPoint();

		auto& N = basis.evaluate_shape(trial_iso);

		double approx = 0.0;
		for (size_t i = 0; i < NUM; ++i){
			approx += N[i] * polyValue(elem_iso_coor[i]);
		}

		double exact = polyValue(trial_iso);
		std::cout << "testIsoIntp===================================\n\n";
		std::cout << exact*exact - approx*approx << "\n\n";
	}

	void testIsoGrad(void){

		PointType trial_iso = getIsoPoint();

		auto& dNdxi = basis.evaluate_dNdxi(trial_iso);

		std::vector<double> approx(dim, 0.0);
		for (size_t i = 0; i < NUM; ++i){
			for (size_t j = 0; j < dim; ++j){
				approx[j] += dNdxi[j][i] * polyValue(elem_iso_coor[i]);
			}

		}

		auto exact = polyGrad(trial_iso);
		std::cout << "testIsoGrad==================================\n\n" << "[";
		for (size_t i = 0; i < dim; ++i){
			std::cout <<
				(exact[i] * exact[i] - approx[i] * approx[i]) << " ";
		}
		std::cout << " ]\n\n";
	}

	void testJacobian(){

		PointType trial_iso = getIsoPoint();
		Eigen::MatrixXd J_ = getJ(tag, trial_iso);

		auto& J = basis.evaluate_Jacobian(trial_iso, elem_phy_coor);
		std::cout << "testJ=====================================\n\n";
		for (size_t i = 0; i < dim; ++i){
			std::cout << "[ ";
			for (size_t j = 0; j < dim; ++j){
				std::cout << J[i][j] - J_(i, j) << " ";
			}
			std::cout << " ]" << "\n\n";
		}

	}

	void testInvJ(){

		PointType trial_iso = getIsoPoint();
		Eigen::MatrixXd J_ = getJ(tag, trial_iso);

		auto& invJ_ = J_.inverse();
		auto& invJ = basis.evaluate_invJacobian(trial_iso, elem_phy_coor);

		std::cout << "testInvJ=====================================\n\n";
		for (size_t i = 0; i < dim; ++i){
			std::cout << "[ ";
			for (size_t j = 0; j < dim; ++j){
				std::cout << invJ_(i, j) - invJ[i][j] << " ";
			}
			std::cout << " ]" << "\n\n";
		}

	}

	void testdetJ(){

		PointType trial_iso = getIsoPoint();
		Eigen::MatrixXd J_ = getJ(tag, trial_iso);

		double detJ_ = J_.determinant();
		double detJ = basis.evaluate_detJacobian(trial_iso, elem_phy_coor);
		std::cout << "testdetJ=====================================\n\n";
		std::cout << detJ_ - detJ << "\n\n";

	}

protected:

	Eigen::MatrixXd getJ(::Dim2D, PointType& trial_iso){

		Eigen::MatrixXd J(dim, dim);
		J.setZero();
		auto& dNdxi = basis.evaluate_dNdxi(trial_iso);

		for (size_t i = 0; i < NUM; ++i){
			J(0, 0) += dNdxi[0][i] * elem_phy_coor[i].getX();
			J(1, 0) += dNdxi[1][i] * elem_phy_coor[i].getX();

			J(0, 1) += dNdxi[0][i] * elem_phy_coor[i].getY();
			J(1, 1) += dNdxi[1][i] * elem_phy_coor[i].getY();
		}
		return J;
	}

	Eigen::MatrixXd getJ(::Dim3D, PointType& trial_iso){

		Eigen::MatrixXd J(dim, dim);
		J.setZero();
		auto& dNdxi = basis.evaluate_dNdxi(trial_iso);

		for (int i = 0; i < NUM; ++i){
			J(0, 0) += dNdxi[0][i] * elem_phy_coor[i].getX();
			J(1, 0) += dNdxi[1][i] * elem_phy_coor[i].getX();
			J(2, 0) += dNdxi[2][i] * elem_phy_coor[i].getX();

			J(0, 1) += dNdxi[0][i] * elem_phy_coor[i].getY();
			J(1, 1) += dNdxi[1][i] * elem_phy_coor[i].getY();
			J(2, 1) += dNdxi[2][i] * elem_phy_coor[i].getY();

			J(0, 2) += dNdxi[0][i] * elem_phy_coor[i].getZ();
			J(1, 2) += dNdxi[1][i] * elem_phy_coor[i].getZ();
			J(2, 2) += dNdxi[2][i] * elem_phy_coor[i].getZ();
		}
		return J;
	}

	virtual PointType getIsoPoint() = 0;
	virtual double polyValue(PointType& coor) = 0;
	virtual std::vector<double> polyGrad(PointType& coor) = 0;

	typename PointT::dimension_tag tag;
	std::size_t dim;
	std::size_t NUM{ 1000 };
	std::vector<double> coeff;
	std::vector<PointType> elem_iso_coor;
	std::vector<PointType> elem_phy_coor;
	BasisT<PointType>& basis;
};

class TestQ4 :public TestBase<Coor2, dim2::BasisFunction>{

public:
	TestQ4() :
		TestBase<Coor2, dim2::BasisFunction>(basis_func2d.getInstance(dim2::ElementType::Q4)){
		init();
	}
protected:
	void init() {
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
		for (size_t i = 0; i < NUM; ++i){
			coeff.push_back(1.0*rand2());
		}
	}
	virtual double polyValue(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * x*y;
	}
	virtual std::vector<double> polyGrad(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();
		std::vector<double> grad(2);
		grad[0] = coeff[1] + coeff[3] * y;
		grad[1] = coeff[2] + coeff[3] * x;
		return grad;
	}
	virtual PointType getIsoPoint() override{
		return PointType(rand2(), rand2());
	}
};

class TestQ9 :public TestBase<Coor2, dim2::BasisFunction>{

public:
	TestQ9() :
		TestBase<Coor2, dim2::BasisFunction>(basis_func2d.getInstance(dim2::ElementType::Q9)){
		init();
	}
protected:
	void init() {
		NUM = 9;
		elem_iso_coor = {
			PointType(-1, -1),
			PointType(1, -1),
			PointType(1, 1),
			PointType(-1, 1),
			PointType(0, -1),
			PointType(1, 0),
			PointType(0, 1),
			PointType(-1, 0),
			PointType(0, 0)
		};
		elem_phy_coor = {
			PointType(0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.5 + 0.2*rand1()),
			PointType(0.2*rand1(), 0.5 + 0.2*rand1()),

			PointType(0.25 + 0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.25 + 0.2*rand1()),
			PointType(0.25 + 0.2*rand1(), 0.25 + 0.2*rand1()),
			PointType(0.2*rand1(), 0.25 + 0.2*rand1()),
			PointType(0.2*rand1(), 0.2*rand1())
		};
		for (size_t i = 0; i < NUM; ++i){
			coeff.push_back(1.0*rand2());
		}
	}
	virtual double polyValue(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * x*y
			+ coeff[4] * x*x
			+ coeff[5] * y*y
			+ coeff[6] * x*x*y
			+ coeff[7] * x*y*y
			+ coeff[8] * x*x*y*y;
	}
	virtual std::vector<double> polyGrad(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();
		std::vector<double> grad(2);
		grad[0] = 
			coeff[1] + coeff[3] * y + 2.*coeff[4] * x 
			+ 2.*coeff[6] * x*y + coeff[7] * y*y + 2.*coeff[8]*x*y*y;
		grad[1] = 
			coeff[2] + coeff[3] * x + 2.*coeff[5]*y+coeff[6]*x*x
			+ 2.*coeff[7] * x*y + 2.*coeff[8]*x*x*y;
		return grad;
	}
	virtual PointType getIsoPoint() override{
		return PointType(rand2(), rand2());
	}
};

class TestQ8 :public TestBase<Coor2, dim2::BasisFunction>{

public:
	TestQ8() :
		TestBase<Coor2, dim2::BasisFunction>(basis_func2d.getInstance(dim2::ElementType::Q8)){
		init();
	}
protected:
	void init() {
		NUM = 8;
		elem_iso_coor = {
			PointType(-1, -1),
			PointType(1, -1),
			PointType(1, 1),
			PointType(-1, 1),
			PointType(0, -1),
			PointType(1, 0),
			PointType(0, 1),
			PointType(-1, 0)
		};
		elem_phy_coor = {
			PointType(0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.5 + 0.2*rand1()),
			PointType(0.2*rand1(), 0.5 + 0.2*rand1()),

			PointType(0.25 + 0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.25 + 0.2*rand1()),
			PointType(0.25 + 0.2*rand1(), 0.25 + 0.2*rand1()),
			PointType(0.2*rand1(), 0.25 + 0.2*rand1())
		};
		for (size_t i = 0; i < NUM; ++i){
			coeff.push_back(1.0*rand2());
		}
	}
	virtual double polyValue(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * x*y
			+ coeff[4] * x*x
			+ coeff[5] * y*y
			+ coeff[6] * x*x*y
			+ coeff[7] * x*y*y;
	}
	virtual std::vector<double> polyGrad(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();
		std::vector<double> grad(2);
		grad[0] =
			coeff[1] + coeff[3] * y + 2.*coeff[4] * x
			+ 2.*coeff[6] * x*y + coeff[7] * y*y;
		grad[1] =
			coeff[2] + coeff[3] * x + 2.*coeff[5] * y + coeff[6] * x*x
			+ 2.*coeff[7] * x*y;
		return grad;
	}
	virtual PointType getIsoPoint() override{
		return PointType(rand2(), rand2());
	}
};

class TestT3 :public TestBase<Coor2, dim2::BasisFunction>{

public:
	TestT3() :
		TestBase<Coor2, dim2::BasisFunction>(basis_func2d.getInstance(dim2::ElementType::T3)){
		init();
	}
protected:
	void init() {
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
		for (size_t i = 0; i < NUM; ++i){
			coeff.push_back(1.0*rand2());
		}
	}
	virtual double polyValue(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y;
	}

	virtual std::vector<double> polyGrad(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();
		std::vector<double> grad(2);
		grad[0] = coeff[1];
		grad[1] = coeff[2];
		return grad;
	}
	virtual PointType getIsoPoint() override{
		return PointType(0.5*rand1(), 0.5*rand1());
	}
};

class TestTET4 :public TestBase<Coor3, dim3::BasisFunction>{

public:
	TestTET4() :
		TestBase<Coor3, dim3::BasisFunction>(basis_func3d.getInstance(dim3::ElementType::Tetra4)){
		init();
	}
protected:
	void init() {
		NUM = 4;
		elem_iso_coor = {
			PointType(1, 0, 0),
			PointType(0, 1, 0),
			PointType(0, 0, 1),
			PointType(0, 0, 0),
		};
		elem_phy_coor = {
			PointType(0.5 + 0.2*rand1(), 0.2*rand1(), 0.2*rand1()),
			PointType(0.2*rand1(), 0.5 + 0.2*rand1(), 0.2*rand1()),
			PointType(0.2*rand1(), 0.2*rand1(), 0.5 + 0.2*rand1()),
			PointType(0.2*rand1(), 0.2*rand1(), 0.2*rand1())
		};
		for (size_t i = 0; i < NUM; ++i){
			coeff.push_back(1.0*rand2());
		}
	}
	virtual double polyValue(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * z;
	}

	virtual std::vector<double> polyGrad(PointType& coor) override{

		std::vector<double> grad(3);
		grad[0] = coeff[1];
		grad[1] = coeff[2];
		grad[2] = coeff[3];
		return grad;
	}
	virtual PointType getIsoPoint() override{
		auto a = PointType(0.5*rand1(), 0.5*rand1(), 0.5*rand1());
		//std::cout << a.getX() << " " << a.getY() << " " << a.getZ() << "\n";
		return a;
	}
};

class TestHEXA8 :public TestBase<Coor3, dim3::BasisFunction>{

public:
	TestHEXA8() :
		TestBase<Coor3, dim3::BasisFunction>(basis_func3d.getInstance(dim3::ElementType::Hexa8)){
		init();
	}
protected:
	void init() {
		NUM = 8;
		elem_iso_coor = {
			PointType(-1, -1, -1),
			PointType(+1, -1, -1),
			PointType(+1, +1, -1),
			PointType(-1, +1, -1),
			PointType(-1, -1, +1),
			PointType(+1, -1, +1),
			PointType(+1, +1, +1),
			PointType(-1, +1, +1)
		};
		elem_phy_coor = {
			PointType(0.2*rand2(), 0.2*rand2(), 0.2*rand2()),
			PointType(0.5 + 0.2*rand2(), 0.2*rand2(), 0.2*rand2()),
			PointType(0.5 + 0.2*rand2(), 0.5 + 0.2*rand2(), 0.2*rand2()),
			PointType(0.2*rand2(), 0.5 + 0.2*rand2(), 0.2*rand2()),
			PointType(0.2*rand2(), 0.2*rand2(), 0.5 + 0.2*rand2()),
			PointType(0.5 + 0.2*rand2(), 0.2*rand2(), 0.5 + 0.2*rand2()),
			PointType(0.5 + 0.2*rand2(), 0.5 + 0.2*rand2(), 0.5 + 0.2*rand2()),
			PointType(0.2*rand2(), 0.5 + 0.2*rand2(), 0.5 + 0.2*rand2())

		};
		for (size_t i = 0; i < NUM; ++i){
			coeff.push_back(1.0*rand2());
		}
	}
	virtual double polyValue(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * z
			+ coeff[4] * x*y
			+ coeff[5] * z*y
			+ coeff[6] * x*z
			+ coeff[7] * x*y*z;
	}

	virtual std::vector<double> polyGrad(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();
		std::vector<double> grad(3);

		grad[0] = coeff[1] + coeff[4] * y + coeff[6] * z + coeff[7] * y*z;
		grad[1] = coeff[2] + coeff[4] * x + coeff[5] * z + coeff[7] * z*x;
		grad[2] = coeff[3] + coeff[5] * y + coeff[6] * x + coeff[7] * x*y;
		return grad;
	}
	virtual PointType getIsoPoint() override{
		return PointType(rand2(), rand2(), rand2());
	}
};

class TestPrism6 :public TestBase<Coor3, dim3::BasisFunction> {

public:
	TestPrism6() :
		TestBase<Coor3, dim3::BasisFunction>(basis_func3d.getInstance(dim3::ElementType::Prism6)) {
		init();
	}
protected:
	void init() {
		NUM = 6;
		elem_iso_coor = {

			PointType(1, 0, 0),
			PointType(1, 1, 0),
			PointType(1, 0, 1),
			PointType(-1, 0, 0),
			PointType(-1, 1, 0),
			PointType(-1, 0, 1)
		};
		elem_phy_coor = {
			PointType(0.5 + 0.2*rand1(), 0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.5 + 0.2*rand1(), 0.2*rand1()),
			PointType(0.5 + 0.2*rand1(), 0.2*rand1(), 0.5 + 0.2*rand1()),
			PointType(-0.5 - 0.2*rand1(), 0.2*rand1(), 0.2*rand1()),
			PointType(-0.5 - 0.2*rand1(), 0.5 + 0.2*rand1(), 0.2*rand1()),
			PointType(-0.5 - 0.2*rand1(), 0.2*rand1(), 0.5 + 0.2*rand1())
		};
		for (size_t i = 0; i < NUM; ++i) {
			coeff.push_back(1.0*rand2());
		}
	}
	virtual double polyValue(PointType& coor) override {
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * z
			+ coeff[4] * x*y
			+ coeff[5] * x*z;
	}

	virtual std::vector<double> polyGrad(PointType& coor) override {
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();
		std::vector<double> grad(3);

		grad[0] = coeff[1] + coeff[4] * y + coeff[5] * z;
		grad[1] = coeff[2] + coeff[4] * x;
		grad[2] = coeff[3] + coeff[5] * x;
		return grad;
	}
	virtual PointType getIsoPoint() override {
		return PointType(rand2(), 0.5*rand1(), 0.5*rand1());
	}
};

class TestPyra5 :public TestBase<Coor3, dim3::BasisFunction> {

public:
	TestPyra5() :
		TestBase<Coor3, dim3::BasisFunction>(basis_func3d.getInstance(dim3::ElementType::Pyramid5)) {
		init();
	}
protected:
	void init() {
		NUM = 5;
		elem_iso_coor = {
			PointType(+1, 0, 0),
			PointType(0, +1, 0),
			PointType(-1, 0, 0),
			PointType(0, -1, 0),
			PointType(0, 0, +1)
		};
		elem_phy_coor = {
			PointType(10 + 2 * rand2(), 2 * rand2(), 2 * rand2()),
			PointType(2 * rand2(), 10 + 2 * rand2(), 2 * rand2()),
			PointType(-10 + 2 * rand2(), 2 * rand2(), 2 * rand2()),
			PointType(2 * rand2(), -10 + 2 * rand2(), 2 * rand2()),
			PointType(2 * rand2(), 2 * rand2(), 10 + 2 * rand2())
		};
		for (size_t i = 0; i < NUM; ++i) {
			coeff.push_back(1.0*rand2());
		}
	}
	virtual double polyValue(PointType& coor) override {
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * z;
	}

	virtual std::vector<double> polyGrad(PointType& coor) override {
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();
		std::vector<double> grad(3);

		grad[0] = coeff[1];
		grad[1] = coeff[2];
		grad[2] = coeff[3];
		return grad;
	}
	virtual PointType getIsoPoint() override {
		return PointType(0.5*rand2(), 0.5*rand2(), 0.5*rand2());
	}
};

class TestTET10 :public TestBase<Coor3, dim3::BasisFunction>{

public:
	TestTET10() :
		TestBase<Coor3, dim3::BasisFunction>(basis_func3d.getInstance(dim3::ElementType::Tetra10)){
		init();
	}
protected:
	void init() {
		NUM = 10;
		elem_iso_coor = {
			PointType(1, 0, 0),//1
			PointType(0, 1, 0),//2
			PointType(0, 0, 1),//3
			PointType(0, 0, 0),//4
			PointType(0.5, 0.5, 0),//5
			PointType(0, 0.5, 0.5),//6
			PointType(0.5, 0, 0.5),//7
			PointType(0.5, 0, 0),//8
			PointType(0, 0.5, 0),//9
			PointType(0, 0, 0.5)//10
		};
		double x1, x2, x3, x4;
		double y1, y2, y3, y4;
		double z1, z2, z3, z4;
		x1 = 0.5 + 0.2*rand1();
		y1 = 0.2*rand1();
		z1 = 0.2*rand1();
		x2 = 0.2*rand1();
		y2 = 0.5 + 0.2*rand1();
		z2 = 0.2*rand1();
		x3 = 0.2*rand1();
		y3 = 0.2*rand1();
		z3 = 0.5 + 0.2*rand1();
		x4 = 0.2*rand1();
		y4 = 0.2*rand1();
		z4 = 0.2*rand1();
		elem_phy_coor = {
			PointType(x1, y1, z1),
			PointType(x2, y2, z2),
			PointType(x3, y3, z3),
			PointType(x4, y4, z4),
			PointType(0.5*(x1 + x2), 0.5*(y1 + y2), 0.5*(z1 + z2)),
			PointType(0.5*(x3 + x2), 0.5*(y3 + y2), 0.5*(z3 + z2)),
			PointType(0.5*(x1 + x3), 0.5*(y1 + y3), 0.5*(z1 + z3)),
			PointType(0.5*(x1 + x4), 0.5*(y1 + y4), 0.5*(z1 + z4)),
			PointType(0.5*(x4 + x2), 0.5*(y4 + y2), 0.5*(z4 + z2)),
			PointType(0.5*(x3 + x4), 0.5*(y3 + y4), 0.5*(z3 + z4))
		};
		for (size_t i = 0; i < NUM; ++i){
			coeff.push_back(1.0);
		}
	}
	virtual double polyValue(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * z
			+ coeff[4] * x*x
			+ coeff[5] * y*y
			+ coeff[6] * z*z
			+ coeff[7] * x*y
			+ coeff[8] * y*z
			+ coeff[9] * x*z;
	}

	virtual std::vector<double> polyGrad(PointType& coor) override{

		std::vector<double> grad(3);
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();
		grad[0] = coeff[1] + 2.0 * coeff[4] * x + coeff[7] * y + coeff[9] * z;
		grad[1] = coeff[2] + 2.0 * coeff[5] * y + coeff[7] * x + coeff[8] * z;
		grad[2] = coeff[3] + 2.0 * coeff[6] * z + coeff[8] * y + coeff[9] * x;
		return grad;
	}
	virtual PointType getIsoPoint() override{
		auto a = PointType(0.5*rand1(), 0.5*rand1(), 0.5*rand1());
		return a;
	}
};

class TestHexa20 :public TestBase<Coor3, dim3::BasisFunction>{

public:
	TestHexa20() :
		TestBase<Coor3, dim3::BasisFunction>(basis_func3d.getInstance(dim3::ElementType::Hexa20)){
		init();
	}
protected:
	void init() {
		NUM = 20;

		elem_iso_coor = {
			PointType(-1, -1, -1),//1
			PointType(+1, -1, -1),//2
			PointType(+1, +1, -1),//3
			PointType(-1, +1, -1),//4
			PointType(-1, -1, +1),//5
			PointType(+1, -1, +1),//6
			PointType(+1, +1, +1),//7
			PointType(-1, +1, +1),//8

			PointType(0, -1, -1),//9
			PointType(1, 0, -1),//10
			PointType(0, 1, -1),//11
			PointType(-1, 0, -1),//12
			PointType(0, -1, 1),//13
			PointType(1, 0, 1),//14
			PointType(0, 1, 1),//15
			PointType(-1, 0, 1),//16
			PointType(-1, -1, 0),//17
			PointType(1, -1, 0),//18
			PointType(1, 1, 0),//19
			PointType(-1, 1, 0)//20

		};

		elem_phy_coor = {
			PointType(-10 - 2 * rand1(), -10 - 2 * rand1(), -10 - 2 * rand1()),
			PointType(+10 + 2 * rand1(), -10 - 2 * rand1(), -10 - 2 * rand1()),
			PointType(+10 + 2 * rand1(), +10 + 2 * rand1(), -10 - 2 * rand1()),
			PointType(-10 - 2 * rand1(), +10 + 2 * rand1(), -10 - 2 * rand1()),
			PointType(-10 - 2 * rand1(), -10 - 2 * rand1(), +10 + 2 * rand1()),
			PointType(+10 + 2 * rand1(), -10 - 2 * rand1(), +10 + 2 * rand1()),
			PointType(+10 + 2 * rand1(), +10 + 2 * rand1(), +10 + 2 * rand1()),
			PointType(-10 - 2 * rand1(), +10 + 2 * rand1(), +10 + 2 * rand1()),

			PointType(2 * rand1(), -10 - 2 * rand1(), -10 - 2 * rand1()),
			PointType(10 + 2 * rand1(), 2 * rand1(), -10 - 2 * rand1()),
			PointType(2 * rand1(), 10 + 2 * rand1(), -10 - 2 * rand1()),
			PointType(-10 - 2 * rand1(), 2 * rand1(), -10 - 2 * rand1()),
			PointType(2 * rand1(), -10 - 2 * rand1(), 10 + 2 * rand1()),
			PointType(10 + 2 * rand1(), 2 * rand1(), 10 + 2 * rand1()),
			PointType(2 * rand1(), 10 + 2 * rand1(), 10 + 2 * rand1()),
			PointType(-10 - 2 * rand1(), 2 * rand1(), 10 + 2 * rand1()),
			PointType(-10 - 2 * rand1(), -10 - 2 * rand1(), 2 * rand1()),
			PointType(10 + 2 * rand1(), -10 - 2 * rand1(), 2 * rand1()),
			PointType(10 + 2 * rand1(), 10 + 2 * rand1(), 2 * rand1()),
			PointType(-10 - 2 * rand1(), 10 + 2 * rand1(), 2 * rand1())
		};
		for (size_t i = 0; i < NUM; ++i){
			coeff.push_back(1.0*rand2());
		}
	}
	virtual double polyValue(PointType& coor) override{
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * z
			+ coeff[4] * x*y
			+ coeff[5] * y*z
			+ coeff[6] * x*z
			+ coeff[7] * x*x
			+ coeff[8] * y*y
			+ coeff[9] * z*z
			+ coeff[10] * x*x*y
			+ coeff[11] * x*y*y
			+ coeff[12] * x*x*z
			+ coeff[13] * z*z*x
			+ coeff[14] * y*y*z
			+ coeff[15] * y*z*z
			+ coeff[16] * x*y*z
			+ coeff[17] * x*x*y*z
			+ coeff[18] * x*y*y*z
			+ coeff[19] * x*y*z*z;
	}

	virtual std::vector<double> polyGrad(PointType& coor) override{

		std::vector<double> grad(3);
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();
		grad[0] = coeff[1]
			+ coeff[4] * y
			+ coeff[6] * z
			+ 2.0*coeff[7] * x
			+ 2.0*coeff[10] * x*y
			+ coeff[11] * y*y
			+ 2.0*coeff[12] * x*z
			+ coeff[13] * z*z
			+ coeff[16] * y*z
			+ 2.0*coeff[17] * x*y*z
			+ coeff[18] * y*y*z
			+ coeff[19] * y*z*z;
		grad[1] = coeff[2]
			+ coeff[4] * x
			+ coeff[5] * z
			+ 2.0*coeff[8] * y
			+ coeff[10] * x*x
			+ 2.0*coeff[11] * x*y
			+ 2.0*coeff[14] * y*z
			+ coeff[15] * z*z
			+ coeff[16] * x*z
			+ coeff[17] * x*x*z
			+ 2.0*coeff[18] * x*y*z
			+ coeff[19] * x*z*z;
		grad[2] = coeff[3]
			+ coeff[5] * y
			+ coeff[6] * x
			+ 2.0*coeff[9] * z
			+ coeff[12] * x*x
			+ 2.0*coeff[13] * x*z
			+ coeff[14] * y*y
			+ 2.0*coeff[15] * y*z
			+ coeff[16] * x*y
			+ coeff[17] * x*x*y
			+ coeff[18] * x*y*y
			+ 2.0*coeff[19] * x*y*z;
		return grad;
	}
	virtual PointType getIsoPoint() override{
		return PointType(rand2(), rand2(), rand2());
	}
};

class TestPrism15 :public TestBase<Coor3, dim3::BasisFunction> {

public:
	TestPrism15() :
		TestBase<Coor3, dim3::BasisFunction>(basis_func3d.getInstance(dim3::ElementType::Prism15)) {
		init();
	}
protected:
	void init() {
		NUM = 15;
		elem_iso_coor = {
			PointType(1, 0, 0),
			PointType(1, 1, 0),
			PointType(1, 0, 1),
			PointType(-1, 0, 0),
			PointType(-1, 1, 0),
			PointType(-1, 0, 1),

			PointType(1, 0.5, 0),
			PointType(1, 0.5, 0.5),
			PointType(1, 0, 0.5),

			PointType(-1, 0.5, 0),
			PointType(-1, 0.5, 0.5),
			PointType(-1, 0, 0.5),

			PointType(0, 0, 0),
			PointType(0, 1, 0),
			PointType(0, 0, 1)
		};
		elem_phy_coor = {
			PointType(10 + 2 * rand2(), 2 * rand2(), 2 * rand2()),
			PointType(10 + 2 * rand2(), 10 + 2 * rand2(), 2 * rand2()),
			PointType(10 + 2 * rand2(), 2 * rand2(), 10 + 2 * rand2()),
			PointType(-10 + 2 * rand2(), 2 * rand2(), 2 * rand2()),
			PointType(-10 + 2 * rand2(), 10 + 2 * rand2(), 2 * rand2()),
			PointType(-10 + 2 * rand2(), 2 * rand2(), 10 + 2 * rand2()),

			PointType(10 + 2 * rand2(), 5 + rand2(), 2 * rand2()),
			PointType(10 + 2 * rand2(), 5 + rand2(), 5 + rand2()),
			PointType(10 + 2 * rand2(), 2 * rand2(), 5 + rand2()),

			PointType(-10 + 2 * rand2(), 5 + rand2(), 2 * rand2()),
			PointType(-10 + 2 * rand2(), 5 + rand2(), 5 + rand2()),
			PointType(-10 + 2 * rand2(), 2 * rand2(), 5 + rand2()),

			PointType(2 * rand2(), 2 * rand2(), 2 * rand2()),
			PointType(2 * rand2(), 10 + 2 * rand2(), 2 * rand2()),
			PointType(2 * rand2(), 2 * rand2(), 10 + 2 * rand2())
		};
		for (size_t i = 0; i < NUM; ++i) {
			coeff.push_back(1.0);
		}
	}
	virtual double polyValue(PointType& coor) override {
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * z
			+ coeff[4] * x*x
			+ coeff[5] * y*y
			+ coeff[6] * z*z
			+ coeff[7] * x*y*z
			+ coeff[8] * x*z
			+ coeff[9] * x*y
			+ coeff[10] * y*z
			+ coeff[11] * x*x*y
			+ coeff[12] * x*x*z
			+ coeff[13] * x*z*z
			+ coeff[14] * x*y*y;
	}

	virtual std::vector<double> polyGrad(PointType& coor) override {
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();
		std::vector<double> grad(3);

		grad[0] = coeff[1]
			+ 2.0*coeff[4] * x
			+ coeff[7] * y*z
			+ coeff[8] * z
			+ coeff[9] * y
			+ 2.0*coeff[11] * x*y
			+ 2.0*coeff[12] * x*z
			+ coeff[13] * z*z
			+ coeff[14] * y*y;
		grad[1] = coeff[2]
			+ 2.0*coeff[5] * y
			+ coeff[7] * x*z
			+ coeff[9] * x
			+ coeff[10] * z
			+ coeff[11] * x*x
			+ 2.0*coeff[14] * x*y;
		grad[2] = coeff[3]
			+ 2.0*coeff[6] * z
			+ coeff[7] * x*y
			+ coeff[8] * x
			+ coeff[10] * y
			+ coeff[12] * x*x
			+ 2.0*coeff[13] * x*z;
		return grad;
	}
	virtual PointType getIsoPoint() override {
		return PointType(rand2(), 0.4*rand1(), 0.4*rand1());
	}
};

class TestPyra13 :public TestBase<Coor3, dim3::BasisFunction> {

public:
	TestPyra13() :
		TestBase<Coor3, dim3::BasisFunction>(basis_func3d.getInstance(dim3::ElementType::Pyramid13)) {
		init();
	}
protected:
	void init() {
		NUM = 13;
		elem_iso_coor = {
			PointType(1, 0, 0),
			PointType(0, 1, 0),
			PointType(-1, 0, 0),
			PointType(0, -1, 0),
			PointType(0, 0, 1),
			PointType(0.5, 0.5, 0),
			PointType(-0.5, 0.5, 0),
			PointType(-0.5, -0.5, 0),
			PointType(0.5, -0.5, 0),
			PointType(0.5, 0, 0.5),
			PointType(0, 0.5, 0.5),
			PointType(-0.5, 0, 0.5),
			PointType(0, -0.5, 0.5),
		};
		elem_phy_coor = {
			PointType(10 + 2 * rand2(), 2 * rand2(), 2 * rand2()),
			PointType(2 * rand2(), 10 + 2 * rand2(), 2 * rand2()),
			PointType(-10 + 2 * rand2(), 2 * rand2(), 2 * rand2()),
			PointType(2 * rand2(), -10 + 2 * rand2(), 2 * rand2()),
			PointType(2 * rand2(), 2 * rand2(), 10 + 2 * rand2()),
			PointType(5 + rand2(), 5 + rand2(), 2 * rand2()),
			PointType(-5 + rand2(), 5 + rand2(), 2 * rand2()),
			PointType(-5 + rand2(), -5 + rand2(), 2 * rand2()),
			PointType(5 + rand2(), -5 + rand2(), 2 * rand2()),
			PointType(5 + rand2(), 2 * rand2(), 5 + rand2()),
			PointType(2 * rand2(), 5 + rand2(), 5 + rand2()),
			PointType(-5 + rand2(), 2 * rand2(), 5 + rand2()),
			PointType(2 * rand2(), -5 + rand2(), 5 + rand2()),
		};
		for (size_t i = 0; i < NUM; ++i) {
			coeff.push_back(1.0*rand2());
		}
	}
	virtual double polyValue(PointType& coor) override {
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();

		return coeff[0]
			+ coeff[1] * x
			+ coeff[2] * y
			+ coeff[3] * z
			+ coeff[4] * x*x
			+ coeff[5] * y*y
			+ coeff[6] * z*z
			+ coeff[7] * x*y
			+ coeff[8] * y*z
			+ coeff[9] * x*z;
	}

	virtual std::vector<double> polyGrad(PointType& coor) override {
		double x = coor.getX();
		double y = coor.getY();
		double z = coor.getZ();
		std::vector<double> grad(3);

		grad[0] = coeff[1] + 2.0 * coeff[4] * x + coeff[7] * y + coeff[9] * z;
		grad[1] = coeff[2] + 2.0 * coeff[5] * y + coeff[7] * x + coeff[8] * z;
		grad[2] = coeff[3] + 2.0 * coeff[6] * z + coeff[8] * y + coeff[9] * x;

		return grad;
	}
	virtual PointType getIsoPoint() override {
		PointType a = PointType(0.5*rand2(), 0.5*rand2(), 0.5*rand2());
		return a;
	}
};

int main() {


	{std::cout << "TestQ4 \n";
	TestQ4 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestQ8 \n";
	TestQ8 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestQ9 \n";
	TestQ9 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestT3 \n";
	TestT3 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestTET4 \n";
	TestTET4 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestHEXA8 \n";
	TestHEXA8 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestPrism6 \n";
	TestPrism6 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestPyra5 \n";
	TestPyra5 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestTET10 \n";
	TestTET10 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestHexa20 \n";
	TestHexa20 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestPrism15 \n";
	TestPrism15 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }

	{std::cout << "TestPyra13 \n";
	TestPyra13 trial;
	trial.testDelta();
	trial.testUnity();
	trial.testIsoInterpolation();
	trial.testIsoGrad();
	trial.testJacobian();
	trial.testInvJ();
	trial.testdetJ(); }




	system("pause");
	return 0;
}
