#include "Point.hpp"

class Interpolation;

template<class NumericalMethodUtility, class Dimension>
class FunctionSpace{

public:

	template<class Dimension>
	using GeometryDataPtr = std::shared_ptr<GeometryData<Dimension>>;

	using InterpolationPtr = std::shared_ptr<Interpolation>;

	FunctionSpace() = delete;

	FunctionSpace(const GeometryDataPtr& geo_data_, const InterpolationPtr& method_) :
		geo_data(geo_data_), method(method_)
	{
		initialize();
	}

	FunctionSpace(const GeometryDataPtr& geo_data_, const InterpolationPtr& method_, const Reorder& reorder) :
		geo_data(geo_data_), method(method_)
	{
		initialize();
		reorder(solution_points);
	}

	double evaluate_shape(const PointList& position){
		return method_.evaluate_shape(position);
	}

	bool evaluate(Dof& interpolant, Dof& data, const PointList& old_position, const PointList& new_position){
		method_.evaluate(data, old_position, new_position);
		return true;
	}

private:

	bool initialize(){
		// new pointlist here
		return true;
	}

	GeometryDataPtr geo_data;
	InterpolationPtr method;
	std::shared_ptr<PointList<Dimension>> solution_points;

};