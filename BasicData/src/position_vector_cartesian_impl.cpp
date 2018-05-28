//
// Created by Chingkai Chou on 5/27/18.
//

template <size_t Dimension, bool Authority, typename T> struct CartesianAPI_FuncHelper;

template<typename T> struct CartesianAPI_FuncHelper<1, false, T> :
        public virtual T{
    virtual const double &getX() const { return this->data->at(0); }
};
template<typename T> struct CartesianAPI_FuncHelper<1, true, T> :
        public virtual CartesianAPI_FuncHelper<1, false, T>{
    virtual void setX(const double &x){ this->data->at(0) = x; }
    virtual void setX(double &x){ this->data->at(0) = x; }
};
template<typename T> struct CartesianAPI_FuncHelper<2, false, T> :
        public virtual CartesianAPI_FuncHelper<1, false, T>{
    virtual const double &getY() const { return this->data->at(1); }
};
template<typename T> struct CartesianAPI_FuncHelper<2, true, T> :
        public virtual CartesianAPI_FuncHelper<1, true, T>,
        public virtual CartesianAPI_FuncHelper<2, false, T> {
    virtual void setY(const double &y){ this->data->at(1) = y; }
    virtual void setY(double &y){ this->data->at(1) = y; }
};
template<typename T> struct CartesianAPI_FuncHelper<3, false, T> :
        public virtual CartesianAPI_FuncHelper<2, false, T>{
    virtual const double &getZ() const { return this->data->at(2); }
};
template<typename T> struct CartesianAPI_FuncHelper<3, true, T> :
        public virtual CartesianAPI_FuncHelper<2, true, T>,
        public virtual CartesianAPI_FuncHelper<3, false, T> {
    virtual void setZ(const double &z){ this->data->at(2) = z; }
    virtual void setZ(double &z){ this->data->at(2) = z; }
};

template <size_t Dimension>
struct CartesianAPI<Dimension, false>:
        public virtual CartesianAPI_FuncHelper<Dimension, false, PointData<Dimension>> {

CartesianAPI<Dimension, false>& operator=(const CartesianAPI<Dimension, false>& other)
{
    PointData<Dimension>::operator=(other);
    return *this;
}
};

template <size_t Dimension>
struct CartesianAPI<Dimension, true>:
        public virtual CartesianAPI_FuncHelper<Dimension, true, PointData<Dimension>>{

CartesianAPI<Dimension, true>& operator=(const CartesianAPI<Dimension, true>& other)
{
    PointData<Dimension>::operator=(other);
    return *this;
}

};

