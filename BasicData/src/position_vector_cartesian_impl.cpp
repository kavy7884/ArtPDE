//
// Created by Chingkai Chou on 5/27/18.
//

template <size_t Dimension, bool Authority> class CartesianAPI;
template <size_t Dimension, bool Authority, typename T> struct CartesianAPI_FuncHelper;
template <size_t Dimension, bool Authority, typename T> struct CartesianAPI_FuncWrap;

// -------- CartesianAPI <Start> -----------
template <size_t Dimension>
struct CartesianAPI<Dimension, true>:
        public CartesianAPI_FuncWrap<Dimension, true, PointData<Dimension>>{

    CartesianAPI<Dimension, true> &operator=(const CartesianAPI<Dimension, true> &other) {
        PointData<Dimension>::operator=(other);
        return *this;
    }

};

template <size_t Dimension>
struct CartesianAPI<Dimension, false>:
        public CartesianAPI_FuncWrap<Dimension, false, PointData<Dimension>>{

    CartesianAPI<Dimension, false>& operator=(const CartesianAPI<Dimension, false>& other)
    {
        PointData<Dimension>::operator=(other);
        return *this;
    }

    const double &getDataById(const size_t id) const{
        return PointData<Dimension>::getDataById(id);
    }
};

// -------- CartesianAPI <End> -----------


// -------- CartesianAPI_FuncHelper <Start> -----------
template<typename T> struct CartesianAPI_FuncHelper<1, false, T> :
        public virtual T{
    const double &getX() const { return this->data->at(0); }
};
template<typename T> struct CartesianAPI_FuncHelper<1, true, T> :
        public virtual T{
    void setX(const double &x){ this->data->at(0) = x; }
    void setX(double &x){ this->data->at(0) = x; }
};
template<typename T> struct CartesianAPI_FuncHelper<2, false, T> :
        public virtual T{
    const double &getY() const { return this->data->at(1); }
};
template<typename T> struct CartesianAPI_FuncHelper<2, true, T> :
        public virtual T{
    void setY(const double &y){ this->data->at(1) = y; }
    void setY(double &y){ this->data->at(1) = y; }
};
template<typename T> struct CartesianAPI_FuncHelper<3, false, T> :
        public virtual T{
    const double &getZ() const { return this->data->at(2); }
};
template<typename T> struct CartesianAPI_FuncHelper<3, true, T> :
        public virtual T{
    void setZ(const double &z){ this->data->at(2) = z; }
    void setZ(double &z){ this->data->at(2) = z; }
};
// -------- CartesianAPI_FuncHelper <End> -----------

// -------- CartesianAPI_FuncWrap <Start> -----------
template <typename T>
struct CartesianAPI_FuncWrap<1, false, T>:
        public CartesianAPI_FuncHelper<1, false, T>{};

template <typename T>
struct CartesianAPI_FuncWrap<1, true, T>:
        public CartesianAPI_FuncHelper<1, false, T>,
        public CartesianAPI_FuncHelper<1, true, T>{};

template <typename T>
struct CartesianAPI_FuncWrap<2, false, T>:
        public CartesianAPI_FuncHelper<1, false, T>,
        public CartesianAPI_FuncHelper<2, false, T>{};

template <typename T>
struct CartesianAPI_FuncWrap<2, true, T>:
        public CartesianAPI_FuncHelper<1, false, T>,
        public CartesianAPI_FuncHelper<2, false, T>,
        public CartesianAPI_FuncHelper<1, true, T>,
        public CartesianAPI_FuncHelper<2, true, T>{};

template <typename T>
struct CartesianAPI_FuncWrap<3, false, T>:
        public CartesianAPI_FuncHelper<1, false, T>,
        public CartesianAPI_FuncHelper<2, false, T>,
        public CartesianAPI_FuncHelper<3, false, T>{};

template <typename T>
struct CartesianAPI_FuncWrap<3, true, T>:
        public CartesianAPI_FuncHelper<1, false, T>,
        public CartesianAPI_FuncHelper<2, false, T>,
        public CartesianAPI_FuncHelper<3, false, T>,
        public CartesianAPI_FuncHelper<1, true, T>,
        public CartesianAPI_FuncHelper<2, true, T>,
        public CartesianAPI_FuncHelper<3, true, T>{};
// -------- CartesianAPI_FuncWrap <End> -----------