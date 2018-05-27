//
// Created by Chingkai Chou on 5/26/18.
//

template<typename T>
class CartesianFuncBase: public virtual T{
public:
    template <typename ...Args>
    CartesianFuncBase(Args&& ...args) : T(std::forward<Args>(args)...){}
};

template <unsigned Dimension, bool Authority, typename T> class CartesianFunc;
template<typename T>
class CartesianFunc<1, false,T> : public virtual T{
public:
    const double &getX() const { return this->data->at(0); }
};

template<typename T>
class CartesianFunc<1, true,T> : public virtual T{
public:
    void setX(const double &x){ this->data->at(0) = x; }
    void setX(double &x){ this->data->at(0) = x; }
};

template <unsigned Dimension, bool Authority, typename T> class CartesianFunc;
template<typename T>
class CartesianFunc<2, false,T> : public virtual T{
public:
    const double &getY() const { return this->data->at(1); }
};

template<typename T>
class CartesianFunc<2, true,T> : public virtual T{
public:
    void setY(const double &y){ this->data->at(1) = y; }
    void setY(double &y){ this->data->at(1) = y; }
};

template <unsigned Dimension, bool Authority, typename T> class CartesianFunc;
template<typename T>
class CartesianFunc<3, false,T> : public virtual T{
public:
    const double &getZ() const { return this->data->at(2); }
};

template<typename T>
class CartesianFunc<3, true,T> : public virtual T{
public:
    void setZ(const double &z){ this->data->at(2) = z; }
    void setZ(double &z){ this->data->at(2) = z; }
};

template <unsigned Dimension, bool Authority, typename T> class CartesianWrap;
template<bool Authority, typename T>
class CartesianWrap<1, Authority, T>
        :public virtual CartesianFuncBase<T>,
         public virtual CartesianFunc<1,Authority, T>{
public:
    template <typename ...Args>
    CartesianWrap<1, Authority, T>(Args&& ...args) : CartesianFuncBase<T>(std::forward<Args>(args)...){}
};

template<bool Authority, typename T>
class CartesianWrap<2, Authority, T>
        :public virtual CartesianFuncBase<T>,
         public virtual CartesianFunc<2,Authority, T>,
         public virtual CartesianFunc<1,Authority, T>{
public:
    template <typename ...Args>
    CartesianWrap<2, Authority, T>(Args&& ...args) : CartesianFuncBase<T>(std::forward<Args>(args)...){}
};

template<bool Authority, typename T>
class CartesianWrap<3, Authority, T>
        :public virtual CartesianFuncBase<T>,
         public virtual CartesianFunc<3,Authority, T>,
         public virtual CartesianFunc<2,Authority, T>,
         public virtual CartesianFunc<1,Authority, T>{
public:
    template <typename ...Args>
    CartesianWrap<3, Authority, T>(Args&& ...args) : CartesianFuncBase<T>(std::forward<Args>(args)...){}
};

template <unsigned Dimension, typename T>
class CartesianReadable : public CartesianWrap<Dimension, false, T>{
public:
    template <typename ...Args>
    CartesianReadable(Args&& ...args) : CartesianWrap<Dimension, false, T>(std::forward<Args>(args)...){}

    template <typename ...Args>
    CartesianReadable(std::initializer_list<double> input_list, Args&& ...args) : CartesianWrap<Dimension, false, T>(std::forward<Args>(args)...){
        this->addDataByList(input_list); }
};

template <unsigned Dimension, typename T>
class CartesianWritable : public CartesianWrap<Dimension, true, T>, public CartesianWrap<Dimension, false, T>{
public:
    template <typename ...Args>
    CartesianWritable(Args&& ...args): CartesianWrap<Dimension, true, T>(std::forward<Args>(args)...) {}

    template <typename ...Args>
    CartesianWritable(std::initializer_list<double> input_list, Args&& ...args): CartesianWrap<Dimension, true, T>(std::forward<Args>(args)...){
        this->addDataByList(input_list); }
};
