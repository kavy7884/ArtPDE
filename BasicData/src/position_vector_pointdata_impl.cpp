//
// Created by Chingkai Chou on 5/27/18.
//

template <size_t Dimension>
void PointData<Dimension>::newData(){
    this->data = nullptr;
    this->data = std::make_shared<ArrayType>();
    std::cout << "PointData new data" << std::endl;
}

template <size_t Dimension>
void PointData<Dimension>::addDataByList(const std::initializer_list<double> &v){
    assert(v.size() <= Dimension);
    std::copy(v.begin(), v.end(), data->begin());
}

template <size_t Dimension>
PointData<Dimension>& PointData<Dimension>::operator=(const PointData<Dimension>& other)
{
    std::cout << "copy" << std::endl;
    if(&other == this)
        return *this;
    if(this->data == nullptr)
        this->newData();
    std::copy(other.data->begin(), other.data->end(), this->data->begin());
    return *this;
}

template <size_t Dimension>
std::ostream &operator<<(std::ostream &os, const PointData<Dimension> &point_data) {
    os << "[ ";
    for (size_t i = 0; i < point_data.data->size() - 1; ++i) {
        os << point_data.data->at(i) << " ";
    }
    os << point_data.data->at(point_data.data->size() - 1);
    os << " ] ";
    return os;
}
