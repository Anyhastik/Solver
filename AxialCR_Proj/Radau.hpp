//WORKING IMPLICIT RADAU IIA



#ifndef RADAU_METHOD
#define RADAU_METHOD

#include "Types.hpp"

const mp_t secondDeriv = -pow(c0, 2) / (m0 * gama);

namespace ublas = boost::numeric::ublas;

using VectorXd = std::vector<mp_t>;


/**
* Оператор для умножения числа на вектор
*/
static VectorXd operator*(mp_t scalar, const VectorXd& vec) noexcept {
    VectorXd result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = scalar * vec[i];
    }
    return result;
}

/**
* Оператор для умножения вектора на число
*/
static VectorXd operator*(const VectorXd& vec, mp_t scalar) noexcept {
    return scalar * vec; // просто обратный порядок аргументов
}

/**
* Оператор для сложения векторов
*/
static VectorXd operator+(const VectorXd& vec1, const VectorXd& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Vectors must have the same size for addition.");
    }
    VectorXd result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

class Radau {
public:
    Radau(){}
    virtual ~Radau(){}
    /**
    * Функции правой части системы уравнений
    */
    VectorXd f(mp_t t, const VectorXd& y);
    VectorXd fun(mp_t t, const VectorXd& y);
    /**
    *  Метод Радау IIА пятого порядка
    */
    void radauIIA(mp_t h, mp_t t_start, mp_t t_end, const VectorXd& y0, VectorXd& times,
        std::vector<std::vector<mp_t>>& answer);
    void radauIIA2(mp_t h, mp_t t_start, mp_t t_end, const VectorXd& y0, VectorXd& times,
        std::vector<std::vector<mp_t>>& answer);

    template<typename T>
    bool solve_linear_system(ublas::matrix<T>& A, ublas::vector<T>& b);

    void rhs(const ublas::vector<mp_t>& u, std::vector<mp_t>& dudt);
    void lockeWartingSolver(std::vector<mp_t>& u, mp_t t0, mp_t t1, mp_t h, std::vector<mp_t>& times,
        std::vector<std::vector<mp_t>>& answer);

};
#endif //RADAU_METHOD