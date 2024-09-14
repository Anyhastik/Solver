#include "Radau.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include "AxialCR.hpp"

/**
* Функции правой части системы уравнений
*/
VectorXd Radau::f(mp_t t, const VectorXd& y) {
    VectorXd result(4);
    ChanRad* chr = new ChanRad;
    result[0] = y[2];
    result[1] = y[3];
    result[2] = chr->DXcrystal(y[0], y[1]) * secondDeriv;   //cos(y[1] * y[0]) + y[0];
    result[3] = chr->DYcrystal(y[0], y[1]) * secondDeriv;   //cos(y[1] * y[0]) - y[1];
    //result[4] = ;
    //result[5] = ;
    delete chr;
    return result;
}

VectorXd Radau::fun(mp_t t, const VectorXd& y) {
    VectorXd result(4);
    ChanRad* chr = new ChanRad;
    result[0] = y[2];
    result[1] = y[3];
    result[2] = cos(y[1] + y[0]) + y[0];   //cos(y[1] * y[0]) + y[0];
    result[3] = sin(y[1] + y[0]) - y[1];   //cos(y[1] * y[0]) - y[1];
    //result[4] = ;
    //result[5] = ;
    delete chr;
    return result;
}

/**
*  Метод Радау IIА пятого порядка
*/
void Radau::radauIIA(mp_t h, mp_t t_start, mp_t t_end, const VectorXd& y0, VectorXd& times,
    std::vector<std::vector<mp_t>>& answer) {
    auto start = chrono::steady_clock::now();
    int n = y0.size(); /// размерность системы

    /// Коэффициенты метода Радау IIА
    mp_t alpha = (6.0 + sqrt(6.0)) / 10.0;
    mp_t beta = 1.0 / (4.0 * alpha * (4.0 * alpha - 1.0));

    mp_t t = t_start;
    VectorXd y = y0;
    //VectorXd axay(6);

    while (t < t_end) {
        mp_t h_actual = (h < t_end - t) ? h : (t_end - t);  //std::min(h, t_end - t);
        //cout << "h_actual = " << h_actual << endl;
        VectorXd k1 = f(t, y);
        VectorXd k2 = f(t + alpha * h_actual, y + alpha * h_actual * k1);
        VectorXd k3 = f(t + 0.5 * h_actual, y + h_actual * ((1.0 - alpha) * k1 + alpha * k2));
        for (int i = 0; i < n; ++i) {
            y.at(i) += h_actual * (beta * (k1.at(i) + k2.at(i)) + (1.0 - 2.0 * beta) * k3.at(i)); // + b[3] * k4[i]);
            //axay.at(i) = y.at(i);
        }
        answer.push_back(y);
        times.push_back(t);
        t += h_actual;

        // Вывод результатов
       /* std::cout << "t = " << t << ", y = [";
        for (int i = 0; i < n; ++i) {
            std::cout << y[i];
            if (i < n - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;*/
    }
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in seconds: "
        << chrono::duration_cast<chrono::seconds>(end - start).count()
        << " sec";
}

void Radau::radauIIA2(mp_t h, mp_t t_start, mp_t t_end, const VectorXd& y0, VectorXd& times,
    std::vector<std::vector<mp_t>>& answer) {

    int n = y0.size(); /// размерность системы

    /// Коэффициенты метода Радау IIА
    mp_t alpha = (6.0 + sqrt(6.0)) / 10.0;
    mp_t beta = 1.0 / (4.0 * alpha * (4.0 * alpha - 1.0));

    mp_t t = t_start;
    VectorXd y = y0;
    //VectorXd axay(6);

    while (t < t_end) {
        mp_t h_actual = (h < t_end - t) ? h : (t_end - t);  //std::min(h, t_end - t);
        //cout << "h_actual = " << h_actual << endl;
        VectorXd k1 = fun(t, y);
        VectorXd k2 = fun(t + alpha * h_actual, y + alpha * h_actual * k1);
        VectorXd k3 = fun(t + 0.5 * h_actual, y + h_actual * ((1.0 - alpha) * k1 + alpha * k2));
        for (int i = 0; i < n; ++i) {
            y.at(i) += h_actual * (beta * (k1.at(i) + k2.at(i)) + (1.0 - 2.0 * beta) * k3.at(i)); // + b[3] * k4[i]);
            //axay.at(i) = y.at(i);
        }
        answer.push_back(y);
        times.push_back(t);
        t += h_actual;

        // Вывод результатов
        /*std::cout << "t = " << t << ", y = [";
        for (int i = 0; i < n; ++i) {
            std::cout << y[i];
            if (i < n - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;*/
    }
}

// Функция для решения системы линейных уравнений методом LU-разложения
template<typename T>
bool Radau::solve_linear_system(ublas::matrix<T>& A, ublas::vector<T>& b) {
    ublas::permutation_matrix<std::size_t> pm(A.size1());
    int res = lu_factorize(A, pm);
    if (res != 0) return false;
    lu_substitute(A, pm, b);
    return true;
}

// Определение правой части системы уравнений
void Radau::rhs(const ublas::vector<mp_t>& u, std::vector<mp_t>& dudt) {
    mp_t x = u[0];
    mp_t y = u[1];
    mp_t vx = u[2];
    mp_t vy = u[3];
    ChanRad* chr = new ChanRad;

    dudt[0] = vx;
    dudt[1] = vy;
    dudt[2] = chr->DXcrystal(x, y) * secondDeriv;  // std::cos(x + y) + x;
    dudt[3] = chr->DYcrystal(x, y) * secondDeriv;  // std::cos(x + y) - y;
}

// Метод Локка-Вартинга для неявной интеграции
void Radau::lockeWartingSolver(std::vector<mp_t>& u, mp_t t0, mp_t t1, mp_t h, std::vector<mp_t>& times,
    std::vector<std::vector<mp_t>>& answer) {
    const int n = u.size(); // Число переменных (4)
    const int s = 2;        // Число стадий для метода Локка-Вартинга

    ublas::matrix<mp_t> A(n, n);
    ublas::vector<mp_t> F(n), U(n), deltaU(n), U_temp(n);
    ublas::matrix<mp_t> Jacobian(n, n);

    while (t0 < t1) {
        // Инициализация U и Jacobian
        for (int i = 0; i < n; ++i) {
            U(i) = u[i];
            U_temp(i) = u[i] + h;
        }

        // Метод Ньютона для решения нелинейной системы
        for (int iter = 0; iter < 10; ++iter) { // Количество итераций может быть увеличено
            std::vector<mp_t> dudt(n);
            rhs(U_temp, dudt);
            for (int i = 0; i < n; ++i) {
                F(i) = U_temp(i) - u[i] - h * dudt[i];
            }

            // Построение Якобиана
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    Jacobian(i, j) = (i == j) ? 1.0 : 0.0;
                    if (i == j) {
                        Jacobian(i, j) -= h * (1.0 + std::cos(static_cast<ld1>(u[0] + u[1])));
                    }
                }
            }

            deltaU = F;
            if (!solve_linear_system(Jacobian, deltaU)) {
                std::cerr << "Ошибка: не удалось решить систему линейных уравнений" << std::endl;
                return;
            }

            // Обновление U_temp
            for (int i = 0; i < n; ++i) {
                U_temp(i) -= deltaU(i);
            }

            // Проверка сходимости
            if (ublas::norm_2(deltaU) < 1e-10) break; // Уменьшен порог сходимости
        }

        // Обновление значений u на следующем шаге
        u = std::vector<mp_t>(U_temp.data().begin(), U_temp.data().end());
        answer.push_back(u);
        times.push_back(t0);
        t0 += h;

        // Отладка: вывод текущих значений
       /* std::cout << "t = " << t0 << ", x = " << u[0] << ", y = " << u[1]
            << ", vx = " << u[2] << ", vy = " << u[3] << std::endl;*/
    }
}
