#ifndef TYPES_H
#define TYPES_H

# include <limits>
# include <iostream>
# include <iomanip>
# include <type_traits>

#include <cmath>
#include <array>
#include <vector>
#include <fstream>

#include <boost/numeric/odeint.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/interpolators/pchip.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace std;
using namespace std::chrono;
using namespace boost::numeric::odeint;
using namespace boost::multiprecision;
using boost::math::interpolators::pchip;
using boost::multiprecision::number;
using boost::multiprecision::cpp_dec_float;

// Re-compute using 5 extra decimal digits precision (22) than double (17).
#define MP_DIGITS10 unsigned (std::numeric_limits<double>::max_digits10 + 5)

typedef cpp_dec_float<MP_DIGITS10> mp_backend;
typedef number<mp_backend> mp_type;

//typedef cpp_dec_float<MP_DIGITS10> mp_backend1;
//typedef number<mp_backend1> mp_t;

#define ld1 long double
#define mp_t double
//#define mp_t cpp_dec_float_50

// Константы в старой системе единиц

//#define c0 299792458
//#define echarge 4.8032 * pow(10, -10) * sqrt(1 / (1.60217733 * pow(10, -6)) * pow(10, -2))
//#define hbar 6.582 * pow(10, -22)        //(*MeV s*)
//#define m0 0.5109996 * pow(10, 6)        //(*ev/c^2*)
//#define e0 0.5109996                    //(*MeV*)
//#define e 800-e0                          // MeV INPUT!!!!!
//#define mom sqrt(pow(e, 2) - pow(e0, 2)) //(*MeV*)
//#define pi 3.14159265358979323846      // GEANT
//#define gama  (e + e0)/ e0
//#define coeff  1.123443973421022 * pow(10, 28) // coeff = - pow(c0, 2) * pow(10, 20) / (gama * m0);
//// #define dp 1.35775 //3.13559
//#define tetta 0
//#define nbeam 5
//#define n_steps 1000
//#define abor 0.529 * pow(10, -10)
//#define angl pow(10, -10)


#define angl (pow(10, -6))
#define c0 (299792.4580)   //299792458 * 1e-3
#define hbar ((6.582) * (pow(10, -13)) )
#define abor ((0.529) * (pow(10, -10)) / (angl))
#define echarge (sqrt(hbar * c0 / 137))
#define m0 ((0.5109996) * (pow(10, 6)) ) 
#define e0 (0.5109996)                   
#define e ((800) - (e0))
#define gama (((e) + (e0))/ (e0))
#define mom (sqrt((pow(e, 2)) - (pow(e0, 2))))
#define tetta (0)
#define nbeam (5)
#define n_steps (10000)
#define pi (3.14159265358979323846) 
#define crystaltickness (20)

#endif //TYPES_H