
#ifndef AXIAL_CR
#define AXIAL_CR

#include "Radau.hpp"
#include "Types.hpp"

/**
 * Класс содержит функции для описания осевого каналирования
 */
class ChanRad {
public:
	const mp_t secondDeriv = -pow(c0, 2) / (m0 * gama);

	const mp_t ampl = 1.23;  //mp_t ampl = 0;
	const mp_t muX = 0;
	const mp_t muY = 0;
	const mp_t sigmaX = 0.21;
	const mp_t sigmaY = 0.21;
	const mp_t a = 88;
	const mp_t T = 2.713 * 2 * 1e-4;
	mp_t dp;
	mp_t t0;
	mp_t SizeTrajectory;
	vector<mp_t> x_values;
	vector<mp_t> y_values;
	vector<mp_t> vx_values;
	vector<mp_t> vy_values;
	vector<mp_t> ax_values;
	vector<mp_t> ay_values;
	vector<mp_t> t_values;

	VectorXd maxElements;
	VectorXd minElements;

	ChanRad() {
		x_values = {};
		y_values = {};
		vx_values = {};
		vy_values = {};
		ax_values = {};
		ay_values = {};
		t_values = {};
		maxElements = {};
		minElements = {};
		dp = 0;
		t0 = 1e-16;
		SizeTrajectory = 0;
	}

	virtual ~ChanRad() {}
	template <typename value_type, typename function_type>
	value_type derivative(const value_type x, const value_type dx, function_type function)
	{
		/*! \brief Compute the derivative of function using a 3-point central difference rule of O(dx^6).
		  \tparam value_type, floating-point type, for example: `double` or `cpp_dec_float_50`
		  \tparam function_type

		  \param x Value at which to evaluate derivative.
		  \param dx Incremental step-size.
		  \param function Function whose derivative is to computed.

		  \return derivative at x.
		*/

		static_assert(false == std::numeric_limits<value_type>::is_integer, "value_type must be a floating-point type!");

		const value_type dx2(dx * 2U);
		const value_type dx3(dx * 3U);
		// Difference terms.
		const value_type m1((function(x + dx) - function(x - dx)) / 2U);
		const value_type m2((function(x + dx2) - function(x - dx2)) / 4U);
		const value_type m3((function(x + dx3) - function(x - dx3)) / 6U);
		const value_type fifteen_m1(m1 * 15U);
		const value_type six_m2(m2 * 6U);
		const value_type ten_dx(dx * 10U);
		return ((fifteen_m1 - six_m2) + m3) / ten_dx;  // Derivative.
	}
	static mp_t pow2(mp_t x);

	mp_t theta(mp_t x);
	mp_t omega(int n, mp_t TT);
	mp_t dWdE(int Nt, mp_t w, mp_t TT, mp_t integral[nbeam]);

	struct betta {
		vector <mp_t>  x;
		vector <mp_t>  y;
		vector <mp_t>  z;
	};

	mp_t trapezoid(vector<mp_t> times, vector<mp_t> zdots, mp_t TT);

	mp_t trapezoid2(vector<mp_t> times, betta betas, int n, mp_t TT);
	
	void print(array<mp_t, 3>gg);

	void cpotss(string crystaltype);

	array<mp_t, 3> cross_vector(array<mp_t, 3> a, array<mp_t, 3> b);

	array<mp_t, 3> multyply_vec_num(mp_t num, array<mp_t, 3> vec);
	
	mp_t scalar_mult_vec(array<mp_t, 3> g1, array<mp_t, 3> g2);

	typedef vector<mp_t> state_type;
	void system_function(const state_type& y, state_type& dydt, mp_t t);
	void observer(const state_type& y, mp_t t);
	
	// Функции для вычисления потенциала и решения системы дифференциальных уравнений
	mp_t periodicPotntial1(mp_t x, mp_t y) ;
	mp_t periodicPotntial2(mp_t x, mp_t y);
	mp_t crystal(mp_t x, mp_t y);
	mp_t DXcrystal(mp_t x, mp_t y);
	mp_t DYcrystal(mp_t x, mp_t y);
	mp_t cosTest(mp_t x);

	void solveOde(size_t i, size_t j);

	void CleanArrays();
	void findMinX(const std::vector<mp_t>& x);
	void findMaxX(const std::vector<mp_t>& x);

	mp_t DcosTest(mp_t y);

	struct radiation
	{
		vector<mp_t> w;
		vector<mp_t> rad;
	};
};

#endif //AXIAL_CR

