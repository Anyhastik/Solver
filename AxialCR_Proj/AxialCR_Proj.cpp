#include "AxialCR.hpp"

/**
 * Класс содержит функции для описания осевого каналирования
 */
mp_t ChanRad::theta(mp_t x)
{
	if (x < 0)
		return 0;
	else
		return 1;
}

mp_t ChanRad::omega(int n, mp_t TT)
{
	return n * 2 * pi / TT;
}

mp_t ChanRad::pow2(mp_t x) {
	return pow(x, 2);
}

/**
* Расчет излучения
*/
mp_t ChanRad::dWdE(int Nt, mp_t w, mp_t TT, mp_t integral[nbeam])
{
	mp_t multiplier = pow(10, -2) * pow2(echarge) / (hbar * pow2(c0) * pow2(TT)) * Nt;
	mp_t sum = 0;
	for (size_t n = 1; n <= 5; n++)    //число гармоник
	{
		sum += theta(2 * pow2(gama) * omega(n, TT) - w) * w * (1 - 2 * w / (2 * pow2(gama) * omega(n, TT))
			+ 2 * pow(w / (2 * pow2(gama) * omega(n, TT)), 2)) * integral[n - 1];
	}
	sum *= multiplier;
	return sum;
}

/**
* Функция интегрирования методом трапеций
*/
mp_t ChanRad::trapezoid(vector<mp_t> times, vector<mp_t> zdots, mp_t TT)
{
	mp_t tmp1 = 0;
	//mp_t tmp2 = 0;
	int i = 1;
	while (times.at(i) <= TT)
	{
		tmp1 += (zdots.at(i) / (c0) + zdots.at(i - 1) / (c0)) / 2 * (times.at(i) - times.at(i - 1));
		i++;
		if (i == times.size())
		{
			throw std::invalid_argument ("error in calculations of TT.");
			break;
		}
	}
	return tmp1;
}

/**
* Функция интегрирования методом трапеций
*/
mp_t ChanRad::trapezoid2(vector<mp_t> times, betta betas, int n, mp_t TT)
{
	mp_t omega1 = omega(n, TT); //число гармоник = 5 n=[0,5]
	mp_t tmp1 = 0;
	mp_t tmp2 = 0;
	mp_t tmp3 = 0;
	mp_t tmp4 = 0;
	int i = 1;
	while (times.at(i) <= TT)
	{
		tmp1 += (-betas.x.at(i - 1) * cos(omega1 * times[i - 1]) + betas.x.at(i) * cos(omega1 * times[i]))
			 / 2 * (times[i] - times[i - 1]) / (c0);
		tmp3 += (-betas.y.at(i - 1) * cos(omega1 * times[i - 1]) + betas.y.at(i) * cos(omega1 * times[i]))
			 / 2 * (times[i] - times[i - 1]) / (c0);
		tmp2 += (betas.x.at(i - 1) * sin(omega1 * times[i - 1]) + betas.x.at(i) * sin(omega1 * times[i]))
			/ 2 * (times[i] - times[i - 1]) / (c0);
		tmp4 += (betas.y.at(i - 1) * sin(omega1 * times[i - 1]) + betas.y.at(i) * sin(omega1 * times[i]))
			/ 2 * (times[i] - times[i - 1]) / (c0);
		i++;
	}
	return pow2(tmp1) + pow2(tmp3);
}

void ChanRad::print(array<mp_t, 3>gg)
{
	cout << "gg = " << gg[0] << " " << gg[1] << " " << gg[2] << endl;
}

/**
* Содержит коэффициенты для расчета формы потенциала
* в зависимости от типа кристалла crystaltype
* При этом сейчас здесь вычисляется dp - межплоскостное расстояние
* Векторы a1 и b1 нужно будет убрать
*/
void ChanRad::cpotss(string crystaltype)
{
	mp_t a0;
	mp_t p, EMinus4 = 1e-4, EMinus8 = 1e-8;
	array<mp_t, 5> a1;
	array<mp_t, 5> b1;

	if (crystaltype == "si")
	{
		a0 = 5.431 * EMinus4;
		p = 0.075 * EMinus4;
		dp = a0 / 2.0;
		a1 = { 2.1186 * EMinus4, 2.4960 * EMinus4, 0.8104 * EMinus4, 0.3365 * EMinus4, 0.0567 * EMinus4 };
		b1 = { 57.6767 * EMinus8, 16.7929 * EMinus8, 3.2522 * EMinus8, 0.6155 * EMinus8, 0.0582 * EMinus8 };
		//return { a0, p, a1, b1 };
	}

	if (crystaltype == "si111")
	{
		a0 = 4 * 1.18 * EMinus4;
		p = 0.075 * EMinus4;
		dp = a0 / 2.0;
		a1 = { 2.1186 * EMinus4, 2.4960 * EMinus4, 0.8104 * EMinus4, 0.3365 * EMinus4, 0.0567 * EMinus4 };
		b1 = { 57.6767 * EMinus8, 16.7929 * EMinus8, 3.2522 * EMinus8, 0.6155 * EMinus8, 0.0582 * EMinus8 };
		//return { a0, p, a1, b1 };
	}

	else if (crystaltype == "c")
	{
		a0 = 3.5668 * EMinus4;
		p = 0.04 * EMinus4;
		dp = a0 / 2.0;
		a1 = { 0.3555 * EMinus4, 1.1420 * EMinus4, 0.7537 * EMinus4, 0.2091 * EMinus4, 0.0489 * EMinus4 };
		b1 = { 51.1341 * EMinus8, 17.8811 * EMinus8, 5.4281 * EMinus8, 1.0825 * EMinus8, 0.1140 * EMinus8 };
		//return { a0, p, a1, b1 };
	}
	else if (crystaltype == "ge")
	{
		a0 = 5.658 * EMinus4;
		p = 0.085 * EMinus4;
		dp = a0 / 2.0;
		a1 = { 1.6356 * EMinus4, 2.8938 * EMinus4, 1.6555 * EMinus4, 0.9761 * EMinus4, 0.2135 * EMinus4 };
		b1 = { 70.3903 * pow(10, -20), 21.5563 * pow(10, -20), 4.5527 * pow(10, -20), 0.9845 * pow(10, -20), 0.0989 * pow(10, -20) };
		//return { a0, p, a1, b1 };
	}
	else
	{
		a0 = 5.431 * EMinus4;
		p = 0.075 * EMinus4;
		dp = a0 / 2.0;
		a1 = { 2.1186 * EMinus4, 2.4960 * EMinus4, 0.8104 * EMinus4, 0.3365 * EMinus4, 0.0567 * EMinus4 };
		b1 = { 57.6767 * EMinus8, 16.7929 * EMinus8, 3.2522 * EMinus8, 0.6155 * EMinus8, 0.0582 * EMinus8 };
		//return { a0, p, a1, b1 };
	}
}

/**
* Векторное произведение векторов
*/
array<mp_t, 3> ChanRad::cross_vector(array<mp_t, 3> a, array<mp_t, 3> b)
{
	array<mp_t, 3> vec;
	vec[0] = (a[1] * b[2] - a[2] * b[1]);
	vec[1] = (a[2] * b[0] - a[0] * b[2]);
	vec[2] = (a[0] * b[1] - a[1] * b[0]);
	return vec;
}

/**
* Произведение вектора на число
*/
array<mp_t, 3> ChanRad::multyply_vec_num(mp_t num, array<mp_t, 3> vec)
{
	for (int i = 0; i < 3; i++)
		vec[i] *= num;
	return vec;
}

/**
* Скалярное произведение векторов
*/
mp_t ChanRad::scalar_mult_vec(array<mp_t, 3> g1, array<mp_t, 3> g2)
{
	if ((g1[0] == 0 && g1[1] == 0 && g1[2] == 0) || (g2[0] == 0 && g2[1] == 0 && g2[2] == 0))
		return 0;
	else
		return g1[0] * g2[0] + g1[1] * g2[1] + g1[2] * g2[2];
}

// here will be new approximation
mp_t ChanRad::periodicPotntial1(mp_t x, mp_t y) {
	/**
	* Функция для описания ямы с распределением Коши по ценру (четные ряды)
	*/
	return ampl + (a / (1 + pow2((sin(2 * pi * (x) / T) - muX) / sigmaX) +
		pow2((sin(2 * pi * (y) / T) - muY) / sigmaY)));
}

mp_t ChanRad::periodicPotntial2(mp_t x, mp_t y) {
	/**
	* Функция для описания ямы с распределением Коши для нечетных рядов
	*/
	return ampl + (a / (1 + pow2((sin(2 * pi * (x - T / 4) / T) - muX) / sigmaX) +
		pow2((sin(2 * pi * (y + T / 4) / T) - muY) / sigmaY))) + T / 2;
}

mp_t ChanRad::crystal(mp_t x, mp_t y) {
	/**
	* Результирующая функция для описания ямы с распределением Коши
	*/
	return - periodicPotntial1(x, y) - periodicPotntial2(x, y);
}

mp_t ChanRad::DXcrystal(mp_t x, mp_t y) {
	/**
	* Производная результирующей функция для описания ямы с распределением Коши по x
	*/
	const mp_type mp =
		derivative(mp_type(mp_type(x)),
			mp_type(mp_type(1) / 100000000U), 
			[&](const mp_type& x) -> mp_type
			{
				return crystal(static_cast<mp_t>(x), y) ;  /// Function
			});
	mp_t res = mp.convert_to<mp_t>(); /// Convert to closest double.
	return res;
}

mp_t ChanRad::DYcrystal(mp_t x, mp_t y) {
	/**
	* Производная результирующей функция для описания ямы с распределением Коши по y
	*/
	const mp_type mp =
		derivative(mp_type(mp_type(y)),
			mp_type(mp_type(1) / 100000000U), 
			[&](const mp_type& y) -> mp_type
			{
				return crystal(x, static_cast<mp_t>(y)) ;  /// Function
			});
	mp_t res = mp.convert_to<mp_t>(); /// Convert to closest double.
	return res;
}


/** 
* Функция, описывающая правую часть системы ODE
*/
void ChanRad::system_function(const state_type& y, state_type& dydt, mp_t t) {
	dydt[0] = y[2];
	dydt[1] = y[3];
	dydt[2] = DXcrystal(y[0], y[1]) * secondDeriv;
	dydt[3] = DYcrystal(y[0], y[1]) * secondDeriv;
}

/** 
* Функция для сохранения данных
*/
void ChanRad::observer(const state_type& y, mp_t t) {
	x_values.push_back(y[0]);
	y_values.push_back(y[1]);
	vx_values.push_back(y[2]);
	vy_values.push_back(y[3]);
	ax_values.push_back(DXcrystal(y[0], y[1]) * secondDeriv);
	ay_values.push_back(DYcrystal(y[0], y[1]) * secondDeriv);
	t_values.push_back(t);
}

/**
* Функция решает систему диффенциальных уравнений
*/
void ChanRad::solveOde(size_t i, size_t j) {

	t0 = crystaltickness / c0;
	mp_t h = t0 / n_steps;  /// Шаг интегрирования
	mp_t t_start = 0.0;
	mp_t t_end = t0;
	VectorXd y0 = { 1.0 / (nbeam + 1.0) * (i + 1.0) * dp / 2.0,
						1.0 / (nbeam + 1.0) * (j + 1.0) * dp / 2.0,
						mom * pow(10, 6) * c0 * sin(tetta) / (m0 * gama),
						mom * pow(10, 6) * c0 * sin(tetta) / (m0 * gama) };  /// Начальные условия [x, y, x', y']
	// Выбор интегратора (RK54)
	typedef vector<mp_t> state_type;
	auto stepper = make_controlled<runge_kutta_cash_karp54<state_type>>(1e-6, 1e-6);
	// Вычисление решения с сохранением данных в массивы
	integrate_adaptive(stepper, std::bind(&ChanRad::system_function, this,
		std::placeholders::_1, std::placeholders::_2, std::placeholders::_3), 
		y0, t_start, t_end, h, std::bind(&ChanRad::observer, this, std::placeholders::_1, std::placeholders::_2));
	SizeTrajectory = t_values.size();
	/*cout << "----------------------------------------" << endl;
	cout << "t\tx\ty\tx'\ty'" << endl;
	for (size_t i = 0; i < t_values.size(); ++i) {
		cout << t_values[i] << "\t"
			<< x_values[i] << "\t"
			<< y_values[i] << "\t"
			<< vx_values[i] << "\t"
			<< vy_values[i] << "\t"
			<< ax_values[i] << "\t"
			<< ay_values[i] << endl;*/
	//}
}

void ChanRad::CleanArrays() {
	x_values = {};
	y_values = {};
	vx_values = {};
	vy_values = {};
	ax_values = {};
	ay_values = {};
	t_values = {};
}

void ChanRad::findMinX(const std::vector<double>& x) {
	minElements.push_back(*std::min_element(x.begin(), x.end()));
}

void ChanRad::findMaxX(const std::vector<double>& x) {
	maxElements.push_back(*std::max_element(x.begin(), x.end()));
}

mp_t ChanRad::cosTest(mp_t x){
	return cos(2 * x);
}

mp_t ChanRad::DcosTest(mp_t y) {
	/**
	* Производная test
	*/
	const mp_type mp =
		derivative(mp_type(mp_type(y)),
			mp_type(mp_type(1) / 100000000U), /// Step size 10^-8.
			[&](const mp_type& y) -> mp_type
			{
				return cosTest(static_cast<mp_t>(y));  /// Function
			});
	ld1 res = mp.convert_to<ld1>(); /// Convert to closest double.
	return res;
}



