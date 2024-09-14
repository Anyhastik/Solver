/**
* This file contains the main function.
*/

#include "Types.hpp"
#include "Radau.hpp"
#include "AxialCR.hpp"


int main()
{
	ChanRad* chr = new ChanRad;
	string crystaltype = "si";  ///"si111";
	
	chr->cpotss(crystaltype);
	
	ChanRad::betta beta;
	ChanRad::radiation emission;

	/*vector <mp_t> t;
	vector <mp_t> x;
	vector <mp_t> vx;
	vector <mp_t> ax;
	vector <mp_t> y;
	vector <mp_t> vy;
	vector <mp_t> ay;*/

	/**
	* Массив сплайнов (интерполяций) функций
	*/
	std::vector<pchip<std::vector<mp_t, std::allocator<mp_t>>>> splines;
	// evaluate at a point:
	//double z = spline2(3.4);
	mp_t Tc = 0, betaave = 1;
	Radau* solver = new Radau;
	/**
	* Циклы по точкам входа
	*/
	for (size_t point_i = 0; point_i < nbeam; point_i++)
	{
		for (size_t point_j = 0; point_j < nbeam; point_j++)
		{
			emission.w.clear();
			emission.rad.clear();
			/*t.clear();
			x.clear();
			vx.clear();
			ax.clear();
			y.clear();
			vy.clear();
			ay.clear();*/
			/**
			* Решение системы уравнений методом RK54
			*/
			chr->solveOde(point_i, point_j);

			/**
			* Решение системы уравнений методом Радау IIV
			*/
			/// Задаем начальные параметры
			//mp_t h = t0 / n_steps;   //0.001; /// Шаг интегрирования
			//mp_t t_start = 0.0;
			//mp_t t_end = t0;
			//VectorXd y0 = { 1.0 / (nbeam + 1.0) * (point_i + 1.0) * dp / 2.0,
			//				1.0 / (nbeam + 1.0) * (point_j + 1.0) * dp / 2.0,
			//				mom * pow(10, 6) * c0 * sin(tetta) / (m0 * gama),
			//				mom * pow(10, 6) * c0 * sin(tetta) / (m0 * gama) };  /// Начальные условия [x, y, x', y']
			//cout << y0[0] << " " << y0[1] << " " << y0[2] << " " << y0[3] << endl;
			//VectorXd times;
			//times.clear();
			//std::vector<std::vector<mp_t>> answer;
			//answer.clear();
			////// Запуск метода Радау IIА
			//solver->radauIIA(h, t_start, t_end, y0, times, answer);

			/**
			* Решение системы уравнений методом lockeWartingSolver
			*/
			/// Начальные условия: x(0) = 1.0, y(0) = 2.0, vx(0) = 2.3, vy(0) = 3.0
			//std::vector<mp_t> u = { 1.0 / (nbeam + 1.0) * (point_i + 1.0) * dp / 2.0,
			//					1.0 / (nbeam + 1.0) * (point_j + 1.0) * dp / 2.0,
			//					mom * pow(10, 6) * c0 * sin(tetta) / (m0 * gama),
			//					mom * pow(10, 6) * c0 * sin(tetta) / (m0 * gama) };  /// Начальные условия [x, y, x', y']
			//mp_t t00 = 0.0, t1 = t0, h = t0 / n_steps;  // Шаг интегрирования может быть увеличен при необходимости
			//vector<mp_t> times;
			//times.clear();
			//std::vector<std::vector<mp_t>> answer;
			//answer.clear();
			//solver->lockeWartingSolver(u, t00, t1, h, times, answer);
			//for (int i = 0; i < n_steps; i++)
			//{
			//	//fx << times.at(i) << " " << answer[i][0] << endl;
			//	x.push_back(answer[i][0]);
			//	//fvx << times.at(i) << " " << answer[i][2] << endl;
			//	vx.push_back(answer[i][2]);
			//	//fy << times.at(i) << " " << answer[i][1] << endl;
			//	y.push_back(answer[i][1]);
			//	//fvy << times.at(i) << " " << answer[i][3] << endl;
			//	vy.push_back(answer[i][3]);
			//	//fay << times.at(i) << " " << chr->DYcrystal(x.back(), y.back()) * secondDeriv << endl;
			//	ay.push_back(chr->DYcrystal(x.back(), y.back()) * secondDeriv);
			//	//fax << times.at(i) << " " << chr->DXcrystal(x.back(), y.back()) * secondDeriv << endl;
			//	ax.push_back(chr->DXcrystal(x.back(), y.back()) * secondDeriv);
			//}
			////надо убрать
			//fx.close();
			//fvx.close();
			//fax.close();
			//fy.close();
			//fvy.close();
			//fay.close();
			

			beta.x.clear();
			beta.y.clear();
			beta.z.clear();
			for (size_t i = 0; i < chr->SizeTrajectory; i++)
			{
				beta.x.push_back(chr->vx_values.at(i) / (c0));
				beta.y.push_back(chr->vy_values.at(i) / (c0));
				beta.z.push_back(0);
			}
			vector <mp_t> timez;
			timez.clear();
			ofstream timeArray("C:/ВУЗ/Магистратура/DIPLOM/ResultTxt/ResBoost/time.txt");
			for (size_t i = 0; i < chr->SizeTrajectory; i++)
			{
				timez.push_back( (2.0 * chr->vx_values.at(i) * chr->ax_values.at(i) + 2.0 * chr->vy_values.at(i) *
					chr->ay_values.at(i)) / pow(c0, 3) );
				timeArray << chr->t_values.at(i) << " " << timez.back() << endl;
			}
			timeArray.close();

			int k = 0;
			mp_t tp[2] = {};
			for (size_t i = 0; i < chr->SizeTrajectory - 1; i++)
			{
				if (timez.at(i) * timez.at(i + 1) < 0)
				{
					tp[k] = (chr->t_values.at(i) + chr->t_values.at(i + 1)) / 2.0;
					k++;
					if (k == 2)
						break;
				}
			}
			mp_t vz0 = c0 * sqrt(1 - 1.0 / chr->pow2(gama)) * cos(tetta);
			vector <mp_t> zdot;
			zdot.clear();
			for (size_t i = 0; i < chr->SizeTrajectory; i++)
			{
				zdot.push_back(vz0 * (1 - 1.0 / (2.0 * chr->pow2(c0)) * (chr->pow2(chr->vx_values.at(i)) +
					chr->pow2(chr->vy_values.at(i)) - chr->pow2(chr->vx_values.at(0)) - chr->pow2(chr->vy_values.at(0)))));
			}
			Tc = 4 * (tp[1] - tp[0]);
			betaave = chr->trapezoid(chr->t_values, zdot, Tc) / Tc;

			int Nt;
			Nt = int(chr->t0 / Tc);
			mp_t integrGxGy[5];
			/**
			* Здесь вычисляется сумма двух интегралов от gx и gy
			* Они используются в расчете излучения dWdE
			* Надо бы проверить насколько правильно сумма квадратов в функции написана
			*/
			for (size_t j = 0; j < 5; j++) //по гармоникам
			{
				integrGxGy[j] = chr->trapezoid2(chr->t_values, beta, j + 1, Tc);
				cout << "int = " << integrGxGy[j] << endl;
			}
			/**
			* Построение таблицы излучения
			*/
			/**
			* Расчитвается rad как массив массивов точек. То есть для каждой омеги у нас свой массив точек
			* Потом каждый из (nbeam^2) них мы отдельно интерполируем и в нужных точках все эти (nbeam^2) проинтерполированные
			* функции складываем по формуле и получаем излучение!!!!
			*/
			mp_t max_w = chr->omega(5, Tc) / (1.0 - betaave);
			mp_t min_w = chr->omega(1, Tc) / (1.0 + betaave);
			mp_t rad_counter_w = 0.01 * (max_w - min_w);

			while (min_w <= max_w)
			{
				emission.w.push_back(min_w * hbar);
				emission.rad.push_back(chr->dWdE(Nt, min_w, Tc, integrGxGy) / (nbeam * min_w * hbar * 10.0));
				//cout<<"back="<<emission.rad.back()<<endl;

				min_w += rad_counter_w;
			}
			chr->findMinX(emission.w);
			chr->findMaxX(emission.w);

			cout << point_i << " " << point_j << endl;
			auto spline = pchip<std::vector<mp_t, std::allocator<mp_t>>>(std::move(emission.w), std::move(emission.rad));
			splines.push_back(spline);
			chr->CleanArrays();
			
		}
	}
	ChanRad::radiation totalRadiation;
	totalRadiation.w.clear();
	totalRadiation.rad.clear();
	mp_t countW = 0, maxW = 0.5 * hbar * chr->omega(5, Tc) / (1.0 - betaave);
	mp_t step_counter = 0.005 * 0.5 * hbar * chr->omega(5, Tc) / (1.0 - betaave);

	ofstream fout89;
	fout89.open("C:/ВУЗ/Магистратура/DIPLOM/ResultTxt/ResBoost/radiation.txt");

	while (countW < maxW)
	{
		mp_t sumRad = 0.0;
		for (size_t ind = 0; ind < nbeam * nbeam; ind++)
		{
			if ((countW < chr->minElements.at(ind)) || (countW > chr->maxElements.at(ind))) {
				sumRad += 0.0;
			} 
			else { 
				sumRad += splines[ind](countW); 
			}
		}
		totalRadiation.w.push_back(countW);
		totalRadiation.rad.push_back(sumRad / nbeam);
		fout89 << countW << " " << totalRadiation.rad.back() << "\n";
		countW += step_counter;
	}

	fout89.close();
	delete chr;
	delete solver;

	return 0;
}
