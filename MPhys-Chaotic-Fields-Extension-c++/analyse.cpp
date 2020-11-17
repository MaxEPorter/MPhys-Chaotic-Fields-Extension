#include"analyse.h"

# define M_PI 3.14159265358979323846

std::array<std::vector<double>, 2> abc_poincare(
	double t_start,
	double t_end,
	double t_step,
	std::array<double, 3> x_ini,
	std::vector<double> params,
	std::string plane,
	double planevalue)
{

	Solution sol{ abc_field(t_start, t_end, t_step, x_ini, params) };
	std::array<std::vector<double>, 3> line{ sol.x, sol.y, sol.z };

	int i, j, k;
	if (plane == "z") { i = 2; j = 0; k = 1; }
	else if (plane == "x") { i = 0; j = 1; k = 2; }
	else if (plane == "y") { i = 1; j = 0; k = 2; }

	std::vector<double> xp;
	std::vector<double> yp;

	for (int m{ 1 }; m < line[0].size(); m++) {  // every point along s

		double x{ line[i][m - 1] / (2 * M_PI) - planevalue }; // get previous point and divide by 2pi, then shift by plane value
		double y{ line[i][m] / (2 * M_PI) - planevalue };  // get current point

		double num1; // integer part
		double num2;
		double rem1; // remainder part
		double rem2;
		rem1 = std::modf(x, &num1);  // remainder after dividing by 2pi
		rem2 = std::modf(y, &num2);
		//if (rem1 < 0) { rem1 += 1; }
		//if (rem2 < 0) { rem2 += 1; }


		if (num1 == 0 && num2 == 0) { // when integer parts are both 0
			if (rem1 < 0 && rem2>0 || rem1 > 0 && rem2 < 0) { // crossing the integer boundary both ways
				xp.push_back(line[j][m]);
				yp.push_back(line[k][m]);
				continue;
			}
			continue;
		}
		if (num1<num2 || num1>num2) { // crossing integer boundary either way
			xp.push_back(line[j][m]);
			yp.push_back(line[k][m]);
			continue;
		}

	}

	return std::array<std::vector<double>, 2>{periodic_projection(xp), periodic_projection(yp)};

}

std::array<std::vector<double>, 2> double_abc_poincare(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params) {

	std::array<std::vector<double>, 2> to_return;
	state_type prev{ x_ini[0], x_ini[1], x_ini[2] };

	fields::DoubleABC abc(params);  // initialise abc field
	boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;  // initialise rkf stepper/method
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };  // initialise initial starting positions

	boost::numeric::odeint::integrate_const(stepper,  // pass rkf stepper
		[&abc](const state_type& x, state_type& dxdt, double s) {abc.ode(x, dxdt, s); },  // give the function to integrate as abc.ode
		x,  // give the x, y, z positions
		t_start,  // starting s
		t_end,  // ending s
		t_step,  // s step
		// lambda to save field trajectory data into solution
		[&to_return, &prev](const state_type& x, const double t) {
			if (prev[0] < 1272 && x[0] > 1272) { to_return[0].push_back(x[0]); to_return[1].push_back(x[1]); } prev = x; });

	return to_return;  // return solution

}

double periodic_projection_once(double x) {

	x = x / (2 * M_PI);

	double w;
	x = std::modf(x, &w);
	if (x < 0) { x += 1; }

	return x;

}
std::vector<double> periodic_projection(std::vector<double> x) {

	std::vector<double> ret;

	for (auto i{ x.begin() }; i < x.end(); i++) {

		ret.push_back(periodic_projection_once(*i));

	}

	return ret;

}

std::array<std::vector<double>, 2> line_variance(double n, double end, double step, std::vector<double> params, std::string method) {

	int vol{ static_cast<int>(pow(n, 3)) };
	std::vector<std::array<double, 3>> ini;
	for (int i{ }; i < n; i++) {
		for (int j{ }; j < n; j++) {
			for (int k{ }; k < n; k++) {

				std::array<double, 3> wow{ -0.02 + (static_cast<double>(i) * 0.01),
					-0.02 + (static_cast<double>(j) * 0.01),
					-0.02 + (static_cast<double>(k) * 0.01) };
				ini.push_back(wow);

			}
		}
	}

	int mid = (vol + 1) / 2;

	double start{ 0 };

	std::vector<Solution> lines;  // for holding field line solutions
	if (method == "abc") { for (auto it{ ini.begin() }; it < ini.end(); it++) { lines.push_back(abc_field(start, end, step, *it, params)); } }
	else if (method == "double") { for (auto it{ ini.begin() }; it < ini.end(); it++) { lines.push_back(double_abc_field(start, end, step, *it, params)); } }
	else { return std::array<std::vector<double>, 2>{std::vector<double>(1), std::vector<double>(1)}; } // if no field is called

	std::vector<double> variance;
	for (int i{}; i < lines[0].s.size(); i++) { // loop over entire path length

		// calculate mu b
		double mub{ 0 };
		for (int j{ }; j < vol; j++) {

			mub += sqrt(
				pow(lines[j].x[i] - lines[mid].x[i], 2) +
				pow(lines[j].y[i] - lines[mid].y[i], 2) +
				pow(lines[j].z[i] - lines[mid].z[i], 2));

		}
		mub = mub / vol;

		// calculate variance
		double v{ 0 };
		for (int j{ }; j < vol; j++) {
			v += pow(
				sqrt(
					pow(lines[j].x[i] - lines[mid].x[i], 2) +
					pow(lines[j].y[i] - lines[mid].y[i], 2) +
					pow(lines[j].z[i] - lines[mid].z[i], 2)) -
				mub, 2);
		}
		v = v / vol;

		variance.push_back(v);

	}

	return std::array<std::vector<double>, 2>{lines[0].s, variance};

}

std::array<std::vector<double>, 2> lyapunov(double end, double step, std::array<double, 3> ini, std::vector<double> params) {

	double z0 = 0.0000001;

	double lyapunov_distance{ 50 };
	int lyapunov_steps{ static_cast<int>(lyapunov_distance / step) };

	int max_steps{ static_cast<int>(end / step) };
	int n_points{ max_steps / lyapunov_steps };

	std::vector<int> points;
	for (int i{}; i < n_points; i++) {
		points.push_back(lyapunov_steps * i);
	}

	Solution line1{ abc_field(0, end, step, ini, params) };

	std::vector<double> s;
	std::vector<double> exponent;

	for (int i{ 1 }; i < points.size(); i++) {

		std::vector<std::array<double, 3>> starts{
			std::array<double, 3>{line1.x[points[i - 1]] + z0, line1.y[points[i - 1]], line1.z[points[i - 1]]},
			std::array<double, 3>{line1.x[points[i - 1]] - z0, line1.y[points[i - 1]], line1.z[points[i - 1]]},
			std::array<double, 3>{line1.x[points[i - 1]], line1.y[points[i - 1]] + z0, line1.z[points[i - 1]]},
			std::array<double, 3>{line1.x[points[i - 1]], line1.y[points[i - 1]] - z0, line1.z[points[i - 1]]},
			std::array<double, 3>{line1.x[points[i - 1]], line1.y[points[i - 1]], line1.z[points[i - 1]] + z0},
			std::array<double, 3>{line1.x[points[i - 1]], line1.y[points[i - 1]], line1.z[points[i - 1]] - z0},
		};

		std::vector<Solution> lines;
		for (auto it{ starts.begin() }; it < starts.end(); it++) {
			lines.push_back(abc_field(0, lyapunov_distance, step, *it, params));
		}

		std::vector<double> zts;
		for (auto it{ lines.begin() }; it < lines.end(); it++) {
			zts.push_back(sqrt(pow(line1.x[points[i] - 1] - it->x.back(), 2) +
				pow(line1.y[points[i] - 1] - it->y.back(), 2) +
				pow(line1.z[points[i] - 1] - it->z.back(), 2)));
		}

		//double av{};
		//for (auto it{ zts.begin() }; it < zts.end(); it++) {
		//	av += *it;
		//}
		//av = av / 6;
		double biggest{};
		for (int i{}; i < 6; i++) {
			if (zts[i] > biggest) { biggest = zts[i]; }
		}


		s.push_back(line1.s[points[i]]);

		exponent.push_back(log(biggest / z0) / lyapunov_distance);

	}

	return std::array<std::vector<double>, 2>{s, exponent};

}