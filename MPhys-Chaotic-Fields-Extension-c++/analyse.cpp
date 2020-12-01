#include"analyse.h"

# define M_PI 3.14159265358979323846


template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
	T h = (b - a) / static_cast<T>(N - 1);
	std::vector<T> xs(N);
	typename std::vector<T>::iterator x;
	T val;
	for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
		*x = val;
	return xs;
}




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

std::array<std::vector<double>, 2> double_abc_poincare(
	double t_start, 
	double t_end, 
	double t_step, 
	std::array<double, 3> x_ini, 
	std::vector<double> params) {

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



double line_distance_once(const Solution& firstline, const Solution& secondline, const int& index) {

	double ret{ sqrt(
				pow(firstline.x[index] - secondline.x[index], 2) +
				pow(firstline.y[index] - secondline.y[index], 2) +
				pow(firstline.z[index] - secondline.z[index], 2)) };

	return ret;

}


std::vector<double> line_distance(const Solution& line1, const Solution& line2) {

	std::vector<double> to_ret;
	for (int i{}; i < line1.s.size(); i++) {
		to_ret.push_back(log(line_distance_once(line1, line2, i)));
	}

	return to_ret;

}




std::array<std::vector<double>, 2> line_variance(
	double n, 
	double end, 
	double step, 
	std::array<double, 3> ini,
	std::vector<double> params, 
	std::string method) {

	int vol{ static_cast<int>(pow(n, 3)) };
	std::vector<std::array<double, 3>> origin;
	for (int i{ }; i < n; i++) {
		for (int j{ }; j < n; j++) {
			for (int k{ }; k < n; k++) {

				std::array<double, 3> wow{ ini[0] + -0.02 + (static_cast<double>(i) * 0.01),
					ini[1] + -0.02 + (static_cast<double>(j) * 0.01),
					ini[2] + -0.02 + (static_cast<double>(k) * 0.01) };
				origin.push_back(wow);

			}
		}
	}

	int mid = (vol + 1) / 2;

	double start{ 0 };

	std::vector<Solution> lines;  // for holding field line solutions
	if (method == "abc") { for (auto it{ origin.begin() }; it < origin.end(); it++) { lines.push_back(abc_field(start, end, step, *it, params)); } }
	else if (method == "double") { for (auto it{ origin.begin() }; it < origin.end(); it++) { lines.push_back(double_abc_field(start, end, step, *it, params)); } }
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


std::array<std::vector<double>, 3> trajectory_split(
	double end, 
	double step, 
	std::array<double, 3> ini, 
	std::vector<double> params,
	double z0) {

	Solution refline{ abc_field(0, end, step, ini, params) };
	Solution divline{ abc_field(0, end, step, std::array<double, 3>{ini[0] + z0, ini[1], ini[2]}, params) };

	std::vector<double> distance_apart;
	std::vector<double> distance_apart_log;
	double dis{};
	for (int i{}; i < refline.s.size(); i++) {

		dis = line_distance_once(refline, divline, i);
		distance_apart_log.push_back(log(dis));
		distance_apart.push_back(dis);

	}

	return std::array<std::vector<double>, 3>{refline.s, distance_apart, distance_apart_log};

}


std::array<std::vector<double>, 2> lyapunov(
	double end, 
	double step, 
	std::array<double, 3> ini, 
	std::vector<double> params) {

	double z0 = 0.0000001;

	double lyapunov_distance{ 60 };
	int n_points{ 500 };
	int index{ static_cast<int>(end / (step * static_cast<double>(n_points))) };
	
	std::vector<double> starts{ linspace(0., end, n_points) };
	std::vector<double> exponent;

	Solution refline;
	Solution divline;

	int count{};
	while (count < n_points) {

		refline = abc_field(0, lyapunov_distance, step, ini, params);
		divline = abc_field(0, lyapunov_distance, step, std::array<double, 3>{ini[0] + z0, ini[1], ini[2]}, params);
		ini = { refline.x[index], refline.y[index], refline.z[index] };

		exponent.push_back(log(line_distance_once(refline, divline, refline.s.size() - 1) / z0) / lyapunov_distance);

		count++;

	}
	
	return std::array<std::vector<double>, 2>{starts, exponent};

}



std::vector<double> recurrence_single(
	double end, 
	double step, 
	std::array<double, 3> ini, 
	std::vector<double> params, 
	std::array<double, 3> spherepos,
	double sphererad) {

	std::vector<double> rec;
	std::vector<double> s;
	bool inside{ false };
	std::array<double, 6> minmax{
		spherepos[0] + sphererad,
		spherepos[0] - sphererad,
		spherepos[1] + sphererad,
		spherepos[1] - sphererad,
		spherepos[2] + sphererad,
		spherepos[2] - sphererad };

	Solution line{ abc_field(0, end, step, ini, params) };
	std::vector<double> x{ periodic_projection(line.x) };
	std::vector<double> y{ periodic_projection(line.y) };
	std::vector<double> z{ periodic_projection(line.z) };
	
	for (int i{}; i < line.s.size(); i++) {

		if (x[i] < minmax[0] && x[i] > minmax[1] &&
			y[i] < minmax[2] && y[i] > minmax[3] &&
			z[i] < minmax[4] && z[i] > minmax[5]) {

			if (!inside) {

				s.push_back(line.s[i]);
				inside = true;
				continue;

			}
			else {

				continue;

			}

		}

		inside = false;

	}

	for (int i{1}; i < s.size(); i++) {

		rec.push_back(s[i] - s[i - 1]);

	}

	return rec;




}

std::vector<double> recurrence(
	double end, 
	double step, 
	std::array<double, 3> ini, 
	std::vector<double> params, 
	std::array<double, 3> spherepos, 
	double sphererad, 
	int n) {

	std::vector<double> x{ linspace(-0.01 + ini[0], 0.01 + ini[0], n) };
	std::vector<double> y{ linspace(-0.01 + ini[1], 0.01 + ini[1], n) };
	std::vector<double> z{ linspace(-0.01 + ini[2], 0.01 + ini[2], n) };
	std::vector < std::array<double, 3>> beginnings;
	
	for (auto i{ x.begin() }; i < x.end(); i++) {
		for (auto j{ y.begin() }; j < y.end(); j++) {
			for (auto k{ z.begin() }; k < z.end(); k++) {

				beginnings.push_back(std::array<double, 3>{*i, *j, *k});

			}
		}
	}

	std::vector<double> rec;
	for (auto i{ beginnings.begin() }; i < beginnings.end(); i++) {

		std::vector<double> temp{ recurrence_single(end, step, *i, params, spherepos, sphererad) };
		rec.insert(rec.end(), temp.begin(), temp.end());

	}

	return rec;

}




double test(double num) {

	return pow(num, 2);

}