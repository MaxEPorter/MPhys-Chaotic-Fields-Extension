#pragma once

#include"fields.h"

#include <array>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <numeric>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#include <boost/array.hpp>
#include <boost/range/combine.hpp>

// pack data to return to python as solution class
class Solution {

public:

	std::vector<double> s;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;

	std::string field;
	std::vector<double> params;
	std::string method;
	std::array<double, 3> start;

	double step;
	double time_taken;

};

Solution abc_field(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params);
Solution abc_field_euler(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params);
Solution abc_field_rk4(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params);

Solution wire_field(double t_start, double t_end, double t_step, std::array<double, 3> x_ini);
Solution wire_field_euler(double t_start, double t_end, double t_step, std::array<double, 3> x_ini);
Solution wire_field_rk4(double t_start, double t_end, double t_step, std::array<double, 3> x_ini);

Solution double_abc_field(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params);

Solution uniform_field(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params);
Solution uniform_field_euler(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params);
Solution uniform_field_rk4(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params);

Solution exp_field(double t_start, double t_end, double t_step, std::array<double, 3> x_ini);

