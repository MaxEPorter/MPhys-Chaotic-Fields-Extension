#pragma once

#include<array>
#include<vector>
#include<string>

#include"fields.h"
#include"solvers.h"

std::array<std::vector<double>, 2> abc_poincare(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params, std::string plane, double planevalue);
std::array<std::vector<double>, 2> double_abc_poincare(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params);

double periodic_projection_once(double x);
std::vector<double> periodic_projection(std::vector<double> x);

std::array<std::vector<double>, 2> line_variance(double n, double end, double step, std::vector<double> params, std::string method);

std::array<std::vector<double>, 2> lyapunov(double end, double step, std::array<double, 3> ini, std::vector<double> params);