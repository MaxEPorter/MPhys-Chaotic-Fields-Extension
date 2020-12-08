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

double line_distance_once(const Solution& firstline, const Solution& secondline, const int& index);
std::vector<double> line_distance(const Solution& line1, const Solution& line2);

std::array<std::vector<double>, 2> line_variance(double n, double end, double step, std::array<double, 3> ini, std::vector<double> params, std::string method);

std::array<std::vector<double>, 3> trajectory_split(double end, double step, std::array<double, 3> ini, std::vector<double> params, double z0);
std::array<std::vector<double>, 2> lyapunov(double end, double step, std::array<double, 3> ini, std::vector<double> params);

std::vector<double> recurrence_single(std::string method, double end, double step, std::array<double, 3> ini, std::vector<double> params, std::array<double, 3> spherepos, double sphererad);
std::vector<double> recurrence(std::string method, double end, double step, std::array<double, 3> ini, std::vector<double> params, std::array<double, 3> spherepos, double sphererad, int n);

std::array<std::vector<double>, 3> coord_frequency(std::string method, double end, double step, std::array<double, 3> ini, std::vector<double> params, int n);

double test(double num);