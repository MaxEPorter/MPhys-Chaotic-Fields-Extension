#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include"fields.h"
#include"solvers.h"
#include"analyse.h"


PYBIND11_MODULE(chaoticfields, m) {

	pybind11::class_<Solution>(m, "Solution")
		//.def(pybind11::init<double&>())
		.def_readonly("s", &Solution::s)
		.def_readonly("x", &Solution::x)
		.def_readonly("y", &Solution::y)
		.def_readonly("z", &Solution::z)
		.def_readonly("step", &Solution::step)
		.def_readonly("params", &Solution::params)
		.def_readonly("start", &Solution::start)
		.def_readonly("method", &Solution::method)
		.def_readonly("field", &Solution::field);

	m.def("abc_field", &abc_field, R"pbdoc(
		abc_field line solver
	)pbdoc");
	m.def("abc_field_euler", &abc_field_euler, R"pbdoc(
		abc_field line solver with euler method
	)pbdoc");
	m.def("abc_field_rk4", &abc_field_rk4, R"pbdoc(
		abc_field line solver with rk4 method
	)pbdoc");
	m.def("wire_field", &wire_field, R"pbdoc(
		wire_field line solver
	)pbdoc");
	m.def("wire_field_euler", &wire_field_euler, R"pbdoc(
		wire_field line solver with euler method
	)pbdoc");
	m.def("wire_field_rk4", &wire_field_rk4, R"pbdoc(
		wire_field line solver with rk4 method
	)pbdoc");
	m.def("double_abc_field", &double_abc_field, R"pbdoc(
		double abc field line solver with rkf method
	)pbdoc");
	m.def("uniform_field", &uniform_field, R"pbdoc(
		uniform field line solver with rkf method
	)pbdoc");
	m.def("uniform_field_euler", &uniform_field_euler, R"pbdoc(
		uniform field line solver with euler method
	)pbdoc");
	m.def("uniform_field_rk4", &uniform_field_rk4, R"pbdoc(
		uniform field line solver with rk4 method
	)pbdoc");
	m.def("exp_field", &exp_field, R"pbdoc(
		exp field line solver with rk4 method
	)pbdoc");

	m.def("abc_poincare", &abc_poincare,
		R"pbdoc(abc field poincare)pbdoc");

	m.def("double_abc_poincare", &double_abc_poincare, R"pbdoc(
		double abc field poincare
	)pbdoc");

	m.def("periodic_projection", &periodic_projection, R"pbdoc(
		periodic projection
	)pbdoc");

	m.def("line_variance", &line_variance);

	m.def("trajectory_split", &trajectory_split);

	m.def("lyapunov", &lyapunov);

	m.def("line_distance", &line_distance);

	m.def("recurrence", &recurrence);

	m.def("coord_frequency", &coord_frequency);

	m.def("test", &test);

#ifdef VERSION_INFO
	m.attr("__version__") = VERSION_INFO;
#else
	m.attr("__version__") = "dev";
#endif
}