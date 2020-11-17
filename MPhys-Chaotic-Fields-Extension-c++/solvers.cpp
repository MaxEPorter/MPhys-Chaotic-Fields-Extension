#include"solvers.h"

# define M_PI 3.14159265358979323846

Solution abc_field(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params) {

	Solution to_return;  // initialise solution to return
	to_return.params = params;
	to_return.step = t_step;
	to_return.method = "runge_kutta_fehlberg78";
	to_return.start = x_ini;
	to_return.field = "abc";

	fields::ABC abc(params);  // initialise abc field
	boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;  // initialise rkf stepper/method
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };  // initialise initial starting positions

	boost::numeric::odeint::integrate_const(stepper,  // pass rkf stepper
		[&abc](const state_type& x, state_type& dxdt, double s) {abc.ode(x, dxdt, s); },  // give the function to integrate as abc.ode
		x,  // give the x, y, z positions
		t_start,  // starting s
		t_end,  // ending s
		t_step,  // s step
		// lambda to save field trajectory data into solution
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;  // return solution

}
Solution abc_field_euler(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params) {

	Solution to_return;
	fields::ABC abc(params);
	boost::numeric::odeint::euler<state_type> stepper;
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };

	boost::numeric::odeint::integrate_const(stepper,
		[&abc](const state_type& x, state_type& dxdt, double s) {abc.ode(x, dxdt, s); },
		x,
		t_start,
		t_end,
		t_step,
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;

}
Solution abc_field_rk4(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params) {

	Solution to_return;
	fields::ABC abc(params);
	boost::numeric::odeint::runge_kutta4<state_type> stepper;
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };

	boost::numeric::odeint::integrate_const(stepper,
		[&abc](const state_type& x, state_type& dxdt, double s) {abc.ode(x, dxdt, s); },
		x,
		t_start,
		t_end,
		t_step,
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;

}


Solution wire_field(double t_start, double t_end, double t_step, std::array<double, 3> x_ini) {

	Solution to_return;
	fields::Wire w;
	boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };

	boost::numeric::odeint::integrate_const(stepper,
		[&w](const state_type& x, state_type& dxdt, double s) {w.ode(x, dxdt, s); },
		x,
		t_start,
		t_end,
		t_step,
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;

}
Solution wire_field_euler(double t_start, double t_end, double t_step, std::array<double, 3> x_ini) {

	Solution to_return;
	fields::Wire w;
	boost::numeric::odeint::euler<state_type> stepper;
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };

	boost::numeric::odeint::integrate_const(stepper,
		[&w](const state_type& x, state_type& dxdt, double s) {w.ode(x, dxdt, s); },
		x,
		t_start,
		t_end,
		t_step,
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;

}
Solution wire_field_rk4(double t_start, double t_end, double t_step, std::array<double, 3> x_ini) {

	Solution to_return;
	fields::Wire w;
	boost::numeric::odeint::runge_kutta4<state_type> stepper;
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };

	boost::numeric::odeint::integrate_const(stepper,
		[&w](const state_type& x, state_type& dxdt, double s) {w.ode(x, dxdt, s); },
		x,
		t_start,
		t_end,
		t_step,
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;

}

Solution double_abc_field(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params) {

	Solution to_return;
	fields::DoubleABC d_abc(params);
	boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };

	boost::numeric::odeint::integrate_const(stepper,
		[&d_abc](const state_type& x, state_type& dxdt, double s) {d_abc.ode(x, dxdt, s); },
		x,
		t_start,
		t_end,
		t_step,
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;

}

Solution uniform_field(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params) {

	Solution to_return;
	fields::Uniform u(params);
	boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };

	boost::numeric::odeint::integrate_const(stepper,
		[&u](const state_type& x, state_type& dxdt, double s) {u.ode(x, dxdt, s); },
		x,
		t_start,
		t_end,
		t_step,
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;

}
Solution uniform_field_euler(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params) {

	Solution to_return;
	fields::Uniform u(params);
	boost::numeric::odeint::euler<state_type> stepper;
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };

	boost::numeric::odeint::integrate_const(stepper,
		[&u](const state_type& x, state_type& dxdt, double s) {u.ode(x, dxdt, s); },
		x,
		t_start,
		t_end,
		t_step,
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;

}
Solution uniform_field_rk4(double t_start, double t_end, double t_step, std::array<double, 3> x_ini, std::vector<double> params) {

	Solution to_return;
	fields::Uniform u(params);
	boost::numeric::odeint::runge_kutta4<state_type> stepper;
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };

	boost::numeric::odeint::integrate_const(stepper,
		[&u](const state_type& x, state_type& dxdt, double s) {u.ode(x, dxdt, s); },
		x,
		t_start,
		t_end,
		t_step,
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;

}

Solution exp_field(double t_start, double t_end, double t_step, std::array<double, 3> x_ini) {

	Solution to_return;
	fields::Exptest e;
	boost::numeric::odeint::runge_kutta_fehlberg78<state_type> stepper;
	state_type x{ x_ini[0], x_ini[1], x_ini[2] };

	boost::numeric::odeint::integrate_const(stepper,
		[&e](const state_type& x, state_type& dxdt, double s) {e.ode(x, dxdt, s); },
		x,
		t_start,
		t_end,
		t_step,
		[&to_return](const state_type& x, const double t) {to_return.s.push_back(t); to_return.x.push_back(x[0]); to_return.y.push_back(x[1]); to_return.z.push_back(x[2]); });

	return to_return;

}

