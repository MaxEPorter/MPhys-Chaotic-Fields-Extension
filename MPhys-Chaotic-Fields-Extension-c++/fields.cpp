#include"fields.h"


namespace fields {

	Field::Field() {}
	Field::Field(std::vector<double> p) { params = p; }  // parametarised constructor to fill params
	double Field::magnitude(const state_type& x) { return sqrt(pow(bx(x), 2) + pow(by(x), 2) + pow(bz(x), 2)); }
	void Field::ode(const state_type& x, state_type& dxdt, double s) {

		// to calculate magnetic field lines
		dxdt[0] = bx(x) / magnitude(x);
		dxdt[1] = by(x) / magnitude(x);
		dxdt[2] = bz(x) / magnitude(x);
	}

	double Wire::r2(const state_type& x) { return pow(x[0], 2) + pow(x[1], 2); }
	double Wire::bx(const state_type& x) { return -x[1] / r2(x); }
	double Wire::by(const state_type& x) { return x[0] / r2(x); }
	double Wire::bz(const state_type& x) { return 0; }

	ABC::ABC(std::vector<double> p) :Field(p) {}
	double ABC::bx(const state_type& x) { return params[0] * sin(params[3] * x[2]) + params[2] * cos(params[3] * x[1]); }
	double ABC::by(const state_type& x) { return params[1] * sin(params[3] * x[0]) + params[0] * cos(params[3] * x[2]); }
	double ABC::bz(const state_type& x) { return params[2] * sin(params[3] * x[1]) + params[1] * cos(params[3] * x[0]); }

	DoubleABC::DoubleABC(std::vector<double> p) : Field(p) {}
	double DoubleABC::bx(const state_type& x) {
		return params[0] * sin(params[6] * x[2]) + params[4] * cos(params[6] * x[1]) +
			params[1] * sin(params[7] * x[2]) + params[5] * cos(params[7] * x[1]);
	}
	double DoubleABC::by(const state_type& x) {
		return params[2] * sin(params[6] * x[0]) + params[0] * cos(params[6] * x[2]) +
			params[3] * sin(params[7] * x[0]) + params[1] * cos(params[7] * x[2]);
	}
	double DoubleABC::bz(const state_type& x) {
		return params[4] * sin(params[6] * x[1]) + params[2] * cos(params[6] * x[0]) +
			params[5] * sin(params[7] * x[1]) + params[3] * cos(params[7] * x[0]);
	}

	Uniform::Uniform(std::vector<double> p) :Field(p) {}
	double Uniform::bx(const state_type& x) { return params[0]; }
	double Uniform::by(const state_type& x) { return params[1]; }
	double Uniform::bz(const state_type& x) { return params[2]; }

	double Exptest::bx(const state_type& x) { return pow(x[0], 2); }
	double Exptest::by(const state_type& x) { return 0; }
	double Exptest::bz(const state_type& x) { return 0; }

}