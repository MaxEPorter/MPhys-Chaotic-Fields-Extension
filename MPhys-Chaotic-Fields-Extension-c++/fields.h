#pragma once
#include<vector>
#include<math.h>

#include<boost/array.hpp>

typedef boost::array< double, 3 > state_type;

namespace fields {

	// create parent class Field to be inherited when creating other fields
	class Field {

	protected:

		// may have parameters to describe the field
		std::vector<double> params;

	public:

		Field();
		Field(std::vector<double> p);

		// 3 cartesian components of the field to be overridden in derived classes
		virtual double bx(const state_type& x) = 0;
		virtual double by(const state_type& x) = 0;
		virtual double bz(const state_type& x) = 0;

		// magnitude of the field at a point
		double magnitude(const state_type& x);

		// to return the field for odeint in the required form
		void ode(const state_type& x, state_type& dxdt, double s);

	};

	// field for straight infinite wire
	class Wire :public Field {

	public:

		double r2(const state_type& x); // radius for x, y
		double bx(const state_type& x);
		double by(const state_type& x);
		double bz(const state_type& x);

	};

	// field for abc field
	class ABC :public Field {

	public:

		ABC(std::vector<double>);
		double bx(const state_type& x);
		double by(const state_type& x);
		double bz(const state_type& x);

	};

	// field for double abc field
	class DoubleABC :public Field {

	public:

		DoubleABC(std::vector<double>);
		double bx(const state_type& x);
		double by(const state_type& x);
		double bz(const state_type& x);

	};

	// uniform field
	class Uniform :public Field {

	public:

		Uniform(std::vector<double>);
		double bx(const state_type& x);
		double by(const state_type& x);
		double bz(const state_type& x);

	};

	class Exptest :public Field {

	public:

		double bx(const state_type& x);
		double by(const state_type& x);
		double bz(const state_type& x);

	};


}