/*
*/
#pragma once

enum bvpType {
	poly_a,
	poly_b,
	sinu_a,
	sinu_b,
	expo_a,
	expo_b,
};

template <class R>
class Bvp {
public:
	typedef Eigen::Matrix<R, 2, 1> Vec;
	Bvp(const bvpType& _type) : type(_type) {
	}
	~Bvp() {}

	const R func(const Vec& v) const {
		switch (type) {
		case poly_a:
			return bvpFunc_poly_a(v);
		case poly_b:
			return bvpFunc_poly_b(v);
		case sinu_a:
			return bvpFunc_sinu_a(v);
		case sinu_b:
			return bvpFunc_sinu_b(v);
		case expo_a:
			return bvpFunc_expo_a(v);
		case expo_b:
			return bvpFunc_expo_b(v);
		default:
			std::cout << " !!! no type found !!! " << std::endl;
			return 0.;
		}
	}
	const R lap(const Vec& v) const {
		switch (type) {
		case poly_a:
			return bvpLap_poly_a(v);
		case poly_b:
			return bvpLap_poly_b(v);
		case sinu_a:
			return bvpLap_sinu_a(v);
		case sinu_b:
			return bvpLap_sinu_b(v);
		case expo_a:
			return bvpLap_expo_a(v);
		case expo_b:
			return bvpLap_expo_b(v);
		default:
			std::cout << " !!! no type found !!! " << std::endl;
			return 0.;
		}
	}
	const Vec grad(const Vec& v) const {
		switch (type) {
		case poly_a:
			return bvpGrad_poly_a(v);
		case poly_b:
			return bvpGrad_poly_b(v);
		case sinu_a:
			return bvpGrad_sinu_a(v);
		case sinu_b:
			return bvpGrad_sinu_b(v);
		case expo_a:
			return bvpGrad_expo_a(v);
		case expo_b:
			return bvpGrad_expo_b(v);
		default:
			std::cout << " !!! no type found !!! " << std::endl;
			return 0.;
		}
	}

public:
	bvpType type;

private:
	/*poly_a*/
	const R bvpFunc_poly_a(const Vec& v) const {
		return pow(v[0], 2) + pow(v[1], 2);
	}
	const R bvpLap_poly_a(const Vec& v) const {
		return 4.;
	}
	const Vec bvpGrad_poly_a(const Vec& v) const {
		return Vec(2.*v[0], 0., 2.*v[1]);
	}

	/*poly_b*/
	const R bvpFunc_poly_b(const Vec& v) const {
		return pow(v[0], 2) + pow(v[1], 2) + v[0]*v[1];
	}
	const R bvpLap_poly_b(const Vec& v) const {
		return 4.;
	}
	const Vec bvpGrad_poly_b(const Vec& v) const {
		return Vec(2.*v[0] + v[1], 0., 2.*v[1] + v[0]);
	}

	/*sinu_a*/
	const R bvpFunc_sinu_a(const Vec& v) const {
		return (cos(M_PI*v[0]) + cos(M_PI*v[1]));
	}
	const R bvpLap_sinu_a(const Vec& v) const {
		return (-M_PI*M_PI *(cos(M_PI*v[0]) + cos(M_PI*v[1])));
	}
	const Vec bvpGrad_sinu_a(const Vec& v) const {
		return Vec(-M_PI*sin(M_PI* v[0]), 0., -M_PI*sin(M_PI* v[1]));
	}

	/*sinu_b*/
	const R bvpFunc_sinu_b(const Vec& v) const {
		return (cos(M_PI*v[0]) + cos(M_PI*v[1]) + cos(M_PI*v[0]*v[1]));
	}
	const R bvpLap_sinu_b(const Vec& v) const {
		return (-M_PI*M_PI *(cos(M_PI*v[0]) + cos(M_PI*v[1]) + cos(M_PI*v[0]*v[1])*(v[0]*v[0] + v[1]*v[1])));
	}
	const Vec bvpGrad_sinu_b(const Vec& v) const {
		return Vec(-M_PI*sin(M_PI* v[0]) - M_PI*v[1]*sin(M_PI*v[0]*v[1]), 0., -M_PI*sin(M_PI* v[1]) - M_PI*v[0]*sin(M_PI*v[0]*v[1]));
	}

	/*expo_a*/
	const R bvpFunc_expo_a(const Vec& v) const {
		return (exp(M_PI*v[0]) + exp(M_PI*v[1]));
	}
	const R bvpLap_expo_a(const Vec& v) const {
		return (M_PI*M_PI *(exp(M_PI*v[0]) + exp(M_PI*v[1])));
	}
	const Vec bvpGrad_expo_a(const Vec& v) const {
		return Vec(M_PI*exp(M_PI* v[0]), 0., M_PI*exp(M_PI* v[1]));
	}

	/*expo_b*/
	const R bvpFunc_expo_b(const Vec& v) const {
		return (exp(M_PI*v[0]) + exp(M_PI*v[1]) + exp(M_PI*v[0]*v[1]));
	}
	const R bvpLap_expo_b(const Vec& v) const {
		return (M_PI*M_PI *(exp(M_PI*v[0]) + exp(M_PI*v[1]) + (v[1]*v[1] + v[0]*v[0])*exp(M_PI*v[0]*v[1])));
	}
	const Vec bvpGrad_expo_b(const Vec& v) const {
		return Vec(M_PI*(exp(M_PI* v[0]) + v[1]*exp(M_PI*v[0]*v[1])), 0., M_PI*(exp(M_PI* v[1]) + v[0]*exp(M_PI*v[0]*v[1])));
	}

};

