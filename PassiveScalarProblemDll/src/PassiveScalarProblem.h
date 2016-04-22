/*
*/
#pragma once
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>

template <typename R, int D = 2>
class PassiveScalarProblem {
	typedef Eigen::Matrix<R, D, 1> Vec;
public:
	PassiveScalarProblem() {}
	~PassiveScalarProblem() {}

	template <typename U>
	__forceinline R scalar(const U& v) const { 
		if(scalarArea(v)) return static_cast<R>(1.);
		else return static_cast<R>(0.);
	}

	template <typename U>
	__forceinline Vec velocity(const U& v) const {
		Vec ret;
		static const R omega = static_cast<R>(M_PI / 5.);
		const auto r = v.norm();
		if (r < 1.e-10) return Vec::Zero();
		const auto cTheta = v[0] / r;
		const auto sTheta = v[1] / r;
		const auto vel = omega * r;
		ret << -vel*sTheta, vel*cTheta;
		return ret;
	}

private:
	template <typename U>
	__forceinline bool scalarArea(const U& v) const {
		if (v[0] >= -0.1 && v[0] <= 0.1 && v[1] >= 0.16 && v[1] <= 0.36) return true;
		else return false;
	}
};
