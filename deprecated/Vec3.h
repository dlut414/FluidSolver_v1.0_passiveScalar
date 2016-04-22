/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
#ifndef VEC3_H
#define VEC3_H

#include <cmath>

template <typename real>
class Vec3 {
public:
	template <typename T>
	Vec3(const Vec3<T>& v) : x(real(v.x)), y(real(v.y)), z(real(v.z)) {}
	Vec3(const real& _v = 0.) : x(_v), y(_v), z(_v) {}
	Vec3(const real& _x, const real& _y, const real& _z) : x(_x), y(_y), z(_z) {}
	~Vec3() {}

	inline const Vec3 operator+ (const Vec3& v) const {
		return Vec3(x + v.x, y + v.y, z + v.z);
	}

	inline const Vec3 operator- (const Vec3& v) const {
		return Vec3(x - v.x, y - v.y, z - v.z);
	}
	inline const Vec3 operator- () const {
		return Vec3(-x, -y, -z);
	}

	inline void operator+= (const Vec3& v) {
		x += v.x;   y += v.y;   z += v.z;
	}
	inline void operator-= (const Vec3& v) {
		x -= v.x;   y -= v.y;   z -= v.z;
	}

	inline const Vec3 operator* (const real& s) const {
		return Vec3(s*x, s*y, s*z);
	}

	inline const real operator* (const Vec3& v) const {
		return v.x*x + v.y*y + v.z*z;
	}

	inline const Vec3 operator/ (const real& s) const {
		return (1. / s) * (*this);
	}

	inline real mag2() const {
		return (x*x + y*y + z*z);
	}

	inline real mag() const {
		return sqrt(mag2());
	}

	inline Vec3 norm() const {
		real s = (*this).mag();
		if (s <= 0.) return Vec3(0., 0., 0.);
		else return (*this) / s;
	}

public:
	real x, y, z;
};

template <typename real>
inline const Vec3<real> operator* (const real& s, const Vec3<real>& p) {
	return Vec3<real>(s*p.x, s*p.y, s*p.z);
}

#endif // VEC3_H
