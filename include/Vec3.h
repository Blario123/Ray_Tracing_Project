#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

class Vec3 {
public:
	Vec3(double x_, double y_, double z_) {
		x = x_;
		y = y_;
		z = z_;
	}
	Vec3() { x = y = z = 0; }
	double norm() const { return sqrt(x * x + y * y + z * z); }
	double norm2() const { return x * x + y * y + z * z; }
	void normalise() {
		double d = norm();
		x /= d;
		y /= d;
		z /= d;
	}
	Vec3 &operator+=(const Vec3 &rhs) {
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}
	Vec3 &operator-=(const Vec3 &rhs) {
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return *this;
	}
	Vec3 &operator*=(const Vec3 &rhs) {
		x *= rhs.x;
		y *= rhs.y;
		z *= rhs.z;
		return *this;
	}
	Vec3 &operator*=(double rhs) {
		x *= rhs;
		y *= rhs;
		z *= rhs;
		return *this;
	}
	Vec3 &operator/=(double rhs) {
		x /= rhs;
		y /= rhs;
		z /= rhs;
		return *this;
	}
    void clear() {
        x = 0;
        y = 0;
        z = 0;
    }
	double x, y, z;
};

inline Vec3 operator+(const Vec3 &lhs, const Vec3 &rhs) { return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z}; }
inline Vec3 operator-(const Vec3 &lhs, const Vec3 &rhs) { return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z}; }
inline Vec3 operator-(const Vec3 &rhs) { return {-rhs.x, -rhs.y, -rhs.z}; }
inline Vec3 operator*(const double d, const Vec3 &rhs) { return {d * rhs.x, d * rhs.y, d * rhs.z}; }
inline Vec3 operator*(const Vec3 &rhs, const double d) { return {d * rhs.x, d * rhs.y, d * rhs.z}; }
inline Vec3 operator*(const Vec3 &lhs, const Vec3 &rhs) { return {lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z}; }
inline Vec3 operator/(const Vec3 &rhs, const double d) { return {rhs.x / d, rhs.y / d, rhs.z / d}; }
inline double dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3 &a, const Vec3 &b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline Vec3 pow(const Vec3 &a, const int &b) { return {pow(a.x, b), pow(a.y, b), pow(a.z, b)}; }
inline Vec3 pow(const Vec3 &a, const double &b) { return {pow(a.x, b), pow(a.y, b), pow(a.z, b)}; }
inline std::ostream &operator<<(std::ostream &os, const Vec3 &v) {
	os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return os;
}
typedef Vec3 Radiance;

#endif
