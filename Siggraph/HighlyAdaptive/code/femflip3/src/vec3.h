/*
 *	vec3.h
 *	
 *	Created by Ryoichi Ando on 11/20/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include <math.h>
#ifndef _VEC3_H
#define _VEC3_H

template <class T> struct vec3 {
	T v[3];
	vec3() {
		v[0] = v[1] = v[2] = 0.0;
	}
	vec3( const vec3<int> &vd) {
		v[0] = vd.v[0];
		v[1] = vd.v[1];
		v[2] = vd.v[2];
	}
	vec3( const vec3<double> &vd) {
		v[0] = vd.v[0];
		v[1] = vd.v[1];
		v[2] = vd.v[2];
	}
	vec3( const vec3<float> &vd) {
		v[0] = vd.v[0];
		v[1] = vd.v[1];
		v[2] = vd.v[2];
	}
	vec3( const vec3<long double> &vd) {
		v[0] = vd.v[0];
		v[1] = vd.v[1];
	}
	vec3( FLOAT64 x, FLOAT64 y, FLOAT64 z ) {
		v[0] = x;
		v[1] = y;
		v[2] = z;
	}
	const T &operator[](uint idx) const {
		return v[idx];
	}
	T &operator[](uint idx) {
		return v[idx];
	}
	vec3 operator+=(const vec3 &vec) {
		v[0]+=vec[0];
		v[1]+=vec[1];
		v[2]+=vec[2];
		return *this;
	}
	vec3 operator+(const vec3 &vec) const {
		return vec3( v[0]+vec[0], v[1]+vec[1], v[2]+vec[2] );
	}
	vec3 operator-(const vec3 &vec) const {
		return vec3( v[0]-vec[0], v[1]-vec[1], v[2]-vec[2] );
	}
	vec3 operator*=(FLOAT64 s) {
		v[0] *= s;
		v[1] *= s;
		v[2] *= s;
		return *this;
	}
	vec3 operator^(const vec3 &vec) const {
		return vec3(v[1]*vec[2]-v[2]*vec[1],v[2]*vec[0]-v[0]*vec[2],v[0]*vec[1]-v[1]*vec[0]);
	}
	vec3 operator-=(const vec3 &vec) {
		v[0] -= vec[0];
		v[1] -= vec[1];
		v[2] -= vec[2];
		return *this;
	}
	vec3 operator*(T s) const {
		return vec3( s*v[0], s*v[1], s*v[2] );
	}
	vec3 operator/(T s) const {
		return vec3( v[0]/s, v[1]/s, v[2]/s );
	}
	T operator*(vec3<T> vec) const {
		return v[0]*vec[0]+v[1]*vec[1]+v[2]*vec[2];
	}
	T len2() const {
		return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	}
	T len() const {
		return sqrtf(len2());
	}
	vec3 normal() const {
		vec3 copy = *this;
		copy.normalize();
		return copy;
	}
	bool normalize() {
		T length = len();
		if( length ) {
			v[0] /= length;
			v[1] /= length;
			v[2] /= length;
		}
		return length > 0.0;
	}
};

template <class T> static inline vec3<T> operator*(FLOAT64 s, const vec3<T> &vec) {
	return vec*s;
}

typedef vec3<FLOAT32>	vec3f;
typedef vec3<FLOAT64>	vec3d;
typedef vec3<int>		vec3i;

#endif