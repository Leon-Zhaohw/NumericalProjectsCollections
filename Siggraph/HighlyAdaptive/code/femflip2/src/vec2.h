/*
 *	vec2d.h
 *	
 *	Created by Ryoichi Ando on 11/3/11
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#include "macros.h"
#include <math.h>
#ifndef _VEC2_H
#define _VEC2_H

template <class T> struct vec2 {
	T v[2];
	vec2() {
		v[0] = v[1] = 0.0;
	}
	vec2( const vec2<int> &vf) {
		v[0] = vf.v[0];
		v[1] = vf.v[1];
	}
	vec2( const vec2<float> &vf) {
		v[0] = vf.v[0];
		v[1] = vf.v[1];
	}
	vec2( const vec2<double> &vd) {
		v[0] = vd.v[0];
		v[1] = vd.v[1];
	}
	vec2( const vec2<long double> &vd) {
		v[0] = vd.v[0];
		v[1] = vd.v[1];
	}
	vec2( FLOAT64 x, FLOAT64 y ) {
		v[0] = x;
		v[1] = y;
	}
	const T &operator[](uint idx) const {
		return v[idx];
	}
	T &operator[](uint idx) {
		return v[idx];
	}
	vec2 operator+=(const vec2 &vec) {
		v[0]+=vec[0];
		v[1]+=vec[1];
		return *this;
	}
	vec2 operator+(const vec2 &vec) const {
		return vec2( v[0]+vec[0], v[1]+vec[1] );
	}
	vec2 operator-(const vec2 &vec) const {
		return vec2( v[0]-vec[0], v[1]-vec[1] );
	}
	vec2 operator*=(T s) {
		v[0] *= s;
		v[1] *= s;
		return *this;
	}
	T operator^(const vec2 &vec) const {
		return v[0]*vec[1]-v[1]*vec[0];
	}
	
	vec2 operator-=(const vec2 &vec) {
		v[0] -= vec[0];
		v[1] -= vec[1];
		return *this;
	}
	vec2 operator*(FLOAT64 s) const {
		return vec2( s*v[0], s*v[1] );
	}
	vec2 operator/(FLOAT64 s) const {
		return vec2( v[0]/s, v[1]/s );
	}
	T operator*(vec2<T> vec) const {
		return v[0]*vec[0]+v[1]*vec[1];
	}
	T len2() const {
		return v[0]*v[0]+v[1]*v[1];
	}
	T len() const {
		return sqrtf(len2());
	}
	vec2 normal() const {
		vec2 copy = *this;
		copy.normalize();
		return copy;
	}
	bool normalize() {
		T length = len();
		if( length ) {
			v[0] /= length;
			v[1] /= length;
		}
		return length > 0.0;
	}
	vec2 rotate() const {
		return vec2(v[1],-v[0]);
	}
};

template <class T> static inline vec2<T> operator*(FLOAT64 s, const vec2<T> &vec) {
	return vec*s;
}

typedef vec2<FLOAT32>	vec2f;
typedef vec2<FLOAT64>	vec2d;
typedef vec2<int>		vec2i;
#endif