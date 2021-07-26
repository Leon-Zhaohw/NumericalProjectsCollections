#ifndef VEC2_H
#define VEC2_H

#include <cassert>
#include <cmath>
#include <iostream>
#include "util.h"

template<class T>
struct Vec2
{
   T v[2];

   Vec2(void)
   {}

   Vec2(T value_for_all)
   { v[0]=v[1]=value_for_all; }

   template<class S>
   Vec2(const Vec2<S> source)
   { v[0]=(T)source.v[0]; v[1]=(T)source.v[1]; }

   template<class S>
   Vec2(const S source[2])
   { v[0]=(T)source[0]; v[1]=(T)source[1]; }

   Vec2(T v0, T v1)
   { v[0]=v0; v[1]=v1; }

   T &operator[](int index)
   {
      assert(0<=index && (unsigned int)index<2);
      return v[index];
   }

   const T &operator[](int index) const
   {
      assert(0<=index && (unsigned int)index<2);
      return v[index];
   }

   Vec2<T> operator+=(const Vec2<T> &w)
   { v[0]+=w.v[0]; v[1]+=w.v[1]; return *this; }

   Vec2<T> operator-=(const Vec2<T> &w)
   { v[0]-=w.v[0]; v[1]-=w.v[1]; return *this; }

   template<class S>
   Vec2<T> operator*=(S scalar)
   { v[0]*=scalar; v[1]*=scalar; return *this; }

   template<class S>
   Vec2<T> operator/=(S scalar)
   { v[0]/=scalar; v[1]/=scalar; return *this; }
};

template<class T> inline T mag2(const Vec2<T> &a)
{ return a.v[0]*a.v[0] + a.v[1]*a.v[1]; }

template<class T> inline T mag(const Vec2<T> &a)
{ return std::sqrt(mag2(a)); }

template<class T> inline T dist2(const Vec2<T> &a, const Vec2<T> &b)
{ return sqr(a.v[0]-b.v[0]) + sqr(a.v[1]-b.v[1]); }

template<class T> inline T dist(const Vec2<T> &a, const Vec2<T> &b)
{ return std::sqrt(dist2(a,b)); }

template<class T> inline bool operator==(const Vec2<T> &a, const Vec2<T> &b)
{ return a.v[0]==b.v[0] && a.v[1]==b.v[1]; }

template<class T> inline bool operator!=(const Vec2<T> &a, const Vec2<T> &b)
{ return a.v[0]!=b.v[0] || a.v[1]!=b.v[1]; }

template<class T> inline Vec2<T> operator-(const Vec2<T> &a)
{ return Vec2<T>(-a.v[0], -a.v[1]); }

template<class T> inline Vec2<T> operator+(const Vec2<T> &a, const Vec2<T> &b)
{ return Vec2<T>(a.v[0]+b.v[0], a.v[1]+b.v[1]); }

template<class T>
inline Vec2<T> operator-(const Vec2<T> &a, const Vec2<T> &b)
{ return Vec2<T>(a.v[0]-b.v[0], a.v[1]-b.v[1]); }

template<class S, class T>
inline Vec2<T> operator*(const Vec2<T> &a, S scalar)
{ return Vec2<T>(scalar*a.v[0], scalar*a.v[1]); }

template<class S, class T>
inline Vec2<T> operator*(S scalar, const Vec2<T> &a)
{ return Vec2<T>(scalar*a.v[0], scalar*a.v[1]); }

template<class S, class T>
inline Vec2<T> operator/(const Vec2<T> &a, S scalar)
{ return Vec2<T>(a.v[0]/scalar, a.v[1]/scalar); }

template<class T> inline T dot(const Vec2<T> &a, const Vec2<T> &b)
{ return a.v[0]*b.v[0] + a.v[1]*b.v[1]; }

template<class T> inline T cross(const Vec2<T> &a, const Vec2<T> &b)
{ return a.v[0]*b.v[1]-a.v[1]*b.v[0]; }

template<class T> inline Vec2<T> perp(const Vec2<T> &a)
{ return Vec2<T>(-a.v[1], a.v[0]); }

template<class T> inline void normalize(Vec2<T> &a)
{ a/=mag(a); }

template<class T> inline Vec2<T> normalized(const Vec2<T> &a)
{ return a/mag(a); }

template<class T>
inline std::ostream &operator<<(std::ostream &out, const Vec2<T> &a)
{ return out<<a.v[0]<<' '<<a.v[1]; }

template<class T>
inline std::istream &operator>>(std::istream &in, Vec2<T> &a)
{ return in>>a.v[0]>>a.v[1]; }

// common types of vectors ===================================================
typedef Vec2<float> Vec2f;
typedef Vec2<double> Vec2d;
typedef Vec2<int> Vec2i;

// type-specific operations ==================================================

/* may raise problems if hash() isn't defined for ints yet
// presupposes a hash() function defined for ints
inline unsigned int hash(const Vec2i &a)
{ return hash(a.v[0]) ^ a.v[1]; }
*/

template<class T> inline Vec2i round(const Vec2<T> &a)
{ return Vec2i(lround(a.v[0]), lround(a.v[1])); }

#endif
