// templated point2 class

#ifndef __POINT2__

#define __POINT2__

#include <iostream>
#include <fstream>

#ifndef ZEROTOLERANCE
#define ZEROTOLERANCE 1.0e-7
#endif

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////

template <class T1>
class point2 {
public:
   // Constructors
   point2(T1 _x = 0, T1 _y = 0) : x(_x), y(_y) {}
   template <class T2>
      point2(const point2<T2> &p) : x((T1)p.x), y((T1)p.y) {}

   // Binary operators
   template <class T2>
      inline point2 &operator=(const point2<T2> &p) { x = (T1)p.x; y = (T1)p.y; return *this; }
   template <class T2>
      inline point2 operator+(const point2<T2> &p) const { return point2(p.x + x, p.y + y); }
   template <class T2>
      inline point2 operator-(const point2<T2> &p) const { return point2(x - p.x, y - p.y); }
   template <class T2>
      inline T1 operator*(const point2<T2> &p) const { return (T1)(p.x * x + p.y * y); }
   inline point2 operator*(const T1 &c) const { return point2((T1)(c * x), (T1)(c * y)); }
   friend inline point2 operator*(const T1 &c, const point2 &p) { return point2((T1)(c * p.x), (T1)(c * p.y)); }
   inline point2 operator/(const T1 c) const { return point2((T1)(x / c), (T1)(y / c)); }

   // Unary operators
   inline point2 &operator+() const { return *this; }
   inline point2 operator-() const { return point2(-x, -y); }

   // Self-modification operators
   template <class T2>
      inline point2 &operator+=(const point2<T2> &p) { x += (T1)p.x; y += (T1)p.y; return *this; }
   template <class T2>
      inline point2 &operator-=(const point2<T2> &p) { x -= (T1)p.x; y -= (T1)p.y; return *this; }
   template <class M>
      inline point2 &operator*=(const M c) { x *= c; y *= c; return *this; }
   template <class M>
      inline point2 &operator/=(const M c) { x /= c; y /= c; return *this; }
   void assign(T1 _x, T1 _y) { x = _x; y = _y; }

   // Comparison operators
   template <class T2>
      inline int operator==(const point2<T2> &p) const { return ((p.x == x) && (p.y == y)); }
   template <class T2>
      inline int operator!=(const point2<T2> &p) const { return ((p.x != x) || (p.y != y)); }

   // Misc math methods
   inline T1 magsq() { return (*this * *this); }
   inline double mag() { return sqrt((double)(*this * *this)); }
   inline double dir() { return atan2((double)y, (double)x); }
   inline double orientation() { return ( (fabs((double)x) > ZEROTOLERANCE) ? atan(double(y) / double(x)) : M_PI_2 ); }
   inline void normalize() { double m = mag(); if (m > 0) { x /= m; y /= m; } }
   inline point2 normal() { return point2(y, -x); }

   // Misc properties
   inline int zero() { return ((x == (T1)0) && (y == (T1)0)); }

   T1 x, y;
};

///////////////////////////////////////////////////////////////////////////////////////////

// Output operators
template <class T1> inline istream &operator>>(istream &fin, point2<T1> &p) { return fin >> p.x >> p.y; }
template <class T1> inline ostream &operator<<(ostream &fout, const point2<T1> &p) { return fout << '(' << p.x << ", " << p.y << ')'; }

///////////////////////////////////////////////////////////////////////////////////////////

#endif // __POINT2__
