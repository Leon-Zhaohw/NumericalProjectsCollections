// templated point3 class

#ifndef __POINT3__

#define __POINT3__

#include <iostream>
#include <fstream>

#ifndef ZEROTOLERANCE
#define ZEROTOLERANCE 1.0e-7
#endif

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////

// Point3 class
template <class T1>
class point3 {
public:
   T1 x, y, z;

   // Constructors
   point3(T1 _x = 0, T1 _y = 0, T1 _z = 0) : x(_x), y(_y), z(_z) {}
   template <class T2>
      point3(const point3<T2> &p) : x((T1)p.x), y((T1)p.y), z((T1)p.z) {}
   void assign(T1 _x = 0, T1 _y = 0, T1 _z = 0) { x = _x; y = _y; z = _z; }

   // Binary operators
   template <class T2>
      inline point3 &operator=(const point3<T2> &p)
      { x = (T1)p.x; y = (T1)p.y; z = (T1)p.z; return *this; }
   template <class T2>
      inline point3 operator+(const point3<T2> &p) const
      { return point3(p.x + x, p.y + y, p.z + z); }
   template <class T2>
      inline point3 operator-(const point3<T2> &p) const
      { return point3(x - p.x, y - p.y, z - p.z); }
   template <class T2>
      inline T1 operator*(const point3<T2> &p) const
      { return (T1)(p.x * x + p.y * y + p.z * z); }
   inline point3 operator*(const T1 &c) const
      { return point3((T1)(c * x), (T1)(c * y), (T1)(c * z)); }
   inline point3 operator/(const T1 &c) const
      { return point3((T1)(x / c), (T1)(y / c), (T1)(z / c)); }

   // Unary operators
   inline point3 &operator+() const { return *this; }
   inline point3 operator-() const { return point3(-x, -y, -z); }

   // Self-modification operators
   template <class T2>
      inline point3 &operator+=(const point3<T2> &p) { x += (T1)p.x; y +=
	 (T1)p.y; z += (T1)p.z; return *this; }
   template <class T2>
      inline point3 &operator-=(const point3<T2> &p) { x -= (T1)p.x; y -=
	 (T1)p.y; z -= (T1)p.z; return *this; }
   template <class M>
      inline point3 &operator*=(const M c) { x *= c; y *= c; z *= c; return *this; }
   template <class M>
      inline point3 &operator/=(const M c) { x /= c; y /= c; z /= c; return *this; }

   // Comparison operators
   template <class T2>
      inline int operator==(const point3<T2> &p) const
      { return ((p.x == x) && (p.y == y) && (p.z == z)); }
   template <class T2>
      inline int operator!=(const point3<T2> &p) const
      { return ((p.x != x) || (p.y != y) || (p.z != z)); }

   // Misc math methods
   inline point3<T1> &normalize() { T1 m = T1(mag()); if (!m) m = 1;
      x /= m; y /= m; z /= m; return(*this); }
   inline T1 magsq() const { return dot(*this, *this); }
   inline double mag() const { return sqrt((double)magsq()); }
   inline double dir() const { return atan2((double)y, (double)x); }
   inline double orientation() const {
      return ( (fabs((double)x) > ZEROTOLERANCE) ? atan(double(y) / double(x)) : M_PI_2 );
   }
   // Computes the q*v*conj(q), where q=a+u is a quaternion.
   inline point3<T1> applyQuaternion(T1 a, const point3<T1> &u) {
      // q*v*conj(q) = -(u.v)*u + 2a*(uxv) + a^2*v - (uxv)xu
      point3<T1> &v = *this;
      point3<T1> crs = cross(u,v);
      point3<T1> vec = a*a*v;
      vec -= (a+a)*crs;
      vec += dot(u,v)*u;
      vec -= cross(crs,u);

      return vec;
   }
   // Rotate via angle in the "positive" direction around the z-axis
   inline point3<T1> &rotateZ(double angle) {
       double c = cos(angle);
       double s = sin(angle);
       T1 px(x), py(y);
       x = px * c - py * s;
       y = px * s + py * c;
       return *this;
   }
};

///////////////////////////////////////////////////////////////////////////////////////////

// Other operations
template <class T1> inline T1 dot(const point3<T1> &p1, const point3<T1> &p2) {
   return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
}
template <class T1> inline point3<T1> cross(const point3<T1> &p1, const point3<T1> &p2) {
   point3<T1> p = point3<T1>(p1.y*p2.z-p2.y*p1.z, p1.z*p2.x-p2.z*p1.x, p1.x*p2.y-p2.x*p1.y);
   return p;
}

///////////////////////////////////////////////////////////////////////////////////////////

// Compute 2 * the volume of the tetrahedron spanned by the four points.
template <class T1>
inline T1
volume2(const point3<T1> &p1, const point3<T1> &p2, const point3<T1> &p3, const point3<T1> &p4) {
    return dot(p2-p1, cross(p3-p2, p4-p3));
}

///////////////////////////////////////////////////////////////////////////////////////////

// Output operators
template <class T1> inline istream &operator>>(istream &fin, point3<T1> &p)
{ return fin >> p.x >> p.y >> p.z; }
template <class T1> inline ostream &operator<<(ostream &fout, const point3<T1> &p)
{ return fout << '(' << p.x << ", " << p.y << ", " << p.z << ')'; }

///////////////////////////////////////////////////////////////////////////////////////////

// Other friend functions
template <class T1, class T2>
   inline point3<T2> operator*(const T1 c, const point3<T2> &p)
   { return point3<T2>((T2)(c * p.x), (T2)(c * p.y), (T2)(c * p.z)); }
//template <class T2>
//   inline point3<T2> operator*(const T2 &c, const point3<T2> &p);

///////////////////////////////////////////////////////////////////////////////////////////

#endif
