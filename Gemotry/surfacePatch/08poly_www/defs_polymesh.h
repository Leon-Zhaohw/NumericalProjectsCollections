
#ifndef DEFS_POLYMESH_H
#define DEFS_POLYMESH_H

#include <complex>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>

typedef  float                      REAL;

#define REAL_format "lf"  // for sprintf

#include "point2.h"
#include "point3.h"

// CGAL types
typedef  CGAL::Cartesian<REAL>       Kernel;
typedef  Kernel::Point_3       KPoint;
typedef  KPoint Point;
typedef  Kernel::Vector_3      Vector;
typedef  point3<REAL>                RealPoint;
typedef  point2<REAL>                RealPoint2;

// Conversion amongst Kpoint and RealPoint
inline RealPoint KToRealPoint(const KPoint &p) {
    return RealPoint(p.x(), p.y(), p.z());
}
inline KPoint RealToKPoint(const RealPoint &p) {
    return KPoint(p.x, p.y, p.z);
}

// Apply a threshold
#define THRESH 1.0e-8
inline REAL thresh(REAL val) {
    return (abs(val) < THRESH) ? 0 : val;
}

////////////////////////////////////////////////////////////////////////////////

typedef point3<complex<REAL> > ComplexPoint;

template <typename T>
inline point3<complex<T> > RealToComplexPoint(const point3<T> &pr, const point3<T> &pi = point3<T>())
{ return point3<complex<T> >(complex<T>(pr.x, pi.x), complex<T>(pr.y, pi.y), complex<T>(pr.z, pi.z)); }
template <typename T>
inline point3<T> real(const point3<complex<T> > &p) { return point3<T>(real(p.x), real(p.y), real(p.z)); }
template <typename T>
inline point3<T> imag(const point3<complex<T> > &p) { return point3<T>(imag(p.x), imag(p.y), imag(p.z)); }
template <typename T>
inline point3<complex<T> > operator*(const point3<T> &p, const complex<T> &c) { return point3<complex<T> >(p.x*c, p.y*c, p.z*c); }

////////////////////////////////////////////////////////////////////////////////

// Print and exit
#define die(text) cout << text << endl, exit(0)

#endif // DEFS_POLYMESH_H

