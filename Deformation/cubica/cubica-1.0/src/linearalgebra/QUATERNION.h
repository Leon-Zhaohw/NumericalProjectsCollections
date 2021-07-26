/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef QUATERNION_H
#define QUATERNION_H

#include <SETTINGS.h>
#include <VEC3.h>
#include <MATRIX3.h>
#include <MATRIX.h>
#include <iostream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// A quaternion class created to support Georgii and Westermann's
// energy minimization based rotation extraction method. This class
// is *not* used extensively in the tet mesh code, and was created
// far, far after the fact.
//
// _w and _entries[3] are the real component
//////////////////////////////////////////////////////////////////////
class QUATERNION {

public:
  QUATERNION();
  QUATERNION(Real x, Real y, Real z, Real w);
  QUATERNION(const VEC3F& vector);
  QUATERNION(const VECTOR& vector);
  QUATERNION(const MATRIX3& matrix);

  // overload operators
  QUATERNION& operator=(const QUATERNION q);
  QUATERNION& operator=(const VEC3F& vec);
  QUATERNION& operator*=(const Real r);
  QUATERNION& operator-=(const QUATERNION& q);
  QUATERNION& operator+=(const QUATERNION& q);
  inline Real& operator[](const int x) { return _entries[x]; };
  inline Real operator[](const int x) const { return _entries[x]; };

  // compute the conjugate
  QUATERNION conjugate();

  // convert to rotation matrix -- should only use this when
  // the quaternion is *known* to be unit length!
  MATRIX3 toRotationMatrix();

  // convert to an explicit 4 x 4 matrix
  MATRIX toExplicitMatrix();
  
  // convert to an explicit 3 x 3 matrix -- (3,3) entry is carved out
  MATRIX3 toExplicitMatrix3x3();

  // populate a VECTOR of size 4
  VECTOR toVector();
  
  // populate a VECTOR of size 5, and put lambda in the last entry
  VECTOR toVector(Real lambda);

  // get the axis-angle representation
  void axisAngle(VEC3F& axis, Real& angle);

  // normalize
  void normalize();

  // negative the imaginary component, effectively retrieving the transpose
  // of the desired rotation
  void negateIm();

  // accessors
  const Real& x() const { return _x; };
  const Real& y() const { return _y; };
  const Real& z() const { return _z; };
  const Real& w() const { return _w; };

  // rotate the passed in vector by this
  QUATERNION rotates(QUATERNION& toRotate);
  VEC3F rotates(VEC3F& toRotate);
  MATRIX rotates(MATRIX& toRotate);

  // return the matrix corresponding to the gradient
  // 'which' refers to the component of the quaternion
  MATRIX gradient(int which);

  // 3x3 version of above
  MATRIX3 gradient3x3(int which);

  // return the matrix corresponding to the hessian
  // with respect to the entries (i,j)
  MATRIX hessian(int i, int j);
  
  // 3x3 version of above
  MATRIX3 hessian3x3(int i, int j);

  // 2 norm of components
  Real magnitude();

  // Shabana quantities: "Dynamics of Multibody Systems", third edition,
  // appears to be partial of quaternion wrt real component
  MATRIX E();
  MATRIX Ebar();
  MATRIX G();
  MATRIX Gbar();

  // convert from axis angle representation to quaternion
  static QUATERNION fromAxisAngle(VEC3F& axisAngle);

  // compute an Open Dynamics Engine (ODE)-style update vector
  static QUATERNION odeUpdate(Real& dt, VEC3F& angularVelocity);

  // in-place equality
  void equals(const QUATERNION& q);

private:
  // _w and _entries[3] are the real component
  union {
     struct { Real _x,_y,_z,_w; };
     Real _entries[4];
  };
};

QUATERNION operator-(const QUATERNION& left, const QUATERNION& right);
QUATERNION operator+(const QUATERNION& left, const QUATERNION& right);
QUATERNION operator*(const QUATERNION& left, const QUATERNION& right);
QUATERNION operator*(const QUATERNION& left, const Real& right);
QUATERNION operator*(const Real& left, const QUATERNION& right);
QUATERNION operator*(const MATRIX& left, const QUATERNION& right);
// assumes the final entry in the quaternion is zero
QUATERNION operator*(const MATRIX3& left, const QUATERNION& right);
Real operator^(const QUATERNION& left, const QUATERNION& right);

ostream &operator<<(ostream &out, const QUATERNION& q);

#endif
