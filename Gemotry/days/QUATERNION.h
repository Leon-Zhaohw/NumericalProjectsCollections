/*
QUIJIBO: Source code for the paper Symposium on Geometry Processing
         2015 paper "Quaternion Julia Set Shape Optimization"
Copyright (C) 2015  Theodore Kim

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _QUATERNION_H
#define _QUATERNION_H

#include <SETTINGS.h>
#include <VEC3.h>
#include <MATRIX3.h>
#include <iostream>

using namespace std;

class QUATERNION {

public:
  QUATERNION();
  QUATERNION(Real w, Real x, Real y, Real z);
  QUATERNION(Real w, Real x);
  QUATERNION(const VEC3F& vector);
  QUATERNION(const QUATERNION& q);
  QUATERNION(const MATRIX3& R);

  // overload operators
  QUATERNION& operator=(const VEC3F& vec);
  inline QUATERNION& operator=(const QUATERNION q)   { _w = q.w(); _x = q.x(); _y = q.y(); _z = q.z(); return *this;};
  inline QUATERNION& operator*=(const Real r)        { _w *= r; _x *= r; _y *= r; _z *= r;  return *this; };
  inline QUATERNION& operator*=(const QUATERNION& q) { 
    const Real& x = this->_y * q._z - this->_z * q._y + q._w * this->_x + this->_w * q._x;
    const Real& y = this->_z * q._x - this->_x * q._z + q._w * this->_y + this->_w * q._y;
    const Real& z = this->_x * q._y - this->_y * q._x + q._w * this->_z + this->_w * q._z;
    const Real& w = this->_w * q._w - this->_x * q._x - q._y * this->_y - this->_z * q._z;
    _x = x; _y = y; _z = z; _w = w;
    return *this;
  };
  inline QUATERNION& operator-=(const QUATERNION& q) { _w -= q.w(); _x -= q.x(); _y -= q.y(); _z -= q.z(); return *this; };
  inline QUATERNION& operator+=(const QUATERNION& q) { _w += q.w(); _x += q.x(); _y += q.y(); _z += q.z(); return *this; };
  inline Real& operator[](const int x) { return _entries[x]; };
  inline Real operator[](const int x) const { return _entries[x]; };

  // compute the conjugate
  QUATERNION conjugate() const;

  // compute the inverse
  QUATERNION inverse() const;

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
  
  Real& x() { return _x; };
  Real& y() { return _y; };
  Real& z() { return _z; };
  Real& w() { return _w; };

  // 2 norm of components
  inline Real magnitude() const { return sqrt(_w * _w + _x * _x + _y * _y + _z * _z); };
  inline Real magnitudeSq() const { return _w * _w + _x * _x + _y * _y + _z * _z; };

  // in-place equality
  void equals(const QUATERNION& q);

  // check if any components are nan
  inline bool anyNans() { return _x != _x || _y != _y || _z != _z || _w != _w; };

  // do a julia set iteration
  void juliaIteration(const QUATERNION& c);

  void write(FILE* file) const;
  void writeTextfile(FILE* file) const;
  void read(FILE* file);
  void readTextfile(FILE* file);

  inline void multiplyAdd(const QUATERNION& point, const QUATERNION& add)
  {
    const Real x = _x;
    const Real y = _y;
    const Real z = _z;
    const Real w = _w;
    _x = y * point._z - z * point._y + point._w * x + w * point._x + add._x;
    _y = z * point._x - x * point._z + point._w * y + w * point._y + add._y;
    _z = x * point._y - y * point._x + point._w * z + w * point._z + add._z;
    _w = w * point._w - x * point._x - point._y * y - z * point._z + add._w;
  };

  static bool wCompare(const QUATERNION& i, const QUATERNION& j) { return i._w < j._w; };
  static bool xCompare(const QUATERNION& i, const QUATERNION& j) { return i._x < j._x; };
  static bool yCompare(const QUATERNION& i, const QUATERNION& j) { return i._y < j._y; };
  static bool zCompare(const QUATERNION& i, const QUATERNION& j) { return i._z < j._z; };

  // take the exponential
  // from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
  QUATERNION exp() const;

  // take the log
  // from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
  QUATERNION log() const;

  // take the power
  // from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
  QUATERNION pow(const Real& exponent) const;

  // dot against another quaternion
  inline Real dot(const QUATERNION& rhs) const { return _w * rhs._w + _x * rhs._x + _y * rhs._y + _z * rhs._z; };

  // convert to rotation matrix -- should only use this when
  // the quaternion is *known* to be unit length!
  MATRIX3 toRotationMatrix() const;

  // get the axis-angle representation
  void axisAngle(VEC3F& axis, Real& angle);

  void debugMultiply(const QUATERNION& right);

  union {
     struct { Real _w,_x,_y,_z; };
     Real _entries[4];
  };
};

QUATERNION operator*(const QUATERNION& left, const QUATERNION& right);
QUATERNION operator/(const QUATERNION& left, const QUATERNION& right);
QUATERNION operator/(const QUATERNION& left, const Real& right);
Real operator^(const QUATERNION& left, const QUATERNION& right);
ostream &operator<<(ostream &out, const QUATERNION& q);

inline QUATERNION operator*(const QUATERNION& left, const Real& right) {
  return QUATERNION(left.w() * right, left.x() * right,
                    left.y() * right, left.z() * right);
};
inline QUATERNION operator*(const Real& left, const QUATERNION& right) {
  return QUATERNION(right.w() * left, right.x() * left,
                    right.y() * left, right.z() * left);
};
inline QUATERNION operator-(const QUATERNION& left, const QUATERNION& right) {
  return QUATERNION(left.w() - right.w(), left.x() - right.x(),
                    left.y() - right.y(), left.z() - right.z());
};
inline QUATERNION operator+(const QUATERNION& left, const QUATERNION& right) {
  return QUATERNION(left.w() + right.w(), left.x() + right.x(),
                    left.y() + right.y(), left.z() + right.z());
};
#endif
