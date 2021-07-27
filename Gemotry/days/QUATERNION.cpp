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
#include "QUATERNION.h"

QUATERNION::QUATERNION() :
  _w(0), _x(0), _y(0), _z(0)
{
}

QUATERNION::QUATERNION(Real w, Real x, Real y, Real z) :
  _w(w), _x(x), _y(y), _z(z)
{
}

QUATERNION::QUATERNION(Real w, Real x) :
  _w(w), _x(x), _y(0), _z(0)
{
}

QUATERNION::QUATERNION(const VEC3F& vector) :
  _w(0), _x(vector[0]), _y(vector[1]), _z(vector[2])
{
}

QUATERNION::QUATERNION(const QUATERNION& q) :
  _w(q._w), _x(q._x), _y(q._y), _z(q._z)
{
}

QUATERNION::QUATERNION(const MATRIX3& matrix)
{
  float trace, u;
  trace = matrix(0,0) + matrix(1,1) + matrix(2,2);

  if (trace >= 0)
  {
    u = (float)sqrt(trace + 1);
    _w = (float)0.5 * u;
    u = (float)0.5 / u;
    _x = (matrix(2,1) - matrix(1,2)) * u;
    _y = (matrix(0,2) - matrix(2,0)) * u;
    _z = (matrix(1,0) - matrix(0,1)) * u;
  }
  else
  {
    int i = 0;
    if (matrix(1,1) > matrix(0,0))
      i = 1;

    //if (matrix(2,2) > (i == 0) ? matrix(0,0) : matrix(1,1))
    if (matrix(2,2) > matrix(i,i))
      i = 2;

    switch (i)
    {
      case 0:
        u = (float)sqrt((matrix(0,0) - (matrix(1,1) + matrix(2,2))) + 1);
        _x = 0.5f * u;
        u = 0.5f / u;
        _y = (matrix(1,0) + matrix(0,1)) * u;
        _z = (matrix(0,2) + matrix(2,0)) * u;
        _w = (matrix(2,1) - matrix(1,2)) * u;
      break;

      case 1:
        u = (float)sqrt((matrix(1,1) - (matrix(2,2) + matrix(0,0))) + 1);
        _y = 0.5f * u;
        u = 0.5f / u;
        _z = (matrix(2,1) + matrix(1,2)) * u;
        _x = (matrix(1,0) + matrix(0,1)) * u;
        _w = (matrix(0,2) - matrix(2,0)) * u;
      break;

      case 2:
        u = (float)sqrt((matrix(2,2) - (matrix(0,0) + matrix(1,1))) + 1);
        _z = 0.5f * u;
        u = 0.5f / u;
        _x = (matrix(0,2) + matrix(2,0)) * u;
        _y = (matrix(2,1) + matrix(1,2)) * u;
        _w = (matrix(1,0) - matrix(0,1)) * u;
      break;
    }
  }
}

QUATERNION QUATERNION::conjugate() const
{
  QUATERNION final(_w, -_x, -_y, -_z);
  return final;
}

QUATERNION& QUATERNION::operator=(const VEC3F& vec)
{
  this->_x = vec[0];
  this->_y = vec[1];
  this->_z = vec[2];
  this->_w = 0.0;
  return *this;
}

void QUATERNION::debugMultiply(const QUATERNION& right)
{
  QUATERNION final;
  cout << " \t\tx: " << final._x << " adding: " << _y * right._z << endl; 
  cout << " \t\tpieces: " << _y << " " << right._z << endl;
  final._x = _y * right._z;
  cout << " \t\tx: " << final._x << " adding: " << -_z * right._y << endl; 
  final._x += -_z * right._y; 
  cout << " \t\tx: " << final._x << " adding: " << right._w * _x << endl; 
  final._x += right._w * _x;
  cout << " \t\tx: " << final._x << " adding: " << _w * right._x << endl; 
  final._x += _w * right._x;
  final._y = _z * right._x - _x * right._z + right._w * _y + _w * right._y;
  final._z = _x * right._y - _y * right._x + right._w * _z + _w * right._z;
  final._w = _w * right._w - _x * right._x - right._y * _y - _z * right._z;
  *this = final;
}

QUATERNION operator*(const QUATERNION& left, const QUATERNION& right)
{
  QUATERNION final;
  final._x = left._y * right._z - left._z * right._y + right._w * left._x + left._w * right._x;
  final._y = left._z * right._x - left._x * right._z + right._w * left._y + left._w * right._y;
  final._z = left._x * right._y - left._y * right._x + right._w * left._z + left._w * right._z;
  final._w = left._w * right._w - left._x * right._x - right._y * left._y - left._z * right._z;
  return final;
}

//////////////////////////////////////////////////////////////////////
// compute the inverse
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::inverse() const
{
  const Real magnitudeSq = (*this)[0] * (*this)[0] + (*this)[1] * (*this)[1] +
                           (*this)[2] * (*this)[2] + (*this)[3] * (*this)[3];
  return QUATERNION(conjugate() / magnitudeSq);
}

QUATERNION operator/(const QUATERNION& left, const QUATERNION& right)
{
  // build the reciprocal of the right hand side
  Real magnitudeSq = right[0] * right[0] + right[1] * right[1] +
                     right[2] * right[2] + right[3] * right[3];
  QUATERNION rightInv = right.conjugate() / magnitudeSq;

  return left * rightInv;
}

Real operator^(const QUATERNION& left, const QUATERNION& right)
{
  return left._x * right._x + left._y * right._y +
         left._z * right._z + left._w * right._w;
}

ostream& operator<<(ostream &out, const QUATERNION& q)
{
  out << "(" << q._w << ", " << q._x << ", " << q._y << ", " << q._z << ")";
  return out;
}

QUATERNION operator/(const QUATERNION& left, const Real& right)
{
  const Real rightInv = 1.0 / right;
  return left * rightInv;
}

void QUATERNION::normalize()
{
  Real invMagnitude = 1.0 / sqrtf(_x * _x + _y * _y + _z * _z + _w * _w);
  _x *= invMagnitude;
  _y *= invMagnitude;
  _z *= invMagnitude;
  _w *= invMagnitude;
}

//////////////////////////////////////////////////////////////////////
// negate the imaginary component -- effetively transpose the rotation
//////////////////////////////////////////////////////////////////////
void QUATERNION::negateIm()
{
  _x *= -1.0;
  _y *= -1.0;
  _z *= -1.0;
}

//////////////////////////////////////////////////////////////////////
// in-place equality
//////////////////////////////////////////////////////////////////////
void QUATERNION::equals(const QUATERNION& q)
{
  _entries[0] = q._entries[0];
  _entries[1] = q._entries[1];
  _entries[2] = q._entries[2];
  _entries[3] = q._entries[3];
}

//////////////////////////////////////////////////////////////////////
// do a julia set iteration
//////////////////////////////////////////////////////////////////////
void QUATERNION::juliaIteration(const QUATERNION& c)
{
  const Real copy[] = {_x, _y, _z, _w};
  
  const Real copy03 = 2.0f * copy[0] * copy[3];
  const Real copy31 = 2.0f * copy[3] * copy[1];
  const Real copy32 = 2.0f * copy[3] * copy[2];

  _x = copy03 + c[0];
  _y = copy31 + c[1];
  _z = copy32 + c[2];
  _w = copy[3] * copy[3] - copy[0] * copy[0] - copy[1] * copy[1] - copy[2] * copy[2] + c[3];
}

void QUATERNION::writeTextfile(FILE* file) const
{
  double entries[] = {(double)_entries[0], 
                      (double)_entries[1], 
                      (double)_entries[2], 
                      (double)_entries[3]};
  fprintf(file, "(%.20e %.20e %.20e %.20e)\n", entries[0], entries[1], entries[2], entries[3]);
}

void QUATERNION::readTextfile(FILE* file)
{
  fscanf(file, "(%lf %lf %lf %lf)\n", &_entries[0], &_entries[1], &_entries[2], &_entries[3]);
}

void QUATERNION::write(FILE* file) const
{
  if (sizeof(Real) == sizeof(double))
  {
    fwrite((void*)&_entries[0], sizeof(Real), 1, file);
    fwrite((void*)&_entries[1], sizeof(Real), 1, file);
    fwrite((void*)&_entries[2], sizeof(Real), 1, file);
    fwrite((void*)&_entries[3], sizeof(Real), 1, file);
  }
  else
  {
    double entries[] = {(double)_entries[0], 
                        (double)_entries[1], 
                        (double)_entries[2], 
                        (double)_entries[3]};
    fwrite((void*)&entries[0], sizeof(double), 1, file);
    fwrite((void*)&entries[1], sizeof(double), 1, file);
    fwrite((void*)&entries[2], sizeof(double), 1, file);
    fwrite((void*)&entries[3], sizeof(double), 1, file);
  }
}

void QUATERNION::read(FILE* file)
{
  if (sizeof(Real) == sizeof(double))
  {
    fread((void*)&_entries[0], sizeof(Real), 1, file);
    fread((void*)&_entries[1], sizeof(Real), 1, file);
    fread((void*)&_entries[2], sizeof(Real), 1, file);
    fread((void*)&_entries[3], sizeof(Real), 1, file);
  }
  else
  {
    double entries[4];
    fread((void*)&entries[0], sizeof(double), 1, file);
    fread((void*)&entries[1], sizeof(double), 1, file);
    fread((void*)&entries[2], sizeof(double), 1, file);
    fread((void*)&entries[3], sizeof(double), 1, file);

    _entries[0] = entries[0];
    _entries[1] = entries[1];
    _entries[2] = entries[2];
    _entries[3] = entries[3];
  }
}

//////////////////////////////////////////////////////////////////////
// take the exponential
// from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::exp() const
{
  const Real& s =  _w;
  const VEC3F v(_x, _y, _z);

  const Real magnitude = v.magnitude();
  const Real exps = std::exp(s);
  const VEC3F vFinal = (v / magnitude) * sin(magnitude);

  return exps * QUATERNION(cos(magnitude), vFinal[0], vFinal[1], vFinal[2]);
}

//////////////////////////////////////////////////////////////////////
// take the log 
// from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::log() const
{
  const Real& s =  _w;
  const VEC3F v(_x, _y, _z);

  const Real qMagnitude = magnitude();
  const Real vMagnitude = v.magnitude();

  const VEC3F vNormalized = (vMagnitude > 0) ? v / vMagnitude : VEC3F(0,0,0);

  const VEC3F vFinal = vNormalized * acos(s / qMagnitude);

  return QUATERNION(std::log(qMagnitude), vFinal[0], vFinal[1], vFinal[2]);
}

//////////////////////////////////////////////////////////////////////
// take the power
// from: http://www.lce.hut.fi/~ssarkka/pub/quat.pdf
//
// some care has been taken to optimize this
// hacking this back to always do single precision doesn't 
// win very much
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::pow(const Real& exponent) const
{
  const Real partial = _x * _x + _y * _y + _z * _z;
  const Real qMagnitude = sqrt(partial + _w * _w);
  const Real vMagnitude = sqrt(partial);
  const Real vMagnitudeInv = (vMagnitude > 0.0) ? 1.0 / vMagnitude : 0.0;

  const Real scale = exponent * acos(_w / qMagnitude) * vMagnitudeInv;

  const Real magnitude = scale * vMagnitude;
  const Real magnitudeInv = (magnitude > 0.0) ? 1.0 / magnitude : 0.0;

  const Real exps = std::exp(exponent * std::log(qMagnitude));
  Real sMag,cMag;
#if __APPLE__
  sMag = sin(magnitude);
  cMag = cos(magnitude);
#else
  // only supported on Linux
  sincos(magnitude, &sMag, &cMag);
#endif
  const Real scale2 = scale * exps * magnitudeInv * sMag;
  return QUATERNION(exps * cMag,
                    scale2 * _x,
                    scale2 * _y,
                    scale2 * _z);
}

MATRIX3 QUATERNION::toRotationMatrix() const
{
  MATRIX3 final;
  final(0,0) = 1 - 2 * _y * _y - 2 * _z * _z; 
  final(0,1) = 2 * _x * _y - 2 * _w * _z;      
  final(0,2) = 2 * _x * _z + 2 * _w * _y;
  
  final(1,0) = 2 * _x * _y + 2 * _w * _z;      
  final(1,1) = 1 - 2 * _x * _x - 2 * _z * _z; 
  final(1,2) = 2 * _y * _z - 2 * _w * _x;

  final(2,0) = 2 * _x * _z - 2 * _w * _y;      
  final(2,1) = 2 * _y * _z + 2 * _w * _x;      
  final(2,2) = 1 - 2 * _x * _x - 2 * _y * _y;

  return final;
}

//////////////////////////////////////////////////////////////////////
// get the axis-angle representation
// as seen in David Baraff's Physically Based Modelling notes
//////////////////////////////////////////////////////////////////////
void QUATERNION::axisAngle(VEC3F& axis, Real& angle)
{
  if ((_w >= ((float)1)) || (_w <= (float)(-1)))
  {
    // identity; this check is necessary to avoid problems with acos if s is 1 + eps
    angle = 0;
    axis[0] = 1;
    axis[1] = 0;
    axis[2] = 0;
    return;
  }

  angle = 2.0 * acos(_w);
  float sin2 = _x*_x + _y*_y + _z*_z; //sin^2(*angle / 2.0)

  if (sin2 == 0)
  {
    // identity rotation; angle is zero, any axis is equally good
    axis[0] = 1;
    axis[1] = 0;
    axis[2] = 0;
  }
  else
  {
    float inv = 1.0 / sqrt(sin2); // note: *angle / 2.0 is on [0,pi], so sin(*angle / 2.0) >= 0, and therefore the sign of sqrt can be safely taken positive
    axis[0] = _x * inv;
    axis[1] = _y * inv;
    axis[2] = _z * inv;
  }

  angle = angle / (2.0 * M_PI) * 360.0f;
}
