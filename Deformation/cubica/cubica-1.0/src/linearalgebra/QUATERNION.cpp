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
#include "QUATERNION.h"

QUATERNION::QUATERNION() :
  _x(0), _y(0), _z(0), _w(1)
{
}

QUATERNION::QUATERNION(Real x, Real y, Real z, Real w) :
  _x(x), _y(y), _z(z), _w(w)
{
}

QUATERNION::QUATERNION(const VEC3F& vector) :
  _x(vector[0]), _y(vector[1]), _z(vector[2]), _w(0)
  //_x(vector[0]), _y(vector[1]), _z(vector[2]), _w(1)
{
}

QUATERNION::QUATERNION(const VECTOR& vector) 
{
  assert(vector.size() >= 4);
  _x = vector[0];
  _y = vector[1];
  _z = vector[2];
  _w = vector[3];
}

//////////////////////////////////////////////////////////////////////
// Stolen form Jernej Barbic's implemenation, which is from David
// Baraff's implementaion
//////////////////////////////////////////////////////////////////////
QUATERNION::QUATERNION(const MATRIX3& matrix)
{
  Real trace, u;
  trace = matrix(0,0) + matrix(1,1) + matrix(2,2);

  if (trace >= 0)
  {
    u = (Real)sqrt(trace + 1);
    _w = (Real)0.5 * u;
    u = (Real)0.5 / u;
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
        u = (Real)sqrt((matrix(0,0) - (matrix(1,1) + matrix(2,2))) + 1);
        _x = 0.5f * u;
        u = 0.5f / u;
        _y = (matrix(1,0) + matrix(0,1)) * u;
        _z = (matrix(0,2) + matrix(2,0)) * u;
        _w = (matrix(2,1) - matrix(1,2)) * u;
      break;

      case 1:
        u = (Real)sqrt((matrix(1,1) - (matrix(2,2) + matrix(0,0))) + 1);
        _y = 0.5f * u;
        u = 0.5f / u;
        _z = (matrix(2,1) + matrix(1,2)) * u;
        _x = (matrix(1,0) + matrix(0,1)) * u;
        _w = (matrix(0,2) - matrix(2,0)) * u;
      break;

      case 2:
        u = (Real)sqrt((matrix(2,2) - (matrix(0,0) + matrix(1,1))) + 1);
        _z = 0.5f * u;
        u = 0.5f / u;
        _x = (matrix(0,2) + matrix(2,0)) * u;
        _y = (matrix(2,1) + matrix(1,2)) * u;
        _w = (matrix(1,0) - matrix(0,1)) * u;
      break;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// return the matrix corresponding to the gradient
// 'which' refers to the component of the quaternion
//////////////////////////////////////////////////////////////////////
MATRIX QUATERNION::gradient(int which)
{
  assert(which >= 0 && which < 4);

  MATRIX dq(4,4);
  QUATERNION& q = *this;
  switch (which)
  {
    case 0:
      dq(0,0) =  q[0]; dq(0,1) =  q[1]; dq(0,2) =  q[2];
      dq(1,0) =  q[1]; dq(1,1) = -q[0]; dq(1,2) = -q[3];
      dq(2,0) =  q[2]; dq(2,1) =  q[3]; dq(2,2) = -q[0];
      dq(3,3) =  q[0];
      break;
    case 1:
      dq(0,0) = -q[1]; dq(0,1) =  q[0]; dq(0,2) =  q[3];
      dq(1,0) =  q[0]; dq(1,1) =  q[1]; dq(1,2) =  q[2];
      dq(2,0) = -q[3]; dq(2,1) =  q[2]; dq(2,2) = -q[1];
      dq(3,3) =  q[1];
      break;
    case 2:
      dq(0,0) = -q[2]; dq(0,1) = -q[3]; dq(0,2) =  q[0];
      dq(1,0) =  q[3]; dq(1,1) = -q[2]; dq(1,2) =  q[1];
      dq(2,0) =  q[0]; dq(2,1) =  q[1]; dq(2,2) =  q[2];
      dq(3,3) =  q[2];
      break;
    case 3:
      dq(0,0) =  q[3]; dq(0,1) = -q[2]; dq(0,2) =  q[1];
      dq(1,0) =  q[2]; dq(1,1) =  q[3]; dq(1,2) = -q[0];
      dq(2,0) = -q[1]; dq(2,1) =  q[0]; dq(2,2) =  q[3];
      dq(3,3) =  q[3];
      break;
  }
  dq *= 2;
  return dq;
}

//////////////////////////////////////////////////////////////////////
// return the matrix corresponding to the gradient
// 'which' refers to the component of the quaternion
//////////////////////////////////////////////////////////////////////
MATRIX3 QUATERNION::gradient3x3(int which)
{
  assert(which >= 0 && which < 4);

  MATRIX3 dq;
  QUATERNION& q = *this;
  switch (which)
  {
    case 0:
      dq(0,0) =  q[0]; dq(0,1) =  q[1]; dq(0,2) =  q[2];
      dq(1,0) =  q[1]; dq(1,1) = -q[0]; dq(1,2) = -q[3];
      dq(2,0) =  q[2]; dq(2,1) =  q[3]; dq(2,2) = -q[0];
      break;
    case 1:
      dq(0,0) = -q[1]; dq(0,1) =  q[0]; dq(0,2) =  q[3];
      dq(1,0) =  q[0]; dq(1,1) =  q[1]; dq(1,2) =  q[2];
      dq(2,0) = -q[3]; dq(2,1) =  q[2]; dq(2,2) = -q[1];
      break;
    case 2:
      dq(0,0) = -q[2]; dq(0,1) = -q[3]; dq(0,2) =  q[0];
      dq(1,0) =  q[3]; dq(1,1) = -q[2]; dq(1,2) =  q[1];
      dq(2,0) =  q[0]; dq(2,1) =  q[1]; dq(2,2) =  q[2];
      break;
    case 3:
      dq(0,0) =  q[3]; dq(0,1) = -q[2]; dq(0,2) =  q[1];
      dq(1,0) =  q[2]; dq(1,1) =  q[3]; dq(1,2) = -q[0];
      dq(2,0) = -q[1]; dq(2,1) =  q[0]; dq(2,2) =  q[3];
      break;
  }
  dq *= 2;
  return dq;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX QUATERNION::hessian(int i, int j)
{
  assert(i >= 0 && i < 4);
  assert(j >= 0 && j < 4);

  MATRIX final(4,4);

  //dxdx =
  if (i == 0 && j == 0)
  {
    final(0,0) = 1;
    final(1,1) = -1;
    final(2,2) = -1;
    final(3,3) = 1;
  }

  //dydy = 
  if (i == 1 && j == 1)
  {
    final(0,0) = -1;
    final(1,1) = 1;
    final(2,2) = -1;
    final(3,3) = 1;
  }

  //dzdz = 
  if (i == 2 && j == 2)
  {
    final(0,0) = -1;
    final(1,1) = -1;
    final(2,2) = 1;
    final(3,3) = 1;
  }

  //dwdw = 
  if (i == 3 && j == 3)
  {
    final(0,0) = 1;
    final(1,1) = 1;
    final(2,2) = 1;
    final(3,3) = 1;
  }

  //dxdy = 
  if ((i == 0 && j == 1) || (i == 1 && j == 0))
  {
    final(0,1) = 1;
    final(1,0) = 1;
  }

  //dxdz = 
  if ((i == 0 && j == 2) || (i == 2 && j == 0))
  {
    final(0,2) = 1;
    final(2,0) = 1;
  }

  //dxdw= 
  if ((i == 0 && j == 3) || (i == 3 && j == 0))
  {
    final(1,2) = -1;
    final(2,1) = 1;
  }

  //dydz= 
  if ((i == 1 && j == 2) || (i == 2 && j == 1))
  {
    final(1,2) = 1;
    final(2,1) = 1;
  }

  //dydw= 
  if ((i == 1 && j == 3) || (i == 3 && j == 1))
  {
    final(0,2) = 1;
    final(2,0) = -1;
  }

  //dzdw= 
  if ((i == 2 && j == 3) || (i == 3 && j == 2))
  {
    final(0,1) = -1;
    final(1,0) = 1;
  }

  final *= 2;
  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX3 QUATERNION::hessian3x3(int i, int j)
{
  assert(i >= 0 && i < 4);
  assert(j >= 0 && j < 4);

  MATRIX3 final;

  //dxdx =
  if (i == 0 && j == 0)
  {
    final(0,0) = 1;
    final(1,1) = -1;
    final(2,2) = -1;
  }

  //dydy = 
  if (i == 1 && j == 1)
  {
    final(0,0) = -1;
    final(1,1) = 1;
    final(2,2) = -1;
  }

  //dzdz = 
  if (i == 2 && j == 2)
  {
    final(0,0) = -1;
    final(1,1) = -1;
    final(2,2) = 1;
  }

  //dwdw = 
  if (i == 3 && j == 3)
  {
    final(0,0) = 1;
    final(1,1) = 1;
    final(2,2) = 1;
  }

  //dxdy = 
  if ((i == 0 && j == 1) || (i == 1 && j == 0))
  {
    final(0,1) = 1;
    final(1,0) = 1;
  }

  //dxdz = 
  if ((i == 0 && j == 2) || (i == 2 && j == 0))
  {
    final(0,2) = 1;
    final(2,0) = 1;
  }

  //dxdw= 
  if ((i == 0 && j == 3) || (i == 3 && j == 0))
  {
    final(1,2) = -1;
    final(2,1) = 1;
  }

  //dydz= 
  if ((i == 1 && j == 2) || (i == 2 && j == 1))
  {
    final(1,2) = 1;
    final(2,1) = 1;
  }

  //dydw= 
  if ((i == 1 && j == 3) || (i == 3 && j == 1))
  {
    final(0,2) = 1;
    final(2,0) = -1;
  }

  //dzdw= 
  if ((i == 2 && j == 3) || (i == 3 && j == 2))
  {
    final(0,1) = -1;
    final(1,0) = 1;
  }

  final *= 2;
  return final;
}

QUATERNION QUATERNION::conjugate()
{
  QUATERNION final(-_x, -_y, -_z, _w);
  return final;
}

QUATERNION& QUATERNION::operator=(const QUATERNION q)
{
  this->_x = q._x;
  this->_y = q._y;
  this->_z = q._z;
  this->_w = q._w;
  return *this;
}

QUATERNION& QUATERNION::operator=(const VEC3F& vec)
{
  this->_x = vec[0];
  this->_y = vec[1];
  this->_z = vec[2];
  this->_w = 0.0;
  return *this;
}

QUATERNION& QUATERNION::operator*=(const Real r)
{
  this->_x *= r; 
  this->_y *= r; 
  this->_z *= r; 
  this->_w *= r; 
  return *this;
}

QUATERNION& QUATERNION::operator+=(const QUATERNION& q)
{
  this->_x += q.x(); 
  this->_y += q.y(); 
  this->_z += q.z(); 
  this->_w += q.w(); 
  return *this;
}

QUATERNION& QUATERNION::operator-=(const QUATERNION& q)
{
  this->_x -= q.x(); 
  this->_y -= q.y(); 
  this->_z -= q.z(); 
  this->_w -= q.w(); 
  return *this;
}

QUATERNION operator-(const QUATERNION& left, const QUATERNION& right)
{
  QUATERNION final(left._x - right._x,
                   left._y - right._y,
                   left._z - right._z,
                   left._w - right._w);
  return final;
}

QUATERNION operator+(const QUATERNION& left, const QUATERNION& right)
{
  QUATERNION final(left._x + right._x,
                   left._y + right._y,
                   left._z + right._z,
                   left._w + right._w);
  return final;
}

QUATERNION operator*(const MATRIX& left, const QUATERNION& right)
{
  assert(left.rows() == 4 && left.cols() == 4);

  QUATERNION final;
  final[0] = left(0,0) * right[0] + left(0,1) * right[1] + left(0,2) * right[2] + left(0,3) * right[3];
  final[1] = left(1,0) * right[0] + left(1,1) * right[1] + left(1,2) * right[2] + left(1,3) * right[3];
  final[2] = left(2,0) * right[0] + left(2,1) * right[1] + left(2,2) * right[2] + left(2,3) * right[3];
  final[3] = left(3,0) * right[0] + left(3,1) * right[1] + left(3,2) * right[2] + left(3,3) * right[3];

  return final;
}

QUATERNION operator*(const MATRIX3& left, const QUATERNION& right)
{
  assert(fabs(right[3]) < 1e-7);

  QUATERNION final;
  final[0] = left(0,0) * right[0] + left(0,1) * right[1] + left(0,2) * right[2];
  final[1] = left(1,0) * right[0] + left(1,1) * right[1] + left(1,2) * right[2];
  final[2] = left(2,0) * right[0] + left(2,1) * right[1] + left(2,2) * right[2];
  final[3] = 0.0;

  return final;
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

Real operator^(const QUATERNION& left, const QUATERNION& right)
{
  return left._x * right._x + left._y * right._y +
         left._z * right._z + left._w * right._w;
}

MATRIX3 QUATERNION::toRotationMatrix()
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

ostream& operator<<(ostream &out, const QUATERNION& q)
{
  out << q._x << " " << q._y << " " << q._z << " " << q._w << " ";
  return out;
}

QUATERNION operator*(const QUATERNION& left, const Real& right)
{
  return QUATERNION(left.x() * right, left.y() * right , left.z() * right, left.w() * right);
}

QUATERNION operator*(const Real& left, const QUATERNION& right)
{
  return QUATERNION(left * right.x(), left * right.y(), left * right.z(), left * right.w());
}

void QUATERNION::normalize()
{
  Real invMagnitude = 1.0 / sqrtf(_x * _x + _y * _y + _z * _z + _w * _w);
  _x *= invMagnitude;
  _y *= invMagnitude;
  _z *= invMagnitude;
  _w *= invMagnitude;
}

MATRIX QUATERNION::toExplicitMatrix()
{
  MATRIX final(4,4);

  final(0,0) = -_z * _z  - _y * _y + _w * _w + _x * _x; 
  final(0,1) = - 2* _z *_w + 2* _y *_x;
  final(0,2) = + 2* _z *_x + 2* _y *_w;

  final(1,0) = + 2* _x *_y + 2* _z *_w;
  final(1,1) = -_x * _x - _z * _z + _w * _w + _y * _y;
  final(1,2) = - 2* _x *_w + 2* _z *_y;

  final(2,0) = - 2* _y *_w + 2* _x *_z;
  final(2,1) = + 2* _y *_z + 2* _x *_w;
  final(2,2) = -_y * _y + _w * _w + _z * _z - _x * _x;

  final(3,3) = _w * _w + _x * _x  + _y * _y  + _z * _z;

  return final; 
}

MATRIX3 QUATERNION::toExplicitMatrix3x3()
{
  MATRIX3 final;

  final(0,0) = -_z * _z  - _y * _y + _w * _w + _x * _x; 
  final(0,1) = - 2.0* _z *_w + 2.0* _y *_x;
  final(0,2) = + 2.0* _z *_x + 2.0* _y *_w;

  final(1,0) = + 2.0* _x *_y + 2.0* _z *_w;
  final(1,1) = -_x * _x - _z * _z + _w * _w + _y * _y;
  final(1,2) = - 2.0* _x *_w + 2.0* _z *_y;

  final(2,0) = - 2.0* _y *_w + 2.0* _x *_z;
  final(2,1) = + 2.0* _y *_z + 2.0* _x *_w;
  final(2,2) = -_y * _y + _w * _w + _z * _z - _x * _x;

  return final; 
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
// get the axis-angle representation
// as seen in David Baraff's Physically Based Modelling notes
//////////////////////////////////////////////////////////////////////
void QUATERNION::axisAngle(VEC3F& axis, Real& angle)
{
  if ((_w >= ((Real)1)) || (_w <= (Real)(-1)))
  {
    // identity; this check is necessary to avoid problems with acos if s is 1 + eps
    angle = 0;
    axis[0] = 1;
    axis[1] = 0;
    axis[2] = 0;
    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //cout << "Quaternion is not a rotation! " << endl;
    return;
  }

  angle = 2.0 * acos(_w);
  Real sin2 = _x*_x + _y*_y + _z*_z; //sin^2(*angle / 2.0)

  if (sin2 == 0)
  {
    // identity rotation; angle is zero, any axis is equally good
    axis[0] = 1;
    axis[1] = 0;
    axis[2] = 0;
  }
  else
  {
    Real inv = 1.0 / sqrt(sin2); // note: *angle / 2.0 is on [0,pi], so sin(*angle / 2.0) >= 0, and therefore the sign of sqrt can be safely taken positive
    axis[0] = _x * inv;
    axis[1] = _y * inv;
    axis[2] = _z * inv;
  }

  angle = angle / (2.0 * M_PI) * 360.0f;
}

//////////////////////////////////////////////////////////////////////
// rotate the passed in vector by this
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::rotates(QUATERNION& toRotate)
{
  QUATERNION inv = this->conjugate();

  // can't assume it's a unit quaternion
  //
  // NOTE: non-unit magnitude problems should be handled outside this
  // class. Secretly normalizing inside here only serves to screw up
  // the gradient descent search when solving for the domain rotations.
  //inv *= 1.0 / (_x * _x + _y * _y + _z * _z + _w * _w);

  return (*this) * toRotate * inv;
}

//////////////////////////////////////////////////////////////////////
// rotate the passed in vector by this
//////////////////////////////////////////////////////////////////////
VEC3F QUATERNION::rotates(VEC3F& toRotate)
{
  QUATERNION inv = this->conjugate();

  // can't assume it's a unit quaternion
  //
  // NOTE: non-unit magnitude problems should be handled outside this
  // class. Secretly normalizing inside here only serves to screw up
  // the gradient descent search when solving for the domain rotations.
  //inv *= 1.0 / (_x * _x + _y * _y + _z * _z + _w * _w);

  QUATERNION qRotate(toRotate);
  QUATERNION rotated = (*this) * qRotate * inv;

  VEC3F final;
  final[0] = rotated[0];
  final[1] = rotated[1];
  final[2] = rotated[2];

  return final;
}

//////////////////////////////////////////////////////////////////////
// rotate the passed in 3x3 matrix
//////////////////////////////////////////////////////////////////////
MATRIX QUATERNION::rotates(MATRIX& toRotate)
{
  assert(toRotate.rows() == 3);
  assert(toRotate.cols() == 3);

  MATRIX final(3,3);
  for (int x = 0; x < 3; x++)
  {
    VEC3F column;
    column[0] = toRotate(0,x);
    column[1] = toRotate(1,x);
    column[2] = toRotate(2,x);

    VEC3F rotated = this->rotates(column);

    final(0,x) = rotated[0];
    final(1,x) = rotated[1];
    final(2,x) = rotated[2];
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// populate a VECTOR of size 4
//////////////////////////////////////////////////////////////////////
VECTOR QUATERNION::toVector()
{
  VECTOR final(4);
  final[0] = _entries[0];
  final[1] = _entries[1];
  final[2] = _entries[2];
  final[3] = _entries[3];
  return final;
}

//////////////////////////////////////////////////////////////////////
// populate a VECTOR of size 5 and put lambda in the last entry
//////////////////////////////////////////////////////////////////////
VECTOR QUATERNION::toVector(Real lambda)
{
  VECTOR final(5);
  final[0] = _entries[0];
  final[1] = _entries[1];
  final[2] = _entries[2];
  final[3] = _entries[3];
  final[4] = lambda;
  return final;
}

//////////////////////////////////////////////////////////////////////
// 2 norm of the entries
//////////////////////////////////////////////////////////////////////
Real QUATERNION::magnitude()
{
  Real final = _entries[0] * _entries[0];
  final += _entries[1] * _entries[1];
  final += _entries[2] * _entries[2];
  final += _entries[3] * _entries[3];

  return sqrt(final);
}

//////////////////////////////////////////////////////////////////////
// Shabana, equation 2.21 (note he packs the real component first)
//////////////////////////////////////////////////////////////////////
MATRIX QUATERNION::E()
{
  MATRIX final(3, 4);
  final(0,0) =  _entries[3];
  final(1,0) =  _entries[2];
  final(2,0) = -_entries[1];

  final(0,1) = -_entries[2];
  final(1,1) =  _entries[3];
  final(2,1) =  _entries[0];

  final(0,2) =  _entries[1];
  final(1,2) = -_entries[0];
  final(2,2) =  _entries[3];

  final(0,3) = -_entries[0];
  final(1,3) = -_entries[1];
  final(2,3) = -_entries[2];
  return final;
}

//////////////////////////////////////////////////////////////////////
// Shabana, equation 2.22 (note he packs the real component first)
//////////////////////////////////////////////////////////////////////
MATRIX QUATERNION::Ebar()
{
  MATRIX final(3, 4);
  final(0,0) =  _entries[3];
  final(1,0) = -_entries[2];
  final(2,0) =  _entries[1];

  final(0,1) =  _entries[2];
  final(1,1) =  _entries[3];
  final(2,1) = -_entries[0];

  final(0,2) = -_entries[1];
  final(1,2) =  _entries[0];
  final(2,2) =  _entries[3];

  final(0,3) = -_entries[0];
  final(1,3) = -_entries[1];
  final(2,3) = -_entries[2];
  return final;
}

//////////////////////////////////////////////////////////////////////
// Shabana, equation 2.59 (note he packs the real component first)
//////////////////////////////////////////////////////////////////////
MATRIX QUATERNION::G()
{
  MATRIX final(3, 4);
  final(0,0) =  _entries[3];
  final(1,0) =  _entries[2];
  final(2,0) = -_entries[1];

  final(0,1) = -_entries[2];
  final(1,1) =  _entries[3];
  final(2,1) =  _entries[0];

  final(0,2) =  _entries[1];
  final(1,2) = -_entries[0];
  final(2,2) =  _entries[3];

  final(0,3) = -_entries[0];
  final(1,3) = -_entries[1];
  final(2,3) = -_entries[2];
  final *= 2.0;
  return final;
}

//////////////////////////////////////////////////////////////////////
// Shabana, equation 2.60 (note he packs the real component first)
//////////////////////////////////////////////////////////////////////
MATRIX QUATERNION::Gbar()
{
  MATRIX final(3, 4);
  final(0,0) =  _entries[3];
  final(1,0) = -_entries[2];
  final(2,0) =  _entries[1];

  final(0,1) =  _entries[2];
  final(1,1) =  _entries[3];
  final(2,1) = -_entries[0];

  final(0,2) = -_entries[1];
  final(1,2) =  _entries[0];
  final(2,2) =  _entries[3];

  final(0,3) = -_entries[0];
  final(1,3) = -_entries[1];
  final(2,3) = -_entries[2];
  final *= 2.0;
  return final;
}

//////////////////////////////////////////////////////////////////////
// compute an Open Dynamics Engine (ODE)-style update vector
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::odeUpdate(Real& dt, VEC3F& angularVelocity)
{
  // Try the ODE way
  Real wlen = norm(angularVelocity);
  Real h = dt * 0.5;

  //dReal wlen = dSqrt (b->avel[0]*b->avel[0] + b->avel[1]*b->avel[1] +
	//    b->avel[2]*b->avel[2]);
  
  Real theta = wlen * h;
  //h *= REAL(0.5);
  //dReal theta = wlen * h;
  
  QUATERNION q;
  q[3] = cos(theta);
  //q[0] = dCos(theta);
  
  Real sinc;
  if (fabs(theta) < 1.0e-4) 
    sinc = Real(1.0) - theta*theta*Real(0.166666666666666666667);
  else 
    sinc = sin(theta)/theta;
  Real s = sinc * h;

  q[0] = angularVelocity[0] * s;
  q[1] = angularVelocity[1] * s;
  q[2] = angularVelocity[2] * s;

  return q;
}

//////////////////////////////////////////////////////////////////////
// Convert from axis-angle representation to quaternion
//
// From Simo and Wong, 
// "Unconditionally Stable Algorithms for Rigid Body Dynamics that
// Exactly Preserve Energy and Momentum", International Journal for
// Numerical Methods in Engineering, 31, 19--52 (1991)
//
// Appendix B
//////////////////////////////////////////////////////////////////////
QUATERNION QUATERNION::fromAxisAngle(VEC3F& axisAngle)
{
  /*
  // Try the ODE way
  Real wlen = norm(axisAngle);
      
  //dReal wlen = dSqrt (b->avel[0]*b->avel[0] + b->avel[1]*b->avel[1] +
	//    b->avel[2]*b->avel[2]);
  
  Real theta = wlen * 0.5;
  //h *= REAL(0.5);
  //dReal theta = wlen * h;
  
  QUATERNION q;
  q[3] = cos(theta);
  //q[0] = dCos(theta);
  
  Real sinc;
  if (fabs(theta) < 1.0e-4) 
    sinc = Real(1.0) - theta*theta*Real(0.166666666666666666667);
  else 
    sinc = sin(theta)/theta;
  //dReal s = sinc(theta) * h;

  q[0] = axisAngle[0] * sinc;
  q[1] = axisAngle[1] * sinc;
  q[2] = axisAngle[2] * sinc;

  //q[1] = b->avel[0] * s;
  //q[2] = b->avel[1] * s;
  //q[3] = b->avel[2] * s;

  return q;
  */

  // if it's too small, use Taylor series
  Real magnitude = norm(axisAngle);
  Real half = 0.5 * magnitude;

  // this appears to be slightly stabler, but the error
  // component in the quaternion is still growing
  VEC3F imaginary = 0.5 * (sin(half) / (half)) * axisAngle;

  //VEC3F normalized = axisAngle;
  //normalized.normalize();
  //VEC3F imaginary = sin(half) * normalized;
 
  //if (magnitude < 1e-15)
  //if (magnitude < 1e-8)
  if (magnitude < 1e-6) // original setting
  //if (magnitude < 1e-4) // original setting
  {
    // do a Taylor approximation
    Real squared = half * half;
    Real fourth = squared * squared;
    Real sixth = fourth * squared;
    Real eight = fourth * fourth;
    Real ten = eight * squared;
    Real taylor = 1.0 - squared / 6.0 + fourth / 120.0 - sixth / 5040.0 + eight / 362880.0 - ten / 39916800.0;

    imaginary = 0.5 * taylor * axisAngle;
    //cout << " Last term: " << eight / 362880.0 << endl;
    //cout << " magnitude: " << magnitude << endl;
    //cout << " imaginary: " << imaginary << endl;
  }

  // Simo and Wong, Appendix B
  QUATERNION final;
  final[0] = imaginary[0];
  final[1] = imaginary[1];
  final[2] = imaginary[2];
  final[3] = cos(0.5 * magnitude);

  return final;

  /*
  // try using the trig form
  Real magnitude = norm(axisAngle);
  VEC3F omega = axisAngle;
  omega.normalize();

  Real real = cos(magnitude * 0.5);
  QUATERNION final(sin(magnitude * 0.5) * omega);
  final[3] = real;

  return final;
  */

  /*
  // try a Cayley map
  MATRIX3 directCayley= MATRIX3::cayley(axisAngle);
  return QUATERNION(directCayley);
  */

  /*
  // try computing the exponential directly
  MATRIX3 directExp = MATRIX3::exp(axisAngle);
  return QUATERNION(directExp);
  */
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
