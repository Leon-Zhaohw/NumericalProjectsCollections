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
#ifndef _MAT3_H_
#define _MAT3_H_

#include <VEC3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//////////////////////////////////////////////////////////////////////
// 3x3 Matrix class - adapted from libgfx by Michael Garland
//////////////////////////////////////////////////////////////////////

class MATRIX3
{
private:
  VEC3F row[3];

public:
  // Standard constructors
  MATRIX3() { *this = 0.0; }
  MATRIX3(const VEC3F& r0,const VEC3F& r1,const VEC3F& r2)
    { row[0]=r0; row[1]=r1; row[2]=r2; }
  MATRIX3(const MATRIX3& m) { *this = m; }
  MATRIX3(Real* data);

  // Access methods
  Real& operator()(int i, int j)       { return row[i][j]; }
  Real  operator()(int i, int j) const { return row[i][j]; }
  VEC3F&       operator[](int i)       { return row[i]; }
  const VEC3F& operator[](int i) const { return row[i]; }
  inline VEC3F col(int i) const {return VEC3F(row[0][i],row[1][i],row[2][i]);}

  operator       Real*()       { return row[0]; }
  operator const Real*()       { return row[0]; }
  operator const Real*() const { return row[0]; }

  // Assignment methods
  inline MATRIX3& operator=(const MATRIX3& m);
  inline MATRIX3& operator=(Real s);

  inline MATRIX3& operator+=(const MATRIX3& m);
  inline MATRIX3& operator-=(const MATRIX3& m);
  inline MATRIX3& operator*=(Real s);
  inline MATRIX3& operator/=(Real s);

  // Construction of standard matrices
  static MATRIX3 I();
  static MATRIX3 outer_product(const VEC3F& u, const VEC3F& v);
  static MATRIX3 outer_product(const VEC3F& v);
  static MATRIX3 cross(const VEC3F& v);

  MATRIX3& diag(Real d);
  MATRIX3& ident() { return diag(1.0); }

  MATRIX3 inverse() {
    MATRIX3 A = adjoint();
    Real d = A[0] * row[0];
    if (d == 0.0f)
      return ident();
    A = A.transpose();
    A /= d;
    return A;
  };
  MATRIX3 adjoint() {
    return MATRIX3(row[1]^row[2], row[2]^row[0], row[0]^row[1]);
  };
  MATRIX3 transpose() const {
    return MATRIX3(this->col(0), this->col(1), this->col(2));
  };
  void clear() {
    row[0].clear(); row[1].clear(); row[2].clear();
  };

  // theta is in RADIANS
  static MATRIX3 rotation(VEC3F axis, Real theta)
  {
    axis *= 1.0 / norm(axis);
    Real cosTheta = cos(theta);
    Real sinTheta = sin(theta);
    Real versine = (1.0 - cosTheta);
   
    MATRIX3 R;
    R(0,0) = axis[0] * axis[0] * versine + cosTheta;
    R(0,1) = axis[0] * axis[1] * versine - axis[2] * sinTheta;
    R(0,2) = axis[0] * axis[2] * versine + axis[1] * sinTheta;
  
    R(1,0) = axis[1] * axis[0] * versine + axis[2] * sinTheta;
    R(1,1) = axis[1] * axis[1] * versine + cosTheta;
    R(1,2) = axis[1] * axis[2] * versine - axis[0] * sinTheta;
  
    R(2,0) = axis[2] * axis[0] * versine - axis[1] * sinTheta;
    R(2,1) = axis[2] * axis[1] * versine + axis[0] * sinTheta;
    R(2,2) = axis[2] * axis[2] * versine + cosTheta;

    return R;
  };

  static MATRIX3 exp(VEC3F omega);
  static MATRIX3 dexp(VEC3F omega);
  static MATRIX3 cayley(VEC3F omega);

  Real contraction(MATRIX3& rhs) {
    Real sum = 0.0;
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
        sum += (*this)(x,y) * rhs(x,y);
    return sum;
  };
  Real squaredSum() {
    Real dot = 0;
    for (int x = 0; x < 3; x++)
      dot += row[x] * row[x];
    return dot;
  };
  Real absSum() {
    Real sum = 0.0;
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
        sum += (*this)(x,y);
    return sum;
  };

  void setToIdentity()
  {
    this->clear();
    (*this)(0,0) = 1.0;
    (*this)(1,1) = 1.0;
    (*this)(2,2) = 1.0;
  }
};

////////////////////////////////////////////////////////////////////////
// Method definitions
////////////////////////////////////////////////////////////////////////

inline MATRIX3& MATRIX3::operator=(const MATRIX3& m)
	{ row[0] = m[0]; row[1] = m[1]; row[2] = m[2];  return *this; }

inline MATRIX3& MATRIX3::operator=(Real s)
	{ row[0]=s;  row[1]=s;  row[2]=s;  return *this; }

inline MATRIX3& MATRIX3::operator+=(const MATRIX3& m)
	{ row[0] += m[0]; row[1] += m[1]; row[2] += m[2]; return *this; }

inline MATRIX3& MATRIX3::operator-=(const MATRIX3& m)
	{ row[0] -= m[0]; row[1] -= m[1]; row[2] -= m[2]; return *this; }

inline MATRIX3& MATRIX3::operator*=(Real s)
	{ row[0] *= s; row[1] *= s; row[2] *= s;  return *this; }

inline MATRIX3& MATRIX3::operator/=(Real s)
	{ row[0] /= s; row[1] /= s; row[2] /= s;  return *this; }

////////////////////////////////////////////////////////////////////////
// Operator definitions
////////////////////////////////////////////////////////////////////////

inline MATRIX3 operator+(const MATRIX3& n, const MATRIX3& m)
	{ return MATRIX3(n[0]+m[0], n[1]+m[1], n[2]+m[2]); }

inline MATRIX3 operator-(const MATRIX3& n, const MATRIX3& m)
	{ return MATRIX3(n[0]-m[0], n[1]-m[1], n[2]-m[2]); }

inline MATRIX3 operator-(const MATRIX3& m)
	{ return MATRIX3(-m[0], -m[1], -m[2]); }

inline MATRIX3 operator*(Real s, const MATRIX3& m)
	{ return MATRIX3(m[0]*s, m[1]*s, m[2]*s); }
inline MATRIX3 operator*(const MATRIX3& m, Real s)
	{ return s*m; }

inline MATRIX3 operator/(const MATRIX3& m, Real s)
	{ return MATRIX3(m[0]/s, m[1]/s, m[2]/s); }

inline VEC3F operator*(const MATRIX3& m, const VEC3F& v)
	{ return VEC3F(m[0]*v, m[1]*v, m[2]*v); }

extern MATRIX3 operator*(const MATRIX3& n, const MATRIX3& m);

inline std::ostream &operator<<(std::ostream &out, const MATRIX3& M)
	{ return out << "[" << std::endl 
               << M[0] << std::endl << M[1] << std::endl << M[2]
               << std::endl << "]" << std::endl; }

inline std::istream &operator>>(std::istream &in, MATRIX3& M)
	{ return in >> M[0] >> M[1] >> M[2]; }

////////////////////////////////////////////////////////////////////////
// Misc. function definitions
////////////////////////////////////////////////////////////////////////

inline Real det(const MATRIX3& m) { return m[0] * (m[1] ^ m[2]); }

inline Real trace(const MATRIX3& m) { return m(0,0) + m(1,1) + m(2,2); }

#endif

