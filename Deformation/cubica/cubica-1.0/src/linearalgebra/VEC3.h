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
#ifndef _VEC3_H_
#define _VEC3_H_

//////////////////////////////////////////////////////////////////////
// Credit: adapted from libgfx by Michael Garland
//////////////////////////////////////////////////////////////////////

#include <SETTINGS.h>
#include <iostream>
#include <math.h>
#include <VECTOR.h>

template<class T>
class TVEC3 {
private:
  T element[3];

public:
  // Standard constructors
  TVEC3(T s=0) { *this = s; }
  TVEC3(T x, T y, T z) { element[0]=x; element[1]=y; element[2]=z; }
  TVEC3(const VECTOR& vec) { assert(vec.size() == 3);
    element[0] = vec[0]; element[1] = vec[1]; element[2] = vec[2]; }

  // Copy constructors & assignment operators
  template<class U> TVEC3(const TVEC3<U>& v) { *this = v; }
  template<class U> TVEC3(const U v[3])
    { element[0]=v[0]; element[1]=v[1]; element[2]=v[2]; }
  template<class U> TVEC3& operator=(const TVEC3<U>& v)
    { element[0]=v[0];  element[1]=v[1];  element[2]=v[2];  return *this; }
  TVEC3& operator=(T s) { element[0]=element[1]=element[2]=s; return *this; }
  TVEC3& operator=(const VECTOR& vec);

  // Descriptive interface
  typedef T value_type;
  static int dim() { return 3; }

  // Access methods
  operator       T*()       { return element; }
  operator const T*() const { return element; }

  T& operator[](int i)       { return element[i]; }
  T  operator[](int i) const { return element[i]; }
  operator const T*()        { return element; }

  // Assignment and in-place arithmetic methods
  inline TVEC3& operator+=(const TVEC3& v);
  inline TVEC3& operator-=(const TVEC3& v);
  inline TVEC3& operator*=(T s);
  inline TVEC3& operator/=(T s);

  void normalize() {
    T l = norm2(*this);
    if( l!=1.0 && l!=0.0 )  *this /= sqrt(l);
  };
  void clear() {
    T zero = 0.0;
    element[0] = zero;
    element[1] = zero;
    element[2] = zero;
  };

  T maxElement() {
    T foundMax = element[0];
    if (element[1] > foundMax) foundMax = element[1];
    if (element[2] > foundMax) foundMax = element[2];
    return foundMax;
  };

  VECTOR toVector() const {
    VECTOR final(3);
    final[0] = element[0]; final[1] = element[1]; final[2] = element[2];
    return final;
  };

  void equals(const TVEC3& v) {
    element[0] = v.element[0]; element[1] = v.element[1]; element[2] = v.element[2];
  };
  void axpy(const T scalar, const TVEC3& v) {
    element[0] += scalar * v.element[0]; element[1] += scalar * v.element[1]; element[2] += scalar * v.element[2]; 
  };
  void clearingAxpy(const T scalar, const TVEC3& v) {
    element[0] = scalar * v.element[0]; element[1] = scalar * v.element[1]; element[2] = scalar * v.element[2]; 
  };
};

////////////////////////////////////////////////////////////////////////
// Method definitions
////////////////////////////////////////////////////////////////////////

template<class T> inline TVEC3<T>& TVEC3<T>::operator+=(const TVEC3<T>& v)
  { element[0] += v[0];   element[1] += v[1];   element[2] += v[2];  return *this; }

template<class T> inline TVEC3<T>& TVEC3<T>::operator-=(const TVEC3<T>& v)
  { element[0] -= v[0];   element[1] -= v[1];   element[2] -= v[2];  return *this; }

template<class T> inline TVEC3<T>& TVEC3<T>::operator*=(T s)
  { element[0] *= s;   element[1] *= s;   element[2] *= s;  return *this; }

template<class T> inline TVEC3<T>& TVEC3<T>::operator/=(T s)
  { element[0] /= s;   element[1] /= s;   element[2] /= s;  return *this; }

template<class T> inline TVEC3<T>& TVEC3<T>::operator=(const VECTOR& vec)
  { element[0] = vec[0];   element[1] = vec[1];   element[2] = vec[2];  return *this; }

////////////////////////////////////////////////////////////////////////
// Operator definitions
////////////////////////////////////////////////////////////////////////

template<class T>
inline TVEC3<T> operator+(const TVEC3<T> &u, const TVEC3<T>& v)
  { return TVEC3<T>(u[0]+v[0], u[1]+v[1], u[2]+v[2]); }

template<class T>
inline TVEC3<T> operator-(const TVEC3<T> &u, const TVEC3<T>& v)
  { return TVEC3<T>(u[0]-v[0], u[1]-v[1], u[2]-v[2]); }

template<class T> inline TVEC3<T> operator-(const TVEC3<T> &v)
  { return TVEC3<T>(-v[0], -v[1], -v[2]); }

template<class T> inline TVEC3<T> operator*(T s, const TVEC3<T> &v)
  { return TVEC3<T>(v[0]*s, v[1]*s, v[2]*s); }
template<class T> inline TVEC3<T> operator*(const TVEC3<T> &v, T s)
  { return s*v; }

template<class T> inline TVEC3<T> operator/(const TVEC3<T> &v, T s)
  { return TVEC3<T>(v[0]/s, v[1]/s, v[2]/s); }

template<class T> inline T operator*(const TVEC3<T> &u, const TVEC3<T>& v)
  { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }

template<class T> inline TVEC3<T> cross(const TVEC3<T>& u, const TVEC3<T>& v)
  { return TVEC3<T>( u[1]*v[2] - v[1]*u[2], -u[0]*v[2] + v[0]*u[2], u[0]*v[1] - v[0]*u[1] ); }

template<class T>
inline TVEC3<T> operator^(const TVEC3<T>& u, const TVEC3<T>& v)
  { return cross(u, v); }


template<class T>
inline std::ostream &operator<<(std::ostream &out, const TVEC3<T>& v)
  { return out << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")"; }

template<class T>
inline std::istream &operator>>(std::istream &in, TVEC3<T>& v)
  { return in >> v[0] >> v[1] >> v[2]; }

////////////////////////////////////////////////////////////////////////
// Misc. function definitions
////////////////////////////////////////////////////////////////////////

template<class T> inline T norm2(const TVEC3<T>& v)  { return v*v; }
template<class T> inline T norm(const TVEC3<T>& v)   { return sqrt(norm2(v)); }

template<class T> inline T infty_norm(const TVEC3<T>& v)   
{
  T result = 0;
  for (int i = 0; i < 3; i++)
  {
    T absVal = (v[i] >= 0) ? v[i] : -v[i];
    result =  (result >= absVal) ? result : absVal;
  }
  return result;
}

template<class T> inline TVEC3<T> fabs(const TVEC3<T>& v)   
{
  TVEC3<T> result;
  result[0] = fabs(v[0]);
  result[1] = fabs(v[1]);
  result[2] = fabs(v[2]);
  return result;
}

////////////////////////////////////////////////////////////////////////
// TODO: wtf is with max/min being taken.  should resolve this issue.
////////////////////////////////////////////////////////////////////////
template<class T> inline TVEC3<T> vec_max (const TVEC3<T>& u, const TVEC3<T>& v)   
{
  TVEC3<T> result;
  for (int i = 0; i < 3; i++)
  {
    result[i] = (u[i] > v[i]) ? u[i] : v[i];
  }
  return result;
}

template<class T> inline TVEC3<T> vec_min (const TVEC3<T>& u, const TVEC3<T>& v)   
{
  TVEC3<T> result;
  for (int i = 0; i < 3; i++)
  {
    result[i] = (u[i] < v[i]) ? u[i] : v[i];
  }
  return result;
}

template<class T> inline void unitize(TVEC3<T>& v)
{
    T l = norm2(v);
    if( l!=1.0 && l!=0.0 )  v /= sqrt(l);
}

template<class T> inline TVEC3<T> project_onto (const TVEC3<T>& vec_to_project, const TVEC3<T> vec_to_project_onto)
{
    return ( vec_to_project_onto * (vec_to_project * vec_to_project_onto / norm2 (vec_to_project_onto) ) );
}

template<class T> inline TVEC3<T> get_orthonormal_vector (TVEC3<T> v)
{
  // make v be a unit vector
  unitize (v);

  // let a be a vector not parallel to v
  TVEC3<T> a (1, 0, 0);
  if (a * v > 0.75)
  {
      a = TVEC3<T> (0, 1, 0);
  }

  // remove from a its projection onto v
  TVEC3<T> aOnV = v * (v * a);
  a -= aOnV;

  // unitize a
  unitize (a);
  return a;
}

typedef TVEC3<Real> VEC3;
typedef TVEC3<Real> VEC3F;

#endif

