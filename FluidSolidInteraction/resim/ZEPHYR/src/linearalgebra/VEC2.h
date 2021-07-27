/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2013 Theodore Kim
 */
#ifndef _VEC2_H_
#define _VEC2_H_

// Credit: adapted from libgfx by Michael Garland

#include <SETTINGS.h>
#include <iostream>

template<class T>
class TVEC2 {
private:
    T elt[2];

public:
    // Standard constructors
    //
    TVEC2(T s=0) { *this = s; }
    TVEC2(T x, T y) { elt[0]=x; elt[1]=y; }

    // Copy constructors & assignment operators
    template<class U> TVEC2(const TVEC2<U>& v) { *this = v; }
    template<class U> TVEC2(const U v[2]) { elt[0]=v[0]; elt[1]=v[1]; }
    template<class U> TVEC2& operator=(const TVEC2<U>& v)
	{ elt[0]=v[0];  elt[1]=v[1];  return *this; }
    TVEC2& operator=(T s) { elt[0]=elt[1]=s; return *this; }


    // Descriptive interface
    //
    typedef T value_type;
    static int dim() { return 2; }


    // Access methods
    //
    operator       T*()       { return elt; }
    operator const T*() const { return elt; }

#ifndef HAVE_CASTING_LIMITS
    T& operator[](int i)       { return elt[i]; }
    T  operator[](int i) const { return elt[i]; }
    operator const T*()       { return elt; }
#endif

    // In-place arithmetic methods
    //
    inline TVEC2& operator+=(const TVEC2& v);
    inline TVEC2& operator-=(const TVEC2& v);
    inline TVEC2& operator*=(T s);
    inline TVEC2& operator/=(T s);
};

////////////////////////////////////////////////////////////////////////
//
// Method definitions
//
template<class T> inline TVEC2<T>& TVEC2<T>::operator+=(const TVEC2<T>& v)
	{ elt[0] += v[0];   elt[1] += v[1];   return *this; }

template<class T> inline TVEC2<T>& TVEC2<T>::operator-=(const TVEC2<T>& v)
	{ elt[0] -= v[0];   elt[1] -= v[1];   return *this; }

template<class T> inline TVEC2<T>& TVEC2<T>::operator*=(T s)
	{ elt[0] *= s;   elt[1] *= s;   return *this; }

template<class T> inline TVEC2<T>& TVEC2<T>::operator/=(T s)
	{ elt[0] /= s;   elt[1] /= s;   return *this; }

////////////////////////////////////////////////////////////////////////
//
// Operator defintions
//

template<class T>
inline TVEC2<T> operator+(const TVEC2<T> &u, const TVEC2<T> &v)
	{ return TVEC2<T>(u[0]+v[0], u[1]+v[1]); }

template<class T>
inline TVEC2<T> operator-(const TVEC2<T> &u, const TVEC2<T> &v)
	{ return TVEC2<T>(u[0]-v[0], u[1]-v[1]); }

template<class T> inline TVEC2<T> operator-(const TVEC2<T> &v)
	{ return TVEC2<T>(-v[0], -v[1]); }

#if _MSC_VER>=1200
// Normally, we use the <class T, class N> construct below to allow the scalar
// argument to be different than the template type.  This, for example, allows
// the user to write things like v/2.  Unfortunately, Microsoft VC6.0 (aka
// v1200) gets confused by this.  We used to include explicit versions for the
// case of int's, but this was causing silent (and incorrect) coercion of
// doubles to ints.
//

  template<class T> inline TVEC2<T> operator*(T s, const TVEC2<T> &v)
	{ return TVEC2<T>(v[0]*s, v[1]*s); }
  template<class T> inline TVEC2<T> operator*(const TVEC2<T> &v, T s)
	{ return s*v; }

  template<class T> inline TVEC2<T> operator/(const TVEC2<T> &v, T s)
	{ return TVEC2<T>(v[0]/s, v[1]/s); }

#else
  template<class T, class N> inline TVEC2<T> operator*(N s, const TVEC2<T> &v)
	{ return TVEC2<T>(v[0]*s, v[1]*s); }
  template<class T, class N> inline TVEC2<T> operator*(const TVEC2<T> &v, N s)
	{ return s*v; }

  template<class T, class N> inline TVEC2<T> operator/(const TVEC2<T> &v, N s)
	{ return TVEC2<T>(v[0]/s, v[1]/s); }
#endif

template<class T> inline T operator*(const TVEC2<T> &u, const TVEC2<T>& v)
	{ return u[0]*v[0] + u[1]*v[1]; }

template<class T> inline TVEC2<T> perp(const TVEC2<T> &v)
	{ return TVEC2<T>(v[1], -v[0]); }

template<class T>
inline std::ostream &operator<<(std::ostream &out, const TVEC2<T> &v)
	{ return out << v[0] << " " << v[1]; }

template<class T>
inline std::istream &operator>>(std::istream &in, TVEC2<T>& v)
	{ return in >> v[0] >> v[1]; }


////////////////////////////////////////////////////////////////////////
//
// Misc. function definitions
//

template<class T> inline T norm2(const TVEC2<T>& v)  { return v*v; }
template<class T> inline T norm(const TVEC2<T>& v)   { return sqrt(norm2(v)); }

template<class T> inline T infty_norm (const TVEC2<T>& v)   
{
	T result = 0;
	for (int i = 0; i < 2; i++)
	{
		T absVal = (v[i] >= 0) ? v[i] : -v[i];
		result =  (result >= absVal) ? result : absVal;
	}
	return result;
}

/*
template<class T> inline TVEC2<T> max (const TVEC2<T>& u, const TVEC2<T>& v)   
{
	TVEC2<T> result;
	for (int i = 0; i < 2; i++)
	{
		result[i] = (u[i] > v[i]) ? u[i] : v[i];
	}
	return result;
}

template<class T> inline TVEC2<T> min (const TVEC2<T>& u, const TVEC2<T>& v)   
{
	TVEC2<T> result;
	for (int i = 0; i < 2; i++)
	{
		result[i] = (u[i] < v[i]) ? u[i] : v[i];
	}
	return result;
}
*/

template<class T> inline T cross (const TVEC2<T>& u, const TVEC2<T>& v)   
	{ return u[0]*v[1] - v[0]*u[1]; }

template<class T> inline void unitize(TVEC2<T>& v)
{
    T l = norm2(v);
    if( l!=1.0 && l!=0.0 )  v /= sqrt(l);
}

typedef TVEC2<Real> VEC2;
typedef TVEC2<Real>  VEC2F;

#endif

