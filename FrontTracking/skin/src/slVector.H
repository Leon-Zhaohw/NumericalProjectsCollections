// Copyright (c) 2011, Regents of the University of Utah
// Copyright (c) 2003-2005, Regents of the University of California.  
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the <organization> nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef SlVECTOR_is_defined
#define SlVECTOR_is_defined

//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Classes from this header:

class SlVector2;
class SlVector3;

#include <cmath>
#include <cstdlib>
#include <iostream>

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

//-------------------------------------------------------------------
//-------------------------------------------------------------------

class SlVector3 {
  // Class to store a three dimensional vector.
  // Expected uses include point in R3, RGB color, etc.
public:

  //----------------------------------------------
  // Constructors

  inline SlVector3();
  inline SlVector3(double d);
  inline SlVector3(double d0,double d1,double d2);

  inline SlVector3(const SlVector3 &da);
  inline SlVector3(const double    *da);
  // da should point to a double[3] that will be copied


  //----------------------------------------------
  // Index operators

  inline double &operator[](unsigned int i)      ;
  inline double  operator[](unsigned int i) const;

  inline double &operator()(unsigned int i)      ;
  inline double  operator()(unsigned int i) const;

  inline double x() const;
  inline double y() const;
  inline double z() const;

  //----------------------------------------------
  // Assignment and set

  inline SlVector3 &set(double d);
  inline SlVector3 &set(double d0, double d1, double d2);

  inline SlVector3 &set(const SlVector3 &da);
  inline SlVector3 &set(const double    *da);
  // da should point to a double[3] that will be copied

  inline SlVector3 &operator=(double d);
  inline SlVector3 &operator=(const SlVector3 &da);
  inline SlVector3 &operator=(const double    *da);


  //----------------------------------------------
  // Comparison operators

  inline int operator==(const SlVector3 &da) const;
  inline int operator!=(const SlVector3 &da) const;

  inline int operator==(double d) const;
  inline int operator!=(double d) const;


  //----------------------------------------------
  // In place arithmetic

  inline SlVector3 &operator+=(double d);
  inline SlVector3 &operator-=(double d);
  inline SlVector3 &operator*=(double d);
  inline SlVector3 &operator/=(double d);

  inline SlVector3 &operator+=(const SlVector3 &da);
  inline SlVector3 &operator-=(const SlVector3 &da);
  inline SlVector3 &operator*=(const SlVector3 &da);
  inline SlVector3 &operator/=(const SlVector3 &da);
  // Componentwise operations

  inline SlVector3 &maxSet(const SlVector3 &da);
  inline SlVector3 &minSet(const SlVector3 &da);
  // Sets data[i] = max(data[i],da[i]) or min


  //----------------------------------------------
  // Static methods

  inline static unsigned int cycleAxis(unsigned int axis, int direction);

  //----------------------------------------------
  // Define Components

  enum Index { X = 0 , Y = 1 , Z = 2 ,
	       U = 0 , V = 1 , W = 2 ,
	       R = 0 , G = 1 , B = 2 };

public:
  //----------------------------------------------
  // Public data members

  double data[3];

  //----------------------------------------------
};


//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Operators for class SlVector3

inline SlVector3 operator-(const SlVector3 &a);

inline SlVector3 operator+(const SlVector3 &a,const SlVector3 &b);
inline SlVector3 operator-(const SlVector3 &a,const SlVector3 &b);
inline SlVector3 operator*(const SlVector3 &a,const SlVector3 &b);
inline SlVector3 operator/(const SlVector3 &a,const SlVector3 &b);

inline SlVector3 operator+(const SlVector3 &a, double b);
inline SlVector3 operator-(const SlVector3 &a, double b);
inline SlVector3 operator*(const SlVector3 &a, double b);
inline SlVector3 operator/(const SlVector3 &a, double b);

inline SlVector3 operator+(double a, const SlVector3 &b);
inline SlVector3 operator-(double a, const SlVector3 &b);
inline SlVector3 operator*(double a, const SlVector3 &b);
inline SlVector3 operator/(double a, const SlVector3 &b);

std::istream &operator>>(std::istream &strm,      SlVector3 &v);
std::ostream &operator<<(std::ostream &strm,const SlVector3 &v);

//-------------------------------------------------------------------
// Norm type functions for SlVector3

inline double   l1Norm(const SlVector3 &a);
inline double   l2Norm(const SlVector3 &a);
inline double lInfNorm(const SlVector3 &a);
// Computes the l1, l2 or lInfinity norm of a

inline double    mag(const SlVector3 &a);
inline double sqrMag(const SlVector3 &a);
// mag is the l2Norm or magnitude of the vector
// sqrMag is mag^2, which is faster to compute

inline double normalize(SlVector3 &a);
// SETS a = a/mag(a)


//-------------------------------------------------------------------
// Other functions for SlVector3

inline unsigned int dominantAxis(const SlVector3 &v);
inline unsigned int subinantAxis(const SlVector3 &v);
inline unsigned int midinantAxis(const SlVector3 &v);
// Returns the index of the component with the largest,
// smallest or middle value.  Note: subinantAxis and
// midinantAxis are nore really words, I made them up.
// If multiple comonents have the same value, then the
// results are not unique.

inline double      dot(const SlVector3 &a,const SlVector3 &b);
inline SlVector3 cross(const SlVector3 &a,const SlVector3 &b);
// Compute the dot and cros product of a and b.

inline double box(const SlVector3 &a,const SlVector3 &b,const SlVector3 &c);
// Compute the box (aka tripple) product of a, b, and d.

inline SlVector3 abs(const SlVector3 &a);
// returns a vector with r[i] = abs(a[i])

inline double sum(const SlVector3 &a);
// return a[0]+a[1]+a[2]

inline double max(const SlVector3 &a);
inline double min(const SlVector3 &a);
// Returns the max or min component of a.

inline SlVector3 max(const SlVector3 &a,const SlVector3 &b);
inline SlVector3 min(const SlVector3 &a,const SlVector3 &b);
// Computes a NEW vector by taking the max component in the
// x,y,z direction from a and b.  ie: r[i] = max(a[i],b[i])
// Note: signed values are used.


//-------------------------------------------------------------------
//-------------------------------------------------------------------


class SlVector2 {
  // Class to store a two dimensional vector.
  // Expected uses include point in R2, etc.
public:

  //----------------------------------------------
  // Constructors

  inline SlVector2();
  inline SlVector2(double d);
  inline SlVector2(double d0,double d1);

  inline SlVector2(const SlVector2 &da);
  inline SlVector2(const double    *da);
  // da should point to a double[2] that will be copied


  //----------------------------------------------
  // Index operators

  inline double &operator[](unsigned int i)      ;
  inline double  operator[](unsigned int i) const;

  inline double &operator()(unsigned int i)      ;
  inline double  operator()(unsigned int i) const;


  //----------------------------------------------
  // Assignment and set

  inline SlVector2 &set(double d);
  inline SlVector2 &set(double d0, double d1);

  inline SlVector2 &set(const SlVector2 &da);
  inline SlVector2 &set(const double    *da);
  // da should point to a double[3] that will be copied

  inline SlVector2 &operator=(double d);
  inline SlVector2 &operator=(const SlVector2 &da);
  inline SlVector2 &operator=(const double    *da);

  //----------------------------------------------
  // Comparison operators

  inline int operator==(const SlVector2 &da) const;
  inline int operator!=(const SlVector2 &da) const;

  inline int operator==(double d) const;
  inline int operator!=(double d) const;


  //----------------------------------------------
  // In place arithmetic

  inline SlVector2 &operator+=(double d);
  inline SlVector2 &operator-=(double d);
  inline SlVector2 &operator*=(double d);
  inline SlVector2 &operator/=(double d);

  inline SlVector2 &operator+=(const SlVector2 &da);
  inline SlVector2 &operator-=(const SlVector2 &da);
  inline SlVector2 &operator*=(const SlVector2 &da);
  inline SlVector2 &operator/=(const SlVector2 &da);
  // Componentwise operations

  inline SlVector2 &maxSet(const SlVector2 &da);
  inline SlVector2 &minSet(const SlVector2 &da);
  // Sets data[i] = max(data[i],da[i]) or min


  //----------------------------------------------
  // Static methods

  inline static unsigned int cycleAxis(unsigned int axis, int direction);

  //----------------------------------------------
  // Define Components

  enum Index { X = 0 , Y = 1 ,
	       U = 0 , V = 1 };

public:
  //----------------------------------------------
  // Public data members

  double data[2];

  //----------------------------------------------
};


//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Operators for class SlVector2

inline SlVector2 operator-(const SlVector2 &a);

inline SlVector2 operator+(const SlVector2 &a, const SlVector2 &b);
inline SlVector2 operator-(const SlVector2 &a, const SlVector2 &b);
inline SlVector2 operator*(const SlVector2 &a, const SlVector2 &b);
inline SlVector2 operator/(const SlVector2 &a, const SlVector2 &b);

inline SlVector2 operator+(const SlVector2 &a, double b);
inline SlVector2 operator-(const SlVector2 &a, double b);
inline SlVector2 operator*(const SlVector2 &a, double b);
inline SlVector2 operator/(const SlVector2 &a, double b);

inline SlVector2 operator+(double a, const SlVector2 &b);
inline SlVector2 operator-(double a, const SlVector2 &b);
inline SlVector2 operator*(double a, const SlVector2 &b);
inline SlVector2 operator/(double a, const SlVector2 &b);

std::istream &operator>>(std::istream &strm,      SlVector2 &v);
std::ostream &operator<<(std::ostream &strm,const SlVector2 &v);

//-------------------------------------------------------------------
// Norm type functions for SlVector2

inline double   l1Norm(const SlVector2 &a);
inline double   l2Norm(const SlVector2 &a);
inline double lInfNorm(const SlVector2 &a);
// Computes the l1, l2 or lInfinity norm of a

inline double    mag(const SlVector2 &a);
inline double sqrMag(const SlVector2 &a);
// mag is the l2Norm or magnitude of the vector
// sqrMag is mag^2, which is faster to compute

inline double normalize(SlVector2 &a);
// SETS a = a/mag(a)


//-------------------------------------------------------------------
// Other functions for SlVector2

inline unsigned int dominantAxis(const SlVector2 &v);
inline unsigned int subinantAxis(const SlVector2 &v);
// Returns the index of the component with the largest,
// smallest value.  Note: subinantAxis and is not really
// a word, I made it up.
// If multiple comonents have the same value, then the
// results are not unique.

inline double   dot(const SlVector2 &a,const SlVector2 &b);
inline double cross(const SlVector2 &a,const SlVector2 &b);
// Compute the dot and cros product of a and b.

inline SlVector2 abs(const SlVector2 &a);
// returns a vector with r[i] = abs(a[i])

inline double sum(const SlVector2 &a);
// return a[0]+a[1]

inline double max(const SlVector2 &a);
inline double min(const SlVector2 &a);
// Returns the max or min component of a.

inline SlVector2 max(const SlVector2 &a,const SlVector2 &b);
inline SlVector2 min(const SlVector2 &a,const SlVector2 &b);
// Computes a NEW vector by taking the max component in the
// x,y,z direction from a and b.  ie: r[i] = max(a[i],b[i])
// Note: signed values are used.




//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Inline implementation below:::

#ifndef DB_CHECK
#ifdef DEBUG
#define DB_CHECK( C ) { if ( ! (C) ) { abort(); } }
#else
#define DB_CHECK( C ) { }
#endif
#endif

//-------------------------------------------------------------------
//-------------------------------------------------------------------

// Inline implementation of SlVector3

inline double &SlVector3::operator[](unsigned int i) {
  DB_CHECK(i<3);
  return data[i];
}
inline double &SlVector3::operator()(unsigned int i) {
  DB_CHECK(i<3);
  return data[i];
}

inline double  SlVector3::operator[](unsigned int i) const {
  DB_CHECK(i<3);
  return data[i];
}

inline double  SlVector3::operator()(unsigned int i) const {
  DB_CHECK(i<3);
  return data[i];
}

inline double SlVector3::x() const {
  return data[0];
}
inline double SlVector3::y() const {
  return data[1];
}
inline double SlVector3::z() const {
  return data[2];
}
//-------------------------------------------------------------------

inline SlVector3::SlVector3() {
  data[0] = data[1] = data[2] = 0.0;
}
inline SlVector3::SlVector3(double d) {
  data[0] = data[1] = data[2] = d;
}

inline SlVector3::SlVector3(double d0,double d1,double d2) {
  data[0] = d0;
  data[1] = d1;
  data[2] = d2;
}

inline SlVector3::SlVector3(const SlVector3 &da) {
  data[0] = da[0];
  data[1] = da[1];
  data[2] = da[2];
}

inline SlVector3::SlVector3(const double *da) {
  data[0] = da[0];
  data[1] = da[1];
  data[2] = da[2];
}

//-------------------------------------------------------------------

inline SlVector3 &SlVector3::set(double d) {
  data[0] = d;
  data[1] = d;
  data[2] = d;
  return (*this);
}

inline SlVector3 &SlVector3::set(double d0, double d1, double d2) {
  data[0] = d0;
  data[1] = d1;
  data[2] = d2;
  return (*this);
}

inline SlVector3 &SlVector3::set(const SlVector3 &da) {
  data[0] = da[0];
  data[1] = da[1];
  data[2] = da[2];
  return (*this);
}

inline SlVector3 &SlVector3::set(const double *da) {
  data[0] = da[0];
  data[1] = da[1];
  data[2] = da[2];
  return (*this);
}

//-------------------------------------------------------------------

inline SlVector3 &SlVector3::operator=(double d) {
  return set(d);
}

inline SlVector3 &SlVector3::operator=(const SlVector3 &da) {
  return set(da);
}

inline SlVector3 &SlVector3::operator=(const double *da) {
  return set(da);
}

//-------------------------------------------------------------------

inline int SlVector3::operator==(const SlVector3 &da) const {
  return ((data[0] == da[0]) &&
	  (data[1] == da[1]) &&
	  (data[2] == da[2]));
}

inline int SlVector3::operator!=(const SlVector3 &da) const {
  return ((data[0] != da[0]) ||
	  (data[1] != da[1]) ||
	  (data[2] != da[2]));
}

inline int SlVector3::operator==(double d) const {
  return ((data[0] == d) &&
	  (data[1] == d) &&
	  (data[2] == d));
}

inline int SlVector3::operator!=(double d) const {
  return ((data[0] != d) ||
	  (data[1] != d) ||
	  (data[2] != d));
}

//-------------------------------------------------------------------

inline SlVector3 &SlVector3::operator+=(double d) {
  data[0] += d;
  data[1] += d;
  data[2] += d;
  return (*this);
}

inline SlVector3 &SlVector3::operator-=(double d) {
  data[0] -= d;
  data[1] -= d;
  data[2] -= d;
  return (*this);
}

inline SlVector3 &SlVector3::operator*=(double d) {
  data[0] *= d;
  data[1] *= d;
  data[2] *= d;
  return (*this);
}

inline SlVector3 &SlVector3::operator/=(double d) {
  data[0] /= d;
  data[1] /= d;
  data[2] /= d;
  return (*this);
}

//-------------------------------------------------------------------

inline SlVector3 &SlVector3::operator+=(const SlVector3 &da) {
  data[0] += da[0];
  data[1] += da[1];
  data[2] += da[2];
  return (*this);
}

inline SlVector3 &SlVector3::operator-=(const SlVector3 &da) {
  data[0] -= da[0];
  data[1] -= da[1];
  data[2] -= da[2];
  return (*this);
}

inline SlVector3 &SlVector3::operator*=(const SlVector3 &da) {
  data[0] *= da[0];
  data[1] *= da[1];
  data[2] *= da[2];
  return (*this);
}

inline SlVector3 &SlVector3::operator/=(const SlVector3 &da) {
  data[0] /= da[0];
  data[1] /= da[1];
  data[2] /= da[2];
  return (*this);
}

//-------------------------------------------------------------------

inline SlVector3 &SlVector3::maxSet(const SlVector3 &da) {
  if (da[0] > data[0]) data[0] = da[0];
  if (da[1] > data[1]) data[1] = da[1];
  if (da[2] > data[2]) data[2] = da[2];
  return (*this);
}

inline SlVector3 &SlVector3::minSet(const SlVector3 &da) {
  if (da[0] < data[0]) data[0] = da[0];
  if (da[1] < data[1]) data[1] = da[1];
  if (da[2] < data[2]) data[2] = da[2];
  return (*this);
}

//-------------------------------------------------------------------

inline unsigned int SlVector3::cycleAxis(unsigned int axis, int direction) {
  switch (axis+direction) {
  case 0: case 3: case 6: return 0;
  case 1: case 4: case 7: return 1;
  case 2: case 5: case 8: return 2;
  default: return (axis+direction)%3;
  }
}

//-------------------------------------------------------------------

inline SlVector3 operator-(const SlVector3 &a) {
  return SlVector3(-a[0],-a[1],-a[2]);
}

//-------------------------------------------------------------------

inline SlVector3 operator+(const SlVector3 &a,const SlVector3 &b) {
  return SlVector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

inline SlVector3 operator-(const SlVector3 &a,const SlVector3 &b) {
  return SlVector3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

inline SlVector3 operator*(const SlVector3 &a,const SlVector3 &b){
  return SlVector3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

inline SlVector3 operator/(const SlVector3 &a,const SlVector3 &b){
  return SlVector3(a[0] / b[0], a[1] / b[1], a[2] / b[2]);
}

//-------------------------------------------------------------------

inline SlVector3 operator+(const SlVector3 &a,double b){
  return SlVector3(a[0] + b, a[1] + b, a[2] + b);
}

inline SlVector3 operator-(const SlVector3 &a,double b){
  return SlVector3(a[0] - b, a[1] - b, a[2] - b);
}

inline SlVector3 operator*(const SlVector3 &a,double b){
  return SlVector3(a[0] * b, a[1] * b, a[2] * b);
}

inline SlVector3 operator/(const SlVector3 &a,double b){
  return SlVector3(a[0] / b, a[1] / b, a[2] / b);
}

//-------------------------------------------------------------------

inline SlVector3 operator+(double a,const SlVector3 &b){
  return SlVector3(a + b[0], a + b[1], a + b[2]);
}

inline SlVector3 operator-(double a,const SlVector3 &b){
  return SlVector3(a - b[0], a - b[1], a - b[2]);
}

inline SlVector3 operator*(double a,const SlVector3 &b){
  return SlVector3(a * b[0], a * b[1], a * b[2]);
}

inline SlVector3 operator/(double a,const SlVector3 &b){
  return SlVector3(a / b[0], a / b[1], a / b[2]);
}

//-------------------------------------------------------------------

inline double l1Norm(const SlVector3 &a) {
  return (((a[0]>0)?a[0]:-a[0])+
	  ((a[1]>0)?a[1]:-a[1])+
	  ((a[2]>0)?a[2]:-a[2]));
}

inline double l2Norm(const SlVector3 &a) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

inline double lInfNorm(const SlVector3 &a) {
  return max(abs(a));
}

inline double mag(const SlVector3 &a) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

inline double sqrMag(const SlVector3 &a) {
  return (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

inline double normalize(SlVector3 &a) {
	double m = mag(a);
  if (m != 0) a /= m;
  return m;
}

//-------------------------------------------------------------------

inline unsigned int dominantAxis(const SlVector3 &v) {
  double x,y,z;
  if (v[0]>0) x = v[0]; else x = -v[0];
  if (v[1]>0) y = v[1]; else y = -v[1];
  if (v[2]>0) z = v[2]; else z = -v[2];
  return ( x > y ) ? (( x > z ) ? 0 : 2) : (( y > z ) ? 1 : 2 );
}

inline unsigned int subinantAxis(const SlVector3 &v) {
  double x,y,z;
  if (v[0]>0) x = v[0]; else x = -v[0];
  if (v[1]>0) y = v[1]; else y = -v[1];
  if (v[2]>0) z = v[2]; else z = -v[2];
  return ( x < y ) ? (( x < z ) ? 0 : 2) : (( y < z ) ? 1 : 2 );
}

inline unsigned int midinantAxis(const SlVector3 &v) {
  double x,y,z;
  if (v[0]>0) x = v[0]; else x = -v[0];
  if (v[1]>0) y = v[1]; else y = -v[1];
  if (v[2]>0) z = v[2]; else z = -v[2];
  unsigned int d = ( x > y ) ? (( x > z ) ? 0 : 2) : (( y > z ) ? 1 : 2 );
  unsigned int s = ( x < y ) ? (( x < z ) ? 0 : 2) : (( y < z ) ? 1 : 2 );
  unsigned int m = 0;
  if (d==0) {
    if (s!= 1) m = 1; else m = 2;
  }else if (d==1) {
    if (s!= 0) m = 0; else m = 2;
  }else if (d==2) {
    if (s!= 0) m = 0; else m = 1;
  }
  return m;
}

//-------------------------------------------------------------------

inline double dot(const SlVector3 &a,const SlVector3 &b) {
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

inline SlVector3 cross(const SlVector3 &a,const SlVector3 &b) {
  return SlVector3(a[1] * b[2] - b[1] * a[2],
		   a[2] * b[0] - b[2] * a[0],
		   a[0] * b[1] - b[0] * a[1]);
}

//-------------------------------------------------------------------

inline double box(const SlVector3 &a,const SlVector3 &b,const SlVector3 &c) {
  return dot(cross(a,b),c);
}

//-------------------------------------------------------------------

inline SlVector3 abs(const SlVector3 &a) {
  return  SlVector3(((a[0]>0)?a[0]:-a[0]),
		    ((a[1]>0)?a[1]:-a[1]),
		    ((a[2]>0)?a[2]:-a[2]));
}

inline double sum(const SlVector3 &a) {
  return a[0]+a[1]+a[2];
}

//-------------------------------------------------------------------

inline double max(const SlVector3 &a) {
  return ((a[0]>a[1])?((a[0]>a[2])?a[0]:a[2]):(a[1]>a[2])?a[1]:a[2]);
}

inline double min(const SlVector3 &a) {
  return ((a[0]<a[1])?((a[0]<a[2])?a[0]:a[2]):(a[1]<a[2])?a[1]:a[2]);
}

//-------------------------------------------------------------------

inline SlVector3 max(const SlVector3 &a,const SlVector3 &b) {
  return SlVector3((a[0]>b[0])?a[0]:b[0],
		   (a[1]>b[1])?a[1]:b[1],
		   (a[2]>b[2])?a[2]:b[2]);
}

inline SlVector3 min(const SlVector3 &a,const SlVector3 &b) {
  return SlVector3((a[0]<b[0])?a[0]:b[0],
		   (a[1]<b[1])?a[1]:b[1],
		   (a[2]<b[2])?a[2]:b[2]);
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------

// Inline implementation of SlVector2

inline double &SlVector2::operator[](unsigned int i) {
  DB_CHECK(i<2);
  return data[i];
}
inline double &SlVector2::operator()(unsigned int i) {
  DB_CHECK(i<2);
  return data[i];
}

inline double  SlVector2::operator[](unsigned int i) const {
  DB_CHECK(i<2);
  return data[i];
}

inline double  SlVector2::operator()(unsigned int i) const {
  DB_CHECK(i<2);
  return data[i];
}

//-------------------------------------------------------------------

inline SlVector2::SlVector2() {
  data[0] = data[1] = 0.0;
}
inline SlVector2::SlVector2(double d) {
  data[0] = data[1] = d;
}

inline SlVector2::SlVector2(double d0,double d1) {
  data[0] = d0;
  data[1] = d1;
}

inline SlVector2::SlVector2(const SlVector2 &da) {
  data[0] = da[0];
  data[1] = da[1];
}

inline SlVector2::SlVector2(const double *da) {
  data[0] = da[0];
  data[1] = da[1];
}

//-------------------------------------------------------------------

inline SlVector2 &SlVector2::set(double d) {
  data[0] = d;
  data[1] = d;
  return (*this);
}

inline SlVector2 &SlVector2::set(double d0, double d1) {
  data[0] = d0;
  data[1] = d1;
  return (*this);
}

inline SlVector2 &SlVector2::set(const SlVector2 &da) {
  data[0] = da[0];
  data[1] = da[1];
  return (*this);
}

inline SlVector2 &SlVector2::set(const double *da) {
  data[0] = da[0];
  data[1] = da[1];
  return (*this);
}

//-------------------------------------------------------------------

inline SlVector2 &SlVector2::operator=(double d) {
  return set(d);
}

inline SlVector2 &SlVector2::operator=(const SlVector2 &da) {
  return set(da);
}

inline SlVector2 &SlVector2::operator=(const double *da) {
  return set(da);
}

//-------------------------------------------------------------------

inline int SlVector2::operator==(const SlVector2 &da) const {
  return ((data[0] == da[0]) &&
	  (data[1] == da[1]));
}

inline int SlVector2::operator!=(const SlVector2 &da) const {
  return ((data[0] != da[0]) ||
	  (data[1] != da[1]));
}

inline int SlVector2::operator==(double d) const {
  return ((data[0] == d) &&
	  (data[1] == d));
}

inline int SlVector2::operator!=(double d) const {
  return ((data[0] != d) ||
	  (data[1] != d));
}

//-------------------------------------------------------------------

inline SlVector2 &SlVector2::operator+=(double d) {
  data[0] += d;
  data[1] += d;
  return (*this);
}

inline SlVector2 &SlVector2::operator-=(double d) {
  data[0] -= d;
  data[1] -= d;
  return (*this);
}

inline SlVector2 &SlVector2::operator*=(double d) {
  data[0] *= d;
  data[1] *= d;
  return (*this);
}

inline SlVector2 &SlVector2::operator/=(double d) {
  data[0] /= d;
  data[1] /= d;
  return (*this);
}

//-------------------------------------------------------------------

inline SlVector2 &SlVector2::operator+=(const SlVector2 &da) {
  data[0] += da[0];
  data[1] += da[1];
  return (*this);
}

inline SlVector2 &SlVector2::operator-=(const SlVector2 &da) {
  data[0] -= da[0];
  data[1] -= da[1];
  return (*this);
}

inline SlVector2 &SlVector2::operator*=(const SlVector2 &da) {
  data[0] *= da[0];
  data[1] *= da[1];
  return (*this);
}

inline SlVector2 &SlVector2::operator/=(const SlVector2 &da) {
  data[0] /= da[0];
  data[1] /= da[1];
  return (*this);
}

//-------------------------------------------------------------------

inline SlVector2 &SlVector2::maxSet(const SlVector2 &da) {
  if (da[0] > data[0]) data[0] = da[0];
  if (da[1] > data[1]) data[1] = da[1];
  return (*this);
}

inline SlVector2 &SlVector2::minSet(const SlVector2 &da) {
  if (da[0] < data[0]) data[0] = da[0];
  if (da[1] < data[1]) data[1] = da[1];
  return (*this);
}

//-------------------------------------------------------------------

inline unsigned int SlVector2::cycleAxis(unsigned int axis, int direction) {
  switch (axis+direction) {
  case 0: case 2: case 4: return 0;
  case 1: case 3: case 5: return 1;
  default: return (axis+direction)%2;
  }
}

//-------------------------------------------------------------------

inline SlVector2 operator-(const SlVector2 &a) {
  return SlVector2(-a[0],-a[1]);
}

//-------------------------------------------------------------------

inline SlVector2 operator+(const SlVector2 &a,const SlVector2 &b) {
  return SlVector2(a[0] + b[0], a[1] + b[1]);
}

inline SlVector2 operator-(const SlVector2 &a,const SlVector2 &b) {
  return SlVector2(a[0] - b[0], a[1] - b[1]);
}

inline SlVector2 operator*(const SlVector2 &a,const SlVector2 &b){
  return SlVector2(a[0] * b[0], a[1] * b[1]);
}

inline SlVector2 operator/(const SlVector2 &a,const SlVector2 &b){
  return SlVector2(a[0] / b[0], a[1] / b[1]);
}

//-------------------------------------------------------------------

inline SlVector2 operator+(const SlVector2 &a,double b){
  return SlVector2(a[0] + b, a[1] + b);
}

inline SlVector2 operator-(const SlVector2 &a,double b){
  return SlVector2(a[0] - b, a[1] - b);
}

inline SlVector2 operator*(const SlVector2 &a,double b){
  return SlVector2(a[0] * b, a[1] * b);
}

inline SlVector2 operator/(const SlVector2 &a,double b){
  return SlVector2(a[0] / b, a[1] / b);
}

//-------------------------------------------------------------------

inline SlVector2 operator+(double a,const SlVector2 &b){
  return SlVector2(a + b[0], a + b[1]);
}

inline SlVector2 operator-(double a,const SlVector2 &b){
  return SlVector2(a - b[0], a - b[1]);
}

inline SlVector2 operator*(double a,const SlVector2 &b){
  return SlVector2(a * b[0], a * b[1]);
}

inline SlVector2 operator/(double a,const SlVector2 &b){
  return SlVector2(a / b[0], a / b[1]);
}

//-------------------------------------------------------------------

inline double l1Norm(const SlVector2 &a) {
  return (((a[0]>0)?a[0]:-a[0])+
	  ((a[1]>0)?a[1]:-a[1]));
}

inline double l2Norm(const SlVector2 &a) {
  return sqrt(a[0] * a[0] + a[1] * a[1]);
}

inline double lInfNorm(const SlVector2 &a) {
  return max(abs(a));
}

inline double mag(const SlVector2 &a) {
  return sqrt(a[0] * a[0] + a[1] * a[1]);
}

inline double sqrMag(const SlVector2 &a) {
  return (a[0] * a[0] + a[1] * a[1]);
}

inline double normalize(SlVector2 &a) {
  double m = mag(a);
  if (m != 0) a /= m;
  return m;
}

//-------------------------------------------------------------------

inline unsigned int dominantAxis(const SlVector2 &v) {
  double x,y;
  if (v[0]>0) x = v[0]; else x = -v[0];
  if (v[1]>0) y = v[1]; else y = -v[1];
  return ( x > y ) ? 0 : 1;
}

inline unsigned int subinantAxis(const SlVector2 &v) {
  double x,y;
  if (v[0]>0) x = v[0]; else x = -v[0];
  if (v[1]>0) y = v[1]; else y = -v[1];
  return ( x < y ) ? 0 : 1;
}

//-------------------------------------------------------------------

inline double dot(const SlVector2 &a,const SlVector2 &b) {
  return (a[0] * b[0] + a[1] * b[1]);
}

inline double cross(const SlVector2 &a,const SlVector2 &b) {
  return (a[0] * b[1] - b[0] * a[1]);
}

//-------------------------------------------------------------------

inline SlVector2 abs(const SlVector2 &a) {
  return  SlVector2(((a[0]>0)?a[0]:-a[0]),
		    ((a[1]>0)?a[1]:-a[1]));
}

inline double sum(const SlVector2 &a) {
  return a[0]+a[1];
}

//-------------------------------------------------------------------

inline double max(const SlVector2 &a) {
  return ((a[0]>a[1])? a[0] : a[1]);
}

inline double min(const SlVector2 &a) {
  return ((a[0]<a[1])? a[0] : a[1]);
}

//-------------------------------------------------------------------

inline SlVector2 max(const SlVector2 &a,const SlVector2 &b) {
  return SlVector2((a[0]>b[0])?a[0]:b[0],
		   (a[1]>b[1])?a[1]:b[1]);
}

inline SlVector2 min(const SlVector2 &a,const SlVector2 &b) {
  return SlVector2((a[0]<b[0])?a[0]:b[0],
		   (a[1]<b[1])?a[1]:b[1]);
}



//-------------------------------------------------------------------
//-------------------------------------------------------------------

#endif

