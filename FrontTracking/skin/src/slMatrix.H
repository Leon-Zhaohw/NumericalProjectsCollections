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

#ifndef SlMATRIX_is_defined
#define SlMATRIX_is_defined


//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Classes from this header:

class SlMatrix3x3;
class SlMatrix2x2;

// External classes

class SlVector3;
class SlVector2;

#include <iostream>

//-------------------------------------------------------------------
//-------------------------------------------------------------------

class SlMatrix3x3 {
  // Stores a 3x3 row ordered matrix. Vectors are assumed to be row
  // vectors and matricies are indexed (row,col).  NOTE: The *
  // operator is matrix multiplication, not a component-wise
  // operation.


public:

  //----------------------------------------------
  // Constructors

  enum ordering { row , col };

  inline SlMatrix3x3();
  inline SlMatrix3x3(double x);
  inline SlMatrix3x3(const SlMatrix3x3 &that);

  inline SlMatrix3x3(double xx, double xy, double xz,
		     double yx, double yy, double yz,
		     double zx, double zy, double zz);

  inline SlMatrix3x3(const double data[3][3], ordering o = row);
  inline SlMatrix3x3(const SlVector3 data[3], ordering o = row);
  inline SlMatrix3x3(const SlVector3 a      ,
		     const SlVector3 b      ,
		     const SlVector3 c      , ordering o = row);

  //----------------------------------------------
  // Index operators

  inline double &operator()(unsigned int r,unsigned int c)      ;
  inline double  operator()(unsigned int r,unsigned int c) const;


  //----------------------------------------------
  // Assignment and set

  inline SlMatrix3x3 &set      (const double d);
  inline SlMatrix3x3 &operator=(const double d);

  inline SlMatrix3x3 &set      (const SlMatrix3x3 &that);  
  inline SlMatrix3x3 &operator=(const SlMatrix3x3 &that); 

  inline SlMatrix3x3 &set      (const double xx, const double xy, const double xz,
																const double yx, const double yy, const double yz,
																const double zx, const double zy, const double zz);


  //----------------------------------------------
  // Comparison operators

  inline int operator==(double d) const;
  inline int operator!=(double d) const;
  
  inline int operator==(const SlMatrix3x3 &that)const; 
  inline int operator!=(const SlMatrix3x3 &that)const; 

  
  //----------------------------------------------
  // In place arithmetic

  inline SlMatrix3x3 &operator+=(double d);
  inline SlMatrix3x3 &operator-=(double d);
  inline SlMatrix3x3 &operator*=(double d);
  inline SlMatrix3x3 &operator/=(double d);


  inline SlMatrix3x3 &operator+=(const SlMatrix3x3 &that);
  inline SlMatrix3x3 &operator-=(const SlMatrix3x3 &that);
  // These are all component-wise operations.

  SlMatrix3x3 &inplaceMultPre (const SlMatrix3x3 &that);
  SlMatrix3x3 &inplaceMultPost(const SlMatrix3x3 &that);
  // Pre : this = that x this  
  // Post: this = this x that

  inline SlMatrix3x3 &componentAdd (const SlMatrix3x3 &that);
  inline SlMatrix3x3 &componentSub (const SlMatrix3x3 &that);
  inline SlMatrix3x3 &componentMult(const SlMatrix3x3 &that);
  inline SlMatrix3x3 &componentDiv (const SlMatrix3x3 &that);

  inline SlMatrix3x3 &setIdentity     ();
  inline SlMatrix3x3 &inplaceTranspose();

  //----------------------------------------------
  // Static methods

  inline static unsigned int cycleAxis(unsigned int axis, int direction);


  //----------------------------------------------
  // Define Components

  enum Index { X = 0 , Y = 1 , Z = 2 ,
	       U = 0 , V = 1 , W = 2 };


public:
  //----------------------------------------------
  // Accessable data members

  double data[3][3];
};
  


//-------------------------------------------------------------------
//-------------------------------------------------------------------

class SlMatrix2x2 {
  // Stores a 2x2 row ordered matrix. Vectors are assumed to be row
  // vectors and matricies are indexed (row,col).  NOTE: The *
  // operator is matrix multiplication, not a component-wise
  // operation.

public:

  //----------------------------------------------
  // Constructors

  enum ordering { row , col };

  inline SlMatrix2x2();
  inline SlMatrix2x2(double x);
  inline SlMatrix2x2(const SlMatrix2x2 &that);

  inline SlMatrix2x2(double xx, double xy, 
		     double yx, double yy);

  inline SlMatrix2x2(const double data[2][2], ordering o = row);
  inline SlMatrix2x2(const SlVector2 data[2], ordering o = row);
  inline SlMatrix2x2(const SlVector2 a      ,
		     const SlVector2 b      , ordering o = row);


  //----------------------------------------------
  // Index operators

  inline double &operator()(unsigned int r,unsigned int c)      ;
  inline double  operator()(unsigned int r,unsigned int c) const;


  //----------------------------------------------
  // Assignment and set

  inline SlMatrix2x2 &set      (const double d);
  inline SlMatrix2x2 &operator=(const double d);

  inline SlMatrix2x2 &set      (const SlMatrix2x2 &that);  
  inline SlMatrix2x2 &operator=(const SlMatrix2x2 &that); 


  //----------------------------------------------
  // Comparison operators

  inline int operator!=(const SlMatrix2x2 &that)const; 
  inline int operator==(const SlMatrix2x2 &that)const; 

  inline int operator==(double d) const;
  inline int operator!=(double d) const;
  
  
  //----------------------------------------------
  // In place arithmetic

  inline SlMatrix2x2 &operator+=(double d);
  inline SlMatrix2x2 &operator-=(double d);
  inline SlMatrix2x2 &operator*=(double d);
  inline SlMatrix2x2 &operator/=(double d);


  inline SlMatrix2x2 &operator+=(const SlMatrix2x2 &that);
  inline SlMatrix2x2 &operator-=(const SlMatrix2x2 &that);
  // These are all component-wise operations.

  SlMatrix2x2 &inplaceMultPre (const SlMatrix2x2 &that);
  SlMatrix2x2 &inplaceMultPost(const SlMatrix2x2 &that);
  // Pre : this = that x this  
  // Post: this = this x that

  inline SlMatrix2x2 &componentAdd (const SlMatrix2x2 &that);
  inline SlMatrix2x2 &componentSub (const SlMatrix2x2 &that);
  inline SlMatrix2x2 &componentDiv (const SlMatrix2x2 &that);
  inline SlMatrix2x2 &componentMult(const SlMatrix2x2 &that);

  inline SlMatrix2x2 &setIdentity     ();
  inline SlMatrix2x2 &inplaceTranspose();

  //----------------------------------------------
  // Static methods

  inline static unsigned int cycleAxis(unsigned int axis, int direction);


  //----------------------------------------------
  // Define Components

  enum Index { X = 0 , Y = 1 ,
	       U = 0 , V = 1 };


public:
  //----------------------------------------------
  // Accessable data members

  double data[2][2];

};
  

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Operators for class SlMatrix3x3

inline SlMatrix3x3 operator+(const SlMatrix3x3 &a, double b);
inline SlMatrix3x3 operator-(const SlMatrix3x3 &a, double b);
inline SlMatrix3x3 operator*(const SlMatrix3x3 &a, double b);
inline SlMatrix3x3 operator/(const SlMatrix3x3 &a, double b);

inline SlMatrix3x3 operator+(double a, const SlMatrix3x3 &b);
inline SlMatrix3x3 operator-(double a, const SlMatrix3x3 &b);
inline SlMatrix3x3 operator*(double a, const SlMatrix3x3 &b);
inline SlMatrix3x3 operator/(double a, const SlMatrix3x3 &b);

inline SlMatrix3x3 operator+(const SlMatrix3x3 &a, const SlMatrix3x3 &b);
inline SlMatrix3x3 operator-(const SlMatrix3x3 &a, const SlMatrix3x3 &b);
// These are component-wise operations. 

inline SlMatrix3x3 operator*(const SlMatrix3x3 &a, const SlMatrix3x3 &b); 
inline SlVector3   operator*(const SlMatrix3x3 &a, const SlVector3   &b);
inline SlVector3   operator*(const SlVector3   &a, const SlMatrix3x3 &b);
// Matrix multiplication.  Vector returning routines treat the vector
// argument like a row or column vector as appropriate.


//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Global funtions for class SlMatrix3x3

inline double determinant(const SlMatrix3x3 &a);

inline double trace(const SlMatrix3x3 &a);


inline SlMatrix3x3 transpose(const SlMatrix3x3 &a);
/*  */ SlMatrix3x3   inverse(const SlMatrix3x3 &a);

inline SlMatrix3x3 diagonal(const SlVector3 &a);


void SlSVDecomp(const SlMatrix3x3 &a,
		/* */ SlMatrix3x3 &u,
		/* */ SlVector3   &w,
		/* */ SlMatrix3x3 &v);
// Computes the singular value decomposition of a matrix, [u,w,v].
// The input matrix, a, is related to the output by a =
// u*diagonal(w)*transpose(v).  Both u and v are orthonormal.

void SlSymetricEigenDecomp(const SlMatrix3x3 &a,
			   /* */ SlVector3   &vals,
			   /* */ SlMatrix3x3 &vecs);
// Given a SYMETRIC matrix, a, it returns the Eigen Values in vals,
// and the Eigen Vectors as the COLUMNS of vecs.  The vectors are
// normalized.  The values/vectors are not sorted in any order.

inline SlMatrix3x3 SlMatrixProduct(const SlVector3 a,const SlVector3 b);
// Returns the matrix given by a * aT  [ m_ij = a_i * b_j ]

inline SlMatrix3x3 SlEigenMatrix(const SlVector3 a);
// Returns SlMatrixProduct(a,a)/|a| [ m_ij = a_i * a_j / (a_k * a_k)].
// Note that the (only) eigen vector of the matrix is a, and the
// matrix is symetric.  (Involves a call to sqrt.)


inline double i2(const SlMatrix3x3 &m);

//Frobenius norm
inline double norm(const SlMatrix3x3 &x); 

std::istream &operator>>(std::istream &strm,      SlMatrix3x3 &m);
std::ostream &operator<<(std::ostream &strm,const SlMatrix3x3 &m);



//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Operators for class SlMatrix2x2

inline SlMatrix2x2 operator+(const SlMatrix2x2 &a, double b);
inline SlMatrix2x2 operator-(const SlMatrix2x2 &a, double b);
inline SlMatrix2x2 operator*(const SlMatrix2x2 &a, double b);
inline SlMatrix2x2 operator/(const SlMatrix2x2 &a, double b);

inline SlMatrix3x3 operator+(double a, const SlMatrix3x3 &b);
inline SlMatrix3x3 operator-(double a, const SlMatrix3x3 &b);
inline SlMatrix3x3 operator*(double a, const SlMatrix3x3 &b);
inline SlMatrix3x3 operator/(double a, const SlMatrix3x3 &b);


inline SlMatrix2x2 operator+(const SlMatrix2x2 &a, const SlMatrix2x2 &b);
inline SlMatrix2x2 operator-(const SlMatrix2x2 &a, const SlMatrix2x2 &b);
// These are component-wise operations.

inline SlMatrix2x2 operator*(const SlMatrix2x2 &a, const SlMatrix2x2 &b); 
inline SlVector2   operator*(const SlMatrix2x2 &a, const SlVector2   &b);
inline SlVector2   operator*(const SlVector2   &a, const SlMatrix2x2 &b);
// Matrix multiplication.  Vector returning routines treat the vector
// argument like a row or column vector as appropriate.

//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Global funtions for class SlMatrix2x2

inline double determinant(const SlMatrix2x2 &a);

inline double trace(const SlMatrix2x2 &a);

inline SlMatrix2x2 transpose(const SlMatrix2x2 &a);
/*  */ SlMatrix2x2   inverse(const SlMatrix2x2 &a);

inline SlMatrix2x2 diagonal(const SlVector2 &a);

void SlSymetricEigenDecomp(const SlMatrix2x2 &a,
			   /* */ SlVector2   &vals,
			   /* */ SlMatrix2x2 &vecs);
// Given a SYMETRIC matrix, a, it returns the Eigen Values in vals,
// and the Eigen Vectors as the COLUMNS of vecs.  The vectors are
// normalized.  The values/vectors are not sorted in any order.


inline SlMatrix2x2 SlMatrixProduct(const SlVector2 a,const SlVector2 b);
// Returns the matrix given by a * aT  [ m_ij = a_i * b_j ]

inline SlMatrix2x2 SlEigenMatrix(const SlVector2 a);
// Returns SlMatrixProduct(a,a)/|a| [ m_ij = a_i * a_j / (a_k * a_k)].
// Note that the (only) eigen vector of the matrix is a, and the
// matrix is symetric.  (Involves a call to sqrt.)


bool decomposeF(const SlMatrix3x3 &F, SlMatrix3x3 &U, SlVector3 &Fhat, SlMatrix3x3 &V);


std::istream &operator>>(std::istream &strm,      SlMatrix2x2 &m);
std::ostream &operator<<(std::ostream &strm,const SlMatrix2x2 &m);


//-------------------------------------------------------------------
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

#include <cmath>
#include <cstdlib>
#include "slVector.H"


//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Inline implementation for Methods from SlMatrix3x3


inline SlMatrix3x3::SlMatrix3x3() {
  data[0][0] = 1;
  data[0][1] = 0;
  data[0][2] = 0;

  data[1][0] = 0;
  data[1][1] = 1;
  data[1][2] = 0;

  data[2][0] = 0;
  data[2][1] = 0;
  data[2][2] = 1;  
}

inline SlMatrix3x3::SlMatrix3x3(double x) {
  data[0][0] = x;
  data[0][1] = x;
  data[0][2] = x;

  data[1][0] = x;
  data[1][1] = x;
  data[1][2] = x;

  data[2][0] = x;
  data[2][1] = x;
  data[2][2] = x;
}

inline SlMatrix3x3::SlMatrix3x3(const SlMatrix3x3 &that) {
  data[0][0] = that.data[0][0];
  data[0][1] = that.data[0][1];
  data[0][2] = that.data[0][2];

  data[1][0] = that.data[1][0];
  data[1][1] = that.data[1][1];
  data[1][2] = that.data[1][2];

  data[2][0] = that.data[2][0];
  data[2][1] = that.data[2][1];
  data[2][2] = that.data[2][2];
};

inline SlMatrix3x3::SlMatrix3x3(double xx, double xy, double xz,
				double yx, double yy, double yz,
				double zx, double zy, double zz) {
  data[0][0] = xx;
  data[0][1] = xy;
  data[0][2] = xz;

  data[1][0] = yx;
  data[1][1] = yy;
  data[1][2] = yz;

  data[2][0] = zx;
  data[2][1] = zy;
  data[2][2] = zz;
}

inline SlMatrix3x3::SlMatrix3x3(const double data[3][3], SlMatrix3x3::ordering o) {
  switch (o) {
  case (SlMatrix3x3::row):
    this->data[0][0] = data[0][0];
    this->data[0][1] = data[0][1];
    this->data[0][2] = data[0][2];

    this->data[1][0] = data[1][0];
    this->data[1][1] = data[1][1];
    this->data[1][2] = data[1][2];

    this->data[2][0] = data[2][0];
    this->data[2][1] = data[2][1];
    this->data[2][2] = data[2][2];
    break;
  case (SlMatrix3x3::col):
    this->data[0][0] = data[0][0];
    this->data[0][1] = data[1][0];
    this->data[0][2] = data[2][0];

    this->data[1][0] = data[0][1];
    this->data[1][1] = data[1][1];
    this->data[1][2] = data[2][1];

    this->data[2][0] = data[0][2];
    this->data[2][1] = data[1][2];
    this->data[2][2] = data[2][2];
    break;
  default:
    abort();
  }
}


inline SlMatrix3x3::SlMatrix3x3(const SlVector3 data[3], SlMatrix3x3::ordering o) {
  switch (o) {
  case (SlMatrix3x3::row):
    this->data[0][0] = data[0][0];
    this->data[0][1] = data[0][1];
    this->data[0][2] = data[0][2];

    this->data[1][0] = data[1][0];
    this->data[1][1] = data[1][1];
    this->data[1][2] = data[1][2];

    this->data[2][0] = data[2][0];
    this->data[2][1] = data[2][1];
    this->data[2][2] = data[2][2];
    break;
  case (SlMatrix3x3::col):
    this->data[0][0] = data[0][0];
    this->data[0][1] = data[1][0];
    this->data[0][2] = data[2][0];

    this->data[1][0] = data[0][1];
    this->data[1][1] = data[1][1];
    this->data[1][2] = data[2][1];

    this->data[2][0] = data[0][2];
    this->data[2][1] = data[1][2];
    this->data[2][2] = data[2][2];
    break;
  default:
    abort();
  }
}


inline SlMatrix3x3::SlMatrix3x3(const SlVector3 a,
				const SlVector3 b,
				const SlVector3 c, SlMatrix3x3::ordering o) {
  switch (o) {
  case (SlMatrix3x3::row):
    this->data[0][0] = a[0];
    this->data[0][1] = a[1];
    this->data[0][2] = a[2];

    this->data[1][0] = b[0];
    this->data[1][1] = b[1];
    this->data[1][2] = b[2];

    this->data[2][0] = c[0];
    this->data[2][1] = c[1];
    this->data[2][2] = c[2];
    break;
  case (SlMatrix3x3::col):
    this->data[0][0] = a[0];
    this->data[0][1] = b[0];
    this->data[0][2] = c[0];

    this->data[1][0] = a[1];
    this->data[1][1] = b[1];
    this->data[1][2] = c[1];

    this->data[2][0] = a[2];
    this->data[2][1] = b[2];
    this->data[2][2] = c[2];
    break;
  default:
    abort();
  }
}


//-------------------------------------------------------------------


inline double &SlMatrix3x3::operator()(unsigned int r,unsigned int c) {
   DB_CHECK(r<3);
   DB_CHECK(c<3);
   return data[r][c];
}

inline double SlMatrix3x3::operator()(unsigned int r,unsigned int c) const { 
   DB_CHECK(r<3);
   DB_CHECK(c<3);
   return data[r][c];
}


//-------------------------------------------------------------------


inline SlMatrix3x3 &SlMatrix3x3::set(const double d) {
  return (*this)=d;
}

inline SlMatrix3x3 &SlMatrix3x3::operator=(const double d) {
  data[0][0] = d;
  data[0][1] = d;
  data[0][2] = d;

  data[1][0] = d;
  data[1][1] = d;
  data[1][2] = d;

  data[2][0] = d;
  data[2][1] = d;
  data[2][2] = d;
  return (*this);
};


//-------------------------------------------------------------------

inline SlMatrix3x3 &SlMatrix3x3::set(const SlMatrix3x3 &that) {
  return (*this)=that;
}

inline SlMatrix3x3 &SlMatrix3x3::operator=(const SlMatrix3x3 &that) {
  data[0][0] = that.data[0][0];
  data[0][1] = that.data[0][1];
  data[0][2] = that.data[0][2];

  data[1][0] = that.data[1][0];
  data[1][1] = that.data[1][1];
  data[1][2] = that.data[1][2];

  data[2][0] = that.data[2][0];
  data[2][1] = that.data[2][1];
  data[2][2] = that.data[2][2];
  return (*this);
};


//--------------------------------------------------------------------

inline SlMatrix3x3 &SlMatrix3x3::set(double xx, double xy, double xz,
																		 double yx, double yy, double yz,
																		 double zx, double zy, double zz) {
  data[0][0] = xx;
  data[0][1] = xy;
  data[0][2] = xz;

  data[1][0] = yx;
  data[1][1] = yy;
  data[1][2] = yz;

  data[2][0] = zx;
  data[2][1] = zy;
  data[2][2] = zz;
	return (*this);
}


//-------------------------------------------------------------------

inline int SlMatrix3x3::operator==(double d) const {
  return  ( (data[0][0] == d) &&
	    (data[0][1] == d) &&
	    (data[0][2] == d) &&
	    (data[1][0] == d) &&
	    (data[1][1] == d) &&
	    (data[1][2] == d) &&
	    (data[2][0] == d) &&
	    (data[2][1] == d) &&
	    (data[2][2] == d) );
}

inline int SlMatrix3x3::operator!=(double d) const {
  return  ( (data[0][0] != d) ||
	    (data[0][1] != d) ||
	    (data[0][2] != d) ||
	    (data[1][0] != d) ||
	    (data[1][1] != d) ||
	    (data[1][2] != d) ||
	    (data[2][0] != d) ||
	    (data[2][1] != d) ||
	    (data[2][2] != d) );
}
  
inline int SlMatrix3x3::operator==(const SlMatrix3x3 &that)const {
  return ( (data[0][0] == that.data[0][0]) &&
	   (data[0][1] == that.data[0][1]) &&
	   (data[0][2] == that.data[0][2]) &&
	   (data[1][0] == that.data[1][0]) &&
	   (data[1][1] == that.data[1][1]) &&
	   (data[1][2] == that.data[1][2]) &&
	   (data[2][0] == that.data[2][0]) &&
	   (data[2][1] == that.data[2][1]) &&
	   (data[2][2] == that.data[2][2]) );
}

inline int SlMatrix3x3::operator!=(const SlMatrix3x3 &that)const {
  return ( (data[0][0] != that.data[0][0]) ||
	   (data[0][1] != that.data[0][1]) ||
	   (data[0][2] != that.data[0][2]) ||
	   (data[1][0] != that.data[1][0]) ||
	   (data[1][1] != that.data[1][1]) ||
	   (data[1][2] != that.data[1][2]) ||
	   (data[2][0] != that.data[2][0]) ||
	   (data[2][1] != that.data[2][1]) ||
	   (data[2][2] != that.data[2][2]) );
}


//-------------------------------------------------------------------

inline SlMatrix3x3 &SlMatrix3x3::operator+=(double d) {
  data[0][0] += d; data[1][0] += d; data[2][0] += d;
  data[0][1] += d; data[1][1] += d; data[2][1] += d;
  data[0][2] += d; data[1][2] += d; data[2][2] += d;
  return (*this);
}

inline SlMatrix3x3 &SlMatrix3x3::operator-=(double d) {
  data[0][0] -= d; data[1][0] -= d; data[2][0] -= d;
  data[0][1] -= d; data[1][1] -= d; data[2][1] -= d;
  data[0][2] -= d; data[1][2] -= d; data[2][2] -= d;
  return (*this);
}

inline SlMatrix3x3 &SlMatrix3x3::operator*=(double d) {
  data[0][0] *= d; data[1][0] *= d; data[2][0] *= d;
  data[0][1] *= d; data[1][1] *= d; data[2][1] *= d;
  data[0][2] *= d; data[1][2] *= d; data[2][2] *= d;
  return (*this);
}

inline SlMatrix3x3 &SlMatrix3x3::operator/=(double d) {
  data[0][0] /= d; data[1][0] /= d; data[2][0] /= d;
  data[0][1] /= d; data[1][1] /= d; data[2][1] /= d;
  data[0][2] /= d; data[1][2] /= d; data[2][2] /= d;
  return (*this);
}

//-------------------------------------------------------------------

inline SlMatrix3x3 &SlMatrix3x3::operator+=(const SlMatrix3x3 &that) {
  return componentAdd(that);
}
  
inline SlMatrix3x3 &SlMatrix3x3::operator-=(const SlMatrix3x3 &that) {
  return componentSub(that);
}

//-------------------------------------------------------------------

inline SlMatrix3x3 &SlMatrix3x3::componentAdd (const SlMatrix3x3 &that) {
  data[0][0] += that(0,0); data[1][0] += that(1,0); data[2][0] += that(2,0);
  data[0][1] += that(0,1); data[1][1] += that(1,1); data[2][1] += that(2,1);
  data[0][2] += that(0,2); data[1][2] += that(1,2); data[2][2] += that(2,2);
  return (*this);
}

inline SlMatrix3x3 &SlMatrix3x3::componentSub (const SlMatrix3x3 &that) {
  data[0][0] -= that(0,0); data[1][0] -= that(1,0); data[2][0] -= that(2,0);
  data[0][1] -= that(0,1); data[1][1] -= that(1,1); data[2][1] -= that(2,1);
  data[0][2] -= that(0,2); data[1][2] -= that(1,2); data[2][2] -= that(2,2);
  return (*this);
}

inline SlMatrix3x3 &SlMatrix3x3::componentMult(const SlMatrix3x3 &that) {
  data[0][0] *= that(0,0); data[1][0] *= that(1,0); data[2][0] *= that(2,0);
  data[0][1] *= that(0,1); data[1][1] *= that(1,1); data[2][1] *= that(2,1);
  data[0][2] *= that(0,2); data[1][2] *= that(1,2); data[2][2] *= that(2,2);
  return (*this);
}

inline SlMatrix3x3 &SlMatrix3x3::componentDiv (const SlMatrix3x3 &that) {
  data[0][0] /= that(0,0); data[1][0] /= that(1,0); data[2][0] /= that(2,0);
  data[0][1] /= that(0,1); data[1][1] /= that(1,1); data[2][1] /= that(2,1);
  data[0][2] /= that(0,2); data[1][2] /= that(1,2); data[2][2] /= that(2,2);
  return (*this);
}


//-------------------------------------------------------------------

inline SlMatrix3x3 &SlMatrix3x3::setIdentity() {
  data[0][0] = 1; data[0][1] = 0; data[0][2] = 0;
  data[1][0] = 0; data[1][1] = 1; data[1][2] = 0;
  data[2][0] = 0; data[2][1] = 0; data[2][2] = 1;
  return (*this);
};

inline SlMatrix3x3 &SlMatrix3x3::inplaceTranspose() {
  double tmp;
  tmp = data[0][1]; data[0][1] = data[1][0]; data[1][0] = tmp;
  tmp = data[0][2]; data[0][2] = data[2][0]; data[2][0] = tmp;
  tmp = data[1][2]; data[1][2] = data[2][1]; data[2][1] = tmp;
  return (*this);
}

//-------------------------------------------------------------------

inline unsigned int SlMatrix3x3::cycleAxis(unsigned int axis, int direction) {
  return SlVector3::cycleAxis(axis,direction);
}

#define sqr(x)(x*x)
inline double i2(const SlMatrix3x3 &m) {
  return ( sqr(m(0,0)) + sqr(m(0,1)) + sqr(m(0,2)) +
					 sqr(m(1,0)) + sqr(m(1,1)) + sqr(m(1,2)) +
					 sqr(m(2,0)) + sqr(m(2,1)) + sqr(m(2,2)) );
}
#undef sqr

//Frobenius norm
inline double norm(const SlMatrix3x3 &x) {
  // square root of the sums of the squares of each entry in the matrix
  return sqrt(i2(x));
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Inline implementation for Methods from SlMatrix2x2


inline SlMatrix2x2::SlMatrix2x2() {
  data[0][0] = 1;
  data[0][1] = 0;

  data[1][0] = 0;
  data[1][1] = 1;
}

inline SlMatrix2x2::SlMatrix2x2(double x) {
  data[0][0] = x;
  data[0][1] = x;

  data[1][0] = x;
  data[1][1] = x;
}


inline SlMatrix2x2::SlMatrix2x2(double xx, double xy, 
				double yx, double yy) {
  data[0][0] = xx;
  data[0][1] = xy;
  data[1][0] = yx;
  data[1][1] = yy;
};



inline SlMatrix2x2::SlMatrix2x2(const SlMatrix2x2 &that) {
  data[0][0] = that.data[0][0];
  data[0][1] = that.data[0][1];

  data[1][0] = that.data[1][0];
  data[1][1] = that.data[1][1];
};

inline SlMatrix2x2::SlMatrix2x2(const double data[2][2], SlMatrix2x2::ordering o) {
  switch (o) {
  case (SlMatrix2x2::row):
    this->data[0][0] = data[0][0];
    this->data[0][1] = data[0][1];

    this->data[1][0] = data[1][0];
    this->data[1][1] = data[1][1];
    break;
  case (SlMatrix2x2::col):
    this->data[0][0] = data[0][0];
    this->data[0][1] = data[1][0];

    this->data[1][0] = data[0][1];
    this->data[1][1] = data[1][1];
  default:
    abort();
  }
}


inline SlMatrix2x2::SlMatrix2x2(const SlVector2 data[2], SlMatrix2x2::ordering o) {
  switch (o) {
  case (SlMatrix2x2::row):
    this->data[0][0] = data[0][0];
    this->data[0][1] = data[0][1];

    this->data[1][0] = data[1][0];
    this->data[1][1] = data[1][1];
    break;
  case (SlMatrix2x2::col):
    this->data[0][0] = data[0][0];
    this->data[0][1] = data[1][0];

    this->data[1][0] = data[0][1];
    this->data[1][1] = data[1][1];
    break;
  default:
    abort();
  }
}


inline SlMatrix2x2::SlMatrix2x2(const SlVector2 a,
				const SlVector2 b, SlMatrix2x2::ordering o) {
  switch (o) {
  case (SlMatrix2x2::row):
    this->data[0][0] = a[0];
    this->data[0][1] = a[1];

    this->data[1][0] = b[0];
    this->data[1][1] = b[1];
    break;
  case (SlMatrix2x2::col):
    this->data[0][0] = a[0];
    this->data[0][1] = b[0];

    this->data[1][0] = a[1];
    this->data[1][1] = b[1];
    break;
  default:
    abort();
  }
}


//-------------------------------------------------------------------


inline double &SlMatrix2x2::operator()(unsigned int r,unsigned int c) {
   DB_CHECK(r<2);
   DB_CHECK(c<2);
   return data[r][c];
}

inline double SlMatrix2x2::operator()(unsigned int r,unsigned int c) const { 
   DB_CHECK(r<2);
   DB_CHECK(c<2);
   return data[r][c];
}


//-------------------------------------------------------------------


inline SlMatrix2x2 &SlMatrix2x2::set(const double d) {
  return (*this)=d;
}

inline SlMatrix2x2 &SlMatrix2x2::operator=(const double d) {
  data[0][0] = d;
  data[0][1] = d;

  data[1][0] = d;
  data[1][1] = d;
  return (*this);
};


//-------------------------------------------------------------------

inline SlMatrix2x2 &SlMatrix2x2::set(const SlMatrix2x2 &that) {
  return (*this)=that;
}

inline SlMatrix2x2 &SlMatrix2x2::operator=(const SlMatrix2x2 &that) {
  data[0][0] = that.data[0][0];
  data[0][1] = that.data[0][1];

  data[1][0] = that.data[1][0];
  data[1][1] = that.data[1][1];
  return (*this);
};



//-------------------------------------------------------------------

inline int SlMatrix2x2::operator==(double d) const {
  return  ( (data[0][0] == d) &&
	    (data[0][1] == d) &&
	    (data[1][0] == d) &&
	    (data[1][1] == d) );
}

inline int SlMatrix2x2::operator!=(double d) const {
  return  ( (data[0][0] != d) ||
	    (data[0][1] != d) ||
	    (data[1][0] != d) ||
	    (data[1][1] != d) );
}
  
inline int SlMatrix2x2::operator==(const SlMatrix2x2 &that)const {
  return ( (data[0][0] == that.data[0][0]) &&
	   (data[0][1] == that.data[0][1]) &&
	   (data[1][0] == that.data[1][0]) &&
	   (data[1][1] == that.data[1][1]) );
}

inline int SlMatrix2x2::operator!=(const SlMatrix2x2 &that)const {
  return ( (data[0][0] != that.data[0][0]) ||
	   (data[0][1] != that.data[0][1]) ||
	   (data[1][0] != that.data[1][0]) ||
	   (data[1][1] != that.data[1][1]) );
}


//-------------------------------------------------------------------

inline SlMatrix2x2 &SlMatrix2x2::operator+=(double d) {
  data[0][0] += d; data[1][0] += d; 
  data[0][1] += d; data[1][1] += d; 
  return (*this);
}

inline SlMatrix2x2 &SlMatrix2x2::operator-=(double d) {
  data[0][0] -= d; data[1][0] -= d; 
  data[0][1] -= d; data[1][1] -= d; 
  return (*this);
}

inline SlMatrix2x2 &SlMatrix2x2::operator*=(double d) {
  data[0][0] *= d; data[1][0] *= d; 
  data[0][1] *= d; data[1][1] *= d; 
  return (*this);
}

inline SlMatrix2x2 &SlMatrix2x2::operator/=(double d) {
  data[0][0] /= d; data[1][0] /= d; 
  data[0][1] /= d; data[1][1] /= d; 
  return (*this);
}

//-------------------------------------------------------------------

inline SlMatrix2x2 &SlMatrix2x2::operator+=(const SlMatrix2x2 &that) {
  return componentAdd(that);
}
  
inline SlMatrix2x2 &SlMatrix2x2::operator-=(const SlMatrix2x2 &that) {
  return componentSub(that);
}

//-------------------------------------------------------------------

inline SlMatrix2x2 &SlMatrix2x2::componentAdd (const SlMatrix2x2 &that) {
  data[0][0] += that(0,0); data[1][0] += that(1,0); 
  data[0][1] += that(0,1); data[1][1] += that(1,1); 
  return (*this);
}

inline SlMatrix2x2 &SlMatrix2x2::componentSub (const SlMatrix2x2 &that) {
  data[0][0] -= that(0,0); data[1][0] -= that(1,0); 
  data[0][1] -= that(0,1); data[1][1] -= that(1,1); 
  return (*this);
}

inline SlMatrix2x2 &SlMatrix2x2::componentMult(const SlMatrix2x2 &that) {
  data[0][0] *= that(0,0); data[1][0] *= that(1,0); 
  data[0][1] *= that(0,1); data[1][1] *= that(1,1); 
  return (*this);
}

inline SlMatrix2x2 &SlMatrix2x2::componentDiv (const SlMatrix2x2 &that) {
  data[0][0] /= that(0,0); data[1][0] /= that(1,0); 
  data[0][1] /= that(0,1); data[1][1] /= that(1,1); 
  return (*this);
}


//-------------------------------------------------------------------

inline SlMatrix2x2 &SlMatrix2x2::setIdentity() {
  data[0][0] = 1; data[0][1] = 0; 
  data[1][0] = 0; data[1][1] = 1; 
  return (*this);
};

inline SlMatrix2x2 &SlMatrix2x2::inplaceTranspose() {
  double tmp;
  tmp = data[0][1]; data[0][1] = data[1][0]; data[1][0] = tmp;
  return (*this);
}

//-------------------------------------------------------------------

inline unsigned int SlMatrix2x2::cycleAxis(unsigned int axis, int direction) {
  return SlVector2::cycleAxis(axis,direction);
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Inline implementation of global functions for SlMatrix3x3

inline SlMatrix3x3 operator+(const SlMatrix3x3 &a,double b) {
  return (SlMatrix3x3(a)+=b);
}

inline SlMatrix3x3 operator-(const SlMatrix3x3 &a,double b) {
  return (SlMatrix3x3(a)-=b);
}

inline SlMatrix3x3 operator*(const SlMatrix3x3 &a,double b) {
  return (SlMatrix3x3(a)*=b);
}
 
inline SlMatrix3x3 operator/(const SlMatrix3x3 &a,double b) {
  return (SlMatrix3x3(a)/=b);  
}


//-------------------------------------------------------------------

inline SlMatrix3x3 operator+(double a, const SlMatrix3x3 &b) {
  return b+a;
}

inline SlMatrix3x3 operator-(double a, const SlMatrix3x3 &b) {
  return SlMatrix3x3(a-b(0,0),a-b(0,1),a-b(0,2),
		     a-b(1,0),a-b(1,1),a-b(1,2),
		     a-b(2,0),a-b(2,1),a-b(2,2));
}

inline SlMatrix3x3 operator*(double a, const SlMatrix3x3 &b) {
  return b*a;
}
 
inline SlMatrix3x3 operator/(double a, const SlMatrix3x3 &b) {
  return SlMatrix3x3(a/b(0,0),a/b(0,1),a/b(0,2),
		     a/b(1,0),a/b(1,1),a/b(1,2),
		     a/b(2,0),a/b(2,1),a/b(2,2));
}


//-------------------------------------------------------------------

inline SlMatrix3x3 operator+(const SlMatrix3x3 &a,const SlMatrix3x3 &b) {
  return (SlMatrix3x3(a)+=b);
}
 
inline SlMatrix3x3 operator-(const SlMatrix3x3 &a,const SlMatrix3x3 &b) {
  return (SlMatrix3x3(a)-=b);
}
  
//-------------------------------------------------------------------

inline SlMatrix3x3 operator*(const SlMatrix3x3 &a,const SlMatrix3x3 &b) {
  SlMatrix3x3 tmp(a.data[0][0] * b.data[0][0] + a.data[0][1] * b.data[1][0] + a.data[0][2] * b.data[2][0],
		     a.data[0][0] * b.data[0][1] + a.data[0][1] * b.data[1][1] + a.data[0][2] * b.data[2][1],
		     a.data[0][0] * b.data[0][2] + a.data[0][1] * b.data[1][2] + a.data[0][2] * b.data[2][2],
		     a.data[1][0] * b.data[0][0] + a.data[1][1] * b.data[1][0] + a.data[1][2] * b.data[2][0],
		     a.data[1][0] * b.data[0][1] + a.data[1][1] * b.data[1][1] + a.data[1][2] * b.data[2][1],
		     a.data[1][0] * b.data[0][2] + a.data[1][1] * b.data[1][2] + a.data[1][2] * b.data[2][2],
		     a.data[2][0] * b.data[0][0] + a.data[2][1] * b.data[1][0] + a.data[2][2] * b.data[2][0],
		     a.data[2][0] * b.data[0][1] + a.data[2][1] * b.data[1][1] + a.data[2][2] * b.data[2][1],
		     a.data[2][0] * b.data[0][2] + a.data[2][1] * b.data[1][2] + a.data[2][2] * b.data[2][2]);
  return tmp;
}

inline SlVector3 operator*(const SlMatrix3x3 &a,const SlVector3 &b) {
  return SlVector3(b.data[0]*a.data[0][0] + b.data[1]*a.data[0][1] + b.data[2]*a.data[0][2],
		   b.data[0]*a.data[1][0] + b.data[1]*a.data[1][1] + b.data[2]*a.data[1][2],
		   b.data[0]*a.data[2][0] + b.data[1]*a.data[2][1] + b.data[2]*a.data[2][2]);
}

inline SlVector3 transmult(const SlMatrix3x3 &a,const SlVector3 &b) {
  return SlVector3(b.data[0]*a.data[0][0] + b.data[1]*a.data[1][0] + b.data[2]*a.data[2][0],
		   b.data[0]*a.data[0][1] + b.data[1]*a.data[1][1] + b.data[2]*a.data[2][1],
		   b.data[0]*a.data[0][2] + b.data[1]*a.data[1][2] + b.data[2]*a.data[2][2]);
}

inline SlVector3 operator*(const SlVector3 &a,const SlMatrix3x3 &b) {
  return SlVector3(a[0]*b(0,0) + a[1]*b(1,0) + a[2]*b(2,0),
		   a[0]*b(0,1) + a[1]*b(1,1) + a[2]*b(2,1),
		   a[0]*b(0,2) + a[1]*b(1,2) + a[2]*b(2,2));
}


//-------------------------------------------------------------------  

inline double determinant(const SlMatrix3x3 &a) {
  return ( a(0,0) * ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) +
	   a(0,1) * ( a(1,2) * a(2,0) - a(1,0) * a(2,2) ) +
	   a(0,2) * ( a(1,0) * a(2,1) - a(1,1) * a(2,0) ) );
}


inline double trace(const SlMatrix3x3 &a) {
  return ( a(0,0) + a(1,1) + a(2,2) );
}



inline SlMatrix3x3 transpose(const SlMatrix3x3 &a) {
  return SlMatrix3x3(a(0,0),a(1,0),a(2,0),a(0,1),a(1,1),a(2,1),a(0,2),a(1,2),a(2,2));
}
  

inline SlMatrix3x3 diagonal(const SlVector3 &a) {
  SlMatrix3x3 tmp(0.0);
  tmp(0,0) = a(0);
  tmp(1,1) = a(1);
  tmp(2,2) = a(2);
  return tmp;
}
  



//-------------------------------------------------------------------

inline SlMatrix3x3 SlMatrixProduct(const SlVector3 a,const SlVector3 b) {
  return SlMatrix3x3( a[0] * b[0] ,  a[0] * b[1] ,  a[0] * b[2] , 
		      a[1] * b[0] ,  a[1] * b[1] ,  a[1] * b[2] , 
		      a[2] * b[0] ,  a[2] * b[1] ,  a[2] * b[2] );
}
		     
inline SlMatrix3x3 SlEigenMatrix(const SlVector3 a) {
  double m = mag(a);
  if (m==0) {
    return 0.0;
  }else{
    return SlMatrixProduct(a,a*(1.0/m));
  }
}



//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
// Inline implementation of global functions for SlMatrix2x2

inline SlMatrix2x2 operator+(const SlMatrix2x2 &a,double b) {
  return (SlMatrix2x2(a)+=b);
}

inline SlMatrix2x2 operator-(const SlMatrix2x2 &a,double b) {
  return (SlMatrix2x2(a)-=b);
}

inline SlMatrix2x2 operator*(const SlMatrix2x2 &a,double b) {
  return (SlMatrix2x2(a)*=b);
}
 
inline SlMatrix2x2 operator/(const SlMatrix2x2 &a,double b) {
  return (SlMatrix2x2(a)/=b);  
}


//-------------------------------------------------------------------

inline SlMatrix2x2 operator+(double a, const SlMatrix2x2 &b) {
return b+a;
}

inline SlMatrix2x2 operator-(double a, const SlMatrix2x2 &b) {
  return SlMatrix2x2(a-b(0,1),a-b(0,1),
		     a-b(1,1),a-b(1,1));
}

inline SlMatrix2x2 operator*(double a, const SlMatrix2x2 &b) {
  return b*a;
}
 
inline SlMatrix2x2 operator/(double a, const SlMatrix2x2 &b) {
  return SlMatrix2x2(a/b(0,1),a/b(0,1),
		     a/b(1,1),a/b(1,1));
}


//-------------------------------------------------------------------

inline SlMatrix2x2 operator+(const SlMatrix2x2 &a,const SlMatrix2x2 &b) {
  return (SlMatrix2x2(a)+=b);
}
 
inline SlMatrix2x2 operator-(const SlMatrix2x2 &a,const SlMatrix2x2 &b) {
  return (SlMatrix2x2(a)-=b);
}
  
  
//-------------------------------------------------------------------

inline SlMatrix2x2 operator*(const SlMatrix2x2 &a,const SlMatrix2x2 &b) {
  SlMatrix2x2 tmp(a);
  tmp.inplaceMultPost(b);
  return tmp;
}

inline SlVector2 operator*(const SlMatrix2x2 &a,const SlVector2 &b) {
  return SlVector2(b[0]*a(0,0) + b[1]*a(0,1),
		   b[0]*a(1,0) + b[1]*a(1,1));
}

inline SlVector2 operator*(const SlVector2 &a,const SlMatrix2x2 &b) {
  return SlVector2(a[0]*b(0,0) + a[1]*b(1,0),
		   a[0]*b(0,1) + a[1]*b(1,1));
}


//-------------------------------------------------------------------  

inline double determinant(const SlMatrix2x2 &a) {
  return ( a(0,0) * a(1,1) - a(0,1) * a(1,0) );
}

inline double trace(const SlMatrix2x2 &a) {
  return ( a(0,0) + a(1,1) );
}



inline SlMatrix2x2 transpose(const SlMatrix2x2 &a) {
  SlMatrix2x2 tmp(a);
  tmp.inplaceTranspose();
  return tmp;
}
 
inline SlMatrix2x2 diagonal(const SlVector2 &a) {
  SlMatrix2x2 tmp(0.0);
  tmp(0,0) = a(0);
  tmp(1,1) = a(1);
  return tmp;
}
  

//-------------------------------------------------------------------

inline SlMatrix2x2 SlMatrixProduct(const SlVector2 a,const SlVector2 b) {
  return SlMatrix2x2( a[0] * b[0] ,  a[0] * b[1] ,
		      a[1] * b[0] ,  a[1] * b[1] );
}
		     
inline SlMatrix2x2 SlEigenMatrix(const SlVector2 a) {
  double m = mag(a);
  if (m==0) {
    return 0.0;
  }else{
    return SlMatrixProduct(a,a*(1.0/m));
  }
}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#endif
