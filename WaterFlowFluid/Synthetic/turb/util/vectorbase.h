/******************************************************************************
 *  DDF Fluid solver with Turbulence extensions
 *
 *  copyright 2009 Nils Thuerey, Tobias Pfaff
 * 
 *  DDF is free software, distributed under the GNU General Public License (GPL v2).
 *  See the file COPYING for more information.
 *
 * Vector class
 *
 *****************************************************************************/


#ifndef DDF_BASICVEC_H
#define DDF_BASICVEC_H

// get rid of windos min/max defines
#ifdef WIN32
#define NOMINMAX
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// if min/max are still around...
#ifdef WIN32
#undef min
#undef max
#endif

// use which fp-precision? 1=float, 2=double
#ifndef FLOATINGPOINT_PRECISION
#if DDF_DEBUG==1
#define FLOATINGPOINT_PRECISION 2
#else // DDF_DEBUG==1
#define FLOATINGPOINT_PRECISION 1
#endif // DDF_DEBUG==1
#endif

// dimension of solver&problem
#ifndef DDF_DIM
#define DDF_DIM 3
#endif

// windos, hardcoded limits for now...
// for e.g. MSVC compiler...
// some of these defines can be needed
// for linux systems as well (e.g. FLT_MAX)
#ifndef __FLT_MAX__
#	ifdef FLT_MAX  // try to use it instead
#		define __FLT_MAX__ FLT_MAX
#	else // FLT_MAX
#		define __FLT_MAX__ 3.402823466e+38f
#	endif // FLT_MAX
#endif // __FLT_MAX__
#ifndef __DBL_MAX__
#	ifdef DBL_MAX // try to use it instead
#		define __DBL_MAX__ DBL_MAX
#	else // DBL_MAX
#		define __DBL_MAX__ 1.7976931348623158e+308
#	endif // DBL_MAX
#endif // __DBL_MAX__

#ifndef M_PI
#	define M_PI 3.1415926536
#	define M_E  2.7182818284
#endif


namespace DDF {


// basic inlined vector class
template<class Scalar>
class ntlVector3Dim
{
public:
  // Constructor
  inline ntlVector3Dim();
  // Copy-Constructor
  inline ntlVector3Dim(const ntlVector3Dim<Scalar> &v );
  inline ntlVector3Dim(const float *);
  inline ntlVector3Dim(const double *);
  // construct a vector from one Scalar
  inline ntlVector3Dim(Scalar);
  // construct a vector from three Scalars
  inline ntlVector3Dim(Scalar, Scalar, Scalar);

	// get address of array for OpenGL
	Scalar *getAddress() { return value; }

  // Assignment operator
  inline const ntlVector3Dim<Scalar>& operator=  (const ntlVector3Dim<Scalar>& v);
  // Assignment operator
  inline const ntlVector3Dim<Scalar>& operator=  (Scalar s);
  // Assign and add operator
  inline const ntlVector3Dim<Scalar>& operator+= (const ntlVector3Dim<Scalar>& v);
  // Assign and add operator
  inline const ntlVector3Dim<Scalar>& operator+= (Scalar s);
  // Assign and sub operator
  inline const ntlVector3Dim<Scalar>& operator-= (const ntlVector3Dim<Scalar>& v);
  // Assign and sub operator
  inline const ntlVector3Dim<Scalar>& operator-= (Scalar s);
  // Assign and mult operator
  inline const ntlVector3Dim<Scalar>& operator*= (const ntlVector3Dim<Scalar>& v);
  // Assign and mult operator
  inline const ntlVector3Dim<Scalar>& operator*= (Scalar s);
  // Assign and div operator
  inline const ntlVector3Dim<Scalar>& operator/= (const ntlVector3Dim<Scalar>& v);
  // Assign and div operator
  inline const ntlVector3Dim<Scalar>& operator/= (Scalar s);


  // unary operator
  inline ntlVector3Dim<Scalar> operator- () const;

  // binary operator add
  inline ntlVector3Dim<Scalar> operator+ (const ntlVector3Dim<Scalar>&) const;
  // binary operator add
  inline ntlVector3Dim<Scalar> operator+ (Scalar) const;
  // binary operator sub
  inline ntlVector3Dim<Scalar> operator- (const ntlVector3Dim<Scalar>&) const;
  // binary operator sub
  inline ntlVector3Dim<Scalar> operator- (Scalar) const;
  // binary operator mult
  inline ntlVector3Dim<Scalar> operator* (const ntlVector3Dim<Scalar>&) const;
  // binary operator mult
  inline ntlVector3Dim<Scalar> operator* (Scalar) const;
  // binary operator div
  inline ntlVector3Dim<Scalar> operator/ (const ntlVector3Dim<Scalar>&) const;
  // binary operator div
  inline ntlVector3Dim<Scalar> operator/ (Scalar) const;

  // Projection normal to a vector
  inline ntlVector3Dim<Scalar>	  getOrthogonalntlVector3Dim() const;
  // Project into a plane
  inline const ntlVector3Dim<Scalar>& projectNormalTo(const ntlVector3Dim<Scalar> &v);
  
  inline Scalar min() const { return (x<y) ? ( (x<z) ? x:z ) : ( (y<z) ? y:z); }
  inline Scalar max() const { return (x>y) ? ( (x>z) ? x:z ) : ( (y>z) ? y:z); }

  // minimize
  inline const ntlVector3Dim<Scalar> &minimize(const ntlVector3Dim<Scalar> &);
  // maximize
  inline const ntlVector3Dim<Scalar> &maximize(const ntlVector3Dim<Scalar> &);
  
  // access operator
  inline Scalar& operator[](unsigned int i);
  // access operator
  inline const Scalar& operator[](unsigned int i) const;

  // return absolutes of all components
  inline ntlVector3Dim<Scalar> getAbsolutes() const { return 
  		ntlVector3Dim<Scalar>(fabs(value[0]), fabs(value[1]), fabs(value[2]) );
  };

  // debug output vector to a string
  std::string toString();

	//! actual values
	union {
		struct {
  		Scalar value[3];  
		};
		struct {
  		Scalar x;
  		Scalar y;
  		Scalar z;
		};
		struct {
  		Scalar X;
  		Scalar Y;
  		Scalar Z;
		};
	};

	// expe compatibility functions
	void makeFloor(const ntlVector3Dim<Scalar>& cmp);
	void makeCeil(const ntlVector3Dim<Scalar>& cmp);
	Scalar squaredDistanceTo(const ntlVector3Dim<Scalar>& vec) const;

    // Returns true if the vector's s components are all greater that the ones of the vector it is compared against.
    inline bool operator < ( const ntlVector3Dim<Scalar>& vec ) const; 
    // Returns true if the vector's s components are all greater or equal that the ones of the vector it is compared against.
    inline bool operator <= ( const ntlVector3Dim<Scalar>& vec ) const; 
    // Returns true if the vector's s components are all smaller that the ones of the vector it is compared against.
    inline bool operator > ( const ntlVector3Dim<Scalar>& vec ) const; 
    // Returns true if the vector's s components are all smaller or equal that the ones of the vector it is compared against.
    inline bool operator >= ( const ntlVector3Dim<Scalar>& vec ) const; 

	// Return the maximal component value.
    inline Scalar maxComponent(void) const; 
    // Return the minimal component value.
    inline Scalar minComponent(void) const; 
    // Return the index of the maximal coordinate value.
    inline int maxComponentId(void) const; 
    // Return the index of the minimal coordinate value.
    inline int minComponentId(void) const;

	// zero element
   static const ntlVector3Dim<Scalar> ZERO;

protected:
  
};




// VECTOR_EPSILON is the minimal vector length
// In order to be able to discriminate floating point values near zero, and
// to be sure not to fail a comparison because of roundoff errors, use this
// value as a threshold.  

#if FLOATINGPOINT_PRECISION==1
typedef float Real;
#define FP_REAL_MAX __FLT_MAX__
#define VECTOR_EPSILON (1e-5f)
#else
typedef double Real;
#define FP_REAL_MAX __DBL_MAX__
#define VECTOR_EPSILON (1e-10)
#endif
// as variable in globals.cpp
extern const Real gVecEpsilon;



//------------------------------------------------------------------------------
// VECTOR inline FUNCTIONS
//------------------------------------------------------------------------------



/*************************************************************************
  Constructor.
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>::ntlVector3Dim( void )
{
  value[0] = value[1] = value[2] = 0;
}



/*************************************************************************
  Copy-Constructor.
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>::ntlVector3Dim( const ntlVector3Dim<Scalar> &v )
{
	value[0] = v.value[0];
	value[1] = v.value[1];
	value[2] = v.value[2];
}
	template<class Scalar>
inline ntlVector3Dim<Scalar>::ntlVector3Dim( const float *fvalue)
{
	value[0] = (Scalar)fvalue[0];
	value[1] = (Scalar)fvalue[1];
	value[2] = (Scalar)fvalue[2];
}
	template<class Scalar>
inline ntlVector3Dim<Scalar>::ntlVector3Dim( const double *fvalue)
{
	value[0] = (Scalar)fvalue[0];
	value[1] = (Scalar)fvalue[1];
	value[2] = (Scalar)fvalue[2];
}



/*************************************************************************
  Constructor for a vector from a single Scalar. All components of
  the vector get the same value.
  \param s The value to set
  \return The new vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>::ntlVector3Dim(Scalar s )
{
  value[0]= s;
  value[1]= s;
  value[2]= s;
}


/*************************************************************************
  Constructor for a vector from three Scalars.
  \param s1 The value for the first vector component
  \param s2 The value for the second vector component
  \param s3 The value for the third vector component
  \return The new vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>::ntlVector3Dim(Scalar s1, Scalar s2, Scalar s3)
{
  value[0]= s1;
  value[1]= s2;
  value[2]= s3;
}


/*************************************************************************
  Compute the vector product of two 3D vectors
  \param v Second vector to compute the product with
  \return A new vector with the product values
  */
/*template<class Scalar>
inline ntlVector3Dim<Scalar> 
ntlVector3Dim<Scalar>::operator^( const ntlVector3Dim<Scalar> &v ) const
{
  return ntlVector3Dim<Scalar>(value[1]*v.value[2] - value[2]*v.value[1],
			value[2]*v.value[0] - value[0]*v.value[2],
			value[0]*v.value[1] - value[1]*v.value[0]);
}*/


/*************************************************************************
  Copy a ntlVector3Dim componentwise.
  \param v vector with values to be copied
  \return Reference to self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::operator=( const ntlVector3Dim<Scalar> &v )
{
  value[0] = v.value[0];
  value[1] = v.value[1];
  value[2] = v.value[2];  
  return *this;
}


/*************************************************************************
  Copy a Scalar to each component.
  \param s The value to copy
  \return Reference to self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::operator=(Scalar s)
{
  value[0] = s;
  value[1] = s;
  value[2] = s;  
  return *this;
}


/*************************************************************************
  Add another ntlVector3Dim componentwise.
  \param v vector with values to be added
  \return Reference to self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::operator+=( const ntlVector3Dim<Scalar> &v )
{
  value[0] += v.value[0];
  value[1] += v.value[1];
  value[2] += v.value[2];  
  return *this;
}


/*************************************************************************
  Add a Scalar value to each component.
  \param s Value to add
  \return Reference to self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::operator+=(Scalar s)
{
  value[0] += s;
  value[1] += s;
  value[2] += s;  
  return *this;
}


/*************************************************************************
  Subtract another vector componentwise.
  \param v vector of values to subtract
  \return Reference to self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::operator-=( const ntlVector3Dim<Scalar> &v )
{
  value[0] -= v.value[0];
  value[1] -= v.value[1];
  value[2] -= v.value[2];  
  return *this;
}


/*************************************************************************
  Subtract a Scalar value from each component.
  \param s Value to subtract
  \return Reference to self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::operator-=(Scalar s)
{
  value[0]-= s;
  value[1]-= s;
  value[2]-= s;  
  return *this;
}


/*************************************************************************
  Multiply with another vector componentwise.
  \param v vector of values to multiply with
  \return Reference to self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::operator*=( const ntlVector3Dim<Scalar> &v )
{
  value[0] *= v.value[0];
  value[1] *= v.value[1];
  value[2] *= v.value[2];  
  return *this;
}


/*************************************************************************
  Multiply each component with a Scalar value.
  \param s Value to multiply with
  \return Reference to self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::operator*=(Scalar s)
{
  value[0] *= s;
  value[1] *= s;
  value[2] *= s;  
  return *this;
}


/*************************************************************************
  Divide by another ntlVector3Dim componentwise.
  \param v vector of values to divide by
  \return Reference to self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::operator/=( const ntlVector3Dim<Scalar> &v )
{
  value[0] /= v.value[0];
  value[1] /= v.value[1];
  value[2] /= v.value[2];  
  return *this;
}


/*************************************************************************
  Divide each component by a Scalar value.
  \param s Value to divide by
  \return Reference to self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::operator/=(Scalar s)
{
  value[0] /= s;
  value[1] /= s;
  value[2] /= s;
  return *this;
}


//------------------------------------------------------------------------------
// unary operators
//------------------------------------------------------------------------------


/*************************************************************************
  Build componentwise the negative this vector.
  \return The new (negative) vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>
ntlVector3Dim<Scalar>::operator-() const
{
  return ntlVector3Dim<Scalar>(-value[0], -value[1], -value[2]);
}



//------------------------------------------------------------------------------
// binary operators
//------------------------------------------------------------------------------


/*************************************************************************
  Build a vector with another vector added componentwise.
  \param v The second vector to add
  \return The sum vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>
ntlVector3Dim<Scalar>::operator+( const ntlVector3Dim<Scalar> &v ) const
{
  return ntlVector3Dim<Scalar>(value[0]+v.value[0],
			value[1]+v.value[1],
			value[2]+v.value[2]);
}


/*************************************************************************
  Build a vector with a Scalar value added to each component.
  \param s The Scalar value to add
  \return The sum vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>
ntlVector3Dim<Scalar>::operator+(Scalar s) const
{
  return ntlVector3Dim<Scalar>(value[0]+s,
			value[1]+s,
			value[2]+s);
}


/*************************************************************************
  Build a vector with another vector subtracted componentwise.
  \param v The second vector to subtract
  \return The difference vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>
ntlVector3Dim<Scalar>::operator-( const ntlVector3Dim<Scalar> &v ) const
{
  return ntlVector3Dim<Scalar>(value[0]-v.value[0],
			value[1]-v.value[1],
			value[2]-v.value[2]);
}


/*************************************************************************
  Build a vector with a Scalar value subtracted componentwise.
  \param s The Scalar value to subtract
  \return The difference vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>
ntlVector3Dim<Scalar>::operator-(Scalar s ) const
{
  return ntlVector3Dim<Scalar>(value[0]-s,
			value[1]-s,
			value[2]-s);
}



/*************************************************************************
  Build a vector with another vector multiplied by componentwise.
  \param v The second vector to muliply with
  \return The product vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>
ntlVector3Dim<Scalar>::operator*( const ntlVector3Dim<Scalar>& v) const
{
  return ntlVector3Dim<Scalar>(value[0]*v.value[0],
			value[1]*v.value[1],
			value[2]*v.value[2]);
}


/*************************************************************************
  Build a ntlVector3Dim with a Scalar value multiplied to each component.
  \param s The Scalar value to multiply with
  \return The product vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>
ntlVector3Dim<Scalar>::operator*(Scalar s) const
{
  return ntlVector3Dim<Scalar>(value[0]*s, value[1]*s, value[2]*s);
}


// allow multiplications of the form: v2 = 3 * v1
template<class Scalar>
inline ntlVector3Dim<Scalar>
operator*(Scalar s, ntlVector3Dim<Scalar> v) 
{
  return ntlVector3Dim<Scalar>(v.value[0]*s, v.value[1]*s, v.value[2]*s);
}



/*************************************************************************
  Build a vector divided componentwise by another vector.
  \param v The second vector to divide by
  \return The ratio vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>
ntlVector3Dim<Scalar>::operator/(const ntlVector3Dim<Scalar>& v) const
{
  return ntlVector3Dim<Scalar>(value[0]/v.value[0],
			value[1]/v.value[1],
			value[2]/v.value[2]);
}



/*************************************************************************
  Build a vector divided componentwise by a Scalar value.
  \param s The Scalar value to divide by
  \return The ratio vector
  */
template<class Scalar>
inline ntlVector3Dim<Scalar>
ntlVector3Dim<Scalar>::operator/(Scalar s) const
{
  return ntlVector3Dim<Scalar>(value[0]/s,
			value[1]/s,
			value[2]/s);
}





/*************************************************************************
  Get a particular component of the vector.
  \param i Number of Scalar to get
  \return Reference to the component
  */
template<class Scalar>
inline Scalar&
ntlVector3Dim<Scalar>::operator[]( unsigned int i )
{
  return value[i];
}


/*************************************************************************
  Get a particular component of a constant vector.
  \param i Number of Scalar to get
  \return Reference to the component
  */
template<class Scalar>
inline const Scalar&
ntlVector3Dim<Scalar>::operator[]( unsigned int i ) const
{
  return value[i];
}



//------------------------------------------------------------------------------
// BLITZ compatibility functions
//------------------------------------------------------------------------------



/*************************************************************************
  Compute the scalar product with another vector.
  \param v The second vector to work with
  \return The value of the scalar product
  */
template<class Scalar>
inline Scalar dot(const ntlVector3Dim<Scalar> &t, const ntlVector3Dim<Scalar> &v )
{
  //return t.value[0]*v.value[0] + t.value[1]*v.value[1] + t.value[2]*v.value[2];
  return ((t[0]*v[0]) + (t[1]*v[1]) + (t[2]*v[2]));
}


/*************************************************************************
  Calculate the cross product of this and another vector
 */
template<class Scalar>
inline ntlVector3Dim<Scalar> cross(const ntlVector3Dim<Scalar> &t, const ntlVector3Dim<Scalar> &v)
{
  ntlVector3Dim<Scalar> cp( 
			((t[1]*v[2]) - (t[2]*v[1])),
		  ((t[2]*v[0]) - (t[0]*v[2])),
		  ((t[0]*v[1]) - (t[1]*v[0])) );
  return cp;
}




/*************************************************************************
  Compute a vector that is orthonormal to self. Nothing else can be assumed
  for the direction of the new vector.
  \return The orthonormal vector
  */
template<class Scalar>
ntlVector3Dim<Scalar>
ntlVector3Dim<Scalar>::getOrthogonalntlVector3Dim() const
{
	// Determine the  component with max. absolute value
	int maxIndex= (fabs(value[0]) > fabs(value[1])) ? 0 : 1;
	maxIndex= (fabs(value[maxIndex]) > fabs(value[2])) ? maxIndex : 2;

	/*************************************************************************
	  Choose another axis than the one with max. component and project
	  orthogonal to self
	 */
	ntlVector3Dim<Scalar> vec(0.0);
	vec[(maxIndex+1)%3]= 1;
	vec.normalize();
	vec.projectNormalTo(this->getNormalized());
	return vec;
}


/*************************************************************************
  Projects the vector into a plane normal to the given vector, which must
  have unit length. Self is modified.
  \param v The plane normal
  \return The projected vector
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar>&
ntlVector3Dim<Scalar>::projectNormalTo(const ntlVector3Dim<Scalar> &v)
{
	Scalar sprod = dot(*this,v);
	value[0]= value[0] - v.value[0] * sprod;
	value[1]= value[1] - v.value[1] * sprod;
	value[2]= value[2] - v.value[2] * sprod;  
	return *this;
}



//------------------------------------------------------------------------------
// Other helper functions
//------------------------------------------------------------------------------



/*************************************************************************
  Minimize the vector, i.e. set each entry of the vector to the minimum
  of both values.
  \param pnt The second vector to compare with
  \return Reference to the modified self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar> &
ntlVector3Dim<Scalar>::minimize(const ntlVector3Dim<Scalar> &pnt)
{
  for (unsigned int i = 0; i < 3; i++)
    value[i] = MIN(value[i],pnt[i]);
  return *this;
}



/*************************************************************************
  Maximize the vector, i.e. set each entry of the vector to the maximum
  of both values.
  \param pnt The second vector to compare with
  \return Reference to the modified self
  */
template<class Scalar>
inline const ntlVector3Dim<Scalar> &
ntlVector3Dim<Scalar>::maximize(const ntlVector3Dim<Scalar> &pnt)
{
  for (unsigned int i = 0; i < 3; i++)
    value[i] = MAX(value[i],pnt[i]);
  return *this;
}






/************************************************************************/
// HELPER FUNCTIONS, independent of implementation
/************************************************************************/

#define VECTOR_TYPE ntlVector3Dim<Scalar>


/*************************************************************************
  Compute the length (norm) of the vector.
  \return The value of the norm
  */
template<class Scalar>
inline Scalar norm( const VECTOR_TYPE &v)
{
  Scalar l = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  return (fabs(l-1.) < VECTOR_EPSILON*VECTOR_EPSILON) ? 1. : sqrt(l);
}

// for e.g. min max operator
inline Real normHelper(const ntlVector3Dim<Real> &v) {
	return norm(v);
}	
inline Real normHelper(const Real &v) {
	return (0. < v) ? v : -v ; 
}	
inline Real normHelper(const int &v) {
	return (0 < v) ? (Real)(v) : (Real)(-v) ; 
}	


/*************************************************************************
  Same as getNorm but doesnt sqrt  
 */
template<class Scalar>
inline Scalar normNoSqrt( const VECTOR_TYPE &v) {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}


/*************************************************************************
  Compute a normalized vector based on this vector.
  \return The new normalized vector
  */
template<class Scalar>
inline VECTOR_TYPE getNormalized( const VECTOR_TYPE &v)
{
	Scalar l = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	if (fabs(l-1.) < VECTOR_EPSILON*VECTOR_EPSILON)
		return v; /* normalized "enough"... */
	else if (l > VECTOR_EPSILON*VECTOR_EPSILON)
	{
		Scalar fac = 1./sqrt(l);
		return VECTOR_TYPE(v[0]*fac, v[1]*fac, v[2]*fac);
	}
	else
		return VECTOR_TYPE((Scalar)0);
}


/*************************************************************************
  Compute the norm of the vector and normalize it.
  \return The value of the norm
 */
	template<class Scalar>
inline Scalar normalize( VECTOR_TYPE &v) 
{
	Scalar norm;
	Scalar l = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];  
	if (fabs(l-1.) < VECTOR_EPSILON*VECTOR_EPSILON) {
		norm = 1.;
	} else if (l > VECTOR_EPSILON*VECTOR_EPSILON) {
		norm = sqrt(l);
		Scalar fac = 1./norm;
		v[0] *= fac;
		v[1] *= fac;
		v[2] *= fac; 
	} else {
		v[0]= v[1]= v[2]= 0;
		norm = 0.;
	}
	return (Scalar)norm;
}

/*************************************************************************
  Stable vector to angle conversion
  \return unique Angles in the of phi=[0,2PI], th=[0,PI]
 */
	template<class Scalar>
inline void vecToAngle(const VECTOR_TYPE &v, Scalar& phi, Scalar& theta) 
{
	if (fabs(v.y) < VECTOR_EPSILON)
		theta = M_PI/2;
	else if (fabs(v.x) < VECTOR_EPSILON && fabs(v.z) < VECTOR_EPSILON )
		theta = (v.y>=0) ? 0:M_PI;
	else
		theta = atan(sqrt(v.x*v.x+v.z*v.z)/v.y);
	if (theta<0) theta+=M_PI;
		
	if (fabs(v.x) < VECTOR_EPSILON)
		phi = M_PI/2;
	else
		phi = atan(v.z/v.x);	
	if (phi<0) phi+=M_PI;
	if (fabs(v.z) < VECTOR_EPSILON)
		phi = (v.x>=0) ? 0 : M_PI;
	else if (v.z < 0)
		phi += M_PI;
}


/*************************************************************************
  Compute a vector, that is self (as an incoming
  vector) reflected at a surface with a distinct normal vector. Note
  that the normal is reversed, if the scalar product with it is positive.
  \param n The surface normal
  \return The new reflected vector
  */
template<class Scalar>
inline VECTOR_TYPE reflectVector(const VECTOR_TYPE &t, const VECTOR_TYPE &n) 
{
  VECTOR_TYPE nn= (dot(t, n) > 0.0) ? (n*-1.0) : n;
  return ( t - nn * (2.0 * dot(nn, t)) );
}



/*************************************************************************
 * My own refraction calculation
 * Taken from Glassner's book, section 5.2 (Heckberts method)
 */
template<class Scalar>
inline VECTOR_TYPE refractVector(const VECTOR_TYPE &t, const VECTOR_TYPE &normal, Scalar nt, Scalar nair, int &refRefl) 
{
	Scalar eta = nair / nt;
	Scalar n = -dot(t, normal);
	Scalar tt = 1.0 + eta*eta* (n*n-1.0);
	if(tt<0.0) {
		// we have total reflection!
		refRefl = 1;
	} else {
		// normal reflection
		tt = eta*n - sqrt(tt);
		return( t*eta + normal*tt );
	}
	return t;
}


/*************************************************************************
  Test two ntlVector3Dims for equality based on the equality of their
  values within a small threshold.
  \param c The second vector to compare
  \return TRUE if both are equal
  \sa getEpsilon()
  */
template<class Scalar>
inline bool equal(const VECTOR_TYPE &v, const VECTOR_TYPE &c)
{
  return (ABS(v[0]-c[0]) + 
	  ABS(v[1]-c[1]) + 
	  ABS(v[2]-c[2]) < VECTOR_EPSILON);
}


/*************************************************************************
 * Assume this vector is an RGB color, and convert it to HSV
 */
template<class Scalar>
inline void rgbToHsv( VECTOR_TYPE &V )
{
	Scalar h=0,s=0,v=0;
	Scalar maxrgb, minrgb, delta;
	// convert to hsv...
	maxrgb = V[0];
	int maxindex = 1;
	if(V[2] > maxrgb){ maxrgb = V[2]; maxindex = 2; }
	if(V[1] > maxrgb){ maxrgb = V[1]; maxindex = 3; }
	minrgb = V[0];
	if(V[2] < minrgb) minrgb = V[2];
	if(V[1] < minrgb) minrgb = V[1];

	v = maxrgb;
	delta = maxrgb-minrgb;

	if(maxrgb > 0) s = delta/maxrgb;
	else s = 0;

	h = 0;
	if(s > 0) {
		if(maxindex == 1) {
			h = ((V[1]-V[2])/delta)  + 0.0; }
		if(maxindex == 2) {
			h = ((V[2]-V[0])/delta)  + 2.0; }
		if(maxindex == 3) {
			h = ((V[0]-V[1])/delta)  + 4.0; }
		h *= 60.0;
		if(h < 0.0) h += 360.0;
	}

	V[0] = h;
	V[1] = s;
	V[2] = v;
}

/*************************************************************************
 * Assume this vector is HSV and convert to RGB
 */
template<class Scalar>
inline void hsvToRgb( VECTOR_TYPE &V )
{
	Scalar h = V[0], s = V[1], v = V[2];
	Scalar r=0,g=0,b=0;
	Scalar p,q,t, fracth;
	int floorh;
	// ...and back to rgb
	if(s == 0) {
		r = g = b = v; }
	else {
		h /= 60.0;
		floorh = (int)h;
		fracth = h - floorh;
		p = v * (1.0 - s);
		q = v * (1.0 - (s * fracth));
		t = v * (1.0 - (s * (1.0 - fracth)));
		switch (floorh) {
		case 0: r = v; g = t; b = p; break;
		case 1: r = q; g = v; b = p; break;
		case 2: r = p; g = v; b = t; break;
		case 3: r = p; g = q; b = v; break;
		case 4: r = t; g = p; b = v; break;
		case 5: r = v; g = p; b = q; break;
		}
	}

	V[0] = r;
	V[1] = g;
	V[2] = b;
}

//------------------------------------------------------------------------------
// STREAM FUNCTIONS
//------------------------------------------------------------------------------



//! global string for formatting vector output in utilities.cpp
extern const char *globVecFormatStr;

/*************************************************************************
  Outputs the object in human readable form using the format
  [x,y,z]
  */
template<class Scalar>
std::ostream&
operator<<( std::ostream& os, const DDF::ntlVector3Dim<Scalar>& i )
{
	char buf[256];
	snprintf(buf,256,globVecFormatStr, (double)i[0],(double)i[1],(double)i[2]);
	os << std::string(buf); 
	//os << '[' << i[0] << ", " << i[1] << ", " << i[2] << ']';
	return os;
}



/*************************************************************************
  Reads the contents of the object from a stream using the same format
  as the output operator.
  */
template<class Scalar>
std::istream&
operator>>( std::istream& is, DDF::ntlVector3Dim<Scalar>& i )
{
	char c;
	char dummy[3];
	is >> c >> i[0] >> dummy >> i[1] >> dummy >> i[2] >> c;
	return is;
}

// helper function for output
template<class Scalar> std::string ntlVector3Dim<Scalar>::toString() {
	char buf[256];
	snprintf(buf,256,globVecFormatStr, (double)(*this)[0],(double)(*this)[1],(double)(*this)[2]);
	return std::string(buf);
}

// helper function for output
template<class Vec> bool intVecIsEqual(Vec a, Vec b) {
	return a[0]==b[0] && a[1]==b[1] && a[2]==b[2];
}

/**************************************************************************/
// typedefs!
/**************************************************************************/

// a 3D vector with double precision
typedef ntlVector3Dim<double>  nVec3d; 

// a 3D vector with single precision
typedef ntlVector3Dim<float>   nVec3f; 

// a 3D integer vector
typedef ntlVector3Dim<int>     nVec3i; 

/* convert int,float and double vectors */
template<class T> inline nVec3i vec2I(T v) { return nVec3i((int)v[0],(int)v[1],(int)v[2]); }
template<class T> inline nVec3i vec2I(T v0, T v1, T v2) { return nVec3i((int)v0,(int)v1,(int)v2); }
template<class T> inline nVec3d vec2D(T v) { return nVec3d(v[0],v[1],v[2]); }
template<class T> inline nVec3i vec2D(T v0, T v1, T v2) { return nVec3d((double)v0,(double)v1,(double)v2); }
template<class T> inline nVec3f vec2F(T v) { return nVec3f(v[0],v[1],v[2]); }
template<class T> inline nVec3i vec2F(T v0, T v1, T v2) { return nVec3f((float)v0,(float)v1,(float)v2); }
template<class T> inline nVec3i vecround(T v) { return nVec3i((int)round(v[0]),(int)round(v[1]),(int)round(v[2])); }



/************************************************************************/
// default vector typing


// a 3D vector for graphics output, typically float?
typedef ntlVector3Dim<Real>  Vec3; 

/* convert to Real vec */
template<class T> inline Vec3   vec2R(T v) { return Vec3(v[0],v[1],v[2]); }


/* get minimal vector length value that can be discriminated.  */
inline Real getVecEpsilon() { return (Real)VECTOR_EPSILON; }

// expe compatibility functions

//template<class Scalar>
//inline XXX ntlVector3Dim<Scalar>::

//--------------------------------------------------------------------------------
template<class Scalar>
inline Scalar ntlVector3Dim<Scalar>::squaredDistanceTo(const ntlVector3Dim<Scalar>& vec) const
{
    Scalar dx,dy,dz;
    dx = x-vec.x;
    dy = y-vec.y;
    dz = z-vec.z;
    return dx*dx + dy*dy + dz*dz;
}
//--------------------------------------------------------------------------------
template<class Scalar>
inline void ntlVector3Dim<Scalar>::makeFloor(const ntlVector3Dim<Scalar>& cmp)
{
    if( cmp.x < x ) x = cmp.x;
    if( cmp.y < y ) y = cmp.y;
    if( cmp.z < z ) z = cmp.z;
}
//--------------------------------------------------------------------------------
template<class Scalar>
inline void ntlVector3Dim<Scalar>::makeCeil(const ntlVector3Dim<Scalar>& cmp)
{
    if( cmp.x > x ) x = cmp.x;
    if( cmp.y > y ) y = cmp.y;
    if( cmp.z > z ) z = cmp.z;
}

//--------------------------------------------------------------------------------
template<class Scalar>
inline bool ntlVector3Dim<Scalar>::operator < ( const ntlVector3Dim<Scalar>& vec ) const
{
    if( x < vec.x && y < vec.y && z < vec.z )
        return true;
    return false;
}
//--------------------------------------------------------------------------------
template<class Scalar>
inline bool ntlVector3Dim<Scalar>::operator <= ( const ntlVector3Dim<Scalar>& vec ) const
{
    if( x <= vec.x && y <= vec.y && z <= vec.z )
        return true;
    return false;
}
//--------------------------------------------------------------------------------
template<class Scalar>
inline bool ntlVector3Dim<Scalar>::operator > ( const ntlVector3Dim<Scalar>& vec ) const
{
    if( x > vec.x && y > vec.y && z > vec.z )
        return true;
    return false;
}
//--------------------------------------------------------------------------------
template<class Scalar>
inline bool ntlVector3Dim<Scalar>::operator >= ( const ntlVector3Dim<Scalar>& vec ) const
{
    if( x >= vec.x && y >= vec.y && z >= vec.z )
        return true;
    return false;
}

//--------------------------------------------------------------------------------
template<class Scalar>
inline Scalar ntlVector3Dim<Scalar>::maxComponent(void) const
{
    return VMAX(*this);
}
//--------------------------------------------------------------------------------
template<class Scalar>
inline Scalar ntlVector3Dim<Scalar>::minComponent(void) const
{
    return VMIN(*this);
}
//--------------------------------------------------------------------------------
template<class Scalar>
inline int ntlVector3Dim<Scalar>::maxComponentId(void) const
{
    if (x<y)
        return (y<z ? 2 : 1);
    else
        return (x<z ? 2 : 0);
}
//--------------------------------------------------------------------------------
template<class Scalar>
inline int ntlVector3Dim<Scalar>::minComponentId(void) const
{
    if (x<y)
        return (x<z ? 0 : 2);
    else
        return (y<Z ? 1 : 2);
}

}; // namespace DDF


#endif /* DDF_BASICVEC_H */
