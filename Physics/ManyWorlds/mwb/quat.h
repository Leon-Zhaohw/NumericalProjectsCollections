
#ifndef __QUAT_H__
#define __QUAT_H__

#ifdef USE_NOVODEX
#include "NxQuat.h"
#endif

#include "twigg/vlutil.h"

#include <boost/array.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <numeric>
#include <cmath>

// order xyzw
template <typename Real>
class Quaternion
{
private:
	typedef boost::array<Real, 4> Array;
	Array elements;

	void init( const vl::Vec3d& scaledAxis )
	{
		double angle = vl::len( scaledAxis );
		if( angle <= std::numeric_limits<double>::epsilon() )
		{
			std::fill( elements.begin(), elements.begin() + 3, 
				boost::numeric_cast<Real>(0) );
			this->elements[3] = boost::numeric_cast<Real>(1);		
		}
		else
		{
			vl::Vec3d normAxis = scaledAxis;
			normAxis /= angle;

			double s = sin( angle/2.0 );
			for( size_t i = 0; i < 3; ++i )
				this->elements[i] = normAxis[i]*s;

			this->elements[3] = cos( angle / 2.0 );
		}
	}

public:
	enum Order
	{
		ORDER_WXYZ,
		ORDER_XYZW
	};

	Quaternion()
	{
		std::fill( elements.begin(), elements.begin() + 3, 
			boost::numeric_cast<Real>(0) );
		this->elements[3] = boost::numeric_cast<Real>(1);
	}

	template <typename Real2>
	Quaternion( const Quaternion<Real2>& q )
	{
		this->elements[0] = q.x();
		this->elements[1] = q.y();
		this->elements[2] = q.z();
		this->elements[3] = q.w();
	}

#ifdef USE_NOVODEX
	Quaternion( const NxQuat& q )
	{
		q.getXYZW( &elements[0] );
	}
#endif

	template <typename Real2>
	Quaternion( const Real2* other, Order order )
	{
		switch( order )
		{
		case ORDER_WXYZ:
			this->elements[3] = other[0];
			std::copy( other + 1, other + 4, this->elements.begin() );
			break;
		case ORDER_XYZW:
			std::copy( other, other + 4, this->elements.begin() );
			break;
		}
	}

	Quaternion( const vl::Vec3d& scaledAxis )
	{
		init( scaledAxis );
	}

	Quaternion( const vl::Vec3f& scaledAxis )
	{
		init( toVec3d(scaledAxis) );
	}

	template <typename Real2>
	Quaternion( const vl::Vec3f& vec, Real2 w )
	{
		std::copy( vec.Ref(), vec.Ref()+3, this->elements.begin() );
		this->elements[3] = w;
	}

	template <typename Real2>
	Quaternion( const vl::Vec3d& vec, Real2 w )
	{
		std::copy( vec.Ref(), vec.Ref()+3, this->elements.begin() );
		this->elements[3] = w;
	}

	Quaternion( Real x, Real y, Real z, Real w )
	{
		this->elements[0] = x;
		this->elements[1] = y;
		this->elements[2] = z;
		this->elements[3] = w;
	}

	// Matrix had better be in SO(3):
	// from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
	template <typename Mat>
	Quaternion( const Mat& mat )
	{
		Real& x = this->elements[0];
		Real& y = this->elements[1];
		Real& z = this->elements[2];
		Real& w = this->elements[3];

		w = sqrt( std::max( Real(0), Real(1) + Real(mat[0][0]) + Real(mat[1][1]) + Real(mat[2][2]) ) ) / Real(2);
		x = sqrt( std::max( Real(0), Real(1) + Real(mat[0][0]) - Real(mat[1][1]) - Real(mat[2][2]) ) ) / Real(2);
		y = sqrt( std::max( Real(0), Real(1) - Real(mat[0][0]) + Real(mat[1][1]) - Real(mat[2][2]) ) ) / Real(2); 
		z = sqrt( std::max( Real(0), Real(1) - Real(mat[0][0]) - Real(mat[1][1]) + Real(mat[2][2]) ) ) / Real(2); 
		x = _copysign( x, mat[2][1] - mat[1][2] );
		y = _copysign( y, mat[0][2] - mat[2][0] );
		z = _copysign( z, mat[1][0] - mat[0][1] );
	}

	template <typename Real2>
    void operator*=( const Quaternion<Real2>& rhs )
	{
		Quaternion<Real> result = (*this) * rhs;
		this->elements.swap( result.elements );
	}

	void operator*=( const Real& rhs )
	{
		for( typename Array::size_type i = 0; i < 4; ++i )
			this->elements[i] *= rhs;
	}

	void operator/=( const Real& rhs )
	{
		for( typename Array::size_type i = 0; i < 4; ++i )
			this->elements[i] /= rhs;
	}

    template <typename Real2>
    Real dot( const Quaternion<Real2>& rhs ) const
    {
        return std::inner_product( elements.begin(), elements.end(), rhs.elements.begin(), Real(0) );
    }

	template <typename Real2>
	Quaternion<Real> operator*( const Quaternion<Real2>& rhs ) const
	{
		//(Q1 * Q2).w = (w1w2 - x1x2 - y1y2 - z1z2)
		//(Q1 * Q2).x = (w1x2 + x1w2 + y1z2 - z1y2)
		//(Q1 * Q2).y = (w1y2 - x1z2 + y1w2 + z1x2)
		//(Q1 * Q2).z = (w1z2 + x1y2 - y1x2 + z1w2)
		Quaternion<Real> result;
		Real x = this->w() * rhs.x()
			+ this->x() * rhs.w()
			+ this->y() * rhs.z()
			- this->z() * rhs.y();
		Real y = this->w() * rhs.y()
			- this->x() * rhs.z()
			+ this->y() * rhs.w()
			+ this->z() * rhs.x();
		Real z = this->w() * rhs.z()
			+ this->x() * rhs.y()
			- this->y() * rhs.x()
			+ this->z() * rhs.w();
		Real w = this->w() * rhs.w()
			- this->x() * rhs.x()
			- this->y() * rhs.y()
			- this->z() * rhs.z();
		return Quaternion<Real>( x, y, z, w );
	}

	Quaternion<Real> operator*( const Real& rhs )
	{
		Quaternion<Real> result( *this );
		result *= rhs;
		return result;
	}

	Quaternion<Real> operator/( const Real& rhs )
	{
		Quaternion<Real> result( *this );
		result /= rhs;
		return result;
	}

	template <typename Real2>
	void operator+=( const Quaternion<Real2>& rhs )
	{
		for( typename Array::size_type i = 0; i < 4; ++i )
			this->elements[i] += rhs.elements[i];
	}

	template <typename Real2>
	Quaternion<Real> operator+( const Quaternion<Real2>& rhs ) const
	{
		Quaternion<Real> result;
		result += rhs;
		return result;
	}

	bool operator==( const Quaternion<Real>& rhs ) const
	{
        return std::equal( this->elements.begin(), this->elements.end(), rhs.elements.begin() );
	}

    bool operator!=( const Quaternion<Real>& rhs ) const
    {
        return !(*this == rhs);
    }

#ifdef USE_NOVODEX
	NxQuat toNxQuat() const
	{
		NxQuat result( 
			NxVec3(this->elements[0], this->elements[1], this->elements[2]),
			this->elements[3] );
		return result;
	}
#endif

	vl::Mat3d toRotMatd() const
	{
		return toRotMat<vl::Mat3d, double>();
	}

	vl::Mat3f toRotMatf() const
	{
		return toRotMat<vl::Mat3f, float>();
	}

	Real w() const
	{
		return elements[3];
	}

	Real x() const
	{
		return elements[0];
	}

	Real y() const
	{
		return elements[1];
	}

	Real z() const
	{
		return elements[2];
	}

	void conjugate()
	{
		for( typename Array::size_type i = 0; i < 3; ++i )
			this->elements[i] = -this->elements[i];
	}

	Real magnitudeSquared() const
	{
		Real result = 
			std::inner_product( this->elements.begin(), this->elements.end(),
				this->elements.begin(), boost::numeric_cast<Real>(0) );
		return result;
	}

	Real magnitude() const
	{
		return sqrt( this->magnitudeSquared() );
	}

	void normalize()
	{
		Real mag = this->magnitude();
		for( typename Array::size_type i = 0; i < 4; ++i )
			this->elements[i] /= mag;
	}

	Real argument() const
	{
		return acos( this->w() / this->magnitude() );
	}

	template <typename Real2, typename VecType>
	void toAxisAngle( VecType& axis, Real2& angle ) const
	{
		// quaternion must be normalized
		Quaternion<Real2> q( *this );
		q.normalize();

		angle = boost::numeric_cast<Real2>(2) * acos(q.w());
		Real2 s = sqrt( boost::numeric_cast<Real2>(1) - q.w()*q.w() );
		if( s <= std::numeric_limits<Real2>::epsilon() )
		{
			axis[0] = boost::numeric_cast<Real2>(0);
			axis[1] = boost::numeric_cast<Real2>(1);
			axis[2] = boost::numeric_cast<Real2>(0);
		}
		else
		{
			axis[0] = q.x() / s;
			axis[1] = q.y() / s;
			axis[2] = q.z() / s;
		}

	}

	boost::array<Real, 4> get( Order order ) const
	{
		switch( order )
		{
		case ORDER_WXYZ:
			{
				boost::array<Real, 4> result = {{ w(), x(), y(), z() }};
				return result;
			}
		case ORDER_XYZW:
			{
				return this->elements;
			}
		}

		assert( false );
		return this->elements;
	}

private:
	template <typename Mat, typename Float>
	Mat toRotMat() const
	{
		Float xx      = x() * x();
		Float xy      = x() * y();
		Float xz      = x() * z();
		Float xw      = x() * w();

		Float yy      = y() * y();
		Float yz      = y() * z();
		Float yw      = y() * w();

		Float zz      = z() * z();
		Float zw      = z() * w();

		Mat mat;
		mat[0][0]  = 1 - 2 * ( yy + zz );
		mat[0][1]  =     2 * ( xy - zw );
		mat[0][2]  =     2 * ( xz + yw );

		mat[1][0]  =     2 * ( xy + zw );
		mat[1][1]  = 1 - 2 * ( xx + zz );
		mat[1][2]  =     2 * ( yz - xw );

		mat[2][0]  =     2 * ( xz - yw );
		mat[2][1]  =     2 * ( yz + xw );
		mat[2][2]  = 1 - 2 * ( xx + yy );
		return mat;
	}

    template <typename Real>
    friend Quaternion<Real> slerp( const Quaternion<Real>& q1, const Quaternion<Real>& q2, const Real t );

    template <typename Real>
    friend Quaternion<Real> lerp( const Quaternion<Real>& q1, const Quaternion<Real>& q2, const Real t );
};

template <typename Real>
Quaternion<Real> conjugate( const Quaternion<Real>& q )
{
	Quaternion<Real> result( q );
	result.conjugate();
	return result;
}

template <typename Real>
Quaternion<Real> inverse( const Quaternion<Real>& q )
{
	Quaternion<Real> result( q );
	result.conjugate();
	result /= result.magnitudeSquared();
	return result;
}

template <typename Real>
Real scalar( const Quaternion<Real>& q )
{
	return q.w();
}

template <typename Real>
Real arg( const Quaternion<Real>& q )
{
	return q.argument();
}

template <typename Real>
bool isBad( const Quaternion<Real>& q )
{
	return (isBad(q.w()) || isBad(q.x()) || isBad(q.y()) || isBad(q.z()));
}

template <typename Real>
Quaternion<Real> lerp( const Quaternion<Real>& q1, const Quaternion<Real>& q2, const Real t )
{
    assert( !isBad(q1) );
    assert( !isBad(q2) );
    assert( !isBad(t) );

    assert( t >= Real(0) && t <= Real(1) );

    Quaternion<Real> q3;
    Real inv_t = Real(1) - t;
    for( size_t i =0; i < 4; ++i )
        q3.elements[i] = t*q2.elements[i] + inv_t*q1.elements[i];
    q3.normalize();

    return q3;
}


template <typename Real>
Quaternion<Real> slerp( const Quaternion<Real>& q1, const Quaternion<Real>& q2, const Real t )
{
    assert( !isBad(q1) );
    assert( !isBad(q2) );
    assert( !isBad(t) );

    // code from here:
    // http://number-none.com/product/Understanding%20Slerp,%20Then%20Not%20Using%20It/
    Real dot = q1.dot( q2 );
    const Real threshold = Real(1) - 1e-5;
    if( dot > threshold )
    {
        // linearly interpolate
        Quaternion<Real> result;
        for( size_t i = 0; i < 4; ++i )
            result.elements[i] = result.elements[i] + t*(q2.elements[i] - q1.elements[i]);

        result.normalize();
        return result;
    }

    dot = std::max( std::min( dot, Real(1) ), Real(-1) );
    Real theta0 = acos( dot );
    Real theta = theta0 * t;

    assert( !isBad(theta0) );
    assert( !isBad(theta) );

    Quaternion<Real> q3;
    for( size_t i = 0; i < 4; ++i )
        q3.elements[i] = q2.elements[i] - dot*q1.elements[i];
    q3.normalize();
    assert( !isBad(q3) );

    Real cosTheta = cos(theta), sinTheta = sin(theta);
    for( size_t i = 0; i < 4; ++i )
        q3.elements[i] = q1.elements[i]*cosTheta + q3.elements[i]*sinTheta;

    assert( !isBad(q3) );
    return q3;
}

template <typename Real>
Real angleBetween( const Quaternion<Real>& left, const Quaternion<Real>& right )
{
    Quaternion<Real> change = inverse(left) * right;
    return arg(change);
}

void testQuat();

#endif

