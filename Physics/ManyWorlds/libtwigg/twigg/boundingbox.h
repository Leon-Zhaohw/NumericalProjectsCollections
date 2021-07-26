#ifndef BOUNDINGBOX_H_INCLUDED
#define BOUNDINGBOX_H_INCLUDED

#include "twigg/util.h"
#include "twigg/vlutil.h"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <algorithm>
#include <boost/array.hpp>

// Borrowed from the ray tracer code.
template <typename VecType, typename RealType>
class Ray
{
public:
	Ray( const VecType& pp, const VecType& dd )
		: position_( pp ), direction_( dd )
	{

		inverseDirection_[0] = static_cast<RealType>(1) / direction_[0];
		inverseDirection_[1] = static_cast<RealType>(1) / direction_[1];
		inverseDirection_[2] = static_cast<RealType>(1) / direction_[2];

		for( vl::Int i = 0; i < 3; ++i )
		{
			if( inverseDirection_[i] < 0.0f )
				sign_[i] = 1;
			else
				sign_[i] = 0;
		}
	}

	VecType at( RealType t ) const
	{
		return position_ + (t*direction_);
	}

	VecType getPosition() const { return position_; }
	VecType getDirection() const { return direction_; }
	VecType getInverseDirection() const { return inverseDirection_; }
	int sign( unsigned int index ) const { return sign_[index]; }

private:
	VecType position_;
	VecType direction_;
	VecType inverseDirection_;
	boost::array<char, 3> sign_;
};

typedef Ray<vl::Vec3f, float> Ray3f;
typedef Ray<vl::Vec3d, double> Ray3d;

template <typename VecType, typename RealType>
class TriangleIntersection
{
public:
	TriangleIntersection( RealType t, const VecType& bary )
		: t_(t), barycentricCoords_(bary) {}

	TriangleIntersection()
		: t_(boost::numeric::bounds<double>::highest()), barycentricCoords_(vl::vl_0) {}

	RealType t() const						{ return t_; }
	VecType barycentricCoordinates() const	{ return barycentricCoords_; }

private:
	RealType t_;
	VecType barycentricCoords_;
};

typedef TriangleIntersection<vl::Vec3f, float> TriangleIntersection3f;
typedef TriangleIntersection<vl::Vec3d, double> TriangleIntersection3d;

template <typename VecType, typename RealType>
std::pair< bool, TriangleIntersection<VecType, RealType> > intersectTriangle(
	const Ray<VecType, RealType>& r, const VecType& a, const VecType& b, const VecType& c )
{
	// @todo check these
	const typename RealType SMALL_EPSILON = std::numeric_limits<RealType>::epsilon();
	const typename RealType epsilon = SMALL_EPSILON;
	const typename RealType NORMAL_EPSILON = std::numeric_limits<RealType>::epsilon() * 10;

	typename VecType p = r.getPosition();
	typename VecType v = r.getDirection();

	typename VecType ab = b - a;
	typename VecType ac = c - a;
	typename VecType ap = p - a;

	typename VecType cv = cross(ab, ac);

	// there exists some bad triangles such that two vertices coincide
	// check this before normalize
	if( sqrlen(cv) < SMALL_EPSILON )
		return std::make_pair(false, TriangleIntersection<VecType, RealType>());
	typename VecType n = vl::norm(cv);

	typename RealType vdotn = dot(v,n);
	if( fabs(vdotn) < NORMAL_EPSILON )
		return std::make_pair(false, TriangleIntersection<VecType, RealType>());

	typename RealType t = -dot(ap,n) / vdotn;

	if( t < epsilon )
		return std::make_pair(false, TriangleIntersection<VecType, RealType>());

	// find k where k is the index of the component
	// of normal vector with greatest absolute value
	typename RealType greatestMag = 0;
	int k = -1;
	for( int j = 0; j < 3; ++j )
	{
		float val = n[j];
		if( val < 0 )
			val *= -1;
		if( val > greatestMag )
		{
			k = j;
			greatestMag = val;
		}
	}

	typename VecType am = ap + t * v;

	typename VecType bary;
	bary[1] = cross(am, ac)[k]/cross(ab, ac)[k];
	bary[2] = cross(ab, am)[k]/cross(ab, ac)[k];
	bary[0] = 1-bary[1]-bary[2];

	if( bary[0] < -SMALL_EPSILON || 
		bary[1] < -SMALL_EPSILON || 
		bary[1] > (1+SMALL_EPSILON) || 
		bary[2] < -SMALL_EPSILON || 
		bary[2] > (1+SMALL_EPSILON) )
		return std::make_pair(false, TriangleIntersection<VecType, RealType>());

	// if we get this far, we have an intersection.  Fill in the info.
	return std::make_pair(true,
		TriangleIntersection<VecType, RealType>( t, bary ));
}

template <typename VecType, typename RealType>
class BoundingBox3
{
public:
	typedef VecType point_type;
	typedef RealType real_type;

	BoundingBox3(const VecType& min, const VecType& max)
	{
		nullify();

		expand( min );
		expand( max );
	}

	BoundingBox3()
	{
		nullify();
	}

	void nullify()
	{
		for( unsigned int i = 0; i < 3; ++i )
			(parameters[0])[i] = boost::numeric::bounds<RealType>::highest();
		for( unsigned int i = 0; i < 3; ++i )
			(parameters[1])[i] = boost::numeric::bounds<RealType>::lowest();
	}

	void expand(const VecType& point)
	{
		assert(!isBad(point));
		parameters[0] = vl::min( point, parameters[0] );
		parameters[1] = vl::max( point, parameters[1] );
	}

	void expand(const BoundingBox3<VecType, RealType>& other)
	{
		parameters[0] = vl::min( parameters[0], other.parameters[0] );
		parameters[1] = vl::max( parameters[1], other.parameters[1] );
	}

	void expand(const RealType amount)
	{
		if( empty() )
			return;

		for( unsigned int iVec = 0; iVec < 3; ++iVec )
			(parameters[0])[iVec] -= amount;

		for( unsigned int iVec = 0; iVec < 3; ++iVec )
			(parameters[1])[iVec] += amount;
	}

	VecType minimum() const
	{
		return parameters[0];
	}

	VecType maximum() const
	{
		return parameters[1];
	}

	RealType area() const
	{
		return 2.0*halfArea();
	}

	RealType halfArea() const
	{
		VecType v3f = parameters[1] - parameters[0];

		return fabs(v3f[0] * v3f[1]) +
			fabs(v3f[1] * v3f[2]) +
			fabs(v3f[2] * v3f[0]);
	}

	bool empty() const
	{
		for( unsigned int i = 0; i < 3; ++i )
		{
			if( (parameters[1])[i] < (parameters[0])[i]  )
				return true;
		}

		return false;
	}

	// if the ray hits the box, put the "t" value of the intersection
	// closest to the origin in tMin and the "t" value of the far intersection
	// in tMax and return true, else return false.
	bool intersect(const Ray<VecType, RealType>& r, RealType& tmin, RealType& tmax,
					RealType t0 = 0.0, RealType t1 = boost::numeric::bounds<RealType>::highest()) const
	{
		// if the ray hits the box, put the "t" value of the intersection
		// closest to the origin in tMin and the "t" value of the far intersection
		// in tMax and return true, else return false.
		// Using code from this paper:
		//   http://www.cs.utah.edu/~awilliam/box/box.pdf
		VecType invDir = r.getInverseDirection();
		VecType origin = r.getPosition();

		tmin = (parameters[r.sign(0)][0] - origin[0]) * invDir[0];
		tmax = (parameters[1-r.sign(0)][0] - origin[0]) * invDir[0];
		float tymin = (parameters[r.sign(1)][1] - origin[1]) * invDir[1];
		float tymax = (parameters[1-r.sign(1)][1] - origin[1]) * invDir[1];
		if ( (tmin > tymax) || (tymin > tmax) )
		return false;
		if (tymin > tmin)
		tmin = tymin;
		if (tymax < tmax)
		tmax = tymax;
		float tzmin = (parameters[r.sign(2)][2] - origin[2]) * invDir[2];
		float tzmax = (parameters[1-r.sign(2)][2] - origin[2]) * invDir[2];
		if ( (tmin > tzmax) || (tzmin > tmax) )
		return false;
		if (tzmin > tmin)
		tmin = tzmin;
		if (tzmax < tmax)
		tmax = tzmax;

		return ( (tmin < t1) && (tmax > t0) );
	}

private:
	boost::array<VecType, 2> parameters;
};

template <typename VecType, typename RealType>
const BoundingBox3<VecType, RealType> unite(const BoundingBox3<VecType, RealType>& bb1, 
											const BoundingBox3<VecType, RealType>& bb2)
{
	VecType v3fMin, v3fMax;

	v3fMin[0] = bb1.minimum()[0] < bb2.minimum()[0] ? bb1.minimum()[0] : bb2.minimum()[0];
	v3fMin[1] = bb1.minimum()[1] < bb2.minimum()[1] ? bb1.minimum()[1] : bb2.minimum()[1];
	v3fMin[2] = bb1.minimum()[2] < bb2.minimum()[2] ? bb1.minimum()[2] : bb2.minimum()[2];

	v3fMax[0] = bb1.maximum()[0] > bb2.maximum()[0] ? bb1.maximum()[0] : bb2.maximum()[0];
	v3fMax[1] = bb1.maximum()[1] > bb2.maximum()[1] ? bb1.maximum()[1] : bb2.maximum()[1];
	v3fMax[2] = bb1.maximum()[2] > bb2.maximum()[2] ? bb1.maximum()[2] : bb2.maximum()[2];

	return BoundingBox3<VecType, RealType>(v3fMin, v3fMax);
}

template <typename VecType, typename RealType>
bool contains(const BoundingBox3<VecType, RealType>& bb3f, 
			  const VecType& v3fPt, 
			  const RealType fEpsilon = 0.0)
{
	return (
		v3fPt[0] >= bb3f.minimum()[0] - fEpsilon &&
		v3fPt[1] >= bb3f.minimum()[1] - fEpsilon &&
		v3fPt[2] >= bb3f.minimum()[2] - fEpsilon &&
		v3fPt[0] <= bb3f.maximum()[0] + fEpsilon &&
		v3fPt[1] <= bb3f.maximum()[1] + fEpsilon &&
		v3fPt[2] <= bb3f.maximum()[2] + fEpsilon);
}

template <typename VecType, typename RealType>
bool intersects(const BoundingBox3<VecType, RealType>& bb1, 
				const BoundingBox3<VecType, RealType>& bb2 )
{
	{
		const VecType bb1Max = bb1.maximum();
		const VecType bb2Min = bb2.minimum();
		if (bb1Max[0] < bb2Min[0] ||
			bb1Max[1] < bb2Min[1] ||
			bb1Max[2] < bb2Min[2])
			return false;
	}
	{
		const VecType bb1Min = bb1.minimum();
		const VecType bb2Max = bb2.maximum();
		if (bb1Min[0] > bb2Max[0] ||
			bb1Min[1] > bb2Max[1] ||
			bb1Min[2] > bb2Max[2] )
			return false;
	}

	return true;
}


template <typename VecType, typename RealType>
bool intersects(const BoundingBox3<VecType, RealType>& bb1, 
				const BoundingBox3<VecType, RealType>& bb2, 
				const RealType fEpsilon )
{
	if (bb1.maximum()[0] + fEpsilon < bb2.minimum()[0] ||
		bb1.maximum()[1] + fEpsilon < bb2.minimum()[1] ||
		bb1.maximum()[2] + fEpsilon < bb2.minimum()[2])
		return false;
	if (bb1.minimum()[0] > bb2.maximum()[0] + fEpsilon ||
		bb1.minimum()[1] > bb2.maximum()[1] + fEpsilon ||
		bb1.minimum()[2] > bb2.maximum()[2] + fEpsilon)
		return false;

	return true;
}

template <typename real_type, typename vec_type>
class BoundingBox2
{
public:
	typedef vec_type point_type;

	BoundingBox2(const vec_type& v2fMin, const vec_type& v2fMax)
		: m_v2fMin(v2fMin), m_v2fMax(v2fMax)
	{
		if (m_v2fMin[0] > m_v2fMax[0])
			m_v2fMin[0] = m_v2fMax[0];
		if (m_v2fMin[1] > m_v2fMax[1])
			m_v2fMin[1] = m_v2fMax[1];
	}

	BoundingBox2()
	{
		nullify();
	}

	void nullify()
	{
		m_v2fMin = point_type(
				boost::numeric::bounds<real_type>::highest(),
				boost::numeric::bounds<real_type>::highest());
		m_v2fMax = point_type(
				boost::numeric::bounds<real_type>::lowest(),
				boost::numeric::bounds<real_type>::lowest());
	}

	bool empty() const
	{
		for( unsigned int i = 0; i < 2; ++i )
		{
			if( m_v2fMax[i] < m_v2fMin[i]  )
				return true;
		}

		return false;
	}

	void expand(const point_type& v2fPt)
	{
		for( unsigned int i = 0; i < 2; ++i )
			m_v2fMin[i] = std::min<real_type>(v2fPt[i], m_v2fMin[i]);

		for( unsigned int i = 0; i < 2; ++i )
			m_v2fMax[i] = std::max<real_type>(v2fPt[i], m_v2fMax[i]);
	}

	void expand(const real_type amount)
	{
		if( empty() )
			return;

		for( unsigned int iVec = 0; iVec < 2; ++iVec )
			m_v2fMin[iVec] -= amount;

		for( unsigned int iVec = 0; iVec < 2; ++iVec )
			m_v2fMax[iVec] += amount;
	}

	void expand(const BoundingBox2<real_type, vec_type>& other)
	{
		m_v2fMin = point_type( std::min<real_type>(m_v2fMin[0], other.m_v2fMin[0]), 
			std::min<real_type>(m_v2fMin[1], other.m_v2fMin[1]) );
		m_v2fMax = point_type( std::max<real_type>(m_v2fMax[0], other.m_v2fMax[0]), 
			std::max<real_type>(m_v2fMax[1], other.m_v2fMax[1]) );

	}

	const point_type& minimum() const
	{
		return this->m_v2fMin;
	}

	const point_type& maximum() const
	{
		return this->m_v2fMax;
	}

	real_type area() const
	{
		if( empty() )
			return 0;

		vec_type v2f = m_v2fMax - m_v2fMin;
		assert( v2f[0] >= 0.0 && v2f[1] >= 0.0 );

		return v2f[0] * v2f[1];
	}

	/*
	???
	real_type circumference() const
	{
		if( empty() )
			return 0.0;

		point_type v2f = m_v2fMax - m_v2fMin;

		return v2f[0] + v2f[1] * 2.0;
	}
	*/

	/*
	???
	real_type atanAspectRatio() const
	{
		point_type v2f = m_v2fMax - m_v2fMin;

		return fabs(atan2(v2f[1], v2f[0]));
	}
	*/


	real_type width() const
	{
		return m_v2fMax[0] - m_v2fMin[0];
	}

	real_type height() const
	{
		return m_v2fMax[1] - m_v2fMin[1];
	}

private:
	point_type m_v2fMin, m_v2fMax;
};

template <typename real_type, typename vec_type>
const BoundingBox2<real_type, vec_type> intersect(
	const BoundingBox2<real_type, vec_type>& bb1, 
	const BoundingBox2<real_type, vec_type>& bb2)
{
	vec_type v2fMin, v2fMax;

	v2fMin[0] = bb1.minimum()[0] > bb2.minimum()[0] ? bb1.minimum()[0] : bb2.minimum()[0];
	v2fMin[1] = bb1.minimum()[1] > bb2.minimum()[1] ? bb1.minimum()[1] : bb2.minimum()[1];

	v2fMax[0] = bb1.maximum()[0] < bb2.maximum()[0] ? bb1.maximum()[0] : bb2.maximum()[0];
	v2fMax[1] = bb1.maximum()[1] < bb2.maximum()[1] ? bb1.maximum()[1] : bb2.maximum()[1];

	return BoundingBox2<real_type, vec_type>(v2fMin, v2fMax);
}

template <typename real_type, typename vec_type>
const BoundingBox2<real_type, vec_type> unite(
	const BoundingBox2<real_type, vec_type>& bb1, 
	const BoundingBox2<real_type, vec_type>& bb2)
{
	vec_type v2fMin, v2fMax;

	v2fMin[0] = bb1.minimum()[0] < bb2.minimum()[0] ? bb1.minimum()[0] : bb2.minimum()[0];
	v2fMin[1] = bb1.minimum()[1] < bb2.minimum()[1] ? bb1.minimum()[1] : bb2.minimum()[1];

	v2fMax[0] = bb1.maximum()[0] > bb2.maximum()[0] ? bb1.maximum()[0] : bb2.maximum()[0];
	v2fMax[1] = bb1.maximum()[1] > bb2.maximum()[1] ? bb1.maximum()[1] : bb2.maximum()[1];

	return BoundingBox2<real_type, vec_type>(v2fMin, v2fMax);
}

template <typename real_type, typename vec_type>
bool contains(const BoundingBox2<real_type, vec_type>& bb2f, 
			  const vec_type& v2fPt, 
			  const real_type fEpsilon = 0.0)
{
	return (
		v2fPt[0] >= bb2f.minimum()[0] - fEpsilon &&
		v2fPt[1] >= bb2f.minimum()[1] - fEpsilon &&
		v2fPt[0] <= bb2f.maximum()[0] + fEpsilon &&
		v2fPt[1] <= bb2f.maximum()[1] + fEpsilon);
}

template <typename real_type, typename vec_type>
bool intersects(const BoundingBox2<real_type, vec_type>& bb1, 
				const BoundingBox2<real_type, vec_type>& bb2, 
				const real_type fEpsilon = 0.0)
{
	if (bb1.maximum()[0] + fEpsilon < bb2.minimum()[0] ||
		bb1.maximum()[1] + fEpsilon < bb2.minimum()[1])
		return false;
	if (bb1.minimum()[0] > bb2.maximum()[0] + fEpsilon ||
		bb1.minimum()[1] > bb2.maximum()[1] + fEpsilon)
		return false;

	return true;
}

template <typename ObjectPtr>
class GenericIntersectionHandler
{
public:
	GenericIntersectionHandler( ObjectPtr obj = ObjectPtr() )
		: intersected_(obj), t_( boost::numeric::bounds<double>::highest() ) {}

	bool operator()( ObjectPtr object, const Ray& ray )
	{
		double t;
		if( object->intersect(ray, t) )
		{
			if( t < t_ )
			{
				t_ = t;
				intersected_ = object;
			}
			return true;
		}
		return false;
	}

	ObjectPtr intersected() const
	{
		return intersected_;
	}

private:
	ObjectPtr intersected_;
	double t_;
};

typedef BoundingBox3<vl::Vec3d, double> BoundingBox3d;
typedef BoundingBox3<vl::Vec3f, float> BoundingBox3f;
typedef BoundingBox2<double, vl::Vec2d> BoundingBox2d;
typedef BoundingBox2<float, vl::Vec2f> BoundingBox2f;

template <typename BoundingBox, typename Transform>
BoundingBox transformBounds( const BoundingBox& box, const Transform& transform )
{
	if( box.empty() )
		return box;

	BoundingBox result;

	const typename BoundingBox::point_type mn = box.minimum();
	const typename BoundingBox::point_type mx = box.maximum();
	for( int i = 0; i < 2; ++i )
		for( int j = 0; j < 2; ++j )
			for( int k = 0; k < 2; ++k )
			{
				typename BoundingBox::point_type point( 
					i == 0 ? mn[0] : mx[0],
					j == 0 ? mn[1] : mx[1],
					k == 0 ? mn[2] : mx[2] );
				result.expand( xform(transform, point) );
			}

	return result;
}


#endif // BOUNDINGBOX_H_INCLUDED

