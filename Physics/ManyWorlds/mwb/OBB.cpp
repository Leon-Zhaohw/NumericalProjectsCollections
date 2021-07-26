#include "stdafx.h"

#include "OBB.h"
#include "sphere_tree.h"

#ifdef NAG
#include <nag.h>
#include <nage04.h>
#include "twigg/nagWrappers.h"
#endif

#include "twigg/linalg.h"
#include "twigg/ioutil.h"
#include "twigg/random.h"

#include "Wm4Ray3.h"
#include "Wm4DistVector2Segment2.h"
#include "Wm4DistVector2Line2.h"
#include "Wm4DistLine3Line3.h"
#include "Wm4DistLine3Circle3.h"
#include "Wm4DistVector3Line3.h"
#include "Wm4DistVector3Segment3.h"
#include "Wm4DistSegment3Segment3.h"

#ifdef SLOW_DIST
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Cartesian.h>
#endif


#ifdef WIN32
#include <crtdbg.h>
#endif

#include <iostream>
#include <fstream>


namespace planning
{

#ifdef NAG
class NAGOptions
{
  public:
  NAGOptions()
  {
    nag_opt_init(&options_);
    options_.inf_bound = ldexp( 1.0, 20*3 );
    options_.list = FALSE;
#ifdef DEBUG_SOLVE
    options_.print_level = Nag_Soln;
#else
    options_.print_level = Nag_NoPrint;
#endif
  }

  Nag_E04_Opt* get()
  {
    return &options_;
  }

  ~NAGOptions()
  {
    NagError fail;
    INIT_FAIL(fail);

    nag_opt_free( &options_, "", &fail );

	assert( fail.code == NE_NOERROR );
  }

  private:
  Nag_E04_Opt options_;
};
#endif

template <typename Vec, typename Real>
Real segmentDistanceImpl( const Vec& p, const Vec& segmentStart, const Vec& segmentEnd )
{
	const Vec d( segmentEnd - segmentStart );
	const Vec pVec( p - segmentStart );

	const Real tOrig = vl::dot( pVec, d ) / sqrlen(d);
	const Real t = std::max<Real>( std::min<Real>( tOrig, (Real) 1 ), (Real) 0 );

	Real dist = vl::len( pVec - t*d );

	/*
	Wm4::Segment2f segment;
	segment.Origin = Wm4::Vector2f(segmentStart[0] + 0.5f*d[0], segmentStart[1] + 0.5f*d[1]);
	const float dLen = vl::len(d);
	segment.Direction = Wm4::Vector2f( d[0] / dLen, d[1] / dLen);
	segment.Extent = 0.5f*dLen;
	Wm4::DistVector2Segment2f wmlDistStruct( Wm4::Vector2f(p.Ref()), segment );
	float wmlDist = wmlDistStruct.Get();

	assert( fabs(wmlDist - dist) < (100.0f*std::max(wmlDist, 1.0f)) * std::numeric_limits<float>::epsilon() );
	*/

	return dist;
}

float segmentDistance( const vl::Vec2f& p, const vl::Vec2f& segmentStart, const vl::Vec2f& segmentEnd )
{
	return segmentDistanceImpl<vl::Vec2f, float>( p, segmentStart, segmentEnd );
}

double segmentDistance( const vl::Vec2d& p, const vl::Vec2d& segmentStart, const vl::Vec2d& segmentEnd )
{
	return segmentDistanceImpl<vl::Vec2d, double>( p, segmentStart, segmentEnd );
}

float segmentDistance( const vl::Vec3f& p, const vl::Vec3f& segmentStart, const vl::Vec3f& segmentEnd )
{
	return segmentDistanceImpl<vl::Vec3f, float>( p, segmentStart, segmentEnd );
}

double segmentDistance( const vl::Vec3d& p, const vl::Vec3d& segmentStart, const vl::Vec3d& segmentEnd )
{
	return segmentDistanceImpl<vl::Vec3d, double>( p, segmentStart, segmentEnd );
}

float segmentDistance( const vl::Vec2f& p, const vl::Vec3f& segmentStart, const vl::Vec3f& segmentEnd, float maxZ )
{
	// toss points outside the frustum:
	if( (segmentEnd[2] < -1.0f && segmentStart[2] < -1.0f) ||
        (segmentStart[2] > maxZ && segmentEnd[2] > maxZ) )
		return boost::numeric::bounds<float>::highest();

	return segmentDistance( p, strip(segmentStart), strip(segmentEnd) );
}

OrientedBoundingBox::OrientedBoundingBox( const std::vector<vl::Vec3f>& points )
{
	/*
	this->axes_ = vl::vl_1;
	BoundingBox3f box;
	for( std::vector<vl::Vec3f>::const_iterator pointItr = points.begin();
		pointItr != points.end(); ++pointItr )
	{
		box.expand( *pointItr );
	}
	this->center_ = 0.5*(box.maximum() + box.minimum());
	vl::Vec3f lengths = 0.5*(box.maximum() - box.minimum());
	for( vl::Int i = 0; i < 3; ++i )
		this->axes_[i] *= lengths[i];
		*/

	if( points.empty() )
	{
		this->axes_ = vl::vl_0;
		this->center_ = vl::vl_0;
		return;
	}

	// Okay, need to figure out axes now
	float epsilon = 1e-3;
	vl::Vec3f diff = points.back() - points.front();
	float diffLen = vl::len( diff );
	if( diffLen < epsilon )
		this->axes_[0] = vl::Vec3f( 0, 1, 0 );
	else
		this->axes_[0] = vl::norm(points.back() - points.front());

	this->axes_[1] = vl::vl_0;
	float bestAxis = 0.0f;
	for( std::vector<vl::Vec3f>::const_iterator pointItr = points.begin();
		pointItr != points.end(); ++pointItr )
	{
		vl::Vec3f vec = *pointItr - points.front();
		float proj = vl::dot( vec, this->axes_[0] );
		vl::Vec3f residual = vec - proj*this->axes_[0];

		if( vl::sqrlen( residual ) > bestAxis )
		{
			this->axes_[1] = residual;
			bestAxis = vl::sqrlen( residual );
		}
	}

	if( bestAxis < epsilon )
	{
		for( size_t i = 0; i < 3; ++i )
		{
			vl::Vec3f v( vl::vl_0 );
			v[i] = 1.0;
			vl::Vec3f xprod = vl::cross( this->axes_[0], v );
			if( vl::sqrlen(xprod) > bestAxis )
			{
				this->axes_[1] = xprod;
				bestAxis = vl::sqrlen(xprod);
			}
		}
	}
	else
	{
		this->axes_ = vl::vl_1;
	}

	// need a perpendicular basis
	this->axes_[2] = vl::norm( vl::cross( this->axes_[0], this->axes_[1] ) );
	this->axes_[1] = vl::norm( vl::cross( this->axes_[2], this->axes_[0] ) );

	BoundingBox3f bounds;
	for( std::vector<vl::Vec3f>::const_iterator pointItr = points.begin();
		pointItr != points.end(); ++pointItr )
	{
		bounds.expand( this->axes_ * *pointItr );
	}

	this->center_ = vl::trans(this->axes_) * (0.5*(bounds.minimum() + bounds.maximum()));
	vl::Vec3f lengths = 0.5*(bounds.maximum() - bounds.minimum());
	for( vl::Int i = 0; i < 3; ++i )
	{
		lengths[i] = std::max( lengths[i], epsilon );
		this->axes_[i] *= lengths[i];
	}
}

OrientedBoundingBox::OrientedBoundingBox( const vl::Vec3f& center, const vl::Vec3f& lengths,
	const vl::Mat3f& axes )
	: center_(center), axes_(axes)
{
	for( vl::Int i = 0; i < 3; ++i )
		this->axes_[i] *= lengths[i];
}

template <typename T>
T signum( const T& value )
{
	return (value < 0) ? -1 : 1;
}

bool OrientedBoundingBox::intersect( const Ray3f& ray ) const
{
	BoundingBox3f localBounds( 
		vl::Vec3f(-1.0f, -1.0f, -1.0f), 
		vl::Vec3f(1.0f, 1.0f, 1.0f) );
	vl::Mat3f invTransform = vl::inv( vl::trans(this->axes_) );

	//only for axis-aligned boxes
#ifdef _DEBUG
	bool myInt;
	{
		BoundingBox3f b( this->center_ - this->axes_*vl::Vec3f(1.0f, 1.0f, 1.0f),
			this->center_ + this->axes_*vl::Vec3f(1.0f, 1.0f, 1.0f) );
		float tMin, tMax;
		myInt = b.intersect( ray, tMin, tMax, 0.0f, boost::numeric::bounds<float>::highest() );
	}
#endif

	vl::Vec3f pos = invTransform * (ray.getPosition() - this->center_);
	vl::Vec3f dir = invTransform * (ray.getPosition() + ray.getDirection() - this->center_) - pos;
	Ray3f localRay( pos, dir );
	float tMin;
	float tMax;

	bool intersected = localBounds.intersect( localRay, tMin, tMax, 0.0f, boost::numeric::bounds<float>::highest() );
#ifdef _DEBUG
	assert( intersected == myInt );
#endif

	if( intersected )
	{
#ifdef _DEBUG
		// we will create a bigger box and make sure the ray misses that box too
		BoundingBox3f alignedBox;
		for( int i = 0; i < 2; ++i )
		{
			for( int j = 0; j < 2; ++j )
			{
				for( int k = 0; k < 2; ++k )
				{
					vl::Vec3f point = this->center_;
					point += (i == 0) ? (this->axes_[0]) : (-this->axes_[0]);
					point += (j == 0) ? (this->axes_[1]) : (-this->axes_[1]);
					point += (k == 0) ? (this->axes_[2]) : (-this->axes_[2]);

					alignedBox.expand( point );
				}
			}
		}

		float tMin2, tMax2;
		assert( alignedBox.intersect( ray, tMin2, tMax2, 0.0f, boost::numeric::bounds<float>::highest() ) );
#endif
		return true;
	}
	else
		return false;
}

#ifdef SLOW_DIST
typedef CGAL::Cartesian<float> K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon;

Point_2 toPoint_2( const vl::Vec3f& value )
{
	return Point_2( value[0], value[1] );
}

vl::Vec2f toVec2f( const Point_2& p )
{
	return vl::Vec2f( p.x(), p.y() );
}

float OrientedBoundingBox::dist( const vl::Vec2f& point, const vl::Mat4f& transform ) const
{

	// the first step is we compute the convex hull of the projection in the 2D plane
	boost::array<Point_2, 8> points = {{
		toPoint_2( vl::xform(transform, center_ - axes_[0] - axes_[1] - axes_[2]) ),
		toPoint_2( vl::xform(transform, center_ - axes_[0] - axes_[1] + axes_[2]) ),
		toPoint_2( vl::xform(transform, center_ - axes_[0] + axes_[1] - axes_[2]) ),
		toPoint_2( vl::xform(transform, center_ - axes_[0] + axes_[1] + axes_[2]) ),
		toPoint_2( vl::xform(transform, center_ + axes_[0] - axes_[1] - axes_[2]) ),
		toPoint_2( vl::xform(transform, center_ + axes_[0] - axes_[1] + axes_[2]) ),
		toPoint_2( vl::xform(transform, center_ + axes_[0] + axes_[1] - axes_[2]) ),
		toPoint_2( vl::xform(transform, center_ + axes_[0] + axes_[1] + axes_[2]) ) }};


	/*
	for( size_t i = 0; i < 2; ++i )
	{
		for( size_t j = 0; j < 2; ++j )
		{
			for( size_t k = 0; k < 2; ++k )
			{
				vl::Vec3f point;
                if( i == 0 )
					point += axes_[0];
				else
					point -= axes_[0];

				if( j == 0 )
					point += axes_[1];
				else
					point -= axes_[1]

				if( k == 0 )
					point += axes_[2];
				else
					point -= axes_[2];

				vl::Vec3f projPoint = vl::xform( transform, point );
				points[ i*1 + j*2 + k*4 ] = Point_2( projPoint[0], projPoint[1] );
			}
		}
	}
	*/
#ifdef _DEBUG
	float bruteForceMin = boost::numeric::bounds<float>::highest();
	for( size_t i = 0; i < 8; ++i )
	{
		for( size_t j = 1; j < 8; ++j )
		{
			bruteForceMin = std::min( bruteForceMin, 
				segmentDistance( point, toVec2f(points[i]), toVec2f(points[j]) ) );
		}
	}
#endif

	std::vector<Point_2> outHull;
	outHull.reserve(8);

	CGAL::ch_graham_andrew( points.begin(), points.end(), std::back_inserter(outHull) );
	CGAL::Bounded_side result = bounded_side_2( outHull.begin(), outHull.end(), 
		Point_2(point[0], point[1]), Polygon::Traits() );

	/*
#ifdef _DEBUG
	{
		fortran_matrix P( 8, 2 );
		for( size_t i = 0; i < 8; ++i )
		{
			P(i, 0) = points[i].x();
			P(i, 1) = points[i].y();
		}

		fortran_matrix H( outHull.size(), 2 );
		for( size_t i = 0; i < outHull.size(); ++i )
		{
			H(i, 0) = outHull[i].x();
			H(i, 1) = outHull[i].y();
		}

		fortran_matrix pointMat( 2, 1 );
		pointMat(0, 0) = point[0];
		pointMat(1, 0) = point[1];

		MATFile hull( "hull.mat", "Convex hull test" );
		hull.add( "P", P );
		hull.add( "H", H );
		hull.add( "point", pointMat );
	}
#endif
	*/

	if( result == CGAL::ON_BOUNDED_SIDE )
		return 0.0f;
	
	float minDist = boost::numeric::bounds<float>::highest();
	for( size_t i = 0; i < outHull.size(); ++i )
	{
		minDist = std::min( minDist, 
			segmentDistance( point, toVec2f(outHull[i]), toVec2f(outHull[(i+1)%outHull.size()]) ) );
	}

#ifdef _DEBUG
	assert( fabs(minDist - bruteForceMin) < 100.0f*std::numeric_limits<float>::epsilon() );
#endif

	return minDist;

	/*
	{
		vl::Vec3f pNear = vl::xform( this->invTransform_, vl::Vec3f(point, 0.0f) );
		vl::Vec3f pFar = vl::xform( this->invTransform_, vl::Vec3f(point, 1.0f) );
		Ray3f r( pNear, pFar - pNear );
		if( this->intersect(r) )
			return 0.0f;
	}

	// iterate through all sides of the box
	float myDist;
	{
		float minDist = boost::numeric::bounds<float>::highest();
		for( vl::Int iDimension = 0; iDimension < 3; ++iDimension )
		{
			vl::Int otherDimension1 = (iDimension + 1)%3;
			vl::Int otherDimension2 = (iDimension + 2)%3;

			vl::Vec3f dimCenter = this->center_ - this->axes_[iDimension];
			for( int j = 0; j < 4; ++j )
			{
				vl::Vec3f startPoint = dimCenter;
				if( j%2 )
					startPoint -= this->axes_[otherDimension1];
				else
					startPoint += this->axes_[otherDimension1];

				if( j/2 )
					startPoint -= this->axes_[otherDimension2];
				else
					startPoint += this->axes_[otherDimension2];

				vl::Vec3f endPoint = startPoint + 2.0f*this->axes_[iDimension];

				float curDist = segmentDistance( point, 
					vl::xform( transform, startPoint ),
					vl::xform( transform, endPoint ) );

				minDist = std::min( minDist, curDist );
			}
		}

		myDist = minDist;
	}
	*/
}
#endif // SLOW_DIST

bool OrientedBoundingBox::query( const ProjectedInBoxQuery& query ) const
{
	// test the n vector, as described in Optimized View Frustum Culling Algorithms for Bounding Boxes
	// by Ulf Assarsson and Tomas Moller
	for( size_t iPlane = 0; iPlane < 6; ++iPlane )
	{
		const plane& p = query.getPlane( iPlane );

		vl::Vec3f nVec( this->center_ );
		for( vl::Int i = 0; i < 3; ++i )
		{
			/*
			nVec += (lengths_[i] * signum( p.dotWithNormal( axes_[i] ) )) * axes_[i];
			*/

			if( p.dotWithNormal( axes_[i] ) > 0.0f )
				nVec += axes_[i];
			else
				nVec -= axes_[i];
		}
        
		if( p.evaluate( nVec ) < 0.0f )
			return false;
		/*
		bool outside = true;
		for( size_t iPoint = 0; iPoint < 8; ++iPoint )
		{
			vl::Vec3f point( this->center_ );
			for( size_t iAxis = 0; iAxis < 3; ++iAxis )
				point += (((iPoint >> iAxis) & 1) ? (lengths_[iAxis]) : (-lengths_[iAxis])) * axes_[iAxis];
			outside = outside && (p.evaluate( point ) < 0.0f);
			if( !outside )
				break;
		}

		if( outside )
			return false;
			*/
	}

	return true;
}

OBBTree::OBBTree( const std::vector<vl::Vec3f>& points )
	: root_( OrientedBoundingBox(points) )
{
	this->constructTree( points, &root_ );
}


OBBTree::Node* OBBTree::constructTree( const std::vector<vl::Vec3f>& points, Node* newNode )
{
	const size_t maxPoints = 30;
	if( newNode == 0 )
		newNode = new Node( OrientedBoundingBox(points) );

	if( points.size() > maxPoints )
	{
		for( size_t i = 0; i < NumberSplitting; ++i )
		{
			size_t low = i * points.size() / NumberSplitting;
			size_t high = (i+1) * points.size() / NumberSplitting;
			assert( low >= 0 && low < points.size() );
			assert( high >= 0 && high <= points.size() );

			std::vector<vl::Vec3f>::const_iterator highItr = 
				points.begin() + high;
			if( highItr != points.end() )
				++highItr;

			std::vector<vl::Vec3f> childrenPoints( points.begin() + low, highItr );
			newNode->children[i].reset( constructTree( childrenPoints ) );
		}
	}
	else if( !points.empty() )
	{
		newNode->points.reset( new vl::Vec3f[points.size()] );
		std::copy( points.begin(), points.end(), 
			&newNode->points[0] );
		newNode->numPoints = points.size();
	}

	return newNode;
}

bool OBBTree::query( const ProjectedInBoxQuery& q ) const
{
	return this->query( q, &this->root_ );
}

bool OBBTree::query( const ProjectedInBoxQuery& q, const Node* node ) const
{
	if( !node->box.query(q) )
		return false;

	if( node->points )
	{
		for( size_t iPoint = 0; (iPoint+1) < node->numPoints; ++iPoint )
		{
			if( q( node->points[iPoint], node->points[iPoint+1] ) )
				return true;
		}

		return false;
	}

	for( size_t iChild = 0; iChild < node->children.size(); ++iChild )
	{
		if( !node->children[iChild] )
			break;

		if( query(q, node->children[iChild].get()) )
			return true;
	}

	return false;
}

std::vector<vl::Vec3f> OBBTree::points() const
{
	std::vector<vl::Vec3f> result;
	this->points( result );
	return result;
}

void OBBTree::points( std::vector<vl::Vec3f>& v ) const
{
	this->points( v, &this->root_ );
}

void OBBTree::points( std::vector<vl::Vec3f>& v, const Node* node ) const
{
	if( node->points )
	{
		std::copy( &node->points[0], &node->points[0] + node->numPoints, 
			std::back_inserter(v) );
		return;
	}

	for( size_t iChild = 0; iChild < node->children.size(); ++iChild )
	{
		if( !node->children[iChild] )
			break;

		this->points( v, node->children[iChild].get() );

		// remove repeated index
		if( (iChild+1) != node->children.size() )
			v.pop_back();
	}
}

float OBBTree::dist( const vl::Vec2f& point, const vl::Mat4f& transform, float maxDist, float maxZ ) const
{
	return this->dist( point, transform, &this->root_, maxDist, maxZ );
}

float OBBTree::dist( const vl::Vec2f& point, const vl::Mat4f& transform, const Node* node, float minDist, float maxZ ) const
{
	/*
	float boxDist = node->box.dist( point, transform );
	if( boxDist > minDist )
	{
#ifdef _DEBUG
		for( size_t iChild = 0; iChild < node->children.size(); ++iChild )
		{
			if( node->children[iChild] )
			{
				float childDist = this->dist(point, transform, node->children[iChild], boost::numeric::bounds<float>::highest(), boxDist);
				assert( (childDist + 100.0f*std::numeric_limits<float>::epsilon()) >= boxDist );
			}
		}
#endif
		return boxDist;
	}
	*/

	float myMinDist = boost::numeric::bounds<float>::highest();

	if( node->points )
	{
		vl::Vec3f prevPoint = vl::xform(transform, node->points[0]);
		const size_t numPoints = node->numPoints;
		for( size_t iPoint = 1; iPoint < numPoints; ++iPoint )
		{
			vl::Vec3f nextPoint = vl::xform(transform, node->points[iPoint]);
			if( vl::len( point - strip(prevPoint) ) - vl::len( strip(prevPoint) - strip(nextPoint) ) > minDist )
			{
				prevPoint = nextPoint;
				continue;
			}

			float dist = planning::segmentDistance( point, prevPoint, nextPoint, maxZ );
			myMinDist = std::min( myMinDist, dist );

			prevPoint = nextPoint;
		}
	}

	for( size_t iChild = 0; iChild < node->children.size(); ++iChild )
	{
		if( !node->children[iChild] )
			break;

		myMinDist = std::min( 
			this->dist(point, transform, node->children[iChild].get(), 
				std::min(minDist, myMinDist), maxZ ), myMinDist );
	}	

	return myMinDist;
}

bool bruteForceQuery( const std::vector<vl::Vec3f>& path, const ProjectedInBoxQuery& query )
{
	const size_t tmp = path.size() - 1;
	for( size_t j = 0; j < tmp; ++j )
	{
		if( query( path.at(j), path.at(j+1) ) )
			return true;
	}

	return false;
}

// Borrowed from the ray tracer code.
/*
class TestRay
{
public:
	TestRay( const vl::Vec3f& pp, const vl::Vec3f& dd )
		: position_( pp ), direction_( dd )
	{
		for( vl::Int i = 0; i < 3; ++i )
			inverseDirection_[i] = 1. / direction_[i];

		for( vl::Int i = 0; i < 3; ++i )
		{
			if( inverseDirection_[i] < 0.0f )
				sign_[i] = 1;
			else
				sign_[i] = 0;
		}
	}

private:
	vl::Vec3f position_;
	vl::Vec3f direction_;
	vl::Vec3f inverseDirection_;
	boost::array<char, 3> sign_;
};
*/


struct PickQuery
{
	vl::Mat4f projectionMatrix;
	vl::Mat4f modelViewMatrix;
	vl::Vec2f selectedPoint;
};

float bruteForcePickQuery( const std::vector<vl::Vec3f>& path, const PickQuery& q )
{
	float minDist = boost::numeric::bounds<float>::highest();
	if( path.empty() )
		return minDist;

	vl::Mat4f mat = q.projectionMatrix*q.modelViewMatrix;
	vl::Vec3f prevPoint = vl::xform( mat, path.front() );

	const size_t pathSize = path.size();
	for( size_t i = 1; i < pathSize; ++i )
	{
		vl::Vec3f nextPoint = vl::xform( mat, path[i] );

		float dist = planning::segmentDistance( q.selectedPoint, 
			prevPoint,
			nextPoint, 
			1.0f );
		minDist = std::min( dist, minDist );
		prevPoint = nextPoint;
	}

	return minDist;
}


void testOBB()
{
	std::deque< std::vector<vl::Vec3f> > paths;

	{
        std::ifstream pathStream( "//volley.graphics.cs.cmu.edu/cdtwigg/planning/testPaths.bin", std::ios::binary );
		assert( pathStream );
		pathStream.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
		std::ifstream::pos_type m = filelen( pathStream );

		while( pathStream.tellg() < m )
		{
			paths.push_back( std::vector<vl::Vec3f>() );
			readArray( paths.back(), pathStream );
		}
	}

	typedef std::pair< vl::Mat4f, BoundingBox2f > ConstraintPr;
	std::deque< ConstraintPr > constraints;

	{
        std::ifstream constraintStream( "//volley.graphics.cs.cmu.edu/cdtwigg/planning/constraints.bin", std::ios::binary );
		assert( constraintStream );
		constraintStream.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
		std::ifstream::pos_type m = filelen( constraintStream );

		while( constraintStream.tellg() < m )
		{
			vl::Mat4f mat;
			constraintStream.read( reinterpret_cast<char*>( &mat ), sizeof( vl::Mat4f ) );
			assert( !isBad(mat) );

			BoundingBox2f box;
			constraintStream.read( reinterpret_cast<char*>( &box ), sizeof( BoundingBox2f ) );
			assert( !isBad(box.minimum()) );
			assert( !isBad(box.maximum()) );

			constraints.push_back( ConstraintPr(mat, box) );
		}
	}

	std::deque<PickQuery> pickQueries;
	{
		std::ifstream pickStream( "//volley.graphics.cs.cmu.edu/cdtwigg/planning/clicks.bin", std::ios::binary );
		assert( pickStream );
		pickStream.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
		std::ifstream::pos_type m = filelen( pickStream );

		while( pickStream.tellg() < m )
		{
			PickQuery q;
			pickStream.read( reinterpret_cast<char*>( q.modelViewMatrix.Ref() ), sizeof( vl::Mat4f ) );
			pickStream.read( reinterpret_cast<char*>( q.projectionMatrix.Ref() ), sizeof( vl::Mat4f ) );
			pickStream.read( reinterpret_cast<char*>( q.selectedPoint.Ref() ), sizeof( vl::Vec2f ) );
			pickQueries.push_back( q );
		}
	}

	boost::ptr_deque< OBBTree > obbTrees;
	for( size_t iPath = 0; iPath < paths.size(); ++iPath )
		obbTrees.push_back( new OBBTree( paths.at(iPath) ) );

	boost::ptr_deque< sphere_tree > sphereTrees;
	for( size_t iPath = 0; iPath < paths.size(); ++iPath )
		sphereTrees.push_back( new sphere_tree( paths.at(iPath), 0.0f ) );

	std::deque<ProjectedInBoxQuery> queries;

	for( size_t iConstraint = 0; iConstraint < constraints.size(); ++iConstraint )
	{
		ConstraintPr con = constraints[iConstraint];
		queries.push_back( ProjectedInBoxQuery(con.first, con.second) );
		const ProjectedInBoxQuery& query = queries.back();

		for( size_t iPath = 0; iPath < paths.size(); ++iPath )
		{
			const std::vector<vl::Vec3f>& path = paths.at(iPath);
			bool bruteForceResult = bruteForceQuery( path, query );

			bool obbResult = obbTrees.at(iPath).query( query );
			assert( obbResult == bruteForceResult );

			/*
			bool sphereTreeResult = sphereTrees.at(iPath).query( con.first, con.second );
			assert( sphereTreeResult == bruteForceResult );
			*/
		}
	}

	// now test pick queries
	{
		for( size_t iQuery = 0; iQuery < pickQueries.size(); ++iQuery )
		{
			float minDist = boost::numeric::bounds<float>::highest();
			for( size_t iPath = 0; iPath < paths.size(); ++iPath )
			{
				float bruteForceDist = bruteForcePickQuery( paths.at(iPath), pickQueries.at(iQuery) );
				float obbTreeDist = obbTrees.at(iPath).dist( pickQueries[iQuery].selectedPoint,
					pickQueries[iQuery].projectionMatrix * pickQueries[iQuery].modelViewMatrix, minDist, 1.0 );
				assert( (bruteForceDist > minDist) || fabs(bruteForceDist - obbTreeDist) < (std::numeric_limits<float>::epsilon()*10.0f) );
				minDist = std::min( minDist, bruteForceDist );
			}
		}
	}

#ifndef _DEBUG // performance testing doesn't make sense in debug mode
	const size_t numTest = 20;

	double obbPickTime = 0, bruteForcePickTime = 0;
	{
		time_t start = clock();
		for( size_t iTest = 0; iTest < numTest; ++iTest )
		// now test pick queries
		{
			for( size_t iQuery = 0; iQuery < pickQueries.size(); ++iQuery )
			{
				float minDist = boost::numeric::bounds<float>::highest();
				for( size_t iPath = 0; iPath < paths.size(); ++iPath )
				{
					float obbTreeDist = obbTrees.at(iPath).dist( pickQueries[iQuery].selectedPoint,
						pickQueries[iQuery].projectionMatrix * pickQueries[iQuery].modelViewMatrix, minDist, 1.0f );
					minDist = std::min( minDist, obbTreeDist );
				}
			}
		}
		time_t end = clock();
		obbPickTime = boost::numeric_cast<double>(end - start) / boost::numeric_cast<double>( CLOCKS_PER_SEC );
	}

	{
		time_t start = clock();
		for( size_t iTest = 0; iTest < numTest; ++iTest )
		// now test pick queries
		{
			for( size_t iPath = 0; iPath < paths.size(); ++iPath )
			{
				for( size_t iQuery = 0; iQuery < pickQueries.size(); ++iQuery )
				{
					float bruteForceDist = bruteForcePickQuery( paths.at(iPath), pickQueries.at(iQuery) );
				}
			}
		}
		time_t end = clock();
		bruteForcePickTime = boost::numeric_cast<double>(end - start) / boost::numeric_cast<double>( CLOCKS_PER_SEC );
	}

	/*
	double obbTime = 0, sphereTreeTime = 0, bruteForceTime = 0;


	{
		time_t start = clock();
		for( size_t iTest = 0; iTest < numTest; ++iTest )
		{
			for( size_t iConstraint = 0; iConstraint < constraints.size(); ++iConstraint )
			{
				const ProjectedInBoxQuery& query = queries.at(iConstraint);

				for( size_t iPath = 0; iPath < paths.size(); ++iPath )
				{
					bool res = obbTrees[iPath].query( query );
				}
			}
		}
		time_t end = clock();
		obbTime = boost::numeric_cast<double>(end - start) / boost::numeric_cast<double>( CLOCKS_PER_SEC );
	}

	{
		time_t start = clock();
		for( size_t iTest = 0; iTest < numTest; ++iTest )
		{
			for( size_t iConstraint = 0; iConstraint < constraints.size(); ++iConstraint )
			{
				const ConstraintPr& constraint = constraints.at(iConstraint);
				const ProjectedInBoxQuery& query = queries.at(iConstraint);

				for( size_t iPath = 0; iPath < paths.size(); ++iPath )
				{
					bool res = sphereTrees[iPath].query( constraint.first, constraint.second );
				}
			}
		}
		time_t end = clock();
		sphereTreeTime = boost::numeric_cast<double>(end - start) / boost::numeric_cast<double>( CLOCKS_PER_SEC );
	}

	{
		time_t start = clock();
		for( size_t iTest = 0; iTest < numTest; ++iTest )
		{
			for( size_t iConstraint = 0; iConstraint < constraints.size(); ++iConstraint )
			{
				const ProjectedInBoxQuery& query = queries.at(iConstraint);

				for( size_t iPath = 0; iPath < paths.size(); ++iPath )
				{
					bool bruteForceResult = bruteForceQuery( paths.at(iPath), query );
				}
			}

		}
		time_t end = clock();
		bruteForceTime = boost::numeric_cast<double>(end - start) / boost::numeric_cast<double>( CLOCKS_PER_SEC );
	}
	*/

	int numQueries = numTest * constraints.size() * paths.size();
	std::ofstream ofs( "results.txt" );
	ofs << "num queries: " << numQueries << "\n";
	ofs << "num pick queries: " << (numTest * paths.size() * pickQueries.size()) << "\n";
	ofs << "obb tree picking: " << obbPickTime << "\n";
	ofs << "brute force picking: " << bruteForcePickTime << "\n";
//	ofs << "obb tree: " << obbTime << "\n";
//	ofs << "sphere tree: " << sphereTreeTime << "\n";
//	ofs << "brute force: " << bruteForceTime << "\n";
//	ofs << "ray build: " << rayTime << "\n";
	std::cout << "test.";
#endif
}

} // namespace planning
