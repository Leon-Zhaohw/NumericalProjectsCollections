#include "stdafx.h"
#include "sphere_tree.h"

#include "twigg/linalg.h"
#include "twigg/random.h"
#include "twigg/boundingbox.h"
#include "twigg/magicWrappers.h"

#include "Wm4DistVector3Line3.h"
#include "Wm4DistVector3Segment3.h"
#include "Wm4DistLine3Line3.h"
#include "Wm4DistSegment3Segment3.h"

#include <boost/spirit/core.hpp>
#include <boost/spirit/attribute.hpp>
#include <boost/spirit/symbols/symbols.hpp>
#include <boost/spirit/utility/chset.hpp>
#include <boost/spirit/utility/escape_char.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <boost/spirit/iterator.hpp>
#include <boost/spirit/error_handling/exceptions.hpp>

#include <algorithm>
#include <numeric>
#include <fstream>

namespace planning {

template <typename IntegerType>
IntegerType log2int( IntegerType v )
{
	BOOST_STATIC_ASSERT( !std::numeric_limits<IntegerType>::is_signed );
	BOOST_STATIC_ASSERT( std::numeric_limits<IntegerType>::is_integer );

	// from http://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
	IntegerType b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000};
	IntegerType S[] = {1, 2, 4, 8, 16};
	int i;

	register IntegerType r = 0; // result of log2(v) will go here
	for (i = 4; i >= 0; --i) // unroll for speed...
	{
		if (v & b[i])
		{
			v >>= S[i];
			r |= S[i];
		} 
	}

	return r;
}

float error( size_t start, size_t end, 
			const std::vector<vl::Vec3f>& points )
{
	float result = 0.0f;

	vl::Vec3f origin = points.at( start );
	vl::Vec3f direction = points.at( end ) - points.at( start );
	float dirLen = vl::len( direction );
	direction /= dirLen;

	Wm4::Segment3f segment;
	segment.Origin = Wm4::Vector3f( origin.Ref() );
	segment.Direction = Wm4::Vector3f( direction.Ref() );
	segment.Extent = dirLen;

	for( size_t iPoint = start; iPoint < end; ++iPoint )
	{
		Wm4::Vector3f point( points.at(iPoint).Ref() );
		Wm4::DistVector3Segment3f dist( point, segment );
		result = std::max( dist.Get(), result );
	}

	return result;
}

sphere_tree::sphere_tree( const std::vector<vl::Vec3f>& points, float epsilon )
{
	if( points.empty() )
		return;

	// need to trim points that are too close together:
	std::vector<vl::Vec3f> newPoints;
	newPoints.reserve( points.size() );

	{
		newPoints.push_back( points.front() );
		size_t start = 0;

		while( start < points.size() - 1 )
		{
			size_t low = start;
			size_t high = start + 1;

			while( high < points.size() )
			{
				if( error( start, high, points ) > epsilon )
					break;

				low = high;
				high += (high - start);
			}

			high = std::min( high, points.size() );
			while( high - low > 1 )
			{
				size_t probe = (high + low) / 2;
				if( error( start, probe, points ) > epsilon )
					high = probe;
				else
					low = probe;
			}

			size_t end = std::max(low, start+1);
            
			newPoints.push_back( points.at(end) );
			start = end;
		}
	}

	this->nodes_.reserve( newPoints.size() );

	// we will build it such that node i's children are at 2*i+1 and 2*i+2
	std::deque< std::vector<vl::Vec3f> > pointsQueue;
	pointsQueue.push_back( newPoints );

	bool noMoreChildren = false;

	while( !pointsQueue.empty() )
	{
		std::vector<vl::Vec3f> curPoints = pointsQueue.front();
		pointsQueue.pop_front();

		// we want everything to be perfectly packed, so we will need
		//   to make sure any 'dangling' nodes end up in the first set, not
		//   the second.
		size_t curSize = curPoints.size();
		size_t depth = 0;
		while( ((1 << depth) - 1) <= curSize )
			++depth;

		--depth;
		size_t leftover = curSize - ((1 << depth) - 1);
		size_t toLeft = std::min<size_t>( leftover, 1u << (depth-1) );
		size_t toRight = leftover - toLeft;
		toLeft += ((1 << (depth-1)) - 1);
		toRight += ((1 << (depth-1)) - 1);
		assert( (toRight + toLeft + 1) == curSize );

		std::vector<vl::Vec3f>::iterator mid = curPoints.begin() + toLeft;

		float radius = 0.0f;
		for( std::vector<vl::Vec3f>::const_iterator iter = curPoints.begin();
			iter != curPoints.end(); ++iter )
		{
			radius = std::max( radius, vl::len(*mid - *iter) );
		}

		vl::Vec3f adjacent = *mid;
		if( mid+1 != curPoints.end() )
			adjacent = *(mid+1);

		this->nodes_.push_back( node(*mid, adjacent, radius) );

		std::vector<vl::Vec3f> left( curPoints.begin(), mid );		size_t leftSize = left.size();
		std::vector<vl::Vec3f> right( mid+1, curPoints.end() );		size_t rightSize = right.size();
		assert( left.size() >= right.size() );

		if( !left.empty() )
		{
			assert( !noMoreChildren );
			pointsQueue.push_back( left );
		}
		else
			noMoreChildren = true;

		if( !right.empty() )
			pointsQueue.push_back( right );
	}
}

sphere_tree::sphere_tree( std::istream& is )
{
	readArray( this->nodes_, is );
}

sphere_tree::~sphere_tree()
{
}

void sphere_tree::write( std::ostream& out ) const
{
	dumpArray( this->nodes_, out );
}

std::vector<vl::Vec3f> sphere_tree::points() const
{
	std::vector<vl::Vec3f> result;
	this->getPoints( 0, result );
	return result;
}

void sphere_tree::getPoints( size_t nodeId, std::vector<vl::Vec3f>& vec ) const
{
	if( nodeId >= this->nodes_.size() )
		return;

	getPoints( 2*nodeId + 1, vec );
	vec.push_back( this->nodes_.at(nodeId).position );
	getPoints( 2*nodeId + 2, vec );
}

ProjectedInBoxQuery::ProjectedInBoxQuery( const vl::Mat4f& projectionMat, const BoundingBox2f& box )
	: projMat_(projectionMat)
{
	this->box_.expand( vl::Vec3f(box.minimum(), -1.0f) );
	this->box_.expand( vl::Vec3f(box.maximum(), 1.0f) );

	// first, need to transform the projection matrix:
	vl::Mat4f projMatXform( vl::vl_1 );
	for( size_t iCoord = 0; iCoord < 2; ++iCoord )
	{
		float x_0 = (box.minimum())[iCoord];
		float x_1 = (box.maximum())[iCoord];

		projMatXform[iCoord][iCoord] = 2.0f / (x_1 - x_0);
		projMatXform[iCoord][3] = 1.0f - 2.0f*x_1 / (x_1 - x_0);
	}

	vl::Mat4f newMat = projMatXform*projectionMat;

	projPlanes_ = computePlanes(newMat);

	/*
	this->box_.expand( vl::Vec3f(-1.0, -1.0f, 0.0f) );
	this->box_.expand( vl::Vec3f(1.0f, 1.0f, 1.0f) );
	this->projMat_ = newMat;
	*/
}

bool ProjectedInBoxQuery::operator()( const vl::Vec3f& point, float radius ) const
{
	return intersects(point, radius, projPlanes_);
}

bool ProjectedInBoxQuery::operator()( const vl::Vec3f& point, const vl::Vec3f& otherPoint ) const
{
	vl::Vec3f projected = vl::xform( projMat_, point );
	vl::Vec3f projectedParent = vl::xform( projMat_, otherPoint );

#ifdef _DEBUG
	assert( contains( box_, projected ) == (*this)(point, 0.0f) );
#endif

	if( contains( box_, projected ) )
		return true;

	// see if we can avoid having to run the expensive line-box test
	{
		const vl::Vec3f boxMin = box_.minimum();
		const vl::Vec3f boxMax = box_.maximum();
		for( vl::Int i = 0; i < 3; ++i )
		{
			if( (projected[i] > boxMax[i] && projectedParent[i] > boxMax[i]) ||
				(projected[i] < boxMin[i] && projectedParent[i] < boxMin[i]) )
				return false;
		}
	}

	vl::Vec3f diff = projectedParent - projected;

#ifdef _DEBUG
	// verify that the ray intersection code works even when ray is not
	// normalized
	{
		double len = vl::len( diff );
		if( len > 1e-4 )
		{
			BoundingBox3d box2( toVec3d(this->box_.minimum() ), 
				toVec3d(this->box_.maximum()) );
			Ray3d ray1( toVec3d(projected), vl::norm( toVec3d(diff) ) );
			double tMin1;
			double tMax1;
			bool intersected1 = box2.intersect( ray1, tMin1, tMax1, 0.0, len );

			Ray3f ray2( projected, diff );
			float tMin2;
			float tMax2;
			bool intersected2 = box_.intersect( ray2, tMin2, tMax2, 0.0f, 1.0f );
			assert( intersected1 == intersected2 );
		}
	}
#endif

	Ray3f ray( projected, diff );
	float tMin;
	float tMax;
	return box_.intersect( ray, tMin, tMax, 0.0f, 1.0f );
}


float sphere_tree::distance( const node& left, const node& right ) const
{
	Wm4::Segment3f leftSegment;
	leftSegment.Origin = toMagicVec( left.position );
	vl::Vec3f leftDir = left.adjacent - left.position;
	float leftDirLen = vl::len( leftDir );
	leftDir /= leftDirLen;
	leftSegment.Direction = toMagicVec( leftDir );
	leftSegment.Extent = leftDirLen;

	Wm4::Segment3f rightSegment;
	rightSegment.Origin = toMagicVec( right.position );
	vl::Vec3f rightDir = right.adjacent - right.position;
	float rightDirLen = vl::len( rightDir );
	rightDir /= rightDirLen;
	rightSegment.Direction = toMagicVec( rightDir );
	rightSegment.Extent = rightDirLen;

	Wm4::DistSegment3Segment3f dist( leftSegment, rightSegment );
	return dist.Get();
}

/*
float sphere_tree::maxDistance( const sphere_tree& right ) const
{
	const sphere_tree& left = *this;
	float maxDistance = 0.0f;

	float maxDistBound = boost::numeric::bounds<float>::max();
	float minDistBound = 0.0f;

	typedef std::pair<size_t, size_t> NodePr;
	std::deque<NodePr> toHandle( 1, NodePr(0, 0) );
	while( minDistBound < maxDistBound && !toHandle.empty() )
	{
		NodePr current = toHandle.front();
		toHandle.pop_front();

		float curDist = this->distance( 
			left.nodes_[ current.first ], right.nodes_[ current.second ] );

		float curMaxDist = curDist + std::min( 
			left.nodes_[ current.first ].radius, 
			right.nodes_[ current.second ].radius );
		float curMinDist = std::max( 0.0f, curDist - 
			left.nodes_[ current.first ].radius -
			right.nodes_[ current.second ].radius );

		if( curMaxDist
	}

	return maxDistance;
}
*/


bool sphere_tree::query( const vl::Mat4f& projectionMat, const BoundingBox2f& box ) const
{
	ProjectedInBoxQuery queryObject( projectionMat, box );

	return this->query( queryObject );
}

vl::Vec3f randomVec( RandomStream& random )
{
	return vl::Vec3f( random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0) );
}

void testSphereTree( const vl::Mat4f& projMatrix, const std::vector<vl::Vec3f>& points, RandomStream& random )
{
	for( size_t i = 0; i < 50; ++i )
	{
		std::vector<vl::Vec3f> points;
		for( size_t j = 0; j < i; ++j )
			points.push_back( randomVec(random) );

		sphere_tree tree( points, 0.0f );
	}

	sphere_tree tree( points, 0.0f );
	for( size_t iTest = 0; iTest < 10; ++iTest )
	{
		BoundingBox2f box;
		for( size_t i = 0; i < 2; ++i )
			box.expand( vl::Vec2f( random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0) ) );

		bool foundInSphereTree = tree.query( projMatrix, box );

		BoundingBox3f allBounds( vl::Vec3f(box.minimum(), -1.0), vl::Vec3f(box.maximum(), 1.0) );
		bool foundOne = false;
		for( size_t iPoint = 0; iPoint < points.size(); ++iPoint )
		{
			vl::Vec3f point = points[iPoint];
			vl::Vec3f xformed = vl::xform( projMatrix, point );
			if( contains(allBounds, xformed) )
			{
				foundOne = true;
				assert( foundInSphereTree );
			}
		}

		assert( foundOne == foundInSphereTree );
	}
}

void testSphereTree( const vl::Mat4f& projMatrix )
{
	RandomStream random;

	{
		boost::array<plane, 6> projPlanes = computePlanes(projMatrix);
		BoundingBox3f allBounds( vl::Vec3f(-1.0, -1.0, -1.0), vl::Vec3f(1.0, 1.0, 1.0) );
		for( size_t iTest = 0; iTest < 40; ++iTest )
		{
			vl::Vec3f vec;

			// test both small values and large ones
			if( iTest > 20 )
				vec = vl::Vec3f( random.uniform(-1000.0, 1000.0), 
					random.uniform(-1000.0, 1000.0), random.uniform(-1000.0, 1000.0) );
			else
				vec = randomVec(random);

			bool planeTest = contained( vec, projPlanes );

			vl::Vec3f projected = vl::xform( projMatrix, vec );
			bool projTest = contains( allBounds, projected );

			assert( projTest == planeTest );
		}
	}


	{
		BoundingBox2f box;
		for( size_t i = 0; i < 2; ++i )
			box.expand( vl::Vec2f( random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0) ) );

		vl::Mat4f projMatXform( vl::vl_1 );
		for( size_t iCoord = 0; iCoord < 2; ++iCoord )
		{
			float x_0 = (box.minimum())[iCoord];
			float x_1 = (box.maximum())[iCoord];

			projMatXform[iCoord][iCoord] = 2.0f / (x_1 - x_0);
			projMatXform[iCoord][3] = 1.0f - 2.0f*x_1 / (x_1 - x_0);
		}

		vl::Mat4f newProjMat = projMatXform*projMatrix;

		float epsilon = 1e-3;
		for( size_t iTest = 0; iTest < 5; ++iTest )
		{
			float z = random.uniform(-1000.0, 1000.0);
			vl::Vec3f pos1 = vl::xform( projMatrix, vl::Vec3f(box.minimum(), z) );
			vl::Vec3f pos2 = vl::xform( projMatXform, pos1 );

			vl::Vec3f lowPos = vl::xform( projMatXform, vl::Vec3f(box.minimum(), z) );
			assert( fabs(lowPos[0] + 1.0) < epsilon && fabs(lowPos[1] + 1.0) < epsilon );

			vl::Vec3f highPos = vl::xform( projMatXform, vl::Vec3f(box.maximum(), z) );
			assert( fabs(highPos[0] - 1.0) < epsilon && fabs(highPos[1] - 1.0) < epsilon );
		}

		boost::array<plane, 6> newPlanes = computePlanes(newProjMat);
		BoundingBox3f allBounds( vl::Vec3f(box.minimum(), -1.0), vl::Vec3f(box.maximum(), 1.0) );

		for( size_t iTest = 0; iTest < 40; ++iTest )
		{
			vl::Vec3f vec;

			// test both small values and large ones
			if( iTest > 20 )
				vec = vl::Vec3f( random.uniform(-1000.0, 1000.0), 
					random.uniform(-1000.0, 1000.0), random.uniform(-1000.0, 1000.0) );
			else
				vec = randomVec(random);

			bool planeTest = contained( vec, newPlanes );

			vl::Vec3f projected = vl::xform( projMatrix, vec );
			vl::Vec3f otherProj = vl::xform( newProjMat, vec );
			bool projTest = contains( allBounds, projected );

			assert( projTest == planeTest );
		}

	}

	for( size_t iNumNodes = 0; iNumNodes < 10; ++iNumNodes )
	{
		std::vector<vl::Vec3f> points;
		for( size_t j = 0; j < iNumNodes; ++j )
			points.push_back( randomVec(random) );

		testSphereTree( projMatrix, points, random );
	}

	for( size_t iTest = 0; iTest < 10; ++iTest )
	{
		size_t numNodes = random.uniform( 10, 1000 );
		std::vector<vl::Vec3f> points;
		for( size_t j = 0; j < numNodes; ++j )
			points.push_back( randomVec(random) );

		testSphereTree( projMatrix, points, random );
	}

	// try all the points on a line:
	vl::Vec3f line = vl::norm( randomVec(random) );
	vl::Vec3f center = randomVec(random);
	size_t numNodes = random.uniform( 10, 500 );
	std::vector<vl::Vec3f> points;
	for( size_t j = 0; j < numNodes; ++j )
		points.push_back( random.uniform(-1.0, 1.0) * line );
	testSphereTree( projMatrix, points, random );
}

// from this document:
// http://www2.ravensoft.com/users/ggribb/plane%20extraction.pdf
boost::array<plane, 6> computePlanes( vl::Mat4f matrix )
{
	matrix = trans(matrix);
	boost::array<plane, 6> planes;
	
	// Left clipping plane
	planes[0] = plane(
		matrix[0][3] + matrix[0][0],
		matrix[1][3] + matrix[1][0],
		matrix[2][3] + matrix[2][0],
		matrix[3][3] + matrix[3][0] );

	// Right clipping plane
	planes[1] = plane(
		matrix[0][3] - matrix[0][0],
		matrix[1][3] - matrix[1][0],
		matrix[2][3] - matrix[2][0],
		matrix[3][3] - matrix[3][0] );

	// Top clipping plane
	planes[2] = plane(
		matrix[0][3] - matrix[0][1],
		matrix[1][3] - matrix[1][1],
		matrix[2][3] - matrix[2][1],
		matrix[3][3] - matrix[3][1] );

	// Bottom clipping plane
	planes[3] = plane(
		matrix[0][3] + matrix[0][1], 
		matrix[1][3] + matrix[1][1], 
		matrix[2][3] + matrix[2][1], 
		matrix[3][3] + matrix[3][1] );

	// Near clipping plane
	planes[4] = plane(
		matrix[0][3] + matrix[0][2],
		matrix[1][3] + matrix[1][2],
		matrix[2][3] + matrix[2][2],
		matrix[3][3] + matrix[3][2] );

	// Far clipping plane
	planes[5] = plane(
		matrix[0][3] - matrix[0][2],
		matrix[1][3] - matrix[1][2],
		matrix[2][3] - matrix[2][2],
		matrix[3][3] - matrix[3][2] );

	return planes;
}

bool contained( const vl::Vec3f& vec, const boost::array<plane, 6>& planes )
{
	for( boost::array<plane, 6>::const_iterator iter = planes.begin();
		iter != planes.end(); ++iter )
	{
		float value = iter->evaluate(vec);
		if( value < 0.0f )
			return false;
	}

	return true;
}

bool intersects( const vl::Vec3f& point, float radius, const boost::array<plane, 6>& planes )
{
	for( boost::array<plane, 6>::const_iterator iter = planes.begin();
		iter != planes.end(); ++iter )
	{
		float value = iter->evaluate(point);
		if( value < -radius )
			return false;
	}

	return true;
}

using namespace boost::spirit;

SphereWrapperTree::SphereWrapperTree( const std::string& filename )
{
	std::ifstream ifs( filename.c_str() );
	if( !ifs )
		throw IOException( "Unable to open file '" + filename + "' for reading." );
	ifs.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

	std::vector<char> line( 1000 );
	ifs.getline( &line[0], line.size() );

	try
	{
		size_t lineNum = 1;

		size_t depth;
		{
			parse_info<char const*> info = parse( &line[0], 
				uint_p[assign(depth)] >> uint_p[assign(branching_)],
				space_p );
			if( !info.full )
				throw ParserException("Expected: [depth] [branching factor]", filename, lineNum, 
					info.stop - &line[0]);
		}

		size_t pow = 1;
		size_t totalSpheres = 0;
		for( size_t i = 0; i < depth; ++i )
		{
			totalSpheres += pow;
			pow *= branching_;
		}
		this->spheres_.reserve( totalSpheres );

		for( size_t iSphere = 0; iSphere < totalSpheres; ++iSphere )
		{
			++lineNum;

			if(!ifs)
			{
				std::ostringstream oss;
				oss << "Not enough spheres: only found " << this->spheres_.size() 
					<< "; expecting " << totalSpheres;
				throw ParserException(oss.str(), filename, lineNum, 1);
			}

			ifs.getline( &line[0], line.size() );
			vl::Vec3f center;
			float radius;
			parse_info<char const*> info = parse( &line[0], 
				real_p[assign(center[0])] >> real_p[assign(center[1])] >> real_p[assign(center[2])]
					>> real_p[assign(radius)] >> real_p, space_p );
			if( !info.full )
				throw ParserException("Expected: [x] [y] [z] [center] [weight]", filename, lineNum, 
					info.stop - &line[0]);

			/*
			if( !(radius > 0.0f) )
			{
				std::ostringstream errStr;
				throw ParserException( "Radius should be strictly positive.", filename, lineNum, 1 );
			}
			*/

			spheres_.push_back( Sphere(center, radius) );
		}
	}
	catch( std::ifstream::failure& e )
	{
		std::ostringstream oss;
		oss << "Error reading file '" << filename << "': " << e.what();
		throw IOException( oss.str() );
	}
}

std::vector<SphereWrapperTree::Sphere> SphereWrapperTree::level( size_t level ) const
{
	size_t current = 0;
	size_t currentSize = 1;
	for( size_t iLevel = 1; iLevel < level; ++iLevel )
	{
		current += currentSize;
		currentSize *= branching_;
	}

	std::vector<Sphere> result;
	result.reserve( currentSize );
	for( size_t i = current; i < current+currentSize; ++i )
	{
		if( !(spheres_[i].radius > 0.0) )
			continue;

		result.push_back( spheres_[i] );
	}

	return result;
}

} // namespace planning

