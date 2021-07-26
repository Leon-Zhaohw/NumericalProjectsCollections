#include "stdafx.h"

#include "twigg/primitive.h"
#include "twigg/boundingbox.h"
#include <cfloat>
#include <boost/numeric/conversion/bounds.hpp>

using namespace std;
using namespace vl;

bool contains(const BoundingBox2d& bb2f, const vl::Vec2& v2fPt, const Real fEpsilon/* = 0.0*/)
{
	return (
		v2fPt[0] >= bb2f.minimum()[0] - fEpsilon &&
		v2fPt[1] >= bb2f.minimum()[1] - fEpsilon &&
		v2fPt[0] <= bb2f.maximum()[0] + fEpsilon &&
		v2fPt[1] <= bb2f.maximum()[1] + fEpsilon);
}

bool intersects(const BoundingBox2d& bb1, const BoundingBox2d& bb2, const Real fEpsilon/* = 0.0*/)
{
	if (bb1.maximum()[0] + fEpsilon < bb2.minimum()[0] ||
		bb1.maximum()[1] + fEpsilon < bb2.minimum()[1])
		return false;
	if (bb1.minimum()[0] > bb2.maximum()[0] + fEpsilon ||
		bb1.minimum()[1] > bb2.maximum()[1] + fEpsilon)
		return false;

	return true;
}


BoundingBox3d transformBounds( const BoundingBox3d& box, const vl::Mat4d& transform )
{
	BoundingBox3d newBox;

	const vl::Vec3d min = box.minimum();
	const vl::Vec3d max = box.maximum();
	for( int i = 0; i < 2; ++i )
		for( int j = 0; j < 2; ++j )
			for( int k = 0; k < 2; ++k )
			{
				vl::Vec3d point( 
					i == 0 ? min[0] : max[0],
					j == 0 ? min[1] : max[1],
					k == 0 ? min[2] : max[2] );
				newBox.expand( xform(transform, point) );
			}

	return newBox;
}

