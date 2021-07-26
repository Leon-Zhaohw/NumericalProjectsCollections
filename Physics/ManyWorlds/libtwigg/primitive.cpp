#include "primitive.h"


std::pair< bool, TriangleIntersection > intersectTriangle( 
	const Ray& r, const Vec3& a, const Vec3& b, const Vec3& c )
{
	const double SMALL_EPSILON = 1e-10;
	const double epsilon = SMALL_EPSILON;
	const double NORMAL_EPSILON = 1e-8;

	Vec3 p = r.getPosition();
	Vec3 v = r.getDirection();

	Vec3 ab = b - a;
	Vec3 ac = c - a;
	Vec3 ap = p - a;

	Vec3 cv = cross(ab, ac);

	// there exists some bad triangles such that two vertices coincide
	// check this before normalize
	if( sqrlen(cv) < SMALL_EPSILON )
		return std::make_pair(false, TriangleIntersection());
	vl::Vec3 n = vl::norm(cv);

	Real vdotn = dot(v,n);
	if( fabs(vdotn) < NORMAL_EPSILON )
		return std::make_pair(false, TriangleIntersection());

	Real t = -dot(ap,n) / vdotn;

	if( t < epsilon )
		return std::make_pair(false, TriangleIntersection());

	// find k where k is the index of the component
	// of normal vector with greatest absolute value
	Real greatestMag = 0;
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

	Vec3 am = ap + t * v;

	Vec3 bary;
	bary[1] = cross(am, ac)[k]/cross(ab, ac)[k];
	bary[2] = cross(ab, am)[k]/cross(ab, ac)[k];
	bary[0] = 1-bary[1]-bary[2];

	if( bary[0] < -SMALL_EPSILON || 
		bary[1] < -SMALL_EPSILON || 
		bary[1] > (1+SMALL_EPSILON) || 
		bary[2] < -SMALL_EPSILON || 
		bary[2] > (1+SMALL_EPSILON) )
		return std::make_pair(false, TriangleIntersection());

	// if we get this far, we have an intersection.  Fill in the info.
	return std::make_pair(true,
		TriangleIntersection( t, bary ));
}




