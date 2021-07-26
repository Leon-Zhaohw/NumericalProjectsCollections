#ifndef __OBB_H__
#define __OBB_H__

#include "twigg/boundingbox.h"

namespace planning
{

class ProjectedInBoxQuery;

class OrientedBoundingBox
{
public:
	OrientedBoundingBox( const std::vector<vl::Vec3f>& points );
	OrientedBoundingBox( const vl::Vec3f& center, const vl::Vec3f& lengths,
		const vl::Mat3f& axes = vl::Mat3f(vl::vl_1) );
	bool query( const ProjectedInBoxQuery& planes ) const;
	float dist( const vl::Vec2f& point, const vl::Mat4f& transform ) const;
	bool intersect( const Ray3f& r ) const;

private:

	vl::Vec3f center_;
	vl::Mat3f axes_;
};

#define NumberSplitting 2

class OBBTree
{
	struct Node
	{
		Node( const OrientedBoundingBox& b )
			: box(b), numPoints(0)
		{
		}

		OrientedBoundingBox box;

		boost::array< boost::scoped_ptr<Node>, NumberSplitting > children;
		boost::scoped_array<vl::Vec3f> points;
		unsigned int numPoints;
	};

public:
	OBBTree( const std::vector<vl::Vec3f>& points );
	bool query( const ProjectedInBoxQuery& planes ) const; 
	void points( std::vector<vl::Vec3f>& v ) const;
	std::vector<vl::Vec3f> points() const;
	float dist( const vl::Vec2f& point, const vl::Mat4f& transform, float maxDist, float maxZ ) const;

private:
	bool query( const ProjectedInBoxQuery& planes, const Node* node ) const;
	float dist( const vl::Vec2f& point, const vl::Mat4f& transform, const Node* node, float minDist, float maxZ ) const;
	void points( std::vector<vl::Vec3f>& v, const Node* node ) const;

	Node* constructTree( const std::vector<vl::Vec3f>& points, Node* node = 0 );

	Node root_;
};

void testOBB();

float segmentDistance( const vl::Vec2f& p, const vl::Vec2f& segmentStart, const vl::Vec2f& segmentEnd );
double segmentDistance( const vl::Vec2d& p, const vl::Vec2d& segmentStart, const vl::Vec2d& segmentEnd );
float segmentDistance( const vl::Vec3f& p, const vl::Vec3f& segmentStart, const vl::Vec3f& segmentEnd );
double segmentDistance( const vl::Vec3d& p, const vl::Vec3d& segmentStart, const vl::Vec3d& segmentEnd );


} // namespace planning

#endif
