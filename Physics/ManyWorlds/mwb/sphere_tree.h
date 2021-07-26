#ifdef WIN32
#pragma once
#endif

#ifndef __SPHERETREE_H__
#define __SPHERETREE_H__

#include "twigg/boundingbox.h"
#include "twigg/random.h"

#include <boost/array.hpp>

namespace planning {

// eventually might want to template all this:
class sphere_tree
{
public:
	sphere_tree( const std::vector<vl::Vec3f>& points, float epsilon );
	sphere_tree( std::istream& is );
	~sphere_tree();

	void write( std::ostream& out ) const;

	bool query( const vl::Mat4f& projectionMat, const BoundingBox2f& box ) const;
//	float maxDistance( const sphere_tree& right ) const;  unimplemented

	std::vector<vl::Vec3f> points() const;

	template <typename QueryObject>
	bool query( QueryObject& queryObject ) const
	{
		std::deque<size_t> toHandle( 1, 0 );
		while( !toHandle.empty() )
		{
			size_t current = toHandle.front();

			toHandle.pop_front();
			if( current >= this->nodes_.size() )
				continue;

			// if queryObject tells us to, we will perform an early exit.
			const node& curNode = this->nodes_.at( current );
			if( queryObject(curNode.position, curNode.adjacent) )
				return true;

			if( queryObject(curNode.position, curNode.radius) )
			{
				toHandle.push_back( 2*current+1 );
				toHandle.push_back( 2*current+2 );
			}
		}

		// no node found
		return false;
	}

private:
	struct node
	{
		node( const vl::Vec3f& p, const vl::Vec3f& adj, float r )
			: position(p), adjacent(adj), radius(r) {}

		node() {}

		vl::Vec3f position;
		vl::Vec3f adjacent;
		float radius;
	};

	std::vector<node> nodes_;

	void construct_sphere_tree( std::vector<vl::Vec3f> points, float epsilon );
	void getPoints( size_t nodeId, std::vector<vl::Vec3f>& vec ) const;
	float distance( const node& left, const node& right ) const;
};

void testSphereTree( const vl::Mat4f& projMatrix );

class plane
{
public:
	plane() {}

	plane( float a, float b, float c, float d )
	{
		float magnitude = a*a + b*b + c*c;
		a_ = a/magnitude;
		b_ = b/magnitude;
		c_ = c/magnitude;
		d_ = d/magnitude;
	}

	float evaluate( const vl::Vec3f& point ) const
	{
		return (a_*point[0] + b_*point[1] + c_*point[2] + d_);
	}

	float dotWithNormal( const vl::Vec3f& vec ) const
	{
		return (a_*vec[0] + b_*vec[1] + c_*vec[2]);
	}

private:
	float a_;
	float b_;
	float c_;
	float d_;
};

bool contained( const vl::Vec3f& vec, const boost::array<plane, 6>& planes );
bool intersects( const vl::Vec3f& point, float radius, const boost::array<plane, 6>& planes );

boost::array<plane, 6> computePlanes( vl::Mat4f matrix );

class SphereWrapperTree
{
public:
	SphereWrapperTree( const std::string& filename );

	struct Sphere
	{
		Sphere( vl::Vec3f c, float r )
			: center(c), radius(r) {}

		vl::Vec3f center;

		// note: if radius <= 0.0, this sphere is not actually in the tree.
		float radius;
	};

	std::vector<Sphere> level( size_t level ) const;

	// single descent query, for objects that aren't sphere trees themselves.
	template <typename Handler>
	bool query( Handler& handler ) const
	{
		std::deque<size_t> active;
		if( handler.intersects( spheres_.front().center, spheres_.front().radius ) )
			active.push_back( 0 );

		while( !active.empty() )
		{
			size_t current = active.front();
			active.pop_front();

			size_t childrenStart = branching_*current + 1;
			if( childrenStart >= spheres_.size() )
			{
				handler.add( spheres_[current].center, spheres_[current].radius );
				continue;
			}

			std::vector<size_t> children;
			bool childrenBad = false;
			for( size_t iChild = childrenStart; iChild < childrenStart+branching_; ++iChild )
			{
				if( !spheres_[iChild].radius > 0.0 )
					continue;

				if( handler.inside( spheres_[iChild].center, spheres_[Child].radius ) )
				{
					childrenBad = true;
					break;
				}
			}

			if( !childrenBad && handler.hasRoom(children.size()) )
				std::copy( children.begin(), children.end(), std::back_inserter(active) );
			else
			{
				handler.add( spheres_[current].center, spheres_[current].radius );
				continue;
			}
		}
	}

private:

	size_t branching_;
	std::vector<Sphere> spheres_;
};

class ProjectedInBoxQuery
{
public:
	ProjectedInBoxQuery( const vl::Mat4f& projectionMat, const BoundingBox2f& box );
	bool operator()( const vl::Vec3f& point, float radius ) const;
	bool operator()( const vl::Vec3f& point, const vl::Vec3f& otherPoint ) const;
	const plane& getPlane( size_t iPlane ) const { return this->projPlanes_[iPlane]; }

private:
	vl::Mat4f projMat_;
	BoundingBox3f box_;

	// projPlanes_ are for modified projection matrix
	boost::array<plane, 6> projPlanes_;
};


float error( size_t start, size_t end, 
			const std::vector<vl::Vec3f>& points );

} // namespace planning

#endif // __SPHERETREE_H__
