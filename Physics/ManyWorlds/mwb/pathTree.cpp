#include "stdafx.h"

#include "pathTree.h"
#include "simulationTree.h"
#include "sphere_tree.h"

#include "twigg/magicWrappers.h"

#include "Wm4DistVector3Segment3.h"
#include "Wm4DistVector3Line3.h"
#include "Wm4DistLine3Line3.h"

namespace planning {

PathTree::PathTree( size_t objectId )
	: object_( objectId ), root_(0)
{
}

PathTree::~PathTree()
{
}

float asymmetricDifference( const std::vector<vl::Vec3f>& leftPoints,
						  const std::vector<vl::Vec3f>& rightPoints )
{
	float maxDist = 0.0f;

	const int dist = 10;
	// heuristic technique
	// we march along the left path, at each point looking up to dist steps in each
	// direction on the right path to find the closest point of approach
	int rightCur = 0;

	for( int iLeft = 0; iLeft < leftPoints.size(); ++iLeft )
	{
		float minDist = boost::numeric::bounds<float>::highest();

		Wm4::Vector3f point = toMagicVec( leftPoints[iLeft] );
		for( int iRight = std::max<int>( 0, rightCur - dist);
			iRight < std::min<int>( rightCur + dist, rightPoints.size()-2 );
			++iRight )
		{
			Wm4::Segment3f rightSegment;
			rightSegment.Origin = toMagicVec( rightPoints[iRight] );
			vl::Vec3f dir = rightPoints[iRight+1] - rightPoints[iRight];
			float dirLen = vl::len(dir);
			rightSegment.Direction = toMagicVec( dir / dirLen );
			rightSegment.Extent = dirLen;

			Wm4::DistVector3Segment3f wmlDist( point, rightSegment );
			float curDist = wmlDist.Get();
			if( curDist < minDist )
			{
				rightCur = iRight;
				minDist = curDist;
			}
		}

		maxDist = std::max( maxDist, minDist );
	}

	return maxDist;
}

float PathTree::distance( const Path* left, const Path* right ) const
{
	const std::vector<vl::Vec3f> leftPoints = left->obbTree( this->object_ ).points();
	const std::vector<vl::Vec3f> rightPoints = right->obbTree( this->object_ ).points();

	return std::max( asymmetricDifference(leftPoints, rightPoints),
		asymmetricDifference(rightPoints, leftPoints) );

}

void PathTree::add( const Path* path )
{
	if( this->root_ == 0 )
	{
		this->root_ = this->nodePool_.construct(path);
		return;
	}

	Node* currentNode = root_;
	std::pair<float, const Path*> currentPath( distance( root_->path, path ), path );

	while( true )
	{
		currentNode->radius = std::max( currentNode->radius, currentPath.first );

		// first try to insert into the current node
		for( Node::OtherPathArray::size_type i = 0; i < currentNode->otherPaths.size(); ++i )
		{
			if( currentPath.first < currentNode->otherPaths[i].first )
				std::swap( currentPath, currentNode->otherPaths[i] );
		}

		if( currentPath.second == 0 )
			return;

		if( currentNode->leftChild == 0 )
		{
			currentNode->leftChild = this->nodePool_.construct(path);
			return;
		}

		if( currentNode->rightChild == 0 )
		{
			currentNode->rightChild = this->nodePool_.construct(path);
			return;
		}

		float leftDist = distance( currentNode->leftChild->path, path );
		float rightDist = distance( currentNode->rightChild->path, path );
		if( leftDist < rightDist )
		{
			currentNode = currentNode->leftChild;
			currentPath.first = leftDist;
		}
		else
		{
			currentNode = currentNode->rightChild;
			currentPath.first = rightDist;
		}
	}
}

std::deque<const Path*> PathTree::query( const vl::Mat4f& matrix, const BoundingBox2f& box ) const
{
	ProjectedInBoxQuery query( matrix, box );

	// first, need to transform the projection matrix:
	vl::Mat4f projMatXform( vl::vl_1 );
	for( size_t iCoord = 0; iCoord < 2; ++iCoord )
	{
		float x_0 = (box.minimum())[iCoord];
		float x_1 = (box.maximum())[iCoord];

		projMatXform[iCoord][iCoord] = 2.0f / (x_1 - x_0);
		projMatXform[iCoord][3] = 1.0f - 2.0f*x_1 / (x_1 - x_0);
	}

	vl::Mat4f newMat = projMatXform*matrix;
	boost::array<plane, 6> planes = computePlanes( newMat );

	std::deque<const Path*> result;

	std::deque<const Node*> toHandle(1, root_);
	while( !toHandle.empty() )
	{
		const Node* current = toHandle.front();
		toHandle.pop_front();

		// @todo this is really going to be painfully slow I think
		const std::vector<vl::Vec3f> points = current->path->obbTree( this->object_ ).points();

		float furthestInside = 0.0f;
		for( size_t iPoint = 0; iPoint < points.size(); ++iPoint )
		{
			float minValue = boost::numeric::bounds<float>::highest();
			for( size_t i = 0; i < planes.size(); ++i )
			{
				minValue = std::min( minValue,
					planes[i].evaluate( points[iPoint] ) );
			}

			furthestInside = std::max( minValue, furthestInside );
		}

		if( furthestInside > 0.0f )
		{
			if( current->radius < furthestInside )
			{
				std::deque< const Node* > dfsQueue(1, current);
				while( !dfsQueue.empty() )
				{
					const Node* cur = dfsQueue.front();
					dfsQueue.pop_front();

					result.push_back( cur->path );
					for( size_t i = 0; (i < cur->otherPaths.size())
						&& (cur->otherPaths[i].second != 0); ++i )
					{
						result.push_back( cur->otherPaths[i].second );
					}

					if( cur->leftChild )
						dfsQueue.push_front( cur->leftChild );

					if( cur->rightChild )
						dfsQueue.push_front( cur->rightChild );
				}
			}
			else
			{
				result.push_back( current->path );

				// passes through query rectangle
				for( size_t i = 0; (i < current->otherPaths.size())
					&& (current->otherPaths[i].second != 0); ++i )
				{
					if( current->otherPaths[i].first < furthestInside )
						result.push_back( current->otherPaths[i].second );
					else if( current->otherPaths[i].second->obbTree( this->object_ ).query( query ) )
						result.push_back( current->otherPaths[i].second );
				}

				if( current->leftChild )
					toHandle.push_front( current->leftChild );

				if( current->rightChild )
					toHandle.push_front( current->rightChild );			
			}
		}
		else
		{
			// doesn't pass through query rectangle
			// do something here.
		}
	}

	return result;
}

} // namespace planning
