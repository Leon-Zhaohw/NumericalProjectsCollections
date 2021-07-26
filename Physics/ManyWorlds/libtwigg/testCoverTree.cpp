// testCoverTree.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "../planning/cover_tree.h"

#include "twigg/random.h"

#include <boost/array.hpp>

const size_t N = 5;
typedef boost::array<float, N> PointType;

struct DistanceType
{
	PointType::value_type operator()( const PointType& left, const PointType& right ) const
	{
		PointType::value_type result = 0.0;
		for( size_t i = 0; i < N; ++i )
		{
			PointType::value_type diff = left[i] - right[i];
			result += diff*diff;
		}

		return std::sqrt( result );
	}
};


class DistanceComparator
{
public:
	DistanceComparator( const PointType& point )
		: point_( point ) {}

	bool operator()( const PointType& lhs, const PointType& rhs ) const
	{
		const DistanceType d;
		return std::less<float>()( d(lhs, point_), d(rhs, point_) );
	}

private:
	PointType point_;
};

class PointTypeComparator
{
public:
	bool operator()( const PointType& left, const PointType& right )
	{
		for( size_t i = 0; i < N; ++i )
		{
			if( left[i] < right[i] )
				return true;

			if( left[i] > right[i] )
				return false;
		}

		return false;
	}
};

PointType randomPoint()
{
	static RandomStream random;
	PointType p;
	random.normal( 0.0, 1.0, N, &p[0] );
	return p;
}

typedef cover_tree<PointType, DistanceType> CoverTree;

void test_knn( size_t k, const CoverTree& tree, std::deque<PointType>& allPoints )
{
	assert( k > 0 );

	PointType p = randomPoint();

	std::vector<PointType> knn = tree.knn( p, k );

	std::deque<PointType>::iterator kth = allPoints.begin() + std::min(allPoints.size(), k) - 1;
	std::nth_element( allPoints.begin(), kth, allPoints.end(), DistanceComparator(p) );
	std::vector<PointType> knn_naive( allPoints.begin(), kth + 1 );

	std::sort( knn.begin(), knn.end(), PointTypeComparator() );
	std::sort( knn_naive.begin(), knn_naive.end(), PointTypeComparator() );
	assert( knn.size() == knn_naive.size() );
	for( size_t i = 0; i < std::min(k, allPoints.size()); ++i )
	{
		const DistanceType d;
		PointType left = knn[i];
		PointType right = knn_naive[i];
		float dLeft = d( left, p );
		float dRight = d( right, p );

		assert( left == right );
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	RandomStream random;

	{
		// test the empty tree case
		CoverTree tree;
		std::vector<PointType> points = tree.knn( randomPoint() );
		assert( points.empty() );
	}

	for( size_t j = 1; j < 20; ++j )
	{
		CoverTree tree;
		std::deque<PointType> allPoints;

		// first, build the tree and make sure there are
		//   no problems
		size_t num = random.uniform( 10, 110 );
		if( j < 10 )
			num = j;

		for( size_t i = 0; i < num; ++i )
		{
			assert( tree.size() == i );
			PointType p = randomPoint();
			allPoints.push_back( p );
			tree.insert( p );
			tree.check_invariants();
		}

		CoverTree batchTree( allPoints );
		assert( batchTree.size() == tree.size() );
		
		// now, run a bunch of k nearest neighbor searches
		for( size_t i = 0; i < 20; ++i )
		{
			test_knn( 1, tree, allPoints );
			test_knn( 2, tree, allPoints );
			test_knn( tree.size(), tree, allPoints );
			test_knn( 2*tree.size(), tree, allPoints );

			if( tree.size() > 1 )
				test_knn( random.uniform( 1, tree.size() ), tree, allPoints );

			test_knn( 1, batchTree, allPoints );
			test_knn( 2, batchTree, allPoints );
			test_knn( tree.size(), batchTree, allPoints );
			if( tree.size() > 1 )
				test_knn( random.uniform( 1, batchTree.size() ), tree, allPoints );
		}

        
	}

	return 0;
}

