#include "stdafx.h"

#include "spatial_tree.h"

#include "twigg/random.h"
#include "twigg/boundingbox.h"

namespace planning {

typedef spatial_tree< float, 2, size_t > TreeType;

void test_tree( const TreeType& tree, const std::deque<vl::Vec2f>& elts, BoundingBox2f bounds )
{
	float boxSize = ldexp( 1.0, -boost::numeric_cast<int>(tree.max_depth()) );

	TreeType::data_set treeRes = tree.query( bounds.minimum().Ref(), bounds.maximum().Ref(), true );
	std::deque<size_t> treeResVec( treeRes.begin(), treeRes.end() );
	std::sort( treeResVec.begin(), treeResVec.end() );

	std::deque<size_t> bruteForceRes;
	std::deque<size_t> bruteForceResBiggest;
	for( size_t i = 0; i < elts.size(); ++i )
	{
		if( contains( bounds, elts[i] ) )
			bruteForceRes.push_back( i );

		if( contains( bounds, elts[i], boxSize ) )
			bruteForceResBiggest.push_back( i );
	}

	for( std::deque<size_t>::const_iterator itr = treeResVec.begin();
		itr != treeResVec.end(); ++itr )
	{
		assert( std::binary_search(bruteForceResBiggest.begin(), bruteForceResBiggest.end(), *itr) );
	}

	for( std::deque<size_t>::const_iterator itr = bruteForceRes.begin();
		itr != bruteForceRes.end(); ++itr )
	{
		assert( std::binary_search(treeResVec.begin(), treeResVec.end(), *itr) );
	}
}

void test_spatial_tree()
{
	RandomStream random;

	for( size_t depth = 2; depth < 8; ++depth )
	{
		for( size_t iNumber = 0; iNumber < 20; ++iNumber )
		{
			size_t actualNumber = iNumber;
			if( iNumber > 10 )
				actualNumber = random.uniform( 10, 1000 );

			TreeType tree( depth );
			std::deque< vl::Vec2f > elts;
			for( size_t i = 0; i < actualNumber; ++i )
			{
				elts.push_back( vl::Vec2f( random.uniform(), random.uniform() ) );
				tree.insert( elts.back(), i );
 
				// test inserting something outside the tree, should just get dropped
				tree.insert( vl::Vec2f(2.0, 2.0), actualNumber+i );
				tree.insert( vl::Vec2f(-0.1, -0.1), actualNumber+i );
			}

			// full box
			test_tree( tree, elts, 
				BoundingBox2f( vl::Vec2f(0.0, 0.0), vl::Vec2f(1.0, 1.0) ) );

			// empty box
			test_tree( tree, elts, BoundingBox2f() );

			for( size_t iTest = 0; iTest < 10; ++iTest )
			{
				BoundingBox2f box;
				for( size_t i = 0; i < 2; ++i )
					box.expand( vl::Vec2f( random.uniform(), random.uniform() ) );

				test_tree( tree, elts, box );

				// test a box that (probably) extends beyond the tree:
				box.expand( vl::Vec2f( (random.uniform()-0.5) * 5.0, (random.uniform()-0.5) * 5.0) );

				test_tree( tree, elts, box );
			}
		}
	}
}


} // namespace planning
