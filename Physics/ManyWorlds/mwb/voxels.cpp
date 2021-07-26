#include "stdafx.h"

#include "voxels.h"
#include "twigg/random.h"
#include "twigg/linalg.h"

namespace planning {

VoxelGrid::VoxelGrid( size_t xDim, size_t yDim, size_t zDim, size_t maxValue )
{
	++maxValue;
	numBits_ = 1;
	while( maxValue >>= 1 )
		++numBits_;

	bits_.resize( xDim*yDim*zDim*numBits_ );
	dimensions_[0] = xDim;
	dimensions_[1] = yDim;
	dimensions_[2] = zDim;
}

VoxelGrid::VoxelGrid( std::istream& is )
{
	is.read( reinterpret_cast<char*>( &numBits_ ), sizeof(size_t) );
	is.read( reinterpret_cast<char*>( &dimensions_ ), sizeof(TinyVec<size_t, 3>) );

	std::vector<BitSet::block_type> blocks;
	readArray( blocks, is );
	bits_.swap( BitSet( blocks.begin(), blocks.end() ) );

	const size_t totalBits = numBits_*dimensions_[0]*dimensions_[1]*dimensions_[2];
	assert( bits_.size() >= totalBits );
	bits_.resize( totalBits );
}

VoxelGrid::VoxelGrid()
	: numBits_(1)
{
	for( size_t i = 0; i < 3; ++i )
		dimensions_[i] = 0;
}

void VoxelGrid::dump( std::ostream& os ) const
{
	os.write( reinterpret_cast<const char*>( &numBits_ ), sizeof(size_t) );
	os.write( reinterpret_cast<const char*>( &dimensions_ ), sizeof(TinyVec<size_t, 3>) );
	std::vector<BitSet::block_type> blocks( bits_.num_blocks() );
	to_block_range( bits_, blocks.begin() );
	dumpArray( blocks, os );
}

void testVoxels(size_t x, size_t y, size_t z, RandomStream& rand)
{
	// different numbers of possible values
	for( size_t numValues = 1; numValues < 10; ++numValues )
	{
		VoxelGrid grid( x, y, z, numValues );
		size_t n = 0;
		for( size_t i = 0; i < x; ++i )
			for( size_t j = 0; j < y; ++j )
				for( size_t k = 0; k < z; ++k )
					grid.set( i, j, k, (n++) % numValues );

		n = 0;
		for( size_t i = 0; i < x; ++i )
		{
			for( size_t j = 0; j < y; ++j )
			{
				for( size_t k = 0; k < z; ++k )
				{
					size_t gridVal = grid( i, j, k );
					assert( gridVal == (n % numValues) );
					++n;
				}
			}
		}
	}
}

void testVoxels()
{
	RandomStream rand;

	// all ones
	testVoxels( 1, 1, 1, rand );

	// all primes
	testVoxels( 3, 7, 23, rand );

	// all squares
	testVoxels( 4, 16, 25, rand );

	// a mix
	testVoxels( 4, 1, 7, rand );
}

} // namespace planning

