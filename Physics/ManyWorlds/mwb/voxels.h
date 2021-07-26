#ifndef __VOXELS_H__
#define __VOXELS_H__

#include "twigg/boundingbox.h"

#include <boost/dynamic_bitset.hpp>

namespace planning {

// This voxel grid uses a bitset as the underlying representation
//   for speed and space reasons.  
class VoxelGrid
{
public:
	VoxelGrid();
	VoxelGrid( size_t xDim, size_t yDim, size_t zDim, size_t maxValue );
	VoxelGrid( std::istream& is );

	void dump( std::ostream& os ) const;

	void set( const TinyVec<size_t, 3>& cell, size_t value )
	{
		this->set( cell[0], cell[1], cell[2], value );
	}

	void set( const size_t x, const size_t y, const size_t z, size_t value )
	{
		assert( value < (1 << numBits_) );
		checkCell( x, y, z );

		const size_t cellStart = this->cell(x, y, z);
		for( size_t i = 0; i < numBits_; ++i )
			bits_.set( cellStart + i, (value >> i) & 1 );
	}

	size_t operator()( const TinyVec<size_t, 3>& cell ) const
	{
		return (*this)( cell[0], cell[1], cell[2] );
	}

	size_t operator()( const size_t x, const size_t y, const size_t z ) const
	{
		checkCell( x, y, z );

		const size_t cellStart = this->cell(x, y, z);
		size_t result = 0;
		for( size_t i = 0; i < numBits_; ++i )
			result |= bits_[ cellStart + i ] << i;

		return result;
	}

	TinyVec<size_t, 3> dimensions() const
	{
		return dimensions_;
	}

	/*
	std::size_t hash() const
	{
		std::size_t result = 0;
		to_block_range( this->bits_,
			boost::make_function_output_iterator(
				boost::bind( &boost::hash_combine, boost::ref(result), _1 ) ) );
		return result;
	}
	*/

private:
	void checkCell( const TinyVec<size_t, 3>& cell ) const 
	{
		assert( cell[0] < dimensions_[0] &&
			cell[1] < dimensions_[1] &&
			cell[2] < dimensions_[2] );
	}

	void checkCell( const size_t x, const size_t y, const size_t z ) const
	{
		assert( x < dimensions_[0] &&
			y < dimensions_[1] &&
			z < dimensions_[2] );
	}

	size_t cell( const TinyVec<size_t, 3>& c ) const
	{
		return cell( c[0], c[1], c[2] );
	}

	size_t cell( const size_t x, const size_t y, const size_t z ) const
	{
		return (z*(dimensions_[0]*dimensions_[1]) + 
			y*(dimensions_[0]) +
			x)*numBits_;
	}

	typedef boost::dynamic_bitset<> BitSet;
	BitSet bits_;
	TinyVec<size_t, 3> dimensions_;
	size_t numBits_;
};

// This places a voxel grid in space
class BoundedVoxelGrid
{
public:
	BoundedVoxelGrid( const BoundingBox3f& bounds, const TinyVec<size_t, 3>& resolution, size_t numMaterials )
		: bounds_( bounds ), grid_( resolution[0], resolution[1], resolution[2], numMaterials )
	{
		vl::Vec3f diff = bounds_.maximum() - bounds_.minimum();
		for( vl::Int i = 0; i < 3; ++i )
			this->increment_[i] = diff[i] / boost::numeric_cast<float>( resolution[i] );
	}

	BoundedVoxelGrid( const VoxelGrid& grid, const BoundingBox3f& bounds )
		: grid_(grid), bounds_(bounds)
	{
		TinyVec<size_t, 3> resolution = grid.dimensions();
		vl::Vec3f diff = bounds_.maximum() - bounds_.minimum();
		for( vl::Int i = 0; i < 3; ++i )
			this->increment_[i] = diff[i] / boost::numeric_cast<float>( resolution[i] );
	}

	vl::Vec3f center( const TinyVec<size_t, 3>& voxel ) const
	{
		return bounds_.minimum() + vl::Vec3f(
			(boost::numeric_cast<float>( voxel[0] ) + 0.5f)*increment_[0],
			(boost::numeric_cast<float>( voxel[1] ) + 0.5f)*increment_[1],
			(boost::numeric_cast<float>( voxel[2] ) + 0.5f)*increment_[2] );
	}

	BoundingBox3f box( const TinyVec<size_t, 3>& voxel ) const
	{
		vl::Vec3f lowerLeft(
			boost::numeric_cast<float>( voxel[0] )*increment_[0],
			boost::numeric_cast<float>( voxel[1] )*increment_[1],
			boost::numeric_cast<float>( voxel[2] )*increment_[2] );
		lowerLeft += bounds_.minimum();
		return BoundingBox3f( lowerLeft, lowerLeft + increment_ );
	}

	TinyVec<size_t, 3> voxel( const vl::Vec3f& pos )
	{
		vl::Vec3f offsetPos = pos - bounds_.minimum();
		for( unsigned i = 0; i < 3; ++i )
			offsetPos[i] = std::max( offsetPos[i], 0.0f );

		TinyVec<size_t, 3> result(
			boost::numeric_cast<size_t>( offsetPos[0] / increment_[0] ),
			boost::numeric_cast<size_t>( offsetPos[1] / increment_[1] ),
			boost::numeric_cast<size_t>( offsetPos[2] / increment_[2] ) );

		const TinyVec<size_t, 3> dims = dimensions();
		// clip to grid
		for( unsigned i = 0; i < 3; ++i )
			result[i] = std::min( result[i], dims[i] - 1 );

		return result;
	}

	TinyVec<size_t, 3> dimensions() const
	{
		return grid_.dimensions();
	}

	void set( const TinyVec<size_t, 3>& voxel, size_t value )
	{
		grid_.set( voxel, value );
	}

	void set( const size_t i, const size_t j, const size_t k, size_t value )
	{
		grid_.set( i, j, k, value );
	}

	size_t operator()( const TinyVec<size_t, 3>& voxel ) const
	{
		return grid_(voxel);
	}

	size_t operator()( const size_t i, const size_t j, const size_t k ) const
	{
		return grid_(i, j, k);
	}

	BoundingBox3f bounds() const
	{
		return this->bounds_;
	}

	const VoxelGrid& voxelGrid() const
	{
		return grid_;
	}

private:
	BoundingBox3f bounds_;
	vl::Vec3f increment_;
	VoxelGrid grid_;
};

void testVoxels();

} // namespace planning

#endif

