#include "stdafx.h"
#include "distanceField.h"

#include <fstream>

namespace planning {

DistanceField::DistanceField( const std::string& filename )
	: lookup_(1)
{
	std::ifstream ifs( filename.c_str(), std::ios::binary );
	if( !ifs )
		throw IOException( "Unable to read file '" + filename + "' for reading." );
	ifs.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

	try
	{
		int resolution[3];
		ifs.read( reinterpret_cast<char*>( resolution ), 3*sizeof(int) );
		assert( (resolution[0] == resolution[1]) && (resolution[1] == resolution[2]) );
		std::copy( resolution, resolution+3, dimensions_.begin() );
		lookup_ = FastTableLookup(resolution[0]);

		vl::Vec3d lowerBound;
		ifs.read( reinterpret_cast<char*>( lowerBound.Ref() ), 3*sizeof(double) );
		this->bounds_.expand( toVec3f(lowerBound) );

		vl::Vec3d upperBound;
		ifs.read( reinterpret_cast<char*>( upperBound.Ref() ), 3*sizeof(double) );
		this->bounds_.expand( toVec3f(upperBound) );

		std::vector<double> data( resolution[0]*resolution[1]*resolution[2] );
		ifs.read( reinterpret_cast<char*>( &data[0] ), data.size()*sizeof(double) );
		data_.resize( data.size() );
		std::copy( data.begin(), data.end(), data_.begin() );
	}
	catch( std::ifstream::failure& e )
	{
		std::ostringstream oss;
		oss << "Unable to read signed distance data from '" << filename << "': " << e.what();
		throw IOException( oss.str() );
	}

	delta_ = (bounds_.maximum() - bounds_.minimum());
	delta_[0] /= static_cast<float>(dimensions_[0] - 1);
	delta_[1] /= static_cast<float>(dimensions_[1] - 1);
	delta_[2] /= static_cast<float>(dimensions_[2] - 1);

	strides_[0] = 1;
	for ( int i = 1; i < 3; ++i )
		strides_[i] = strides_[i-1] * dimensions_[i-1];

}

DistanceField::~DistanceField()
{
}

float DistanceField::value( unsigned int i0, unsigned int i1, unsigned int i2 ) const
{
	return data_.at( i0 * strides_[0] + i1 * strides_[1] + i2 * strides_[2] );
}

std::pair<float, vl::Vec3f> DistanceField::value( const vl::Vec3f& position ) const
{
	boost::array<unsigned int, 3> indices;
	vl::Vec3f offsets;

	vl::Vec3f normalizedPos = (position - this->bounds_.minimum()) / delta_;
	// convert into a usable position for trilinear interp
	for( size_t i = 0; i < 3; ++i )
		this->lookup_.tableLookupIndices( normalizedPos[i], indices[i], offsets[i] );

	float dist = 0.0f;
	dist += (1.0f-offsets[0])*(1.0f-offsets[1])*(1.0f-offsets[2])*value(indices[0], indices[1], indices[2]);
	dist += (offsets[0])*(1.0f-offsets[1])*(1.0f-offsets[2])*value(indices[0] + 1, indices[1], indices[2]);
	dist += (1.0f-offsets[0])*(offsets[1])*(1.0f-offsets[2])*value(indices[0], indices[1] + 1, indices[2]);
	dist += (1.0f-offsets[0])*(1.0f-offsets[1])*(offsets[2])*value(indices[0], indices[1], indices[2] + 1);
	dist += (offsets[0])*(1.0f-offsets[1])*(offsets[2])*value(indices[0] + 1, indices[1], indices[2] + 1);
	dist += (1.0f-offsets[0])*(offsets[1])*(offsets[2])*value(indices[0], indices[1] + 1, indices[2] + 1);
	dist += (offsets[0])*(offsets[1])*(1.0f-offsets[2])*value(indices[0] + 1, indices[1] + 1, indices[2]);
	dist += (offsets[0])*(offsets[1])*(offsets[2])*value(indices[0] + 1, indices[1] + 1, indices[2] + 1);

	vl::Vec3f gradient;
	gradient[0] = 
		(1.0f-offsets[1])*(1.0f-offsets[2])*
				(value(indices[0]+1, indices[1], indices[2]) - value(indices[0], indices[1], indices[2])) +
		(offsets[1])*(1.0f-offsets[2])*
				(value(indices[0]+1, indices[1]+1, indices[2]) - value(indices[0], indices[1]+1, indices[2])) +
		(1.0f-offsets[1])*(offsets[2])*
				(value(indices[0]+1, indices[1], indices[2]+1) - value(indices[0], indices[1], indices[2]+1)) +
		(offsets[1])*(offsets[2])*
				(value(indices[0]+1, indices[1]+1, indices[2]+1) - value(indices[0], indices[1]+1, indices[2]+1));
	gradient[1] = 
		(1.0f-offsets[0])*(1.0f-offsets[2])*
				(value(indices[0], indices[1]+1, indices[2]) - value(indices[0], indices[1], indices[2])) +
		(offsets[0])*(1.0f-offsets[2])*
				(value(indices[0]+1, indices[1]+1, indices[2]) - value(indices[0]+1, indices[1], indices[2])) +
		(1.0f-offsets[0])*(offsets[2])*
				(value(indices[0], indices[1]+1, indices[2]+1) - value(indices[0], indices[1], indices[2]+1)) +
		(offsets[0])*(offsets[2])*
				(value(indices[0]+1, indices[1]+1, indices[2]+1) - value(indices[0]+1, indices[1], indices[2]+1));
	gradient[2] = 
		(1.0f-offsets[0])*(1.0f-offsets[1])*
				(value(indices[0], indices[1], indices[2]+1) - value(indices[0], indices[1], indices[2])) +
		(offsets[0])*(1.0f-offsets[1])*
				(value(indices[0]+1, indices[1], indices[2]+1) - value(indices[0]+1, indices[1], indices[2])) +
		(1.0f-offsets[0])*(offsets[1])*
				(value(indices[0], indices[1]+1, indices[2]+1) - value(indices[0], indices[1]+1, indices[2])) +
		(offsets[0])*(offsets[1])*
				(value(indices[0]+1, indices[1]+1, indices[2]+1) - value(indices[0]+1, indices[1]+1, indices[2]));

	return std::make_pair( dist, gradient );
}

} // namespace planning

