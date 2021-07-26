#ifdef WIN32
#pragma once
#endif

#include "twigg/boundingbox.h"
#include "twigg/fastTableLookup.h"

namespace planning {

class DistanceField
{
public:
	DistanceField( const std::string& filename );
	~DistanceField();

	// might as well return both the value and the gradient in the 
	//   same call:
	std::pair<float, vl::Vec3f> value( const vl::Vec3f& position ) const;

private:
	float value( unsigned int i, unsigned int j, unsigned int k ) const;

	std::vector<float> data_;
	boost::array<size_t, 3> dimensions_;
	boost::array<size_t, 3> strides_;
	BoundingBox3f bounds_;
	vl::Vec3f delta_;

	FastTableLookup lookup_;
};

} // namespace planning

