#ifndef __COMPRESS_H__
#define __COMPRESS_H__

#include "scene.h"
#include "quat.h"

#include "twigg/random.h"

namespace planning {

vl::Vec3d differentialRotation( const Quaternion<double>& left, const Quaternion<double>& right );


struct CompressedPath
{
    int startFrame;

	std::vector<char> positionTimes;
	boost::array< std::vector<char>, 3 > positionLinearCoeffs;
	boost::array< std::vector<char>, 3 > positionQuadraticCoeffs;

	std::vector<char> rotationTimes;
	boost::array< std::vector<char>, 3 > rotationCoeffs;

	size_t size() const
	{
		size_t result = 0;
		result += positionTimes.size();
		for( size_t i = 0; i < 3; ++i )
			result += positionLinearCoeffs[i].size();
		for( size_t i = 0; i < 3; ++i )
			result += positionQuadraticCoeffs[i].size();

		result += rotationTimes.size();
		for( size_t i = 0; i < 3; ++i )
			result += rotationCoeffs[i].size();

		return result;
	}
};

// This is a path that has been decompressed
//   but is still in piecewise form
class PiecewisePath
{
public:
    typedef int FrameType;
    
    PiecewisePath();
	PiecewisePath( const CompressedPath& path, const PiecewisePath& other = PiecewisePath(), bool backwards = false );
	PiecewisePath( const std::vector<FrameType>& positionTimes,
		const std::vector<float>& positionCoeffs,
		const std::vector<FrameType>& rotationTimes,
		const std::vector<float>& rotations,
		const std::vector<float>& rotationCoeffs );

	vl::Vec3f position( const float time ) const;
	Quaternion<double> rotation( const float time ) const;

	vl::Vec3f linearVelocity( const float time ) const;
	vl::Vec3f angularVelocity( const float time ) const;

	FrameType startFrame() const;
	FrameType endFrame() const;

	// Dumps out as Maya animation
	void toMayaAsciiFile( std::ostream& ofs, const std::string& name, size_t frameRate, 
        FrameType startFrame = boost::numeric::bounds<FrameType>::lowest(), 
        FrameType endFrame = boost::numeric::bounds<FrameType>::highest() ) const;

	// dumps the path out as a NURBS curve
	void pathToMayaAsciiFile( std::ostream& ofs, const std::string& name, const vl::Vec3d& pointInLocalFrame ) const;

	void dumpToArrays( 
		std::vector<FrameType>& positionTimes,
		std::vector<float>& positionCoeffs,
		std::vector<FrameType>& rotationTimes,
		std::vector<float>& rotations,
		std::vector<float>& rotationCoeffs ) const;

	float angularEnergy() const;
private:

	typedef vl::Vec3f CoeffList;
	typedef std::pair<FrameType, boost::array<CoeffList, 3> > TimedCoeffList;
	std::vector<TimedCoeffList> positionCoeffs_;

	typedef boost::tuple<FrameType, Quaternion<double>, vl::Vec3f> TimedDifferentialRotation;
	std::vector<TimedDifferentialRotation> rotationCoeffs_;

	struct CompareTimedRot
	{
		bool operator()( const TimedDifferentialRotation& left, const TimedDifferentialRotation& right ) const
		{
			return (left.get<0>() < right.get<0>());
		}
	};

};

CompressedPath compress( 
	const std::deque<RigidDynamicState>& states, 
	float quadConstant,
	size_t startOffset,
    size_t frameRate,
	const float maxPositionError = 1e-3,
	const float maxRotationError = 1e-3 );

std::vector<CompressedPath> compress( const std::deque<Simulation::State>& states, 
	const float quadConstant,
    size_t frameRate,
	const float maxPositionError = 1e-3,
	const float maxRotationError = 1e-3 );

void testCompression();

void computeRateDistortionCurve( const std::deque<Simulation::State>& states, 
	const float quadConstant, const std::string& matFilename );

} // namespace planning

#endif

