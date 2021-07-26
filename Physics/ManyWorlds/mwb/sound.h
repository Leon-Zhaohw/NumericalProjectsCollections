#ifndef __SOUND_H__
#define __SOUND_H__

#include "physicsFwd.h"
#include "barbicUtils.h"
#include "voxels.h"
#include "threads.h"

#include "twigg/linalg.h"

#ifdef SOUND

namespace planning {

struct SquashingCubesModes;
// the Event is so we can exit early if necessary
boost::shared_ptr<SquashingCubesModes> computeModes( PhysicsObjectPtr object, Event& event );

// for i/j/k indices, the numbers will actually be
//   pretty small (<<1000), so we can save space by packing
//   them into a smaller integer type:
typedef unsigned short SmallIndex;
typedef boost::tuple<SmallIndex, SmallIndex, SmallIndex> SmallIndexTriple;
size_t findVertex( size_t i, size_t j, size_t k, 
				const std::deque<SmallIndexTriple>& vertices );

void runCalculix( const std::string& inputFile, Event& event );
std::pair< std::vector<double>, fortran_matrix > parseFRDFile( const std::string& filename );

struct SquashingCubesModes
{
	struct ShortMaterial
	{
		ShortMaterial( float ym, float pr )
			: youngsModulus(ym), poissonRatio(pr) {}

		float youngsModulus;
		float poissonRatio;
	};

	SquashingCubesModes( 
		const BoundedVoxelGrid& vox )
		: voxels(vox) {}

	BoundedVoxelGrid voxels;
	std::deque<SmallIndexTriple> vertices;

	std::vector<double> frequencies;
	fortran_matrix modes;

	std::vector<ShortMaterial> materials;
};

class AudioGenerator
{
public:
	AudioGenerator( size_t rate, 
		const std::vector<double>& frequencies, 
		double alpha, double beta );
	virtual ~AudioGenerator();

	virtual void addImpulse( const vl::Vec3f& position, const vl::Vec3f& force, float dt ) = 0;
	std::vector<double> integrate( size_t numSteps );

protected:
	void addImpulses( const std::vector<double>& impulses, double width );
	double dt() const { return this->dt_; }

	struct GaussianImpulse
	{
		GaussianImpulse( const std::vector<double>& imp, double sd )
			: impulses(imp), stdDev(sd), currentPosition( -3.0*sd ) {}

		std::vector<double> impulses;
		double currentPosition;
		double stdDev;
	};

	typedef std::list<GaussianImpulse> ImpulseList;
	ImpulseList impulses_;

private:
	size_t rate_;
	double dt_;
	std::vector<barbic::HarmonicOscillator_IIR> oscillators_;
};

class SquashingCubesAudioGenerator
	: public AudioGenerator
{
public:
	// alpha and beta are parameters for Raleigh damping
	// see http://www.cs.cornell.edu/~djames/papers/DyRT.pdf
	SquashingCubesAudioGenerator( 
		size_t rate,
		boost::shared_ptr<SquashingCubesModes> modes,
		double alpha,
		double beta );
	void addImpulse( const vl::Vec3f& position, const vl::Vec3f& force, float dt );

private:
	boost::shared_ptr<SquashingCubesModes> modes_;
};

template <typename IntegerType>
void createWavFile( const std::string& filename, const std::vector<IntegerType>& samples );

typedef std::pair< std::string, boost::array<float, 6> > AudioHash;

struct CompareAudioHash
{
	bool operator()( const AudioHash& left, const AudioHash& right ) const;
};

// We want to be able ot automatically detect if two objects in the scene are identical
//   for the purposes of computing audio, as this will save us a lot of time in the 
//   long run.  
AudioHash hashForAudio( ConstPhysicsObjectPtr physicsObject );

} // namespace planning

#endif // SOUND

#endif
