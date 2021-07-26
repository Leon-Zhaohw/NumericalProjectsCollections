#ifdef WIN32
#pragma once
#endif

#ifndef __SIMULATIONTREE_H__
#define __SIMULATIONTREE_H__

#include "physicsObject.h"
#include "scene.h"
#include "compression.h"
#include "pathTree.h"
#include "OBB.h"
#include "threads.h"

#include "twigg/cover_tree.h"
#include "twigg/linalg.h"
#include "twigg/random.h"

#include <list>

#ifdef CONDOR_SAMPLING
class MWRMComm;
#endif

namespace planning {

class sphere_tree;

vl::Vec3f randomDirection(const vl::Vec3f& source, double xi_1, double xi_2, double n);


// a list of objects that you need to check out from:
template <typename T>
class CheckoutList
{
public:
	template <typename List>
	CheckoutList( const List& list )
	{
		std::copy( list.begin(), list.end(), 
			std::back_inserter( this->list_ ) );
	}

	CheckoutList()
	{
	}

	template <typename List>
	void add( const List& list )
	{
		boost::mutex::scoped_lock lock( this->mutex_ );
		std::copy( list.begin(), list.end(), 
			std::back_inserter(this->objects_) );
	}

	~CheckoutList()
	{
	}

	friend class Acquire;
	class Acquire
	{
	public:
		Acquire(CheckoutList<T>& list)
			: list_(list)
		{
			boost::mutex::scoped_lock lock( list.mutex_ );
			while( list_.objects_.empty() )
				list.nonemptyCondition_.wait( lock );

			this->object_ = list_.objects_.front();
			list_.objects_.pop_front();
		}

		~Acquire()
		{
			boost::mutex::scoped_lock lock( list_.mutex_ );
			list_.objects_.push_back( object_ );
			list_.nonemptyCondition_.notify_one();
		}

		SimulationPtr operator()()
		{
			return object_;
		}

	private:
		T object_;
		CheckoutList<T>& list_;
	};

private:
	std::list<SimulationPtr> objects_;
	boost::mutex mutex_;
	boost::condition nonemptyCondition_;
};

class Path
{
public:
	typedef std::vector<RigidStaticState> StateSnapshot;
	typedef std::pair<float, StateSnapshot> TimedSnapshot;

	Path( const std::vector<CompressedPath>& paths, 
		const std::vector<TimedSnapshot>& states, 
		const std::vector< std::vector<vl::Vec3f> >& cmPaths,
		short frameRate, 
		const std::vector<float>& scores);
	Path( const std::vector<CompressedPath>& paths, 
		const std::deque<Simulation::State>& states, 
		const std::vector< std::vector<vl::Vec3f> >& cmPaths,
		short frameRate, 
		const std::vector<float>& scores,
		float epsilon = 1e-3 );

	Path( std::istream& is, size_t numPaths, short frameRate, const std::string& filename );
	Path( const Path& other );
	~Path();

	size_t frameRate() const;

	std::vector<Simulation::State> samples(size_t numSamples) const;
	std::vector<Simulation::State> samples() const;

	const CompressedPath& path(size_t iObject) const;

	void write( std::ostream& os ) const;

	size_t treeCount() const;
	const OBBTree& obbTree( size_t iObject ) const;

	size_t snapshotCount() const                 { return this->snapshots_.size(); }
	float snapshotTime( size_t iSnapshot ) const { return this->snapshots_.at(iSnapshot).first; }
	std::vector<RigidStaticState> snapshot( size_t iSnapshot ) const
		{ return this->snapshots_.at(iSnapshot).second; }

	float score( size_t iObject, size_t iMetric ) const;
	size_t metricCount() const;
	std::vector<float> allScores() const { return this->scores_; }
private:
	CompressedPath readPath( std::istream& is ) const;

	mutable std::vector<CompressedPath> paths_;

	// we will store a sampling of states to be able to quickly
	//   render a snapshot
	std::vector<TimedSnapshot> snapshots_;
	short frameRate_;

	boost::ptr_vector<OBBTree> obbTrees_;
	std::vector<float> scores_;

	std::string filename_;
	std::vector< std::streampos > pathPositions_;
};


class SimulationNodeCallback;

class SimulationTree
	: public Named
{
public:
	SimulationTree( const std::string& name,
		boost::shared_ptr<Scene> scene, 
		size_t frameRate,
		boost::shared_ptr<SimulationFactory> simulationFactory,
        bool timeReversed = false );

	// refinement when all we have is the compressed state:
	SimulationTree( const std::string& name,
		ScenePtr scene,
		boost::shared_ptr<SimulationTree> parent,
		std::deque<ConstPhysicsObjectPtr> activeObjects, 
		const std::vector<PiecewisePath>& piecewisePaths,
		size_t piecewisePathsFrameRate,
		float startTime,
		bool timeReversed = false );

	~SimulationTree();

	// we will need this to be able to interpret what is in the 
	//   SimulationNodes:
	std::vector<ConstPhysicsObjectPtr> dynamicObjects() const;

	boost::shared_ptr<const Scene> scene() const;
	boost::shared_ptr<const Scene> parentScene() const;

	float dt() const
	{
		return this->dt_;
	}

	float startTime() const
	{
		return this->startState_.time();
	}

	void read( const std::string& filename );

	const CameraSource& cameraSource() const;

	void sample( RandomStream& random, bool backwards );

	// @todo: make this private
	// expand a node until we find another choice point
	Path* sample( 
		size_t maxTimesteps,                          // max. number of timesteps to run for
		size_t minTimesteps,
		RandomStream& random );             // says if we're allowed to add impulses to the passed in node

    std::vector<Simulation::State> sampleForwardToState( 
        SimulationPtr simulation,
        size_t numTimesteps,
        RandomStream& random );

	std::vector<ObjectivePtr> objectives() const
	{
		return this->objectives_;
	}

	boost::dynamic_bitset<> addConstraint( ConstraintPtr constraint, const boost::dynamic_bitset<>& hint = boost::dynamic_bitset<>() );
	void removeConstraint( ConstraintPtr constraint );

	boost::dynamic_bitset<> activePaths() const;
	std::deque<ConstraintPtr> constraints() const;
	const Path* path( size_t pathId ) const;
	size_t numPaths() const;
	size_t idForPath( const Path* path ) const;
	size_t constraintCount() const;
	bool satisfiesConstraint( size_t pathId, size_t constraintId ) const;

	void addCallback( SimulationNodeCallback* callback );
	void removeCallback( SimulationNodeCallback* callback );

    void setTimeRange( std::pair<float, float> timeRange );

	size_t samplingRate() const { return this->samplingRate_; }

	// should be called in application's idle loop
	bool checkForNewPaths();
	std::vector<PiecewisePath> decode( const Path* path );

	std::vector<size_t> fixedObjects() const;
	std::vector<PiecewisePath> fixedObjectsHistory() const
	{
		return this->fixedObjectsHistory_;
	}

	void addPath( 
		const std::vector<CompressedPath>& paths, 
		const std::vector<Path::TimedSnapshot>& states, 
		const std::vector< std::vector<vl::Vec3f> >& cmPaths,
		short frameRate,
		const std::vector<float>& objectives );
	void addPath( const Path& path );

	std::vector<vl::Vec3f> pathRoot(size_t iObject) const;

    void setSamplingProperties( size_t iDynamicObject, SamplingProperties props );
    SamplingProperties getSamplingProperties( size_t iDynamicObject ) const;

private:
	std::vector<ConstPhysicsObjectPtr> dynamicObjects_;
	std::vector<ConstPhysicsObjectPtr> parentDynamicObjects_;
	std::vector< std::vector<vl::Vec3f> > pathRoots_;

    std::vector<SamplingProperties> samplingProperties_;
    mutable boost::mutex samplingPropertiesMutex_;

	// we will frequently need to figure out where in the parent's
	//  set of dynamic objects a particular object lies:
	typedef STLEXT hash_map<size_t, size_t> DynamicObjectIdMap;
	DynamicObjectIdMap dynamicObjectToParentDynamicObject_;
	// we will use this to prime the set of parent objects that are "active"
	std::vector< ConstPhysicsObjectPtr > onlyInParentObjects_;
	std::vector< ConstPhysicsObjectPtr > onlyInChildObjects_;
	// Need to know an object's 'doppelganger' so we can activate it at the 
	//   appropriate time:
	typedef STLEXT hash_map<const PhysicsObject*, const PhysicsObject*, 
		hash_ptr<const PhysicsObject*> > ObjectMap;
	ObjectMap childObjectToParentObject_;
	ObjectMap parentObjectToChildObject_;

	vl::Vec3f randomImpulse(RandomStream& random);
	vl::Vec3f randomImpulse(const vl::Vec3f& source, RandomStream& random);
	vl::Vec3f randomImpulse(const vl::Vec3f& source, double xi_1, double xi_2, double nu);

	boost::shared_ptr<Scene> scene_;

	// this is the actual root of the tree
	Simulation::State startState_;

	float dt_;

	// we will refuse to expand past this time, will be useful for constraining
	//   animations without endings
    // note that this is a range so that we know where to stop if
    //   we're simulating backwards as well
    std::pair<float, float> timeRange_;

	// the minimum impulses we will consider "consequential," in terms
	//   of making something a choice point
	float minLinearImpulse_;
	float minAngularImpulse_;

	std::vector<ObjectivePtr> objectives_;

	// a body is disabled if the _square_ of its linear/angular velocity
	//   drops below this threshold for two consecutive timesteps
	float linearVelocityThreshold_;
	float angularVelocityThreshold_;

	boost::shared_ptr<CameraSource> cameraSource_;
	size_t samplingRate_;
	size_t simulationRate_;
	size_t multiplier_;

	// parent scene
	boost::shared_ptr<Scene> parentScene_;

	std::vector<size_t> fixedObjects_;
	std::vector<PiecewisePath> fixedObjectsHistory_;
	size_t fixedObjectsHistoryFramerate_;

	// we are going to maintain sets of _all_ possible nodes (anything that
	//   runs all the way out) and _active_ nodes, which are nodes that 
	//   satisfy all current constraints
	boost::dynamic_bitset<> activePaths_;
	mutable ReaderWriterLock activePathsLock_;

	std::deque<const Path*> allPaths_;
	mutable ReaderWriterLock allPathsLock_;

	typedef STLEXT hash_map< const Path*, size_t, hash_ptr<const Path*> > PathToIdMap;
	PathToIdMap pathToId_;

	// all constraints
	std::deque<ConstraintPtr> constraints_;
	mutable ReaderWriterLock constraintsLock_;

	std::deque< boost::dynamic_bitset<> > constraintResults_;

	/*
	boost::ptr_vector<PathTree> pathTrees_;
	mutable ReaderWriterLock pathTreesLock_;
	*/

	boost::mutex nodesToAddMutex_;
	std::deque<const Path*> pathsToAdd_;

	boost::object_pool<Path> pathPool_;
	mutable boost::mutex pathPoolMutex_;

	std::list<SimulationNodeCallback*> callbacks_;
	boost::mutex callbackMutex_;

	boost::shared_ptr<SimulationFactory> simulationFactory_;

	// Determines if the tree is simulating backwards
	bool timeReversed_;

    std::vector<SamplingProperties> backwardsSamplingProperties_;

    size_t numSamples_;
    
#ifdef STREAMING
	std::ofstream ofs_;
	size_t counter_;

	// gonna dump out paths in a .mat file for analysis
	std::deque< std::vector<vl::Vec3f> > pathsForMatFile_;
#endif
};

// We will use this callback to let others know when nodes have been added or removed.
//   For now the "delete" interface is very simple which helps us elsewhere.
class SimulationNodeCallback
{
public:
	virtual ~SimulationNodeCallback() {}
	virtual void addNodes( const SimulationTree* tree, const std::deque<const Path*>& nodes, const boost::dynamic_bitset<>& active ) = 0;
	virtual void setNodes( const SimulationTree* tree, const boost::dynamic_bitset<>& nodes ) = 0;
};

class SimulationNodeCallbackRegistrar
{
public:
	SimulationNodeCallbackRegistrar( SimulationTreePtr tree, SimulationNodeCallback* callback )
		: callback_(callback), tree_(tree)
	{
		tree_->addCallback( callback_ );
	}

	~SimulationNodeCallbackRegistrar()
	{
		tree_->removeCallback( callback_ );
	}

	ConstSimulationTreePtr tree() const
	{
		return tree_;
	}

private:
	SimulationTreePtr tree_;
	SimulationNodeCallback* callback_;
};

} // namespace planning


#endif
