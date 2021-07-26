#include "stdafx.h" // NOT_ME
#include "simulationTree.h"
#include "constraints.h"
#include "sphere_tree.h"
#include "quat.h"
#include "compression.h"

#include "twigg/AABBTree.h"
#include "twigg/renderUtils.h"
#include "twigg/linalg.h"
#include "twigg/ioutil.h"
#include "twigg/stlext.h"


#include "Wm4DistVector3Line3.h"
#include "Wm4DistLine3Line3.h"

#include <numeric>
#include <iostream>
#include <mkl_lapack.h>

// we need to use file locking for some output stuff
#include <cstdio>
#ifndef _WINDOWS
#include <cunistd>
#endif

#ifdef GUI
#include <wx/log.h>
#endif


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/connected_components.hpp>

#ifdef WIN32
#include <crtdbg.h>
#endif

namespace planning {

template <typename Pr>
class Compare1st
{
public:
	bool operator()( const Pr& left, const Pr& right ) const
	{
		return std::less<typename Pr::first_type>()( left.first, right.first );
	}
};

class PrecomputedCameraSource
	: public CameraSource
{
public:
	PrecomputedCameraSource( CameraWrapperPtr camera, size_t sampling )
	{
		// back it up first
		HasKeyable<float>::AttributeList attributes = 
			camera->camera()->keyable();
		std::vector<float> backup; backup.reserve( attributes.size() );
		std::transform( attributes.begin(), attributes.end(), std::back_inserter(backup), 
			boost::bind( &KeyableAttribute<float>::get, _1 ) );

		camera->camera()->setRatio( 1.5 );

		size_t iFrame = 0;
		float firstFrame = camera->firstKey();
		float lastFrame = camera->lastKey();

		while( (boost::numeric_cast<float>(iFrame) / boost::numeric_cast<float>(sampling)) < firstFrame )
			++iFrame;

		while( true )
		{
			float time = (boost::numeric_cast<float>(iFrame) / boost::numeric_cast<float>(sampling));
			camera->setTime( time );

			vl::Mat4f projectiveTrans = camera->camera()->projectiveTransform();
			vl::Mat4f modelViewTrans = camera->camera()->modelViewTransform();

			projectionMatrices_.push_back( MatWithTime(time, projectiveTrans) );
			modelViewMatrices_.push_back( MatWithTime(time, modelViewTrans) );
			invModelViewMatrices_.push_back( MatWithTime(time, vl::inv(modelViewTrans)) );

			if( time > lastFrame )
				break;

			++iFrame;
		}

		for( size_t i = 0; i < attributes.size(); ++i )
			attributes[i]->set( backup[i] );
	}

	virtual ~PrecomputedCameraSource() {}

	vl::Mat4f projectionMatrix( float time ) const        { return findMatrix( time, projectionMatrices_ );   }
	vl::Mat4f modelViewMatrix( float time ) const         { return findMatrix( time, modelViewMatrices_ );    }
	vl::Mat4f inverseModelViewMatrix( float time ) const  { return findMatrix( time, invModelViewMatrices_ ); }

private:
	typedef std::pair<float, vl::Mat4f> MatWithTime;
	typedef std::vector<MatWithTime> TimedMatList;

	TimedMatList projectionMatrices_;
	TimedMatList modelViewMatrices_;
	TimedMatList invModelViewMatrices_;

	vl::Mat4f findMatrix( float time, const TimedMatList& list ) const
	{
		typedef TimedMatList::const_iterator Iterator;
		Iterator iter = std::lower_bound( list.begin(), list.end(), MatWithTime(time, vl::Mat4f()), Compare1st<MatWithTime>() );
		if( iter == list.end() )
			return list.back().second;
		else
			return iter->second;
	}
};

Simulation::State startState( const std::vector<ConstPhysicsObjectPtr>& objects )
{
	// need to generate a start state for all those objects
	std::vector<RigidDynamicState> states;
	for( size_t iDynOb = 0; iDynOb < objects.size(); ++iDynOb )
	{
		vl::Mat4f xf = objects[iDynOb]->rigidTransform();
		vl::Vec3f position = vl::xform( xf, vl::Vec3f(vl::vl_0) );
		Quaternion<float> rotation( xf );
		states.push_back( 
			RigidDynamicState(
				position,
				rotation,
				objects[iDynOb]->linearVelocity(),
				objects[iDynOb]->angularVelocity() ) );
	}

	return Simulation::State( states, 0.0f );
}

SimulationTree::SimulationTree( const std::string& name, boost::shared_ptr<Scene> scene,
							   size_t frameRate, boost::shared_ptr<SimulationFactory> simulationFactory,
                               bool timeReversed )
	:	Named(name),
		scene_( scene ),
		samplingRate_(frameRate),
		timeRange_( 0.0f, 30.0f ),
		simulationFactory_( simulationFactory ),
        timeReversed_(timeReversed),
        numSamples_(0)
{
	this->fixedObjectsHistoryFramerate_ = frameRate;

	{
		this->simulationRate_ = this->samplingRate_;
		this->multiplier_ = 1;
		while( this->multiplier_ * this->simulationRate_ < 240 )
			++this->multiplier_;

		this->simulationRate_ = multiplier_ * this->samplingRate_;
		this->dt_ = 1.0 / boost::numeric_cast<double>( this->samplingRate_ );
	}

	angularVelocityThreshold_ = 0.001;
	linearVelocityThreshold_ = 0.001;

	objectives_.push_back(
		ObjectivePtr( new AngularEnergyObjective() ) );
	objectives_.push_back(
		ObjectivePtr( new FinalOrientationObjective() ) );	/*
	objectives_.push_back(
		ObjectivePtr( new ViewDependentAngularEnergyObjective() ) );
	objectives_.push_back(
		ObjectivePtr( new RunTimeObjective() ) );
	objectives_.push_back(
		ObjectivePtr( new CountCollisionsObjective() ) );
	objectives_.push_back(
		ObjectivePtr( new CollisionCoverageObjective() ) );
	objectives_.push_back(
		ObjectivePtr( new StackingObjective() ) );
	objectives_.push_back(
		ObjectivePtr( new YPositionObjective() ) );
		*/

	cameraSource_.reset( new PrecomputedCameraSource(scene_->camera(), this->samplingRate_) );

	/*
	// minimum impulse will correspond to the integration of 
	//   gravity's acceleration over the maxImpulseTime time period:
	minLinearImpulse_ = maxImpulseTime_ * 9.81;
	// not sure what to put here; I've arbitrarily chosen maxImpulseTime * 1 rev/s/s
	minAngularImpulse_ = 2.0 * maxImpulseTime_ * M_PI;
	*/
	// minimum impulse will correspond to the integration of 
	//   gravity's acceleration over the maxImpulseTime time period:
	minLinearImpulse_ = 2.0 * dt_ * 9.81;
	// not sure what to put here; I've arbitrarily chosen maxImpulseTime * 1 rev/s/s
	minAngularImpulse_ = 2.0 * dt_ * M_PI;

	this->dynamicObjects_ = physicsObjects( *this->scene_ );
    this->samplingProperties_.resize( dynamicObjects_.size() );
	this->startState_ = startState( this->dynamicObjects_ );

	/*
	this->pathTrees_.reserve( dynamicObjects_.size() );
	for( size_t i = 0; i < dynamicObjects_.size(); ++i )
		this->pathTrees_.push_back( new PathTree(i) );
		*/

#ifdef STREAMING
	ofs_.open( name.c_str(), std::ios::binary );
	if( !ofs_ )
		throw IOException( "Unable to open file '" + name + "' for writing." );

	// types:
	//   1: old-style tree dump out, dumps out individual nodes
	//   2: old-style 'sompressed' output; dumps back-to-back ShortStates
	//   3: zlibbed version of (2)
	//   101: new-style compressed output
	//   102: new-style compressed output, zlibbed

	size_t type = 101;

	std::cout << "Writing type..." << std::flush;
	ofs_.write( reinterpret_cast<const char*>(&type), sizeof(size_t) );
	std::cout << "done." << std::endl;

	size_t numObjects = this->dynamicObjects_.size();
	size_t fr = this->samplingRate_;

	std::cout << "Writing header..." << std::flush;
	ofs_.write( reinterpret_cast<const char*>( &numObjects ), sizeof(size_t) );
	ofs_.write( reinterpret_cast<const char*>( &fr ), sizeof(size_t) );
	std::cout << "done." << std::endl;

	ofs_.flush();
#endif // STREAMING
}

SimulationTree::SimulationTree( const std::string& name, 
	ScenePtr scene,
	boost::shared_ptr<SimulationTree> parent, 
	std::deque<ConstPhysicsObjectPtr> activeObjects, 
	const std::vector<PiecewisePath>& piecewisePaths,
	size_t piecewisePathsFrameRate,
	float startTime,
	bool timeReversed )
	:	Named(name), 
		timeRange_( parent->timeRange_ ),
		fixedObjectsHistory_(piecewisePaths),
		fixedObjectsHistoryFramerate_(piecewisePathsFrameRate),
		timeReversed_( timeReversed ),
        numSamples_(0)
{
	this->samplingRate_ = parent->samplingRate_;
	this->simulationRate_ = parent->simulationRate_;
	this->multiplier_ = parent->multiplier_;

	this->scene_ = scene;
	this->parentScene_.reset( new Scene( *parent->scene() ) );

	this->dt_ = parent->dt();
	this->angularVelocityThreshold_ = parent->angularVelocityThreshold_;
	this->linearVelocityThreshold_ = parent->linearVelocityThreshold_;
	this->objectives_ = parent->objectives_;

	this->cameraSource_ = parent->cameraSource_;
	this->minLinearImpulse_ = parent->minLinearImpulse_;
	this->minAngularImpulse_ = parent->minAngularImpulse_;

	/*
	this->pathTrees_.reserve( dynamicObjects_.size() );
	for( size_t i = 0; i < dynamicObjects_.size(); ++i )
		this->pathTrees_.push_back( new PathTree(i) );
		*/

	this->simulationFactory_ = parent->simulationFactory_;
	// first step is to figure out the dynamic object correspondence
	this->dynamicObjects_ = physicsObjects( *this->scene_ );
    this->samplingProperties_.resize( dynamicObjects_.size() );

	{
		ODEDefaultSimulationType parentSim( this->parentScene_ );
		this->parentDynamicObjects_ = parentSim.dynamicObjects();
	}

	// Figure out what has changed between scenes; the symmetric_diff routine
	//   is gauranteed to only take one pass on the scenes so it is significantly
	//   faster than always calling equals() on the objects, plus it captures
	//   relations between objects and joints that we'd otherwise miss
	STLEXT hash_set<const SceneGraphElement*> dirtyObjects;
	std::deque<ConstSceneGraphElementPtr> changed = 
		symmetricDiffScene( *this->scene_, *this->parentScene_ );
	for( std::deque<ConstSceneGraphElementPtr>::const_iterator itr = changed.begin();
		itr != changed.end(); ++itr )
		dirtyObjects.insert( itr->get() );

	{
		for( size_t iDynOb = 0; iDynOb < parentDynamicObjects_.size(); ++iDynOb )
		{
			// If the object is changed from the previous scene, storing the mapping
			//   back doesn't make any sense, really.
			// since the difference is symmetric, it doesn't matter if we look for the
			//   child object or the parent object, but it is cheaper to look for the parent
			//   since we already have a pointer
			if( dirtyObjects.find( parentDynamicObjects_.at(iDynOb).get() ) != dirtyObjects.end() )
				continue;

			std::vector<ConstPhysicsObjectPtr>::iterator myDynItr = 
				std::find_if( this->dynamicObjects_.begin(), this->dynamicObjects_.end(), 
					NameEquals( parentDynamicObjects_.at(iDynOb)->name() ) );
			if( myDynItr != this->dynamicObjects_.end() )
			{
				dynamicObjectToParentDynamicObject_[ 
					std::distance( this->dynamicObjects_.begin(), myDynItr ) ] = iDynOb;
			}
		}
	}

	{
		// go through objects in the old scene and find out if they have no correspondence
		// in the new scene:
		std::deque<ConstSceneGraphElementPtr> toHandle( 1, parentScene_->root() );
		while( !toHandle.empty() )
		{
			ConstSceneGraphElementPtr current = toHandle.front();
			toHandle.pop_front();
			
			ConstPhysicsObjectPtr physicsObject = 
				boost::dynamic_pointer_cast<const PhysicsObject>( current );
			boost::shared_ptr<const CombinedObject> combined = 
				boost::dynamic_pointer_cast<const CombinedObject>( physicsObject );
			if( !combined )
			{
				std::copy( current->children_begin(), current->children_end(),
					std::front_inserter(toHandle) );
			}

			if( !physicsObject )
				continue;

			ConstPhysicsObjectPtr myObject = 
				boost::dynamic_pointer_cast<const PhysicsObject>( this->scene_->object(
					physicsObject->name() ) );
			if( myObject && (dirtyObjects.find(physicsObject.get()) == dirtyObjects.end()) )
			{
				parentObjectToChildObject_[physicsObject.get()] = myObject.get();
			}
			else
			{
				onlyInParentObjects_.push_back( physicsObject );
			}
		}
	}

	{
		// the reverse of the previous:
		// go through objects in the new scene and find out if they have no correspondence
		// in the old scene:
		std::deque<ConstSceneGraphElementPtr> toHandle( 1, scene_->root() );
		while( !toHandle.empty() )
		{
			ConstSceneGraphElementPtr current = toHandle.front();
			toHandle.pop_front();
			
			ConstPhysicsObjectPtr physicsObject = 
				boost::dynamic_pointer_cast<const PhysicsObject>( current );
			boost::shared_ptr<const CombinedObject> combined = 
				boost::dynamic_pointer_cast<const CombinedObject>( physicsObject );
			if( !combined )
			{
				std::copy( current->children_begin(), current->children_end(),
					std::front_inserter(toHandle) );
			}

			if( !physicsObject )
				continue;

			ConstPhysicsObjectPtr parentObject = 
				boost::dynamic_pointer_cast<const PhysicsObject>( parentScene_->object(
					physicsObject->name() ) );
			if( parentObject && (dirtyObjects.find(physicsObject.get()) == dirtyObjects.end()) )
			{
				childObjectToParentObject_[ physicsObject.get() ] =
					parentObject.get();
			}
			else
			{
				onlyInChildObjects_.push_back( physicsObject );
			}
		}
	}

	for( size_t iDynObj = 0; iDynObj < dynamicObjects_.size(); ++iDynObj )
	{
		// objects are automatically active if we don't have a simulation history for them
		if( (std::find_if( activeObjects.begin(), activeObjects.end(), 
					NameEquals(dynamicObjects_[iDynObj]->name()) ) == activeObjects.end())
			&& std::find( onlyInChildObjects_.begin(), onlyInChildObjects_.end(), 
				dynamicObjects_[iDynObj] ) == onlyInChildObjects_.end() )
			fixedObjects_.push_back(iDynObj);
	}

	RandomStream random;

	pathRoots_.resize( dynamicObjects_.size() );
	{
		Simulation::State origState = startState( this->dynamicObjects_ );

		// have to start from the beginning for any new objects
		if( !onlyInChildObjects_.empty() || !onlyInParentObjects_.empty() )
			startTime = 0.0f;

		float samplingRate = boost::numeric_cast<float>( parent->samplingRate_ );
		float scaledTime = startTime * samplingRate;
		std::vector<RigidDynamicState> states;
		for( size_t iDynObj = 0; iDynObj < dynamicObjects_.size(); ++iDynObj )
		{
			DynamicObjectIdMap::const_iterator dynIdItr = 
				dynamicObjectToParentDynamicObject_.find( iDynObj );
			if( dynIdItr == dynamicObjectToParentDynamicObject_.end() )
			{
				// if we don't have state for an object, we'd better be starting
				//   at the beginning:
				assert( startTime == 0.0f );
				states.push_back( origState.state(iDynObj) );
			}
			else
			{
				const PiecewisePath& path = piecewisePaths.at( dynIdItr->second );
				RigidDynamicState state(
					path.position( scaledTime ),
					path.rotation( scaledTime ),
					path.linearVelocity( scaledTime ) * static_cast<float>(fixedObjectsHistoryFramerate_),
					path.angularVelocity( scaledTime ) * static_cast<float>(fixedObjectsHistoryFramerate_) );
				states.push_back( state );

				for( size_t j = path.startFrame(); j < path.endFrame(); j += 4 )
					pathRoots_.at(iDynObj).push_back( path.position(j) );
			}
		}

		Simulation::State state( states, startTime );
		this->startState_ = state;
	}
}

void SimulationTree::setTimeRange( std::pair<float, float> timeRange )
{
	this->timeRange_ = timeRange;
}

void SimulationTree::setSamplingProperties( size_t iDynamicObject, SamplingProperties props )
{
    boost::mutex::scoped_lock lock( samplingPropertiesMutex_ );
    this->samplingProperties_.at( iDynamicObject ) = props;
}

SamplingProperties SimulationTree::getSamplingProperties( size_t iDynamicObject ) const
{
    boost::mutex::scoped_lock lock( samplingPropertiesMutex_ );
    return this->samplingProperties_.at( iDynamicObject );
}

size_t SimulationTree::idForPath( const Path* path ) const
{
	PathToIdMap::const_iterator itr = pathToId_.find( path );
	assert( itr != pathToId_.end() );
	return itr->second;
}

bool SimulationTree::satisfiesConstraint( size_t pathId, size_t constraintId ) const
{
	return (this->constraintResults_[constraintId]).test(pathId);
}

size_t SimulationTree::constraintCount() const
{
	return this->constraints_.size();
}

bool SimulationTree::checkForNewPaths()
{
#ifdef DUMP_PATHS
	static std::ofstream pathsDumpOFS( "dumpedPaths.bin", std::ofstream::binary );
#endif
	const size_t maxNodesToAdd = 50;
	size_t nodesAdded = 0;
	std::deque<const Path*> newPaths;
	boost::dynamic_bitset<> active;

	while( nodesAdded < maxNodesToAdd )
	{
		boost::mutex::scoped_lock lock( this->nodesToAddMutex_ );
		if( !pathsToAdd_.empty() )
		{
			const Path* path = pathsToAdd_.front();
			pathsToAdd_.pop_front();
			++nodesAdded;

#ifdef DUMP_PATHS
			// this is not threadsafe but only gets called from the main
			//   thread so we don't need to worry.
			{
				std::deque< std::vector<vl::Vec3f> > paths;
				for( size_t iTree = 0; iTree < path->treeCount(); ++iTree )
				{
					std::vector<vl::Vec3f> points;
					path->obbTree( iTree ).points( points );
					if( points.size() > 1 )
						paths.push_back( points );
				}

				int numPaths = paths.size();
				pathsDumpOFS.write( reinterpret_cast<const char*>(&numPaths), sizeof(int) );
				for( size_t iPath = 0; iPath < paths.size(); ++iPath )
				{
					const std::vector<vl::Vec3f>& path = paths.at( iPath );
					dumpArray( path, pathsDumpOFS );
				}

				int numSnaps = path->snapshotCount();
				pathsDumpOFS.write( reinterpret_cast<const char*>(&numSnaps), sizeof(int) );
				for( int iSnap = 0; iSnap < numSnaps; ++iSnap )
				{
					float snapshotTime = path->snapshotTime(iSnap);
					std::vector<RigidStaticState> snapshot = path->snapshot( iSnap );
					pathsDumpOFS.write( reinterpret_cast<const char*>( &snapshotTime ), sizeof(float) );
					dumpArray( snapshot, pathsDumpOFS );
				}

				pathsDumpOFS.flush();
			}
#endif

			{
				ReaderWriterLock::ScopedWriterLock lock( this->allPathsLock_ );
				pathToId_.insert( PathToIdMap::value_type(path, allPaths_.size()) );
				allPaths_.push_back( path );
			}

			ReaderWriterLock::ScopedReaderLock lock( this->constraintsLock_ );
			bool nodeOkay = true;
			for( size_t iConstraint = 0; iConstraint < constraints_.size(); ++iConstraint )
			{
				bool constraintResult = constraints_[iConstraint]->satisfies( *path );
				constraintResults_[iConstraint].push_back( constraintResult );
				nodeOkay = nodeOkay && constraintResult;
			}

			/*
			{
				ReaderWriterLock::ScopedWriterLock pathTreesLock( this->pathTreesLock_ );
				assert( pathTrees_.size() == this->dynamicObjects_.size() );
				boost::mutex::scoped_lock lock( this->pathTreesMutex_ );
				for( size_t i = 0; i < pathTrees_.size(); ++i )
					pathTrees_.at(i).add( path );
			}
			*/

			// we need to add it immediately; if there were paths "in flight"
			// this would tend to screw up the constraint handling
			{
				ReaderWriterLock::ScopedWriterLock lock( this->activePathsLock_ );
				activePaths_.push_back( nodeOkay );

				newPaths.push_back( path );
				active.push_back( nodeOkay );
			}
		}
		else
		{
			break;
		}
	}

	if( !newPaths.empty() )
	{
		boost::mutex::scoped_lock lock( this->callbackMutex_ );
		std::for_each( this->callbacks_.begin(),
			this->callbacks_.end(),
			boost::bind( &SimulationNodeCallback::addNodes, _1, this, 
				newPaths, active ) );
		return true;
	}
	else
	{
		return false;
	}
}

boost::dynamic_bitset<> SimulationTree::addConstraint( ConstraintPtr constraint,
								   const boost::dynamic_bitset<>& hint )
{
	ReaderWriterLock::ScopedWriterLock constraintLock( this->constraintsLock_ );
	this->constraints_.push_back( constraint );
	this->constraintResults_.push_back( hint );

	this->constraintResults_.back().resize( this->activePaths_.size() );

	for( size_t iPath = hint.size(); iPath < constraintResults_.back().size(); ++iPath )
		(constraintResults_.back())[iPath] = constraint->satisfies( *this->path(iPath) );

	this->activePaths_ &= constraintResults_.back();

	{
		boost::mutex::scoped_lock lock( this->callbackMutex_ );
		std::for_each( this->callbacks_.begin(),
			this->callbacks_.end(),
			boost::bind( &SimulationNodeCallback::setNodes, _1, this, 
				this->activePaths_ ) );
	}

	return constraintResults_.back();
}

void SimulationTree::removeConstraint( ConstraintPtr constraint )
{
	ReaderWriterLock::ScopedWriterLock constraintLock( this->constraintsLock_ );

	// remove constraint
	{
		std::deque<ConstraintPtr>::iterator constraintItr = 
			std::find( constraints_.begin(), constraints_.end(), constraint );
		size_t iConstraint = std::distance( constraints_.begin(), constraintItr );

		this->constraintResults_.erase( constraintResults_.begin() + iConstraint );
		this->constraints_.erase( constraintItr );
	}

	boost::dynamic_bitset<> newActive( allPaths_.size() );
	newActive.set();
	for( size_t iConstraint = 0; iConstraint < constraints_.size(); ++iConstraint )
		newActive &= constraintResults_[iConstraint];

	activePaths_.swap( newActive );

	{
		boost::mutex::scoped_lock lock( this->callbackMutex_ );
		std::for_each( this->callbacks_.begin(),
			this->callbacks_.end(),
			boost::bind( &SimulationNodeCallback::setNodes, _1, this, 
				this->activePaths_ ) );
	}
}

std::vector<ConstPhysicsObjectPtr> SimulationTree::dynamicObjects() const
{
	return dynamicObjects_;
}

boost::shared_ptr<const Scene> SimulationTree::scene() const
{
	return scene_;
}

boost::shared_ptr<const Scene> SimulationTree::parentScene() const
{
	return this->parentScene_;
}

SimulationTree::~SimulationTree()
{
}

vl::Vec3f SimulationTree::randomImpulse(RandomStream& random)
{
	std::vector<double> impulse( random.uniformSpherical(3) );
	return vl::Vec3f( impulse[0], impulse[1], impulse[2] );
}

// Finds an arbitrary basis, with N in the Z direction
vl::Mat3f computeBasis( const vl::Vec3f& N )
{
	vl::Mat3f basis;
	basis[2] = vl::norm(N);

	vl::Vec3f test1 = cross( vl::Vec3f( 0.0, 1.0, 0.0 ), basis[2] );
	vl::Vec3f test2 = cross( vl::Vec3f( 1.0, 0.0, 0.0 ), basis[2] );
	if( sqrlen( test2 ) > sqrlen( test1 ) )
		basis[1] = vl::norm( test2 );
	else
		basis[1] = vl::norm( test1 );

	basis[0] = norm( cross( basis[1], basis[2] ) );
	basis[1] = norm( cross( basis[2], basis[0] ) );

	return basis;
}

vl::Vec3f SimulationTree::randomImpulse(const vl::Vec3f& source, RandomStream& random)
{
	return randomImpulse( source, 
		random.uniform(), 
		random.uniform(), 
		random.uniform() );
}

vl::Vec3f SimulationTree::randomImpulse(const vl::Vec3f& source, double xi_1, double xi_2, double nu)
{
	if( sqrlen(source) < 1e-5 ) 
		return vl::Vec3f( 0.0, 0.0, 0.0 );

	vl::Mat3f basis = computeBasis( source );

	// this is the n from the cosine lobe:
	const double n = 5.0;

	double shared_radical = sqrt( 1 - pow(xi_1, 2.0/(n+1.0)) );
	vl::Vec3f d( 
		shared_radical * cos( 2*M_PI*xi_2 ),
		shared_radical * sin( 2*M_PI*xi_2 ),
		pow( xi_1, 1.0/(n+1.0) ) );
	d = trans(basis) * vl::norm(d);

	//double lenSource = len(source);
	//double cosTheta = lenSource*pow( static_cast<double>(dot( d, source )) / lenSource, n );
	//double cosTheta = dot( d, source );

	double length = (-0.01 + 0.02*nu) * len(source);
	assert( length <= len(source) );

	vl::Vec3f result = length * d;
	assert( !isBad(result) );
	assert( len(result) < len(source) );
	return result;
}

vl::Vec3f randomDirection(const vl::Vec3f& source, double xi_1, double xi_2, double n)
{
	if( sqrlen(source) < 1e-5 ) 
		return vl::Vec3f( 0.0, 0.0, 0.0 );

	vl::Mat3f basis = computeBasis( source );

	double shared_radical = sqrt( 1 - pow(xi_1, 2.0/(n+1.0)) );
	vl::Vec3f d( 
		shared_radical * cos( 2*M_PI*xi_2 ),
		shared_radical * sin( 2*M_PI*xi_2 ),
		pow( xi_1, 1.0/(n+1.0) ) );
	d = trans(basis) * vl::norm(d);

	return d;
}

class ImpulseWrapper
{
public:
	ImpulseWrapper( float minImpulse )
		: minImpulse_(minImpulse) {}

	float minImpulse() const { return minImpulse_; }

	virtual void perturb( Simulation& simulation, size_t iDynOb, const vl::Vec3f& perturbation ) = 0;
	virtual vl::Vec3f get( const RigidDynamicState& state ) const = 0;
	virtual vl::Vec3f get( const Simulation& simulation, size_t iDynOb ) const = 0;

    virtual bool linear() const = 0;

private:
	float minImpulse_;
};

class AngularImpulseWrapper
	: public ImpulseWrapper
{
public:
	AngularImpulseWrapper( float minImpulse )
		: ImpulseWrapper( minImpulse ) {}

	void perturb( Simulation& simulation, size_t iDynOb, const vl::Vec3f& perturbation )
	{
		simulation.perturbAngularVelocity( iDynOb, perturbation );
	}

	vl::Vec3f get( const RigidDynamicState& state ) const
	{
		return toVec3f(state.angularVelocity());
	}

	vl::Vec3f get( const Simulation& simulation, size_t iDynOb ) const
	{
		return simulation.getAngularVelocity( iDynOb );
	}

    bool linear() const { return false; }
};

class LinearImpulseWrapper
	: public ImpulseWrapper
{
public:
	LinearImpulseWrapper( float minImpulse )
		: ImpulseWrapper( minImpulse ) {}

	void perturb( Simulation& simulation, size_t iDynOb, const vl::Vec3f& perturbation )
	{
		simulation.perturbLinearVelocity( iDynOb, perturbation );
	}

	vl::Vec3f get( const RigidDynamicState& state ) const
	{
		return toVec3f(state.linearVelocity());
	}

	vl::Vec3f get( const Simulation& simulation, size_t iDynOb ) const
	{
		return simulation.getLinearVelocity( iDynOb );
	}

    bool linear() const { return true; }
};


template <typename T>
class FirstEquals
{
public:
	FirstEquals( const typename T::first_type& val )
		: value_(val) {}

	bool operator()( const T& rhs )
	{
		return rhs.first == value_;
	}

private:
	typename T::first_type value_;
};

// We will want to keep track of simulations, their scores, their collision
// events, etc. 
struct ScoredSimulation
{
    ScoredSimulation( size_t numDynamicObjects )
        :   states( numDynamicObjects ), 
            endState( numDynamicObjects ), 
            score( boost::numeric::bounds<float>::highest() ) {}

    std::vector< std::deque<RigidDynamicState> > states;
    std::vector<RigidDynamicState> endState;
    float score;

    void swap( ScoredSimulation& other )
    {
        std::swap( score, other.score );
        states.swap( other.states );
        endState.swap( other.endState );
    }
};

float dist( const RigidDynamicState& state1, const RigidDynamicState& state2 )
{
    float res = vl::len(state1.position() - state2.position());
    res += vl::len(state1.linearVelocity() - state2.linearVelocity());
    res += vl::len(state1.angularVelocity() - state2.angularVelocity());

    {
        vl::Mat3f rot1 = state1.orientation().toRotMatf();
        vl::Mat3f rot2 = state2.orientation().toRotMatf();
        float dist = 0.0f;
        for( vl::Int i = 0; i < 3; ++i )
        {
            for( vl::Int j = 0; j < 3; ++j )
            {
                float diff = rot1[i][j] - rot2[i][j];
                dist += diff*diff;
            }
        }
        res += sqrt(dist);
    }

    return res;
}

vl::Vec3f randomVector( RandomStream& random, float maxLength )
{
    std::vector<double> v = random.uniformSpherical(3);
    vl::Vec3f result;
    std::copy( v.begin(), v.end(), result.Ref() );
    result *= random.uniform( 0.0, maxLength );
    return result;
}

typedef std::deque<RigidDynamicState> ObjectStateList;

// We'd like to use forward motion instead of backward motion here, to try to
// avoid some of ths problems (energy growth, etc.) that are associated with
// backward simulation.  So we'll back the simulation up n timesteps and run
// a whole bunch of forward simulations to try to find ones that have the 
// characteristics we want.  

class IntegerUnionFind
{
public:
    IntegerUnionFind( size_t n )
        : rank_(n, 0)
    {
        rank_.reserve( n );
        for( size_t i = 0; i < n; ++i )
            parent_.push_back( i );
    }

    size_t findElement( size_t element )
    {
        size_t parent = parent_[element];
        if( parent == element )
            return element;
        else
        {
            parent = findElement( parent );
            parent_[element] = parent;
            return parent;
        }
    }

    void unionElements( size_t x, size_t y )
    {
        size_t xRoot = findElement( x );
        size_t yRoot = findElement( y );

        if( rank_[xRoot] > rank_[yRoot] )
            parent_[yRoot] = xRoot;
        else if( rank_[xRoot] < rank_[yRoot] )
            parent_[xRoot] = yRoot;
        else if( xRoot != yRoot )
        {
            parent_[yRoot] = xRoot;
            rank_[xRoot]++;
        }
    }

private:
    std::vector<size_t> rank_;
    std::vector<size_t> parent_;
};

// @todo: for now we won't deal with activating objects, we'll just
// assume they don't interact
std::vector<Simulation::State> SimulationTree::sampleForwardToState( 
    SimulationPtr simulation,
    size_t numTimesteps,
    RandomStream& random )
{
	boost::ptr_vector<ImpulseWrapper> impulseWrappers;
	impulseWrappers.push_back( new AngularImpulseWrapper(minAngularImpulse_) );
	impulseWrappers.push_back( new LinearImpulseWrapper(minLinearImpulse_) );

    // todo: this should eventually back the state up to machine precision
    std::vector<Simulation::State> backupState;
    backupState.reserve( numTimesteps*multiplier_ );

    // @todo this is going to be super-slow, but oh well
    // fix it some other time
    std::vector<ConstPhysicsObjectPtr> dynamicObjects = 
        simulation->dynamicObjects();

    // This is a union-find data structure for keep track of who is connected together:
    IntegerUnionFind connected( dynamicObjects.size() );

    for( size_t i = 0; i < dynamicObjects.size(); ++i )
        for( size_t j = i+1; j < dynamicObjects.size(); ++j )
            if( simulation->connected( dynamicObjects[i].get(), dynamicObjects[j].get() ) )
                connected.unionElements( i, j );

    for( size_t iObject = 0; iObject < dynamicObjects_.size(); ++iObject )
        simulation->setDisabled( iObject, false );

    for( size_t i = 0; i < numTimesteps*multiplier_; ++i )
    {
        simulation->step( -1.0 / boost::numeric_cast<double>(simulationRate_) );
        backupState.push_back( simulation->getSimulationState() );

#ifdef _DEBUG
        if( backupState.size() > 2 )
        {
            RigidDynamicState prev = backupState[ backupState.size() - 2 ].state(0);
            RigidDynamicState next = backupState[ backupState.size() - 1 ].state(0);
            assert( vl::len(prev.linearVelocity()) < 1e-2 ||
                prev.linearVelocity() != next.linearVelocity() ); 
        }
#endif
    }
#ifdef _DEBUG
    {
        RigidDynamicState last = backupState.back().state(0);
        RigidDynamicState first = backupState.front().state(0);

        assert( last.linearVelocity() != first.linearVelocity() );
    }
#endif

    std::reverse( backupState.begin(), backupState.end() );

    std::vector<RigidDynamicState> goalState; goalState.reserve( dynamicObjects_.size() );
    for( size_t iObject = 0; iObject < dynamicObjects_.size(); ++iObject )
        goalState.push_back( backupState.back().state( iObject ) );

    std::vector<ObjectStateList> result( dynamicObjects_.size() );
    
    // Now, start running the sim
    for( size_t curObject = 0; curObject < dynamicObjects_.size(); ++curObject )
    {
startOfLoop:
        // this means we've already dealt with this object:
        if( !result.at(curObject).empty() )
            continue;

        // @todo this is also going to be slow
        std::vector<size_t> currentObjects;
        for( size_t otherObject = 0; otherObject < dynamicObjects_.size(); ++otherObject )
            if( connected.findElement(otherObject) == connected.findElement(curObject) )
                currentObjects.push_back( otherObject );

        std::vector< std::deque<vl::Vec3f> > linearImpulses( dynamicObjects_.size() );
        std::vector< std::deque<vl::Vec3f> > angularImpulses( dynamicObjects_.size() );

        ScoredSimulation currentSim( dynamicObjects_.size() );
        ScoredSimulation bestSim( dynamicObjects_.size() );

        for( size_t iIter = 0, refinementIter = 0; iIter < 500; ++iIter )
        {
            ScoredSimulation newSim( dynamicObjects_.size() );

            if( currentSim.states.front().empty() || refinementIter > 5 )
            {
                for( std::vector<size_t>::const_iterator obIter = currentObjects.begin();
                    obIter != currentObjects.end(); ++obIter )
                {
                    linearImpulses[*obIter].clear();
                    angularImpulses[*obIter].clear();

                    newSim.states[*obIter].push_back( backupState.front().state(*obIter) );

                    // Randomize the state a little
                    newSim.states[*obIter].back().setPosition( newSim.states[*obIter].back().position() + randomVector(random, 0.1) );
                    newSim.states[*obIter].back().setOrientation( Quaternion<float>( randomVector(random, 0.1) ) 
                        * newSim.states[*obIter].back().orientation() );
                    newSim.states[*obIter].back().setLinearVelocity( newSim.states[*obIter].back().linearVelocity() + randomVector(random, 0.5) );
                    newSim.states[*obIter].back().setAngularVelocity( newSim.states[*obIter].back().angularVelocity() + randomVector(random, 0.2) );
                }

                refinementIter = 0;
            }
            else
            {
                for( std::vector<size_t>::const_iterator obIter = currentObjects.begin();
                    obIter != currentObjects.end(); ++obIter )
                {
                    assert( !currentSim.states[*obIter].empty() );
                    newSim.states[*obIter].push_back( currentSim.states[*obIter].front() );
                }

                ++refinementIter;
            }

            // wake up any objects that were asleep in previous sims
            for( size_t iObject = 0; iObject < dynamicObjects_.size(); ++iObject )
                simulation->setDisabled( iObject, false );

            // Initialize just the newly-simulated objects:
            simulation->setTime( backupState.front().time() );
            for( std::vector<size_t>::const_iterator obIter = currentObjects.begin();
                obIter != currentObjects.end(); ++obIter )
            {
                simulation->setSimulationState( *obIter, newSim.states[*obIter].front() );
            }

            for( size_t iStep = 0; iStep + 1 < backupState.size(); ++iStep )
            {
                // Update the state of all the non-simulated objects
                for( size_t iObject = 0; iObject < dynamicObjects_.size(); ++iObject )
                {
                    if( std::binary_search(currentObjects.begin(), currentObjects.end(), iObject) )
                        continue;

                    if( result[iObject].empty() )
                        simulation->setSimulationState( iObject, backupState[iStep].state(iObject) );
                    else
                        simulation->setSimulationState( iObject, result[iObject][iStep] );
                }

                assert(simulation->time() >= backupState.front().time() - 1e-4 &&
                       simulation->time() <= backupState.back().time() + 1e-4 );
                simulation->step( 1.0 / boost::numeric_cast<double>(simulationRate_) );

                // check all interactions
                std::vector<Simulation::ContactPoint> contacts = simulation->contactPoints();
                for( std::vector<Simulation::ContactPoint>::const_iterator iter = contacts.begin();
                    iter != contacts.end(); ++iter )
                {
                    int first = iter->objects.first, second = iter->objects.second;
                    if( first < 0 || second < 0 )
                        continue;

                    bool firstInSet = std::binary_search( currentObjects.begin(), currentObjects.end(), 
                                                          first );
                    bool secondInSet = std::binary_search( currentObjects.begin(), currentObjects.end(), 
                                                          second );
                    if( !firstInSet && secondInSet )
                    {
                        std::swap( first, second );
                        std::swap( firstInSet, secondInSet );
                    }

                    if( firstInSet && !secondInSet )
                    {
                        connected.unionElements( first, second );

                        // need to wipe the simulations for all the connected objects
                        for( size_t iObject = 0; iObject < dynamicObjects_.size(); ++iObject )
                        {
                            if( connected.findElement(iObject) == connected.findElement(first) )
                                result[iObject].clear();
                        }

                        // now go back to start
                        goto startOfLoop; 
                    }
                }

                for( std::vector<size_t>::const_iterator obIter = currentObjects.begin();
                    obIter != currentObjects.end(); ++obIter )
                {
			        const RigidDynamicState& prevObjectState = newSim.states[*obIter].back();

                    assert( linearImpulses.size() == angularImpulses.size() );

                    if( linearImpulses[*obIter].size() < newSim.states[*obIter].size() )
                    {
			            for( boost::ptr_vector<ImpulseWrapper>::iterator wrapperItr = impulseWrappers.begin();
				            wrapperItr != impulseWrappers.end(); ++wrapperItr )
			            {
				            vl::Vec3f prevValue = wrapperItr->get( prevObjectState );
				            vl::Vec3f nextValue = wrapperItr->get( *simulation, *obIter );
				            assert( !isBad( prevValue ) );
				            assert( !isBad( nextValue ) );

				            vl::Vec3f origImpulse = nextValue - prevValue;
				            assert( !isBad(origImpulse) );

                            vl::Vec3f perturbation( vl::vl_0 );
				            if( vl::len(origImpulse) > wrapperItr->minImpulse() )
                            {
				                perturbation = randomImpulse( origImpulse, random );
    				            wrapperItr->perturb( *simulation, *obIter, perturbation );
                            }

                            if( wrapperItr->linear() )
                                linearImpulses[*obIter].push_back( perturbation );
                            else
                                angularImpulses[*obIter].push_back( perturbation );

			            }
                    }
                    else
                    {
                        // apply the force to the simulation
                        simulation->perturbLinearVelocity( *obIter, linearImpulses[*obIter][newSim.states.size()-1] );
                        simulation->perturbAngularVelocity( *obIter, angularImpulses[*obIter][newSim.states.size()-1] );
                    }

                    newSim.states[*obIter].push_back( simulation->getSimulationState(*obIter) );
                }
            }

            // for visualization purposes.
            /*
            {
                std::deque<Simulation::State> fullStates;
                assert( backupState.size() == newSim.states.size() );
                for( size_t i = 0; i < backupState.size(); i += multiplier_ )
                {
                    fullStates.push_back( backupState[i] );
                    fullStates.back().state(iObject) = newSim.states[i];
                }

                std::vector<CompressedPath> compressedPaths( this->dynamicObjects_.size() );
    	        float quadraticCoeff = -0.5 * 9.81 * boost::numeric_cast<double>(this->dt_) * boost::numeric_cast<double>(this->dt_);
                compressedPaths.at( iObject ) = 
                    compress( newSim.states, quadraticCoeff, fullStates.front().time(), this->samplingRate_*this->multiplier_, 1e-4, 1e-4 );

                PiecewisePath compressedPath( compressedPaths[iObject] );
                assert( vl::len( compressedPath.position( compressedPath.startFrame() ) - 
                    compressedPath.position( compressedPath.endFrame() ) ) > 1e-1 );

		        float epsilon = 1e-5;

	            std::vector<float> objectives( dynamicObjects_.size() * objectives_.size(), 1000 + iIter );

	            std::vector< std::vector<vl::Vec3f> > cmPaths( dynamicObjects_.size() );
		        cmPaths[iObject].reserve( newSim.states.size() );
		        for( size_t i = 0; i < newSim.states.size(); ++i )
			        cmPaths[iObject].push_back( toVec3f(newSim.states[i].position()) );

                Path* path;
                {
		            boost::mutex::scoped_lock lock( this->pathPoolMutex_ );
		            path = new (pathPool_.malloc())
			            Path( compressedPaths, fullStates, cmPaths, this->samplingRate_, objectives, epsilon );
                }

		        boost::mutex::scoped_lock lock( this->nodesToAddMutex_ );
		        pathsToAdd_.push_back( path );
	        }
            */

            newSim.score = 0.0;
            for( std::vector<size_t>::const_iterator obIter = currentObjects.begin();
                obIter != currentObjects.end(); ++obIter )
            {
                newSim.endState[*obIter] = newSim.states[*obIter].back();
                newSim.score += 
                    angleBetween( newSim.endState[*obIter].orientation(), goalState[*obIter].orientation() ) +
                    vl::len( newSim.endState[*obIter].position() - goalState[*obIter].position() ) +
                    vl::len( newSim.endState[*obIter].linearVelocity() - goalState[*obIter].linearVelocity() ) +
                    vl::len( newSim.endState[*obIter].angularVelocity() - goalState[*obIter].angularVelocity() );

                // need to warp the one motion into the other
                if( numTimesteps > 1 )
                {
                    float delta_t = 1.0f / static_cast<float>(simulationRate_);

                    size_t maxBlendFrames = this->simulationRate_*1; // 1 seconds maximum
                    size_t blendFrames = std::min( maxBlendFrames, newSim.states[*obIter].size()-1 );
                    typedef generic_matrix< float, boost::numeric::ublas::column_major > float_matrix;
                    float_matrix T( 2, blendFrames );

                    boost::scoped_array<float> weights( new float[T.ncols()] );
                    size_t startBlendFrame = newSim.states[*obIter].size() - blendFrames - 1;
                    for( size_t i = 0; i < blendFrames; ++i )
                    {
                        const float epsilon = 1e-2;
                        vl::Vec3f impulse = newSim.states[*obIter][startBlendFrame+i+1].linearVelocity() - 
                            newSim.states[*obIter][startBlendFrame+i].linearVelocity();
                        weights[i] = sqrt(epsilon + vl::len(impulse));
                    }

                    for( size_t i = 0; i < blendFrames; ++i )
                    {
                        T( 0, i ) = delta_t * weights[i];
                        T( 1, i ) = (delta_t * (backupState.back().time() - 
                            0.5f*(backupState[startBlendFrame + i + 1].time() + backupState[startBlendFrame+i].time()))) * weights[i];
                    }

                    int m = T.nrows();
                    int n = T.ncols();
                    int lda = T.nrows();
                    char trans = 'N';
                    int nrhs = 3;

                    float_matrix rhs( n, 3 );
                    vl::Vec3f posDiff = backupState.back().state(*obIter).position() - 
                        newSim.endState[*obIter].position();
                    vl::Vec3f velDiff = backupState.back().state(*obIter).linearVelocity() - 
                        newSim.endState[*obIter].linearVelocity(); 
                    for( size_t i = 0; i < 3; ++i )
                    {
                        rhs( 0, i ) = velDiff[i];
                        rhs( 1, i ) = posDiff[i];
                    }
                    int ldb = rhs.nrows();

                    // todo apply weights

                    int blocksize = 64;
                    int lwork = std::min( m, n ) + std::max( std::max( std::max( 1, m ), n), nrhs ) * blocksize;

                    int info = 0;

                    boost::scoped_array<float> work( new float[lwork] );

                    sgels( &trans,
                        &m,
                        &n,
                        &nrhs,
                        T.data(),
                        &lda,
                        rhs.data(),
                        &ldb,
                        &work[0],
                        &lwork,
                        &info );

                    vl::Vec3f posAccum = vl::vl_0;
                    vl::Vec3f velAccum = vl::vl_0;
                    assert( newSim.states[*obIter].size() == numTimesteps*multiplier_ );
                    linearImpulses[*obIter].resize( newSim.states[*obIter].size() - 1, vl::vl_0 );
                    for( size_t i = 1; i <= blendFrames; ++i )
                    {
                        vl::Vec3f accel( rhs(i-1, 0), rhs(i-1, 1), rhs(i-1, 2) );
                        assert( !isBad(accel) );

                        accel *= weights[i-1];

                        posAccum += delta_t*velAccum + 0.5f*delta_t*delta_t*accel;
                        velAccum += delta_t*accel;
                        linearImpulses[*obIter][ startBlendFrame+i-1 ] += delta_t*accel;

                        newSim.states[*obIter][startBlendFrame+i].setPosition( newSim.states[*obIter][startBlendFrame+i].position() + posAccum );
                        newSim.states[*obIter][startBlendFrame+i].setLinearVelocity( newSim.states[*obIter][startBlendFrame+i].linearVelocity() + velAccum );
                    }

                    vl::Vec3f posDiffNew = newSim.states[*obIter].back().position() - goalState[*obIter].position();
                    vl::Vec3f velDiffNew = newSim.states[*obIter].back().linearVelocity() - goalState[*obIter].linearVelocity();
                    //assert( vl::len(posDiffNew) < 1e-3 && vl::len(velDiffNew) < 1e-3 );

                    // Need to fix rotations, too
                    angularImpulses[*obIter].resize( newSim.states[*obIter].size() - 1, vl::vl_0 );

                    vl::Vec3f prevAngularVel;
                    const size_t rotBlendFrames = this->simulationRate_ / 4;  // half a second total (1/4s on either side)
                    const size_t maxFrame = this->simulationRate_ * 1;     // 1 seconds
                    const double velocityWeighting = 1.0;

                    assert( newSim.states[*obIter].size() > 2*rotBlendFrames );
                    assert( maxFrame > 2*rotBlendFrames );
                    size_t iBestFrame = 0;
                    double bestScore = boost::numeric::bounds<double>::highest();
                    for( size_t iFrame = (newSim.states[*obIter].size() > maxFrame) ? 
                                (newSim.states[*obIter].size() - maxFrame + rotBlendFrames + 1) : (rotBlendFrames + 1); 
                         iFrame < (newSim.states[*obIter].size() - rotBlendFrames);
                         ++iFrame )
                    {
                        const vl::Mat3d newSimRot  = newSim.states[*obIter][iFrame].orientation().toRotMatd();
                        const vl::Mat3d origSimRot = backupState[iFrame].state(*obIter).orientation().toRotMatd();

                        const vl::Vec3f newSimAngVel = newSim.states[*obIter][iFrame].angularVelocity();
                        const vl::Vec3f origAngVel = backupState[iFrame].state(*obIter).angularVelocity();

                        double error = 0.0;
                        for( vl::Int i = 0; i < 3; ++i )
                            for( vl::Int j = 0; j < 3; ++j )
                                error += newSimRot[i][j] - origSimRot[i][j];

                        error += velocityWeighting*vl::len(newSimAngVel - origAngVel);

                        if( error < bestScore )
                        {
                            iBestFrame = iFrame;
                            bestScore = error;
                        }
                    }

                    vl::Vec3d curVel = toVec3d(newSim.states[*obIter][ iBestFrame - rotBlendFrames ].angularVelocity());
                    const Quaternion<double> startRot = newSim.states[*obIter][ iBestFrame - rotBlendFrames ].orientation();
                    Quaternion<double> endRot = backupState[ iBestFrame + rotBlendFrames ].state(*obIter).orientation();
                    if( startRot.dot( endRot ) < 0 )
                        endRot *= -1.0;
                    Quaternion<double> prevRot = startRot;

                    assert( iBestFrame > rotBlendFrames && iBestFrame < newSim.states[*obIter].size() - rotBlendFrames );
                    for( size_t iBlendFrame = iBestFrame - rotBlendFrames + 1; iBlendFrame <= iBestFrame + rotBlendFrames; ++iBlendFrame )
                    {
                        double weight = static_cast<double>( iBlendFrame + rotBlendFrames - iBestFrame ) / static_cast<double>( 2*rotBlendFrames );
                        assert( weight >= 0.0 && weight <= 1.0 );

                        Quaternion<double> newRot = slerp( startRot, endRot, weight );

                        vl::Vec3d desiredAngularVel = differentialRotation( prevRot, newRot ) / delta_t;
                        vl::Vec3d delta = desiredAngularVel - curVel;

                        if( refinementIter >= 2 )
                            angularImpulses[*obIter][iBlendFrame - 1] += toVec3f(delta);

                        curVel = desiredAngularVel;
                        prevRot = newRot;

                        newSim.states[*obIter][iBlendFrame].setOrientation( newRot );
                        newSim.states[*obIter][iBlendFrame - 1].setAngularVelocity( toVec3f(desiredAngularVel) );
                    }

                    for( size_t iFrame = iBestFrame + rotBlendFrames + 1; iFrame < newSim.states[*obIter].size(); ++iFrame )
                    {
                        newSim.states[*obIter][iFrame].setOrientation( backupState[ iFrame ].state(*obIter).orientation() );
                        newSim.states[*obIter][iFrame].setAngularVelocity( backupState[ iFrame ].state(*obIter).angularVelocity() );
                    }
                }

                assert( !newSim.states.empty() );
#ifdef _DEBUG
                {
                    RigidDynamicState front = newSim.states[*obIter].front();
                    RigidDynamicState back = newSim.states[*obIter].back();
                    assert( front.position() != back.position() );
                    assert( front.linearVelocity() != back.linearVelocity() );
                    assert( front.orientation() != back.orientation() );
                }
#endif
            }

            if( newSim.score < bestSim.score )
                bestSim = newSim;

            currentSim.swap( newSim );
        }

        bool bad = false;
        for( std::vector<size_t>::const_iterator obIter = currentObjects.begin();
            obIter != currentObjects.end(); ++obIter )
        {
            if(    vl::len(bestSim.endState[*obIter].position() - goalState[*obIter].position()) > 0.5 
                || angleBetween(bestSim.endState[*obIter].orientation(), goalState[*obIter].orientation()) > 0.5
                /*
                || vl::len(bestSim.endState[*obIter].linearVelocity() - goalState[*obIter].linearVelocity()) > 2.0
                || vl::len(bestSim.endState[*obIter].angularVelocity() - goalState[*obIter].angularVelocity()) > 2.0 */ )
                bad = true;
        }

        if( !bad )
        {
            // sim is just too far away!
            for( std::vector<size_t>::const_iterator obIter = currentObjects.begin();
                obIter != currentObjects.end(); ++obIter )
            {
                result[*obIter] = bestSim.states[*obIter];
            }
        }
    }

    for( size_t iObject = 0; iObject < dynamicObjects_.size(); ++iObject )
    {
        if( result[iObject].empty() )
            continue;

        assert( backupState.size() == result[iObject].size() );
        for( size_t i = 0; i < backupState.size(); ++i )
            backupState[i].state(iObject) = result[iObject][i];
    }

    // We simulated at twice the rate of playback, so strip the extras out:
    assert( (backupState.size() % multiplier_) == 0 );
    for( size_t i = 0; i < backupState.size(); i += multiplier_ )
        backupState[i / multiplier_] = backupState[i];
    backupState.resize(backupState.size() / multiplier_);

    simulation->setSimulationState( backupState.front() );

    return backupState;
}

Path* SimulationTree::sample( 
	size_t maxTimesteps,
	size_t minTimesteps,
	RandomStream& random )
{
	assert( maxTimesteps > 0 );

    // Need to randomize all the random properties in the scene
	// make a local copy of the scene before we change anything
	ScenePtr localScene( new Scene(*scene_) );
	for( Scene::object_iter itr = localScene->begin_objects();
		itr != localScene->end_objects(); ++itr )
	{
		std::vector<RandomizedProperty> props = 
			(*itr)->randomizedProperties();
		for( std::vector<RandomizedProperty>::const_iterator propItr = props.begin();
			propItr != props.end(); ++propItr )
		{
			assert( propItr->low <= propItr->high );
			if( propItr->low == propItr->high )
			{
				(*itr)->setAttribute( propItr->name, propItr->low );
			}
			else
			{
				double val = random.uniform( 
					boost::numeric_cast<double>( propItr->low ),
					boost::numeric_cast<double>( propItr->high ) );
				(*itr)->setAttribute( propItr->name, val );
			}
		}
	}

	ObjectMap childObjectToParentObjectLocal;
	for( ObjectMap::const_iterator itr = this->childObjectToParentObject_.begin();
		itr != childObjectToParentObject_.end(); ++itr )
	{
		const PhysicsObject* childObject = itr->first;
		ConstSceneGraphElementPtr localObject = localScene->object( childObject->name() );
		const PhysicsObject* localPhysicsObject = dynamic_cast<const PhysicsObject*>( localObject.get() );
		assert( localPhysicsObject != 0 );
		childObjectToParentObjectLocal.insert( ObjectMap::value_type(localPhysicsObject, itr->second) );
	}

	SimulationPtr simulation( simulationFactory_->construct(localScene) );

    {
        for( size_t iDynOb = 0; iDynOb < this->samplingProperties_.size(); ++iDynOb )
        {
            SamplingProperties props;
            {
                boost::mutex::scoped_lock lock( samplingPropertiesMutex_ );
                props = this->samplingProperties_.at( iDynOb );
            }

            if( props.frictionDirection != vl::vl_0 )
                props.frictionDirection = vl::norm( randomDirection( 
                    props.frictionDirection, random.uniform(), random.uniform(), 50.0 ) );

            simulation->setSamplingProperties( iDynOb, props );
        }
    }

    if( this->parentScene_ )
    	simulation->setSimulationState( this->startState_ );

	std::vector<ConstPhysicsObjectPtr> dynamicObjects = simulation->dynamicObjects();
	dynamicObjects.resize( this->dynamicObjects_.size() );
	for( size_t iDynOb = 0; iDynOb < dynamicObjects.size(); ++iDynOb )
		assert( dynamicObjects[iDynOb]->name() == this->dynamicObjects_[iDynOb]->name() );

	// for disabled objects, we don't want to store their paths.  So when an object
	//   becomes UNdisabled, we will store the starting time in startTimes and the 
	//   states in states.
	std::vector<ObjectStateList> objectStates( dynamicObjects.size() );
	std::vector<float> startTimes( dynamicObjects.size(), 
		this->startState_.time() );

	// we will use these to abstract away the concepts of
	//   'linearVelocity' and 'angularVelocity' so we don't
	//   have to dup all our code.
	boost::ptr_vector<ImpulseWrapper> impulseWrappers;
	impulseWrappers.push_back( new AngularImpulseWrapper(minAngularImpulse_) );
	impulseWrappers.push_back( new LinearImpulseWrapper(minLinearImpulse_) );

    PiecewisePath::FrameType lastFixedObjectsFrame = boost::numeric::bounds<PiecewisePath::FrameType>::lowest();
    PiecewisePath::FrameType firstFixedObjectsFrame = boost::numeric::bounds<PiecewisePath::FrameType>::highest();
	for( size_t iObject = 0; iObject < fixedObjectsHistory_.size(); ++iObject )
	{
		lastFixedObjectsFrame = std::max( lastFixedObjectsFrame, 
			fixedObjectsHistory_[iObject].endFrame() );
		firstFixedObjectsFrame = std::min( lastFixedObjectsFrame, 
			fixedObjectsHistory_[iObject].startFrame() );
	}
	float lastFixedObjectsTime = boost::numeric_cast<float>(lastFixedObjectsFrame) / 
		boost::numeric_cast<float>(fixedObjectsHistoryFramerate_);
	float firstFixedObjectsTime = boost::numeric_cast<float>(firstFixedObjectsFrame) / 
		boost::numeric_cast<float>(fixedObjectsHistoryFramerate_);

	// this is all the objects that have been activated in either the regular set
	//   or the "ghost" set
	Simulation::PhysicsObjectSet activeObjects;


	// We will be inserting ghost objects as we go, and we need to be able to 
	// map their dynamic ids in the new scene to the id of the stored history path
	// we have in fixedObjectsHistory_ 
	typedef std::pair<size_t, size_t> GhostObjectIdToHistoryPathId;
	std::deque<GhostObjectIdToHistoryPathId> ghostObjectIdToHistoryPathIdMap;
	for( std::vector<ConstPhysicsObjectPtr>::const_iterator iter = onlyInParentObjects_.begin();
		iter != onlyInParentObjects_.end(); ++iter )
	{
		activeObjects.insert( iter->get() );
		int iDynamicId = simulation->addObject( *iter, true );
		if( iDynamicId >= 0 )
		{
			std::vector<ConstPhysicsObjectPtr>::iterator dynItr = 
				std::find( parentDynamicObjects_.begin(), parentDynamicObjects_.end(), *iter );
			// since iDynamicId is positive, this must be the case:
			assert( dynItr != parentDynamicObjects_.end() );

			ghostObjectIdToHistoryPathIdMap.push_back( 
				GhostObjectIdToHistoryPathId( 
					boost::numeric_cast<size_t>( iDynamicId ),
					std::distance( parentDynamicObjects_.begin(), dynItr ) ) );
		}
	}

	// disable all the fixedObjects_ (that is, objects which the user has marked
	//   as not needing refinement)
	boost::dynamic_bitset<> disabledObjects( dynamicObjects.size() );

	// prep the list of active objects with the ones the user wanted simulated
	for( size_t iObject = 0; iObject < dynamicObjects.size(); ++iObject )
	{
		if( std::binary_search(this->fixedObjects_.begin(), this->fixedObjects_.end(), iObject) )
		{
			simulation->setDisabled( iObject, true );
			assert( simulation->getDisabled( iObject ) );
			disabledObjects.set( iObject );
		}
    }

    // if all the objects are active, there's no point doing any tracking
    if( disabledObjects.any() )
    {
	    for( size_t iObject = 0; iObject < dynamicObjects.size(); ++iObject )
        {
            if( !disabledObjects[iObject] )
                continue;

            activeObjects.insert( dynamicObjects.at(iObject).get() );

			// we also need to activate its corresponding "ghost" object:
			ObjectMap::const_iterator itr = childObjectToParentObjectLocal.find( 
				dynamicObjects.at(iObject).get() );

			// if it's a dynamic object, we need to create a ghost object for it:
			if( itr != childObjectToParentObjectLocal.end() )
			{
				// mark it as active
				activeObjects.insert( itr->second );

				// now add to the simulation
				// need to find the corresponding dynamic object id
				// normally I would use stl::find here, but shared_ptr doesn't have an
				//   operator== that takes a plain vanilla pointer
				std::vector<ConstPhysicsObjectPtr>::iterator dynItr = parentDynamicObjects_.begin();
				while( dynItr != parentDynamicObjects_.end() )
				{
					if( dynItr->get() == itr->second )
						break;

					++dynItr;
				}

				// If there is a corresponding object in the parent scene, insert a ghost object
				assert( dynItr != parentDynamicObjects_.end() );

				int iDynamicId = simulation->addObject( *dynItr, true );
				assert( iDynamicId >= 0 );

				ghostObjectIdToHistoryPathIdMap.push_back( 
					GhostObjectIdToHistoryPathId( 
						boost::numeric_cast<size_t>( iDynamicId ),
						std::distance( parentDynamicObjects_.begin(), dynItr ) ) );
			}
		}
	}

	// Since we have dup'd them into a new scene, we have to translate the
	//   pointers across
	for( std::vector<ConstPhysicsObjectPtr>::const_iterator iter = onlyInChildObjects_.begin();
		iter != onlyInChildObjects_.end(); ++iter )
	{
		// names are guaranteed to be unique, so this is okay
		ConstSceneGraphElementPtr elt = localScene->object( (*iter)->name() );
		const PhysicsObject* localPhysicsObject = dynamic_cast<const PhysicsObject*>( elt.get() );
		assert( localPhysicsObject != 0 );
		activeObjects.insert( localPhysicsObject );
	}

	Simulation::State prevState = simulation->getSimulationState();
	std::deque< Simulation::State > fullStates;

    // Keep track of the contact count for statistics
    std::vector<size_t> numContacts;
    std::vector<size_t> prevContactsPerObject;

	assert( minTimesteps <= maxTimesteps );
	for( size_t iStep = 0; (iStep < maxTimesteps) || (iStep < minTimesteps); ++iStep )
	{
		// store the full state twice per second
#ifndef _DEBUG
		if( (iStep % (this->samplingRate_/2)) == 0 )
#endif
        {
            if( this->timeReversed_ )
			    fullStates.push_front( prevState );
            else
			    fullStates.push_back( prevState );
        }

		for( size_t iObject = 0; iObject < dynamicObjects.size(); ++iObject )
		{
			if( !disabledObjects[iObject] )
			{
				if( objectStates[iObject].empty() )
					startTimes[iObject] = prevState.time();

				if( this->timeReversed_ )
					objectStates[iObject].push_front( prevState.state(iObject) );
				else
					objectStates[iObject].push_back( prevState.state(iObject) );
			}
		}

		if( simulation->asleep() && iStep >= minTimesteps && 
            ((!this->timeReversed_ && prevState.time() >= lastFixedObjectsTime) ||
            (this->timeReversed_ && prevState.time() <= firstFixedObjectsTime)) )
			break;

        if( ((this->timeReversed_ && prevState.time() < this->timeRange_.first) ||
            (!this->timeReversed_ && prevState.time() > this->timeRange_.second)) && iStep >= minTimesteps )
			break;

        /*
        if( this->timeReversed_ )
        {
            // We need to check whether or not we have gathered too much energy
            // if so, we terminate the sample
            BoundingBox3f bounds = this->scene_->bounds();
            double yMin = (bounds.minimum())[1];
            double yDiff = (bounds.maximum())[1] - yMin;
            assert( yDiff > 0.0f );

            // Check each object individually; would it be better to check them all together?

            double sumMaxPotentialEnergy = 0.0;
            double sumCurrentEnergy = 0.0;
            bool objectIsBad = false;

            for( size_t iDynObject = 0; iDynObject < dynamicObjects.size(); ++iDynObject )
            {
                const RigidDynamicState& state = prevState.state( iDynObject );

                // todo this could be maybe be cached
                dMass mass = dynamicObjects.at(iDynObject)->getOdeMassProps();
                double maxPotentialEnergy = 9.81*mass.mass*yDiff;
                sumMaxPotentialEnergy += maxPotentialEnergy;

                // linear part
                vl::Vec3f linearVel = state.linearVelocity();
                double linearKineticEnergy = 
                    0.5 * mass.mass * vl::dot( linearVel, linearVel );

                vl::Vec3f angularVel = state.angularVelocity();

                double angularKineticEnergy = 0.0;
                // angular part
                for( int i = 0; i < 3; ++i )
                    for( int j = 0; j < 3; ++j )
                        angularKineticEnergy += 0.5 * mass.I[4*i + j] * angularVel[i]*angularVel[j];

                // note that this doesn't include possible additional sources of energy, like
                // motors
                double currentPotentialEnergy = 0.5 * mass.mass * ((state.position())[1] - yMin);

                double currentTotalEnergy = angularKineticEnergy + linearKineticEnergy + currentPotentialEnergy;

                sumCurrentEnergy += currentTotalEnergy;

                // TODO fudge factor here?
                if( currentTotalEnergy > maxPotentialEnergy )
                    objectIsBad = true;
            }

            if( objectIsBad || (sumCurrentEnergy > sumMaxPotentialEnergy) )
                break;
        }
        */

#ifdef _DEBUG
		// check that the backup state is okay before we do anything
		for( size_t iObject = 0; iObject < dynamicObjects.size(); ++iObject )
		{
			const RigidDynamicState& state = prevState.state( iObject );
			assert( !isBad(state.linearVelocity()) && !isBad(state.angularVelocity()) );
			assert( !isBad(state.position()) && !isBad(state.orientation()) );
		}
#endif // _DEBUG

		// disabled objects need to have their state updated appropriately
		for( size_t iObject = 0; iObject < disabledObjects.size(); ++iObject )
		{
			if( !disabledObjects[iObject] )
				continue;

			// regardless of whether we activate it or not, we should set its state 
			//   correctly for future simulation
			float scaledTime = simulation->time() * fixedObjectsHistoryFramerate_;

			DynamicObjectIdMap::const_iterator dynamicIdItr = 
				dynamicObjectToParentDynamicObject_.find(iObject);
			assert( dynamicIdItr != dynamicObjectToParentDynamicObject_.end() );
			const PiecewisePath& path = fixedObjectsHistory_.at( dynamicIdItr->second );

            if( scaledTime < path.startFrame() || scaledTime > path.endFrame() )
                simulation->setDisabled(iObject, false);
            else
            {
			    RigidDynamicState state(
				    path.position( scaledTime ),
				    path.rotation( scaledTime ),
				    path.linearVelocity( scaledTime ) * static_cast<float>(fixedObjectsHistoryFramerate_),
				    path.angularVelocity( scaledTime ) * static_cast<float>(fixedObjectsHistoryFramerate_) );
			    simulation->setSimulationState( iObject, state );
            }

            if( !simulation->getDisabled(iObject) )
			{
				// in this case, the object has been newly enabled during simulation
				disabledObjects.set( iObject, false );

				// we set the object as active
				activeObjects.insert( dynamicObjects.at(iObject).get() );

				// we also need to activate its corresponding "ghost" object:
				ObjectMap::const_iterator itr = childObjectToParentObjectLocal.find( 
					dynamicObjects.at(iObject).get() );
				// if it's a dynamic object, we need to create a ghost object for it:
				if( itr != childObjectToParentObjectLocal.end() )
				{
					// mark it as active
					activeObjects.insert( itr->second );

					// now add to the simulation
					// need to find the corresponding dynamic object id
					// normally I would use stl::find here, but shared_ptr doesn't have an
					//   operator== that takes a plain vanilla pointer
					std::vector<ConstPhysicsObjectPtr>::iterator dynItr = parentDynamicObjects_.begin();
					while( dynItr != parentDynamicObjects_.end() )
					{
						if( dynItr->get() == itr->second )
							break;

						++dynItr;
					}

					// here this should always be the case because if the object doesn't
					//   have a corresponding version in the parent it should be in the
					//   onlyInChildObjects_ set.
					assert( dynItr != parentDynamicObjects_.end() );
					int iDynamicId = simulation->addObject( *dynItr, true );
					assert( iDynamicId >= 0 );

					ghostObjectIdToHistoryPathIdMap.push_back( 
						GhostObjectIdToHistoryPathId( 
							boost::numeric_cast<size_t>( iDynamicId ),
							std::distance( parentDynamicObjects_.begin(), dynItr ) ) );
				}
			}
		}

        // If every object is active, there's no point tracking the ghost objects any more
        if( disabledObjects.none() )
        {
            for( std::deque<GhostObjectIdToHistoryPathId>::const_iterator ghostItr = ghostObjectIdToHistoryPathIdMap.begin();
                ghostItr != ghostObjectIdToHistoryPathIdMap.end(); ++ghostItr )
            {
                simulation->removeObject( ghostItr->first );
            }

            ghostObjectIdToHistoryPathIdMap.clear();
        }

		// update state of all the "ghost" objects
		for( std::deque<GhostObjectIdToHistoryPathId>::iterator ghostItr = 
				ghostObjectIdToHistoryPathIdMap.begin();
			ghostItr != ghostObjectIdToHistoryPathIdMap.end(); )
		{
			float scaledTime = simulation->time() * fixedObjectsHistoryFramerate_;
			const PiecewisePath& path = fixedObjectsHistory_.at( ghostItr->second );

            // Ghost objects can't exist outside of the time that we have motion stored for them:
            if( (scaledTime < path.startFrame() && this->timeReversed_) || 
                (scaledTime > path.endFrame() && !this->timeReversed_) )
            {
                simulation->removeObject( ghostItr->first );
                ghostItr = ghostObjectIdToHistoryPathIdMap.erase( ghostItr );

                continue;
            }

			// velocity doesn't matter for ghost objects since they'll never be
			//   simulated:
			RigidDynamicState state(
				path.position( scaledTime ),
				path.rotation( scaledTime ),
				vl::vl_0,
				vl::vl_0 );

			// ghost objects should always be disabled
			simulation->setDisabled(ghostItr->first, true);
			simulation->setSimulationState( ghostItr->first, state );

            ++ghostItr;
		}

		simulation->setActiveObjects( activeObjects );

		for( size_t i = 0; i < this->multiplier_; ++i )
			simulation->step( (this->timeReversed_ ? -1.0 : 1.0) / boost::numeric_cast<double>(simulationRate_) );


        /*
        std::vector<Simulation::ContactPoint> contactPoints = 
            simulation->contactPoints();
        if( vl::len(simulation->getLinearVelocity(0)) > 1e-3 )
            numContacts.push_back( contactPoints.size() );

        std::vector<size_t> nextContactsPerObject( dynamicObjects.size(), 0 );
        for( std::vector<Simulation::ContactPoint>::const_iterator itr = contactPoints.begin();
            itr != contactPoints.end(); ++itr )
        {
            if( itr->objects.first >= 0 )
                nextContactsPerObject[ itr->objects.first ]++;
            if( itr->objects.second >= 0 )
                nextContactsPerObject[ itr->objects.second ]++;
        }
        */

        /*
        if( this->timeReversed_ && !prevContactsPerObject.empty() )
        {
            std::vector<size_t> contactsChanged;

            for( size_t iObject = 0; iObject < dynamicObjects_.size(); ++iObject )
            {
                if( nextContactsPerObject[iObject] != prevContactsPerObject[iObject] )
                    contactsChanged.push_back(iObject);
            }

            if( !contactsChanged.empty() )
            {
                std::vector<Simulation::State> states =
                    this->sampleForwardToState( simulation, 240, random );

                // Now, need to insert it into the stream:
                for( size_t iFrame = 0; iFrame < states.size(); ++iFrame )
                {
	                for( size_t iObject = 0; iObject < dynamicObjects.size(); ++iObject )
	                {
		                if( disabledObjects[iObject] )
                            continue;

                        objectStates[iObject].push_front( states[states.size() - iFrame - 1].state(iObject) );
                    }
                }

                prevState = states.front();
                simulation->setSimulationState( prevState );

                iStep += states.size() - 1;
                prevContactsPerObject.clear();

                //if( iStep > 350 )
                //    break;
                //break;
                continue;
            }
        }
        */

        //nextContactsPerObject.swap( prevContactsPerObject );

        /*
        // For some debugging purposes we want to dump out lots of statistics about where the simulation
        // ends up at when the vertical position is 0.  
        {
            Simulation::State nextState = simulation->getSimulationState();
            for( size_t iObject = 0; iObject < dynamicObjects.size(); ++iObject )
            {
                // Not threadsafe to open this file in multiple threads
                static boost::mutex mutexLock;
                boost::mutex::scoped_lock lock( mutexLock );

                // write start and end position, orientation, linear, and angular velocity
                // for object
                vl::Vec3f prevPos = prevState.state(iObject).position();
                vl::Vec3f nextPos = nextState.state(iObject).position();
                // Only dump if the object passed through the plane y=0 at the given time
                if( prevPos[1] < 0.0 || nextPos[1] >= 0.0 )
                    continue;

                RigidDynamicState startState = this->startState_.state(iObject);
                RigidDynamicState endState = nextState.state(iObject);

                ConstPhysicsObjectPtr object = dynamicObjects.at(iObject);
                std::string outputFile( object->name() + ".bin" );

                // Need to open the file
                int fd = open( outputFile.c_str(), O_WRONLY | O_CREAT );
                int err = errno;

                assert( fd != -1 );

                int res;
#ifndef _WINDOWS
                struct flock fl = { F_WRLCK, SEEK_END, 0,       2*sizeof(RigidDynamicState),     0 };
                struct timespec ts = { 0, 1000 };

                // Try to acquire a lock
                res = fcntl(fd, F_SETLKW, &fl);
                while( fcntl(fd, F_SETLKW, &fl) != -1 )
                    assert( errno == EINTR ); // Only signal interruptions are ok
#endif

                lseek( fd, 0, SEEK_END );

                res = write( fd, &startState, sizeof( RigidDynamicState ) );
                assert( res != -1 );
                res = write( fd, &endState, sizeof( RigidDynamicState ) );
                assert( res != -1 );

#ifndef _WINDOWS
                // release lock
                fl.l_type = F_UNLCK;
                res = fcntl(fd, F_SETLK, &fl);
                assert( res != -1 );
#endif

                res = close(fd);
                assert( res != -1 );
            }
        }
        */

        for( size_t iObject = 0; iObject < dynamicObjects.size(); ++iObject )
		{
			// perturbing disabled objects could activate them in some simulators
			if( disabledObjects[iObject] )
				continue;

			const RigidDynamicState& prevObjectState = prevState.state( iObject );

			for( boost::ptr_vector<ImpulseWrapper>::iterator wrapperItr = impulseWrappers.begin();
				wrapperItr != impulseWrappers.end(); ++wrapperItr )
			{
				vl::Vec3f prevValue = wrapperItr->get( prevObjectState );
				vl::Vec3f nextValue = wrapperItr->get( *simulation, iObject );
				assert( !isBad( prevValue ) );
				assert( !isBad( nextValue ) );

				vl::Vec3f origImpulse = nextValue - prevValue;
				assert( !isBad(origImpulse) );
				if( vl::len(origImpulse) < wrapperItr->minImpulse() )
					continue;

				vl::Vec3f perturbation = randomImpulse( origImpulse, random );

				wrapperItr->perturb( *simulation, iObject, perturbation );
			}
		}

		// after applying impulses inside of simulation, need to re-get the state:
		prevState = simulation->getSimulationState();
	}

	// always save the end state
	fullStates.push_back( simulation->getSimulationState() );


	std::vector<CompressedPath> compressedPaths( dynamicObjects.size() );
	std::vector< std::vector<vl::Vec3f> > cmPaths( dynamicObjects.size() );

	// don't know how well this is going to work
	// @todo check this
	float quadraticCoeff = -0.5 * 9.81 * boost::numeric_cast<double>(this->dt_) * boost::numeric_cast<double>(this->dt_);
	std::vector<float> objectives;
	objectives.reserve( dynamicObjects.size() * objectives_.size() );

	size_t originalBytes = 0;
	size_t compressedBytes = 0;

	for( size_t iObject = 0; iObject < dynamicObjects.size(); ++iObject )
	{
		DynamicObjectIdMap::const_iterator dynamicIdItr = 
			dynamicObjectToParentDynamicObject_.find(iObject);

		if( !objectStates.at(iObject).empty() )
		{
			long startTime = boost::numeric_cast<long>(startTimes[iObject] * this->samplingRate_);

			// If we simulated backwards, startTime actually equals the time of the _last_ frame,
			// so we need to subtract off the number of frames
			if( this->timeReversed_ )
				startTime -= objectStates[iObject].size();

			compressedPaths.at(iObject) = 
				compress( objectStates[iObject], quadraticCoeff, startTime, this->samplingRate_, 1e-3, 1e-3 );

			size_t originalSize = objectStates.at(iObject).size() * sizeof(RigidStaticState);
			size_t newSize = compressedPaths.at(iObject).size();

			originalBytes += originalSize;
			compressedBytes += newSize;
		}

		boost::scoped_ptr<PiecewisePath> path;
        if( (dynamicIdItr != dynamicObjectToParentDynamicObject_.end()) )
		{
			const PiecewisePath& fixedObjectPath = 
				fixedObjectsHistory_.at( dynamicIdItr->second );
			path.reset( new PiecewisePath(compressedPaths.at(iObject), fixedObjectPath, timeReversed_) );
		}
		else
		{
			path.reset( new PiecewisePath(compressedPaths.at(iObject)) );
		}

		for( std::vector<ObjectivePtr>::const_iterator itr = objectives_.begin();
			itr != objectives_.end(); ++itr )
		{
			objectives.push_back( 
				(*itr)->evaluate( 
					*path, 
					this->samplingRate_,
					dynamicObjects.at(iObject),
					*cameraSource_ ) );
		}

		cmPaths[iObject].reserve( objectStates[iObject].size() );
		for( size_t i = 0; i < objectStates[iObject].size(); ++i )
			cmPaths[iObject].push_back( toVec3f(objectStates[iObject].at(i).position()) );
	}
    /*
	{
		std::cout << "Compressed: original size: " << originalBytes << " bytes; new size " 
			<< compressedBytes << "; ratio: " << 
				(boost::numeric_cast<double>(originalBytes) / boost::numeric_cast<double>(compressedBytes) ) << std::endl;
	}
    */

    /*
    {
        static boost::mutex numContactsMutex;
        boost::mutex::scoped_lock lock( numContactsMutex );
        std::ofstream ofs( "countContacts.bin", std::ios::binary | std::ios::app );
        assert( ofs );
        dumpArray( numContacts, ofs );
    }
    */

    /*
    // dump end position
    {
        static boost::mutex mutexLock;
        boost::mutex::scoped_lock lock( mutexLock );

        static std::ofstream ofs( "endPositions.bin", std::ios::binary | std::ios::app );
        assert( ofs );
        assert( dynamicObjects.size() == 1 );
        Simulation::State lastState = simulation->getSimulationState();
        RigidDynamicState endState = lastState.state(0);
        vl::Vec3f pos = endState.position();
        ofs.write( reinterpret_cast<const char*>( &pos ), sizeof( vl::Vec3f ) );
    }
    */

    {
        ++numSamples_;

        BoundingBox3f box = localScene->bounds();
		float epsilon = 5e-4 * vl::len(box.maximum() - box.minimum());

		boost::mutex::scoped_lock lock( this->pathPoolMutex_ );
		Path* path = new (pathPool_.malloc())
			Path( compressedPaths, fullStates, cmPaths, this->samplingRate_, objectives, epsilon );
		return path;
	}
}

std::vector<PiecewisePath> SimulationTree::decode( const Path* path )
{
	std::vector<PiecewisePath> result;
	result.reserve( this->dynamicObjects_.size() );
	for( size_t iObject = 0; iObject < dynamicObjects_.size(); ++iObject )
	{
		const CompressedPath& compressedPath = path->path(iObject);

		DynamicObjectIdMap::const_iterator dynIdItr = 
			dynamicObjectToParentDynamicObject_.find( iObject );
		if( dynIdItr != dynamicObjectToParentDynamicObject_.end() )
		{
			result.push_back( PiecewisePath(compressedPath, 
				fixedObjectsHistory_.at(dynIdItr->second),
                this->timeReversed_) );
		}
		else
		{
			result.push_back( PiecewisePath(compressedPath) );
		}
	}

	return result;
}

boost::dynamic_bitset<> SimulationTree::activePaths() const
{
	ReaderWriterLock::ScopedReaderLock lock( this->activePathsLock_ );
	boost::dynamic_bitset<> result( this->activePaths_ );
	return result;
}

std::deque<ConstraintPtr> SimulationTree::constraints() const
{
	ReaderWriterLock::ScopedReaderLock lock( this->constraintsLock_ );
	std::deque<ConstraintPtr> result( this->constraints_ );
	return result;
}

const Path* SimulationTree::path( size_t pathId ) const
{
	ReaderWriterLock::ScopedReaderLock lock( this->allPathsLock_ );
	return this->allPaths_.at( pathId );
}

size_t SimulationTree::numPaths() const
{
	ReaderWriterLock::ScopedReaderLock lock( this->allPathsLock_ );
	return this->allPaths_.size();
}

const vl::Vec3f toPosition( const char* state )
{
	const RigidStaticState* rigidState = 
		reinterpret_cast<const RigidStaticState*>( state );
	return toVec3f( rigidState->position() );
}

void SimulationTree::addPath( const Path& path )
{
	Path* p;
	{
		boost::mutex::scoped_lock lock( this->pathPoolMutex_ );
		p = new (pathPool_.malloc())	Path( path );
	}
	{
		boost::mutex::scoped_lock lock( this->nodesToAddMutex_ );
		pathsToAdd_.push_back( p );
	}
}

void SimulationTree::addPath( 
	const std::vector<CompressedPath>& compressedPaths, 
	const std::vector<Path::TimedSnapshot>& states, 
	const std::vector< std::vector<vl::Vec3f> >& cmPaths,
	short frameRate,
	const std::vector<float>& objectives )
{
	Path* path;
	{
		boost::mutex::scoped_lock lock( this->pathPoolMutex_ );
		path = new (pathPool_.malloc())
			Path( compressedPaths, states, cmPaths, frameRate, objectives );
	}
	{
		boost::mutex::scoped_lock lock( this->nodesToAddMutex_ );
		pathsToAdd_.push_back( path );
	}
}

std::vector<size_t> SimulationTree::fixedObjects() const
{
	return this->fixedObjects_;
}

std::vector<vl::Vec3f> SimulationTree::pathRoot(size_t iObject) const
{
	if( iObject < pathRoots_.size() )
		return pathRoots_.at( iObject );
	else
		return std::vector<vl::Vec3f>();
}

void SimulationTree::sample( RandomStream& random, bool backwards )
{
#ifndef _WIN32
	std::cout << "Computing sample... ";
#endif

	Path* path = this->sample(
		boost::numeric::bounds<size_t>::highest(),
		14,
		random );

    if( path == 0 )
        return;

#ifdef STREAMING
	/*
	size_t compressedSize = 0;
	for( std::vector<CompressedPath>::const_iterator itr = compressed.begin();
		itr != compressed.end(); ++itr )
	{
		compressedSize += itr->positionTimes.size();

		for( size_t i = 0; i < 3; ++i )
			compressedSize += itr->positionLinearCoeffs[i].size();

		compressedSize += itr->positionQuadraticCoeffs[1].size();

		compressedSize += itr->rotationTimes.size();
		for( size_t i = 0; i < 3; ++i )
			compressedSize += itr->rotationCoeffs[i].size();
	}

	// 4 float * (3 position components + 3 orientation components) *
	//   num. of objects * num. of states
	size_t uncompressedSize = 4*6*compressed.size()*states.size();
	double ratio = boost::numeric_cast<double>( uncompressedSize ) / 
		boost::numeric_cast<double>( compressedSize );
	std::cout << "Compressed " << states.size() << " frame animation.\n"
		<< "  original size:   " << uncompressedSize << " bytes; \n"
		<< "  compressed size: " << compressedSize << " bytes; \n"
		<< "  ratio:           " << ratio << std::endl;
		*/

	path->write( ofs_ );
	ofs_.flush();

	{
		boost::mutex::scoped_lock lock( this->pathPoolMutex_ );
		pathPool_.free( path );
	}
#else
	{
		boost::mutex::scoped_lock lock( this->nodesToAddMutex_ );
		pathsToAdd_.push_back( path );
	}
#endif
}

void writeState( const Simulation::State& state, std::ostream& os )
{
	float time = state.time();
	os.write( reinterpret_cast<const char*>(&time), sizeof( float ) );

	for( size_t iState = 0; iState < state.stateCount(); ++iState )
	{
		const RigidDynamicState& s = state.state( iState );
		vl::Vec3f position = s.position();
		Quaternion<float> orientation = s.orientation();
		vl::Vec3f linearVel = s.linearVelocity();
		vl::Vec3f angularVel = s.angularVelocity();

		os.write( reinterpret_cast<const char*>( position.Ref() ), sizeof( vl::Vec3f ) );
		os.write( reinterpret_cast<const char*>( &orientation ), sizeof( Quaternion<float> ) );
		os.write( reinterpret_cast<const char*>( linearVel.Ref() ), sizeof( vl::Vec3f ) );
		os.write( reinterpret_cast<const char*>( angularVel.Ref() ), sizeof( vl::Vec3f ) );
	}
}

Simulation::State readState( std::istream& is, size_t stateCount )
{
	float time;
	is.read( reinterpret_cast<char*>(&time), sizeof( float ) );

	std::vector<RigidDynamicState> states;
	states.reserve( stateCount );
	for( size_t iState = 0; iState < stateCount; ++iState )
	{
		vl::Vec3f position;
		Quaternion<float> orientation;
		vl::Vec3f linearVel;
		vl::Vec3f angularVel;

		is.read( reinterpret_cast<char*>( position.Ref() ), sizeof( vl::Vec3f ) );
		is.read( reinterpret_cast<char*>( &orientation ), sizeof( Quaternion<float> ) );
		is.read( reinterpret_cast<char*>( linearVel.Ref() ), sizeof( vl::Vec3f ) );
		is.read( reinterpret_cast<char*>( angularVel.Ref() ), sizeof( vl::Vec3f ) );
		
		states.push_back( 
			RigidDynamicState(position, orientation, linearVel, angularVel) );
	}

	return Simulation::State( states, time );
}

/*
class MemoryStateComparator
{
public:
	MemoryStateComparator( std::istream& i )
		: ifs(i)
	{
		_CrtMemCheckpoint( &s1 );

		this->startOffset = ifs.tellg();
	}

	void difference()
	{
		_CrtMemState s2, s3;
		_CrtMemCheckpoint( &s2 );
		_CrtMemDifference( &s3, &s1, &s2);
        _CrtMemDumpStatistics( &s3 );

		size_t endOffset = ifs.tellg();
		size_t diff = endOffset - startOffset;
		std::ostringstream oss;
		oss << "Bytes in file: " << diff << "\n";
		_CrtDbgReport( _CRT_WARN, NULL, NULL, "tree", oss.str().c_str() );
		std::cout << diff;
	}

private:
	_CrtMemState s1;
	std::istream& ifs;
	size_t startOffset;
};
*/

void SimulationTree::read( const std::string& filename )
{
	std::ifstream ifs( filename.c_str(), std::ios::binary );
	if( !ifs )
		throw IOException( "Unable to open file '" + filename + "' for reading." );

	ifs.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

	try
	{
		std::ifstream::pos_type m = filelen( ifs );

		// opens with a "type" field
		size_t type;
		ifs.read( reinterpret_cast<char*>(&type), sizeof( size_t ) );

		// first word tells how many objects, this will be used to determine
		//   the state size throughout the tree:
		size_t numDynamicObjects;
		ifs.read( reinterpret_cast<char*>(&numDynamicObjects), sizeof( size_t ) );

		if( numDynamicObjects != this->dynamicObjects_.size() )
		{
			size_t myDynOb = this->dynamicObjects().size();
			std::ostringstream oss;
			oss << "Mismatch in object counts; file has " << numDynamicObjects 
				<< " but actual number is " << myDynOb;
			throw IOException( oss.str() );
		}

		if( type == 101 )
		{
			size_t frameRate;
			ifs.read( reinterpret_cast<char*>( &frameRate ), sizeof( size_t ) );

			std::deque<Path*> paths;
			while( ifs.tellg() < m )
			{
				boost::mutex::scoped_lock lock( this->pathPoolMutex_ );
				Path* path = new (pathPool_.malloc())
					Path( ifs, numDynamicObjects, frameRate, filename );
				paths.push_back( path );

				/*
				dumping paths for our constraint testing code

				static std::ofstream ofs( "testPaths.bin", std::ios::binary );
				for( size_t i = 0; i < path->treeCount(); ++i )
				{
					std::vector<vl::Vec3f> points = path->tree(i).points();
					dumpArray( points, ofs );
				}
				*/
			}

			{
				boost::mutex::scoped_lock lock( this->nodesToAddMutex_ );
				std::copy( paths.begin(), paths.end(), 
					std::back_inserter(pathsToAdd_) );
			}
		}
		else
		{
			throw IOException( "Unknown type" );
		}
	}
	catch( std::ifstream::failure& )
	{
		throw IOException( "Caught exception reading from file '" + filename + "'" );
	}
}

const CameraSource& SimulationTree::cameraSource() const
{
	return *this->cameraSource_;
}

void SimulationTree::addCallback( SimulationNodeCallback* callback )
{
	{
		boost::mutex::scoped_lock lock( this->callbackMutex_ );
		this->callbacks_.push_back( callback );
	}
}

void SimulationTree::removeCallback( SimulationNodeCallback* callback )
{
	boost::mutex::scoped_lock lock( this->callbackMutex_ );
	this->callbacks_.erase( std::remove( callbacks_.begin(), callbacks_.end(), callback ) );
}

std::vector<vl::Vec3f> approximatePath( const std::vector<vl::Vec3f>& points, float epsilon )
{
	if( points.empty() )
		return points;

	// need to trim points that are too close together:
	std::vector<vl::Vec3f> result;
	result.reserve( points.size()/10 );

	{
		result.push_back( points.front() );
		size_t start = 0;

		while( start < points.size() - 1 )
		{
			size_t low = start;
			size_t high = start + 1;

			while( high < points.size() )
			{
				if( error( start, high, points ) > epsilon )
					break;

				low = high;
				high += (high - start);
			}

			high = std::min( high, points.size() );
			while( high - low > 1 )
			{
				size_t probe = (high + low) / 2;
				if( error( start, probe, points ) > epsilon )
					high = probe;
				else
					low = probe;
			}

			size_t end = std::max(low, start+1);
            
			result.push_back( points.at(end) );
			start = end;
		}
	}

	return result;
}

Path::Path( 
	const std::vector<CompressedPath>& paths, 
	const std::deque<Simulation::State>& states, 
	const std::vector< std::vector<vl::Vec3f> >& cmPaths,
	short frameRate, 
	const std::vector<float>& scores,
	float epsilon )
	: scores_( scores ), frameRate_(frameRate), paths_(paths)
{
	std::cout << "Got epsilon: " << epsilon << std::endl;
	const size_t nSamples = 10;

	assert( paths.size() <= states.front().stateCount() );
	for( size_t iObject = 0; iObject < paths.size(); ++iObject )
	{
		obbTrees_.push_back(
			new OBBTree( approximatePath(cmPaths[iObject], epsilon) ) );
	}

	// need to strip the states down to just the static portion for
	//   space reasons
	this->snapshots_.reserve( nSamples );
	size_t skip = states.size() / nSamples;
	for( size_t i = 0; i < nSamples; ++i )
	{
		Simulation::State state = states.at(skip*i);
		std::vector<RigidStaticState> states;
		assert( state.stateCount() >= obbTrees_.size() );
		states.reserve( obbTrees_.size() );
		for( size_t j = 0; j < obbTrees_.size(); ++j )
			states.push_back( state.state(j) );

		this->snapshots_.push_back( 
			std::make_pair( state.time(), states ) );
	}
}

Path::Path( const std::vector<CompressedPath>& paths, 
	const std::vector<TimedSnapshot>& snapshots, 
	const std::vector< std::vector<vl::Vec3f> >& cmPaths,
	short frameRate, 
	const std::vector<float>& scores )
	: scores_( scores ), frameRate_(frameRate), paths_(paths), snapshots_( snapshots )
{
	for( size_t iObject = 0; iObject < paths.size(); ++iObject )
	{
		obbTrees_.push_back(
			new OBBTree( cmPaths[iObject] ) );
	}
}

Path::Path( const Path& other )
	:	paths_(other.paths_),
		snapshots_(other.snapshots_),
		frameRate_(other.frameRate_),
		scores_(other.scores_),
		filename_(other.filename_),
		pathPositions_(other.pathPositions_)
{
	for( size_t i = 0; i < other.obbTrees_.size(); ++i )
		obbTrees_.push_back( new OBBTree(other.obbTrees_[i].points()) );
}

CompressedPath Path::readPath( std::istream& ifs ) const
{
	CompressedPath compressedPath;

    ifs.read( reinterpret_cast<char*>( &compressedPath.startFrame ), sizeof(int) );

	readArray( compressedPath.positionTimes, ifs );
	assert( compressedPath.positionTimes.capacity() == compressedPath.positionTimes.size() );
	for( size_t i = 0; i < 3; ++i )
	{
		readArray( compressedPath.positionLinearCoeffs[i], ifs );
		assert( compressedPath.positionLinearCoeffs[i].capacity() == 
			compressedPath.positionLinearCoeffs[i].size() );
	}

	readArray( compressedPath.positionQuadraticCoeffs[1], ifs );
	assert( compressedPath.positionQuadraticCoeffs[1].capacity() == 
		compressedPath.positionQuadraticCoeffs[1].size() );

	readArray( compressedPath.rotationTimes, ifs );
	assert( compressedPath.rotationTimes.capacity() == 
		compressedPath.rotationTimes.size() );
	for( size_t i = 0; i < 3; ++i )
	{
		readArray( compressedPath.rotationCoeffs[i], ifs );
		assert( compressedPath.rotationCoeffs[i].capacity() == 
			compressedPath.rotationCoeffs[i].size() );
	}

	return compressedPath;
}

Path::Path( std::istream& ifs, size_t numPaths, short frameRate, const std::string& filename )
	: frameRate_( frameRate ), filename_( filename )
{
	this->pathPositions_.reserve( numPaths );
	for( size_t i = 0; i < numPaths; ++i )
	{
		pathPositions_.push_back( ifs.tellg() );
		readPath(ifs);
		//paths_.push_back( readPath(ifs) );

		{
			// read tree from input stream
			std::vector<vl::Vec3f> points;
			readArray( points, ifs );
			size_t beforeSize = points.size();
			points.swap( approximatePath(points, 1e-3) );
			size_t afterSize = points.size();
			this->obbTrees_.push_back( new OBBTree(points) );
		}
	}

	{
		size_t nSamples;
		ifs.read( reinterpret_cast<char*>(&nSamples), sizeof(size_t) );
		this->snapshots_.reserve( nSamples );
		for( size_t i = 0; i < nSamples; ++i )
		{
			float time;
			ifs.read( reinterpret_cast<char*>(&time), sizeof(float) );

			std::vector<RigidStaticState> state;
			readArray( state, ifs, numPaths );

			this->snapshots_.push_back( std::make_pair(time, state) );
		}
	}

	{
		readArray( scores_, ifs );
	}
}

Path::~Path()
{
}

size_t Path::treeCount() const
{
	return this->obbTrees_.size();
}

const OBBTree& Path::obbTree( size_t iObject ) const
{
	return this->obbTrees_.at(iObject);
}

float Path::score( size_t iObject, size_t iMetric ) const
{
	assert( iObject < this->obbTrees_.size() );
	assert( iMetric < this->metricCount() );
	return this->scores_.at( iObject * this->metricCount() + iMetric );
}

size_t Path::metricCount() const
{
	assert( (scores_.size() % obbTrees_.size()) == 0 );
	return (this->scores_.size() / this->obbTrees_.size());
}

void Path::write( std::ostream& os ) const
{
	for( size_t iObject = 0; iObject < this->paths_.size(); ++iObject )
	{
		const CompressedPath& path = this->paths_.at( iObject );

        os.write( reinterpret_cast<const char*>( &path.startFrame ), sizeof(int) );

        dumpArray( path.positionTimes, os );

		for( size_t i = 0; i < 3; ++i )
			dumpArray( path.positionLinearCoeffs[i], os );

		dumpArray( path.positionQuadraticCoeffs[1], os );

		dumpArray( path.rotationTimes, os );
		for( size_t i = 0; i < 3; ++i )
			dumpArray( path.rotationCoeffs[i], os );

		std::vector<vl::Vec3f> points = this->obbTrees_.at(iObject).points();
		dumpArray( points, os );
	}

	size_t nSamples = this->snapshots_.size();
	os.write( reinterpret_cast<const char*>(&nSamples), sizeof(size_t) );

	for( size_t i = 0; i < nSamples; ++i )
	{
		float time = snapshots_.at(i).first;
		os.write( reinterpret_cast<const char*>(&time), sizeof(float) );
		dumpArray( snapshots_.at(i).second, os, snapshots_.at(i).second.size() );
	}

	dumpArray( scores_, os );
}

size_t Path::frameRate() const
{
	return this->frameRate_;
}

std::vector<Simulation::State> Path::samples(size_t numSamples) const
{
	std::vector<Simulation::State> result;
	result.reserve( numSamples );

	assert( numSamples <= this->snapshots_.size() );
	size_t skip = this->snapshots_.size() / numSamples;
	for( size_t i = 0; i < this->snapshots_.size(); i += skip )
	{
		std::vector<RigidDynamicState> states( this->snapshots_[i].second.begin(),
			this->snapshots_[i].second.end() );
		result.push_back( Simulation::State(states, this->snapshots_[i].first) );
	}

	return result;
}

std::vector<Simulation::State> Path::samples() const
{
	return this->samples( this->snapshots_.size() );
}

const CompressedPath& Path::path(size_t iObject) const
{
	if( this->paths_.empty() )
	{
		// load it off disk
		assert( !this->filename_.empty() );
		assert( pathPositions_.size() == this->obbTrees_.size() );
		std::ifstream ifs( this->filename_.c_str(), std::ios::binary );
		this->paths_.reserve( pathPositions_.size() );
		for( size_t i = 0; i < pathPositions_.size(); ++i )
		{
			ifs.seekg( pathPositions_[i] );
			this->paths_.push_back( readPath(ifs) );
		}
	}

	return this->paths_.at( iObject );
}

} // namespace planning
