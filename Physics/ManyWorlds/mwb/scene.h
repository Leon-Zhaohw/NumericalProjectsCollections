#ifndef __SCENE_H__

#define __SCENE_H__

#include "physicsFwd.h"
#include "physicsObject.h"

#include "twigg/stlext.h"
#include "twigg/keyable.h"
#include "twigg/camera.h"
#include "twigg/random.h"

#include <ode/ode.h>

#include <boost/type_traits.hpp> 
#include <boost/graph/adjacency_list.hpp>
#include <boost/logic/tribool.hpp>

#include <set>
#include <map>

namespace planning {

class PiecewisePath;

// These will be used to bias the sampling
struct SamplingProperties
{
    SamplingProperties()
        : bounce( boost::indeterminate ), slide( boost::indeterminate ), frictionDirection( vl::vl_0 ) {}

    boost::tribool bounce;
    boost::tribool slide;
    vl::Vec3f frictionDirection;
};


// this will wrap up a camera together with its keyframes and stick it in 
//   one big Named package
class CameraWrapper
	: public Named
{
public:
	typedef boost::shared_ptr<AttributeWithKeys> AttributeWithKeysPtr;

	CameraWrapper( const std::string& name, boost::shared_ptr<Camera> camera );
	~CameraWrapper();

	boost::shared_ptr<Camera> camera();
	boost::shared_ptr<const Camera> camera() const;

	MechModelObjectPtr toMechModel();

	float nextKey(float time) const;
	float prevKey(float time) const;
	void deleteKey(float time);
	void setKey(float time);
	void setTime(float time);
	float firstKey() const;
	float lastKey() const;

	typedef std::deque<AttributeWithKeysPtr>::iterator attribute_iterator;
	attribute_iterator keyable_begin();
	attribute_iterator keyable_end();
	CameraWrapper* clone() const;

private:
	boost::shared_ptr<Camera> camera_;

	std::deque<AttributeWithKeysPtr> keyable_;
};

class Simulation
	: boost::noncopyable
{
public:
	// okay, so different types of simulators have different
	//   types of stored state.  Here we will treat state as
	//   just a collection of bytes and return chunks of 
	//   it according to the stored stride
	class State
	{
	public:
		State()
			:	time_(0.0), stride_( sizeof(RigidStaticState) ) {}

		// note that this requires that a byte copy is sufficient for
		//   StateType.
		template <typename StateType>
		State( const std::vector<StateType>& states, float time )
			:	time_(time)
		{
			BOOST_STATIC_ASSERT( sizeof(StateType) >= sizeof(RigidDynamicState) );

			this->stride_ = sizeof( StateType );
			if( !states.empty() )
			{
				this->states_.resize( stride_ * states.size() );
				const char* start = reinterpret_cast<const char*>( &states[0] );
				std::copy( start, start + states_.size(), states_.begin() );
			}
		}

		float time() const { return time_; }
		void setTime( float time ) { time_ = time; }

		const char* extendedState( size_t iState ) const
		{
			return &states_[iState*stride_];
		}

		char* extendedState( size_t iState )
		{
			return &states_[iState*stride_];
		}
		
		RigidDynamicState& state( size_t iState )
		{
			char* s = this->extendedState( iState );
			RigidDynamicState* rs = reinterpret_cast<RigidDynamicState*>( s );
			return *rs;
		}

		const RigidDynamicState& state( size_t iState ) const
		{
			const char* s = this->extendedState( iState );
			const RigidDynamicState* rs = reinterpret_cast<const RigidDynamicState*>( s );
			return *rs;
		}

		size_t stateCount() const
		{
			return states_.size() / stride_;
		}

		size_t stride() const
		{
			return stride_;
		}

		bool operator<( const State& other ) const
		{
			std::less<float> comparator;
			return comparator( this->time(), other.time() );
		}

		bool operator==( const State& other )
		{
			assert( stateCount() == other.stateCount() );
			assert( stride() == other.stride() );

			if( time() != other.time() )
				return false;

			// bitwise comparison
			for( unsigned int i = 0; i < states_.size(); ++i )
			{
				if( !(this->states_[i] == other.states_[i]) )
					return false;
			}

			return true;
		}

	private:
		std::vector<char> states_;
		size_t stride_;
		float time_;
	};

	Simulation( boost::shared_ptr<const Scene> scene );
	virtual ~Simulation();

	// we need to be able to add objects one at a time in order to represent
	//   'ghost' objects: those that are there only to detect interactions, 
	//   but which should never be simulated.
	virtual int addObject( ConstPhysicsObjectPtr object, bool ghost ) = 0;

    // For now you can only remove dynamic objects, since this is the capability that
    // we need during sampling
	virtual void removeObject( int dynamicId ) = 0;

	std::vector<ConstPhysicsObjectPtr> dynamicObjects() const;

	// for brevity, these state functions will only include
	//   dynamic objects:
	virtual Simulation::State getSimulationState() const = 0;
	virtual RigidDynamicState getSimulationState( size_t iObject ) const = 0;
	virtual void setSimulationState( const Simulation::State& state ) = 0;
	virtual void setSimulationState( size_t iObject, const RigidDynamicState& state ) = 0;

	virtual void step( const float stepSize ) = 0;
	virtual void multiStep( const float endTime, const float stepSize );

	virtual bool asleep() const;

	virtual void setDisabled( size_t iDynamicObject, bool disabled = true ) = 0;
	virtual bool getDisabled( size_t iDynamicObject ) const = 0;

	float time() const;
	void setTime( float time );

	virtual void perturbLinearVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation ) = 0;
	virtual void perturbAngularVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation ) = 0;
	virtual vl::Vec3f getLinearVelocity( size_t iDynamicObject ) const = 0;
	virtual vl::Vec3f getAngularVelocity( size_t iDynamicObject ) const = 0;

    struct ContactPoint
    {
        ContactPoint( const vl::Vec3f& p, const vl::Vec3f& n, const size_t ob1, const size_t ob2 )
            : position(p), normal(n), objects(ob1, ob2) {}

        vl::Vec3f position;
        vl::Vec3f normal;
        std::pair<int, int> objects;
    };
	virtual std::vector<ContactPoint> contactPoints() const;

	// These will be the objects that activate other objects:
	typedef STLEXT hash_set< const PhysicsObject*, hash_ptr<const PhysicsObject*> > PhysicsObjectSet;
	void setActiveObjects( PhysicsObjectSet activeObjects );
	bool isActive( const PhysicsObject* object ) const;

	bool ignoreCollisions( const PhysicsObject* left, const PhysicsObject* right ) const;
    bool connected( const PhysicsObject* left, const PhysicsObject* right ) const;

    void setSamplingProperties( size_t iDynamicObject, SamplingProperties props );
    SamplingProperties getSamplingProperties( size_t iDynamicObject ) const;

    virtual void dumpState( const char* filename ) const = 0;
    virtual void loadState( const char* filename ) = 0;

protected:
    boost::shared_ptr<RandomStream> randomStream_;
    
	// ignore collisions for object and all its children
	void addIgnoreCollisions( ConstSceneGraphElementPtr root );
	void addIgnoreCollisions( const PhysicsObject* left, const PhysicsObject* right );

    size_t addDynamicObject( ConstPhysicsObjectPtr object );
    int dynamicId( ConstPhysicsObjectPtr object ) const;

private:
	std::vector<ConstPhysicsObjectPtr> dynamicObjects_;
    
    typedef std::pair<const PhysicsObject*, const PhysicsObject*> PhysicsObjectPair;
	struct HashPhysicsObjectPair
		: stdext::hash_compare<PhysicsObjectPair>
	{
		size_t operator()( const PhysicsObjectPair& pr ) const;
		bool operator()( const PhysicsObjectPair& left, const PhysicsObjectPair& right ) const;
	};

	// we will want to ignore certain collisions, this is how we'll keep
	//   track of them:
	typedef STLEXT hash_set< PhysicsObjectPair, HashPhysicsObjectPair > PhysicsObjectPairSet;
	PhysicsObjectPairSet ignoreCollisionsSet_;

	PhysicsObjectSet activeObjects_;
	float time_;

    // We want to keep track of objects that are connected to each other
    // through a path of joints.  
    typedef STLEXT hash_map< const PhysicsObject*, size_t > ConnectedComponentMap;
    ConnectedComponentMap connectedComponents_;

    std::vector<SamplingProperties> samplingProperties_;
};

class SimulationFactory
{
public:
	virtual ~SimulationFactory();

	virtual SimulationPtr construct( boost::shared_ptr<const Scene> scene ) const = 0;
};

template <typename SimulationType>
class SimulationFactoryInst
	: public SimulationFactory
{
public:
	SimulationPtr construct( boost::shared_ptr<const Scene> scene ) const
	{
		SimulationPtr result( new SimulationType(scene) );
		return result;
	}
};

class Scene
{
	typedef STLEXT hash_map< std::string, SceneGraphElementPtr > ObjectNameMap;
	

public:
	Scene();
	~Scene();

	Scene( const Scene& other );

	// generate a unique name; takes as input a prefix of the form "blah" and finds the
	//   first name 'blah*' that is available.
	std::string uniqueObjectName( const std::string& prefix ) const;
	std::string uniqueCameraName( const std::string& prefix ) const;

	// this only adds the objects to the name lookup table.
	//   they must be made children of some other scene element
	//   to actually appear in the scene graph
	void addObject( SceneGraphElementPtr object );
	void addObjects( std::deque<SceneGraphElementPtr> objects );
	void removeObject( SceneGraphElementPtr object );
	void removeObjects( std::deque<SceneGraphElementPtr> objects );

	// retrieve an object by name, useful in certain cases
	SceneGraphElementPtr object( const std::string& name );
	ConstSceneGraphElementPtr object( const std::string& name ) const;

	// iterate over _all_ objects in the scene
	typedef boost::transform_iterator<
		std::select2nd<ObjectNameMap::value_type>, 
		ObjectNameMap::const_iterator> object_iter;
	object_iter begin_objects() const;
	object_iter end_objects() const;

	SceneGraphElementPtr root();
	ConstSceneGraphElementPtr root() const;

	size_t materialCount() const;
	MaterialPtr material( size_t iMat );
	ConstMaterialPtr material( size_t iMat ) const;
	MaterialPtr material( const std::string& name );
	ConstMaterialPtr material( const std::string& name ) const;
	void addMaterial( MaterialPtr material );
	void removeMaterial( MaterialPtr material );
	MaterialPtr defaultMaterial();

	BoundingBox3f bounds() const;

	// Export scene to file:
	MechModelObjectPtr toMechModel() const;

	TriangleMeshPtr triMesh( const std::string& filename );

	typedef std::deque<CameraWrapperPtr>::const_iterator camera_iterator;
	camera_iterator begin_cameras() const			{ return cameras_.begin(); }
	camera_iterator end_cameras() const				{ return cameras_.end(); }
	void addCameras( std::deque<CameraWrapperPtr> cameras );
	void addCamera( CameraWrapperPtr cameras );
	void removeCameras( std::deque<CameraWrapperPtr> cameras );
	void removeCamera( CameraWrapperPtr camera );
	CameraWrapperPtr camera( const std::string& name );
	CameraWrapperPtr camera();
	void setCamera(CameraWrapperPtr cam);

	void writeRibFile( const std::string& filename, 
		const std::deque<Simulation::State>& states, 
		const std::vector<ConstSceneGraphElementPtr>& stateObjects );

	// @todo figure out how to do animation
	void writeMayaAsciiFile( const std::string& filename,
		const std::vector<PiecewisePath>& paths, 
		size_t frameRate,
		const std::vector<ConstPhysicsObjectPtr>& stateObjects) const;

private:
	typedef NamedPtrComparator ObjectComparator;
	typedef NamedPtrComparator MaterialComparator;
	typedef NamedPtrComparator CameraComparator;

	boost::shared_ptr<TransformGroup> root_;

	ObjectNameMap objectNameMap_;

	std::deque< CameraWrapperPtr > cameras_;

	// this is the one that we'll use in evaluating all the metrics
	CameraWrapperPtr sceneCamera_;

	std::deque< MaterialPtr > materials_;
	MaterialPtr defaultMaterial_;

	typedef STLEXT hash_map< std::string, TriangleMeshPtr > TriMeshMap;
	TriMeshMap meshes_;
};

// utility functions:

// we will need to keep track of the copies of things so we can
//   run through and fix up all the joints, etc.
typedef STLEXT hash_map<const SceneGraphElement*, SceneGraphElementPtr, 
	hash_ptr<const SceneGraphElement*> > SceneGraphElementMap;

// clones the object plus its whole tree:
SceneGraphElementPtr cloneObject( ConstSceneGraphElementPtr object, Scene* scene, 
								 SceneGraphElementMap& oldToNewMap );

// fix up joints so they use the new objects rather than the old ones
void fixJoints( SceneGraphElementPtr current, const SceneGraphElementMap& oldToNewMap );

class ODESimulation
	: public Simulation
{
public:
	ODESimulation(boost::shared_ptr<const Scene> scene);
	~ODESimulation();

	int addObject( ConstPhysicsObjectPtr object, bool ghost );
    void removeObject( int dynamicId );

	typedef std::pair<vl::Vec3f, vl::Vec3f> PointNormal;
	std::vector<PointNormal> findEpsilonContacts( size_t iDynamicObject, float maxEpsilon );

	Simulation::State getSimulationState() const;
	RigidDynamicState getSimulationState( size_t iObject ) const;
	void setSimulationState( const Simulation::State& state );
	void setSimulationState( size_t iObject, const RigidDynamicState& state );

	void step(const float);
	void setDisabled( size_t iDynamicObject, bool disabled );
	bool getDisabled( size_t iDynamicObject ) const;

#ifdef SOUND
	std::vector<double> stepAudio(const Simulation::State& previousState, 
		float previousTimestep, size_t numSteps);

	// this one is for the precomputed examples; we will need to infer the
	// collision impulses from the instantaneous velocity changes and
	// use that do comptue some plausible impulses:
	std::vector<double> stepAudio(const Simulation::State& previousState, 
		const Simulation::State& nextState,
		size_t stateFrameRate,
		float epsilon,
		size_t numTimesteps );
#endif

	dBodyID odeBodyID( ConstPhysicsObjectPtr object ) const;

	void perturbLinearVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation );
	void perturbAngularVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation );

	vl::Vec3f getAngularVelocity( size_t iDynamicObject ) const;
	vl::Vec3f getLinearVelocity( size_t iDynamicObject ) const;

	std::vector<Simulation::ContactPoint> contactPoints() const;

    void dumpState( const char* filename ) const;
    void loadState( const char* filename );

protected:
	// we won't be able to build collision objects until the derived class is fully
	//   initialized, so we push scene-building into this function and require derived
	//   to call it.
	void initialize(boost::shared_ptr<const Scene> scene);
	void addTree( ConstSceneGraphElementPtr root );

	// in order to allow us to allow handling collisions using either ODE or Bullet Dynamics,
	//   we allow overloading here.  This function is required to fill up our contactGeoms_
	//   array with appropriate contact constraints.
	virtual void handleCollisions(bool backwards) = 0;

	// add a collision object and tie it to the specified bodyId
	virtual void addCollisionObject( 
		boost::shared_ptr<const PhysicsObject> object,
		int iDynamicId,
		dBodyID bodyId,
		const vl::Vec3f& centerOfMass,
		float padding,
		bool isGhost ) = 0;

	virtual void removeCollisionObject( int iDynamicId,
        dBodyID bodyId ) = 0;

	void addContactJoint( dBodyID b1, dBodyID b2,
		const dContact& geom );
private:
	std::vector<double> stepAudio( size_t numTimesteps );

	// check for disabled joints:
	void updateJoints(float time);
	RigidDynamicState state( size_t iObject ) const;
	void setState( size_t iObject, const RigidDynamicState& s );

	dWorldID worldId_;

	// we will keep track of contacts ourselves
	// rather than entrusting them to a joint group,
	// since ODE doesn't seem to provide any way to 
	// iterate over all the joints in a joint group
	std::deque<dJointID> contactJoints_;
	std::deque<dContactGeom> contactGeoms_;
	std::deque<dJointFeedback> contactJointFeedback_;

	std::deque<dBodyID> dynamicBodies_;
    boost::dynamic_bitset<> ghostBodies_;

    // These bodies need to get fixed up at the end of every timestep
    // Note that this is a bit of a hack
    // @todo make this better.
    std::deque<dBodyID> planeJointBodies_;

#ifdef SOUND
	typedef STLEXT hash_map< const PhysicsObject*, 
		boost::shared_ptr<AudioGenerator>,
		hash_ptr<const PhysicsObject*> > ObjectToAudioHash;
	ObjectToAudioHash audioGenerators_;
#endif

	Simulation::State initialState_;

	// center of mass offsets for various objects:
	// We need this here because ODE does not support having a nonzero 
	// center of mass, so it needs to be corrected outside the solver.
	std::deque<vl::Vec3f> offsets_;

	typedef boost::shared_ptr< const Joint > ConstJointPtr;
	typedef std::pair< ConstJointPtr, dJointID > JointWithId;
	std::deque<JointWithId> activeJoints_;
	std::deque<ConstJointPtr> inactiveJoints_;

protected:
	dBodyID odeBodyId( size_t iObject );

	// Want to be able to quickly locate the material
	//   given a geom ID.  Also this will prevent materials
	//   from getting changed "under the simulation's nose"
	//   so to speak.
	struct GeomInfo
	{
		GeomInfo(const Material& mat)
            : originalObject(0), iDynamicId(-1), material(mat), ghost(false), 
              lastCollisionTime(boost::numeric::bounds<float>::lowest()), 
              lastSlidingTime(boost::numeric::bounds<float>::highest()), 
              lastRollingTime(boost::numeric::bounds<float>::highest()) 
        {}

		const PhysicsObject* originalObject;
		int iDynamicId;
		Material material;
		bool ghost;

        mutable float lastCollisionTime;
        mutable float lastSlidingTime;
        mutable float lastRollingTime;
        mutable vl::Vec3d lastFrictionDir;
	};

	friend class MyBulletCollisionDispatcher;

    // just need to keep track of this:
    double stepSize_;
    bool backwards_;

    typedef STLEXT hash_map<dBodyID, size_t> BodyToIntHash;
    BodyToIntHash bodyToId_;
};

class ODENativeSimulation
	: public ODESimulation
{
public:
	ODENativeSimulation(boost::shared_ptr<const Scene> scene);
	~ODENativeSimulation();

protected:
	void handleCollisions(bool backwards);
	void addCollisionObject( 
		boost::shared_ptr<const PhysicsObject> object,
		int iDynamicId,
		dBodyID bodyId,
		const vl::Vec3f& centerOfMass,
		float padding,
		bool isGhost );
    void removeCollisionObject( int iDynamicId, dBodyID bodyId );

protected:
    // HACK HACK HACK should make this private
	dSpaceID spaceId_;

private:
	friend void potentialCollision( void* data, dGeomID o1, dGeomID o2 );
	void checkContact( dGeomID o1, dGeomID o2 );

	boost::object_pool<GeomInfo> geomInfoPool_;

    // bug in the quadSpace implementation means that we need to remove all these
    // from their space BEFORE deleting the space.
	std::deque<dGeomID> allGeoms_;
	std::deque<dSpaceID> allSpaces_;

#if dTRIMESH_ENABLED
	std::deque<ODETriMeshWrapperPtr> meshes_;
#endif

	// we need this so we can go through and set the prev. transform
	//   for ODE's trimesh handling:
	std::deque<dGeomID> triMeshGeoms_;
};

class PenaltyForceSimulation
	: public ODENativeSimulation
{
public:
	PenaltyForceSimulation(boost::shared_ptr<const Scene> scene);
	~PenaltyForceSimulation();

protected:
	void handleCollisions(bool backwards);

private:
	friend void potentialCollisionPenalty( void* data, dGeomID o1, dGeomID o2 );
	void checkContact( dGeomID o1, dGeomID o2, bool backwards );
};

#ifdef USE_BULLET
// uses ODE for dynamics but the Bullet library for collision handling
class ODEBulletSimulation
	: public ODESimulation
{
public:
	ODEBulletSimulation(boost::shared_ptr<const Scene> scene);
	~ODEBulletSimulation();

protected:
	void handleCollisions(bool backwards);
	void addCollisionObject( 
		boost::shared_ptr<const PhysicsObject> object,
		int iDynamicId,
		dBodyID bodyId,
		const vl::Vec3f& centerOfMass,
		float padding,
		bool isGhost );
    void removeCollisionObject( int dynamicId, dBodyID bodyId );

private:
	// bullet stuff
	boost::scoped_ptr<btCollisionDispatcher> collisionDispatcher_;
	boost::scoped_ptr<btBroadphaseInterface> broadPhase_;

	// keep the collision shapes, for deletion/cleanup
	boost::ptr_deque<btCollisionAlgorithmCreateFunc> collisionDispatcherFunctions_;

	boost::ptr_deque<btCollisionObject> staticRigidBodies_;
	boost::ptr_deque<btCollisionObject> dynamicRigidBodies_;

	boost::ptr_deque<btCollisionShape> collisionShapes_;

	// it is important that this be the first thing deleted
	boost::scoped_ptr<btCollisionWorld> collisionWorld_;

	boost::object_pool<GeomInfo> geomInfoPool_;

    btDefaultCollisionConfiguration collisionConfig_;
};

typedef ODEBulletSimulation ODEDefaultSimulationType;
#else
typedef ODENativeSimulation ODEDefaultSimulationType;
#endif


#ifdef USE_NEWTON
class NewtonSimulation
	: public Simulation
{
public:
	NewtonSimulation(boost::shared_ptr<const Scene> scene);
	~NewtonSimulation();

	Simulation::State getSimulationState() const;
	void setSimulationState( const Simulation::State& state );

	void step( const float stepSize );

	void setDisabled( size_t iDynamicObject, bool disabled = true );
	bool getDisabled( size_t iDynamicObject ) const;


private:
	NewtonWorld* world_;
	std::vector<NewtonBody*> dynamicBodies_;

	// since Newton dynamic requires a diagonalized mass matrix, we
	//   will need to apply some transformations to get things to act
	//   correctly
	std::vector<dFloatMat4> offsetTransformations_;
};
#endif

#ifdef USE_NOVODEX
class NxSimulation
	: public Simulation, public NxUserTriggerReport
{
public:
	NxSimulation(NxPhysicsSDKPtr sdk, boost::shared_ptr<const Scene> scene);
	~NxSimulation();

	NxScenePtr nxScene();
	NxPhysicsSDKPtr nxSDK();

	Simulation::State getSimulationState() const;
	RigidDynamicState getSimulationState( size_t iObject ) const;
	void setSimulationState( const Simulation::State& state );
	void setSimulationState( size_t iObject, const RigidDynamicState& state );

	void perturbLinearVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation );
	void perturbAngularVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation );

	vl::Vec3f getAngularVelocity( size_t iDynamicObject ) const;
	vl::Vec3f getLinearVelocity( size_t iDynamicObject ) const;

	void step( const float stepSize );
	void multiStep( const float endTime, const float stepSize );
	void setDisabled( size_t iDynamicObject, bool disabled );
	bool getDisabled( size_t iDynamicObject ) const;

	NxActorPtr nxActor( ConstPhysicsObjectPtr object ) const;

	int addObject( ConstPhysicsObjectPtr object, bool ghost );

	void onTrigger(NxShape& triggerShape, NxShape& otherShape, NxTriggerFlag status);

private:
	NxScenePtr nxScene_;
	NxPhysicsSDKPtr physicsSDK_;

	typedef boost::shared_ptr<NxMaterial> NxMaterialPtr;
	typedef std::map< ConstMaterialPtr, NxMaterialPtr > NxMaterialMap;
	NxMaterialMap nxMaterials_;

	typedef boost::shared_ptr<NxActor> NxActorPtr;
	std::deque<NxActorPtr> staticActors_;
	std::deque<NxActorPtr> dynamicActors_;

	std::deque<NxJointPtr> joints_;

	typedef STLEXT hash_map<const PhysicsObject*, NxActorPtr> PhysicsObjectToNxActorMap;
	PhysicsObjectToNxActorMap physicsObjectToActorMap_;

	void addTree( ConstSceneGraphElementPtr root, bool selfCollisions, std::deque<NxActorPtr>& actors = std::deque<NxActorPtr>() );

	// Since we can't use the machinery built into Novodex for tracking
	//   which objects have been activated, we will have to maintain our
	//   own list.
	std::vector<bool> active_;

	// We will store information about an object in this structure which
	//   will then get assigned to an NxActor's user pointer:
	struct GeomInfo
	{
		GeomInfo()
			: originalObject(0), iDynamicId(-1), ghost(false) {}

		const PhysicsObject* originalObject;
		int iDynamicId;
		bool ghost;
	};

	boost::object_pool<GeomInfo> geomInfoPool_;

	friend class ContactReport;


	boost::ptr_deque<std::string> nameCache_;
};

class NovodexSimulationFactory
	: public SimulationFactory
{
public:
	NovodexSimulationFactory(NxPhysicsSDKPtr sdk);
	~NovodexSimulationFactory();

	SimulationPtr construct( boost::shared_ptr<const Scene> scene ) const;

private:
	NxPhysicsSDKPtr sdk_;
};

#endif

#ifdef USE_BULLET
class BulletSimulation
	: public Simulation
{
public:
	BulletSimulation( boost::shared_ptr<const Scene> scene );
	~BulletSimulation();

	int addObject( ConstPhysicsObjectPtr object, bool ghost );
    void removeObject( int dynamicId );

	Simulation::State getSimulationState() const;
	RigidDynamicState getSimulationState( size_t iObject ) const;
	void setSimulationState( const Simulation::State& state );
	void setSimulationState( size_t iObject, const RigidDynamicState& state );

	void step( const float stepSize );
	void multiStep( const float endTime, const float stepSize );

	void setDisabled( size_t iDynamicObject, bool disabled );
	bool getDisabled( size_t iDynamicObject ) const;

	void perturbLinearVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation );
	void perturbAngularVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation );

	vl::Vec3f getAngularVelocity( size_t iDynamicObject ) const;
	vl::Vec3f getLinearVelocity( size_t iDynamicObject ) const;

    virtual void dumpState( const char* filename ) const { assert(false); }
    virtual void loadState( const char* filename ) { assert(false); }

private:
	// bullet stuff
	boost::scoped_ptr<btCollisionDispatcher> collisionDispatcher_;
	boost::scoped_ptr<btConstraintSolver> constraintSolver_;
	boost::scoped_ptr<btBroadphaseInterface> broadPhase_;

	// keep the collision shapes, for deletion/cleanup
	boost::ptr_deque<btCollisionAlgorithmCreateFunc> collisionDispatcherFunctions_;

	boost::ptr_deque<btRigidBody> staticRigidBodies_;
	boost::ptr_deque<btRigidBody> dynamicRigidBodies_;

	boost::ptr_deque<btMotionState> motionStates_;
	boost::ptr_deque<btCollisionShape> collisionShapes_;
	std::vector<btMat4> offsetTransformations_;

	// it is important that this be the first thing deleted
	boost::scoped_ptr<btDynamicsWorld> dynamicsWorld_;

    btDefaultCollisionConfiguration collisionConfig_;
};
#endif

// This kind of simulation only ever returns the same SimulationState and 
// contact points.  This is so we can load in a particular simulation w/ 
// contacts and visualize it easily.  
class DummySimulation
    : public Simulation
{
public:
    DummySimulation( boost::shared_ptr<const Scene> scene, 
        const Simulation::State& state,
        const std::vector<Simulation::ContactPoint>& contactPoints );
	virtual ~DummySimulation();

    virtual int addObject( ConstPhysicsObjectPtr object, bool ghost ) { assert( false ); return 0; }
    virtual void removeObject( int dynamicId ) { assert( false ); }

	// for brevity, these state functions will only include
	//   dynamic objects:
    virtual Simulation::State getSimulationState() const                              { return this->state_; }
    virtual RigidDynamicState getSimulationState(size_t iObject) const                { return this->state_.state(iObject); }
    virtual void setSimulationState( const Simulation::State& state )                 { assert( false ); }
    virtual void setSimulationState( size_t iObject, const RigidDynamicState& state ) { assert( false ); }

    virtual std::vector<ContactPoint> contactPoints() const                           { return this->contactPoints_; }

    virtual void step( const float stepSize )                                         {}

    virtual void setDisabled( size_t iDynamicObject, bool disabled = true )           { assert( false ); }
    virtual bool getDisabled( size_t iDynamicObject ) const                           { return false; }

    virtual void perturbLinearVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation )  { assert( false ); }
    virtual void perturbAngularVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation ) { assert( false ); }
    virtual vl::Vec3f getLinearVelocity( size_t iDynamicObject ) const;
    virtual vl::Vec3f getAngularVelocity( size_t iDynamicObject ) const;

    virtual void dumpState( const char* filename ) const {}
    virtual void loadState( const char* filename ) {}

private:
    Simulation::State state_;
    std::vector<Simulation::ContactPoint> contactPoints_;
};

template <typename T>
void setAttr( std::ostream& ofs, const char* attributeName, const T& value );

std::deque<ConstSceneGraphElementPtr> symmetricDiffScene(
	const Scene& origScene,
	const Scene& changedScene );

std::vector<ConstPhysicsObjectPtr> physicsObjects( const Scene& scene );

} // namespace planning


#endif

