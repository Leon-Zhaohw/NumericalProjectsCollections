#ifndef __PHYSICSFWD_H__

#define __PHYSICSFWD_H__

#include "twigg/primitive.h"
#include "twigg/AABBTree.h"

#include <vector>

class ObjFile;

#ifdef USE_NEWTON
#include "newton.h"
#endif

#ifdef USE_BULLET
#include "btBulletDynamicsCommon.h"
#include "btBulletCollisionCommon.h"
#endif

namespace planning {

// Forward declarations

typedef boost::shared_ptr<ObjFile> ObjFilePtr;

class Simulation;
class SimulationFactory;
class NxSimulation;
class ODESimulation;
typedef boost::shared_ptr<Simulation> SimulationPtr;

class Constraint;
typedef boost::shared_ptr<const Constraint> ConstraintPtr;

class CameraWrapper;
typedef boost::shared_ptr<CameraWrapper> CameraWrapperPtr;
typedef boost::shared_ptr<const CameraWrapper> ConstCameraWrapperPtr;

class SceneGraphElement;
typedef boost::shared_ptr<SceneGraphElement> SceneGraphElementPtr;
typedef boost::shared_ptr<const SceneGraphElement> ConstSceneGraphElementPtr;

class TransformedSceneGraphElement;
typedef boost::shared_ptr<TransformedSceneGraphElement> TransformedElementPtr;
typedef boost::shared_ptr<const TransformedSceneGraphElement> ConstTransformedElementPtr;

class PhysicsObject;
typedef boost::shared_ptr<PhysicsObject> PhysicsObjectPtr;
typedef boost::shared_ptr<const PhysicsObject> ConstPhysicsObjectPtr;

class Joint;
typedef boost::shared_ptr<const Joint> JointPtr;

class Material;
typedef boost::shared_ptr<Material> MaterialPtr;
typedef boost::shared_ptr<const Material> ConstMaterialPtr;

class Objective;
typedef boost::shared_ptr<const Objective> ObjectivePtr;

class Plug;
typedef boost::shared_ptr<Plug> PlugPtr;
typedef std::pair< PhysicsObjectPtr, PlugPtr > ObjectPlug;

class MechModelObject;
typedef boost::shared_ptr<MechModelObject> MechModelObjectPtr;

class Scene;
typedef boost::shared_ptr<Scene> ScenePtr;
typedef boost::shared_ptr<const Scene> ConstScenePtr;

class PMapWrapper;
typedef boost::shared_ptr<PMapWrapper> PMapPtr;

class TriangleMesh;
typedef boost::shared_ptr<const TriangleMesh> TriangleMeshPtr;

class SimulationTree;
typedef boost::shared_ptr<SimulationTree> SimulationTreePtr;
typedef boost::shared_ptr<const SimulationTree> ConstSimulationTreePtr;

class RibObject;
class RibNamedBlock;
typedef boost::shared_ptr<const RibObject> RibObjectPtr;

#ifdef SOUND
struct SquashingCubesModes;
class AudioGenerator;
#endif

// This is the part of the state vector that we can use a 
//   standard-type metric on:
typedef boost::array<float, 9> EuclideanState;

#ifdef USE_NEWTON
typedef vl::Mat4f dFloatMat4;
typedef vl::Mat3f dFloatMat3;
typedef vl::Vec4f dFloatVec4;
typedef vl::Vec3f dFloatVec3;
#endif

#ifdef BT_USE_DOUBLE_PRECISION
typedef vl::Mat4d btMat4;
typedef vl::Mat3d btMat3;
typedef vl::Vec4d btVec4;
typedef vl::Vec3d btVec3;
#else
typedef vl::Mat4f btMat4;
typedef vl::Mat3f btMat3;
typedef vl::Vec4f btVec4;
typedef vl::Vec3f btVec3;
#endif

// We will pass one of these in to the various metrics to compute
//   view-dependent metrics.  
class CameraSource
{
public:
	virtual ~CameraSource() {}

	virtual vl::Mat4f projectionMatrix( float time ) const = 0;
	virtual vl::Mat4f modelViewMatrix( float time ) const = 0;
	virtual vl::Mat4f inverseModelViewMatrix( float time ) const = 0;
};

/*
struct PrimGroup
{
	GLenum glPrimitiveType() const
	{
		switch(type)
		{
		case LIST:		return GL_TRIANGLES;
		case STRIP:		return GL_TRIANGLE_STRIP;
		case FAN:		return GL_TRIANGLE_FAN;
		}

		throw Exception( "Shouldn't reach here." );
	}

	enum Type
	{
		LIST = 0,
		STRIP = 1,
		FAN = 2,
	} type;

	std::vector<unsigned short> indices;
};
*/

// will be used for things that take a while to run:
class ProgressUpdater
{
public:
	virtual ~ProgressUpdater() {}
	virtual void create( const std::string& text, int maxValue ) = 0;
	virtual void update( int value ) = 0;
};

typedef boost::shared_ptr<ProgressUpdater> ProgressUpdaterPtr;

// An action is something that is undoable
class Action
{
public:
	virtual ~Action();
	virtual void doIt() = 0;
	virtual void undoIt() = 0;
	virtual bool changesScene() const = 0;
	virtual bool changesHierarchy() const = 0;
};

typedef boost::shared_ptr<Action> ActionPtr;

} // namespace planning

#endif
