#ifdef WIN32
#pragma once
#endif

#ifndef __PHYSICSOBJECT_H__
#define __PHYSICSOBJECT_H__

#include "novodexWrappers.h"
#include "odeWrappers.h"
#include "physicsFwd.h"
#include "quat.h"

#include "twigg/volumeIntegrals.h"
#include "twigg/boundingbox.h"
#include "twigg/AABBTree.h"
#include "twigg/primitive.h"

#include "loki/AssocVector.h"

#include <string>

#ifdef CONDOR_SAMPLING
class MWRMComm;
#endif

class GLTexture;
class GLDisplayList;

namespace planning {

class BoundedVoxelGrid;

// For the sound stuff, we want to be able to dump out
//   a special mech model that only includes properties
//   important for audio
enum MechModelType
{
	MECH_MODEL_FULL,
	MECH_MODEL_AUDIO_ROOT,
	MECH_MODEL_AUDIO,
};

class AttributeException
	: public Exception
{
public:
	AttributeException( const std::string& attribute, char* message )
		: Exception(message), message_( message ) {}

	std::string attribute() const
	{
		return this->attribute_;
	}

	std::string message() const
	{
		return this->message_;
	}

private:
	std::string attribute_;
	const char* message_;
};

class Material
	: public Named
{
public:
	Material( const std::string& name );
	Material( const std::string& name, float density, float dynamicFriction, float staticFriction, float restitution );
	~Material();

	bool equals( const Material& other ) const;

    float density() const;
	float dynamicFriction() const;
	float staticFriction() const;
	float restitution() const;

	float youngsModulus() const;
	float poissonRatio() const;
	float raleighAlpha() const;
	float raleighBeta() const;

	vl::Vec3f color() const;
	bool coulombFriction() const;

	void setDensity( float density );
	void setDynamicFriction( float friction );
	void setStaticFriction( float friction );
	void setRestitution( float restitution );
	void setColor( const vl::Vec3f& color );
	void setCoulombFriction( bool value );

	void setRaleighAlpha( float value );
	void setRaleighBeta( float value );
	void setYoungsModulus( float value );
	void setPoissonRatio( float value );

	MechModelObjectPtr toMechModel(MechModelType type = MECH_MODEL_FULL) const;

private:
	float density_;
	float dynamicFriction_;
	float staticFriction_;
	float restitution_;
	bool coulombFriction_;

	// sound params
	float youngsModulus_;
	float poissonRatio_;
	float raleighAlpha_;
	float raleighBeta_;

	vl::Vec3f color_;
};

// An ordering that preserves only the sound-sensitive properties
//   of Young's modulus and Poisson ratio
struct MaterialSoundPropertiesOrder
{
	bool operator()( const Material& left, const Material& right )
	{
		float propsLeft[] = { left.youngsModulus(), left.poissonRatio() };
		float propsRight[] = { right.youngsModulus(), right.poissonRatio() };
		return std::lexicographical_compare( propsLeft, propsLeft+2,
			propsRight, propsRight+2 );
	}

	bool operator()( const boost::shared_ptr<const Material>& left, 
		const boost::shared_ptr<const Material>& right )
	{
		return (*this)(*left, *right);
	}

	bool operator()( const Material* left, const Material* right )
	{
		return (*this)(*left, *right);
	}
};

struct MaterialSoundPropertiesEqual
{
	bool operator()( const Material& left, const Material& right )
	{
		return left.youngsModulus() == right.youngsModulus()
			&& left.poissonRatio() == right.poissonRatio();
	}

	bool operator()( const boost::shared_ptr<const Material>& left, 
		const boost::shared_ptr<const Material>& right )
	{
		return (*this)(*left, *right);
	}

	bool operator()( const Material* left, const Material* right )
	{
		return (*this)(*left, *right);
	}
};

class MechModelDeclaredObject;

struct RandomizedProperty
{
	RandomizedProperty() {}

	RandomizedProperty( const std::string& nm )
		: name(nm) {}

	std::string name;
	float low;
	float high;
};

class SceneGraphElement
	: public Named, public boost::enable_shared_from_this<SceneGraphElement>
{
public:
	struct RenderState
	{
		RenderState()
			:	polygonMode( FILL ), 
				wireframeColor( vl::vl_0 ), 
				alpha( 1.0 ), 
				overrideColor(false), 
				newColor( vl::vl_1 ),
				renderHulls(false)
		{
		}

		enum PolygonMode
		{
			POINTS,
			LINES,
			FILL,
			NONE
		} polygonMode;

		vl::Vec4f wireframeColor;
		float alpha;

		bool overrideColor;
		vl::Vec3f newColor;

		bool renderHulls;
	};

private:
	template <class SceneGraphElementPtrType>
	class child_iterator_base
		:	public boost::iterator_facade<
				child_iterator_base<SceneGraphElementPtrType>,
				SceneGraphElementPtrType,
				boost::bidirectional_traversal_tag,
				SceneGraphElementPtrType>
	{
		struct enabler {};

	public:
		child_iterator_base()
			: child_() {}

		explicit child_iterator_base(SceneGraphElementPtrType p)
			: child_(p) {}

		template <class OtherSceneGraphElementPtrType>
		child_iterator_base( child_iterator_base<OtherSceneGraphElementPtrType> const& other,
			typename boost::enable_if<
				boost::is_convertible<OtherSceneGraphElementPtrType, SceneGraphElementPtrType>,
				enabler >::type = enabler() )
			: child_(other.child_) {}

	private:
		friend class boost::iterator_core_access;
		template <class> friend class child_iterator_base;

		template <class OtherValue>
		bool equal(child_iterator_base<OtherValue> const& other) const
		{ 
			return this->child_ == other.child_;
		}

		void increment()
		{
			child_ = child_->nextChild_;
		}

		void decrement()
		{
			child_ = child_->prevChild_;
		}

		SceneGraphElementPtrType dereference() const
		{
			return this->child_;
		}

		SceneGraphElementPtrType child_;
	};

public:
	typedef child_iterator_base<SceneGraphElementPtr> child_iterator;
	typedef child_iterator_base<ConstSceneGraphElementPtr> const_child_iterator;

public:
	SceneGraphElement( const std::string& name );

	virtual ~SceneGraphElement();

	virtual bool allowUniformScale() const		{ return true; }
	virtual bool allowNonUniformScale() const	{ return true; }
	virtual bool allowShear() const				{ return false; }

	void setParent( SceneGraphElementPtr parent );
	boost::weak_ptr<SceneGraphElement> parent();
	boost::weak_ptr<const SceneGraphElement> parent() const;

	// adds the given child just before the element specified, or at
	//   end if not specified
	void addChild( SceneGraphElementPtr child, SceneGraphElementPtr after = SceneGraphElementPtr() );

	// returns the object just after the removed object; that is,
	//   addChild( child, removeChild(child) ) is the same as a noop
	SceneGraphElementPtr removeChild( SceneGraphElementPtr child );

	child_iterator children_begin();
	child_iterator children_end();

	const_child_iterator children_begin() const;
	const_child_iterator children_end() const;

	vl::Mat4f transform( const char* state = 0 ) const;

	bool visible() const;
	void setVisible( bool visible );

	bool selfCollisions() const;
	void setSelfCollisions( bool value );

	virtual size_t stateSize() const            { return 0; }

	boost::shared_ptr<MechModelDeclaredObject> toMechModel(MechModelType type = MECH_MODEL_FULL) const;
	virtual boost::shared_ptr<MechModelDeclaredObject> toMechModelLocal(MechModelType type) const;

	boost::shared_ptr<RibNamedBlock> toRibObject(const std::string& prefix, 
		const std::string& timesFilename,
		const std::vector<const char*>& states) const;
	void toMayaAscii( std::ofstream& ofs, const std::string& materialName, const std::string& prefix = std::string() ) const;

	// this one clones the whole tree:
	SceneGraphElementPtr clone() const;

	// This returns the bounding box for _just this object_,
	//   not its children
	virtual BoundingBox3f bounds(const vl::Mat4f& transform) const = 0;
	BoundingBox3f bounds(const char* state = 0) const;

	bool intersect(const Ray3d& ray, 
		double& t, 
		const vl::Mat4f& transform, 
		const char* state = 0) const;
	bool intersect(const BoundingBox2d& box, 
		const ScreenSpaceConverter& converter, 
		const vl::Mat4f& transform,
		const char* state = 0) const;

	void render(const RenderState& rs, 
		const vl::Mat4f& transform, 
		const char* state = 0) const;

	virtual void renderPoints( const std::deque<unsigned int>& selected, 
		const vl::Mat4f& transform, 
		const char* state = 0 ) const;
	virtual std::pair<bool, unsigned int> selectPoint( const vl::Vec2d& screenPos, 
		const ScreenSpaceConverter& converter, 
		const vl::Mat4f& transform,
		const char* state = 0) const;

	typedef SceneGraphElement* TransformChangedTarget;
	void addTransformChangedTarget( TransformChangedTarget target );
	void removeTransformChangedTarget( TransformChangedTarget target );

	virtual std::vector<std::string> floatAttributes() const;
	virtual void setAttribute( const std::string& attributeName, float value );
	virtual float getAttribute( const std::string& attributeName ) const;

	virtual std::vector<std::string> boolAttributes() const;
	virtual void setBoolAttribute( const std::string& attributeName, bool value );
	virtual bool getBoolAttribute( const std::string& attributeName ) const;

	bool equals( const SceneGraphElement& other, const SceneGraphElement* stopRecursingAt = 0 ) const;

	std::vector<RandomizedProperty> randomizedProperties() const;
	bool isRandomized( const std::string& propertyName ) const;
	void setRandomized( const RandomizedProperty& property );
	RandomizedProperty getRandomized( const std::string& propertyName ) const;
	void clearRandomized( const std::string& propertyName );

	virtual bool equalsLocal( const SceneGraphElement& other ) const;

	virtual vl::Mat4f localTransform( const char* state ) const = 0;

	virtual vl::Vec3f furthestPoint( const vl::Vec3f& point ) const;
protected:
	SceneGraphElement( const SceneGraphElement& other );
	SceneGraphElement& operator=( const SceneGraphElement& other );

	void updateScale();
	void updateTransform();
	void updateChildren();

	virtual void updateScaleLocal();
	virtual void updateTransformLocal();
	virtual void updateChildrenLocal();

	virtual boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const = 0;

	virtual vl::Vec3f color() const;

	virtual void renderLocalFill(
		float alpha, 
		const char* state) const = 0;
	virtual void renderLocalWire(
		const vl::Vec4f& wireframeColor, 
		const char* state) const = 0;
	virtual void renderLocalHulls(
		const float alpha, 
		const char* statePtr) const;


	virtual RibObjectPtr ribObject() const = 0;
	virtual void toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const;

	virtual bool intersectLocal(
		const Ray3d& ray, 
		double& t, 
		const char* state) const = 0;
	virtual bool intersectLocal(
		const BoundingBox2d& box, 
		const ScreenSpaceConverter& converter, 
		const char* state) const = 0;

	// this is just to clone me, ignoring the rest of 
	//   the tree
	virtual SceneGraphElement* cloneLocal() const = 0;

private:
	void checkConsistency() const;

	// we'll only keep a weak pointer to our parent to avoid
	//   resource leaks due to loops
	boost::weak_ptr<SceneGraphElement> parent_;

	SceneGraphElementPtr firstChild_;
	SceneGraphElementPtr lastChild_;

	// we'll make this a doubly-linked list to
	//   make removals very fast:
	SceneGraphElementPtr nextChild_;

	// however, back pointers are always weak_ptrs to 
	//   not screw up reference counts
	boost::weak_ptr<SceneGraphElement> prevChild_;

	std::deque< TransformChangedTarget > transformChangedRegistry_;
	std::deque<RandomizedProperty> randomizedProperties_;

	bool visible_;
	bool selfCollisions_;

	friend class ReparentAction;
};

Quaternion<float> eulerToQuat( const vl::Vec3f& angles );
vl::Vec3f quatToEuler( const Quaternion<float>& quat );
vl::Mat3f eulerToMat( const vl::Vec3f& angles );
vl::Vec3f matToEuler( const vl::Mat3f& mat );


// this uses a basic transform: Euler angles for simplicity
class TransformedSceneGraphElement
	: public SceneGraphElement
{
public:
	TransformedSceneGraphElement(const std::string& name);
	virtual ~TransformedSceneGraphElement();

	vl::Vec3f scale() const;
	void setScale( const vl::Vec3f& scale, bool temporary = false );

	vl::Vec3f translate() const;
	void setTranslate( const vl::Vec3f& translate, bool temporary = false );

	vl::Vec3f rotate() const;
	void setRotate( const vl::Vec3f& rotate, bool temporary = false );

	vl::Vec3f pivot() const;
	void setPivot( const vl::Vec3f& pivot, bool temporary = false );

	boost::shared_ptr<MechModelDeclaredObject> toMechModelLocal(MechModelType type) const;

	std::vector<std::string> floatAttributes() const;
	void setAttribute( const std::string& attributeName, float value );
	float getAttribute( const std::string& attributeName ) const;

	vl::Mat4f localTransform( const char* state ) const;

protected:
	bool equalsLocal( const SceneGraphElement& other ) const;

private:
	// order of transform is scale, then rotate, then translate

	vl::Vec3f translate_;

	// euler angles:
	vl::Vec3f rotate_;

	vl::Vec3f scale_;

	vl::Vec3f pivot_;
};

// a basic group
class TransformGroup
	: public TransformedSceneGraphElement
{
public:
	TransformGroup( const std::string& name );
	virtual ~TransformGroup();

protected:
	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;

	void renderLocalFill(
		float alpha, 
		const char* state) const;
	void renderLocalWire(
		const vl::Vec4f& wireframeColor, 
		const char* state) const;

	RibObjectPtr ribObject() const;

	bool intersectLocal(
		const Ray3d& ray, 
		double& t, 
		const char* state) const;
	bool intersectLocal(
		const BoundingBox2d& box, 
		const ScreenSpaceConverter& converter, 
		const char* state) const;

	BoundingBox3f bounds( const vl::Mat4f& transform ) const;

	TransformGroup* cloneLocal() const;
};

class RigidStaticState
{
public:
	typedef Quaternion<float>		OrientationType;

	RigidStaticState( const vl::Vec3f& position, const OrientationType& orientation );
	RigidStaticState();

	vl::Vec3f position() const				{ return this->position_; }
	OrientationType orientation() const	{ return this->orientation_; }

	void setPosition( const vl::Vec3f& position )
	{
		this->position_ = position;
	}

	void setOrientation( const OrientationType& orientation )
	{
		this->orientation_ = orientation;
	}

#ifdef CONDOR_SAMPLING
	RigidStaticState( MWRMComm * RMC );
	void write( MWRMComm * RMC ) const;
#endif

private:
	vl::Vec3f		position_;
	OrientationType	orientation_;
};

class RigidDynamicState
	: public RigidStaticState
{
public:
	RigidDynamicState( const vl::Vec3f& position, 
		const OrientationType& orientation, 
		const vl::Vec3f& linearVelocity, 
		const vl::Vec3f& angularVelocity );
	RigidDynamicState(const RigidStaticState& state);
	RigidDynamicState();

	vl::Vec3f linearVelocity() const		{ return linearVelocity_; }
	vl::Vec3f angularVelocity() const		{ return angularVelocity_; }

	void setLinearVelocity( const vl::Vec3f& linearVelocity )
	{
		linearVelocity_ = linearVelocity;
	}

	void setAngularVelocity( const vl::Vec3f& angularVelocity )
	{
		angularVelocity_ = angularVelocity;
	}

private:
	vl::Vec3f		linearVelocity_;
	vl::Vec3f		angularVelocity_;
};

struct ODEGeomResult
{
	ODEGeomResult() {}
	ODEGeomResult( dGeomID id )
		: geoms(1, id) {}

	std::vector<dGeomID> geoms;
#if dTRIMESH_ENABLED
	std::vector<ODETriMeshWrapperPtr> triMeshes;
#endif
};


class PhysicsObject
	: public TransformedSceneGraphElement
{
public:
	PhysicsObject( const std::string& name, boost::shared_ptr<Material> material );
	virtual ~PhysicsObject();

	// create a geom in the given space:
	ODEGeomResult odeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const;
#ifdef USE_NEWTON
	NewtonCollision* newtonCollision( const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const;
#endif
#ifdef USE_BULLET
	virtual btCollisionShape* bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const;
#endif

	bool isStatic() const;
	void setStatic( bool isStatic = true );

	bool collides() const;
	void setCollides( bool collides );

	bool hasMass() const;
	void setHasMass( bool hasMass );

	bool hasAudio() const;
	void setHasAudio( bool hasAudio );

	boost::shared_ptr<Material> material();
	ConstMaterialPtr material() const;
	void setMaterial( boost::shared_ptr<Material> material );

	boost::shared_ptr<MechModelDeclaredObject> toMechModelLocal(MechModelType type) const;
	virtual dMass getOdeMassProps() const = 0;
#ifdef USE_NOVODEX
	virtual std::vector< boost::shared_ptr<NxShapeDesc> > getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const = 0;
#endif

	vl::Vec3f color() const;

	vl::Vec3f linearVelocity() const;
	void setLinearVelocity( const vl::Vec3f& linVel );

	vl::Vec3f angularVelocity() const;
	void setAngularVelocity( const vl::Vec3f& angVel );

	vl::Mat4f localTransform( const char* state ) const;

	// this returns just the rigid part of the transform, excluding any scales, etc.
	vl::Mat4f rigidTransform( const char* state = 0 ) const;

	std::vector<std::string> floatAttributes() const;
	void setAttribute( const std::string& attributeName, float value );
	float getAttribute( const std::string& attributeName ) const;

	std::vector<std::string> boolAttributes() const;
	void setBoolAttribute( const std::string& attributeName, bool value );
	bool getBoolAttribute( const std::string& attributeName ) const;

	void addJoint( JointPtr joint );
	void removeJoint( JointPtr joint );

	vl::Vec3f fullScale() const;

	// for sound, we will need to be able to splat the object down onto a voxel grid,
	//   which we will use to generate a squashing cubes approximation
	virtual void voxelize( BoundedVoxelGrid& grid, 
		const std::vector<ConstMaterialPtr>& availableMaterials,
		const vl::Mat4f& worldToGrid ) const;

#ifdef SOUND
	void setSquashingCubesModes( boost::shared_ptr<SquashingCubesModes> modes );
	boost::shared_ptr<AudioGenerator> getAudioGenerator(size_t rate) const;
#endif

protected:
	bool equalsLocal( const SceneGraphElement& other ) const;

	virtual ODEGeomResult getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const = 0;
#ifdef USE_NEWTON
	virtual NewtonCollision* getNewtonCollision( const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const = 0;
#endif

private:
	bool isStatic_;

	// these will be used to distinguish between collision-only and
	//   render-only geometry:
	bool hasMass_;
	bool collides_;
	bool hasAudio_;

	vl::Vec3f angularVelocity_;
	vl::Vec3f linearVelocity_;

	boost::shared_ptr<Material> material_;

#ifdef SOUND
	boost::shared_ptr<SquashingCubesModes> modes_;
#endif
};

// represents an object that is a fused collection
//   of primitives
class CombinedObject
	: public PhysicsObject
{
public:
	CombinedObject( const std::string& name, boost::shared_ptr<Material> material );
	virtual ~CombinedObject();

	dMass getOdeMassProps() const;
#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxShapeDesc> > getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const;
#endif
#ifdef USE_BULLET
	btCollisionShape* bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const;
#endif

protected:
	CombinedObject* cloneLocal() const;
	BoundingBox3f bounds( const vl::Mat4f& transform ) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
	RibObjectPtr ribObject() const;

	void renderLocalFill(float alpha, const char* statePtr) const;
	void renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const;
	ODEGeomResult getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const;
#ifdef USE_NEWTON
	NewtonCollision* getNewtonCollision( const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const;
#endif

	bool intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const;
	bool intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const;

	void updateTransformLocal();
	void updateChildrenLocal();

	bool equalsLocal( const SceneGraphElement& other ) const;

	vl::Vec3f furthestPoint( const vl::Vec3f& point ) const;

private:
	std::vector< boost::weak_ptr<PhysicsObject> > objects_;
};

// particular object types:
class BoxObject
	: public PhysicsObject
{
public:
	BoxObject( const std::string& name, boost::shared_ptr<Material> material );
	virtual ~BoxObject();

	dMass getOdeMassProps() const;
#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxShapeDesc> > getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const;
#endif
#ifdef USE_BULLET
	btCollisionShape* bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const;
#endif

	void voxelize( BoundedVoxelGrid& grid, 
		const std::vector<ConstMaterialPtr>& availableMaterials,
		const vl::Mat4f& worldToGrid ) const;

protected:
	BoxObject* cloneLocal() const;
	BoundingBox3f bounds( const vl::Mat4f& transform ) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
	RibObjectPtr ribObject() const;
	void toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const;

	void renderLocalFill(float alpha, const char* statePtr) const;
	void renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const;
	ODEGeomResult getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const;
#ifdef USE_NEWTON
	NewtonCollision* getNewtonCollision( const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const;
#endif

	bool intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const;
	bool intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const;
};

class PlaneObject
	: public PhysicsObject
{
public:
	PlaneObject( const std::string& name, boost::shared_ptr<Material> material );
	virtual ~PlaneObject();

	bool allowUniformScale() const		{ return true; }
	bool allowNonUniformScale() const	{ return true; }

#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxShapeDesc> > getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const;
#endif
#ifdef USE_BULLET
	btCollisionShape* bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const;
#endif

protected:
	PlaneObject* cloneLocal() const;
	BoundingBox3f bounds( const vl::Mat4f& transform ) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
	RibObjectPtr ribObject() const;

	void renderLocalFill(float alpha, const char* statePtr) const;
	void renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const;
	bool intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const;
	bool intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const;

	ODEGeomResult getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const;
#ifdef USE_NEWTON
	NewtonCollision* getNewtonCollision( const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const;
#endif
	dMass getOdeMassProps() const;
};

class SphereObject
	: public PhysicsObject
{
public:
	SphereObject( const std::string& name, boost::shared_ptr<Material> material );
	virtual ~SphereObject();

	bool allowUniformScale() const		{ return true; }
	bool allowNonUniformScale() const	{ return false; }

#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxShapeDesc> > getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const;
#endif
#ifdef USE_BULLET
	btCollisionShape* bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const;
#endif

	vl::Vec3f furthestPoint( const vl::Vec3f& point ) const;


protected:
	void toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const;

	SphereObject* cloneLocal() const;
	BoundingBox3f bounds( const vl::Mat4f& transform ) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
	RibObjectPtr ribObject() const;

	void renderLocalFill(float alpha, const char* statePtr) const;
	void renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const;
	bool intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const;
	bool intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const;

	ODEGeomResult getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const;
#ifdef USE_NEWTON
	NewtonCollision* getNewtonCollision( const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const;
#endif
	dMass getOdeMassProps() const;
};

class ConeObject
	: public PhysicsObject
{
public:
	ConeObject( const std::string& name, boost::shared_ptr<Material> material );
	virtual ~ConeObject();

	bool allowUniformScale() const		{ return true; }
	bool allowNonUniformScale() const	{ return false; }

#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxShapeDesc> > getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const;
#endif
#ifdef USE_BULLET
	btCollisionShape* bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const;
#endif

	RibObjectPtr ribObject() const;

	std::vector<std::string> floatAttributes() const;
	void setAttribute( const std::string& attributeName, float value );
	float getAttribute( const std::string& attributeName ) const;

	void setHeight( float height );
	float height() const;

protected:
	bool equalsLocal( const SceneGraphElement& other ) const;

	ConeObject* cloneLocal() const;
	BoundingBox3f bounds( const vl::Mat4f& transform ) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;

	void renderLocalFill(float alpha, const char* statePtr) const;
	void renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const;
	bool intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const;
	bool intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const;

	ODEGeomResult getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const;
	dMass getOdeMassProps() const;

private:
	float height_;
};

class CappedCylinderObject
	: public PhysicsObject
{
public:
	CappedCylinderObject( const std::string& name, boost::shared_ptr<Material> material );
	virtual ~CappedCylinderObject();

	bool allowUniformScale() const		{ return true; }
	bool allowNonUniformScale() const	{ return false; }
	bool allowShear() const				{ return false; }

	std::vector<std::string> floatAttributes() const;
	void setAttribute( const std::string& attributeName, float value );
	float getAttribute( const std::string& attributeName ) const;

	float lengthRatio() const			{ return this->lengthRatio_; }
	void setLengthRatio(float value) 	{ this->lengthRatio_ = value; }

	void toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const;

#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxShapeDesc> > getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const;
#endif
#ifdef USE_BULLET
	btCollisionShape* bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const;
#endif

protected:
	bool equalsLocal( const SceneGraphElement& other ) const;

	CappedCylinderObject* cloneLocal() const;
	BoundingBox3f bounds( const vl::Mat4f& transform ) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
	RibObjectPtr ribObject() const;

	void renderLocalFill(float alpha, const char* statePtr) const;
	void renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const;
	bool intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const;
	bool intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const;

	ODEGeomResult getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const;
#ifdef USE_NEWTON
	NewtonCollision* getNewtonCollision( const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const;
#endif
	dMass getOdeMassProps() const;

private:
	// ratio of length to radius
	float lengthRatio_;
};

class CylinderObject
	: public PhysicsObject
{
public:
	CylinderObject( const std::string& name, boost::shared_ptr<Material> material );
	virtual ~CylinderObject();

	bool allowUniformScale() const		{ return true; }
	bool allowNonUniformScale() const	{ return false; }
	bool allowShear() const				{ return false; }

	std::vector<std::string> floatAttributes() const;
	void setAttribute( const std::string& attributeName, float value );
	float getAttribute( const std::string& attributeName ) const;

	float lengthRatio() const			{ return this->lengthRatio_; }
	void setLengthRatio(float value) 	{ this->lengthRatio_ = value; }

	void toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const;

#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxShapeDesc> > getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const;
#endif
#ifdef USE_BULLET
	btCollisionShape* bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const;
#endif

protected:
	bool equalsLocal( const SceneGraphElement& other ) const;

	CylinderObject* cloneLocal() const;
	BoundingBox3f bounds( const vl::Mat4f& transform ) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
	RibObjectPtr ribObject() const;

	void renderLocalFill(float alpha, const char* statePtr) const;
	void renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const;
	bool intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const;
	bool intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const;

	ODEGeomResult getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const;
#ifdef USE_NEWTON
	NewtonCollision* getNewtonCollision( const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const;
#endif
	dMass getOdeMassProps() const;

private:
	// ratio of length to radius
	float lengthRatio_;
};


struct ConvexHull
{
	typedef TinyVec<int, 3> Triangle;

	std::vector<vl::Vec3f> vertices;
	std::vector<vl::Vec3f> normals;
	std::vector<Triangle> triangles;
};

typedef boost::shared_ptr<const ConvexHull> ConvexHullPtr;
typedef std::deque<ConvexHullPtr> ConvexHullList;

class MeshObject
	: public PhysicsObject
{
public:
	MeshObject( const std::string& name,
		TriangleMeshPtr mesh,
		boost::shared_ptr<Material> material );
	MeshObject( const std::string& name,
		boost::shared_ptr<Material> material );

	void setMesh( TriangleMeshPtr mesh );
	TriangleMeshPtr getMesh() const;

	virtual ~MeshObject();

	dMass getOdeMassProps() const;

	bool allowNonUniformScale() const	{ return true; }
#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxShapeDesc> > getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const;
#endif
#ifdef USE_BULLET
	btCollisionShape* bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const;
#endif

	void toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const;

	vl::Vec3f furthestPoint( const vl::Vec3f& point ) const;

	void setConvexHulls(const ConvexHullList& hulls);
	ConvexHullList getConvexHulls() const;
	ConvexHullPtr convexHullForPoints( const std::deque<unsigned int>& vertices );
	ConvexHullPtr convexHullForPoints( const std::vector<vl::Vec3d>& vertices );

	void renderPoints( const std::deque<unsigned int>& selected, 
		const vl::Mat4f& transform, 
		const char* state ) const;
	std::pair<bool, unsigned int> selectPoint( const vl::Vec2d& screenPos, 
		const ScreenSpaceConverter& converter, 
		const vl::Mat4f& transform,
		const char* state ) const;

protected:
	MeshObject* cloneLocal() const;
	BoundingBox3f bounds( const vl::Mat4f& transform ) const;

	RibObjectPtr ribObject() const;
	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;

	void renderLocalFill(float alpha, const char* statePtr) const;
	void renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const;
	void renderLocalHulls(const float alpha, const char* statePtr) const;
	bool intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const;
	bool intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const;

	ODEGeomResult getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const;
#ifdef USE_NEWTON
	NewtonCollision* getNewtonCollision( const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const;
#endif
	void updateScaleLocal();

	bool equalsLocal( const SceneGraphElement& other ) const;

private:
	ConvexHullList convexHulls_;

	TriangleMeshPtr mesh_;
#if dTRIMESH_ENABLED
	mutable ODETriMeshDataWrapperPtr triMeshData_;
#endif
};

class sphere_tree;
class SphereWrapperTree;

class TriangleMesh
{
public:
	TriangleMesh( ObjFilePtr objFile );

	const std::vector<vl::Vec3f>& positions() const;
	const std::vector<vl::Vec3f>& normals() const;
	const std::vector<Triangle>& triangles() const;

	const AABBTree<Triangle>& triangleTree() const;
	ObjFilePtr objFile() const;
	BoundingBox3f bounds() const;

#ifdef GUI
	boost::shared_ptr<GLDisplayList> displayList() const;
#endif

	ConvexHullList computeConvexHulls( int depth = 6, int maxVertsPerHull = 32 ) const;

private:
	void init( ObjFilePtr objFile );

	ObjFilePtr objFile_;

	std::vector<vl::Vec3f> positions_;
	std::vector<vl::Vec3f> normals_;
	std::vector<Triangle> triangles_;
	BoundingBox3f bounds_;

//	std::deque<PrimGroup> strips;

	boost::shared_ptr< AABBTree<Triangle> > triangleTree_;

	mutable boost::shared_ptr<GLDisplayList> displayList_;
};

class StaticTriangleMesh
{
public:
	StaticTriangleMesh( 
		const float* positions, 
		size_t numVertices,
		const float* normals, 
		const unsigned int* triangles,
		size_t numTriangles,
		const unsigned int* boundarySegments,
		size_t numSegments );

	const AABBTree<Triangle>& triangleTree() const;
	BoundingBox3f bounds() const;

	const unsigned int* segments() const;
	size_t numSegments() const;

	const vl::Vec3f* positions() const;
	size_t numVertices() const;

#ifdef GUI
	boost::shared_ptr<GLDisplayList> displayList() const;
#endif


private:
	const vl::Vec3f* positions_;
	const vl::Vec3f* normals_;
	size_t numVertices_;

	const Triangle* triangles_;
	size_t numTriangles_;

	const unsigned int* segments_;
	size_t numSegments_;

	BoundingBox3f bounds_;
	boost::shared_ptr< AABBTree<Triangle> > triangleTree_;

	mutable boost::shared_ptr<GLDisplayList> displayList_;
};

typedef std::pair< PhysicsObjectPtr, PhysicsObjectPtr > PhysicsObjectPair;

class Joint
	: public TransformedSceneGraphElement
{
public:
	Joint( const std::string& name, 
		const StaticTriangleMesh* leftMeshObj = 0, 
		const StaticTriangleMesh* rightMeshObj = 0 );
	virtual ~Joint();

	void setObjects( PhysicsObjectPair objects );
	void setObjects( PhysicsObjectPtr object );
    void setObjects( const std::vector<PhysicsObjectPtr>& object );

	bool allowShear() const				{ return true; }

	dJointID odeJoint(dWorldID worldId, dJointGroupID jointGroupId, const ODESimulation& sim) const;
#ifdef USE_NOVODEX
	virtual std::vector< boost::shared_ptr<NxJointDesc> > getNxJointDesc(const NxSimulation& sim) const = 0;
#endif

	virtual void updateControllers( dJointID odeJoint, float time, float timestep ) const;

	float getLowStop( size_t iAxis ) const;
	float getHighStop( size_t iAxis ) const;
	void setLowStop( size_t iAxis, float value ) const;
	void setHighStop( size_t iAxis, float value ) const;

	RibObjectPtr ribObject() const;

	std::vector<ConstPhysicsObjectPtr> objects() const;
	boost::shared_ptr<MechModelDeclaredObject> toMechModelLocal(MechModelType type) const;

	// want to be able to extract only the local attributes:
	std::vector<std::string> jointAttributes() const;

	std::vector<std::string> floatAttributes() const;
	void setAttribute( const std::string& attributeName, float value );
	float getAttribute( const std::string& attributeName ) const;

	vl::Mat4f localTransform( const char* state ) const;

	virtual size_t numObjects() const = 0;
	virtual size_t numControllableAxes() const = 0;

	float disableAfterTime() const;
	void setDisableAfterTime( float value );

	void setCFM( float value );
	void setJointFriction( float value );

protected:
    void clearObjects();

	bool equalsLocal( const SceneGraphElement& other ) const;

	typedef void(*ODEParamFunction)(dJointID, int, dReal);
	typedef dReal(*ODEGetAngleFunction)(dJointID);
	typedef dReal(*ODEGetAngleRateFunction)(dJointID);

	// these functions are essentially the same for every joint
	// so it is much easier to just be able to get the function
	// so we can set most parameters in Joint rather than in the
	// subclasses
	virtual ODEParamFunction odeParamFunction() const = 0;
	virtual ODEGetAngleFunction odeGetAngleFunction(size_t iAxis) const = 0;
	virtual ODEGetAngleRateFunction odeGetAngleRateFunction(size_t iAxis) const = 0;

	void init();

	// find the closest rotation from positive z to direction that doesn't 
	//   violate the joint axis constraints
	virtual vl::Mat4f jointTransform( const vl::Vec3f& direction ) const = 0;

	virtual dJointID createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const = 0;
	virtual void setJointParameters(dJointID id) const = 0;

	std::pair<vl::Vec3f, vl::Vec3f> objectPositionsLocal(const char* state = 0) const;

	BoundingBox3f bounds( const vl::Mat4f& transform ) const;

	void renderLocalFill(float alpha, const char* statePtr) const;
	void renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const;
	bool intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const;
	bool intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const;

	void updateTransformLocal();

#ifdef USE_NOVODEX
	void initNxJointDesc(NxJointDesc& joint, const NxSimulation& sim) const;
#endif

protected:
	vl::Mat4f localJointTransform( size_t iObject,
		const vl::Mat4f& myTransform, const char* statePtr = 0 ) const;

private:
	bool intersectMesh(const BoundingBox2d& box, const ScreenSpaceConverter& converter, 
		const StaticTriangleMesh* mesh, vl::Mat4f transform ) const;
	bool intersectMesh(const Ray3d& ray, double& t, 
		const StaticTriangleMesh* mesh, vl::Mat4f transform ) const;
	BoundingBox3f boundsMesh(const StaticTriangleMesh* mesh, vl::Mat4f transform) const;

	vl::Vec3f objectPositionInLocalFrame( size_t iObject, 
		const vl::Mat4f& myTransform, const char* statePtr = 0 ) const;

	std::vector<dReal> minAngle_;
	std::vector<dReal> maxAngle_;
	std::vector<float> cfm_;
	std::vector<float> jointFriction_;

	// for controllers:
	std::vector<float> goalAngle_;
	std::vector<float> maxForce_;

	float disableAfterTime_;

	std::vector< boost::weak_ptr<PhysicsObject> > objects_;
	std::vector<const StaticTriangleMesh*> meshes_;
	std::vector<const StaticTriangleMesh*> renderedMeshes_;
	std::vector<vl::Mat4f> objectFrameToLocal_;
};

class BallAndSocketJoint
	: public Joint
{
public:
	BallAndSocketJoint( const std::string& name );

#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxJointDesc> > getNxJointDesc(const NxSimulation& sim) const;
#endif

	size_t numObjects() const             { return 2; }
	size_t numControllableAxes() const    { return 0; }

protected:
	ODEParamFunction odeParamFunction() const                             { assert(false); return 0; }
	ODEGetAngleFunction odeGetAngleFunction(size_t iAxis) const           { assert(false); return 0; }
	ODEGetAngleRateFunction odeGetAngleRateFunction(size_t iAxis) const   { assert(false); return 0; }

	BallAndSocketJoint* cloneLocal() const;

	vl::Mat4f jointTransform( const vl::Vec3f& direction ) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;

	dJointID createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const;
	void setJointParameters(dJointID id) const;
};

// hinge joint axis defaults to vertical
class HingeJoint
	: public Joint
{
public:
	HingeJoint( const std::string& name );

#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxJointDesc> > getNxJointDesc(const NxSimulation& sim) const;
#endif

	size_t numObjects() const             { return 2; }
	size_t numControllableAxes() const    { return 1; }

protected:
	ODEParamFunction odeParamFunction() const                             { return &dJointSetHingeParam; }
	ODEGetAngleFunction odeGetAngleFunction(size_t iAxis) const           { return &dJointGetHingeAngle; }
	ODEGetAngleRateFunction odeGetAngleRateFunction(size_t iAxis) const   { return &dJointGetHingeAngleRate; }

	vl::Mat4f jointTransform( const vl::Vec3f& direction ) const;

	HingeJoint* cloneLocal() const;

	dJointID createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const;
	void setJointParameters(dJointID id) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
};

class UniversalJoint
	: public Joint
{
public:
	UniversalJoint( const std::string& name );

#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxJointDesc> > getNxJointDesc(const NxSimulation& sim) const;
#endif

	size_t numObjects() const             { return 2; }
	size_t numControllableAxes() const    { return 2; }

protected:
	ODEParamFunction odeParamFunction() const                             { return &dJointSetUniversalParam; }
	ODEGetAngleFunction odeGetAngleFunction(size_t iAxis) const
	{
		switch( iAxis ) {
			case 0: return &dJointGetUniversalAngle1;
			case 1: return &dJointGetUniversalAngle2;
			default: assert(false); return 0;
		}
	}

	ODEGetAngleRateFunction odeGetAngleRateFunction(size_t iAxis) const
	{
		switch( iAxis ) {
			case 0: return &dJointGetUniversalAngle1Rate;
			case 1: return &dJointGetUniversalAngle2Rate;
			default: assert(false); return 0;
		}
	}

	vl::Mat4f jointTransform( const vl::Vec3f& direction ) const;
	UniversalJoint* cloneLocal() const;

	dJointID createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const;
	void setJointParameters(dJointID id) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
};

class Hinge2Joint
	: public Joint
{
public:
	Hinge2Joint( const std::string& name );

	size_t numObjects() const             { return 2; }
	size_t numControllableAxes() const    { return 1; }

	ODEParamFunction odeParamFunction() const                             { return &dJointSetHinge2Param; }
	ODEGetAngleFunction odeGetAngleFunction(size_t iAxis) const           { return &dJointGetHinge2Angle1; }
	ODEGetAngleRateFunction odeGetAngleRateFunction(size_t iAxis) const   { return &dJointGetHinge2Angle1Rate; }

#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxJointDesc> > getNxJointDesc(const NxSimulation& sim) const;
#endif

protected:
	vl::Mat4f jointTransform( const vl::Vec3f& direction ) const;
	Hinge2Joint* cloneLocal() const;

	dJointID createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const;
	void setJointParameters(dJointID id) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;

	void renderLocalFill(float alpha, const char* statePtr) const;
	void renderLocalWire(
		const vl::Vec4f& wireframeColor, 
		const char* state) const;

private:
	vl::Mat4f middleTransform(const char* statePtr) const;
};

class SliderJoint
	: public Joint
{
public:
	SliderJoint( const std::string& name );

	size_t numObjects() const             { return 2; }
	size_t numControllableAxes() const    { return 1; }

	ODEParamFunction odeParamFunction() const                             { return &dJointSetSliderParam; }
	ODEGetAngleFunction odeGetAngleFunction(size_t iAxis) const           { return &dJointGetSliderPosition; }
	ODEGetAngleRateFunction odeGetAngleRateFunction(size_t iAxis) const   { return &dJointGetSliderPositionRate; }

#ifdef USE_NOVODEX
	std::vector< boost::shared_ptr<NxJointDesc> > getNxJointDesc(const NxSimulation& sim) const;
#endif

protected:
	vl::Mat4f jointTransform( const vl::Vec3f& direction ) const;
	SliderJoint* cloneLocal() const;

	dJointID createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const;
	void setJointParameters(dJointID id) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
};

// This wraps an ODE plane joint, which keeps the object in the xy plane.  
// Unfortunately, this doesn't allow us to hold objects to other planes, 
// although this might be possible in other simulators.
class PlaneJoint
    : public Joint
{
public:
    PlaneJoint( const std::string& name );

	size_t numObjects() const             { return 1; }
	size_t numControllableAxes() const    { return 0; }

	ODEParamFunction odeParamFunction() const                             { assert(false); return 0; }
	ODEGetAngleFunction odeGetAngleFunction(size_t iAxis) const           { assert(false); return 0; }
	ODEGetAngleRateFunction odeGetAngleRateFunction(size_t iAxis) const   { assert(false); return 0; }

    vl::Mat4f jointTransform( const vl::Vec3f& direction ) const { return vl::vl_1; }

protected:
	PlaneJoint* cloneLocal() const;

	dJointID createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const;
	void setJointParameters(dJointID id) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
};

/*
// We'll use these to keep two object rotations identical w.r.t. 
// the y axis:
class OneAxisRotationalMotor
	: public Joint
{
public:
	OneAxisRotationalMotor( const std::string& name );

protected:
	vl::Mat4f jointTransform( const vl::Vec3f& direction ) const;
	OneAxisRotationalMotor* cloneLocal() const;

	dJointID createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const;
	void setJointParameters(dJointID id) const;

	boost::shared_ptr<MechModelDeclaredObject> mechModelObject() const;
};
*/


struct NameEquals
{
	NameEquals( const std::string& name )
		: name_(name) {}

	bool operator()( ConstSceneGraphElementPtr obj ) const
	{
		return (obj->name() == this->name_);
	}

	bool operator()( const SceneGraphElement* obj ) const
	{
		return (obj->name() == this->name_);
	}

	bool operator()( const SceneGraphElement& obj ) const
	{
		return (obj.name() == this->name_);
	}
private:
	std::string name_;
};

vl::Mat4f rigidInverse( const vl::Mat4f& rigidTransform );
vl::Mat4d rigidInverse( const vl::Mat4d& rigidTransform );
bool withinEpsilon( float value1, float value2 );
bool withinEpsilon( const vl::Vec3f& value1, const vl::Vec3f& value2 );
vl::Vec3d toMayaRotation( const Quaternion<double>& q );
vl::Mat4f rigidTransform(const RigidStaticState& state);
bool intersectSphere( const Ray3d& ray, double& t, const vl::Vec3d& center = vl::Vec3d(vl::vl_0), double radius = 1.0 );

} // namespace planning


#endif

