#ifndef __ACTION_H__
#define __ACTION_H__

#ifdef GUI

#include "scene.h"
#include "physicsFwd.h"
#include "physicsObject.h"
#include "simulationTree.h"

class Camera;
class ProjectiveFixedYCamera;

namespace planning {

class MainFrame;
class GLView;
class SimulationTree;
class SimulationNode;

class CompoundAction
	: public Action
{
public:
	CompoundAction( const std::deque< ActionPtr > actions );
	~CompoundAction();

	void doIt();
	void undoIt();

	bool changesScene() const;
	bool changesHierarchy() const;

private:
	std::deque< ActionPtr > actions_;
};

class NameChangeAction
	: public Action
{
public:
	NameChangeAction( const SceneGraphElementPtr object, boost::shared_ptr<Scene> scene, const std::string& oldName, const std::string& newName );
	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return true; }

private:
	SceneGraphElementPtr object_;
	boost::shared_ptr<Scene> scene_;
	const std::string oldName_;
	const std::string newName_;
};

class CameraPropertyChangeAction
	: public Action
{
public:
	CameraPropertyChangeAction( boost::shared_ptr<Camera> cam, float zNear, float zFar );
	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	boost::shared_ptr<Camera> camera_;
	std::pair<float, float> zOld_;
	std::pair<float, float> zNew_;
};

class ProjectiveCameraChangeFOVAction
	: public Action
{
public:
	ProjectiveCameraChangeFOVAction( boost::shared_ptr<ProjectiveFixedYCamera> cam, double fov );
	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	double fovOld_;
	double fovNew_;
	boost::shared_ptr<ProjectiveFixedYCamera> camera_;
};

class CameraNameChangeAction
	: public Action
{
public:
	CameraNameChangeAction( 
		const CameraWrapperPtr camera, 
		boost::shared_ptr<Scene> scene, 
		const std::string& oldName, 
		const std::string& newName,
		MainFrame* frame );
	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return true; }

private:
	CameraWrapperPtr camera_;
	boost::shared_ptr<Scene> scene_;
	const std::string oldName_;
	const std::string newName_;
	MainFrame* frame_;
};

class SetFloatAttributeAction
	: public Action
{
public:
	SetFloatAttributeAction( const std::string& attName, 
		float prevValue, 
		float nextValue, 
		SceneGraphElementPtr object );
	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::string attributeName_;
	float prevValue_;
	float nextValue_;
	SceneGraphElementPtr object_;
};


class SetBoolAttributeAction
	: public Action
{
public:
	SetBoolAttributeAction( const std::string& attName, 
		bool prevValue, 
		bool nextValue, 
		SceneGraphElementPtr object );
	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::string attributeName_;
	bool prevValue_;
	bool nextValue_;
	SceneGraphElementPtr object_;
};

class MaterialNameChangeAction
	: public Action
{
public:
	MaterialNameChangeAction( const MaterialPtr material, boost::shared_ptr<Scene> scene, const std::string& oldName, const std::string& newName );
	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	MaterialPtr material_;
	boost::shared_ptr<Scene> scene_;
	const std::string oldName_;
	const std::string newName_;
};

class ChangeMaterialAction
	: public Action
{
public:
	ChangeMaterialAction( const std::deque<PhysicsObjectPtr>& objects, 
		const std::deque<MaterialPtr>& oldMaterials, 
		MaterialPtr newMaterial );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::deque<PhysicsObjectPtr> objects_;
	std::deque<MaterialPtr> oldMaterials_;
	MaterialPtr newMaterial_;
};

class MaterialPropertyAction
	: public Action
{
public:
	MaterialPropertyAction( MaterialPtr material, 
		Material oldMat, Material newMat );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	MaterialPtr material_;
	Material newProperties_;
	Material oldProperties_;
};

class SetVisibleAction
	: public Action
{
public:
	SetVisibleAction( const std::deque<SceneGraphElementPtr>& objects, 
		const std::vector<bool>& oldState, const std::vector<bool>& newState );

	void doIt();
	void undoIt();

	bool changesScene() const { return false; }
	bool changesHierarchy() const { return false; }
private:
	std::vector<bool> oldState_;
	std::vector<bool> newState_;
	std::deque<SceneGraphElementPtr> objects_;
};

class ChangeSelectionAction
	: public Action
{
public:
	ChangeSelectionAction( MainFrame* view, 
		const std::deque<SceneGraphElementPtr>& oldSelected,
		const std::deque<SceneGraphElementPtr>& newSelected,
		ConstraintPtr oldSelectedConstraint,
		ConstraintPtr newSelectedConstraint,
		const std::pair<SceneGraphElementPtr, std::deque<unsigned int> >& oldSelectedPoints,
		const std::pair<SceneGraphElementPtr, std::deque<unsigned int> >& newSelectedPoints );
	ChangeSelectionAction( MainFrame* view, 
		const std::deque<SceneGraphElementPtr>& oldSelected,
		const std::deque<SceneGraphElementPtr>& newSelected,
		ConstraintPtr oldSelectedConstraint,
		ConstraintPtr newSelectedConstraint );

	void doIt();
	void undoIt();

	bool changesScene() const { return false; }
	bool changesHierarchy() const { return false; }
private:
	std::deque<SceneGraphElementPtr> oldSelected_;
	std::deque<SceneGraphElementPtr> newSelected_;
	ConstraintPtr oldSelectedConstraint_;
	ConstraintPtr newSelectedConstraint_;
	std::pair<SceneGraphElementPtr, std::deque<unsigned int> > oldSelectedPoints_;
	std::pair<SceneGraphElementPtr, std::deque<unsigned int> > newSelectedPoints_;
	MainFrame* view_;
};


class CreateObjectAction
	: public Action
{
public:
	CreateObjectAction( SceneGraphElementPtr object, boost::shared_ptr<Scene> scene );
	CreateObjectAction( std::deque<SceneGraphElementPtr> objects, boost::shared_ptr<Scene> scene );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return true; }

private:
	std::deque<SceneGraphElementPtr> objects_;
	boost::shared_ptr<Scene> scene_;
};

class DeleteObjectAction
	: public Action
{
public:
	DeleteObjectAction( std::deque< SceneGraphElementPtr > objects, boost::shared_ptr<Scene> scene );
	DeleteObjectAction( SceneGraphElementPtr object, boost::shared_ptr<Scene> scene );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return true; }

private:
	std::deque< SceneGraphElementPtr > objects_;
	boost::shared_ptr<Scene> scene_;
};

// @todo: make this change the used camera in any views as well
class CreateCameraAction
	: public Action
{
public:
	CreateCameraAction( CameraWrapperPtr camera, boost::shared_ptr<Scene> scene, MainFrame* frame );
	
	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return true; }

private:
	std::deque<CameraWrapperPtr> cameras_;
	boost::shared_ptr<Scene> scene_;
	MainFrame* frame_;
};

class ScaleAction
	: public Action
{
public:
	ScaleAction( const std::vector<TransformedElementPtr>& objects, 
		const std::vector<vl::Vec3f>& backupScales, const vl::Vec3f& newScale );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::vector<TransformedElementPtr> objects_;
	std::vector< vl::Vec3f > backupScales_;
	vl::Vec3f newScale_;
};

class SetScaleAction
	: public Action
{
public:
	SetScaleAction( const std::vector<TransformedElementPtr>& objects, 
		const std::vector<vl::Vec3f>& backupScales, const std::vector<vl::Vec3f>& newScales );
	SetScaleAction( TransformedElementPtr object, 
		const vl::Vec3f& backupScale, const vl::Vec3f& newScale );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::vector<vl::Vec3f> backupScales_;
	std::vector<vl::Vec3f> newScales_;
	std::vector<TransformedElementPtr> objects_;
};

class SetRotateAction
	: public Action
{
public:
	SetRotateAction( const std::vector<TransformedElementPtr>& objects, 
		const std::vector<vl::Vec3f>& backupRotations, const std::vector<vl::Vec3f>& newRotations );
	SetRotateAction( TransformedElementPtr object, 
		const vl::Vec3f& backupRotation, const vl::Vec3f& newRotation );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::vector<vl::Vec3f> backupRotations_;
	std::vector<vl::Vec3f> newRotations_;
	std::vector<TransformedElementPtr> objects_;
};

class SetTranslateAction
	: public Action
{
public:
	SetTranslateAction( const std::vector<TransformedElementPtr>& objects, 
		const std::vector<vl::Vec3f>& backupTranslate, const std::vector<vl::Vec3f>& newTranslate );
	SetTranslateAction( TransformedElementPtr object, 
		const vl::Vec3f& backupTranslate, const vl::Vec3f& newTranslate );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::vector<vl::Vec3f> backupTranslations_;
	std::vector<vl::Vec3f> newTranslations_;
	std::vector<TransformedElementPtr> objects_;
};

class SetPivotAction
	: public Action
{
public:
	SetPivotAction( const std::vector<TransformedElementPtr>& objects, 
		const std::vector<vl::Vec3f>& backupPivot, const std::vector<vl::Vec3f>& newPivot );
	SetPivotAction( TransformedElementPtr object, 
		const vl::Vec3f& backupPivot, const vl::Vec3f& newPivot );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::vector<vl::Vec3f> backupPivots_;
	std::vector<vl::Vec3f> newPivots_;
	std::vector<TransformedElementPtr> objects_;
};

class TranslateAction
	: public Action
{
public:
	TranslateAction( const std::vector<TransformedElementPtr>& objects, 
		const std::vector<vl::Vec3f>& backupTranslations, 
		const std::vector<vl::Vec3f>& newTranslations );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::vector<TransformedElementPtr> objects_;
	std::vector<vl::Vec3f> backupTranslations_;
	std::vector<vl::Vec3f> newTranslations_;
};

class RotateAction
	: public Action
{
public:
	RotateAction( const std::vector<TransformedElementPtr>& objects, 
		const std::vector<vl::Vec3f>& backupRotations, 
		const std::vector<vl::Vec3f>& newRotations );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::vector<TransformedElementPtr> objects_;
	std::vector<vl::Vec3f> backupRotations_;
	std::vector<vl::Vec3f> newRotations_;
};

class CreateMaterialAction
	: public Action
{
public:
	CreateMaterialAction( MaterialPtr material, boost::shared_ptr<Scene> scene );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	MaterialPtr material_;
	boost::shared_ptr<Scene> scene_;
};

class DeleteMaterialAction
	: public Action
{
public:
	DeleteMaterialAction( MaterialPtr material, boost::shared_ptr<Scene> scene );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	MaterialPtr material_;
	boost::shared_ptr<Scene> scene_;
};

class SetStaticAction
	: public Action
{
public:
	SetStaticAction( bool isStatic, const std::deque< PhysicsObjectPtr >& objects );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	std::deque< PhysicsObjectPtr > objects_;
	std::deque< bool > prevStatic_;
	bool newStatic_;
};

class AddConstraintAction
	: public Action
{

public:
	AddConstraintAction( ConstraintPtr constraint, SimulationTreePtr tree );

	void doIt();
	void undoIt();

	bool changesScene() const { return false; }
	bool changesHierarchy() const { return false; }
private:
	SimulationTreePtr tree_;
	ConstraintPtr constraint_;

	boost::dynamic_bitset<> constraintResult_;
};

class RemoveConstraintAction
	: public Action
{

public:
	RemoveConstraintAction( ConstraintPtr constraint, SimulationTreePtr tree );

	void doIt();
	void undoIt();

	bool changesScene() const { return false; }
	bool changesHierarchy() const { return false; }
private:
	SimulationTreePtr tree_;
	ConstraintPtr constraint_;
};

class SetPathSelectionAction
	: public Action
{
public:
	typedef std::pair<SimulationTreePtr, const Path*> PathWithTree;
	SetPathSelectionAction( PathWithTree newNode, PathWithTree oldNode, MainFrame* frame );

	void doIt();
	void undoIt();

	bool changesScene() const { return false; }
	bool changesHierarchy() const { return false; }
private:
	MainFrame* frame_;
	PathWithTree newPath_;
	PathWithTree oldPath_;
};

class AddTreeAction
	: public Action
{
public:
	AddTreeAction( MainFrame* frame, SimulationTreePtr tree );
	
	void doIt();
	void undoIt();
	bool changesScene() const { return false; }
	bool changesHierarchy() const { return false; }

private:
	MainFrame* frame_;
	SimulationTreePtr tree_;
};

class SetTimeRangeAction
	: public Action
{
public:
    SetTimeRangeAction( SimulationTreePtr tree, 
        std::pair<float, float> oldEndTime, 
        std::pair<float, float> newEndTime );

	void doIt();
	void undoIt();

	bool changesScene() const { return false; }
	bool changesHierarchy() const { return false; }
private:
	SimulationTreePtr tree_;
	std::pair<float, float> oldTimeRange_;
	std::pair<float, float> newTimeRange_;
};

class SetTimeAction
	: public Action
{
public:
	SetTimeAction( MainFrame* frame, float oldTime, float newTime );

	void doIt();
	void undoIt();

	bool changesScene() const { return false; }
	bool changesHierarchy() const { return false; }
private:
	MainFrame* frame_;
	float oldTime_;
	float newTime_;
};

class ReparentAction
	: public Action
{
public:
	ReparentAction( SceneGraphElementPtr child, SceneGraphElementPtr newParent );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return true; }

private:
	SceneGraphElementPtr object_;

	SceneGraphElementPtr oldParent_;
	SceneGraphElementPtr newParent_;

	SceneGraphElementPtr oldNextChild_;
	SceneGraphElementPtr newNextChild_;
};

class RandomizePropertyAction
	: public Action
{
public:
	RandomizePropertyAction( SceneGraphElementPtr elt, RandomizedProperty prop );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	SceneGraphElementPtr object_;
	RandomizedProperty origProperty_;
	RandomizedProperty newProperty_;
};

class UnrandomizePropertyAction
	: public Action
{
public:
	UnrandomizePropertyAction( SceneGraphElementPtr elt, const std::string& name );

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	SceneGraphElementPtr object_;
	RandomizedProperty origProperty_;
};

class SetHullsAction
	: public Action
{
public:
	SetHullsAction( boost::shared_ptr<MeshObject> meshObject, 
		const ConvexHullList& origHulls, 
		const ConvexHullList& newHulls );
	~SetHullsAction();

	void doIt();
	void undoIt();
	bool changesScene() const { return true; }
	bool changesHierarchy() const { return false; }

private:
	boost::shared_ptr<MeshObject> meshObject_;
	ConvexHullList origHulls_;
	ConvexHullList newHulls_;
};

} // namespace planning

#endif // GUI

#endif
