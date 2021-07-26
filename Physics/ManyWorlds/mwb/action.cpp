#include "stdafx.h"

#ifdef GUI

#include "action.h"
#include "planningUI.h"
#include "GLView.h"

#include "twigg/camera.h"

namespace planning {

Action::~Action()
{
}


CompoundAction::CompoundAction( const std::deque< ActionPtr > actions )
	: actions_(actions)
{
}

CompoundAction::~CompoundAction()
{
}


bool CompoundAction::changesScene() const
{
	for( std::deque< ActionPtr >::const_iterator itr = actions_.begin();
		itr != actions_.end(); ++itr )
	{
		if( (*itr)->changesScene() )
			return true;
	}

	return false;
}

bool CompoundAction::changesHierarchy() const
{
	for( std::deque< ActionPtr >::const_iterator itr = actions_.begin();
		itr != actions_.end(); ++itr )
	{
		if( (*itr)->changesHierarchy() )
			return true;
	}

	return false;
}


void CompoundAction::doIt()
{
	std::for_each( actions_.begin(), actions_.end(), boost::mem_fn( &Action::doIt ) );
}

void CompoundAction::undoIt()
{
	// need to be in reverse order:
	std::for_each( actions_.rbegin(), actions_.rend(), boost::mem_fn( &Action::undoIt ) );
}

SetScaleAction::SetScaleAction( const std::vector<TransformedElementPtr>& objects, 
	const std::vector<vl::Vec3f>& backupScales, const std::vector<vl::Vec3f>& newScales )
	: backupScales_(backupScales), newScales_(newScales), objects_(objects)
{
}

SetScaleAction::SetScaleAction( TransformedElementPtr object, 
	const vl::Vec3f& backupScale, const vl::Vec3f& newScale )
	: backupScales_(1, backupScale), newScales_(1, newScale), objects_(1, object)
{
}

void SetScaleAction::doIt()
{
	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
		objects_[iObject]->setScale( newScales_[iObject] );
}

void SetScaleAction::undoIt()
{
	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
		objects_[iObject]->setScale( backupScales_[iObject] );
}


ScaleAction::ScaleAction( const std::vector<TransformedElementPtr>& objects, 
	const std::vector<vl::Vec3f>& backupScales, const vl::Vec3f& newScale )
	: objects_(objects), backupScales_(backupScales), newScale_(newScale)
{
}

void ScaleAction::doIt()
{
	for( unsigned int iObject = 0; iObject < objects_.size(); ++iObject )
		objects_[iObject]->setScale( newScale_ * backupScales_[iObject] );
}

void ScaleAction::undoIt()
{
	for( unsigned int iObject = 0; iObject < objects_.size(); ++iObject )
		objects_[iObject]->setScale( backupScales_[iObject] );
}

TranslateAction::TranslateAction( const std::vector<TransformedElementPtr>& objects, 
	const std::vector<vl::Vec3f>& backupTranslations, 
	const std::vector<vl::Vec3f>& newTranslations )
	: objects_(objects), backupTranslations_(backupTranslations), newTranslations_(newTranslations)
{
	assert( backupTranslations_.size() == objects_.size() );
	assert( newTranslations_.size() == objects_.size() );
}

void TranslateAction::doIt()
{
	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
	{
		objects_[iObject]->setTranslate( 
			this->backupTranslations_[iObject] + this->newTranslations_[iObject] );
	}
}

void TranslateAction::undoIt()
{
	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
	{
		objects_[iObject]->setTranslate( 
			this->backupTranslations_[iObject] );
	}
}

RotateAction::RotateAction( const std::vector<TransformedElementPtr>& objects, 
	const std::vector<vl::Vec3f>& backupRotations, 
	const std::vector<vl::Vec3f>& newRotations )
	: objects_(objects), backupRotations_(backupRotations), newRotations_(newRotations)
{
}

void RotateAction::doIt()
{
	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
	{
		objects_[iObject]->setRotate( newRotations_[iObject] );
	}
}

void RotateAction::undoIt()
{
	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
	{
		objects_[iObject]->setRotate( 
			this->backupRotations_[iObject] );
	}
}


CameraNameChangeAction::CameraNameChangeAction( 
	const CameraWrapperPtr camera, 
	boost::shared_ptr<Scene> scene, 
	const std::string& oldName, 
	const std::string& newName,
	MainFrame* frame )
	:	camera_(camera),
		scene_(scene),
		oldName_(oldName),
		newName_(newName),
		frame_(frame)
{
}

void CameraNameChangeAction::doIt()
{
	scene_->removeCamera( camera_ );
	camera_->setName( newName_ );
	scene_->addCamera( camera_ );
	frame_->updateCameraMenu();
}

void CameraNameChangeAction::undoIt()
{
	scene_->removeCamera( camera_ );
	camera_->setName( oldName_ );
	scene_->addCamera( camera_ );
	frame_->updateCameraMenu();
}

CameraPropertyChangeAction::CameraPropertyChangeAction( boost::shared_ptr<Camera> cam, float zNear, float zFar )
	: camera_(cam), zNew_(zNear, zFar), zOld_(cam->zMin(), cam->zMax())
{
}

void CameraPropertyChangeAction::doIt()
{
	camera_->setZRange( zNew_.first, zNew_.second );
}

void CameraPropertyChangeAction::undoIt()
{
	camera_->setZRange( zOld_.first, zOld_.second );
}

ProjectiveCameraChangeFOVAction::ProjectiveCameraChangeFOVAction( 
	boost::shared_ptr<ProjectiveFixedYCamera> cam, 
	double fov )
	: camera_(cam), fovOld_(cam->FOV()), fovNew_(fov)
{
}

void ProjectiveCameraChangeFOVAction::doIt()
{
	camera_->setFOV(fovNew_);
}

void ProjectiveCameraChangeFOVAction::undoIt()
{
	camera_->setFOV(fovOld_);
}

NameChangeAction::NameChangeAction( const SceneGraphElementPtr object, 
								   boost::shared_ptr<Scene> scene, 
								   const std::string& oldName, 
								   const std::string& newName )
	: object_(object), scene_(scene), oldName_(oldName), newName_(newName)
{
}

void NameChangeAction::doIt()
{
	scene_->removeObject( object_ );
	object_->setName( newName_ );
	scene_->addObject( object_ );
}

void NameChangeAction::undoIt()
{
	scene_->removeObject( object_ );
	object_->setName( oldName_ );
	scene_->addObject( object_ );
}

MaterialNameChangeAction::MaterialNameChangeAction( const MaterialPtr material, 
								   boost::shared_ptr<Scene> scene, 
								   const std::string& oldName, 
								   const std::string& newName )
	: material_(material), scene_(scene), oldName_(oldName), newName_(newName)
{
}

void MaterialNameChangeAction::doIt()
{
	scene_->removeMaterial( material_ );
	material_->setName( newName_ );
	scene_->addMaterial( material_ );
}

void MaterialNameChangeAction::undoIt()
{
	scene_->removeMaterial( material_ );
	material_->setName( oldName_ );
	scene_->addMaterial( material_ );
}

template <typename T>
std::deque<T> toDeque( const T& t )
{
	return std::deque<T>( 1, t );
}

ChangeSelectionAction::ChangeSelectionAction( MainFrame* view, 
	const std::deque<SceneGraphElementPtr>& oldSelected,
	const std::deque<SceneGraphElementPtr>& newSelected,
	ConstraintPtr oldSelectedConstraint,
	ConstraintPtr newSelectedConstraint,
	const std::pair<SceneGraphElementPtr, std::deque<unsigned int> >& oldSelectedPoints,
	const std::pair<SceneGraphElementPtr, std::deque<unsigned int> >& newSelectedPoints )
	:	view_(view), 
		oldSelected_(oldSelected), newSelected_(newSelected), 
		oldSelectedConstraint_(oldSelectedConstraint), newSelectedConstraint_(newSelectedConstraint),
		oldSelectedPoints_(oldSelectedPoints), newSelectedPoints_(newSelectedPoints)
{
}

ChangeSelectionAction::ChangeSelectionAction( MainFrame* view, 
	const std::deque<SceneGraphElementPtr>& oldSelected,
	const std::deque<SceneGraphElementPtr>& newSelected,
	ConstraintPtr oldSelectedConstraint,
	ConstraintPtr newSelectedConstraint )
	:	view_(view), 
		oldSelected_(oldSelected), newSelected_(newSelected), 
		oldSelectedConstraint_(oldSelectedConstraint), newSelectedConstraint_(newSelectedConstraint)
{
	oldSelectedPoints_ = view_->getSelectedPoints();
	if( !newSelected.empty() )
		newSelectedPoints_ = std::make_pair( newSelected.back(), std::deque<unsigned int>() );
}

void ChangeSelectionAction::doIt()
{
	view_->setSelectedConstraint( newSelectedConstraint_ );
	view_->setSelected( newSelected_ );
	view_->setSelectedPoints( newSelectedPoints_ );
}

void ChangeSelectionAction::undoIt()
{
	view_->setSelectedConstraint( oldSelectedConstraint_ );
	view_->setSelected( oldSelected_ );
	view_->setSelectedPoints( oldSelectedPoints_ );
}

CreateCameraAction::CreateCameraAction( CameraWrapperPtr camera, boost::shared_ptr<Scene> scene, MainFrame* frame )
	:	cameras_( 1, camera ),
		scene_(scene),
		frame_(frame)
{
}

void CreateCameraAction::doIt()
{
	scene_->addCameras( cameras_ );
	frame_->updateCameraMenu();
}

void CreateCameraAction::undoIt()
{
	scene_->removeCameras( cameras_ );
	frame_->updateCameraMenu();
}

CreateObjectAction::CreateObjectAction( SceneGraphElementPtr object, boost::shared_ptr<Scene> scene )
	: objects_(1, object), scene_(scene)
{
}

CreateObjectAction::CreateObjectAction( std::deque<SceneGraphElementPtr> objects, boost::shared_ptr<Scene> scene )
	: objects_(objects), scene_(scene)
{
}

void CreateObjectAction::doIt()
{
	scene_->addObjects( objects_ );
}

void CreateObjectAction::undoIt()
{
	scene_->removeObjects( objects_ );
}

DeleteObjectAction::DeleteObjectAction( std::deque< SceneGraphElementPtr > objects, boost::shared_ptr<Scene> scene )
	: objects_(objects), scene_(scene)
{
}

DeleteObjectAction::DeleteObjectAction( SceneGraphElementPtr object, boost::shared_ptr<Scene> scene )
	: objects_(1, object), scene_(scene)
{
}

void DeleteObjectAction::doIt()
{
	scene_->removeObjects( this->objects_ );
}

void DeleteObjectAction::undoIt()
{
	scene_->addObjects( objects_ );
}

MaterialPropertyAction::MaterialPropertyAction( MaterialPtr material, 
	Material oldMat, Material newMat )
	: material_(material), newProperties_(newMat), oldProperties_(oldMat)
{
	assert( newProperties_.name() == material->name() );
	assert( oldProperties_.name() == material->name() );
}

void MaterialPropertyAction::doIt()
{
	*material_ = newProperties_;
}

void MaterialPropertyAction::undoIt()
{
	*material_ = oldProperties_;
}

CreateMaterialAction::CreateMaterialAction( MaterialPtr material, boost::shared_ptr<Scene> scene )
	: material_(material), scene_(scene)
{
}

void CreateMaterialAction::doIt()
{
	scene_->addMaterial( material_ );
}

void CreateMaterialAction::undoIt()
{
	scene_->removeMaterial( material_ );
}

DeleteMaterialAction::DeleteMaterialAction( MaterialPtr material, boost::shared_ptr<Scene> scene )
	: material_(material), scene_(scene)
{
}

void DeleteMaterialAction::doIt()
{
	scene_->removeMaterial( material_ );
}

void DeleteMaterialAction::undoIt()
{
	scene_->addMaterial( material_ );
}

ChangeMaterialAction::ChangeMaterialAction( const std::deque<PhysicsObjectPtr>& objects, 
	const std::deque<MaterialPtr>& oldMaterials, 
	MaterialPtr newMaterial )
	: objects_(objects), oldMaterials_(oldMaterials), newMaterial_(newMaterial)
{
}

void ChangeMaterialAction::doIt()
{
	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
		objects_[iObject]->setMaterial( newMaterial_ );
}

void ChangeMaterialAction::undoIt()
{
	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
		objects_[iObject]->setMaterial( oldMaterials_[iObject] );
}

SetStaticAction::SetStaticAction( bool isStatic, const std::deque< PhysicsObjectPtr >& objects )
	: objects_(objects), newStatic_(isStatic)
{
	std::transform( objects_.begin(), objects_.end(), std::back_inserter(prevStatic_), 
		boost::mem_fn( &PhysicsObject::isStatic ) );
}

void SetStaticAction::doIt()
{
	std::for_each( objects_.begin(), objects_.end(), 
		boost::bind( &PhysicsObject::setStatic, _1, newStatic_ ) );
}

void SetStaticAction::undoIt()
{
	for( unsigned int iObject = 0; iObject < objects_.size(); ++iObject )
		objects_[iObject]->setStatic( prevStatic_[iObject] );
}

AddConstraintAction::AddConstraintAction( ConstraintPtr constraint, SimulationTreePtr tree )
	: constraint_(constraint), tree_(tree)
{
}

void AddConstraintAction::doIt()
{
	constraintResult_ = tree_->addConstraint( constraint_, constraintResult_ );
}

void AddConstraintAction::undoIt()
{
	tree_->removeConstraint( constraint_ );
}

RemoveConstraintAction::RemoveConstraintAction( ConstraintPtr constraint, SimulationTreePtr tree )
	: constraint_(constraint), tree_(tree)
{
}

void RemoveConstraintAction::doIt()
{
	tree_->removeConstraint( constraint_ );
}

void RemoveConstraintAction::undoIt()
{
	tree_->addConstraint( constraint_ );
}

SetPathSelectionAction::SetPathSelectionAction( PathWithTree newPath, 
	PathWithTree oldPath, MainFrame* frame )
	: frame_(frame), oldPath_(oldPath), newPath_(newPath)
{
}

void SetPathSelectionAction::doIt()
{
	this->frame_->setPath( newPath_.first, newPath_.second );
}

void SetPathSelectionAction::undoIt()
{
	this->frame_->setPath( oldPath_.first, oldPath_.second );
}

AddTreeAction::AddTreeAction( MainFrame* frame, SimulationTreePtr tree )
	: frame_(frame), tree_(tree)
{
}

void AddTreeAction::doIt()
{
	frame_->addTree( tree_ );
}

void AddTreeAction::undoIt()
{
	frame_->removeTree( tree_ );
}

SetVisibleAction::SetVisibleAction( const std::deque<SceneGraphElementPtr>& objects, 
	const std::vector<bool>& oldState, const std::vector<bool>& newState )
	: objects_(objects), oldState_(oldState), newState_(newState)
{
	assert( newState.size() == oldState.size() );
	assert( objects.size() == newState.size() );
}

void SetVisibleAction::doIt()
{
	for( size_t iObject = 0; iObject < this->objects_.size(); ++iObject )
		objects_[iObject]->setVisible( newState_[iObject] );
}

void SetVisibleAction::undoIt()
{
	for( size_t iObject = 0; iObject < this->objects_.size(); ++iObject )
		objects_[iObject]->setVisible( oldState_[iObject] );
}

SetTimeRangeAction::SetTimeRangeAction( SimulationTreePtr tree, 
                                       std::pair<float, float> oldEndTime, 
                                       std::pair<float, float> newTimeRange_ )
	: tree_(tree), oldTimeRange_(oldEndTime), newTimeRange_(newTimeRange_)
{
}

void SetTimeRangeAction::doIt()
{
	tree_->setTimeRange( newTimeRange_ );
}

void SetTimeRangeAction::undoIt()
{
	tree_->setTimeRange( oldTimeRange_ );
}

SetTimeAction::SetTimeAction( MainFrame* frame, float oldTime, float newTime )
	: frame_(frame), oldTime_(oldTime), newTime_(newTime)
{
}

void SetTimeAction::doIt()
{
	frame_->setTime( newTime_, true );
}

void SetTimeAction::undoIt()
{
	frame_->setTime( oldTime_, true );
}

ReparentAction::ReparentAction( SceneGraphElementPtr child, SceneGraphElementPtr newParent )
{
	this->object_ = child;

	this->newParent_ = newParent;
	this->oldParent_ = child->parent().lock();
}

void ReparentAction::doIt()
{
	if( this->oldParent_ )
		this->oldNextChild_ = oldParent_->removeChild( this->object_ );

	if( this->newParent_ )
		this->newParent_->addChild( this->object_, this->newNextChild_ );

	this->object_->setParent( newParent_ );
}

void ReparentAction::undoIt()
{
	if( this->newParent_ )
		this->newNextChild_ = this->newParent_->removeChild( this->object_ );

	if( this->oldParent_ )
		this->oldParent_->addChild( this->object_, this->oldNextChild_ );

	this->object_->setParent( oldParent_ );
}

SetRotateAction::SetRotateAction( const std::vector<TransformedElementPtr>& objects, 
	const std::vector<vl::Vec3f>& backupRotations, const std::vector<vl::Vec3f>& newRotations )
	: backupRotations_(backupRotations), newRotations_(newRotations), objects_(objects)
{
	assert( objects_.size() == backupRotations_.size() );
	assert( objects_.size() == newRotations_.size() );
}

SetRotateAction::SetRotateAction( TransformedElementPtr object, 
	const vl::Vec3f& backupRotation, const vl::Vec3f& newRotation )
	: backupRotations_(1, backupRotation), newRotations_(1, newRotation), objects_(1, object)
{
}

void SetRotateAction::doIt()
{
	for( size_t i = 0; i < objects_.size(); ++i )
		objects_.at(i)->setRotate( newRotations_.at(i) );
}

void SetRotateAction::undoIt()
{
	for( size_t i = 0; i < objects_.size(); ++i )
		objects_.at(i)->setRotate( backupRotations_.at(i) );
}

SetTranslateAction::SetTranslateAction( const std::vector<TransformedElementPtr>& objects, 
	const std::vector<vl::Vec3f>& backupTranslations, const std::vector<vl::Vec3f>& newTranslations )
	: backupTranslations_(backupTranslations), newTranslations_(newTranslations), objects_(objects)
{
	assert( objects_.size() == backupTranslations_.size() );
	assert( objects_.size() == newTranslations_.size() );
}

SetTranslateAction::SetTranslateAction( TransformedElementPtr object, 
	const vl::Vec3f& backupTranslation, const vl::Vec3f& newTranslation )
	: backupTranslations_(1, backupTranslation), newTranslations_(1, newTranslation), objects_(1, object)
{
}

void SetTranslateAction::doIt()
{
	for( size_t i = 0; i < objects_.size(); ++i )
		objects_.at(i)->setTranslate( newTranslations_.at(i) );
}

void SetTranslateAction::undoIt()
{
	for( size_t i = 0; i < objects_.size(); ++i )
		objects_.at(i)->setTranslate( backupTranslations_.at(i) );
}

SetPivotAction::SetPivotAction( const std::vector<TransformedElementPtr>& objects, 
	const std::vector<vl::Vec3f>& backupPivots, const std::vector<vl::Vec3f>& newPivots )
	: backupPivots_(backupPivots), newPivots_(newPivots), objects_(objects)
{
	assert( objects_.size() == backupPivots_.size() );
	assert( objects_.size() == newPivots_.size() );
}

SetPivotAction::SetPivotAction( TransformedElementPtr object, 
	const vl::Vec3f& backupPivot, const vl::Vec3f& newPivot )
	: backupPivots_(1, backupPivot), newPivots_(1, newPivot), objects_(1, object)
{
}

void SetPivotAction::doIt()
{
	for( size_t i = 0; i < objects_.size(); ++i )
		objects_.at(i)->setPivot( newPivots_.at(i) );
}

void SetPivotAction::undoIt()
{
	for( size_t i = 0; i < objects_.size(); ++i )
		objects_.at(i)->setPivot( backupPivots_.at(i) );
}

SetFloatAttributeAction::SetFloatAttributeAction( const std::string& attName, 
	float prevValue, 
	float nextValue, 
	SceneGraphElementPtr object )
	: attributeName_(attName), prevValue_(prevValue), nextValue_(nextValue), object_(object)
{
}

void SetFloatAttributeAction::doIt()
{
	object_->setAttribute( attributeName_, nextValue_ );
}

void SetFloatAttributeAction::undoIt()
{
	object_->setAttribute( attributeName_, prevValue_ );
}

SetBoolAttributeAction::SetBoolAttributeAction( const std::string& attName, 
	bool prevValue, 
	bool nextValue, 
	SceneGraphElementPtr object )
	: attributeName_(attName), prevValue_(prevValue), nextValue_(nextValue), object_(object)
{
}

void SetBoolAttributeAction::doIt()
{
	object_->setBoolAttribute( attributeName_, nextValue_ );
}

void SetBoolAttributeAction::undoIt()
{
	object_->setBoolAttribute( attributeName_, prevValue_ );
}

RandomizePropertyAction::RandomizePropertyAction( SceneGraphElementPtr elt, RandomizedProperty prop )
	: object_(elt), newProperty_(prop)
{
	if( object_->isRandomized( newProperty_.name ) )
		this->origProperty_ = object_->getRandomized( newProperty_.name );
}

void RandomizePropertyAction::doIt()
{
	object_->setRandomized( newProperty_ );
}

void RandomizePropertyAction::undoIt()
{
	if( origProperty_.name.empty() )
		object_->clearRandomized( newProperty_.name );
	else
		object_->setRandomized( this->origProperty_ );
}

UnrandomizePropertyAction::UnrandomizePropertyAction( SceneGraphElementPtr elt, const std::string& name )
	: object_(elt)
{
	assert( object_->isRandomized( name ) );
	this->origProperty_ = object_->getRandomized( name );
}


void UnrandomizePropertyAction::doIt()
{
	object_->clearRandomized( origProperty_.name );
}

void UnrandomizePropertyAction::undoIt()
{
	object_->setRandomized( origProperty_ );
}

SetHullsAction::SetHullsAction( boost::shared_ptr<MeshObject> meshObject, 
	const ConvexHullList& origHulls, 
	const ConvexHullList& newHulls )
	: meshObject_(meshObject), origHulls_(origHulls), newHulls_(newHulls)
{
}

SetHullsAction::~SetHullsAction()
{
}

void SetHullsAction::doIt()
{
	meshObject_->setConvexHulls( newHulls_ );
}

void SetHullsAction::undoIt()
{
	meshObject_->setConvexHulls( origHulls_ );
}

} // namespace planning

#endif // GUI

