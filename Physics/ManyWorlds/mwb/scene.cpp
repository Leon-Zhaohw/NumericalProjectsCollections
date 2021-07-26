#include "stdafx.h"

#include "scene.h"
#include "constraints.h"
#include "mechModel.h"
#include "physicsObject.h"
#include "sphere_tree.h"
#include "ribFile.h"
#include "compression.h"
#include "sound.h"

#ifdef USE_BULLET
#include "bulletWrappers.h"
#endif
#include "novodexWrappers.h"

#include "twigg/ioutil.h"
#include "twigg/camera.h"
#include "twigg/objfile.h"
#include "twigg/linalg.h"
#include "twigg/EulerAngles.h"

#include <boost/functional/hash/hash.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <mkl_lapack.h>

#include <fstream>

namespace planning {


template <>
void setAttr<vl::Vec3f>( std::ostream& ofs, const char* attributeName, const vl::Vec3f& value )
{
	ofs.precision( 10 );
	ofs << "\tsetAttr \"" << attributeName << "\" -type \"double3\" " 
		<< value[0] << " " << value[1] << " " << value[2] << ";\n";
}

template <>
void setAttr<vl::Vec3d>( std::ostream& ofs, const char* attributeName, const vl::Vec3d& value )
{
	ofs.precision( 10 );
	ofs << "\tsetAttr \"" << attributeName << "\" -type \"double3\" " 
		<< value[0] << " " << value[1] << " " << value[2] << ";\n";
}

template <>
void setAttr<float>( std::ostream& ofs, const char* attributeName, const float& value )
{
	ofs.precision( 10 );
	ofs << "\tsetAttr \"" << attributeName << "\" " 
		<< value << ";\n";
}

template <>
void setAttr<int>( std::ostream& ofs, const char* attributeName, const int& value )
{
	ofs.precision( 10 );
	ofs << "\tsetAttr \"" << attributeName << "\" " 
		<< value << ";\n";
}

template <>
void setAttr< std::vector<vl::Vec3f> >( std::ostream& ofs, const char* attributeName, const std::vector<vl::Vec3f>& values )
{
	ofs.precision( 10 );
	ofs << "\tsetAttr -s " << values.size() << " \"" << attributeName << "\";\n";
	const size_t maxInSet = 165;
	for( size_t start = 0; start < values.size(); start += maxInSet )
	{
		const size_t end = std::min( values.size(), start+maxInSet );
		ofs << "\tsetAttr \"" << attributeName << "[" << start << ":" << (end-1) << "]\" -type \"float3\" ";
		for( size_t i = start; i < end; ++i )
		{
			if( (i % 3) == 0 && i != start )
				ofs << "\n\t\t";
			for( vl::Int j = 0; j < 3; ++j )
			{
				if( fabs( (values[i])[j] - 1e20) < 1.0f )
					ofs << "1e20 ";
				else
					ofs << (values[i])[j] << " ";
			}
		}
		ofs << ";\n";
	}
}

template <>
void setAttr< std::vector< TinyVec<size_t, 3> > >( std::ostream& ofs, const char* attributeName, const std::vector< TinyVec<size_t, 3> >& values )
{
	ofs.precision( 10 );
	ofs << "\tsetAttr -s " << values.size() << " \"" << attributeName << "\";\n";
	const size_t maxInSet = 165;
	for( size_t start = 0; start < values.size(); start += maxInSet )
	{
		const size_t end = std::min( values.size(), start+maxInSet );
		ofs << "\tsetAttr \"" << attributeName << "[" << start << ":" << (end-1) << "]\" ";
		for( size_t i = start; i < end; ++i )
		{
			if( (i % 3) == 0 && i != start )
				ofs << "\n\t\t";
			for( vl::Int j = 0; j < 3; ++j )
				ofs << (values[i])[j] << " ";
		}
		ofs << ";\n";
	}
}

Scene::Scene()
	: defaultMaterial_( new Material("default") )
{
	{
		boost::shared_ptr<Camera> perspectiveCamera( new ProjectiveFixedYCamera );
		CameraWrapperPtr wrapper( new CameraWrapper("perspectiveCam1", perspectiveCamera) );
		this->sceneCamera_ = wrapper;
		addCamera( wrapper );
	}

	{
		boost::shared_ptr<OrthoCamera> cam( new OrthoCamera );
		cam->setUpVector( vl::Vec3f(1.0, 0.0, 0.0) );
		CameraWrapperPtr wrapper( 
			new CameraWrapper("top", cam) );
		addCamera( wrapper );
	}

	{
		boost::shared_ptr<OrthoCamera> cam( new OrthoCamera );
		CameraWrapperPtr wrapper( 
			new CameraWrapper("front", cam) );
		addCamera( wrapper );
	}

	{
		boost::shared_ptr<OrthoCamera> cam( new OrthoCamera );
		CameraWrapperPtr wrapper( 
			new CameraWrapper("side", cam) );
		addCamera( wrapper );
	}

	root_.reset( new TransformGroup("root") );
	this->addObject( root_ );

	addMaterial( defaultMaterial_ );
}

Scene::~Scene()
{
    // We need to clear out the root first which means that objects
    // get deleted when they get removed from objectNameMap_; otherwise,
    // the scene gets deleted in a giant chain which causes a stack overflow in Windows
    root_.reset();
}

// we will need to keep track of the copies of things so we can
//   run through and fix up all the joints, etc.
typedef STLEXT hash_map<const SceneGraphElement*, SceneGraphElementPtr, 
	hash_ptr<const SceneGraphElement*> > SceneGraphElementMap;

// clones the object plus its whole tree:
SceneGraphElementPtr cloneObject( ConstSceneGraphElementPtr object, Scene* scene, 
								 SceneGraphElementMap& oldToNewMap )
{
	SceneGraphElementPtr result = object->clone();
	std::pair<SceneGraphElementMap::iterator, bool> inserted = 
		oldToNewMap.insert( SceneGraphElementMap::value_type(object.get(), result) );
	// should not have any duplicate entries!
	assert( inserted.second );

	for( SceneGraphElement::const_child_iterator childItr = object->children_begin();
		childItr != object->children_end(); ++childItr )
	{
		SceneGraphElementPtr newChild = cloneObject( *childItr, scene, oldToNewMap );
		newChild->setParent( result );
		result->addChild( newChild );
	}
	scene->addObject( result );
	return result;
}

void fixJoints( SceneGraphElementPtr current, const SceneGraphElementMap& oldToNewMap )
{
	boost::shared_ptr<Joint> joint = 
		boost::dynamic_pointer_cast<Joint>( current );
	if( !joint )
		return;

	std::vector<ConstPhysicsObjectPtr> jointObjects = joint->objects();
	std::vector<PhysicsObjectPtr> newObjects; newObjects.reserve( jointObjects.size() );
	for( std::vector<ConstPhysicsObjectPtr>::const_iterator objectItr = jointObjects.begin();
		objectItr != jointObjects.end(); ++objectItr )
	{
		SceneGraphElementMap::const_iterator mapItr = 
			oldToNewMap.find( objectItr->get() );
		assert( mapItr != oldToNewMap.end() );
		PhysicsObjectPtr ptr = 
			boost::dynamic_pointer_cast<PhysicsObject>( mapItr->second );
		assert( ptr );
		newObjects.push_back( ptr );
	}

	// currently only deal with this
    joint->setObjects( newObjects );
}

Scene::Scene(const Scene& other)
{
	SceneGraphElementMap oldToNewMap;
	this->root_ = boost::dynamic_pointer_cast<TransformGroup>( cloneObject( other.root_, this, oldToNewMap ) );
	assert( this->root_ );

	std::for_each( this->begin_objects(), this->end_objects(), 
		boost::bind( &fixJoints, _1, oldToNewMap ) );

	for( std::deque< MaterialPtr >::const_iterator materialItr = other.materials_.begin();
		materialItr != other.materials_.end(); ++materialItr )
	{
		MaterialPtr newMat( new Material( *(*materialItr) ) );
		this->materials_.push_back( newMat );
		if( *materialItr == other.defaultMaterial_ )
			this->defaultMaterial_ = newMat;
	}

	for( std::deque< CameraWrapperPtr >::const_iterator cameraItr = other.cameras_.begin();
		cameraItr != other.cameras_.end(); ++cameraItr )
	{
		CameraWrapperPtr newCam( (*cameraItr)->clone() );
		this->cameras_.push_back( newCam );
		if( *cameraItr == other.sceneCamera_ )
			this->sceneCamera_ = newCam;
	}
}

SceneGraphElementPtr Scene::root()
{
	return this->root_;
}

ConstSceneGraphElementPtr Scene::root() const
{
	return this->root_;
}

Scene::object_iter Scene::begin_objects() const
{
	return object_iter( this->objectNameMap_.begin() );
}

Scene::object_iter Scene::end_objects() const
{
	return object_iter( this->objectNameMap_.end() );
}


BoundingBox3f Scene::bounds() const
{
	BoundingBox3f result;

	for( object_iter iter = this->begin_objects();
		iter != this->end_objects(); ++iter )
	{
		result.expand( (*iter)->bounds() );
	}

	return result;
}

MechModelObjectPtr Scene::toMechModel() const
{
	boost::shared_ptr<MechModelList> result( new MechModelList );

	// save all the cameras
	std::for_each( cameras_.begin(), cameras_.end(),
		boost::bind( &MechModelList::add, result, 
			boost::bind( &CameraWrapper::toMechModel, _1 ) ) );

	result->add( defaultMaterial_->toMechModel() );
	for( std::deque< MaterialPtr >::const_iterator materialItr = materials_.begin();
		materialItr != materials_.end(); ++materialItr )
	{
		if( *materialItr == defaultMaterial_ )
			continue;

		result->add( (*materialItr)->toMechModel() );
	}

	for( SceneGraphElement::const_child_iterator iter = root_->children_begin();
		iter != root_->children_end(); ++iter )
	{
		boost::shared_ptr<MechModelDeclaredObject> childObject = 
			(*iter)->toMechModel();
		result->add( childObject );
	}

	return result;
}

CameraWrapperPtr Scene::camera()
{
	return this->sceneCamera_;
}

void Scene::setCamera( CameraWrapperPtr cam )
{
	this->sceneCamera_ = cam;
}

template <typename Collection>
typename Collection::value_type namedObject( const std::string& name, const Collection& objects )
{
	typedef typename Collection::value_type ObjectType;
	typedef typename std::deque<ObjectType>::const_iterator Iter;
	typedef std::pair<Iter, Iter> IterPr;

	boost::shared_ptr<Named> named( new Named(name) );
	IterPr pr = std::equal_range( objects.begin(), objects.end(), named, NamedPtrComparator() );
	if( pr.first == pr.second )
		return ObjectType();

	return *pr.first;
}

CameraWrapperPtr Scene::camera( const std::string& name )
{
	return namedObject(name, this->cameras_);
}

void Scene::addCamera( CameraWrapperPtr cameras )
{
	this->addCameras( std::deque<CameraWrapperPtr>(1, cameras) );
}

void Scene::addCameras( std::deque<CameraWrapperPtr> cameras )
{
	std::sort( cameras.begin(), cameras.end(), CameraComparator() );
	{
		std::deque<CameraWrapperPtr> result;
		std::set_union( 
			cameras.begin(), cameras.end(), 
			cameras_.begin(), cameras_.end(),
			std::back_inserter(result),
			CameraComparator() );
		std::swap( cameras_, result );
	}
}

void Scene::removeCamera( CameraWrapperPtr camera )
{
	this->removeCameras( std::deque<CameraWrapperPtr>(1, camera) );
}

void Scene::removeCameras( std::deque<CameraWrapperPtr> cameras )
{
	std::sort( cameras.begin(), cameras.end(), CameraComparator() );
	{
		std::deque<CameraWrapperPtr> result;
		std::set_difference(
			cameras_.begin(), cameras_.end(), 
			cameras.begin(), cameras.end(),
			std::back_inserter( result ),
			CameraComparator() );
		std::swap( cameras_, result );
	}
}

// caller is responsible for ensuring no duplicate objects
void Scene::addObject( SceneGraphElementPtr object )
{
	std::pair<ObjectNameMap::iterator, bool> pr = 
		this->objectNameMap_.insert( 
			ObjectNameMap::value_type(object->name(), object) );
	assert( pr.second );
}

void Scene::addObjects( std::deque<SceneGraphElementPtr> objects )
{
	std::for_each( objects.begin(), objects.end(), 
		boost::bind( &Scene::addObject, this, _1 ) );
}

void Scene::removeObject( SceneGraphElementPtr object )
{
	ObjectNameMap::iterator itr = 
		this->objectNameMap_.find( 
			object->name() );

	if( itr == this->objectNameMap_.end() || 
		itr->second != object )
	{
		std::ostringstream oss;
		oss << "Object '" << object->name() << "' not found in system.";
		throw CreationException( oss.str() );
	}

	this->objectNameMap_.erase( itr );
}


void Scene::removeObjects( std::deque<SceneGraphElementPtr> objects )
{
	std::for_each( objects.begin(), objects.end(), 
		boost::bind( &Scene::removeObject, this, _1 ) );
}

ConstSceneGraphElementPtr Scene::object( const std::string& name ) const
{
	ObjectNameMap::const_iterator itr = 
		this->objectNameMap_.find( name );

	if( itr == this->objectNameMap_.end() )
		return SceneGraphElementPtr();
	else
		return itr->second;
}


SceneGraphElementPtr Scene::object( const std::string& name )
{
	ObjectNameMap::const_iterator itr = 
		this->objectNameMap_.find( name );

	if( itr == this->objectNameMap_.end() )
		return SceneGraphElementPtr();
	else
		return itr->second;
}

MaterialPtr Scene::defaultMaterial()
{
	return this->defaultMaterial_;
}

ConstMaterialPtr Scene::material( const std::string& name ) const
{
	return namedObject( name, materials_ );
}

MaterialPtr Scene::material( const std::string& name )
{
	return namedObject( name, materials_ );
}

void Scene::addMaterial( MaterialPtr material )
{
	typedef std::deque<MaterialPtr>::iterator Iter;
	typedef std::pair< Iter, Iter > IterPr;
	IterPr pr = std::equal_range( materials_.begin(), materials_.end(), material, MaterialComparator() );
	if( pr.first != pr.second )
	{
		std::ostringstream oss;
		oss << "Material '" << material->name() << "' already exists in system!";
		throw CreationException( oss.str() );
	}

	materials_.insert( pr.first, material );
}

void Scene::removeMaterial( MaterialPtr material )
{
	typedef std::deque<MaterialPtr>::iterator Iter;
	typedef std::pair< Iter, Iter > IterPr;
	IterPr pr = std::equal_range( materials_.begin(), materials_.end(), material, MaterialComparator() );
	if( pr.first == pr.second )
	{
		std::ostringstream oss;
		oss << "Material '" << material->name() << "' not found in system.";
		throw CreationException( oss.str() );
	}

	materials_.erase( pr.first, pr.second );
}

size_t Scene::materialCount() const
{
	return materials_.size();
}

MaterialPtr Scene::material( size_t iMat )
{
	return materials_.at( iMat );
}

ConstMaterialPtr Scene::material( size_t iMat ) const
{
	return materials_.at( iMat );
}

TriangleMeshPtr Scene::triMesh( const std::string& filename )
{
	TriMeshMap::iterator iter = meshes_.find( filename );
	if( iter == meshes_.end() )
	{
		ObjFilePtr objFile( new ObjFile(filename) );
/*
		std::string prefix( filename );
		std::string::size_type dotPos = filename.rfind(".");
		if( dotPos != std::string::npos )
			prefix = filename.substr(0, dotPos);
		std::string sphereTreeFilename( prefix + ".sph" );
		boost::shared_ptr<SphereWrapperTree> sphereWrapperTree(
			new SphereWrapperTree(sphereTreeFilename) );
			*/

		iter = meshes_.insert( std::make_pair(filename, 
			TriangleMeshPtr( new TriangleMesh(objFile) ) ) ).first;
	}

	return iter->second;
}

template <typename ObjectType>
std::string uniqueName( const std::string& prefix, const std::deque<ObjectType>& objects )
{
	assert( !prefix.empty() );

	typedef typename std::deque< ObjectType >::const_iterator Iter;
	typedef std::pair<Iter, Iter> IterPr;

	std::string after = prefix;
	++(*after.rbegin());

	boost::shared_ptr<Named> named( new Named(prefix) );
	Iter begin = std::lower_bound( objects.begin(), objects.end(), named, NamedPtrComparator() );

	named.reset( new Named( after ) );
	Iter end = std::upper_bound( objects.begin(), objects.end(), named, NamedPtrComparator() );
	for( int i = 1; ; ++i )
	{
		std::ostringstream oss;
		oss << prefix << i;

		named.reset( new Named(oss.str()) );
		if( !std::binary_search( begin, end, named, NamedPtrComparator() ) )
			return oss.str();

		// in this case, we're probably screwed:
		assert( i < 100000 );
	}
}

std::string Scene::uniqueObjectName( const std::string& prefix ) const
{
	for( int i = 1; ; ++i )
	{
		std::ostringstream oss;
		oss << prefix << i;
		
		ObjectNameMap::const_iterator itr = this->objectNameMap_.find( oss.str() );
		if( itr == this->objectNameMap_.end() )
			return oss.str();

		// in this case, we're probably screwed:
		assert( i < 100000 );
	}
}

std::string Scene::uniqueCameraName( const std::string& prefix ) const
{
	return uniqueName( prefix, cameras_ );
}

void Scene::writeMayaAsciiFile( const std::string& filename,
	const std::vector<PiecewisePath>& paths, 
	size_t frameRate,
	const std::vector<ConstPhysicsObjectPtr>& stateObjects ) const
{
	std::ofstream ofs( filename.c_str() );
	if( !ofs )
		throw IOException( "Unable to open file '" + filename + "' for writing." );

	// boilerplate:
	ofs << "//Maya ASCII 7.0 scene\n";
	ofs << "//Name: " << filename <<";\n";
	ofs << "requires maya \"7.0\";\n";
	ofs << "currentUnit -l centimeter -a degree -t ntsc;\n";
	ofs << "fileInfo \"application\" \"maya\";\n";
	ofs << "fileInfo \"product\" \"Maya Unlimited 7.0\";\n";
	ofs << "fileInfo \"version\" \"7.0\";\n";
	ofs << "fileInfo \"osv\" \"Microsoft Windows XP Service Pack 2 (Build 2600)\";\n";

	ofs << "createNode lightLinker -n \"lightLinker1\";\n";
	ofs << "\tsetAttr -s 3 \".lnk\";\n";
	ofs << "createNode displayLayerManager -n \"layerManager\";\n";
	ofs << "createNode displayLayer -n \"defaultLayer\";\n";
	ofs << "createNode renderLayerManager -n \"renderLayerManager\";\n";
	ofs << "createNode renderLayer -n \"defaultRenderLayer\";\n";
	ofs << "\tsetAttr \".g\" yes;\n";

	ofs << "connectAttr \":defaultLightSet.msg\" \"lightLinker1.lnk[0].llnk\";\n";
	ofs << "connectAttr \":initialShadingGroup.msg\" \"lightLinker1.lnk[0].olnk\";\n";
	ofs << "connectAttr \":defaultLightSet.msg\" \"lightLinker1.lnk[1].llnk\";\n";
	ofs << "connectAttr \":initialParticleSE.msg\" \"lightLinker1.lnk[1].olnk\";\n";
	ofs << "connectAttr \":defaultLightSet.msg\" \"lightLinker1.lnk[2].llnk\";\n";
	ofs << "connectAttr \"layerManager.dli[0]\" \"defaultLayer.id\";\n";
	ofs << "connectAttr \"renderLayerManager.rlmi[0]\" \"defaultRenderLayer.rlid\";\n";
	ofs << "connectAttr \"lightLinker1.msg\" \":lightList1.ln\" -na;\n";


	size_t iMatInfo = 1;
	for( std::deque< MaterialPtr >::const_iterator materialItr = materials_.begin();
		materialItr != materials_.end(); ++materialItr )
	{
		std::string name = (*materialItr)->name();
		if( name == "default" )
			name = "defaultMat";
		std::string sgName = name + "SG";
		std::string matInfoName = "materialInfo" + boost::lexical_cast<std::string>( iMatInfo );
		ofs << "createNode lambert -n \"" << name << "\";\n";
		setAttr( ofs, ".c", (*materialItr)->color() );
		ofs << "createNode shadingEngine -n \"" << sgName << "\";\n";
		ofs << "createNode materialInfo -n \"" << matInfoName << "\";\n";
		ofs << "connectAttr \"" << sgName << ".msg\" \"lightLinker1.lnk[2].olnk\";\n";
		ofs << "connectAttr \"" << name << ".oc\" \"" << sgName << ".ss\";\n";
		ofs << "connectAttr \"" << sgName << ".msg\" \"" << matInfoName << ".sg\";\n";
		ofs << "connectAttr \"" << name << ".msg\" \"" << matInfoName << ".m\";\n";
		ofs << "connectAttr \"" << sgName << ".pa\" \":renderPartition.st\" -na;\n";
		ofs << "connectAttr \"" << name << ".msg\" \":defaultShaderList1.s\" -na;\n";
		++iMatInfo;
	}

	std::deque<ConstPhysicsObjectPtr> dynamicObjects;

	typedef std::pair<ConstSceneGraphElementPtr, std::string> ElementWithParentName;
	std::deque<ElementWithParentName> toHandle(1, ElementWithParentName(this->root_, std::string()));
	while( !toHandle.empty() )
	{
		ConstSceneGraphElementPtr current = toHandle.front().first;
		std::string name = current->name();
		std::string parentName = toHandle.front().second;
		toHandle.pop_front();

		ofs << "createNode transform -n \"" << name << "\"";
		if( !parentName.empty() )
			ofs << " -p \"" << parentName << "\"";
		ofs << ";\n";

		// okay, we could do this with inheritance and virtual functions
		// but we are going to need to perform some major surgery when we
		// add dynamic objects in here
		boost::shared_ptr<const TransformedSceneGraphElement> transformed = 
			boost::dynamic_pointer_cast<const TransformedSceneGraphElement>( current );

		// treat dynamic objects specially
		std::vector<ConstPhysicsObjectPtr>::const_iterator obIter = 
            std::find_if(stateObjects.begin(), stateObjects.end(), NameEquals(current->name()));
		if( obIter != stateObjects.end() && !paths.empty() )
		{
			dynamicObjects.push_back( *obIter );

			// since we'll be using positions in global coordinates, 
			//  don't inherit transform
			ofs << "\tsetAttr \".it\" no;\n";
			setAttr( ofs, ".s", (*obIter)->fullScale() );

			size_t iDynOb = std::distance( stateObjects.begin(), obIter );
			paths.at( iDynOb ).toMayaAsciiFile( ofs, name, frameRate );

            /*
			dMass mass = (*obIter)->getOdeMassProps();
			paths.at( iDynOb ).pathToMayaAsciiFile( ofs, name + "_path1", vl::toVec3d(mass.c) );

			vl::Vec3d farthestPoint = toVec3d( (*obIter)->furthestPoint( vl::toVec3f(mass.c) ) );
			paths.at( iDynOb ).pathToMayaAsciiFile( ofs, name + "_path2", farthestPoint );
            */
 		}
		else if( transformed )
		{
			vl::Vec3d ang = toMayaRotation( eulerToQuat(transformed->rotate()) );

			setAttr( ofs, ".t", transformed->translate() );
			// @todo check if rotation order is correct
			setAttr( ofs, ".r", ang );
			setAttr( ofs, ".s", transformed->scale() );
		}

		std::string matName( ":initialShadingGroup" );
		ConstPhysicsObjectPtr physicsObject = boost::dynamic_pointer_cast<const PhysicsObject>(current);
		if( physicsObject )
		{
			ConstMaterialPtr mat = physicsObject->material();
			matName = mat->name() + "SG";
			if( matName == "defaultSG" )
				matName = "defaultMatSG";
		}

		current->toMayaAscii(ofs, matName);

		for( SceneGraphElement::const_child_iterator itr = current->children_begin();
			itr != current->children_end(); ++itr )
		{
			toHandle.push_front( ElementWithParentName(*itr, name) );
		}
	}

    /*
	if( !paths.empty() )
	{
		// one per second
		const size_t jump = frameRate;
		const size_t startFrame = paths.front().startFrame();
		const size_t endFrame = paths.front().endFrame();
		for( size_t iFrame = startFrame; iFrame < endFrame; iFrame += frameRate )
		{
			std::string framePrefix( "frame" + boost::lexical_cast<std::string>(iFrame) + "_" );
			std::string rootName = framePrefix + "group";
			ofs << "createNode transform -n \"" << rootName << "\";\n";

			for( std::deque<ConstPhysicsObjectPtr>::const_iterator obIter = dynamicObjects.begin();
				obIter != dynamicObjects.end(); ++obIter )
			{
				std::deque<ElementWithParentName> toHandle(1, ElementWithParentName(*obIter, rootName));
				while( !toHandle.empty() )
				{
					ConstSceneGraphElementPtr current = toHandle.front().first;
					std::string name = framePrefix + current->name();
					std::string parentName = toHandle.front().second;
					toHandle.pop_front();

					ofs << "createNode transform -n \"" << name << "\"";
					if( !parentName.empty() )
						ofs << " -p \"" << parentName << "\"";
					ofs << ";\n";

					// okay, we could do this with inheritance and virtual functions
					// but we are going to need to perform some major surgery when we
					// add dynamic objects in here
					boost::shared_ptr<const TransformedSceneGraphElement> transformed = 
						boost::dynamic_pointer_cast<const TransformedSceneGraphElement>( current );

					// treat dynamic objects specially
					std::vector<ConstPhysicsObjectPtr>::const_iterator obIter = 
						std::find_if(stateObjects.begin(), stateObjects.end(), NameEquals(current->name()));
					if( obIter != stateObjects.end() && !paths.empty() )
					{
						// since we'll be using positions in global coordinates, 
						//  don't inherit transform
						ofs << "\tsetAttr \".it\" no;\n";
						setAttr( ofs, ".s", (*obIter)->fullScale() );

						size_t iDynOb = std::distance( stateObjects.begin(), obIter );
						paths.at( iDynOb ).toMayaAsciiFile( ofs, name, frameRate, iFrame, iFrame+jump );
 					}
					else if( transformed )
					{
						vl::Vec3d ang = toMayaRotation( eulerToQuat(transformed->rotate()) );

						setAttr( ofs, ".t", transformed->translate() );
						// @todo check if rotation order is correct
						setAttr( ofs, ".r", ang );
						setAttr( ofs, ".s", transformed->scale() );
					}

					std::string matName( ":initialShadingGroup" );
					ConstPhysicsObjectPtr physicsObject = boost::dynamic_pointer_cast<const PhysicsObject>(current);
					if( physicsObject )
					{
						ConstMaterialPtr mat = physicsObject->material();
						matName = mat->name() + "SG";
						if( matName == "defaultSG" )
							matName = "defaultMatSG";
					}

					current->toMayaAscii(ofs, matName, framePrefix);

					for( SceneGraphElement::const_child_iterator itr = current->children_begin();
						itr != current->children_end(); ++itr )
					{
						toHandle.push_front( ElementWithParentName(*itr, name) );
					}
				}
			}
		}
	}
    */
}

void Scene::writeRibFile( const std::string& filename, 
	const std::deque<Simulation::State>& states, 
	const std::vector<ConstSceneGraphElementPtr>& stateObjects )
{
	std::string prefix( filename );
	std::string::size_type dotPos = filename.rfind(".");
	if( dotPos != std::string::npos )
		prefix = filename.substr(0, dotPos);

	boost::shared_ptr<RibBlock> result( new RibBlock );

	// add camera transforms
	boost::shared_ptr<ProjectiveFixedYCamera> camera = 
		boost::dynamic_pointer_cast<ProjectiveFixedYCamera>( this->camera()->camera() );
	if( camera )
	{
		boost::shared_ptr<RibObjectWithAttributes> object( 
			new RibObjectWithAttributes("Projection") );
		object->addAttribute( std::string("perspective") );
		object->addAttribute( std::string("fov") );
		object->addAttribute( camera->FOV() );
		result->addChild( object );

		boost::shared_ptr<RibObjectWithAttributes> scaleObject(
			new RibObjectWithAttributes( "Scale" ) );
		scaleObject->addAttribute( 1 );
		scaleObject->addAttribute( 1 );
		scaleObject->addAttribute( -1 );
		result->addChild( scaleObject );

		vl::Mat4f xform = camera->modelViewTransform();
		result->addChild( 
			ribObjectWithAttribute( 
				"ConcatTransform", xform ) );

		vl::Mat4f invXform = vl::inv(xform);
		StringForRibValue<vl::Mat4f> converter;
		boost::shared_ptr<RibComment> comment( 
			new RibComment( "inverse: " + converter(invXform) ) );
		result->addChild( comment );
	}

	std::string timesFile;
	if( states.size() > 1 )
	{
		timesFile = prefix + ".times";
		std::ofstream ofs( timesFile.c_str() );
		if( !ofs )
			throw IOException( "Unable to open file '" + timesFile + "' for writing." );
		for( std::deque<Simulation::State>::const_iterator stateIter = states.begin();
			stateIter != states.end(); ++stateIter )
		{
			ofs << stateIter->time() << " ";
		}
		ofs << std::endl;
	}
	boost::filesystem::path timesPath( timesFile );

	typedef STLEXT hash_map<const SceneGraphElement*, size_t, 
		hash_ptr<const SceneGraphElement*> > ObjectToIntMap;
	ObjectToIntMap objectMap;
	for( size_t i = 0; i < stateObjects.size(); ++i )
		objectMap.insert( ObjectToIntMap::value_type( stateObjects.at(i).get(), i ) );

	for( Scene::object_iter iter = this->begin_objects();
		iter != this->end_objects(); ++iter )
	{
		std::vector<const char*> objectStates;

		ObjectToIntMap::const_iterator mapItr = objectMap.find( iter->get() );
		if( mapItr != objectMap.end() )
		{
			size_t iObject = mapItr->second;
			objectStates.reserve( states.size() );
			for( std::deque<Simulation::State>::const_iterator stateIter = states.begin();
				stateIter != states.end(); ++stateIter )
			{
				objectStates.push_back( stateIter->extendedState(iObject) );
			}
		}

		result->addChild( RibObjectPtr( new RibComment( "Object '" + (*iter)->name() + "':" ) ) );
		boost::shared_ptr<RibNamedBlock> attributeBlock =
			(*iter)->toRibObject(prefix, timesPath.leaf(), objectStates);
		if( attributeBlock )
			result->addChild( attributeBlock );
	}

	std::ofstream ofs( filename.c_str() );
	if( !ofs )
		throw IOException( "Unable to open file '" + filename + "' for writing." );
	result->dump( ofs, "" );
}

class RandomPool
{
public:
    typedef boost::shared_ptr<RandomStream> RandomStreamPtr;

    RandomStreamPtr get()
    {
        boost::mutex::scoped_lock lock( this->availableMutex_ );
        if( available_.empty() )
        {
            RandomStreamPtr result( new RandomStream );
            return result;
        }
        else
        {
            RandomStreamPtr result = available_.back();
            available_.pop_back();
            return result;
        }
    }

    void put( RandomStreamPtr stream )
    {
        boost::mutex::scoped_lock lock( this->availableMutex_ );
        available_.push_front( stream );
    }
    
private:
    std::deque< RandomStreamPtr > available_;
    boost::mutex availableMutex_;
};

RandomPool randomStreamPool;


Simulation::Simulation( boost::shared_ptr<const Scene> scene )
	: time_(0.0f)
{

    this->randomStream_ = randomStreamPool.get();

    // initialize all objects connected by joints
    // @todo right now we're only doing objects directly connected by joints
    // but we should do the transitive closure eventually.
    std::vector<const PhysicsObject*> physicsObjects;

	for( Scene::object_iter iter = scene->begin_objects(); 
		iter != scene->end_objects(); ++iter )
	{
		boost::shared_ptr< const PhysicsObject > ob = 
			boost::dynamic_pointer_cast<const PhysicsObject>( *iter );
		if( ob )
            physicsObjects.push_back( ob.get() );
    }

    std::sort( physicsObjects.begin(), physicsObjects.end() );

    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
    Graph G( physicsObjects.size() );

	for( Scene::object_iter iter = scene->begin_objects(); 
		iter != scene->end_objects(); ++iter )
	{
		boost::shared_ptr< const Joint > joint = 
			boost::dynamic_pointer_cast<const Joint>( *iter );
		if( !joint )
            continue;

        std::vector<ConstPhysicsObjectPtr> jointObjects = joint->objects();
        std::vector<int> ids; ids.reserve( jointObjects.size() );
        for( std::vector<ConstPhysicsObjectPtr>::const_iterator itr = jointObjects.begin();
            itr != jointObjects.end(); ++itr )
        {
            typedef std::vector<const PhysicsObject*>::iterator Iter;
            typedef std::pair<Iter, Iter> IterPr;
            IterPr pr = std::equal_range( 
                physicsObjects.begin(), physicsObjects.end(), 
                itr->get() );
            if( pr.first == pr.second )
                continue;

            ids.push_back( std::distance( physicsObjects.begin(), pr.first ) );
        }

        if( ids.empty() )
            continue;

        for( std::vector<int>::const_iterator itr1 = ids.begin(); itr1 != ids.end(); ++itr1 )
            for( std::vector<int>::const_iterator itr2 = itr1+1; itr2 != ids.end(); ++itr2 )
                add_edge( *itr1, *itr2, G );
    }

    std::vector<int> component(num_vertices(G));
    int num = connected_components(G, &component[0]);

    for( size_t i = 0; i < component.size(); ++i )
    {
        connectedComponents_.insert( 
            std::make_pair( physicsObjects.at(i), component.at(i) ) );
    }
}

Simulation::~Simulation()
{
    randomStreamPool.put( this->randomStream_ );
}

void Simulation::setSamplingProperties( size_t iDynamicObject, SamplingProperties props )
{
    this->samplingProperties_.at( iDynamicObject ) = props;
}

SamplingProperties Simulation::getSamplingProperties( size_t iDynamicObject ) const
{
    return this->samplingProperties_.at( iDynamicObject );
}

bool Simulation::asleep() const
{
	for( size_t i = 0; i < this->dynamicObjects_.size(); ++i )
	{
		if( !this->getDisabled(i) )
			return false;
	}

	return true;
}

size_t Simulation::addDynamicObject( ConstPhysicsObjectPtr object )
{
    size_t result = this->dynamicObjects_.size();
    this->dynamicObjects_.push_back( object );
    this->samplingProperties_.push_back( SamplingProperties() );
    return result;
}

int Simulation::dynamicId( ConstPhysicsObjectPtr object ) const
{
	std::vector<ConstPhysicsObjectPtr>::const_iterator bodyItr = 
		std::find( dynamicObjects_.begin(), dynamicObjects_.end(), object );
	if( bodyItr == dynamicObjects_.end() )
	{
		return -1;
	}
	else
	{
        return std::distance( dynamicObjects_.begin(), bodyItr );
	}
}

size_t Simulation::HashPhysicsObjectPair::operator()( 
	const Simulation::PhysicsObjectPair& pr ) const
{
	std::size_t seed = 0;
	boost::hash_combine(seed, pr.first);
	boost::hash_combine(seed, pr.second);
	return seed;
}

bool Simulation::HashPhysicsObjectPair::operator()( 
	const Simulation::PhysicsObjectPair& left, const Simulation::PhysicsObjectPair& right ) const
{
	return left < right;
}

void Simulation::addIgnoreCollisions( const PhysicsObject* left, const PhysicsObject* right )
{
	// symmetric
	PhysicsObjectPair pr(left, right);
	if( left > right )
		std::swap( pr.first, pr.second );

	ignoreCollisionsSet_.insert( pr );
}

void Simulation::addIgnoreCollisions( ConstSceneGraphElementPtr root )
{
	std::deque< ConstPhysicsObjectPtr > allPhysicsObjects;

	std::deque< ConstSceneGraphElementPtr > queue( 1, root );
	while( !queue.empty() )
	{
		ConstSceneGraphElementPtr current = queue.front();
		queue.pop_front();

		boost::shared_ptr<const CombinedObject> fused = 
			boost::dynamic_pointer_cast<const CombinedObject>(current);
		if( !fused )
		{
			std::copy( current->children_begin(), current->children_end(),
				std::front_inserter( queue ) );
		}

		boost::shared_ptr< const PhysicsObject > object = 
			boost::dynamic_pointer_cast<const PhysicsObject>( current );
		if( object )
			allPhysicsObjects.push_back( object );
	}

	for( size_t i = 0; i < allPhysicsObjects.size(); ++i )
		for( size_t j = i+1; j < allPhysicsObjects.size(); ++j )
			this->addIgnoreCollisions( allPhysicsObjects.at(i).get(), allPhysicsObjects.at(j).get() );
}

bool Simulation::ignoreCollisions( const PhysicsObject* left, const PhysicsObject* right ) const
{
	PhysicsObjectPair pr(left, right);
	if( left > right )
		std::swap( pr.first, pr.second );
	return ( ignoreCollisionsSet_.find( pr ) != ignoreCollisionsSet_.end() );
}

bool Simulation::connected( const PhysicsObject* left, const PhysicsObject* right ) const
{
    ConnectedComponentMap::const_iterator leftItr = connectedComponents_.find( left );
    ConnectedComponentMap::const_iterator rightItr = connectedComponents_.find( right );
    return (leftItr != connectedComponents_.end()) &&
           (rightItr != connectedComponents_.end()) &&
           (leftItr->second == rightItr->second);
}


void Simulation::setActiveObjects( PhysicsObjectSet activeObjects )
{
	this->activeObjects_ = activeObjects;
}

bool Simulation::isActive( const PhysicsObject* object ) const
{
	return activeObjects_.find( object ) != activeObjects_.end();
}

std::vector<ConstPhysicsObjectPtr> Simulation::dynamicObjects() const
{
	return dynamicObjects_;
}

std::vector<Simulation::ContactPoint> Simulation::contactPoints() const
{
	return std::vector<Simulation::ContactPoint>();
}

float Simulation::time() const
{
	return time_;
}

void Simulation::setTime(float time)
{
	time_ = time;
}

void Simulation::multiStep( const float endTime, const float stepSize )
{
	float actualStepSize = (endTime < 0) ? -fabs(stepSize) : fabs(stepSize);
	float timeDiff = endTime - time_;
	size_t numSteps = boost::numeric_cast<size_t>(timeDiff / actualStepSize);
	for( size_t i = 0; i < numSteps; ++i )
		step( actualStepSize );
}

bool placeable( dGeomID id )
{
	int geomClass = dGeomGetClass( id );
	if( geomClass == dPlaneClass )
		return false;
	else
		return true;
}


ODESimulation::ODESimulation(boost::shared_ptr<const Scene> scene)
	: Simulation(scene)
{
	// not thread safe
	static ODEInitializer odeInit;

	worldId_ = dWorldCreate();
	//dWorldSetERP( worldId_, 0.2 );
	dWorldSetQuickStepNumIterations( worldId_, 50 );
	dWorldSetGravity( worldId_, 0.0, -9.81, 0.0 );
	dWorldSetAutoDisableFlag( worldId_, 0 );

    dWorldSetContactSurfaceLayer( worldId_, 0.001 );
}

void ODESimulation::initialize(boost::shared_ptr<const Scene> scene)
{
	this->addTree( scene->root() );

	// must do joints after other objects
	for( Scene::object_iter iter = scene->begin_objects(); 
		iter != scene->end_objects(); ++iter )
	{
		boost::shared_ptr< const Joint > joint = 
			boost::dynamic_pointer_cast<const Joint>( *iter );
		if( !joint )
			continue;

		dJointID id = joint->odeJoint( 
			worldId_, 0, *this );

		// WARNING: this introduces a dependence between the simulation and
		//   the Scene that wasn't there before.
		// @todo remove this
		activeJoints_.push_back( JointWithId( joint, id ) );

        boost::shared_ptr<const PlaneJoint> planeJoint = 
            boost::dynamic_pointer_cast<const PlaneJoint>( joint );
        if( planeJoint )
        {
            std::vector<ConstPhysicsObjectPtr> physicsObjects = 
                planeJoint->objects();

            for( std::vector<ConstPhysicsObjectPtr>::const_iterator physObItr = physicsObjects.begin();
                physObItr != physicsObjects.end(); ++physObItr )
            {
                // Need to add the corresponding dynamic objects for all these
                // physics objects to the set of "plane constrained" objects
                dBodyID bodyId = this->odeBodyID( *physObItr );
                if( bodyId )
                    this->planeJointBodies_.push_back( bodyId );
            }
        }
	}

	this->initialState_ = this->getSimulationState();
}

void setPlaneAngularVel( dBodyID bodyID )
{
   const dReal *rot = dBodyGetAngularVel( bodyID );
   const dReal *quat_ptr;
   dReal quat[4], quat_len;
   quat_ptr = dBodyGetQuaternion( bodyID );
   quat[0] = quat_ptr[0];
   quat[1] = 0;
   quat[2] = 0; 
   quat[3] = quat_ptr[3]; 
   quat_len = sqrt( quat[0] * quat[0] + quat[3] * quat[3] );
   quat[0] /= quat_len;
   quat[3] /= quat_len;
   dBodySetQuaternion( bodyID, quat );
   dBodySetAngularVel( bodyID, 0, 0, rot[2] );
}

void ODESimulation::addTree( ConstSceneGraphElementPtr root )
{
	if( !root->selfCollisions() )
		addIgnoreCollisions( root );

	boost::shared_ptr< const PhysicsObject > object = 
		boost::dynamic_pointer_cast<const PhysicsObject>(root);
	if( object )
		this->addObject( object, false );

	boost::shared_ptr<const CombinedObject> fused = 
		boost::dynamic_pointer_cast<const CombinedObject>(root);
	if( !fused )
	{
		for( SceneGraphElement::const_child_iterator childItr = root->children_begin();
			childItr != root->children_end(); ++childItr )
		{
			this->addTree( *childItr );
		}
	}
}

int ODESimulation::addObject( ConstPhysicsObjectPtr object, bool ghost )
{
	float padding = 0.0f;
	if( ghost )
		padding = 1e-1;

	vl::Mat4f transform = object->rigidTransform();
	vl::Vec3f position = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	vl::Mat4f odeRotation = transform;
	for( vl::Int i = 0; i < 3; ++i )
		odeRotation[i][3] = 0.0f;

	dBodyID bodyId = 0;
	vl::Vec3f centerOfMass( vl::vl_0 );

#ifdef SOUND
	boost::shared_ptr<AudioGenerator> audioGenerator =
		object->getAudioGenerator( 44100 );
	if( audioGenerator )
		audioGenerators_.insert( std::make_pair(object.get(), audioGenerator) );
#endif

	int iDynamicId = -1;
	if( !object->isStatic() )
	{
		dMass mass = object->getOdeMassProps();
		if( !dMassCheck( &mass ) )
		{
			// invalid object
			// issue some sort of warning here
		}
		else
		{
			iDynamicId = this->addDynamicObject( object );

			bodyId = dBodyCreate( worldId_ );

			centerOfMass = vl::Vec3f( mass.c[0], mass.c[1], mass.c[2] );
			dMassTranslate( &mass, -mass.c[0], -mass.c[1], -mass.c[2] );
			dBodySetMass( bodyId, &mass );

			// @todo not sure about this one
			position += strip( odeRotation * vl::Vec4f( centerOfMass, 0.0f ) );

			//dBodySetAutoDisableFlag( bodyId, 0 );
			dynamicBodies_.push_back( bodyId );
            ghostBodies_.push_back( ghost );
			offsets_.push_back( centerOfMass );
            bodyToId_.insert( std::make_pair(bodyId, dynamicBodies_.size()-1) );

			std::vector<dReal> rotation( odeRotation.Ref(), odeRotation.Ref() + 12 );
			dBodySetPosition( bodyId, position[0], position[1], position[2] );
			dBodySetRotation( bodyId, &rotation[0] );

			vl::Vec3f linearVel = object->linearVelocity();
			dBodySetLinearVel( bodyId, linearVel[0], linearVel[1], linearVel[2] );
			vl::Vec3f angularVel = object->angularVelocity();
			dBodySetAngularVel( bodyId, angularVel[0], angularVel[1], angularVel[2] );

			if( ghost )
			{
				// don't waste resources simulating "ghost" objects
				dBodyDisable( bodyId );
			}
		}
	}

	this->addCollisionObject( object, iDynamicId, bodyId, centerOfMass, padding, ghost );

	return iDynamicId;
}

void ODESimulation::removeObject( int dynamicId )
{
    assert( dynamicId >= 0 );
    this->removeCollisionObject( dynamicId, this->dynamicBodies_.at(dynamicId) );
}

class GeomExpander
{
public:
	GeomExpander( dGeomID id )
		: id_(id) {}
	virtual ~GeomExpander() {}

	virtual void setExpansion( float expansionDist ) = 0;

	dGeomID geomID() { return this->id_; }

private:
	dGeomID id_;
};

class BoxGeomExpander
	: public GeomExpander
{
public:
	BoxGeomExpander( dGeomID id )
		: GeomExpander( id )
	{
		dGeomBoxGetLengths( id, axes_ );
	}

	void setExpansion( float expansionDist )
	{
		dGeomBoxSetLengths( this->geomID(),
			axes_[0] + 2*expansionDist,
			axes_[1] + 2*expansionDist,
			axes_[2] + 2*expansionDist );
	}

	~BoxGeomExpander()
	{
		dGeomBoxSetLengths( this->geomID(),
			axes_[0],
			axes_[1],
			axes_[2] );
	}

private:
	dReal axes_[4];
};

class SphereGeomExpander
	: public GeomExpander
{
public:
	SphereGeomExpander( dGeomID id )
		: GeomExpander( id )
	{
		this->radius_ = dGeomSphereGetRadius( id );
	}

	void setExpansion( float expansionDist )
	{
		dGeomSphereSetRadius( this->geomID(),
			radius_ + expansionDist );
	}

	~SphereGeomExpander()
	{
		dGeomSphereSetRadius( this->geomID(),
			radius_ );	
	}

private:
	dReal radius_;
};

class PlaneGeomExpander
	: public GeomExpander
{
public:
	PlaneGeomExpander( dGeomID id )
		: GeomExpander( id )
	{
		dGeomPlaneGetParams( id, params_ );
	}

	void setExpansion( float expansionDist )
	{
		dGeomPlaneSetParams( this->geomID(),
			params_[0], params_[1], params_[2], params_[3] + expansionDist );
	}

	~PlaneGeomExpander()
	{
		dGeomPlaneSetParams( this->geomID(),
			params_[0], params_[1], params_[2], params_[3] );
	}

private:
	dReal params_[4];
};

class CappedCylinderGeomExpander
	: public GeomExpander
{
public:
	CappedCylinderGeomExpander( dGeomID id )
		: GeomExpander( id )
	{
		dGeomCCylinderGetParams( this->geomID(), &radius_, &length_ );
	}

	void setExpansion( float expansionDist )
	{
		dGeomCCylinderSetParams( this->geomID(), 
			radius_ + expansionDist, length_ + expansionDist );
	}

	~CappedCylinderGeomExpander()
	{
		dGeomCCylinderSetParams( this->geomID(), radius_, length_ );
	}

private:
	dReal radius_;
	dReal length_;
};

class CylinderGeomExpander
	: public GeomExpander
{
public:
public:
	CylinderGeomExpander( dGeomID id )
		: GeomExpander( id )
	{
		dGeomCylinderGetParams( this->geomID(), &radius_, &length_ );
	}

	void setExpansion( float expansionDist )
	{
		dGeomCylinderSetParams( this->geomID(), 
			radius_ + expansionDist, length_ + expansionDist );
	}

	~CylinderGeomExpander()
	{
		dGeomCylinderSetParams( this->geomID(), radius_, length_ );
	}

private:
	dReal radius_;
	dReal length_;
};

void contactGeomsOnly( void* data, dGeomID o1, dGeomID o2 )
{
	std::vector<dContactGeom>& allGeoms = *reinterpret_cast< std::vector<dContactGeom>* >( data );

	bool isSpace1 = dGeomIsSpace(o1);
	bool isSpace2 = dGeomIsSpace(o2);

	if( isSpace1 || isSpace2 )
	{
		dSpaceCollide2( o1, o2, data, contactGeomsOnly );
	}
	else
	{
		const int maxContacts = 128;

		dContactGeom contacts[ maxContacts ];

		int numc = dCollide( o1, o2, maxContacts, &contacts[0], sizeof(dContactGeom) );
		assert( numc >= 0 );
		if( numc == 0 )
			return;

		std::copy( contacts, contacts + numc, std::back_inserter(allGeoms) );
	}
}

#ifdef SOUND
std::vector<double> ODESimulation::stepAudio(const Simulation::State& prevState, 
							  const Simulation::State& nextState,
							  size_t stateFrameRate,
							  float epsilon,
							  size_t numTimesteps )
{
	this->setSimulationState( prevState );

/*
	// this is unnecessary; there is no reason we can't just create the 
	// scene with all the geometry expanded out a bit

	std::deque<dSpaceID> spaceQueue( 1, this->spaceId_ );
	boost::ptr_deque<GeomExpander> geomExpanders;
	while( !spaceQueue.empty() )
	{
		dSpaceID spaceId = spaceQueue.front();
		spaceQueue.pop_front();
		int numGeoms = dSpaceGetNumGeoms( spaceId );
		for( int iGeom = 0; iGeom < numGeoms; ++iGeom )
		{
			dGeomID geomId = dSpaceGetGeom(spaceId, iGeom);
			if( dGeomGetClass(geomId) == dGeomTransformClass )
				geomId = dGeomTransformGetGeom( geomId );

			if( dGeomIsSpace(geomId) )
			{
				spaceQueue.push_front( reinterpret_cast<dSpaceID>(geomId) );
			}
			else
			{
				switch( dGeomGetClass(geomId) )
				{
				case dBoxClass:
					geomExpanders.push_back( new BoxGeomExpander(geomId) );
					break;
				case dSphereClass:
					geomExpanders.push_back( new SphereGeomExpander(geomId) );
					break;
				case dPlaneClass:
					geomExpanders.push_back( new PlaneGeomExpander(geomId) );
					break;
				case dCCylinderClass:
					geomExpanders.push_back( new CappedCylinderGeomExpander(geomId) );
					break;
				case dCylinderClass:
					geomExpanders.push_back( new CylinderGeomExpander(geomId) );
					break;
				case dTriMeshClass:
					assert( false );
					break;
				case dGeomTransformClass:
				case dRayClass:
					assert( false );
					break;
				case dHashSpaceClass:
				case dSimpleSpaceClass:
					// we should have taken care of spaces with the if statement above
					assert( false );
					break;
				default:
					assert(false);
				}
			}
		}
	}

	std::for_each( geomExpanders.begin(), geomExpanders.end(),
		boost::bind( &GeomExpander::setExpansion, _1, epsilon ) );

	std::vector<dContactGeom> contactGeoms;

	// perform collision detection
	dSpaceCollide( this->spaceId_, &contactGeoms, contactGeomsOnly );

	if( contactGeoms.empty() )
	{
		return this->stepAudio( numTimesteps );
	}

	// So we need the forces to be positive.  The correct way to do this is
	// with some kind of active set algorithm.  For now, we'll fake this using
	// a very simple one which just keeps tossing out negative forces until they're
	// all positive.

	// rhs only needs to get set up once
	double dt = nextState.time() - prevState.time();
	std::vector<double> rhs( 6*this->dynamicBodies_.size() );
	for( size_t iBody = 0; iBody < this->dynamicBodies_.size(); ++iBody )
	{
		vl::Vec3f linVelDiff = 
			toVec3f( nextState.state(iBody).linearVelocity() -
				prevState.state(iBody).linearVelocity() ) 
					* static_cast<float>(stateFrameRate); 
		vl::Vec3f angVelDiff = 
			toVec3f( nextState.state(iBody).angularVelocity() -
				prevState.state(iBody).angularVelocity() )
					* static_cast<float>(stateFrameRate); 

		dMass mass;
		dBodyGetMass( this->dynamicBodies_.at(iBody), &mass );
		for( vl::Int i = 0; i < 3; ++i )
			rhs[6*iBody + i] = (mass.mass/dt) * linVelDiff[i];

		vl::Mat3f bodyInertiaTensor;
		for( vl::Int i = 0; i < 3; ++i )
			for( vl::Int j = 0; j < 3; ++j )
				bodyInertiaTensor[i][j] = mass.I[ i*4 + j ];

		vl::Mat3f rotMat;
		{
			NxMat33 rmTemp;
			rmTemp.fromQuat( nextState.state(iBody).orientation() );
			rmTemp.getRowMajor( rotMat.Ref() );
		}

		vl::Mat3f transformedInertia = 
			rotMat * bodyInertiaTensor * trans(rotMat);
		vl::Vec3f rotImpulse = transformedInertia * angVelDiff;

		for( vl::Int i = 0; i < 3; ++i )
			rhs[6*iBody + 3 + i] = rotImpulse[i] / dt;
	}

	while( !contactGeoms.empty() )
	{
		fortran_matrix A( 6*this->dynamicBodies_.size(), contactGeoms.size() );
		std::fill( A.data(), A.data() + A.nrows()*A.ncols(), 0.0 );

		for( size_t iContact = 0; iContact < contactGeoms.size(); ++iContact )
		{
			dGeomID geoms[2] = { contactGeoms[iContact].g1,
				contactGeoms[iContact].g2 };

			vl::Vec3d normal = toVec3d( contactGeoms[iContact].normal );
			vl::Vec3d pos = toVec3d( contactGeoms[iContact].pos );
			for( size_t i = 0; i < 2; ++i )
			{
				dBodyID body = dGeomGetBody( geoms[i] );
				if( body == 0 )
					continue;

				// todo make this not linear
				std::deque<dBodyID>::iterator itr = 
					std::find( this->dynamicBodies_.begin(), this->dynamicBodies_.end(), body );
				assert( itr != dynamicBodies_.end() );
				size_t iBody = std::distance( dynamicBodies_.begin(), itr );

				vl::Vec3d position = toVec3d( dBodyGetPosition(body) );
				vl::Vec3d r = position - pos;
				vl::Vec3d torqueDir = vl::cross( r, normal );
				for( vl::Int i = 0; i < 3; ++i )
					A( 6*iBody + i, iContact ) = normal[i];
				for( vl::Int i = 0; i < 3; ++i )
					A( 6*iBody + 3 + i, iContact ) = torqueDir[i];

				normal *= -1.0;
			}
		}

#ifdef _DEBUG
		boost::shared_ptr<MATFile> test( new MATFile("//volley/cdtwigg/test.mat", "test") );
		test->add( "rhs", rhs );
		test->add( "A", A );
#endif

		// Okay, matrix has been set up
		// perform least squares solve
		std::vector<double> f = leastSquares_QR( A, &rhs[0] );
		std::vector<dContactGeom> newGeoms;
		newGeoms.reserve( contactGeoms.size() );

#ifdef _DEBUG
		test->add( "f", f );
		test.reset();
#endif

		for( size_t iContact = 0; iContact < contactGeoms.size(); ++iContact )
		{
			if( f[iContact] >= 0.0 )
				newGeoms.push_back( contactGeoms[iContact] );
		}

		if( newGeoms.size() == contactGeoms.size() )
		{
			// we're done!
			// add impulses
			for( size_t iContact = 0; iContact < contactGeoms.size(); ++iContact )
			{
				const dContactGeom& contactGeom = contactGeoms.at(iContact);
				boost::array<dGeomID, 2> geoms = {{ contactGeom.g1, contactGeom.g2 }};

				vl::Vec3f normal = toVec3f( contactGeom.normal );
				vl::Vec3f pos = toVec3f( contactGeom.pos );

				for( size_t i = 0; i < 2; ++i )
				{
					const GeomInfo *info = reinterpret_cast<const GeomInfo*>( 
						dGeomGetData( geoms[i] ) );

					ObjectToAudioHash::const_iterator iter = 
						audioGenerators_.find( info->originalObject );
					if( iter == audioGenerators_.end() )
						continue;

					int iDynId = info->iDynamicId;
					vl::Mat4f transform( vl::vl_1 );
					if( iDynId >= 0 )
						transform = info->originalObject->rigidTransform( prevState.extendedState(iDynId) );

					vl::Mat3f rotate = toMat3f( transform );
					vl::Vec3f transformedPos = vl::xform( rigidInverse( transform ), pos );
					// for rigid transform, inv(trans(mat)) == mat
					vl::Vec3f normalForce = f[iContact] * normal;
					vl::Vec3f transformedForce = rotate * normalForce;

					iter->second->addImpulse( transformedPos, transformedForce, dt );
				}
			}

			return this->stepAudio( numTimesteps );
		}

		newGeoms.swap( contactGeoms );
	}
	*/

	// all impulses removed from set
	return this->stepAudio( numTimesteps );
}
#endif

void ODESimulation::perturbLinearVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation )
{
	dBodyID bodyId = this->dynamicBodies_.at(iDynamicObject);

	dReal linearVel[] = { 0.0f, 0.0f, 0.0f, 0.0f };
	const dReal* internalVel = dBodyGetLinearVel( bodyId );
	std::copy( internalVel, internalVel+4, linearVel );
	for( vl::Int i = 0; i < 3; ++i )
		linearVel[i] += perturbation[i];

	dBodySetLinearVel( bodyId, linearVel[0], linearVel[1], linearVel[2] );
}

void ODESimulation::perturbAngularVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation )
{
	dBodyID bodyId = this->dynamicBodies_.at(iDynamicObject);

	dReal angularVel[] = { 0.0f, 0.0f, 0.0f, 0.0f };
	const dReal* internalVel = dBodyGetAngularVel( bodyId );
	std::copy( internalVel, internalVel+4, angularVel );
	for( vl::Int i = 0; i < 3; ++i )
		angularVel[i] += perturbation[i];

	dBodySetAngularVel( bodyId, angularVel[0], angularVel[1], angularVel[2] );
}

vl::Vec3f ODESimulation::getLinearVelocity( size_t iDynamicObject ) const
{
	dBodyID bodyId = this->dynamicBodies_.at(iDynamicObject);
	const dReal* internalVel = dBodyGetLinearVel( bodyId );
	vl::Vec3f result;
	std::copy( internalVel, internalVel+3, result.Ref() );
	return result;
}

vl::Vec3f ODESimulation::getAngularVelocity( size_t iDynamicObject ) const
{
	dBodyID bodyId = this->dynamicBodies_.at(iDynamicObject);
	const dReal* internalVel = dBodyGetAngularVel( bodyId );
	vl::Vec3f result;
	std::copy( internalVel, internalVel+3, result.Ref() );
	return result;
}

void ODESimulation::dumpState( const char* filename ) const
{
    std::ofstream ofs( filename, std::ios::binary );
    assert( ofs );

    for( size_t iBody = 0; iBody < this->dynamicBodies_.size(); ++iBody )
    {
        dBodyID id = this->dynamicBodies_.at( iBody );
        ofs.write( reinterpret_cast<const char*>( dBodyGetPosition(id) ), 4*sizeof(dReal) );
        ofs.write( reinterpret_cast<const char*>( dBodyGetQuaternion(id) ), 4*sizeof(dReal) );
        ofs.write( reinterpret_cast<const char*>( dBodyGetLinearVel(id) ), 4*sizeof(dReal) );
        ofs.write( reinterpret_cast<const char*>( dBodyGetAngularVel(id) ), 4*sizeof(dReal) );
    }
}

void ODESimulation::loadState( const char* filename )
{
    std::ifstream ifs( filename, std::ios::binary );
    assert( this->dynamicBodies_.size()*16*sizeof(dReal) == filelen(ifs) );

    for( size_t iBody = 0; iBody < this->dynamicBodies_.size(); ++iBody )
    {
        dBodyID id = this->dynamicBodies_.at( iBody );
        dReal pos[4], rot[4], lvel[4], avel[4];
        ifs.read( reinterpret_cast<char*>( pos ), 4*sizeof(dReal) );
        ifs.read( reinterpret_cast<char*>( rot ), 4*sizeof(dReal) );
        ifs.read( reinterpret_cast<char*>( lvel ), 4*sizeof(dReal) );
        ifs.read( reinterpret_cast<char*>( avel ), 4*sizeof(dReal) );

        dBodySetPosition( id, pos[0], pos[1], pos[2] );
        dBodySetQuaternion( id, rot );
        dBodySetLinearVel( id, lvel[0], lvel[1], lvel[2] );
        dBodySetAngularVel( id, avel[0], avel[1], avel[2] );
    }
}

void ODESimulation::updateJoints(float time)
{
	for( std::deque<JointWithId>::iterator itr = activeJoints_.begin();
		itr != activeJoints_.end(); )
	{
		if( itr->first->disableAfterTime() < time )
		{
			dJointDestroy( itr->second );
			inactiveJoints_.push_back( itr->first );
			itr = activeJoints_.erase( itr );
		}
		else
		{
			++itr;
		}
	}

	for( std::deque<ConstJointPtr>::iterator itr = inactiveJoints_.begin();
		itr != inactiveJoints_.end(); )
	{
		if( (*itr)->disableAfterTime() >= time )
		{
			// need to restore initial states of all the involved bodies because
			//  Joint::odeJoint assumes that objects are in initial position.
			std::vector<ConstPhysicsObjectPtr> jointObjects = (*itr)->objects();
			std::vector<size_t> objectIds; objectIds.reserve( jointObjects.size() );
			for( std::vector<ConstPhysicsObjectPtr>::iterator objectItr = jointObjects.begin();
				objectItr != jointObjects.end(); )
			{
                int iDynamicId = this->dynamicId( *objectItr );
                if( iDynamicId < 0 )
                {
					objectItr = jointObjects.erase(objectItr);
				}
				else
				{
					objectIds.push_back( iDynamicId );
					++objectItr;
				}
			}

			std::vector<RigidDynamicState> backupStates;
			std::transform( objectIds.begin(), objectIds.end(), 
				std::back_inserter(backupStates), boost::bind(&ODESimulation::state, this, _1) );

			for( std::vector<size_t>::const_iterator bodyItr = objectIds.begin();
				bodyItr != objectIds.end(); ++bodyItr )
			{
				this->setState( *bodyItr, this->initialState_.state(*bodyItr) );
			}

			// Now create joint
			dJointID id = (*itr)->odeJoint( 
				worldId_, 0, *this );
			activeJoints_.push_back( JointWithId( (*itr), id ) );
			itr = inactiveJoints_.erase( itr );

			// restore old state
			for( std::vector<size_t>::const_iterator bodyItr = objectIds.begin();
				bodyItr != objectIds.end(); ++bodyItr )
			{
				this->setState( *bodyItr, backupStates.at(*bodyItr) );
			}
		}
		else
		{
			++itr;
		}
	}
}

dBodyID ODESimulation::odeBodyID( ConstPhysicsObjectPtr object ) const
{
	// @todo make this faster!
	std::vector<ConstPhysicsObjectPtr> dynamicObjects = this->dynamicObjects();
	for( size_t iBody = 0; iBody < dynamicObjects.size(); ++iBody )
	{
		if( dynamicObjects[iBody] == object )
			return this->dynamicBodies_.at(iBody);
	}

	assert( false );
	throw Exception( "No such dynamic body!" );
}

RigidDynamicState ODESimulation::state( size_t iObject ) const
{
	dBodyID bodyId = dynamicBodies_[iObject];
	const dReal* odeQuat = dBodyGetQuaternion(bodyId);
	Quaternion<float> quat( odeQuat, Quaternion<float>::ORDER_WXYZ );

	vl::Mat3f rotMat = quat.toRotMatf();

	RigidDynamicState s(
		vl::toVec3f(dBodyGetPosition(bodyId)) - rotMat*vl::toVec3f(offsets_[iObject]),
		quat,
		vl::toVec3f(dBodyGetLinearVel(bodyId)),
		vl::toVec3f(dBodyGetAngularVel(bodyId)) );
	return s;
}

void ODESimulation::setState( size_t iObject, const RigidDynamicState& s )
{
	dBodyID id = this->dynamicBodies_[iObject];
	{
		Quaternion<float> quat = s.orientation();
		dQuaternion q = { quat.w(), quat.x(), quat.y(), quat.z() };
		dBodySetQuaternion( id, q );

		vl::Mat3f rotMat = quat.toRotMatf();
		vl::Vec3f pos = s.position() + rotMat*offsets_[iObject];
		dBodySetPosition( id, pos[0], pos[1], pos[2] );
	}
	{
		vl::Vec3f vel = s.linearVelocity();
		dBodySetLinearVel( id, vel[0], vel[1], vel[2] );
	}

	{
		vl::Vec3f vel = s.angularVelocity();
		dBodySetAngularVel( id, vel[0], vel[1], vel[2] );
	}
}

RigidDynamicState ODESimulation::getSimulationState(size_t iBody) const
{
    return this->state( iBody );
}

Simulation::State ODESimulation::getSimulationState() const
{
	std::vector<RigidDynamicState> state;
	state.reserve( dynamicBodies_.size() );
	for( size_t i = 0; i < dynamicBodies_.size(); ++i )
	{
		state.push_back( this->state(i) );
	}

	return Simulation::State(state, time());
}

void ODESimulation::setSimulationState( const Simulation::State& state )
{
	assert( state.stateCount() == this->dynamicBodies_.size() );

	updateJoints( state.time() );

	setTime( state.time() );
	for( size_t i = 0; i < dynamicBodies_.size(); ++i )
		this->setState( i, state.state(i) );
}

void ODESimulation::setSimulationState( size_t iObject, const RigidDynamicState& state )
{
	this->setState( iObject, state );
}

void ODESimulation::setDisabled( size_t iDynamicObject, bool disabled = true )
{
	if( disabled )
		dBodyDisable( this->dynamicBodies_[iDynamicObject] );
	else
		dBodyEnable( this->dynamicBodies_[iDynamicObject] );
}

bool ODESimulation::getDisabled( size_t iDynamicObject ) const
{
	return (dBodyIsEnabled( this->dynamicBodies_[iDynamicObject] ) == 0);
}

ODESimulation::~ODESimulation()
{
	// should destroy all active joints:
	std::for_each( contactJoints_.begin(), contactJoints_.end(),
		&dJointDestroy );

	for( std::deque<JointWithId>::const_iterator itr = activeJoints_.begin();
		itr != activeJoints_.end(); ++itr )
	{
		dJointDestroy( itr->second );
	}

	// this should destroy all bodies too:
	dWorldDestroy( worldId_ );
}

/*
void ODESimulation::addContact(const dContactGeom& contactGeom)
{
	MaterialMap::const_iterator itr1 = materialMap_.find( contactGeom.g1 );
	assert( itr1 != materialMap_.end() );
	const Material& mat1 = itr1->second;

	MaterialMap::const_iterator itr2 = materialMap_.find( contactGeom.g2 );
	assert( itr2 != materialMap_.end() );
	const Material& mat2 = itr2->second;

	dContact contact;
	contact.surface.mode = dContactBounce | dContactApprox1;

	// is max appropriate here?  also should switch between dynamic and static friction
	contact.surface.mu = std::max( mat1.dynamicFriction(), mat2.dynamicFriction() );

	contact.surface.bounce = std::max( mat1.restitution(), mat2.restitution() );
	contact.surface.bounce_vel = 0.01f;

	dJointID jointId = dJointCreateContact( this->worldId_, this->collisionJoints_, &contact );

	dBodyID body1 = 0;
	if( placeable(contactGeom.g1) )
		body1 = dGeomGetBody(contactGeom.g1);

	dBodyID body2 = 0;
	if( placeable(contactGeom.g2) )
		body2 = dGeomGetBody(contactGeom.g2);

	dJointAttach( jointId, body1, body2 );

	dJointAttach( jointId, dGeomGetBody(contactGeom.g1), dGeomGetBody(contactGeom.g2) );
}
*/

void ODESimulation::step(const float step)
{
    assert( step != 0 );

    stepSize_ = std::abs( step );
    backwards_ = (step < 0.0);
	if( step < 0.0 )
	{
        // first, need to back all the positions up.
		for( std::deque<dBodyID>::const_iterator bodyItr = this->dynamicBodies_.begin();
			bodyItr != this->dynamicBodies_.end(); ++bodyItr )
		{
            const dReal* linearVel = dBodyGetLinearVel(*bodyItr);
            const dReal* position = dBodyGetPosition(*bodyItr);
			dBodySetPosition( *bodyItr, 
                position[0] + step*linearVel[0], 
                position[1] + step*linearVel[1], 
                position[2] + step*linearVel[2] );

            const dReal* angularVelocity = dBodyGetAngularVel(*bodyItr);
            vl::Vec3d diffRot;
            for( vl::Int i = 0; i < 3; ++i )
                diffRot[i] = step*angularVelocity[i];

			Quaternion<dReal> quat( 
				dBodyGetQuaternion(*bodyItr), 
				Quaternion<double>::ORDER_WXYZ );
			Quaternion<dReal> diffRotQuat( diffRot );
			Quaternion<dReal> rotation(
				diffRotQuat	* Quaternion<dReal>( quat ) );
			boost::array<dReal, 4> q( rotation.get(Quaternion<dReal>::ORDER_WXYZ) );
			dBodySetQuaternion( *bodyItr, &q[0] );

            double dampingCoeff = 0;
            dBodyAddTorque( *bodyItr,
                            dampingCoeff*angularVelocity[0], 
                            dampingCoeff*angularVelocity[1], 
                            dampingCoeff*angularVelocity[2] );
            /*
#ifdef _DEBUG
            if( !ghostBodies_[std::distance(std::deque<dBodyID>::const_iterator(dynamicBodies_.begin()), bodyItr)] )
            {
                if( vl::len(linearVel) > 1e-3 || vl::len(angularVelocity) > 1e-3 )
                    assert( dBodyIsEnabled(*bodyItr) );
            }
#endif
            */
        }

		for( std::deque<JointWithId>::const_iterator jointItr = activeJoints_.begin();
			jointItr != activeJoints_.end(); ++jointItr )
		{
			(jointItr->first)->updateControllers( jointItr->second, time(), step );
		}

		// clear out all contact joints:
		std::for_each( contactJoints_.begin(), contactJoints_.end(),
			&dJointDestroy );
		contactJoints_.clear();
		contactGeoms_.clear();

		// perform collision detection
		this->handleCollisions(true);

		// now need to mark contact geoms as needing feedback
		contactJointFeedback_.resize( contactJoints_.size() );
		for( size_t iJoint = 0; iJoint < contactJoints_.size(); ++iJoint )
			dJointSetFeedback( contactJoints_[iJoint], &contactJointFeedback_[iJoint] );

        // Need to traverse the collision graph and figure out which collisions aren't
        // allowed to bounce!
        /*
        {
            // always start with the static bodies in the first position
            typedef STLEXT hash_map<dBodyID, size_t> BodyToIntHash;
            size_t iUnique = 0;
            BodyToIntHash bodies;
            bodies.insert( std::make_pair( dBodyID(0), iUnique++ ) );

            for( std::deque<dJointID>::const_iterator jointItr = contactJoints_.begin();
                jointItr != contactJoints_.end(); ++jointItr )
            {
                for( int i = 0; i < 2; ++i )
                {
                    dBodyID id = dJointGetBody( *jointItr, i );

                    BodyToIntHash::iterator itr = bodies.find( id );
                    if( itr == bodies.end() )
                        bodies.insert( std::make_pair(id, iUnique++) );
                }
            }

            // build the interaction graph
            typedef boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS> InteractionGraph;
            InteractionGraph G( iUnique );

            for( std::deque<dJointID>::const_iterator jointItr = contactJoints_.begin();
                jointItr != contactJoints_.end(); ++jointItr )
            {
                boost::array<size_t, 2> ids;
                for( int i = 0; i < 2; ++i )
                {
                    dBodyID id = dJointGetBody( *jointItr, i );

                    BodyToIntHash::iterator itr = bodies.find( id );
                    assert( itr != bodies.end() );
                    ids[i] = itr->second;
                }

                add_edge( ids[0], ids[1], G );
            }

            // now ensure that for any joints not connected to leaves, ReverseBounce is false
            for( std::deque<dJointID>::const_iterator jointItr = contactJoints_.begin();
                jointItr != contactJoints_.end(); ++jointItr )
            {
                boost::array<size_t, 2> ids;
                for( int i = 0; i < 2; ++i )
                {
                    dBodyID id = dJointGetBody( *jointItr, i );

                    BodyToIntHash::iterator itr = bodies.find( id );
                    assert( itr != bodies.end() );
                    ids[i] = itr->second;
                }

                size_t od1 = out_degree(ids[0], G);
                size_t od2 = out_degree(ids[1], G);

                if( (ids[0] == 0 || od1 > 1) && (ids[1] == 0 || od2 > 1) )
                {
                    dContact* contact = dJointGetContact(*jointItr);
                    contact->surface.mode &= ~dContactReverseBounce;
                    contact->surface.mode &= ~dContactReverseSlide;
                    contact->surface.mode &= ~dContactReverseRoll;
                    contact->surface.mode |= dContactReverseHold;
                }
            }
        }
        */

        std::deque< std::vector<vl::Vec3d> > contactNormals( this->dynamicBodies_.size() );
        for( std::deque<dJointID>::const_iterator jointItr = contactJoints_.begin();
            jointItr != contactJoints_.end(); ++jointItr )
        {
            dContact* contact = dJointGetContact(*jointItr);
            vl::Vec3d normal;
            for( int i = 0; i < 3; ++i )
                normal[i] = contact->geom.normal[i];

            for( int i = 0; i < 2; ++i )
            {
                dBodyID id = dJointGetBody( *jointItr, i );
                if( id == 0 )
                    continue;

                BodyToIntHash::iterator itr = bodyToId_.find( id );
                assert( itr != bodyToId_.end() );
                contactNormals[ itr->second ].push_back( normal );

                normal *= -1.0;
            }

        }

        boost::dynamic_bitset<> constrainedBodies( contactNormals.size() );
        for( size_t iBody = 0; iBody < contactNormals.size(); ++iBody )
        {
            const std::vector<vl::Vec3d>& normals = contactNormals.at(iBody);
            for( std::vector<vl::Vec3d>::const_iterator nItr = normals.begin();
                nItr != normals.end(); ++nItr )
            {
                for( std::vector<vl::Vec3d>::const_iterator nItr2 = normals.begin()+1;
                    nItr2 != normals.end(); ++nItr2 )
                {
                    if( vl::dot( *nItr, *nItr2 ) < 0 )
                        constrainedBodies[iBody] = true;
                }
            }
        }

        for( std::deque<dJointID>::const_iterator jointItr = contactJoints_.begin();
            jointItr != contactJoints_.end(); ++jointItr )
        {
            bool constrained = true;
            for( int i = 0; i < 2; ++i )
            {
                dBodyID id = dJointGetBody( *jointItr, i );
                if( id == 0 )
                    continue;

                BodyToIntHash::iterator itr = bodyToId_.find( id );
                assert( itr != bodyToId_.end() );
                constrained = constrained && constrainedBodies[ itr->second ];
            }

            if( constrained )
            {
                dContact* contact = dJointGetContact(*jointItr);
                contact->surface.mode &= ~dContactReverseBounce;
                contact->surface.mode &= ~dContactReverseSlide;
                contact->surface.mode &= ~dContactReverseRoll;
                contact->surface.mode |= dContactReverseHold;
            }
        }

#if 0
        {
            std::ofstream ofs( "debugContacts.bin", std::ios::binary );

            // We need to be able to figure out what is going on with the solver in some cases
            // we will dump out lots of useful information which we can use in visualization

            // First dump out all simulation state
            Simulation::State state = this->getSimulationState();
            float time = state.time();
    		ofs.write( reinterpret_cast<const char*>(&time), sizeof(float) );

            std::vector<RigidDynamicState> dynamicStates;
            for( size_t i = 0; i < state.stateCount(); ++i )
                dynamicStates.push_back( state.state(i) );
            dumpArray( dynamicStates, ofs );

            for( std::deque<dContactGeom>::const_iterator contactItr = contactGeoms_.begin();
                contactItr != contactGeoms_.end(); ++contactItr )
            {
                vl::Vec3f pos = vl::toVec3f( contactItr->pos );
                vl::Vec3f normal = vl::toVec3f( contactItr->normal );
                ofs.write( reinterpret_cast<const char*>( pos.Ref() ), sizeof( vl::Vec3f) );
                ofs.write( reinterpret_cast<const char*>( normal.Ref() ), sizeof( vl::Vec3f) );
            }
        }
#endif

		// now, "advance" the velocities
		dWorldStepBackwards( this->worldId_, -step );
	}
	else
	{
#if (0)
		for( std::deque<dGeomID>::const_iterator geomItr = triMeshGeoms_.begin();
			geomItr != triMeshGeoms_.end(); ++geomItr )
		{
			dBodyID body = dGeomGetBody( *geomItr );
			assert( body != 0 ); // we shouldn't have stored any static objects

			NxVec3 Pos( toNxVec3(dBodyGetPosition(body)) );
			NxVec3 linearVel( toNxVec3(dBodyGetLinearVel(body)) );
			Pos -= step*linearVel;

			NxVec3 angularVelocity( toNxVec3(dBodyGetAngularVel(body)) );
			NxVec3 diffRot = -step * angularVelocity;

			Quaternion<double> quat( 
				dBodyGetQuaternion(body), 
				Quaternion<double>::ORDER_WXYZ );
			Quaternion<double> diffRotQuat( toVec3f(diffRot) );
			Quaternion<dReal> rotation(
				diffRotQuat	* Quaternion<double>( quat ) );

			vl::Mat3d Rot = rotation.toRotMat();

			dReal p_matrix[16];
			p_matrix[ 0 ] = Rot[0][0];  p_matrix[ 1 ] = Rot[1][0];  p_matrix[ 2 ] = Rot[2][0];  p_matrix[ 3 ] = 0;
			p_matrix[ 4 ] = Rot[0][1];  p_matrix[ 5 ] = Rot[1][1];  p_matrix[ 6 ] = Rot[2][1];  p_matrix[ 7 ] = 0;
			p_matrix[ 8 ] = Rot[0][2];  p_matrix[ 9 ] = Rot[1][2];  p_matrix[10 ] = Rot[2][2];  p_matrix[11 ] = 0;
			p_matrix[12 ] = Pos[0];     p_matrix[13 ] = Pos[1];     p_matrix[14 ] = Pos[2];     p_matrix[15 ] = 1;
			dGeomTriMeshSetLastTransform( *geomItr, p_matrix );
		}
#endif

		for( std::deque<dBodyID>::const_iterator bodyItr = this->dynamicBodies_.begin();
			bodyItr != this->dynamicBodies_.end(); ++bodyItr )
		{
            const dReal* linearVel = dBodyGetLinearVel(*bodyItr);
            const dReal dampingCoeff = -0.1;
            dBodyAddForce( *bodyItr,
                dampingCoeff*linearVel[0],
                dampingCoeff*linearVel[1],
                dampingCoeff*linearVel[2] );
        }

		for( std::deque<JointWithId>::const_iterator jointItr = activeJoints_.begin();
			jointItr != activeJoints_.end(); ++jointItr )
		{
			(jointItr->first)->updateControllers( jointItr->second, time(), step );
		}

		// clear out all contact joints:
		std::for_each( contactJoints_.begin(), contactJoints_.end(),
			&dJointDestroy );
		contactJoints_.clear();
		contactGeoms_.clear();

		// perform collision detection
		this->handleCollisions(false);

		// now need to mark contact geoms as needing feedback
		contactJointFeedback_.resize( contactJoints_.size() );
		for( size_t iJoint = 0; iJoint < contactJoints_.size(); ++iJoint )
			dJointSetFeedback( contactJoints_[iJoint], &contactJointFeedback_[iJoint] );

		// take a step:
		//dRandSetSeed(0);
		dWorldQuickStep( this->worldId_, step );
		//dWorldStep( this->worldId_, step );
	}

	setTime( time() + step );

    std::for_each( this->planeJointBodies_.begin(), this->planeJointBodies_.end(), 
        &setPlaneAngularVel );

	updateJoints( time() );

    /*
    {
        static size_t step = 0;
        std::ostringstream oss;
        oss << "stateDump" << std::setfill('0') << std::setw( 4 ) << step++ << ".bin";
        this->dumpState( oss.str().c_str() );
    }

    size_t numDynOb = dynamicObjects().size();
    for( size_t iObject = 0; iObject < numDynOb; ++iObject )
    {
        vl::Vec3f vel = this->getLinearVelocity(iObject);
        assert( fabs(vel[1]) < 1.0 );
    }
    */
}

#ifdef SOUND
std::vector<double> ODESimulation::stepAudio(const Simulation::State& previousState, 
											 float previousTimestep, size_t numSteps)
{
	const float epsilon = 1e-4;

	// first apply impulses from the previous timestep
	for( size_t iContact = 0; iContact < contactGeoms_.size(); ++iContact )
	{
		const dContactGeom& contactGeom = contactGeoms_.at(iContact);
		const dJointFeedback& feedback = contactJointFeedback_.at(iContact);
		boost::array<dGeomID, 2> geoms = {{ contactGeom.g1, contactGeom.g2 }};
		boost::array<dBodyID, 2> bodies = {{ dJointGetBody(contactJoints_[iContact], 0),
			dJointGetBody(contactJoints_[iContact], 1) }};

		// need to figure out if there is relative motion between the bodies to determine
		// whether it is necessary to apply an impulse
		boost::array<vl::Vec3f, 2> contactPointVelocities;
		for( size_t i = 0; i < 2; ++i )
		{
			dBodyID body = dJointGetBody(contactJoints_[iContact], i);
			if( body == 0 )
			{
				contactPointVelocities[i] = vl::vl_0;
				continue;
			}

			vl::Vec3f cmPosition = vl::toVec3f( dBodyGetPosition(bodies[i]) );
			vl::Vec3f linVel = vl::toVec3f( dBodyGetLinearVel(bodies[i]) );
			vl::Vec3f angVel = vl::toVec3f( dBodyGetAngularVel(bodies[i]) );
			vl::Vec3f contactPoint = vl::toVec3f( contactGeom.pos );

			contactPointVelocities[i] = linVel + 
				vl::cross(angVel, contactPoint - cmPosition);
		}
		vl::Vec3f diff = contactPointVelocities[1] - contactPointVelocities[0];
		vl::Vec3f normal = vl::toVec3f( contactGeom.normal );
		float normalVel = vl::dot( normal, diff );
		if( fabs(normalVel) < epsilon )
			continue;

		vl::Vec3f pos = vl::toVec3f( contactGeom.pos );

		boost::array<vl::Vec3f, 2> forces = {{ 
			vl::Vec3f( feedback.f1[0], feedback.f1[1], feedback.f1[2] ),
			vl::Vec3f( feedback.f2[0], feedback.f2[1], feedback.f2[2] ) }};

		for( size_t i = 0; i < 2; ++i )
		{
			const GeomInfo *info = reinterpret_cast<const GeomInfo*>( 
				dGeomGetData( geoms[i] ) );

			ObjectToAudioHash::const_iterator iter = 
				audioGenerators_.find( info->originalObject );
			if( iter == audioGenerators_.end() )
				continue;
			
			// @todo this is clearly not the right thing to do about static
			//   objects
			int iDynId = info->iDynamicId;
			vl::Mat4f transform( vl::vl_1 );
			if( iDynId >= 0 )
				transform = info->originalObject->rigidTransform( previousState.extendedState(iDynId) );

			vl::Mat3f rotate = toMat3f( transform );
			vl::Vec3f transformedPos = vl::xform( rigidInverse( transform ), pos );
			// for rigid transform, inv(trans(mat)) == mat
			vl::Vec3f normalForce = vl::dot( forces[i], normal ) * normal;
			vl::Vec3f transformedForce = rotate * normalForce;

			iter->second->addImpulse( transformedPos, transformedForce, previousTimestep );
		}
	}

	return this->stepAudio( numSteps );
}

std::vector<double> ODESimulation::stepAudio( size_t numSteps )
{
	std::vector<double> result( numSteps, 0.0 );
	for( ObjectToAudioHash::iterator itr = audioGenerators_.begin();
		itr != audioGenerators_.end(); ++itr )
	{
		boost::shared_ptr<AudioGenerator> audioGenerator = itr->second;
		assert( audioGenerator );
		std::vector<double> current = audioGenerator->integrate( numSteps );
		assert( current.size() == result.size() );
		std::transform( current.begin(), current.end(), result.begin(), result.begin(),
			std::plus<double>() );
	}

	return result;
}
#endif

void ODESimulation::addContactJoint(dBodyID b1, dBodyID b2, const dContact& contact)
{
    if( (contact.surface.mode & dContactReverseBounce) 
     || (contact.surface.mode & dContactReverseSlide)
     || (contact.surface.mode & dContactReverseRoll) )
    {
        if( b1 ) dBodyEnable( b1 );
        if( b2 ) dBodyEnable( b2 );
    }

	dJointID c = dJointCreateContact(this->worldId_, 0, &contact);
	dJointAttach(c, b1, b2);

	contactJoints_.push_back( c );
	contactGeoms_.push_back( contact.geom );
}

std::vector<Simulation::ContactPoint> ODESimulation::contactPoints() const
{
	std::vector<ContactPoint> result;
	result.reserve( contactGeoms_.size() );
	for( std::deque<dContactGeom>::const_iterator iter = this->contactGeoms_.begin();
		iter != this->contactGeoms_.end(); ++iter )
	{
        const GeomInfo* info1 = reinterpret_cast<const GeomInfo*>( dGeomGetData(iter->g1) );
        const GeomInfo* info2 = reinterpret_cast<const GeomInfo*>( dGeomGetData(iter->g2) );
        result.push_back( ContactPoint(vl::toVec3f(iter->pos), vl::toVec3f(iter->normal), 
            info1->iDynamicId, info2->iDynamicId) );
	}

	return result;
}

dBodyID ODESimulation::odeBodyId( size_t iObject )
{
	return this->dynamicBodies_.at( iObject );
}

ODENativeSimulation::ODENativeSimulation(boost::shared_ptr<const Scene> scene)
	: ODESimulation(scene)
{
    vl::Vec3d center( 0, 0, 0 );
    vl::Vec3d extents( toVec3d(scene->bounds().maximum() - scene->bounds().minimum()) );
	spaceId_ = dHashSpaceCreate(0);
	//spaceId_ = dQuadTreeSpaceCreate(0, &center[0], &extents[0], 4);

	this->initialize(scene);
}

ODENativeSimulation::~ODENativeSimulation()
{
    std::for_each( allGeoms_.begin(), allGeoms_.end(), &dGeomDestroy );
    std::for_each( allSpaces_.begin(), allSpaces_.end(), &dSpaceDestroy );

	dSpaceDestroy( spaceId_ );
}

void removeGeomsFromSpace( dSpaceID space, dBodyID body )
{
    std::deque<dGeomID> geomsToRemove;
    std::deque<dSpaceID> spaces;

    for( int i = 0; i < dSpaceGetNumGeoms(space); ++i )
    {
        dGeomID geom = dSpaceGetGeom( space, i );
        if( dGeomIsSpace(geom) )
            spaces.push_back( reinterpret_cast<dSpaceID>(geom) );
        else if( dGeomGetBody(geom) == body )
            geomsToRemove.push_back( geom );
    }

    std::for_each( geomsToRemove.begin(), geomsToRemove.end(), 
        boost::bind( &dGeomDestroy, _1 ) );
    for( std::deque<dSpaceID>::const_iterator spaceItr = spaces.begin();
        spaceItr != spaces.end(); ++spaceItr )
    {
        removeGeomsFromSpace( *spaceItr, body );
        if( dSpaceGetNumGeoms(*spaceItr) == 0 )
            dSpaceDestroy( *spaceItr );
    }
}

void ODENativeSimulation::removeCollisionObject( int dynamicId, dBodyID bodyId )
{
    assert( dynamicId >= 0 );

    // For now, we're not sure what will happen if we remove the body itself
    // so we'll just remove all the geoms associated with it
    removeGeomsFromSpace( this->spaceId_, bodyId );
}

void ODENativeSimulation::addCollisionObject( 
	boost::shared_ptr<const PhysicsObject> object,
	int iDynamicId,
	dBodyID bodyId,
	const vl::Vec3f& centerOfMass,
	float padding,
	bool isGhost )
{
	ODEGeomResult geomResult = object->odeGeom(0, centerOfMass, padding);
#if dTRIMESH_ENABLED
	std::copy( geomResult.triMeshes.begin(), geomResult.triMeshes.end(),
		std::back_inserter(this->meshes_) );
#endif

	vl::Mat4f transform = object->rigidTransform();
	vl::Vec3f position = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	vl::Mat4f odeRotation = transform;
	for( vl::Int i = 0; i < 3; ++i )
		odeRotation[i][3] = 0.0f;

    std::copy( geomResult.geoms.begin(), geomResult.geoms.end(), std::back_inserter(allGeoms_) );

    const size_t minGeomsForSpace = 5;
	if( geomResult.geoms.size() > minGeomsForSpace )
	{
		dSpaceID currentSpaceId = dHashSpaceCreate( this->spaceId_ );
		std::for_each( geomResult.geoms.begin(), geomResult.geoms.end(),
			boost::bind( &dSpaceAdd, currentSpaceId, _1 ) );
        allSpaces_.push_back( currentSpaceId );
	}
	else
	{
		std::for_each( geomResult.geoms.begin(), geomResult.geoms.end(),
			boost::bind( &dSpaceAdd, this->spaceId_, _1 ) );
	}

	GeomInfo* info = geomInfoPool_.construct( *(object->material()) );
	info->ghost = isGhost;
	info->originalObject = object.get();
	info->iDynamicId = iDynamicId;
	info->material = *object->material();

	for( std::vector<dGeomID>::const_iterator geomIter = geomResult.geoms.begin(); 
		geomIter != geomResult.geoms.end(); ++geomIter )
	{
		dGeomSetData( *geomIter, info );

		if( bodyId == 0 ) // not dynamic
		{
			if( !placeable(*geomIter) )
				continue;

			std::vector<dReal> rotation( odeRotation.Ref(), odeRotation.Ref() + 12 );
			dGeomSetPosition( *geomIter, position[0], position[1], position[2] );
			// should be in the correct form; ODE expects row-major 3x4 matrix
			dGeomSetRotation( *geomIter, &rotation[0] );

#if (0)
			if( dGeomGetClass(*geomIter) == dTriMeshClass )
			{
				// for static objects, we only need to set the previous transform once:
				const dReal* Pos = dGeomGetPosition(*geomIter);
				const dReal* Rot = dGeomGetRotation(*geomIter);

				dReal p_matrix[16];
				p_matrix[ 0 ] = Rot[ 0 ];	p_matrix[ 1 ] = Rot[ 1 ];	p_matrix[ 2 ] = Rot[ 2 ];	p_matrix[ 3 ] = 0;
				p_matrix[ 4 ] = Rot[ 4 ];	p_matrix[ 5 ] = Rot[ 5 ];	p_matrix[ 6 ] = Rot[ 6 ];	p_matrix[ 7 ] = 0;
				p_matrix[ 8 ] = Rot[ 8 ];	p_matrix[ 9 ] = Rot[ 9 ];	p_matrix[10 ] = Rot[10 ];	p_matrix[11 ] = 0;
				p_matrix[12 ] = Pos[ 0 ];	p_matrix[13 ] = Pos[ 1 ];	p_matrix[14 ] = Pos[ 2 ];	p_matrix[15 ] = 1;

				dGeomTriMeshSetLastTransform( *geomIter, p_matrix );
			}
#endif
		}
		else
		{
			dGeomSetBody( *geomIter, bodyId );

			// we need to remember which ones are trimeshes because later on we'll
			//   need to set the prevTransform at each time step: 
			if( dGeomGetClass(*geomIter) == dTriMeshClass )
				triMeshGeoms_.push_back( *geomIter );
		}
	}
}

void potentialCollision( void* data, dGeomID o1, dGeomID o2 )
{
	ODENativeSimulation* simulation = reinterpret_cast<ODENativeSimulation*>( data );
	assert( simulation != 0 );

	simulation->checkContact( o1, o2 );
}

void ODENativeSimulation::checkContact( dGeomID o1, dGeomID o2 )
{
	bool isSpace1 = dGeomIsSpace(o1);
	bool isSpace2 = dGeomIsSpace(o2);

	const GeomInfo *info1, *info2;
	if( isSpace1 )
		info1 = reinterpret_cast<const GeomInfo*>( 
			dGeomGetData( dSpaceGetGeom(reinterpret_cast<dSpaceID>(o1), 0) ) );
	else
		info1 = reinterpret_cast<const GeomInfo*>( dGeomGetData(o1) );

	if( isSpace2 )
		info2 = reinterpret_cast<const GeomInfo*>( 
			dGeomGetData( dSpaceGetGeom(reinterpret_cast<dSpaceID>(o2), 0) ) );
	else
		info2 = reinterpret_cast<const GeomInfo*>( dGeomGetData(o2) );

	if( this->ignoreCollisions( info1->originalObject, info2->originalObject ) )
		return;

	bool active1 = isActive( info1->originalObject );
	bool active2 = isActive( info2->originalObject );

#ifdef _DEBUG
	// If the objects are different, they'd better both be ACTIVE
	// If they are in correspondence, then they are supposed to get added to the ACTIVE
	//   set at the same time
	if( info1->originalObject->name() == info2->originalObject->name() )
		assert( active1 == active2 );
#endif

	// If they're both either active or inactive, we can skip colliding them if
	//   one of them is a ghost
	if( (active1 == active2) && (info1->ghost || info2->ghost) )
		return;

	// No point colliding ghost objects against each other
	if( info1->ghost && info2->ghost )
		return;

	if( isSpace1 || isSpace2 )
	{
		dSpaceCollide2( o1, o2, this, potentialCollision );
	}
	else
	{
		dBodyID b1 = dGeomGetBody(o1);
		dBodyID b2 = dGeomGetBody(o2);

		const Material& mat1 = info1->material;
		const Material& mat2 = info2->material;

		const int maxContacts = 128;

		dContact contacts[ maxContacts ];

		// not worth getting multiple contacts for ghost objects, since
		//   1 is enough to activate them.
		int maxContactsTemp = maxContacts;
		if( info1->ghost || info2->ghost )
			maxContactsTemp = 1;

		int numc = dCollide( o1, o2, maxContactsTemp, &contacts[0].geom, sizeof(dContact) );
		assert( numc >= 0 );
		if( numc == 0 )
			return;

		// non-ghost objects don't need to activate ghost objects
		if( active1 && !active2 && !info2->ghost && (b2 != 0) )
		{
			dBodyEnable(b2);
		}
		else if( active2 && !active1 && !info1->ghost && (b1 != 0) )
		{
			dBodyEnable(b1);
		}

        const dReal epsilon = 1e-1;
        bool moving = false;
        if( b1 && (dDOT( dBodyGetLinearVel( b1 ), dBodyGetLinearVel( b1 ) ) > epsilon
                || dDOT( dBodyGetAngularVel( b1 ), dBodyGetAngularVel( b1 ) ) > epsilon) )
               moving = true;

        if( b2 && (dDOT( dBodyGetLinearVel( b2 ), dBodyGetLinearVel( b2 ) ) > epsilon
                || dDOT( dBodyGetAngularVel( b2 ), dBodyGetAngularVel( b2 ) ) > epsilon) )
               moving = true;


		// won't collide ghost objects
		if( info1->ghost || info2->ghost )
		{
			return;
		}

        // need to check various properties for sampling
        SamplingProperties props1, props2;
        if( info1->iDynamicId >= 0 )
            props1 = this->getSamplingProperties( info1->iDynamicId );
        if( info2->iDynamicId >= 0 )
            props2 = this->getSamplingProperties( info2->iDynamicId );

        bool bounce = props1.bounce || props2.bounce;
        bool slide = props1.slide || props2.slide;
        bool roll = false;

        vl::Vec3d frictionDir;
        if( props1.frictionDirection != vl::vl_0 && props2.frictionDirection != vl::vl_0 )
        {
            // randomly pick one
            frictionDir = vl::norm( toVec3d(
                (randomStream_->uniform(0, 2) == 0) ? -props1.frictionDirection : props2.frictionDirection ) ); 
        }
        if( props1.frictionDirection != vl::vl_0 )
        {
            frictionDir = vl::norm( toVec3d(props2.frictionDirection) );
        }
        else if( props2.frictionDirection != vl::vl_0 )
        {
            frictionDir = -vl::norm( toVec3d(props1.frictionDirection) );
        }
        else
        {
            frictionDir = vl::vl_0;
        }

        double slidingLength = 0.2;
        if( backwards_ && (info1->lastSlidingTime - this->time() < slidingLength)
                       && (info2->lastSlidingTime - this->time() < slidingLength) )
        {
            slide = true;
            if( info1->iDynamicId >= -1 )
                frictionDir = info1->lastFrictionDir;
            else
                frictionDir = info2->lastFrictionDir;
        }


        double rollingLength = 2.0;
        if( backwards_ && (info1->lastRollingTime - this->time() < rollingLength)
                       && (info2->lastRollingTime - this->time() < rollingLength) )
        {
            roll = true;
            if( info1->iDynamicId >= -1 )
                frictionDir = info1->lastFrictionDir;
            else
                frictionDir = info2->lastFrictionDir;
        }

        double timeBeforeSlide = 3.0;
        if( !moving && randomStream_->poisson( stepSize_/timeBeforeSlide ) > 0 )
        {
            // randomize the slide direction
            std::vector<double> dir = randomStream_->uniformSpherical(3);
            std::copy( dir.begin(), dir.end(), frictionDir.Ref() );
            frictionDir = vl::norm( frictionDir );
            slide = true;

            info1->lastSlidingTime = this->time();
            info2->lastSlidingTime = this->time();

            info1->lastFrictionDir = frictionDir;
            info2->lastFrictionDir = frictionDir;
        }

        /*
        double timeBeforeRoll = 5.0;
        if( !moving && randomStream_->poisson( stepSize_/timeBeforeRoll ) > 0 )
        {
            // randomize the slide direction
            std::vector<double> dir = randomStream_->uniformSpherical(3);
            std::copy( dir.begin(), dir.end(), frictionDir.Ref() );
            frictionDir = vl::norm( frictionDir );
            roll = true;

            info1->lastRollingTime = this->time();
            info2->lastRollingTime = this->time();

            info1->lastFrictionDir = frictionDir;
            info2->lastFrictionDir = frictionDir;
        }
        */

        size_t bounceEvents = 0;
        if( this->backwards_ /* && moving*/ && (randomStream_->uniform() < this->stepSize_*120.0) )
        {
            fortran_matrix tm( 4, 4 );
            tm(0, 0) = 0.9669; tm(0, 1) = 0.4276; tm(0, 2) = 0.0911; tm(0, 3) = 0.0004;
            tm(1, 0) = 0.0295; tm(1, 1) = 0.4903; tm(1, 2) = 0.1877; tm(1, 3) = 0.0013;
            tm(2, 0) = 0.0036; tm(2, 1) = 0.0796; tm(2, 2) = 0.6619; tm(2, 3) = 0.0107;
            tm(3, 0) = 0.0000; tm(3, 1) = 0.0025; tm(3, 2) = 0.0594; tm(3, 3) = 0.9875;
            double randomVal = randomStream_->uniform();

            size_t nextContacts = 0;
            int origContacts = std::min(numc, 3);

            double curSum = 0;
            while( nextContacts <= 3 )
            {
                curSum += tm(nextContacts, origContacts);
                if( randomVal < curSum )
                    break;

                nextContacts += 1;
            }

            if( nextContacts < origContacts  )
                bounceEvents = origContacts - nextContacts;
        }
        const size_t origBounceEvents = bounceEvents;

        double r1 = mat1.restitution();
        double r2 = mat2.restitution();

        const double collisionLife = 0.2;
        const double recentCollisionRest = 1.0;

        /*
        if( backwards_ && info1->lastCollisionTime > this->time() )
        {
            double sinceLastCollision = -this->time() + info1->lastCollisionTime;
            if( sinceLastCollision < collisionLife )
            {
                double weight = 1.0 - sinceLastCollision/collisionLife;
                r1 = weight*recentCollisionRest + (1.0-weight)*r1;
            }
        }

        if( backwards_ && info2->lastCollisionTime > this->time() )
        {
            double sinceLastCollision = -this->time() + info2->lastCollisionTime;
            if( sinceLastCollision < collisionLife )
            {
                double weight = 1.0 - sinceLastCollision/collisionLife;
                r2 = weight*recentCollisionRest + (1.0-weight)*r2;
            }
        }
        */

        info1->lastCollisionTime = this->time();
        info2->lastCollisionTime = this->time();

        double restitution = 0.5*(r1+r2);
        if( backwards_ && this->connected( info1->originalObject, info2->originalObject ) )
            restitution = 1.0;

        for( int i = 0; i < numc; ++i )
		{
			contacts[i].surface.mode = dContactBounce;
			if( mat1.coulombFriction() || mat2.coulombFriction() )
            {
				contacts[i].surface.mode |= dContactApprox1;

			    //contacts[i].surface.mode |= dContactRolling;
                contacts[i].surface.rollingFriction = 0.02;
            }

			// is max appropriate here?  also should switch between dynamic and static friction
			contacts[i].surface.mu = std::max( mat1.dynamicFriction(), mat2.dynamicFriction() );

			contacts[i].surface.bounce = restitution;
			contacts[i].surface.bounce_vel = 0.05f;

            /*
            vl::Vec3f vel_b1 = vl::vl_0;
            if( b1 )
            {
                dVector3 pointLocal;
                dBodyGetPosRelPoint( b1, contacts[i].geom.pos[0], contacts[i].geom.pos[1], contacts[i].geom.pos[2], pointLocal );
                dVector3 result;
                dBodyGetPointVel( b1, pointLocal[0], pointLocal[1], pointLocal[2], result );
                for( vl::Int i = 0; i < 3; ++i )
                    vel_b1[i] = result[i];
            }

            vl::Vec3f vel_b2 = vl::vl_0;
            if( b2 )
            {
                dVector3 pointLocal;
                dBodyGetPosRelPoint( b2, contacts[i].geom.pos[0], contacts[i].geom.pos[1], contacts[i].geom.pos[2], pointLocal );
                dVector3 result;
                dBodyGetPointVel( b2, pointLocal[0], pointLocal[1], pointLocal[2], result );
                for( vl::Int i = 0; i < 3; ++i )
                    vel_b2[i] = result[i];
            }

            vl::Vec3f relVel = vel_b1 - vel_b2;
            if( vl::len(relVel) > 1e-3 )
                frictionDir = vl::norm( relVel );
            */

            if( backwards_ )
            {
                if( bounce )
                {
                    contacts[i].surface.mode |= dContactReverseBounce;
                }
                if( randomStream_->uniform() < static_cast<double>(bounceEvents) / static_cast<double>( numc-i ) )
                {
                    contacts[i].surface.mode |= dContactReverseBounce;
                    --bounceEvents;
                }

                if( slide )
                    contacts[i].surface.mode |= dContactReverseSlide;

                /*
                if( roll )
                    contacts[i].surface.mode |= dContactReverseRoll;
                    */

                if( frictionDir != vl::vl_0 )
                {
                    vl::Vec3d fdir = toVec3d( randomDirection( toVec3f(frictionDir), randomStream_->uniform(), randomStream_->uniform(), 20 ) );

                    vl::Vec3d norm( 
                        contacts[i].geom.normal[0], 
                        contacts[i].geom.normal[1], 
                        contacts[i].geom.normal[2] );
                    norm = vl::norm( norm );
                    double dotProd = vl::dot( fdir, norm );
                    fdir -= dotProd*norm;

                    double length = vl::len( fdir );
                    if( length > 1e-3 )
                    {
                        for( vl::Int j = 0; j < 3; ++j )
                            contacts[i].fdir1[j] = fdir[j] / length;
                        contacts[i].surface.mode |= dContactFDir1;
                    }
                }
            }

			this->addContactJoint( b1, b2, contacts[i] );
		}
	}
}

void ODENativeSimulation::handleCollisions(bool backwards)
{
	// perform collision detection
	dSpaceCollide( this->spaceId_, this, potentialCollision );
}

// For Newton and Bullet simulations, the inertia tensor must be in diagonal form
std::pair< vl::Vec3d, vl::Mat4d > diagonalizeInertiaTensor( dMass mass )
{
	vl::Vec3d centerOfMass( mass.c[0], mass.c[1], mass.c[2] );
	dMassTranslate( &mass, -mass.c[0], -mass.c[1], -mass.c[2] );

	// check if the inertia is already diagonalized:
	{
		const dReal epsilon = 1e-4;

		bool diagonal = true;
		for( int i = 0; i < 3; ++i )
			for( int j = 0; j < 3; ++j )
				diagonal = diagonal && ((i == j) || (fabs(mass.I[ i*4 + j ]) < epsilon));

		if( diagonal )
		{
			vl::Mat4d transform = vl::vl_1;
			for( vl::Int i = 0; i < 3; ++i )
				transform[i][3] = centerOfMass[i];

			vl::Vec3d inertia;
			for( vl::Int i = 0; i < 3; ++i )
				inertia[i] = mass.I[ i*4 + i ];

			return std::make_pair( inertia, transform );
		}
	}



	const int workSize = 512;
	fortran_matrix a( 3, 3 );
	for( int i = 0; i < 3; ++i )
		for( int j = 0; j < 3; ++j )
			a(i, j) = mass.I[ i*4 + j ];

	std::vector<double> eigs;
	fortran_matrix V;
	boost::tie( V, eigs ) = symmetricJacobi( a );
	vl::Vec3d I( eigs[0], eigs[1], eigs[2] );

	vl::Vec3d a1( V(0, 0), V(1, 0), V(2, 0) );
	vl::Vec3d a2( V(0, 1), V(1, 1), V(2, 1) );
	vl::Vec3d a3( V(0, 2), V(1, 2), V(2, 2) );

	// need to ensure a right-handed coordinate system, which is
	// not guaranteed by dsyev:
	vl::Vec3d new_a3 = vl::cross( a1, a2 );
#ifdef _DEBUG
	vl::Vec3d tmp = vl::cross( a3, new_a3 );
	assert( vl::len(tmp) < 1e-4 );
	// check normalized
	assert( fabs(vl::len(new_a3) - 1.0f) < 1e-4 );
#endif


#if 0
	{
		// check whether we match lapack:
		char jobz = 'V';
		char uplo = 'U';
		int n = 3;
		int lda = 3;
		int lwork = workSize;
		double work[workSize];
		std::vector<double> eigsLapack(3);
		int info;

		dsyev(
			&jobz,
			&uplo,
			&n,
			a.data(),
			&lda,
			&eigsLapack[0],
			work,
			&lwork,
			&info );

		std::sort( eigs.begin(), eigs.end() );
		std::sort( eigsLapack.begin(), eigsLapack.end() );

		for( size_t i = 0; i < 3; ++i )
			assert( withinEpsilon( eigs[0], eigsLapack[0] ) );
	}
#endif

	vl::Mat4d geomTransform( vl::vl_1 );
	for( vl::Int i = 0; i < 3; ++i )
		geomTransform[i][0] = a1[i];

	for( vl::Int i = 0; i < 3; ++i )
		geomTransform[i][1] = a2[i];

	for( vl::Int i = 0; i < 3; ++i )
		geomTransform[i][2] = new_a3[i];

	for( vl::Int i = 0; i < 3; ++i )
		geomTransform[i][3] = centerOfMass[i];

	return std::make_pair( I, geomTransform );
}

#ifdef USE_BULLET
class MyBulletCollisionDispatcher
	: public btCollisionDispatcher
{
	typedef ODESimulation::GeomInfo GeomInfo;

public:
	MyBulletCollisionDispatcher(btCollisionConfiguration* collisionConfig, const Simulation* simulation)
		: btCollisionDispatcher(collisionConfig), simulation_(simulation) {}

	~MyBulletCollisionDispatcher() {}

	bool needsCollision(btCollisionObject* body1, btCollisionObject* body2)
	{
		const GeomInfo* info1 = 
			reinterpret_cast<const GeomInfo*>(body1->getUserPointer());
		const GeomInfo* info2 = 
			reinterpret_cast<const GeomInfo*>(body2->getUserPointer());

        // don't collide static objects with each other
        if( info1->iDynamicId < 0 && info2->iDynamicId < 0 )
            return false;

		bool active1 = simulation_->isActive( info1->originalObject );
		bool active2 = simulation_->isActive( info2->originalObject );

		if( simulation_->ignoreCollisions( info1->originalObject, info2->originalObject ) )
			return false;

#ifdef _DEBUG
		// If the objects are different, they'd better both be ACTIVE
		// If they are in correspondence, then they are supposed to get added to the ACTIVE
		//   set at the same time
		if( info1->originalObject->name() == info2->originalObject->name() )
			assert( active1 == active2 );
#endif

		// If they're both either active or inactive, we can skip colliding them if
		//   one of them is a ghost
		if( (active1 == active2) && (info1->ghost || info2->ghost) )
			return false;

		// No point colliding ghost objects against each other
		if( info1->ghost && info2->ghost )
			return false;

		return true;
	}

private:
	const Simulation* simulation_;
};


ODEBulletSimulation::ODEBulletSimulation(boost::shared_ptr<const Scene> scene)
	: ODESimulation( scene )
{
	collisionDispatcher_.reset( new MyBulletCollisionDispatcher(&collisionConfig_, this) );

	// todo figure out what these should be
	const size_t maxProxies = 8192;
	btVector3 worldAabbMin(-10000,-10000,-10000);
	btVector3 worldAabbMax(10000,10000,10000);

	broadPhase_.reset( new btAxisSweep3(worldAabbMin,worldAabbMax,maxProxies) );

	collisionWorld_.reset( new btCollisionWorld(
		collisionDispatcher_.get(), 
		broadPhase_.get()) );

	this->initialize( scene );
}

ODEBulletSimulation::~ODEBulletSimulation()
{
}

void ODEBulletSimulation::addCollisionObject( 
	boost::shared_ptr<const PhysicsObject> object,
	int iDynamicId,
	dBodyID bodyId,
	const vl::Vec3f& centerOfMass,
	float padding,
	bool isGhost )
{
	vl::Mat4f geomTransform( vl::vl_1 );
	for( vl::Int i = 0; i < 3; ++i )
		geomTransform[i][3] = centerOfMass[i];

	btCollisionShape* shape = object->bulletCollisionShape(geomTransform, 0.0f);

	if( shape == 0 )
		return;

	{
		// need to make sure all the child shapes get destroyed when the class is destroyed
		btCompoundShape* compoundShape = 
			dynamic_cast<btCompoundShape*>( shape );
		if( compoundShape )
		{
			for( int i = 0; i < compoundShape->getNumChildShapes(); ++i )
				collisionShapes_.push_back( compoundShape->getChildShape(i) );
		}
	}
	collisionShapes_.push_back( shape );

	btCollisionObject* collObj = new btCollisionObject;
	collObj->setCollisionShape( shape );

	// need to set the appropriate world transform
	// for static objects this will be permanent; for dynamic objects
	//   it will get re-set in the collision handling code
	vl::Mat4f transform = object->rigidTransform();
	collObj->setWorldTransform( toBtTransform(transform) );

	// We will use this structure to get information about the object
	//   during collision processing
	GeomInfo* info = geomInfoPool_.construct( *(object->material()) );
	info->ghost = isGhost;
	info->originalObject = object.get();
	info->iDynamicId = iDynamicId;
	info->material = *object->material();
	collObj->setUserPointer( info );

	if( iDynamicId < 0 )
	{
		this->staticRigidBodies_.push_back( collObj );
	}
	else
	{
		assert( this->dynamicRigidBodies_.size() == iDynamicId );
		this->dynamicRigidBodies_.push_back( collObj );
	}

	collisionWorld_->addCollisionObject( collObj );
}

void ODEBulletSimulation::removeCollisionObject( int dynamicId, dBodyID bodyId )
{
    assert( dynamicId >= 0 );

    // For now, we're not sure what will happen if we remove the body itself
    // so we'll just remove all the geoms associated with it
    collisionWorld_->removeCollisionObject( &this->dynamicRigidBodies_.at( dynamicId ) );
}

// borrowed from test_BulletOde.cpp:
btTransform	GetTransformFromOde(const dReal* pos, const dReal* rot)
{
	btTransform trans;
	trans.setIdentity();

// rot is pointer to object's rotation matrix, 4*3 format!
	
	btMatrix3x3 orn(rot[0],rot[1],rot[2],
		rot[4],rot[5],rot[6],
		rot[8],rot[9],rot[10]);

	trans.setOrigin(btVector3(pos[0],pos[1],pos[2]));
	trans.setBasis(orn);

  return trans;
}

void ODEBulletSimulation::handleCollisions(bool backwards)
{
	// first step is to update all the Bullet collision shapes with transforms
	//   from their ODE versions:
	for( size_t iDynObj = 0; iDynObj < this->dynamicRigidBodies_.size(); ++iDynObj )
	{
		dBodyID id = this->odeBodyId(iDynObj);
		dynamicRigidBodies_.at(iDynObj).setWorldTransform( 
			GetTransformFromOde(dBodyGetPosition(id),dBodyGetRotation(id)) );
	}

	// Now, run collision detection
	this->collisionWorld_->performDiscreteCollisionDetection();

	// copy the various contact points over to ODE
	// based on code in test_BulletODE.cpp
	for( int i=0; i < collisionWorld_->getDispatcher()->getNumManifolds(); ++i)
	{
		btPersistentManifold* manifold = collisionWorld_->getDispatcher()->getManifoldByIndexInternal(i);
		btCollisionObject* obj1 = static_cast<btCollisionObject*>(manifold->getBody0());
		btCollisionObject* obj2 = static_cast<btCollisionObject*>(manifold->getBody1());

		const GeomInfo* info1 = reinterpret_cast<const GeomInfo*>( obj1->getUserPointer() );
		const GeomInfo* info2 = reinterpret_cast<const GeomInfo*>( obj2->getUserPointer() );

		bool active1 = this->isActive( info1->originalObject );
		bool active2 = this->isActive( info2->originalObject );

		const Material& mat1 = info1->material;
		const Material& mat2 = info2->material;

        double restitution = std::max( mat1.restitution(), mat2.restitution() );
        if( this->connected( info1->originalObject, info2->originalObject ) )
            restitution = 1.0;

		dBodyID b1 = 0;
		dBodyID b2 = 0;
		if( info1->iDynamicId >= 0 )
			b1 = this->odeBodyId( info1->iDynamicId );
		if( info2->iDynamicId >= 0 )
			b2 = this->odeBodyId( info2->iDynamicId );

		// non-ghost objects don't need to activate ghost objects
		if( active1 && !active2 && !info2->ghost && (b2 != 0) )
			dBodyEnable(b2);
		else if( active2 && !active1 && !info1->ghost && (b1 != 0) )
			dBodyEnable(b1);

		// won't collide ghost objects
		if( info1->ghost || info2->ghost )
			continue;

		//refreshContactPoints will update and/or remove existing contactpoints from previous frames
		manifold->refreshContactPoints(obj1->getWorldTransform(), obj2->getWorldTransform());
		for( int j = 0; j < manifold->getNumContacts(); ++j )
		{
			btManifoldPoint& pt = manifold->getContactPoint(j);
			if( pt.getDistance() < 0.0 )
			{
				dContact contact;

				contact.surface.mode = dContactBounce;
				if( mat1.coulombFriction() || mat2.coulombFriction() )
					contact.surface.mode |= dContactApprox1;

				// is max appropriate here?  also should switch between dynamic and static friction
				contact.surface.mu = std::max( mat1.dynamicFriction(), mat2.dynamicFriction() );

				contact.surface.bounce = restitution;
				contact.surface.bounce_vel = 0.01;

				contact.geom.depth = -pt.getDistance();
				assert( !isBad( contact.geom.depth ) );
				
				contact.geom.normal[0] = pt.m_normalWorldOnB.x();
				contact.geom.normal[1] = pt.m_normalWorldOnB.y();
				contact.geom.normal[2] = pt.m_normalWorldOnB.z();
				for( int i = 0; i < 3; ++i )
					assert( !isBad(contact.geom.normal[i]) );
				
				//contact.geom.g1 does it really need this?
				contact.geom.g1 = 0;
				contact.geom.g2 = 0;
				contact.geom.pos[0] = pt.getPositionWorldOnB().x();
				contact.geom.pos[1] = pt.getPositionWorldOnB().y();
				contact.geom.pos[2] = pt.getPositionWorldOnB().z();
				for( int i = 0; i < 3; ++i )
					assert( !isBad(contact.geom.pos[i]) );

				this->addContactJoint( b1, b2, contact );
			}
		}
	}
}
#endif

#ifdef USE_NEWTON
BOOST_STATIC_ASSERT( sizeof(dFloatMat4) == 4*4*sizeof(dFloat) );
BOOST_STATIC_ASSERT( sizeof(dFloatMat3) == 3*3*sizeof(dFloat) );
BOOST_STATIC_ASSERT( sizeof(dFloatVec4) == 4*sizeof(dFloat) );
BOOST_STATIC_ASSERT( sizeof(dFloatVec3) == 3*sizeof(dFloat) );

void newtonApplyGravity( const NewtonBody* bodyPtr )
{
	dFloat mass;
	dFloat I[3];
	NewtonBodyGetMassMatrix(
		bodyPtr,
		&mass,
		&I[0],
		&I[1],
		&I[2] );
    
	const dFloatVec3 G( 0.0, -9.81, 0.0 );
	dFloatVec3 gravForce = mass*G;
	NewtonBodyAddForce( bodyPtr, gravForce.Ref() );
}

NewtonSimulation::NewtonSimulation(boost::shared_ptr<const Scene> scene)
	: Simulation( scene )
{
	this->world_ = NewtonCreate( NULL, NULL );

	std::deque< ConstSceneGraphElementPtr > queue( 1, scene->root() );
	while( !queue.empty() )
	{
		ConstSceneGraphElementPtr current = queue.front();
		queue.pop_front();

		boost::shared_ptr<const CombinedObject> fused = 
			boost::dynamic_pointer_cast<const CombinedObject>(current);
		if( !fused )
		{
			std::copy( current->children_begin(), current->children_end(),
				std::front_inserter( queue ) );
		}

		boost::shared_ptr< const PhysicsObject > object = 
			boost::dynamic_pointer_cast<const PhysicsObject>( current );
		if( !object )
			continue;

		// Objects where the object is displaced w.r.t. its center of mass
		//   will need to have their collision geometry transformed appropriately.
		//   Similarly, objects that are static will need their positions
		//   translated w.r.t the world.
		dFloatMat4 geomTransform( vl::vl_1 );

		dFloat I[3];
		dFloat m;

		if( !object->isStatic() )
		{
			dMass mass = object->getOdeMassProps();
			m = mass.mass;

			std::pair< vl::Vec3d, vl::Mat4d > pr = 
				diagonalizeInertiaTensor( mass );

			std::copy( pr.first.Ref(), pr.first.Ref() + 3, I );
			for( vl::Int i = 0; i < 4; ++i )
				for( vl::Int j = 0; j < 4; ++j )
					geomTransform[i][j] = pr.second[i][j];
		}

		NewtonCollision* coll = object->newtonCollision( this->world_, 
			geomTransform, 0.0f );
		if( coll == 0 )
			continue;

		NewtonBody* bodyId = NewtonCreateBody( world_, coll );
		//NewtonReleaseCollision( world_, coll );

		if( !object->isStatic() )
		{
			dynamicObjects_.push_back( object );
			dynamicBodies_.push_back( bodyId );
			offsetTransformations_.push_back( geomTransform );

			NewtonBodySetForceAndTorqueCallback( bodyId, &newtonApplyGravity );
			NewtonBodySetMassMatrix( bodyId, m, I[0], I[1], I[2] );

			vl::Vec3f linearVel = object->linearVelocity();
			NewtonBodySetVelocity ( bodyId, linearVel.Ref() );
			vl::Vec3f angularVel = object->angularVelocity();
			NewtonBodySetOmega( bodyId, angularVel.Ref() );
		}

		dFloatMat4 newtonTransform = vl::trans(object->rigidTransform() * geomTransform);
		NewtonBodySetMatrix( bodyId, newtonTransform.Ref() );
	}

#ifdef _DEBUG
	this->setSimulationState( this->getSimulationState() );
#endif
}

NewtonSimulation::~NewtonSimulation()
{
	NewtonDestroy( this->world_ );
}

Simulation::State NewtonSimulation::getSimulationState() const
{
	std::vector<RigidDynamicState> states;
	states.resize( this->dynamicBodies_.size() );
	for( size_t iBody = 0; iBody < dynamicBodies_.size(); ++iBody )
	{
		dFloatMat4 mat;

		NewtonBodyGetMatrix( dynamicBodies_[iBody], mat.Ref() );
#ifdef _DEBUG
		assert( withinEpsilon( det(trans(mat)), 1.0 ) );
#endif
		mat = trans( mat ) * rigidInverse( this->offsetTransformations_.at(iBody) );
#ifdef _DEBUG
		assert( withinEpsilon( det(trans(mat)), 1.0 ) );
#endif

		states[iBody].setPosition( NxVec3( mat[0][3], mat[1][3], mat[2][3] ) );

		Quaternion<dReal> q( mat );
		states[iBody].setOrientation( q.toNxQuat() );

		dFloatVec3 velocity;
		NewtonBodyGetVelocity( dynamicBodies_[iBody], velocity.Ref() );
		states[iBody].setLinearVelocity( toNxVec(velocity) );

		dFloatVec3 omega;
		NewtonBodyGetOmega( dynamicBodies_[iBody], omega.Ref() );
		states[iBody].setAngularVelocity( toNxVec(omega) );
	}

	return Simulation::State( states, this->time() );
}

void NewtonSimulation::setSimulationState( const Simulation::State& state )
{
	assert( state.stateCount() == this->dynamicBodies_.size() );
	for( size_t iBody = 0; iBody < this->dynamicBodies_.size(); ++iBody )
	{
		const RigidDynamicState& bodyState = state.state(iBody);
		NewtonBody* body = this->dynamicBodies_.at(iBody);

		NxMat33 rotMat( bodyState.orientation() );
		NxVec3 pos = bodyState.position();

		dFloatMat4 newtonState( vl::vl_1 );
		rotMat.getRowMajorStride4( newtonState.Ref() );

		for( vl::Int i = 0; i < 3; ++i )
			newtonState[i][3] = pos[i];

		newtonState = trans(newtonState * this->offsetTransformations_.at(iBody));

		NewtonBodySetMatrix( body, newtonState.Ref() );

		NxVec3 vel = bodyState.linearVelocity();
		dFloatVec3 v( vel[0], vel[1], vel[2] );
		NewtonBodySetVelocity( body, v.Ref() );

		NxVec3 omega = bodyState.angularVelocity();
		dFloatVec3 o( omega[0], omega[1], omega[2] );
		NewtonBodySetVelocity( body, o.Ref() );
	}

	setTime( state.time() );
}

void NewtonSimulation::step( const float stepSize )
{
	NewtonUpdate( this->world_, stepSize );
}

void NewtonSimulation::setDisabled( size_t iDynamicObject, bool disabled )
{
	if( disabled )
		NewtonWorldFreezeBody( this->world_, dynamicBodies_.at(iDynamicObject) );
	else
		NewtonWorldUnfreezeBody( this->world_, dynamicBodies_.at(iDynamicObject) );
}

bool NewtonSimulation::getDisabled( size_t iDynamicObject ) const
{
	// @todo not sure if this is doing the right thing
	int sleeping = NewtonBodyGetSleepingState( dynamicBodies_.at(iDynamicObject) );
	return (sleeping == 0);
}
#endif

#ifdef USE_NOVODEX
// For some reason, occasionally when I try to create a new Novodex scene it 
//   fails silently, which naturally causes major problems for me.  To help
//   alleviate this, we'll maintain a pool of scenes that we can draw on
class NxScenePool
{
public:
	~NxScenePool() {}

	static NxScenePool& instance()
	{
		{
			boost::mutex::scoped_lock lock( singletonLock_ );
			if( !instance_ )
				instance_.reset( new NxScenePool );
		}

		return *instance_;
	}

	NxScenePtr get(NxPhysicsSDKPtr sdk)
	{
		{
			// We believe that it may not be threadsafe to create
			//   more than one scene at the same time, so we'll keep this lock
			//  the entire time:
			boost::mutex::scoped_lock lock( this->poolLock_ );
			if( !pool_.empty() )
			{
				NxScenePtr result = pool_.front().second;
				pool_.pop_front();
				return result;
			}

			assert( sdk );
			NxSceneDesc sceneDesc;
			sceneDesc.gravity.set(0, -9.81f, 0);
			sceneDesc.groundPlane = false;

			boost::recursive_mutex::scoped_lock createSceneLock( novodexCreateDeleteMutex );
			// I'm not sure why but sometimes createScene returns 0
			NxScenePtr result;
			while( !result )
			{
				result.reset(
					sdk->createScene(sceneDesc),
					NxSceneDeleter(sdk) );
				assert( result );
			}
			return result;
		}
	}

	void put(NxScenePtr scene, NxPhysicsSDKPtr sdk)
	{
		boost::mutex::scoped_lock lock( this->poolLock_ );
		pool_.push_back( ScenePair(sdk, scene) );
	}

protected:
	NxScenePool() {}

private:
	static boost::scoped_ptr<NxScenePool> instance_;
	static boost::mutex singletonLock_;

	boost::mutex poolLock_;

	typedef std::pair<NxPhysicsSDKPtr, NxScenePtr> ScenePair;
	std::deque<ScenePair> pool_;
};

boost::scoped_ptr<NxScenePool> NxScenePool::instance_;
boost::mutex NxScenePool::singletonLock_;


NxSimulation::NxSimulation(NxPhysicsSDKPtr sdk, boost::shared_ptr<const Scene> scene)
	: Simulation( scene ), physicsSDK_(sdk)
{
	// Let's just totally block since Novodex seems to have so much trouble
	// with multiple threads creating scenes at the same time
	static boost::mutex createSimulationMutex;
	boost::mutex::scoped_lock createSceneLock( createSimulationMutex );

	// set up Novodex scene
	this->nxScene_ = NxScenePool::instance().get(sdk);
	this->nxScene_->setTiming(
		1.0f/480.0f,
		64,
		NX_TIMESTEP_FIXED );
	this->nxScene_->setUserTriggerReport(this);
	this->nxScene_->setActorGroupPairFlags(0, 0, NX_NOTIFY_ON_START_TOUCH);

	this->addTree( scene->root(), true );

	// must do joints after other objects
	for( Scene::object_iter iter = scene->begin_objects(); 
		iter != scene->end_objects(); ++iter )
	{
		boost::shared_ptr< const Joint > joint = 
			boost::dynamic_pointer_cast<const Joint>( *iter );
		if( !joint )
			continue;

		std::vector< boost::shared_ptr<NxJointDesc> > jointDescs = 
			joint->getNxJointDesc(*this);
		for( std::vector< boost::shared_ptr<NxJointDesc> >::const_iterator jointIter = jointDescs.begin();
			jointIter != jointDescs.end(); ++jointIter )
		{
			boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
			NxJointPtr joint( nxScene_->createJoint( **jointIter ), 
				NxJointDeleter(this->nxScene_) );
			assert( joint );
			joints_.push_back( joint );
		}
	}

	// the first call to simulate does nothing, we call it once to make
	//   sure we aren't going to screw everything else up.
	Simulation::State state = this->getSimulationState();
	this->step( 1.0/240.0 );
	this->setSimulationState( state );
	assert( !asleep() );
}

void NxSimulation::addTree( ConstSceneGraphElementPtr root, bool selfCollisions, std::deque<NxActorPtr>& actors )
{
	if( selfCollisions )
		assert( actors.empty() );

	bool selfCollisionsLocal = selfCollisions && root->selfCollisions();

	ConstPhysicsObjectPtr physicsObject = boost::dynamic_pointer_cast<const PhysicsObject>( root );
	if( physicsObject )
	{
		int id = this->addObject( physicsObject, false );

		// if self collisions were disabled somewhere in the parent, we need
		//   to disable collisions with all the others that share the same parent
		if( !selfCollisionsLocal )
		{
			NxActorPtr actor;
			if( id >= 0 )
				actor = this->dynamicActors_.back();
			else
				actor = this->staticActors_.back();

			for( std::deque<NxActorPtr>::const_iterator actorItr = actors.begin();
				actorItr != actors.end(); ++actorItr )
			{
				this->nxScene_->setActorPairFlags( *(*actorItr), *actor, NX_IGNORE_PAIR );
			}

			actors.push_back( actor );
		}
	}

	boost::shared_ptr<const CombinedObject> fused = 
		boost::dynamic_pointer_cast<const CombinedObject>(root);
	if( !fused )
	{
		for( SceneGraphElement::const_child_iterator childItr = root->children_begin();
			childItr != root->children_end(); ++childItr )
		{
			this->addTree( *childItr, selfCollisionsLocal, actors );
		}
	}

	if( selfCollisions && !selfCollisionsLocal )
		actors.clear();
}


int NxSimulation::addObject( ConstPhysicsObjectPtr object, bool ghost )
{
	float padding = 0.0f;
	if( ghost )
		padding = 1e-2;

	int objectId = -1;

	NxMat33 nxRotation;
	NxVec3 position;

	{
		vl::Mat4f transform = object->transform();
		vl::Mat3f rotate = polarDecomposition( toMat3f(transform) ).first;
		for( vl::Int i = 0; i < 3; ++i )
			for( vl::Int j = 0; j < 3; ++j )
				nxRotation(i, j) = rotate[i][j];

		position = toNxVec( vl::xform( transform, vl::Vec3f(vl::vl_0) ) );
	}

	NxActorDesc actorDesc;

	// initialize actor description
	nameCache_.push_back( new std::string( object->name() ) );
	actorDesc.name = nameCache_.back().c_str();
	actorDesc.globalPose = NxMat34(nxRotation, position);
	actorDesc.density = 0.0f;

	NxBodyDesc desc;
	if( !object->isStatic() )
	{
		desc.angularVelocity = toNxVec(object->angularVelocity());
		desc.linearVelocity = toNxVec(object->linearVelocity());

		desc.solverIterationCount = 10;

		dMass mass = object->getOdeMassProps();
		std::pair< vl::Vec3d, vl::Mat4d > pr = 
			diagonalizeInertiaTensor( mass );
		NxMat33 rot;
		for( size_t i = 0; i < 3; ++i )
			for( size_t j = 0; j < 3; ++j )
				rot(i, j) = pr.second[i][j];

		NxVec3 trans;
		for( size_t i = 0; i < 3; ++i )
			trans[i] = pr.second[i][3];
		desc.massLocalPose = NxMat34(rot, trans);

		desc.massSpaceInertia = toNxVec( pr.first );
		desc.mass = mass.mass;

		actorDesc.body = &desc;

		// objects default to active
		this->active_.push_back( true );
	}

	ConstMaterialPtr mat = object->material();
	NxMaterialMap::const_iterator matItr = this->nxMaterials_.find( mat );
	if( matItr == nxMaterials_.end() )
	{
		NxMaterialDesc nxMat;
		nxMat.dynamicFriction = mat->dynamicFriction();
		nxMat.staticFriction = mat->staticFriction();
		nxMat.restitution = mat->restitution();

		boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
		NxMaterialPtr newMat( nxScene_->createMaterial(nxMat), 
			NxMaterialDeleter(nxScene_) );
		bool inserted;
		boost::tie(matItr, inserted) = nxMaterials_.insert( std::make_pair(mat, newMat) );
		assert( inserted );
	}

	assert( matItr != nxMaterials_.end() );
	NxMaterialPtr nxMat = matItr->second;

	std::vector< boost::shared_ptr<NxShapeDesc> > shapes = 
		object->getNxShapeDesc(this->physicsSDK_, padding);
	for( std::vector< boost::shared_ptr<NxShapeDesc> >::const_iterator shapeIter = shapes.begin();
		shapeIter != shapes.end(); ++shapeIter )
	{
		(*shapeIter)->materialIndex = nxMat->getMaterialIndex();
		if( ghost )
			(*shapeIter)->shapeFlags |= NX_TRIGGER_ENABLE;

		actorDesc.shapes.pushBack( shapeIter->get() );
		if( !actorDesc.isValid() )
			throw CreationException( "Unable to create object; invalid description." );
	}

	// apparently not threadsafe?
	NxActorPtr actor;
	{
		boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );

		assert( actorDesc.isValid() );
		actor.reset( nxScene_->createActor( actorDesc ), 
			NxActorDeleter(this->nxScene_) );
		this->physicsObjectToActorMap_.insert(
			PhysicsObjectToNxActorMap::value_type( object.get(), actor ) );
	}

	/*
	actor->setSleepLinearVelocity( 0.01f );
	actor->setSleepAngularVelocity( 0.01f );
	*/

	assert( actor );

	if( ghost )
		actor->setGroup(1);
	else
		actor->setGroup(0);

	if( object->isStatic() )
		this->staticActors_.push_back( actor );
	else
	{
		objectId = dynamicActors_.size();
		this->dynamicActors_.push_back( actor );
		this->dynamicObjects_.push_back( object );
	}

	// store some useful information in the object
	GeomInfo* info = this->geomInfoPool_.construct();
	info->iDynamicId = objectId;
	info->originalObject = object.get();
	info->ghost = ghost;

	actor->userData = reinterpret_cast<void*>( info );


	return objectId;
}

NxSimulation::~NxSimulation()
{
	NxScenePool::instance().put( this->nxScene_, this->physicsSDK_ );
}

void NxSimulation::onTrigger(NxShape& triggerShape, NxShape& otherShape, NxTriggerFlag status)
{
	const GeomInfo* triggerInfo = 
		reinterpret_cast<const GeomInfo*>( triggerShape.getActor().userData );
	const GeomInfo* otherInfo = 
		reinterpret_cast<const GeomInfo*>( otherShape.getActor().userData );

	bool active = this->isActive( otherInfo->originalObject );
	if( !active && otherInfo->iDynamicId >= 0 )
	{
		this->active_[ otherInfo->iDynamicId ] = true;
		otherShape.getActor().wakeUp();
	}
}

Simulation::State NxSimulation::getSimulationState() const
{
	std::vector<RigidDynamicState> state;
	state.reserve( dynamicActors_.size() );
	for( size_t i = 0; i < dynamicActors_.size(); ++i )
	{
		state.push_back(
			RigidDynamicState( 
				vl::toVec3f(dynamicActors_[i]->getGlobalPosition()),
				Quaternion<float>( dynamicActors_[i]->getGlobalOrientationQuat() ),
				vl::toVec3f(dynamicActors_[i]->getLinearVelocity()),
				vl::toVec3f(dynamicActors_[i]->getAngularVelocity()) ) );
	}

	return Simulation::State(state, time());
}

void NxSimulation::setSimulationState( size_t iObject, const RigidDynamicState& s )
{
	dynamicActors_[iObject]->setGlobalPosition( toNxVec(s.position()) );
	dynamicActors_[iObject]->setGlobalOrientationQuat( s.orientation().toNxQuat() );
	dynamicActors_[iObject]->setLinearVelocity( toNxVec(s.linearVelocity()) );
	dynamicActors_[iObject]->setAngularVelocity( toNxVec(s.angularVelocity()) );
}

void NxSimulation::setSimulationState( const Simulation::State& state )
{
	assert( state.stateCount() == this->dynamicActors_.size() );

	setTime( state.time() );
	for( size_t i = 0; i < dynamicActors_.size(); ++i )
	{
		this->setSimulationState(i, state.state(i));
	}
}

void NxSimulation::perturbLinearVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation )
{
	NxVec3 vel = dynamicActors_.at(iDynamicObject)->getLinearVelocity();
	dynamicActors_.at(iDynamicObject)->setLinearVelocity( vel + toNxVec(perturbation) );
}

void NxSimulation::perturbAngularVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation )
{
	NxVec3 vel = dynamicActors_.at(iDynamicObject)->getAngularVelocity();
	dynamicActors_.at(iDynamicObject)->setAngularVelocity( vel + toNxVec(perturbation) );
}

vl::Vec3f NxSimulation::getAngularVelocity( size_t iDynamicObject ) const
{
	return vl::toVec3f(dynamicActors_.at(iDynamicObject)->getAngularVelocity());
}

vl::Vec3f NxSimulation::getLinearVelocity( size_t iDynamicObject ) const
{
	return vl::toVec3f(dynamicActors_.at(iDynamicObject)->getLinearVelocity());
}


bool NxSimulation::getDisabled( unsigned int iObject ) const
{
	return !this->active_.at(iObject);
}

void NxSimulation::setDisabled( unsigned int iObject, bool disabled )
{
	if( disabled )
	{
		this->active_.at(iObject) = false;
	}
	else
	{
		this->active_.at(iObject) = true;
		this->dynamicActors_.at(iObject)->wakeUp();
	}
}

void NxSimulation::multiStep(const float endTime, const float stepSize)
{
	float timestep = endTime - time();
	this->step( timestep );
}


typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> InteractionGraph;

class ContactReport : public NxUserContactReport
{
public:
	ContactReport( InteractionGraph& g, NxSimulation* simulation )
		: G_(g), nxSimulation_(simulation) {}

	void onContactNotify(NxContactPair& pair, NxU32 events)
	{
		const NxSimulation::GeomInfo* info1 = 
			reinterpret_cast<const NxSimulation::GeomInfo*>( pair.actors[0]->userData );
		const NxSimulation::GeomInfo* info2 = 
			reinterpret_cast<const NxSimulation::GeomInfo*>( pair.actors[1]->userData );

		bool active1 = nxSimulation_->isActive( info1->originalObject );
		bool active2 = nxSimulation_->isActive( info2->originalObject );

		// ghost should be triggers:
		assert( !info1->ghost && !info2->ghost );

		if( info1->iDynamicId >= 0 && info2->iDynamicId >= 0 )
			add_edge( info1->iDynamicId, info2->iDynamicId, G_ );

		if( active1 && !active2 && info2->iDynamicId >= 0 )
		{
			nxSimulation_->active_[info2->iDynamicId] = true;
			pair.actors[1]->wakeUp();
		}
		else if( active2 && !active1 && info1->iDynamicId >= 0 )
		{
			nxSimulation_->active_[info1->iDynamicId] = true;
			pair.actors[0]->wakeUp();
		}
	}

private:
	InteractionGraph& G_;
	NxSimulation* nxSimulation_;
};

void NxSimulation::step(const float simTime)
{
	using namespace boost;

	// initialize the graph with joints
	InteractionGraph G( dynamicActors_.size() );

	this->nxScene_->resetJointIterator();
	NxJoint* curJoint = this->nxScene_->getNextJoint();
	for( ; curJoint != NULL; curJoint = this->nxScene_->getNextJoint() )
	{
		NxActor* actors[2];
		curJoint->getActors( &actors[0], &actors[1] );
		if( actors[0] == 0 || actors[1] == 0 )
			continue;

		const NxSimulation::GeomInfo* info1 = 
			reinterpret_cast<const NxSimulation::GeomInfo*>( actors[0]->userData );
		const NxSimulation::GeomInfo* info2 = 
			reinterpret_cast<const NxSimulation::GeomInfo*>( actors[1]->userData );

		if( info1->iDynamicId >= 0 && info2->iDynamicId >= 0 )
			add_edge( info1->iDynamicId, info2->iDynamicId, G );
	}

	ContactReport userContactReport( G, this );
	this->nxScene_->setUserContactReport(&userContactReport);
	this->nxScene_->simulate(simTime);
	this->nxScene_->flushStream();
	this->nxScene_->fetchResults(NX_RIGID_BODY_FINISHED, true);
	this->nxScene_->setUserContactReport(0);

	setTime( time() + simTime );

	// need to update the active set by propagating through the graph
	if( !dynamicActors_.empty() )
	{
		// start out by setting anything that the simulator thinks has
		//   fallen asleep to not active
		for( size_t i = 0; i < dynamicActors_.size(); ++i )
		{
			if( dynamicActors_[i]->isSleeping() )
				active_[i] = false;
		}

		std::vector<int> component(num_vertices(G));
		int num = connected_components(G, &component[0]);

		// first we decide whether a particular component is active:
		std::vector<bool> activeComponent( num, false );
		for( size_t i = 0; i < component.size(); ++i )
			activeComponent[ component[i] ] = 
				activeComponent[ component[i] ] || this->active_[i];

		// now we set everything to active exactly when its component is active:
		for( size_t i = 0; i < component.size(); ++i )
			this->active_[i] = activeComponent[ component[i] ];
	}
}

NxActorPtr NxSimulation::nxActor( ConstPhysicsObjectPtr object ) const
{
	PhysicsObjectToNxActorMap::const_iterator itr = 
		physicsObjectToActorMap_.find( object.get() );
	assert( itr != physicsObjectToActorMap_.end() );
	return itr->second;
}

NovodexSimulationFactory::NovodexSimulationFactory(NxPhysicsSDKPtr sdk)
	: sdk_(sdk)
{
}

NovodexSimulationFactory::~NovodexSimulationFactory()
{
}

SimulationPtr NovodexSimulationFactory::construct( boost::shared_ptr<const Scene> scene ) const
{
	SimulationPtr result( new NxSimulation(sdk_, scene) );
	return result;
}
#endif // USE_NOVODEX

#ifdef USE_BULLET
BulletSimulation::BulletSimulation( boost::shared_ptr<const Scene> scene )
	: Simulation( scene )
{
	collisionDispatcher_.reset( new btCollisionDispatcher(&collisionConfig_) );

	// todo figure out what these should be
	const size_t maxProxies = 8192;
	btVector3 worldAabbMin(-10000,-10000,-10000);
	btVector3 worldAabbMax(10000,10000,10000);

	broadPhase_.reset( new btAxisSweep3(worldAabbMin,worldAabbMax,maxProxies) );

	constraintSolver_.reset( new btSequentialImpulseConstraintSolver );

	dynamicsWorld_.reset( new btDiscreteDynamicsWorld(collisionDispatcher_.get(),
		broadPhase_.get(),
		constraintSolver_.get() ) );
	dynamicsWorld_->setGravity( btVector3(0, -9.81, 0) );
	dynamicsWorld_->setConstraintSolver( constraintSolver_.get() );

	std::deque< ConstSceneGraphElementPtr > queue( 1, scene->root() );
	while( !queue.empty() )
	{
		ConstSceneGraphElementPtr current = queue.front();
		queue.pop_front();

		boost::shared_ptr<const CombinedObject> fused = 
			boost::dynamic_pointer_cast<const CombinedObject>(current);
		if( !fused )
		{
			std::copy( current->children_begin(), current->children_end(),
				std::front_inserter( queue ) );
		}

		boost::shared_ptr< const PhysicsObject > object = 
			boost::dynamic_pointer_cast<const PhysicsObject>( current );
		if( object )
			this->addObject( object, false );
	}

#ifdef _DEBUG
	this->setSimulationState( this->getSimulationState() );
#endif
}

int BulletSimulation::addObject( ConstPhysicsObjectPtr object, bool ghost )
{
	// Objects where the object is displaced w.r.t. its center of mass
	//   will need to have their collision geometry transformed appropriately.
	//   Similarly, objects that are static will need their positions
	//   translated w.r.t the world.
	btMat4 geomTransform( vl::vl_1 );

	btVector3 I( 0, 0, 0 );
	btScalar m = 0;

	if( !object->isStatic() )
	{
		dMass mass = object->getOdeMassProps();
		m = mass.mass;

		std::pair< vl::Vec3d, vl::Mat4d > pr = 
			diagonalizeInertiaTensor( mass );

		for( vl::Int i = 0; i < 3; ++i )
			I[i] = pr.first[i];

		for( vl::Int i = 0; i < 4; ++i )
			for( vl::Int j = 0; j < 4; ++j )
				geomTransform[i][j] = pr.second[i][j];
	}

	// need to convert into bullet's transformation format
    vl::Mat4f rigidTransform = object->rigidTransform();
    btMat4 btRigidTransform;
    for( vl::Int i = 0; i < 4; ++i )
        for( vl::Int j = 0; j < 4; ++j )
            btRigidTransform[i][j] = rigidTransform[i][j];
	btMat4 localTransform = btRigidTransform * geomTransform;

	btCollisionShape* shape = object->bulletCollisionShape( toMat4f(geomTransform), 0.0);

	if( shape == 0 )
		return -1;

	collisionShapes_.push_back( shape );

	{
		btCompoundShape* compoundShape = 
			dynamic_cast<btCompoundShape*>( shape );
		if( compoundShape )
		{
			for( int i = 0; i < compoundShape->getNumChildShapes(); ++i )
				collisionShapes_.push_back( compoundShape->getChildShape(i) );
		}
	}

	btMotionState* motionState = new btDefaultMotionState( toBtTransform(localTransform) );
	motionStates_.push_back( motionState );

	btRigidBody* rigidBody = new btRigidBody(m, motionState, shape, I);

	int iDynamicId = -1;
	if( !object->isStatic() )
	{
		iDynamicId = this->addDynamicObject( object );
		dynamicRigidBodies_.push_back( rigidBody );
		offsetTransformations_.push_back( geomTransform );

		vl::Vec3f linearVel = object->linearVelocity();
		rigidBody->setLinearVelocity( toBtVector3(linearVel) );
		vl::Vec3f angularVel = object->angularVelocity();
		rigidBody->setAngularVelocity( toBtVector3(angularVel) );
	}
	else
	{
		staticRigidBodies_.push_back( rigidBody );
	}

	dynamicsWorld_->addRigidBody( rigidBody );

	return iDynamicId;
}

void BulletSimulation::removeObject( int dynamicId )
{
    assert( dynamicId >= 0 );

    // For now, we're not sure what will happen if we remove the body itself
    // so we'll just remove all the geoms associated with it
    dynamicsWorld_->removeRigidBody( &this->dynamicRigidBodies_.at( dynamicId ) );
}

BulletSimulation::~BulletSimulation()
{
}

RigidDynamicState BulletSimulation::getSimulationState(size_t iBody) const
{
	const btRigidBody& body = this->dynamicRigidBodies_.at( iBody );
	btTransform transform = body.getCenterOfMassTransform();

	btMat4 mat;
	transform.getOpenGLMatrix( mat.Ref() );
#ifdef _DEBUG
	assert( withinEpsilon( det(trans(mat)), 1.0 ) );
#endif
	mat = trans(mat) * rigidInverse( this->offsetTransformations_.at(iBody) );
#ifdef _DEBUG
	assert( withinEpsilon( det(trans(mat)), 1.0 ) );
#endif


	Quaternion<float> q( mat );

    return RigidDynamicState( vl::Vec3f( mat[0][3], mat[1][3], mat[2][3] ),
        q,
        vl::toVec3f( body.getLinearVelocity() ),
        vl::toVec3f( body.getAngularVelocity() ) );
}

Simulation::State BulletSimulation::getSimulationState() const
{
	std::vector<RigidDynamicState> states;
	states.reserve( this->dynamicRigidBodies_.size() );
	for( size_t iBody = 0; iBody < dynamicRigidBodies_.size(); ++iBody )
        states.push_back( this->getSimulationState(iBody) );

	return Simulation::State( states, this->time() );
}

void BulletSimulation::setSimulationState( size_t iBody, const RigidDynamicState& bodyState )
{
	btRigidBody& body = this->dynamicRigidBodies_.at(iBody);

	const vl::Mat3f rotMat = bodyState.orientation().toRotMatf();
	const vl::Vec3f pos = bodyState.position();

	btMat4 bulletState( vl::vl_1 );
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			bulletState[i][j] = rotMat[i][j];

	for( vl::Int i = 0; i < 3; ++i )
		bulletState[i][3] = pos[i];

	bulletState = trans(bulletState * this->offsetTransformations_.at(iBody));
	btTransform transform;
	transform.setFromOpenGLMatrix( bulletState.Ref() );
	body.setCenterOfMassTransform( transform );

	body.setLinearVelocity( toBtVector3( bodyState.linearVelocity() ) );
	body.setAngularVelocity( toBtVector3( bodyState.angularVelocity() ) );
}

void BulletSimulation::setSimulationState( const Simulation::State& state )
{
	assert( state.stateCount() == this->dynamicRigidBodies_.size() );
	for( size_t iBody = 0; iBody < this->dynamicRigidBodies_.size(); ++iBody )
	{
		this->setSimulationState( iBody, state.state(iBody) );
	}

	setTime( state.time() );
}

void BulletSimulation::step( const float stepSize )
{
	this->dynamicsWorld_->stepSimulation( 
		stepSize,
		1,
		stepSize );
	setTime( time() + stepSize );
}

void BulletSimulation::multiStep( const float endTime, const float stepSize )
{
	float timeDiff = endTime - time();
	this->dynamicsWorld_->stepSimulation( 
		timeDiff
		/*
		,
		boost::numeric::bounds<int>::highest(),
		stepSize
		*/
		);
	setTime( endTime );
}

void BulletSimulation::setDisabled( size_t iDynamicObject, bool disabled )
{
	// @todo fix this
}

bool BulletSimulation::getDisabled( size_t iDynamicObject ) const
{
	// @todo fix this
	return false;
}

void BulletSimulation::perturbLinearVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation )
{
	btRigidBody& body = this->dynamicRigidBodies_.at( iDynamicObject );
	btVector3 v( body.getLinearVelocity() );
	v += toBtVector3( perturbation );
	body.setLinearVelocity( v );
}

void BulletSimulation::perturbAngularVelocity( size_t iDynamicObject, const vl::Vec3f& perturbation )
{
	btRigidBody& body = this->dynamicRigidBodies_.at( iDynamicObject );
	btVector3 v( body.getAngularVelocity() );
	v += toBtVector3( perturbation );
	body.setAngularVelocity( v );
}

vl::Vec3f BulletSimulation::getLinearVelocity( size_t iDynamicObject ) const
{
	const btRigidBody& body = this->dynamicRigidBodies_.at( iDynamicObject );
	return vl::toVec3f( body.getLinearVelocity() );
}

vl::Vec3f BulletSimulation::getAngularVelocity( size_t iDynamicObject ) const
{
	const btRigidBody& body = this->dynamicRigidBodies_.at( iDynamicObject );
	return vl::toVec3f( body.getAngularVelocity() );
}

#endif // USE_BULLET

PenaltyForceSimulation::PenaltyForceSimulation(boost::shared_ptr<const Scene> scene)
    : ODENativeSimulation( scene )
{

}

PenaltyForceSimulation::~PenaltyForceSimulation()
{
}

typedef std::pair< PenaltyForceSimulation*, bool > PenaltyBackwardsPr;

void PenaltyForceSimulation::handleCollisions(bool backwards)
{
    // need to pass the "backwards" bit back
    PenaltyBackwardsPr pr( this, backwards );

	// perform collision detection
	dSpaceCollide( this->spaceId_, &pr, potentialCollisionPenalty );
}

void potentialCollisionPenalty( void* data, dGeomID o1, dGeomID o2 )
{
    PenaltyBackwardsPr* pr = reinterpret_cast<PenaltyBackwardsPr*>( data );
	PenaltyForceSimulation* simulation = pr->first;
	assert( simulation != 0 );

    bool backwards = pr->second;

	simulation->checkContact( o1, o2, backwards );
}

void PenaltyForceSimulation::checkContact( dGeomID o1, dGeomID o2, bool backwards )
{
	bool isSpace1 = dGeomIsSpace(o1);
	bool isSpace2 = dGeomIsSpace(o2);

	const GeomInfo *info1, *info2;
	if( isSpace1 )
		info1 = reinterpret_cast<const GeomInfo*>( 
			dGeomGetData( dSpaceGetGeom(reinterpret_cast<dSpaceID>(o1), 0) ) );
	else
		info1 = reinterpret_cast<const GeomInfo*>( dGeomGetData(o1) );

	if( isSpace2 )
		info2 = reinterpret_cast<const GeomInfo*>( 
			dGeomGetData( dSpaceGetGeom(reinterpret_cast<dSpaceID>(o2), 0) ) );
	else
		info2 = reinterpret_cast<const GeomInfo*>( dGeomGetData(o2) );

	if( this->ignoreCollisions( info1->originalObject, info2->originalObject ) )
		return;

	bool active1 = isActive( info1->originalObject );
	bool active2 = isActive( info2->originalObject );

#ifdef _DEBUG
	// If the objects are different, they'd better both be ACTIVE
	// If they are in correspondence, then they are supposed to get added to the ACTIVE
	//   set at the same time
	if( info1->originalObject->name() == info2->originalObject->name() )
		assert( active1 == active2 );
#endif

	// If they're both either active or inactive, we can skip colliding them if
	//   one of them is a ghost
	if( (active1 == active2) && (info1->ghost || info2->ghost) )
		return;

	// No point colliding ghost objects against each other
	if( info1->ghost && info2->ghost )
		return;

	if( isSpace1 || isSpace2 )
	{
		dSpaceCollide2( o1, o2, this, potentialCollision );
	}
	else
	{
		dBodyID b1 = dGeomGetBody(o1);
		dBodyID b2 = dGeomGetBody(o2);

		const Material& mat1 = info1->material;
		const Material& mat2 = info2->material;

		const int maxContacts = 128;

		dContactGeom contacts[ maxContacts ];

		// not worth getting multiple contacts for ghost objects, since
		//   1 is enough to activate them.
		int maxContactsTemp = maxContacts;
		if( info1->ghost || info2->ghost )
			maxContactsTemp = 1;

		int numc = dCollide( o1, o2, maxContactsTemp, &contacts[0], sizeof(dContactGeom) );
		assert( numc >= 0 );
		if( numc == 0 )
			return;

		// non-ghost objects don't need to activate ghost objects
		if( active1 && !active2 && !info2->ghost && (b2 != 0) )
		{
			dBodyEnable(b2);
		}
		else if( active2 && !active1 && !info1->ghost && (b1 != 0) )
		{
			dBodyEnable(b1);
		}

		// won't collide ghost objects
		if( info1->ghost || info2->ghost )
		{
			return;
		}

		for( int i = 0; i < numc; ++i )
		{
            const dContactGeom& cg = contacts[i];

            double c = 2000.0;
            double k = 1000000.0;

            dVector3 v1 = { 0.0, 0.0, 0.0, 0.0 };
            if( b1 )
                dBodyGetPointVel( b1, cg.pos[0], cg.pos[1], cg.pos[2], v1 );

            dVector3 v2 = { 0.0, 0.0, 0.0, 0.0 };
            if( b2 )
                dBodyGetPointVel( b2, cg.pos[0], cg.pos[1], cg.pos[2], v2 );

            double relVel = -vl::dot( vl::toVec3d( v1 ) - vl::toVec3d( v2 ), vl::toVec3d(cg.normal) );

            // apparently the fancy name for this is "Kelvin-Voigt penalty model"
            // but basically it's just a damped spring
            double forceVal = k*cg.depth + c*relVel;

            if( b1 )
                dBodyAddForceAtPos( b1, 
                    forceVal*cg.normal[0], 
                    forceVal*cg.normal[1], 
                    forceVal*cg.normal[2], 
                    cg.pos[0], cg.pos[1], cg.pos[2] );

            if( b2 )
                dBodyAddForceAtPos( b2, 
                    -forceVal*cg.normal[0], 
                    -forceVal*cg.normal[1], 
                    -forceVal*cg.normal[2], 
                    cg.pos[0], cg.pos[1], cg.pos[2] );
		}
	}
}

CameraWrapper::CameraWrapper( const std::string& name, boost::shared_ptr<Camera> camera )
	: Named(name), camera_(camera)
{
	HasKeyable<float>::AttributeList attributes = 
		camera_->keyable();
	for( HasKeyable<float>::AttributeList::iterator attributeItr = attributes.begin();
		attributeItr != attributes.end(); ++attributeItr )
	{
		keyable_.push_back( 
			AttributeWithKeysPtr(new AttributeWithKeys(*attributeItr) ) );
	}
}

CameraWrapper::~CameraWrapper()
{
}

CameraWrapper* CameraWrapper::clone() const
{
	boost::shared_ptr<Camera> newCam;
	boost::shared_ptr<const ProjectiveFixedYCamera> projCam = 
		boost::dynamic_pointer_cast<const ProjectiveFixedYCamera>( this->camera_ );
	boost::shared_ptr<const OrthoCamera> orthoCam = 
		boost::dynamic_pointer_cast<const OrthoCamera>( this->camera_ );

	// this could be solved by adding a clone() method to Camera
	if( projCam )
		newCam.reset( new ProjectiveFixedYCamera(*projCam) );
	else if( orthoCam )
		newCam.reset( new OrthoCamera(*orthoCam) );
	else
		assert( false );

	CameraWrapper* newCamWrapper = new CameraWrapper( this->name(), newCam );

	// now should copy all the camera keys over
	assert( newCamWrapper->keyable_.size() == this->keyable_.size() );
	for( size_t i = 0; i < this->keyable_.size(); ++i )
	{
		newCamWrapper->keyable_[i]->setKeys( 
			this->keyable_[i]->keys() );
	}

	return newCamWrapper;
}

boost::shared_ptr<Camera> CameraWrapper::camera()
{
	return this->camera_;
}

boost::shared_ptr<const Camera> CameraWrapper::camera() const
{
	return this->camera_;
}

MechModelObjectPtr CameraWrapper::toMechModel()
{
	boost::shared_ptr<MechModelDeclaredObject> object(
		new MechModelDeclaredObject( "camera", this->name() ) );

	boost::shared_ptr<ProjectiveFixedYCamera> projCam = 
		boost::dynamic_pointer_cast<ProjectiveFixedYCamera>( this->camera_ );
	boost::shared_ptr<OrthoCamera> orthoCam =
		boost::dynamic_pointer_cast<OrthoCamera>( this->camera_ );
	// one or the other, please
	assert( (orthoCam || projCam) && !(orthoCam && projCam) );

	if( projCam )
		object->add( mechModelValuePair("type", "perspective") );
	else
		object->add( mechModelValuePair("type", "orthographic") );

    object->add( mechModelValuePair("near_z", camera_->zMin()) );
    object->add( mechModelValuePair("far_z", camera_->zMax()) );
	object->add( mechModelValuePair("up", camera_->upVector()) );

	if( projCam )
	{
		object->add( mechModelValuePair("elevation", projCam->getElevation()) );
		object->add( mechModelValuePair("azimuth", projCam->getAzimuth()) );
		object->add( mechModelValuePair("dolly", projCam->getDolly()) );
		object->add( mechModelValuePair("twist", projCam->getTwist()) );
		object->add( mechModelValuePair("look", projCam->getLookAt()) );
		object->add( mechModelValuePair("fov", projCam->FOV()) );
	}
	else
	{
		assert( orthoCam );
		object->add( mechModelValuePair("port", orthoCam->getPort()) );
		object->add( mechModelValuePair("rotation", orthoCam->getAzimRot()) );
		object->add( mechModelValuePair("center", orthoCam->getCenter()) );
	}

	// now save all the keys:
	for( std::deque<AttributeWithKeysPtr>::const_iterator attributeItr = this->keyable_.begin();
		attributeItr != this->keyable_.end(); ++attributeItr )
	{
		const AttributeWithKeys& attribute = *(*attributeItr);
		const KeyableAttribute<float>& keyable = attribute.keyable();

		// don't need to dump to file if no keys
		if( attribute.begin_keys() == attribute.end_keys() )
			continue;

		boost::shared_ptr<MechModelDeclaredObject> keys(
			new MechModelDeclaredObject( "keys", keyable.name() ) );

		for( AttributeWithKeys::keyframe_iterator keyIter = attribute.begin_keys();
			keyIter != attribute.end_keys(); ++keyIter )
		{
			keys->add( mechModelValuePair("key", vl::Vec2f(keyIter->time(), keyIter->value())) );
		}

		object->add( keys );
	}

    return object;
}

float CameraWrapper::nextKey(float time) const
{
	bool found = false;
	float best = boost::numeric::bounds<float>::highest();
	for( std::deque<AttributeWithKeysPtr>::const_iterator attItr = this->keyable_.begin();
		attItr != this->keyable_.end(); ++attItr )
	{
		float next = (*attItr)->nextKey(time);
		if( next > time )
		{
			best = std::min( best, next );
			found = true;
		}
	}

	if( found )
		return best;
	else
		return time;
}

float CameraWrapper::prevKey(float time) const
{
	bool found = false;
	float best = boost::numeric::bounds<float>::lowest();
	for( std::deque<AttributeWithKeysPtr>::const_iterator attItr = this->keyable_.begin();
		attItr != this->keyable_.end(); ++attItr )
	{
		float prev = (*attItr)->prevKey(time);
		if( prev < time )
		{
			best = std::max( best, prev );
			found = true;
		}
	}

	if( found )
		return best;
	else
		return time;
}

void CameraWrapper::deleteKey(float time)
{
	std::for_each( keyable_.begin(), keyable_.end(), 
		boost::bind(&AttributeWithKeys::deleteKey, _1, time) );
}

void CameraWrapper::setKey(float time)
{
	std::for_each( keyable_.begin(), keyable_.end(), 
		boost::bind(&AttributeWithKeys::setKey, _1, time) );
}

void CameraWrapper::setTime(float time)
{
	std::for_each( keyable_.begin(), keyable_.end(), 
		boost::bind(&AttributeWithKeys::setTime, _1, time) );
}

CameraWrapper::attribute_iterator CameraWrapper::keyable_begin()
{
	return keyable_.begin();
}

CameraWrapper::attribute_iterator CameraWrapper::keyable_end()
{
	return keyable_.end();
}

float CameraWrapper::firstKey() const
{
	float best = boost::numeric::bounds<float>::highest();
	for( std::deque<AttributeWithKeysPtr>::const_iterator attItr = this->keyable_.begin();
		attItr != this->keyable_.end(); ++attItr )
	{
		best = std::min( best, (*attItr)->firstKey() );
	}

	return best;
}

float CameraWrapper::lastKey() const
{
	float best = boost::numeric::bounds<float>::lowest();
	for( std::deque<AttributeWithKeysPtr>::const_iterator attItr = this->keyable_.begin();
		attItr != this->keyable_.end(); ++attItr )
	{
		best = std::max( best, (*attItr)->lastKey() );
	}

	return best;
}

SimulationFactory::~SimulationFactory()
{
}

// Helper function for symmetricDiffScene
void addTree( ConstSceneGraphElementPtr root, std::deque<ConstSceneGraphElementPtr>& obs )
{
	std::deque<ConstSceneGraphElementPtr> toHandle( 1, root );
	while( !toHandle.empty() )
	{
		ConstSceneGraphElementPtr current = toHandle.back();
		toHandle.pop_back();
		
		std::copy( current->children_begin(), current->children_end(), 
			std::back_inserter( toHandle ) );
		obs.push_back( current );
	}
}

struct CompareSceneGraphElementNames
{
	bool operator()( ConstSceneGraphElementPtr left, ConstSceneGraphElementPtr right ) const
	{
		std::less<std::string> comparator;
		return comparator( left->name(), right->name() );
	}
};


// return all the objects in changedScene that are different
//   from the original:
std::deque<ConstSceneGraphElementPtr> symmetricDiffScene(
	const Scene& origScene,
	const Scene& changedScene )
{
	typedef std::pair< ConstSceneGraphElementPtr, ConstSceneGraphElementPtr > ElementToHandle;
	std::deque<ElementToHandle> toHandle( 1, 
		ElementToHandle( origScene.root(), changedScene.root() ) );

	std::deque<ConstSceneGraphElementPtr> changed;
	while( !toHandle.empty() )
	{
		ElementToHandle current = toHandle.back();
		toHandle.pop_back();

		if( !current.first->equalsLocal( *current.second ) )
		{
			addTree( current.first, changed );
			addTree( current.second, changed );
		}

		std::vector<ConstSceneGraphElementPtr> origSceneElts(
			current.first->children_begin(), current.first->children_end() );
		std::vector<ConstSceneGraphElementPtr> changedSceneElts(
			current.second->children_begin(), current.second->children_end() );
		std::sort( origSceneElts.begin(), origSceneElts.end(), 
			CompareSceneGraphElementNames() );
		std::sort( changedSceneElts.begin(), changedSceneElts.end(), 
			CompareSceneGraphElementNames() );
		std::vector<ConstSceneGraphElementPtr>::iterator origItr = origSceneElts.begin();
		std::vector<ConstSceneGraphElementPtr>::iterator changedItr = changedSceneElts.begin();

		while( !(origItr == origSceneElts.end() && changedItr == changedSceneElts.end()) )
		{
			if( changedItr == changedSceneElts.end() )
				addTree( *origItr++, changed );
			else if( origItr == origSceneElts.end() )
				addTree( *changedItr++, changed );
			else if( (*origItr)->name() < (*changedItr)->name() )
				addTree( *origItr++, changed );
			else if( (*changedItr)->name() < (*origItr)->name() )
				addTree( *changedItr++, changed );
			else
				toHandle.push_back( ElementToHandle( *origItr++, *changedItr++ ) );
		}
	}

	// Now, if there are any joints in the dirty object set, we need to make sure that
	//   the objects they join get added, too
	std::deque<ConstSceneGraphElementPtr> changedByJoints;
	for( std::deque<ConstSceneGraphElementPtr>::const_iterator dirtyItr = changed.begin();
		dirtyItr != changed.end(); ++dirtyItr )
	{
		boost::shared_ptr<const Joint> joint = 
			boost::dynamic_pointer_cast<const Joint>( *dirtyItr );
		if( !joint )
			continue;

		std::vector<ConstPhysicsObjectPtr> objects = joint->objects();
		for( std::vector<ConstPhysicsObjectPtr>::const_iterator obItr = objects.begin();
			obItr != objects.end(); ++obItr )
		{
			// This is an okay thing to do because if the name of an object changes, it
			//   will get marked as "changed" regardless
			// This may generate dups in the set, fyi
			ConstSceneGraphElementPtr curElt = origScene.object( (*obItr)->name() );
			if( curElt )
				changedByJoints.push_back( curElt );

			ConstSceneGraphElementPtr treeElt = changedScene.object( (*obItr)->name() );
			if( treeElt )
				changedByJoints.push_back( treeElt );
		}
	}

	// To guarantee uniqueness
	std::sort( changed.begin(), changed.end() );
	std::sort( changedByJoints.begin(), changedByJoints.end() );

	std::deque<ConstSceneGraphElementPtr> result;
	std::merge( changed.begin(), changed.end(), changedByJoints.begin(), changedByJoints.end(), 
		std::back_inserter(result) );

	return result;
}

template <typename PhysicsObjectType>
void addPhysicsObjects( ConstSceneGraphElementPtr current, std::vector< boost::shared_ptr<PhysicsObjectType> >& objects )
{
	boost::shared_ptr< PhysicsObjectType > object = 
		boost::dynamic_pointer_cast<PhysicsObjectType>(current);
	if( object && !object->isStatic() )
		objects.push_back( object );

	boost::shared_ptr<const CombinedObject> fused = 
		boost::dynamic_pointer_cast<const CombinedObject>(current);
	if( !fused )
	{
		std::for_each( current->children_begin(), current->children_end(), 
			boost::bind( &addPhysicsObjects<PhysicsObjectType>, _1, boost::ref(objects) ) );
	}
}

std::vector<ConstPhysicsObjectPtr> physicsObjects( const Scene& scene )
{
	std::vector<ConstPhysicsObjectPtr> result;
	addPhysicsObjects( scene.root(), result );
	return result;
}

DummySimulation::DummySimulation( boost::shared_ptr<const Scene> scene, 
    const Simulation::State& state,
    const std::vector<Simulation::ContactPoint>& contactPoints )
    : Simulation( scene ), state_( state ), contactPoints_( contactPoints )
{
    std::vector<ConstPhysicsObjectPtr> objects = physicsObjects( *scene );
    for( std::vector<ConstPhysicsObjectPtr>::const_iterator iter = objects.begin();
        iter != objects.end(); ++iter )
        this->addDynamicObject( *iter );
}

DummySimulation::~DummySimulation()
{
}

vl::Vec3f DummySimulation::getLinearVelocity( size_t iDynamicObject ) const
{
    return this->state_.state( iDynamicObject ).linearVelocity();
}

vl::Vec3f DummySimulation::getAngularVelocity( size_t iDynamicObject ) const
{
    return this->state_.state( iDynamicObject ).angularVelocity();
}


} // namespace planning
