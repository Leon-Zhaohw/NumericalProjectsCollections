#include "stdafx.h"

#include "physicsObject.h"
#include "scene.h"
#include "constraints.h"
#include "mechModel.h"
#include "sphere_tree.h"
#include "ribFile.h"
#include "voxels.h"
#include "sound.h"

#ifdef USE_BULLET
#include "bulletWrappers.h"

#include "LinearMath/btGeometryUtil.h"
#endif

#ifdef CONDOR_SAMPLING
#include "MWRMComm.h"
#endif

#ifdef USE_NOVODEX
#include <NxCooking.h> 
#endif

#include "twigg/vlutil.h"

#include "ConvexDecomposition/ConvexDecomposition.h"
#include "ConvexDecomposition/cd_hull.h"

#ifdef GUI
#include "checkerboard.h"
#include "twigg/renderUtils.h"

#include <wx/image.h>
#include <wx/mstream.h>
#endif

#include "twigg/magicWrappers.h"
#include "twigg/linalg.h"
#include "twigg/boundingbox.h"
#include "twigg/objfile.h"
#include "twigg/EulerAngles.h"

#include "Wm4Box3.h"
#include "Wm4IntrBox3Box3.h"


#include <mkl_cblas.h>

#include <GL/glut.h>

#include <boost/filesystem/operations.hpp>

#include <iomanip>
#include <fstream>

namespace planning {

#include "meshes/meshes.h"

#define STATICMESH_INIT( _prefix ) StaticTriangleMesh  _prefix##_mesh( \
	_prefix##_positions,                                               \
	sizeof(_prefix##_positions) / (3*sizeof(float)),                   \
	_prefix ##_normals,                                                \
	_prefix ##_triangles,                                              \
	sizeof(_prefix##_triangles) / (3*sizeof(unsigned int)),            \
	_prefix##_boundarySegments,                                        \
	sizeof(_prefix##_boundarySegments) / (2*sizeof(unsigned int)) )

STATICMESH_INIT( ballJointLeft );
STATICMESH_INIT( ballJointRight );
STATICMESH_INIT( univJointLeft );
STATICMESH_INIT( univJointRight );
STATICMESH_INIT( hingeJointLeft );
STATICMESH_INIT( hingeJointRight );
STATICMESH_INIT( sliderJointLeft );
STATICMESH_INIT( sliderJointRight );
STATICMESH_INIT( hinge2JointLeft );
STATICMESH_INIT( hinge2JointRight );
STATICMESH_INIT( hinge2JointMid );
STATICMESH_INIT( rotMotorLeft );
STATICMESH_INIT( rotMotorRight );
STATICMESH_INIT( rotMotorMid );
STATICMESH_INIT( planeJoint );


const GLTexture& checkerboardTextureMap()
{
    static boost::scoped_ptr<GLTexture> texture;
    static boost::mutex mapLock;
    boost::mutex::scoped_lock lock( mapLock );

	if( !texture )
	{
		wxMemoryInputStream pngStream( checkerboard_png, sizeof(checkerboard_png) );
		wxImage checkerboardImage( pngStream, wxBITMAP_TYPE_PNG );

        texture.reset( new GLTexture() );
		GLBindTextureHandler bindTexture( *texture );

		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		gluBuild2DMipmaps( GL_TEXTURE_2D, 3, 
			checkerboardImage.GetWidth(), checkerboardImage.GetHeight(), GL_RGB,
			GL_UNSIGNED_BYTE, checkerboardImage.GetData() );
	}

    return *texture;
}

const int sphereSubdivisions = 20;

static int my_checkMass (dMass *m)
{
  int i;

  if (m->mass <= 0) {
    dDEBUGMSG ("mass must be > 0");
    return 0;
  }
  if (!dIsPositiveDefinite (m->I,3)) {
    dDEBUGMSG ("inertia must be positive definite");
    return 0;
  }

  // verify that the center of mass position is consistent with the mass
  // and inertia matrix. this is done by checking that the inertia around
  // the center of mass is also positive definite. from the comment in
  // dMassTranslate(), if the body is translated so that its center of mass
  // is at the point of reference, then the new inertia is:
  //   I + mass*crossmat(c)^2
  // note that requiring this to be positive definite is exactly equivalent
  // to requiring that the spatial inertia matrix
  //   [ mass*eye(3,3)   M*crossmat(c)^T ]
  //   [ M*crossmat(c)   I               ]
  // is positive definite, given that I is PD and mass>0. see the theorem
  // about partitioned PD matrices for proof.

  dMatrix3 I2,chat;
  dSetZero (chat,12);
  dCROSSMAT (chat,m->c,4,+,-);
  dMULTIPLY0_333 (I2,chat,chat);
  for (i=0; i<3; i++) I2[i] = m->I[i] + m->mass*I2[i];
  for (i=4; i<7; i++) I2[i] = m->I[i] + m->mass*I2[i];
  for (i=8; i<11; i++) I2[i] = m->I[i] + m->mass*I2[i];
  if (!dIsPositiveDefinite (I2,3)) {
    dDEBUGMSG ("center of mass inconsistent with mass parameters");
    return 0;
  }
  return 1;
}

typedef TinyVec<unsigned char, 3> Color;
const std::vector<Color>& mapColors()
{
	static std::vector<Color> mapColors;
	if( mapColors.empty() )
	{
		// Color scheme from ColorBrewer
		// http://www.personal.psu.edu/faculty/c/a/cab38/ColorBrewerBeta.html
		mapColors.push_back( Color(166, 206, 227) );
		mapColors.push_back( Color(31, 120, 180) );
		mapColors.push_back( Color(178, 223, 138) );
		mapColors.push_back( Color(51, 160, 44) );
		mapColors.push_back( Color(251, 154, 153) );
		mapColors.push_back( Color(227, 26, 28) );
		mapColors.push_back( Color(253, 191, 111) );
		mapColors.push_back( Color(255, 127, 0) );
		mapColors.push_back( Color(202, 178, 214) );
		mapColors.push_back( Color(106, 61, 154) );
		mapColors.push_back( Color(255, 255, 153) );
	}
	return mapColors;
}

// if we know that the matrix is orthonormal, inverting it becomes much faster
//   and better conditioned:
template <typename MatType>
MatType rigidInverseImp( const MatType& rigidTransform )
{
#ifdef _DEBUG
	{
		// preconditions:
		// check that the matrix is truly orthogonal
		vl::Mat3d upperLeft = toMat3d( rigidTransform );
		double det = vl::det( rigidTransform );
		assert( withinEpsilon( det, 1.0 ) );
		assert( withinEpsilon( rigidTransform[3][3], 1.0 ) );
		for( vl::Int i = 0; i < 3; ++i )
			assert( withinEpsilon( rigidTransform[3][i], 0.0 ) );
	}
#endif

	MatType result;
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			result[i][j] = rigidTransform[j][i];

	for( vl::Int i = 0; i < 3; ++i )
	{
		result[i][3] = 0.0f;
		for( vl::Int j = 0; j < 3; ++j )
			result[i][3] -= rigidTransform[j][i]*rigidTransform[j][3];
	}

	for( vl::Int i = 0; i < 3; ++i )
		result[3][i] = 0.0f;

	result[3][3] = 1.0f;

#ifdef _DEBUG
	// postcondition: result is inverse
	{
		MatType zero = result*rigidTransform - MatType(vl::vl_1);
		for( vl::Int i = 0; i < 4; ++i )
			for( vl::Int j = 0; j < 4; ++j )
				assert( withinEpsilon( zero[i][j], 0.0 ) );
	}
#endif

	return result;
}

vl::Mat4f rigidInverse( const vl::Mat4f& rigidTransform )
{
	return rigidInverseImp( rigidTransform );
}

vl::Mat4d rigidInverse( const vl::Mat4d& rigidTransform )
{
	return rigidInverseImp( rigidTransform );
}

float epsilon = 1e-3;
bool withinEpsilon( float value1, float value2 )
{
	// The first conditional is to check for stuff like value1 == value2 == -Inf
	return (value1 == value2) || (fabs( value1 - value2 ) < epsilon);
}

bool withinEpsilon( const vl::Vec3f& value1, const vl::Vec3f& value2 )
{
	return withinEpsilon( value1[0], value2[0] ) &&
		withinEpsilon( value1[1], value2[1] ) && 
		withinEpsilon( value1[2], value2[2] );
}

const int order = EulOrdXYZr;

Quaternion<float> eulerToQuat( const vl::Vec3f& angles )
{
	EulerAngles eulAngles( Eul_(angles[0], angles[1], angles[2], order) );
	Quat shoemake_quat = Eul_ToQuat( eulAngles );
	Quaternion<float> quat( 
		vl::Vec3f(shoemake_quat.x, shoemake_quat.y, shoemake_quat.z),
		shoemake_quat.w );
	assert( !isBad(quat) );
	return quat;
}

vl::Vec3f quatToEuler( const Quaternion<float>& quat )
{
	Quat shoemake_quat;
	shoemake_quat.x = quat.x();
	shoemake_quat.y = quat.y();
	shoemake_quat.z = quat.z();
	shoemake_quat.w = quat.w();
	EulerAngles angles = Eul_FromQuat( shoemake_quat, order );

	return vl::Vec3f( angles.x, angles.y, angles.z );
}

vl::Vec3d toMayaRotation( const Quaternion<double>& quat )
{
	Quat q;
	q.x = quat.x();
	q.y = quat.y();
	q.z = quat.z();
	q.w = quat.w();

	EulerAngles angles = Eul_FromQuat(q, EulOrdXYZs);
	vl::Vec3d ang( toDegrees(angles.x), toDegrees(angles.y), toDegrees(angles.z) );
	return ang;
}

vl::Vec3f matToEuler( const vl::Mat3f& mat )
{
	HMatrix matrix;
	matrix[3][0] = matrix[3][1] = matrix[3][2] = 0.0f;
	matrix[0][3] = matrix[1][3] = matrix[2][3] = 0.0f;
	matrix[3][3] = 1.0f;
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			matrix[i][j] = mat[i][j];
	EulerAngles angles = Eul_FromHMatrix( matrix, order );

	return vl::Vec3f( angles.x, angles.y, angles.z );
}

vl::Mat3f eulerToMat( const vl::Vec3f& angles )
{
	EulerAngles eulAngles;
	eulAngles.x = angles[0];
	eulAngles.y = angles[1];
	eulAngles.z = angles[2];

	HMatrix matrix;
	Eul_ToHMatrix( eulAngles, matrix );
	vl::Mat3f result;
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			result[i][j] = matrix[i][j];

	return result;
}


SceneGraphElement::SceneGraphElement( const std::string& name )
	: Named(name), visible_(true), selfCollisions_( true )
{
}

SceneGraphElement::SceneGraphElement( const SceneGraphElement& other )
	: Named(other.name()), 
		visible_(other.visible_), 
		selfCollisions_(other.selfCollisions_), 
		randomizedProperties_(other.randomizedProperties_)
{
}

SceneGraphElement::~SceneGraphElement()
{
}

SceneGraphElement::child_iterator SceneGraphElement::children_begin()
{
	return child_iterator( this->firstChild_ );
}

SceneGraphElement::child_iterator SceneGraphElement::children_end()
{
	return child_iterator();
}

SceneGraphElement::const_child_iterator SceneGraphElement::children_begin() const
{
	return const_child_iterator( this->firstChild_ );
}

SceneGraphElement::const_child_iterator SceneGraphElement::children_end() const
{
	return const_child_iterator();
}


vl::Vec3f SceneGraphElement::furthestPoint( const vl::Vec3f& point ) const
{
	return point;
}

vl::Mat4f SceneGraphElement::transform( const char* state ) const
{
	if( state )
	{
		vl::Mat4f result = localTransform( state );
		assert( !isBad(result) );
		return result;
	}
	else
	{
		vl::Mat4f localTrans = this->localTransform(0);
		assert( !isBad(localTrans) );

		ConstSceneGraphElementPtr parent = this->parent().lock();
		if( parent )
			return parent->transform() * localTrans;
		else
			return localTrans;
	}
}

void SceneGraphElement::updateChildren()
{
	SceneGraphElementPtr parent = this->parent().lock();
	if( parent )
	{
		parent->updateChildren();
	}

	this->updateChildrenLocal();
}

void SceneGraphElement::updateTransform()
{
	for( SceneGraphElement::child_iterator childItr = this->children_begin();
		childItr != this->children_end(); ++childItr )
	{
		(*childItr)->updateTransform();
	}

	std::for_each( this->transformChangedRegistry_.begin(),
		this->transformChangedRegistry_.end(),
		boost::bind( &SceneGraphElement::updateTransformLocal, _1 ) );

	this->updateTransformLocal();
}

void SceneGraphElement::updateScale()
{
	for( SceneGraphElement::child_iterator childItr = this->children_begin();
		childItr != this->children_end(); ++childItr )
	{
		(*childItr)->updateScale();
	}

	this->updateScaleLocal();
}

void SceneGraphElement::updateChildrenLocal()
{
}

void SceneGraphElement::updateTransformLocal()
{
}

void SceneGraphElement::updateScaleLocal()
{
}

void SceneGraphElement::setParent( SceneGraphElementPtr parent )
{
	this->parent_ = parent;
}

boost::weak_ptr<SceneGraphElement> SceneGraphElement::parent()
{
	return this->parent_;
}

boost::weak_ptr<const SceneGraphElement> SceneGraphElement::parent() const
{
	return this->parent_;
}

// adds the given child just before the element specified, or at
//   end if not specified
void SceneGraphElement::addChild( SceneGraphElementPtr child, SceneGraphElementPtr after )
{
	this->checkConsistency();

	SceneGraphElementPtr before;
	if( after )
	{
		SceneGraphElementPtr prevChild = after->prevChild_.lock();
		before = prevChild;
		after->prevChild_ = child;
	}
	else
	{
		before = this->lastChild_;
		this->lastChild_ = child;
	}

	if( before )
		before->nextChild_ = child;
	else
		this->firstChild_ = child;

	child->prevChild_ = before;
	child->nextChild_ = after;

	this->updateChildren();
	this->checkConsistency();
}

// returns the object just after the removed object; that is,
//   addChild( child, removeChild(child) ) is the same as a noop
SceneGraphElementPtr SceneGraphElement::removeChild( SceneGraphElementPtr child )
{
	assert( child->parent().lock().get() == this );
	this->checkConsistency();

	assert( child );
	SceneGraphElementPtr before = child->prevChild_.lock();
	SceneGraphElementPtr after = child->nextChild_;

	if( child == this->firstChild_ )
		this->firstChild_ = after;

	if( child == this->lastChild_ )
		this->lastChild_ = before;

	if( before )
		before->nextChild_ = after;

	if( after )
		after->prevChild_ = before;

	child->prevChild_.reset();
	child->nextChild_.reset();

	this->updateChildren();
	this->checkConsistency();

	return after;
}

void SceneGraphElement::checkConsistency() const
{
#ifdef _DEBUG
	SceneGraphElementPtr current = this->firstChild_;
	if( this->firstChild_ )
		assert( !this->firstChild_->prevChild_.lock() );

	while( current )
	{
		if( current->nextChild_ )
			assert( current->nextChild_->prevChild_.lock() == current );
		else
			assert( current == this->lastChild_ );

		current = current->nextChild_;
	}
#endif
}

void SceneGraphElement::addTransformChangedTarget( SceneGraphElement::TransformChangedTarget target )
{
	this->transformChangedRegistry_.push_back( target );
}

void SceneGraphElement::removeTransformChangedTarget( SceneGraphElement::TransformChangedTarget target )
{
	this->transformChangedRegistry_.erase( 
		std::remove( transformChangedRegistry_.begin(), transformChangedRegistry_.end(), target ),
		transformChangedRegistry_.end() );
}

bool SceneGraphElement::visible() const
{
	return this->visible_;
}

void SceneGraphElement::setVisible( bool visible )
{
	this->visible_ = visible;
}

bool SceneGraphElement::selfCollisions() const
{
	return this->selfCollisions_;
}

void SceneGraphElement::setSelfCollisions( bool value )
{
	this->selfCollisions_ = value;
}

vl::Vec3f SceneGraphElement::color() const
{
	return vl::Vec3f( 0.6f, 0.6f, 0.6f );
}

SceneGraphElementPtr SceneGraphElement::clone() const
{
	SceneGraphElementPtr myClone( this->cloneLocal() );
	return myClone;
}

struct RandomizedPropertyComparator
{
	bool operator()( const RandomizedProperty& left, const RandomizedProperty& right )
	{
		std::less<std::string> comparator;
		return comparator(left.name, right.name);
	}
};

std::vector<RandomizedProperty> SceneGraphElement::randomizedProperties() const
{
	return std::vector<RandomizedProperty>( this->randomizedProperties_.begin(),
		this->randomizedProperties_.end() );
}

bool SceneGraphElement::isRandomized( const std::string& propertyName ) const
{
	typedef std::deque<RandomizedProperty>::const_iterator Iter;
	typedef std::pair<Iter, Iter> IterPr;
	IterPr pr = std::equal_range( this->randomizedProperties_.begin(),
		this->randomizedProperties_.end(), 
		RandomizedProperty(propertyName), 
		RandomizedPropertyComparator() );
	return (pr.first != pr.second);
}

void SceneGraphElement::setRandomized( const RandomizedProperty& property )
{
	typedef std::deque<RandomizedProperty>::iterator Iter;
	typedef std::pair<Iter, Iter> IterPr;
	IterPr pr = std::equal_range( this->randomizedProperties_.begin(),
		this->randomizedProperties_.end(), 
		property, 
		RandomizedPropertyComparator() );
	if( pr.first != pr.second )
		*pr.first = property;
	else
		this->randomizedProperties_.insert( pr.first, property );
}

RandomizedProperty SceneGraphElement::getRandomized( const std::string& propertyName ) const
{
	typedef std::deque<RandomizedProperty>::const_iterator Iter;
	typedef std::pair<Iter, Iter> IterPr;
	IterPr pr = std::equal_range( this->randomizedProperties_.begin(),
		this->randomizedProperties_.end(), 
		RandomizedProperty(propertyName), 
		RandomizedPropertyComparator() );
	assert( pr.first != pr.second );
	return *pr.first;
}

void SceneGraphElement::clearRandomized( const std::string& propertyName )
{
	typedef std::deque<RandomizedProperty>::iterator Iter;
	typedef std::pair<Iter, Iter> IterPr;
	IterPr pr = std::equal_range( this->randomizedProperties_.begin(),
		this->randomizedProperties_.end(), 
		RandomizedProperty(propertyName), 
		RandomizedPropertyComparator() );
	assert( pr.first != pr.second );
	this->randomizedProperties_.erase( pr.first, pr.second );
}

bool SceneGraphElement::equals( const SceneGraphElement& other, const SceneGraphElement* stopRecursingAt ) const
{
	if( !this->equalsLocal(other) )
		return false;

	ConstSceneGraphElementPtr myParent = this->parent().lock();
	ConstSceneGraphElementPtr otherParent = other.parent().lock();
	if( (myParent && !otherParent) || (otherParent && !myParent) )
		return false;

	if( myParent && (myParent.get() != stopRecursingAt) && !myParent->equals(*otherParent, stopRecursingAt) )
		return false;

	return true;
}

struct RandomizedNameEquals
{
	RandomizedNameEquals( const std::string& name )
		: name_(name) {}

	bool operator()( const RandomizedProperty& other )
	{
		return other.name == name_;
	}

private:
	std::string name_;
};

bool SceneGraphElement::equalsLocal( const SceneGraphElement& other ) const
{
	if( typeid( *this ) != typeid( other ) )
		return false;

	if( this->randomizedProperties_.size() != other.randomizedProperties_.size() )
		return false;

	for( std::deque<RandomizedProperty>::const_iterator myRandomized = randomizedProperties_.begin();
		myRandomized != randomizedProperties_.end(); ++myRandomized )
	{
		std::deque<RandomizedProperty>::const_iterator otherRandomized = std::find_if(
			other.randomizedProperties_.begin(), other.randomizedProperties_.end(), 
			RandomizedNameEquals( myRandomized->name ) );
		if( otherRandomized == other.randomizedProperties_.end() )
			return false;

		if( !withinEpsilon(myRandomized->low, otherRandomized->low) ||
			!withinEpsilon(myRandomized->high, otherRandomized->high) )
			return false;
	}

	return true;
}

void SceneGraphElement::render(const RenderState& rs, const vl::Mat4f& transform, const char* statePtr) const
{
#ifdef GUI
	if( !visible_ )
		return;

	vl::Vec4f color( this->color(), rs.alpha );
	if( rs.overrideColor )
		color = vl::Vec4f( rs.newColor, rs.alpha );

	vl::Mat4f transpose = trans(transform);
	GLMatrixStackHandler handler;
	glMultMatrixf( transpose.Ref() );

	if( rs.polygonMode != RenderState::NONE )
	{
		switch( rs.polygonMode )
		{
		case RenderState::POINTS:
			glPolygonMode( GL_FRONT, GL_POINT );
			break;
		case RenderState::LINES:
			glPolygonMode( GL_FRONT, GL_LINE );
			break;
		case RenderState::FILL:
			glPolygonMode( GL_FRONT, GL_FILL );
			break;
		}

		glEnable( GL_LIGHTING );
		glEnable( GL_NORMALIZE );

		boost::shared_ptr<GLEnableHandler> polygonOffset;
		if( rs.wireframeColor != vl::vl_0 )
		{
			polygonOffset.reset( new GLEnableHandler(GL_POLYGON_OFFSET_FILL) );
			glPolygonOffset (1., 1.);
		}

		const GLfloat no_mat[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		glMaterialfv(GL_FRONT, GL_AMBIENT, color.Ref());
		glMaterialfv(GL_FRONT, GL_DIFFUSE, color.Ref());
		glMaterialfv(GL_FRONT, GL_SPECULAR, no_mat);
		glMaterialfv(GL_FRONT, GL_SHININESS, no_mat);
		glMaterialfv(GL_FRONT, GL_EMISSION, no_mat );

		this->renderLocalFill( rs.alpha, statePtr );
	}

	if( rs.wireframeColor != vl::vl_0 )
	{
		glDisable( GL_LIGHTING );
		glColor4fv( rs.wireframeColor.Ref() );
		glLineWidth( 1.0 );

		renderLocalWire( rs.wireframeColor, statePtr );
	}

	if( rs.renderHulls )
	{
		glPolygonMode( GL_FRONT, GL_FILL );
		glEnable( GL_LIGHTING );
		glEnable( GL_NORMALIZE );
		renderLocalHulls( 1.0f, statePtr );
	}
#endif
}

void SceneGraphElement::renderPoints( const std::deque<unsigned int>& selected, 
	const vl::Mat4f& transform, 
	const char* state ) const
{
	// default: do nothing
}

std::pair<bool, unsigned int> SceneGraphElement::selectPoint( const vl::Vec2d& screenPos, 
	const ScreenSpaceConverter& converter, 
	const vl::Mat4f& transform,
	const char* state ) const
{
	return std::make_pair( false, 0 );
}

void SceneGraphElement::renderLocalHulls(
	const float alpha, 
	const char* statePtr) const
{
	// default: nothing
}

void SceneGraphElement::toMayaAscii( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const
{
	this->toMayaAsciiLocal( ofs, materialName, prefix );
}

void SceneGraphElement::toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const
{
	// default: do nothing for now
	// @todo fix this
}

boost::shared_ptr<RibNamedBlock> SceneGraphElement::toRibObject(
	const std::string& prefix,
	const std::string& timesFilename,
	const std::vector<const char*>& states) const
{
	boost::shared_ptr<RibNamedBlock> result( new RibNamedBlock("Attribute") );

	if( states.empty() )
	{
		result->addChild( ribObjectWithAttribute( "ConcatTransform", this->transform() ) );
	}
	else if( states.size() == 1 || timesFilename.empty() )
	{
		result->addChild( ribObjectWithAttribute( "ConcatTransform", this->transform(states.front()) ) );
	}
	else
	{
		std::string timelineFilename( prefix + "." + this->name() + ".timeline" );
		boost::filesystem::path timelinePath( timelineFilename );

		{
			RibBlock block;
			for( std::vector<const char*>::const_iterator stateIter = states.begin();
				stateIter != states.end(); ++stateIter )
			{
				block.addChild( ribObjectWithAttribute( "ConcatTransform", this->transform(*stateIter) ) );
			}

			std::ofstream ofs( timelineFilename.c_str() );
			if( !ofs )
				throw IOException( "Unable to open file '" + timelineFilename + "' for writing." );
			ofs << "1\n";
			block.dump( ofs, "" );
		}

		boost::shared_ptr<RibMetaComment> comment( new RibMetaComment(
			"METARIB-MOTION-TIMELINE " + timesFilename + " " + timelinePath.leaf() ) );
		result->addChild( comment );
	}

	result->addChild( ribObjectWithAttribute( "Color", this->color() ) );

	result->addChild( this->ribObject() );
	return result;
}

boost::shared_ptr<MechModelDeclaredObject> SceneGraphElement::toMechModel(MechModelType type) const
{
	boost::shared_ptr<MechModelDeclaredObject> object = this->toMechModelLocal(type);

	if( type == MECH_MODEL_AUDIO_ROOT )
		type = MECH_MODEL_AUDIO;

	for( SceneGraphElement::const_child_iterator iter = this->children_begin();
		iter != this->children_end(); ++iter )
	{
		boost::shared_ptr<MechModelDeclaredObject> childObject = 
			(*iter)->toMechModel(type);
		object->add( childObject );
	}

	return object;
}

boost::shared_ptr<MechModelDeclaredObject> SceneGraphElement::toMechModelLocal(MechModelType type) const
{
	boost::shared_ptr<MechModelDeclaredObject> object = this->mechModelObject();
	// for audio, we don't care about names.
	if( type == MECH_MODEL_AUDIO || type == MECH_MODEL_AUDIO_ROOT )
		object->setName( std::string() );

	if( !this->visible() && type != MECH_MODEL_AUDIO && type != MECH_MODEL_AUDIO_ROOT )
		object->add( mechModelValuePair("visible", this->visible_) );

	if( !this->selfCollisions() && type != MECH_MODEL_AUDIO && type != MECH_MODEL_AUDIO_ROOT )
		object->add( mechModelValuePair("selfCollisions", false) );

	if( type != MECH_MODEL_AUDIO && type != MECH_MODEL_AUDIO_ROOT )
	{
		for( std::deque<RandomizedProperty>::const_iterator itr = randomizedProperties_.begin();
			itr != randomizedProperties_.end(); ++itr )
		{
			boost::shared_ptr< MechModelValueList > mechModelValue( 
				new MechModelValueList("randomize") );
			mechModelValue->add( itr->name );
			mechModelValue->add( itr->low );
			mechModelValue->add( itr->high );
			object->add( mechModelValue );
		}
	}

	return object;
}

std::vector<std::string> SceneGraphElement::floatAttributes() const
{
	std::vector<std::string> result;
	return result;
}

void SceneGraphElement::setAttribute( const std::string& attributeName, float value )
{
	throw AttributeException( attributeName, "No such attribute" );
}

float SceneGraphElement::getAttribute( const std::string& attributeName ) const
{
	throw AttributeException( attributeName, "No such attribute" );
}

std::vector<std::string> SceneGraphElement::boolAttributes() const
{
	std::vector<std::string> result( 1, "visible" );
	result.push_back( "selfCollisions" );
	return result;
}

void SceneGraphElement::setBoolAttribute( const std::string& attributeName, bool value )
{
	if( attributeName == "visible" )
		this->visible_ = value;
	else if( attributeName == "selfCollisions" )
		this->selfCollisions_ = value;
	else
		throw AttributeException( attributeName, "No such attribute" );
}

bool SceneGraphElement::getBoolAttribute( const std::string& attributeName ) const
{
	if( attributeName == "visible" )
		return this->visible_;
	else if( attributeName == "selfCollisions" )
		return this->selfCollisions_;
	else
		throw AttributeException( attributeName, "No such attribute" );
}

class WrappedScreenSpaceConverter
	: public ScreenSpaceConverter
{
public:
	WrappedScreenSpaceConverter( const ScreenSpaceConverter& converter, const vl::Mat4d& transform )
		: transform_(transform), invTrans_( inv(transform) ), converter_(converter) {}

	vl::Vec3d toScreenSpace( const vl::Vec3d& worldSpace ) const
	{
		return converter_.toScreenSpace( xform(transform_, worldSpace) );
	}

	vl::Vec3d toWorldSpace( const vl::Vec3d& screenSpace ) const
	{
		return xform(invTrans_, converter_.toWorldSpace(screenSpace) );
	}

private:
	vl::Mat4d transform_;
	vl::Mat4d invTrans_;
	const ScreenSpaceConverter& converter_;
};

bool SceneGraphElement::intersect(const Ray3d& ray, double& t, const vl::Mat4f& transform, const char* statePtr) const
{
	if( !visible() )
		return false;

	vl::Mat4d localToGlobal = toMat4d( transform );
	vl::Mat4d globalToLocal = inv( localToGlobal );

	vl::Vec3d pos = xform(globalToLocal, ray.getPosition());
	vl::Vec3d dir = xform(globalToLocal, ray.getPosition() + ray.getDirection()) - pos;
	float length = len( dir );
	dir /= length;

    Ray3d localRay( pos, dir );
	if( intersectLocal(localRay, t, statePtr) )
	{
		// Transform the intersection point & normal returned back into global space.
		t /= length;

		return true;
	}
	else
	{
		return false;
	}
}

bool SceneGraphElement::intersect(
	const BoundingBox2d& box, 
	const ScreenSpaceConverter& converter, 
	const vl::Mat4f& transform, 
	const char* statePtr) const
{
	if( !visible() )
		return false;

	vl::Mat4d localToGlobal = toMat4d( transform );
	WrappedScreenSpaceConverter wrappedConverter( converter, localToGlobal );

    return intersectLocal( box, wrappedConverter, statePtr );
}

BoundingBox3f SceneGraphElement::bounds(const char* statePtr) const
{
	return this->bounds( this->transform(statePtr) );
}


TransformedSceneGraphElement::TransformedSceneGraphElement(const std::string& name)
	:	SceneGraphElement(name),
		translate_(0.0f, 0.0f, 0.0f),
		rotate_(0.0f, 0.0f, 0.0f),
		scale_(1.0f, 1.0f, 1.0f),
		pivot_(0.0f, 0.0f, 0.0f)
{
}

TransformedSceneGraphElement::~TransformedSceneGraphElement()
{
}

bool TransformedSceneGraphElement::equalsLocal( const SceneGraphElement& other ) const
{
	const TransformedSceneGraphElement* transformed = 
		dynamic_cast<const TransformedSceneGraphElement*>( &other );
	assert( transformed != 0 );

	if( !withinEpsilon( this->translate_, transformed->translate_ ) ||
		!withinEpsilon( this->rotate_, transformed->rotate_ ) ||
		!withinEpsilon( this->scale_, transformed->scale_ ) ||
		!withinEpsilon( this->pivot_, transformed->pivot_ ) )
		return false;

	return SceneGraphElement::equalsLocal( other );
}

vl::Vec3f TransformedSceneGraphElement::translate() const
{
	return this->translate_;
}

void TransformedSceneGraphElement::setTranslate( const vl::Vec3f& translate, bool temporary )
{
	assert( !isBad(translate) );
	this->translate_ = translate;

	if( !temporary )
	{
		this->updateTransform();
	}
}

vl::Vec3f TransformedSceneGraphElement::rotate() const
{
	return this->rotate_;
}

void TransformedSceneGraphElement::setRotate( const vl::Vec3f& rotate, bool temporary )
{
	assert( !isBad(rotate) );
	this->rotate_ = rotate;

	if( !temporary )
	{
		this->updateTransform();
	}
}

vl::Vec3f TransformedSceneGraphElement::pivot() const
{
	return this->pivot_;
}

void TransformedSceneGraphElement::setPivot( const vl::Vec3f& pivot, bool temporary )
{
	assert( !isBad(pivot) );
	this->pivot_ = pivot;

	if( !temporary )
	{
		this->updateTransform();
	}
}

vl::Vec3f TransformedSceneGraphElement::scale() const
{
	return scale_;
}

void TransformedSceneGraphElement::setScale( const vl::Vec3f& scale, bool temporary )
{
	assert( !isBad(scale) );

	// simulator is not going to deal well with negative scales
	//   I suspect since that will invert the triangle windings.
	for( vl::Int i = 0; i < 3; ++i )
		assert( scale[i] > 0.0 );

	assert( this->allowUniformScale() );

	if( !this->allowNonUniformScale() )
	{
		for( size_t i = 1; i < 3; ++i )
			assert( scale[i] == scale[0] );
	}

	scale_ = scale;

	if( !temporary )
	{
		this->updateTransform();
		this->updateScale();
	}
}

void addTripleAttributes( std::vector<std::string>& result, const char* prefix )
{
	for( size_t i = 0; i < 3; ++i )
	{
		std::string tmp( prefix );
		tmp.push_back( '.' );
		tmp.push_back( 'x' + i );
		result.push_back( tmp );
	}
}

std::vector<std::string> TransformedSceneGraphElement::floatAttributes() const
{
	std::vector<std::string> result = SceneGraphElement::floatAttributes();
	result.reserve( result.size()+12 );
	addTripleAttributes( result, "translate" );
	addTripleAttributes( result, "rotate" );
	addTripleAttributes( result, "scale" );
	addTripleAttributes( result, "pivot" );
	return result;
}

void TransformedSceneGraphElement::setAttribute( const std::string& attributeName, float value )
{
	std::string::size_type dotPos = attributeName.rfind(".");
	if( dotPos == std::string::npos )
		SceneGraphElement::setAttribute( attributeName, value );

	std::string prefix = attributeName.substr(0, dotPos);
	int component = (attributeName[attributeName.size() - 1] - 'x');
	if( component < 0 || component > 3 )
		SceneGraphElement::setAttribute( attributeName, value );

	if( prefix == "translate" )
		this->translate_[component] = value;
	else if( prefix == "rotate" )
	{
		this->rotate_[component] = value * (M_PI/180.0);
	}
	else if( prefix == "scale" )
		this->scale_[component] = value;
	else if( prefix == "pivot" )
		this->pivot_[component] = value;
	else
		SceneGraphElement::setAttribute( attributeName, value );
}

float TransformedSceneGraphElement::getAttribute( const std::string& attributeName ) const
{
	std::string::size_type dotPos = attributeName.rfind(".");
	if( dotPos == std::string::npos )
		return SceneGraphElement::getAttribute( attributeName );

	std::string prefix = attributeName.substr(0, dotPos);
	int component = (attributeName[attributeName.size() - 1] - 'x');
	if( component < 0 || component > 3 )
		return SceneGraphElement::getAttribute( attributeName );

	if( prefix == "translate" )
		return this->translate_[component];
	else if( prefix == "rotate" )
		return this->rotate_[component] * (180.0/M_PI);
	else if( prefix == "scale" )
		return this->scale_[component];
	else if( prefix == "pivot" )
		return this->pivot_[component];
	else
		return SceneGraphElement::getAttribute( attributeName );
}

boost::shared_ptr<MechModelDeclaredObject> TransformedSceneGraphElement::toMechModelLocal(MechModelType type) const
{
	boost::shared_ptr<MechModelDeclaredObject> object = SceneGraphElement::toMechModelLocal(type);

	if( type != MECH_MODEL_AUDIO_ROOT )
	{
		object->add( mechModelValuePair("translate", this->translate_) );
		object->add( mechModelValuePair("rotate", this->rotate_) );
		object->add( mechModelValuePair("pivot", this->pivot_) );
	}

	object->add( mechModelValuePair("scale", this->scale_) );

	return object;
}

vl::Mat4f TransformedSceneGraphElement::localTransform( const char* state ) const
{
	// I don't know what to do with State
	assert( state == 0 );

	vl::Mat4f pivotMatrix = vl::HTrans4f( -this->pivot_ );

	vl::Mat4f scaleMatrix = vl::HScale4f( this->scale_ );
	vl::Mat4f translateMatrix = vl::HTrans4f( this->translate_ );

	Quaternion<float> quat = eulerToQuat( this->rotate_ );
	vl::Mat4f rotateMatrix( vl::vl_1 );
	{
		vl::Mat3f rotMat = quat.toRotMatf();
		for( vl::Int i = 0; i < 3; ++i )
			for( vl::Int j = 0; j < 3; ++j )
				rotateMatrix[i][j] = rotMat[i][j];
	}

	assert( !isBad(scaleMatrix) );
	assert( !isBad(translateMatrix) );
	assert( !isBad(rotateMatrix) );
	assert( !isBad(pivotMatrix) );

	return translateMatrix*rotateMatrix*scaleMatrix*pivotMatrix;
}

TransformGroup::TransformGroup( const std::string& name )
	: TransformedSceneGraphElement( name )
{
}

TransformGroup::~TransformGroup()
{
}

void TransformGroup::renderLocalFill(
	float alpha, 
	const char* state) const
{
}

void TransformGroup::renderLocalWire(
	const vl::Vec4f& wireframeColor, 
	const char* state) const
{
}

boost::shared_ptr<MechModelDeclaredObject> TransformGroup::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "group", this->name() ) );
	return result;
}

bool TransformGroup::intersectLocal(
	const Ray3d& ray, 
	double& t, 
	const char* state) const
{
	return false;
}

bool TransformGroup::intersectLocal(
	const BoundingBox2d& box, 
	const ScreenSpaceConverter& converter, 
	const char* state) const
{
	return false;
}

BoundingBox3f TransformGroup::bounds( const vl::Mat4f& transform ) const
{
	return BoundingBox3f();
}

RibObjectPtr TransformGroup::ribObject() const
{
	return RibObjectPtr();
}

TransformGroup* TransformGroup::cloneLocal() const
{
	return new TransformGroup( *this );
}

vl::Mat4f rigidTransform(const RigidStaticState& state)
{
	/*
	vl::Vec3d orientation = toVec3d( state.orientation() );
	double len = vl::len( orientation );
	vl::Mat4d rotMat( vl::vl_1 );
	if( len > 1e-4 )
	{
		orientation /= len;
		rotMat = vl::HRot4d( orientation, len );
	}
	*/

	vl::Mat4f mat( vl::vl_1 );
	{
		vl::Mat3f rotMat = state.orientation().toRotMatf();
		for( vl::Int i = 0; i < 3; ++i )
			for( vl::Int j = 0; j < 3; ++j )
				mat[i][j] = rotMat[i][j];
	}

	for( vl::Int i = 0; i < 3; ++i )
		mat[i][3] = (state.position())[i];

	assert( !isBad(mat) );

	return mat;
}

RigidStaticState::RigidStaticState( const vl::Vec3f& position, const OrientationType& orientation )
	:	position_(position),
		orientation_( orientation)
{
	assert( !isBad(position_) );
	assert( !isBad(orientation_) );
}

RigidStaticState::RigidStaticState()
	:	position_(vl::vl_0), 
		orientation_( vl::Vec3f(vl::vl_0), 1.0f)
{
}

#ifdef CONDOR_SAMPLING
RigidStaticState::RigidStaticState( MWRMComm * RMC )
{
	float pos[3];
	RMC->unpack( pos, 3 );
	this->position_.set( pos );

	float quat[4];
	RMC->unpack( quat, 4 );
	this->orientation_.setWXYZ( quat );
}

void RigidStaticState::write( MWRMComm * RMC ) const
{
	RMC->pack( this->position_.get(), 3 );
	float quat[4];
	this->orientation_.getWXYZ( quat );
	RMC->pack( quat, 4 );
}
#endif


RigidDynamicState::RigidDynamicState( 
	const vl::Vec3f& position, 
	const OrientationType& orientation, 
	const vl::Vec3f& linearVelocity, 
	const vl::Vec3f& angularVelocity )
	:	RigidStaticState(position, orientation),
		linearVelocity_(linearVelocity),
		angularVelocity_(angularVelocity)
{
	assert( !isBad(linearVelocity_) );
	assert( !isBad(angularVelocity_) );
}

RigidDynamicState::RigidDynamicState(const RigidStaticState& state)
	:	RigidStaticState(state), 
		linearVelocity_(vl::vl_0),
		angularVelocity_(vl::vl_0)
{
}

RigidDynamicState::RigidDynamicState()
	:	linearVelocity_(vl::vl_0),
		angularVelocity_(vl::vl_0)
{
}


vl::Vec3f PhysicsObject::color() const
{
	return this->material()->color();
}

#ifdef SOUND
void PhysicsObject::setSquashingCubesModes( boost::shared_ptr<SquashingCubesModes> modes )
{
	this->modes_ = modes;
}
#endif

void PhysicsObject::voxelize( BoundedVoxelGrid& grid, 
	const std::vector<ConstMaterialPtr>& availableMaterials,
	const vl::Mat4f& worldToGrid ) const
{
	assert( this->hasAudio() );
}

#ifdef SOUND
boost::shared_ptr<AudioGenerator> PhysicsObject::getAudioGenerator(size_t rate) const
{
	if( modes_ && this->hasAudio() )
	{
		assert( !modes_->frequencies.empty() );
		double baseFrequency = modes_->frequencies.front();
		double realAlpha = this->material()->raleighAlpha() * baseFrequency;
		double realBeta = this->material()->raleighBeta() / baseFrequency;

		return boost::shared_ptr<AudioGenerator>( new SquashingCubesAudioGenerator(rate, modes_, realAlpha, realBeta) );
	}
	else
		return boost::shared_ptr<AudioGenerator>();
}
#endif

boost::shared_ptr<MechModelDeclaredObject> PhysicsObject::toMechModelLocal(MechModelType type) const
{
	boost::shared_ptr<MechModelDeclaredObject> object = 
		TransformedSceneGraphElement::toMechModelLocal(type);

	if( this->isStatic() && type != MECH_MODEL_AUDIO && type != MECH_MODEL_AUDIO_ROOT )
		object->add( mechModelValuePair("static", true ) );

	if( !this->hasMass() && type != MECH_MODEL_AUDIO && type != MECH_MODEL_AUDIO_ROOT )
		object->add( mechModelValuePair("hasMass", false ) );

	if( !this->hasAudio() && type != MECH_MODEL_AUDIO && type != MECH_MODEL_AUDIO_ROOT )
		object->add( mechModelValuePair("hasAudio", false ) );

	if( !this->collides() && type != MECH_MODEL_AUDIO && type != MECH_MODEL_AUDIO_ROOT )
		object->add( mechModelValuePair("collides", false ) );

	if( sqrlen( this->linearVelocity_ ) != 0.0f && type != MECH_MODEL_AUDIO && type != MECH_MODEL_AUDIO_ROOT )
		object->add( mechModelValuePair("linearVelocity", this->linearVelocity_) );

	if( sqrlen( this->angularVelocity_ ) != 0.0f && type != MECH_MODEL_AUDIO && type != MECH_MODEL_AUDIO_ROOT )
		object->add( mechModelValuePair("angularVelocity", this->angularVelocity_) );

	// for audio mech models, we embed the material properties
	if( type == MECH_MODEL_AUDIO || type == MECH_MODEL_AUDIO_ROOT )
		object->add( material()->toMechModel(type) );
	else
		object->add( mechModelValuePair("material", material()->name()) );

	return object;
}

vl::Vec3f PhysicsObject::linearVelocity() const
{
	return this->linearVelocity_;
}

void PhysicsObject::setLinearVelocity( const vl::Vec3f& linVel )
{
	this->linearVelocity_ = linVel;
}

vl::Vec3f PhysicsObject::angularVelocity() const
{
	return this->angularVelocity_;
}

void PhysicsObject::setAngularVelocity( const vl::Vec3f& angVel )
{
	this->angularVelocity_ = angVel;
}

PhysicsObject::PhysicsObject( const std::string& name, boost::shared_ptr<Material> material )
	:	TransformedSceneGraphElement(name), 
		isStatic_(false), 
		hasMass_(true),
		collides_(true),
		hasAudio_(true),
		material_(material),
		linearVelocity_(vl::vl_0),
		angularVelocity_(vl::vl_0)
{
}

PhysicsObject::~PhysicsObject()
{
}

bool PhysicsObject::equalsLocal( const SceneGraphElement& other ) const
{
	const PhysicsObject* po = dynamic_cast<const PhysicsObject*>(&other);
	assert( po != 0 );
	if( !this->material_->equals( *po->material_ ) )
		return false;

	if( (this->isStatic_ != po->isStatic_) || 
		(this->hasMass_ != po->hasMass_) ||
		(this->collides_ != po->collides_) ||
		!withinEpsilon(this->angularVelocity_, po->angularVelocity_) ||
		!withinEpsilon(this->linearVelocity_, po->linearVelocity_) )
		return false;

	return TransformedSceneGraphElement::equalsLocal( other );

}

void PhysicsObject::setStatic( bool isStatic )
{
	this->isStatic_ = isStatic;
}

bool PhysicsObject::isStatic() const
{
	return this->isStatic_;
}

bool PhysicsObject::collides() const
{
	return this->collides_;
}

void PhysicsObject::setCollides( bool collides )
{
	this->collides_ = collides;
}

bool PhysicsObject::hasMass() const
{
	return this->hasMass_;
}

void PhysicsObject::setHasMass( bool hasMass )
{
	this->hasMass_ = hasMass;
}

bool PhysicsObject::hasAudio() const
{
	return this->hasAudio_;
}

void PhysicsObject::setHasAudio( bool hasAudio )
{
	this->hasAudio_ = hasAudio;
}


std::vector<std::string> PhysicsObject::boolAttributes() const
{
	std::vector<std::string> result = TransformedSceneGraphElement::boolAttributes();
	result.push_back( "static" );
	result.push_back( "hasMass" );
	result.push_back( "hasAudio" );
	result.push_back( "collides" );
	return result;
}

void PhysicsObject::setBoolAttribute( const std::string& attributeName, bool value )
{
	if( attributeName == "static" )
		this->isStatic_ = value;
	else if( attributeName == "hasMass" )
		this->hasMass_ = value;
	else if( attributeName == "hasAudio" )
		this->hasAudio_ = value;
	else if( attributeName == "collides" )
		this->collides_ = value;
	else
		TransformedSceneGraphElement::setBoolAttribute( attributeName, value );

}

bool PhysicsObject::getBoolAttribute( const std::string& attributeName ) const
{
	if( attributeName == "static" )
		return this->isStatic_;
	else if( attributeName == "hasMass" )
		return this->hasMass_;
	else if( attributeName == "hasAudio" )
		return this->hasAudio_;
	else if( attributeName == "collides" )
		return this->collides_;
	else
		return TransformedSceneGraphElement::getBoolAttribute( attributeName );
}

std::vector<std::string> PhysicsObject::floatAttributes() const
{
	std::vector<std::string> result = TransformedSceneGraphElement::floatAttributes();
	result.reserve( result.size()+6 );
	addTripleAttributes( result, "linearVelocity" );
	addTripleAttributes( result, "angularVelocity" );\
	return result;
}

void PhysicsObject::setAttribute( const std::string& attributeName, float value )
{
	std::string::size_type dotPos = attributeName.rfind(".");
	if( dotPos == std::string::npos )
		TransformedSceneGraphElement::setAttribute(attributeName, value);

	std::string prefix = attributeName.substr(0, dotPos);
	int component = (attributeName[attributeName.size() - 1] - 'x');
	if( component < 0 || component > 3 )
		TransformedSceneGraphElement::setAttribute(attributeName, value);

	if( prefix == "linearVelocity" )
		this->linearVelocity_[component] = value;
	else if( prefix == "angularVelocity" )
		this->angularVelocity_[component] = value;
	else
		TransformedSceneGraphElement::setAttribute(attributeName, value);
}

float PhysicsObject::getAttribute( const std::string& attributeName ) const
{
	std::string::size_type dotPos = attributeName.rfind(".");
	if( dotPos == std::string::npos )
		return TransformedSceneGraphElement::getAttribute(attributeName);

	std::string prefix = attributeName.substr(0, dotPos);
	int component = (attributeName[attributeName.size() - 1] - 'x');
	if( component < 0 || component > 3 )
		return TransformedSceneGraphElement::getAttribute(attributeName);

	if( prefix == "linearVelocity" )
		return this->linearVelocity_[component];
	else if( prefix == "angularVelocity" )
		return this->angularVelocity_[component];
	else
		return TransformedSceneGraphElement::getAttribute(attributeName);

}

ODEGeomResult PhysicsObject::odeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const
{
	return getOdeGeom(space, centerOfMass, padding);
}

#ifdef USE_NEWTON
NewtonCollision* PhysicsObject::newtonCollision( const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const
{
	return this->getNewtonCollision( newtonWorld, vl::trans(offsetMatrix), padding );
}
#endif

#ifdef USE_BULLET
btCollisionShape* PhysicsObject::bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const
{
	return 0;
}
#endif


ConstMaterialPtr PhysicsObject::material() const
{
	return material_;
}

boost::shared_ptr<Material> PhysicsObject::material()
{
	return material_;
}

void PhysicsObject::setMaterial( boost::shared_ptr<Material> material )
{
	material_ = material;
}

vl::Vec3f PhysicsObject::fullScale() const
{
	vl::Mat3f S = polarDecomposition( toMat3f(this->transform()) ).second;
	vl::Vec3f result;
	for( size_t i = 0; i < 3; ++i )
		result[i] = S[i][i];

	return result;
}

vl::Mat4f PhysicsObject::localTransform( const char* state ) const
{
	if( state != 0 )
	{
		const RigidStaticState* rigidState = 
			reinterpret_cast<const RigidStaticState*>(state);

		vl::Mat4f scaleMatrix = vl::HScale4f( this->fullScale() );
		assert( !isBad(scaleMatrix) );

		return planning::rigidTransform( *rigidState ) * scaleMatrix;
	}

	return TransformedSceneGraphElement::localTransform(0);
}

vl::Mat4f PhysicsObject::rigidTransform( const char* state ) const
{
	vl::Mat4f rotation( vl::vl_1 );

	vl::Mat4f transform = this->transform(state);
	vl::Mat3f rotate = polarDecomposition( toMat3f(transform) ).first;
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			rotation[i][j] = rotate[i][j];
	assert( !isBad(rotation) );

	vl::Vec3f position = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	assert( !isBad(position) );

	return vl::HTrans4f( position ) * rotation;
}

CombinedObject::CombinedObject( const std::string& name, boost::shared_ptr<Material> material )
	: PhysicsObject(name, material)
{
}

CombinedObject::~CombinedObject()
{
}

bool CombinedObject::equalsLocal( const SceneGraphElement& other ) const
{
	// if we have a CombinedObject as a parent, it's already checking all our
	// children so don't bother.
	ConstSceneGraphElementPtr parent = this->parent().lock();
	while( parent )
	{
		if( boost::dynamic_pointer_cast<const CombinedObject>( parent ) )
			return true;

		parent = parent->parent().lock();
	}

	std::deque<ConstSceneGraphElementPtr> myChildren( this->children_begin(), this->children_end() );
	std::deque<ConstSceneGraphElementPtr> otherChildren( other.children_begin(), other.children_end() );

	// need to check my children
	while( !myChildren.empty() && !otherChildren.empty() )
	{
		ConstSceneGraphElementPtr myChild = myChildren.front();         myChildren.pop_front();
		ConstSceneGraphElementPtr otherChild = otherChildren.front();   otherChildren.pop_front();
		if( !myChild->equalsLocal( *otherChild ) )
			return false;

		std::copy( myChild->children_begin(), myChild->children_end(), 
			std::front_inserter(myChildren) );
		std::copy( otherChild->children_begin(), otherChild->children_end(), 
			std::front_inserter(otherChildren) );
	}

	if( !myChildren.empty() || !otherChildren.empty() )
		return false;

	return PhysicsObject::equalsLocal(other);
}

void CombinedObject::updateTransformLocal()
{
}

void CombinedObject::updateChildrenLocal()
{
	for( std::vector< boost::weak_ptr<PhysicsObject> >::iterator itr = this->objects_.begin();
		itr != this->objects_.end(); ++itr )
	{
		itr->lock()->removeTransformChangedTarget( this );
	}
	objects_.clear();

	std::deque<SceneGraphElementPtr> queue( this->children_begin(), this->children_end() );
	while( !queue.empty() )
	{
		SceneGraphElementPtr current = queue.front();
		queue.pop_front();

		std::copy( current->children_begin(), current->children_end(),
			std::front_inserter(queue) );

		PhysicsObjectPtr object = boost::dynamic_pointer_cast<PhysicsObject>( current );
		boost::shared_ptr<CombinedObject> combined = boost::dynamic_pointer_cast<CombinedObject>( object );
		if( object && !combined )
		{
			object->addTransformChangedTarget( this );
			objects_.push_back( object );
		}
	}
}

vl::Vec3f CombinedObject::furthestPoint( const vl::Vec3f& point ) const
{
	vl::Vec3f furthest = point;
	for( std::vector< boost::weak_ptr<PhysicsObject> >::const_iterator iter = this->objects_.begin();
		iter != this->objects_.end(); ++iter )
	{
		ConstPhysicsObjectPtr object = iter->lock();

		// need to transform point into object's local coordinates
		vl::Mat4f transform = rigidInverse(this->rigidTransform()) * object->rigidTransform();
		vl::Vec3f curPt = vl::xform( transform, 
			object->furthestPoint( vl::xform(rigidInverse(transform), point) ) );
		if( vl::len( curPt - point ) > vl::len( furthest - point ) )
			furthest = curPt;
	}

	return furthest;
}

dMass CombinedObject::getOdeMassProps() const
{
	dMass result;
	dMassSetZero( &result );

	for( std::vector< boost::weak_ptr<PhysicsObject> >::const_iterator iter = this->objects_.begin();
		iter != this->objects_.end(); ++iter )
	{
		ConstPhysicsObjectPtr object = iter->lock();
		assert( object );
		if( !object->hasMass() )
			continue;

		dMass childMass = object->getOdeMassProps();
		assert( my_checkMass(&childMass) );
		vl::Vec3f childMassCenterOfMass( childMass.c[0], childMass.c[1], childMass.c[2] );

		vl::Mat4f odeRotation;
		vl::Vec3f position;

		{
			vl::Mat4f transform = rigidInverse(this->rigidTransform()) * object->rigidTransform();
			vl::Mat3f rotate = polarDecomposition( toMat3f(transform) ).first;
			for( vl::Int i = 0; i < 3; ++i )
				for( vl::Int j = 0; j < 3; ++j )
					odeRotation[i][j] = rotate[i][j];

			position = vl::xform( transform, childMassCenterOfMass );
		}
 
		std::vector<dReal> rotation( odeRotation.Ref(), odeRotation.Ref()+12 );
		dMassTranslate( &childMass, -childMass.c[0], -childMass.c[1], -childMass.c[2] );
		dMassRotate( &childMass, &rotation[0] );
		dMassTranslate( &childMass, position[0], position[1], position[2] );
		dMassAdd( &result, &childMass );
		assert( my_checkMass(&result) );
	}

	return result;
}

bool CombinedObject::intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const
{
	return false;
}

bool CombinedObject::intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const
{
	return false;
}

boost::shared_ptr<MechModelDeclaredObject> CombinedObject::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "fused", this->name() ) );
	return result;
}

RibObjectPtr CombinedObject::ribObject() const
{
	return RibObjectPtr();
}

void CombinedObject::renderLocalFill(float alpha, const char* statePtr) const
{
}

void CombinedObject::renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const
{
}

CombinedObject* CombinedObject::cloneLocal() const
{
	return new CombinedObject( *this );
}

BoundingBox3f CombinedObject::bounds( const vl::Mat4f& transform ) const
{
	return BoundingBox3f();
}

ODEGeomResult CombinedObject::getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const
{
	ODEGeomResult result;

	for( std::vector< boost::weak_ptr<PhysicsObject> >::const_iterator iter = this->objects_.begin();
		iter != this->objects_.end(); ++iter )
	{
		ConstPhysicsObjectPtr object = iter->lock();
		assert( object );
		if( !object->collides() )
			continue;

		ODEGeomResult geoms = object->odeGeom(0, vl::Vec3f(vl::vl_0), padding);
#if dTRIMESH_ENABLED
		std::copy( geoms.triMeshes.begin(), geoms.triMeshes.end(), 
			std::back_inserter(result.triMeshes) );
#endif

		for( std::vector<dGeomID>::const_iterator geomItr = geoms.geoms.begin();
			geomItr != geoms.geoms.end(); ++geomItr )
		{
			vl::Mat4f odeRotation( vl::vl_1 );
			vl::Vec3f position;

			{
				vl::Mat4f transform = rigidInverse(this->rigidTransform() * vl::HTrans4f(centerOfMass)) * object->rigidTransform();
				for( vl::Int i = 0; i < 3; ++i )
					for( vl::Int j = 0; j < 3; ++j )
						odeRotation[i][j] = transform[i][j];

				position = vl::xform( transform, vl::Vec3f(vl::vl_0) );
			}
	 
			dGeomID transformId = dCreateGeomTransform(space);
			dGeomTransformSetGeom( transformId, *geomItr );
			dGeomTransformSetCleanup( transformId, 1 );
			dGeomTransformSetInfo( transformId, 1 );

			std::vector<dReal> rotation( odeRotation.Ref(), odeRotation.Ref() + 12 );
			dGeomSetRotation ( *geomItr, &rotation[0] );
			dGeomSetPosition( *geomItr, position[0], position[1], position[2] );

			result.geoms.push_back( transformId );
		}
	}

	return result;
}

#ifdef USE_BULLET
btCollisionShape* CombinedObject::bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const
{
	btCompoundShape* result = new btCompoundShape();
	for( std::vector< boost::weak_ptr<PhysicsObject> >::const_iterator iter = this->objects_.begin();
		iter != this->objects_.end(); ++iter )
	{
		ConstPhysicsObjectPtr object = iter->lock();
		assert( object );
		if( !object->collides() )
			continue;

		vl::Mat4f localTransform = 
			rigidInverse(this->rigidTransform() * offsetMatrix) * object->rigidTransform();
		btTransform btLocalTransform( toBtMatrix3x3(localTransform), 
			btVector3(localTransform[0][3], localTransform[1][3], localTransform[2][3]) );

		// Anything that's currently in a CombinedObject needs to be broken out
		vl::Mat4f identity( vl::vl_1 );
		btCollisionShape* shape = object->bulletCollisionShape(identity, padding);

		btCompoundShape* compoundShape = 
			dynamic_cast<btCompoundShape*>( shape );
		if( compoundShape )
		{
			for( int i = 0; i < compoundShape->getNumChildShapes(); ++i )
			{
				btCollisionShape* childShape = compoundShape->getChildShape(i);
				btTransform childTransform = compoundShape->getChildTransform(i);
				btTransform newChildTransform = btLocalTransform * childTransform;
				result->addChildShape( newChildTransform, childShape );
			}

			delete compoundShape;
		}
		else if( shape )
		{
			result->addChildShape( btLocalTransform, shape );
		}
	}

	if( result->getNumChildShapes() == 0 )
	{
		delete result;
		return 0;
	}

	return result;
}
#endif

#ifdef USE_NEWTON
NewtonCollision* CombinedObject::getNewtonCollision( 
	const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const
{
	std::vector<NewtonCollision*> childCollisions;
	childCollisions.reserve( this->objects_.size() );

	for( std::vector< boost::weak_ptr<PhysicsObject> >::const_iterator iter = this->objects_.begin();
		iter != this->objects_.end(); ++iter )
	{
		ConstPhysicsObjectPtr object = iter->lock();
		assert( object );
		if( !object->collides() )
			continue;

		dFloatMat4 childOffsetMat = rigidInverse(this->rigidTransform() * trans(offsetMatrix)) * object->rigidTransform();

		NewtonCollision* geom = object->newtonCollision( newtonWorld, childOffsetMat, padding );
		if( !geom )
			continue;

		childCollisions.push_back( geom );
	}

	NewtonCollision* result = NewtonCreateCompoundCollision(
		newtonWorld,
		childCollisions.size(), 
		&childCollisions[0] );

	//std::for_each( childCollisions.begin(), childCollisions.end(), 
	//	boost::bind( &NewtonReleaseCollision, newtonWorld, _1 ) );

	return result;
}
#endif

#ifdef USE_NOVODEX
std::vector< boost::shared_ptr<NxShapeDesc> > CombinedObject::getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const
{
	std::vector< boost::shared_ptr<NxShapeDesc> > result;

	std::deque<ConstSceneGraphElementPtr> queue( this->children_begin(), this->children_end() );
	while( !queue.empty() )
	{
		ConstSceneGraphElementPtr current = queue.front();
		queue.pop_front();

		std::copy( current->children_begin(), current->children_end(),
			std::front_inserter(queue) );

		ConstPhysicsObjectPtr object = boost::dynamic_pointer_cast<const PhysicsObject>( current );
		if( !object )
			continue;

		const std::vector< boost::shared_ptr<NxShapeDesc> > descs = 
			object->getNxShapeDesc(sdk, padding);


		NxMat33 nxRotation;
		NxVec3 position;

		{
			vl::Mat4f transform = rigidInverse(this->rigidTransform()) * object->transform();
			vl::Mat3f rotate = polarDecomposition( toMat3f(transform) ).first;
			for( vl::Int i = 0; i < 3; ++i )
				for( vl::Int j = 0; j < 3; ++j )
					nxRotation(i, j) = rotate[i][j];

			position = toNxVec( vl::xform( transform, vl::Vec3f(vl::vl_0) ) );
		}

		for( std::vector< boost::shared_ptr<NxShapeDesc> >::const_iterator shapeItr = descs.begin();
			shapeItr != descs.end(); ++shapeItr )
		{
			(*shapeItr)->localPose = NxMat34(nxRotation, position)*(*shapeItr)->localPose;

			result.push_back( *shapeItr );
		}
	}

	return result;
}
#endif

BoxObject::BoxObject( const std::string& name, boost::shared_ptr<Material> material )
	: PhysicsObject( name, material )
{
}


BoundingBox3f BoxObject::bounds( const vl::Mat4f& transform ) const
{
	// in this case, transformBounds is exactly the right thing to do:
	BoundingBox3f localBounds( vl::Vec3f(-0.5, -0.5, -0.5), vl::Vec3f(0.5, 0.5, 0.5) );
	return transformBounds( localBounds, transform );
}

void BoxObject::voxelize( BoundedVoxelGrid& grid, 
	const std::vector<ConstMaterialPtr>& availableMaterials,
	const vl::Mat4f& worldToGrid ) const
{
	assert( this->hasAudio() );

	// find our material in the available material list:
	typedef std::vector<ConstMaterialPtr>::const_iterator Itr;
	typedef std::pair<Itr, Itr> ItrPr;
	ItrPr pr = std::equal_range( availableMaterials.begin(), availableMaterials.end(), 
		this->material(), MaterialSoundPropertiesOrder() );
	assert( pr.first != pr.second );

	// plus 1 because 0 is reserved for empty voxel
	size_t iMat = std::distance( availableMaterials.begin(), pr.first ) + 1;

	// transform box into the grid's local coordinate frame
	vl::Mat4f boxToGrid = worldToGrid * this->transform();

	Wm4::Box3f wmlMyBox;
	vl::Vec3f center = vl::xform( boxToGrid, vl::Vec3f(vl::vl_0) );
	wmlMyBox.Center = toMagicVec( center );
	for( vl::Int i = 0; i < 3; ++i )
	{
		vl::Vec3f endPointLocal( vl::vl_0 );
		endPointLocal[i] = 0.5f;
		vl::Vec3f endPoint = vl::xform( boxToGrid, endPointLocal );
		vl::Vec3f direction = endPoint - center;
		float len = vl::len(direction);
		direction /= len;

		wmlMyBox.Axis[i] = toMagicVec( direction );
		wmlMyBox.Extent[i] = len;
	}

	// for each voxel grid cell, test against the box:
	BoundingBox3f gridBox = grid.box( TinyVec<size_t, 3>(0, 0, 0) );
	Wm4::Box3<float> wmlGridBox;
	wmlGridBox.Axis[0] = Wm4::Vector3f( 1.0f, 0.0f, 0.0f );
	wmlGridBox.Axis[1] = Wm4::Vector3f( 0.0f, 1.0f, 0.0f );
	wmlGridBox.Axis[2] = Wm4::Vector3f( 0.0f, 0.0f, 1.0f );
	for( size_t i = 0; i < 3; ++i )
		wmlGridBox.Extent[i] = 0.5*((gridBox.maximum())[i] - (gridBox.minimum())[i]);

	TinyVec<size_t, 3> dimensions = grid.dimensions();
	for( size_t i = 0; i < dimensions[0]; ++i )
	{
		for( size_t j = 0; j < dimensions[1]; ++j )
		{
			for( size_t k = 0; k < dimensions[2]; ++k )
			{
				TinyVec<size_t, 3> cell(i, j, k);
				wmlGridBox.Center = toMagicVec( grid.center( cell ) );
				Wm4::IntrBox3Box3f intersect( wmlMyBox, wmlGridBox );
				if( intersect.Test() )
					grid.set( cell, iMat );
			}
		}
	}
}

BoxObject* BoxObject::cloneLocal() const
{
	return new BoxObject( *this );
}

RibObjectPtr BoxObject::ribObject() const
{
	return ribObjectWithAttribute("ObjectInstance", 1);
}

void BoxObject::toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const
{
	std::string name = prefix + this->name();
	std::string shapeName = name + "Shape";
	std::string cubeShapeName = name + "PolyCube";
	ofs << "createNode mesh -n \"" << shapeName << "\" -p \"" << name << "\";\n";
	ofs << "createNode polyCube -n \"" << cubeShapeName << "\";\n";
    //setAttr( ofs, ".createUVs", 3 );
	ofs << "connectAttr \"" << cubeShapeName << ".out\" \"" << shapeName << ".i\";\n";
	ofs << "connectAttr \"" << shapeName << ".iog\" \"" << materialName << ".dsm\" -na;\n";
    ofs << "setAttr \"" << cubeShapeName << ".createUVs\" 3;\n";
}

boost::shared_ptr<MechModelDeclaredObject> BoxObject::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "box", this->name() ) );

	return result;
}

ODEGeomResult BoxObject::getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const
{
	assert( centerOfMass == vl::vl_0 );
	vl::Vec3f scale = this->fullScale();
	return ODEGeomResult( dCreateBox( space, scale[0] + padding, scale[1] + padding, scale[2] + padding ) );
}

#ifdef USE_NEWTON
NewtonCollision* BoxObject::getNewtonCollision( 
	const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const
{
	vl::Vec3f scale = this->fullScale();
	return NewtonCreateBox(	newtonWorld, scale[0] + padding, scale[1] + padding, scale[2] + padding, offsetMatrix.Ref() );
}
#endif

#ifdef USE_BULLET
btCollisionShape* BoxObject::bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const
{
	vl::Vec3f scale = 0.5*this->fullScale();
	vl::Mat3f rotMat = vl::trans( toMat3f( offsetMatrix ) );

	const float epsilon = 1e-4;
	{
		// dont know how to apply a translation offset
		vl::Vec3f translate( offsetMatrix[0][3], offsetMatrix[1][3], offsetMatrix[2][3] );
		for( vl::Int i = 0; i < 3; ++i )
			assert( fabs(translate[i]) < epsilon );
	}

	{
		// needs to be a permutation matrix of some sort
		for( vl::Int i = 0; i < 3; ++i )
			for( vl::Int j = 0; j < 3; ++j )
				assert( (i == j) || (fabs(rotMat[i][j]) < epsilon) );
	}

	btVector3 localScale;
	for( vl::Int i = 0; i < 3; ++i )
		localScale[i] = scale[i] + padding;
	return new btBoxShape( localScale );
}
#endif

#ifdef USE_NOVODEX
std::vector< boost::shared_ptr<NxShapeDesc> > BoxObject::getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const
{
	vl::Vec3f scale = this->fullScale();
	boost::shared_ptr<NxBoxShapeDesc> shapeDesc( new NxBoxShapeDesc );
	shapeDesc->dimensions.set( (0.5*scale).Ref() );
	for( size_t i = 0; i < 3; ++i )
		shapeDesc->dimensions[i] += padding + defaultSkinWidth();

	return std::vector< boost::shared_ptr<NxShapeDesc> >(1, shapeDesc);
}
#endif

dMass BoxObject::getOdeMassProps() const
{
	dMass result;
	vl::Vec3f scale = this->fullScale();
	dMassSetBox( &result, this->material()->density(), scale[0], scale[1], scale[2] );
	return result;
}

bool BoxObject::intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const
{
	BoundingBox3d box( vl::Vec3d(-0.5, -0.5, -0.5), vl::Vec3d(0.5, 0.5, 0.5) );
	double tmax;
	bool intersect = box.intersect(ray, t, tmax);
	if( intersect && t >= 0.0 )
		return true;
	else
		return false;
}

bool BoxObject::intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const
{
	for( unsigned int i = 0; i < 2; ++i )
		for( unsigned int j = 0; j < 2; ++j )
			for( unsigned int k = 0; k < 2; ++k )
			{
				vl::Vec3d corner( 
						i == 0 ? -0.5 : 0.5,
						j == 0 ? -0.5 : 0.5,
						k == 0 ? -0.5 : 0.5 );
				if( contains(box, strip( converter.toScreenSpace( corner ) ) ) )
					return true;
			}

	return false;
}

BoxObject::~BoxObject()
{
}

void BoxObject::renderLocalWire( const vl::Vec4f& wireColor, const char* statePtr ) const
{
#ifdef GUI
	glutWireCube(1.0);
#endif
}

void BoxObject::renderLocalFill( float alpha, const char* statePtr ) const
{
#ifdef GUI
    static vl::Vec3f positions[3*24];
    static vl::Vec3f normals[3*24];
    static vl::Vec3f colors[3*24];

    static vl::Vec3f colors6[] = { 
        vl::Vec3f(228.0f, 26.0f,  28.0f ) / 255.0f,
        vl::Vec3f(55.0f,  126.0f, 184.0f) / 255.0f,
        vl::Vec3f(77.0f,  175.0f, 74.0f ) / 255.0f,
        vl::Vec3f(152.0f, 78.0f,  163.0f) / 255.0f,
        vl::Vec3f(255.0f, 127.0f, 0.0f  ) / 255.0f,
        vl::Vec3f(255.0f, 255.0f, 51.0f ) / 255.0f };

    static bool initialized = false;
            
    if( !initialized )
    {
        initialized = true;

        for( size_t iAxis = 0; iAxis < 3; ++iAxis )
        {
            size_t x = (iAxis + 1)%3;
            size_t y = (iAxis + 2)%3;
            size_t z = iAxis;

            // positive version
            vl::Vec3f normal( vl::vl_0 );
            normal[z] = 1.0;
            for( size_t i = 0; i < 4; ++i )
                normals[8*iAxis + i] = normal;

            for( size_t i = 0; i < 4; ++i )
                colors[8*iAxis + i] = colors6[ iAxis ];

            vl::Vec3f position;
            position[x] = 0.5;
            position[y] = 0.5;
            position[z] = 0.5;
            positions[8*iAxis + 0] = position;

            position[x] = -0.5;
            position[y] = 0.5;
            positions[8*iAxis + 1] = position;

            position[x] = -0.5;
            position[y] = -0.5;
            positions[8*iAxis + 2] = position;

            position[x] = 0.5;
            position[y] = -0.5;
            positions[8*iAxis + 3] = position;

            normal[z] = -1.0;
            for( size_t i = 0; i < 4; ++i )
                normals[8*iAxis + 4 + i] = normal;

            for( size_t i = 0; i < 4; ++i )
                colors[8*iAxis + 4 + i] = colors6[ 3 + iAxis ];

            position[x] = 0.5;
            position[y] = -0.5;
            position[z] = -0.5;
            positions[8*iAxis + 4] = position;

            position[x] = -0.5;
            position[y] = -0.5;
            positions[8*iAxis + 5] = position;

            position[x] = -0.5;
            position[y] = 0.5;
            positions[8*iAxis + 6] = position;

            position[x] = 0.5;
            position[y] = 0.5;
            positions[8*iAxis + 7] = position;
        }
    }

	GLEnableClientStateHandler vertexArray(GL_VERTEX_ARRAY);
    glVertexPointer( 3, GL_FLOAT, 0, positions );

    GLEnableClientStateHandler normalArray(GL_NORMAL_ARRAY);
    glNormalPointer( GL_FLOAT, 0, normals );

    GLEnableClientStateHandler colorsArray(GL_COLOR_ARRAY);
    glColorPointer( 3, GL_FLOAT, 0, colors );

    GLEnableHandler colorMaterial( GL_COLOR_MATERIAL );
    glColorMaterial( GL_FRONT, GL_DIFFUSE );
    glDrawArrays( GL_QUADS, 0, 4*6 );

#endif
}

SphereObject::SphereObject( const std::string& name, boost::shared_ptr<Material> material )
	: PhysicsObject( name, material )
{
}

SphereObject* SphereObject::cloneLocal() const
{
	return new SphereObject( *this );
}

BoundingBox3f SphereObject::bounds( const vl::Mat4f& transform ) const
{
	vl::Vec3f center = vl::xform( transform, vl::Vec3f(0.0f, 0.0f, 0.0f) );
	float radius = 
		vl::len( vl::Vec3f( transform[0][0], transform[1][0], transform[2][0] ) );

	BoundingBox3f result( center, center );
	result.expand( radius );
	return result;
}

RibObjectPtr SphereObject::ribObject() const
{
	boost::shared_ptr<RibObjectWithAttributes> sphereObject(
		new RibObjectWithAttributes( "Sphere" ) );
	sphereObject->addAttribute( 1.0 );
	sphereObject->addAttribute( -1.0 );
	sphereObject->addAttribute( 1.0 );

	sphereObject->addAttribute( 360.0 );

	return sphereObject;
}

void SphereObject::toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const
{
	std::string name = prefix + this->name();
	std::string surfaceName = name + "Shape";
	ofs << "createNode nurbsSurface -n \"" 
		<< surfaceName << "\" -p \"" << name << "\";\n";

	std::string makeName = "make" + name;
	ofs << "createNode makeNurbSphere -n \"" << makeName << "\";\n";
	setAttr( ofs, ".ax", vl::Vec3f(0, 1, 0) );
	ofs << "connectAttr \"" << makeName << ".os\" \"" << surfaceName << ".cr\";\n";
	ofs << "connectAttr \"" << surfaceName << ".iog\" \"" << materialName << ".dsm\" -na;\n";
}

boost::shared_ptr<MechModelDeclaredObject> SphereObject::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "sphere", this->name() ) );

	return result;
}

ODEGeomResult SphereObject::getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const
{
	assert( centerOfMass == vl::vl_0 );
	vl::Vec3f scale = this->fullScale();
	return ODEGeomResult( dCreateSphere( space, scale[0] + padding) );
}

vl::Vec3f SphereObject::furthestPoint( const vl::Vec3f& point ) const
{
	vl::Vec3f scale = this->fullScale();
	return vl::Vec3f( 0, scale[0], 0 );
}

#ifdef USE_NEWTON
NewtonCollision* SphereObject::getNewtonCollision( 
	const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const
{
	vl::Vec3f scale = this->fullScale();
	return NewtonCreateSphere( newtonWorld, scale[0] + padding, scale[1] + padding, scale[2] + padding, offsetMatrix.Ref() );
}
#endif

#ifdef USE_NOVODEX
std::vector< boost::shared_ptr<NxShapeDesc> > SphereObject::getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const
{
	vl::Vec3f scale = this->fullScale();
	boost::shared_ptr<NxSphereShapeDesc> shapeDesc( new NxSphereShapeDesc );
	shapeDesc->radius = scale[0] + padding + defaultSkinWidth();
	return std::vector< boost::shared_ptr<NxShapeDesc> > (1, shapeDesc);
}
#endif

#ifdef USE_BULLET
btCollisionShape* SphereObject::bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const
{
	vl::Vec3f scale = this->fullScale();
	return new btSphereShape( scale[0] + padding );
}
#endif

dMass SphereObject::getOdeMassProps() const
{
	dMass result;
	vl::Vec3f scale = this->fullScale();
	dMassSetSphere( &result, this->material()->density(), scale[0] );
	return result;
}

bool intersectSphere( const Ray3d& ray, double& t, const vl::Vec3d& center, double radius )
{
	vl::Vec3d dst = ray.getPosition() - center;
	double A = dot( ray.getDirection(), ray.getDirection() );
	double B = 2.0*dot( dst, ray.getDirection() );
	double C = dot( dst, dst ) - (radius*radius);
	double D = B*B - 4.0*A*C;

	if( D > 0.0f )
	{
		double t0 = (-B - sqrt(D)) / (2*A);
		double t1 = (-B + sqrt(D)) / (2*A);
		if( t0 >= 0 )
			t = t0;
		else if( t1 >= 0 )
			t = t1;
		else
			return false;

		return true;
	}
	else
	{
		return false;
	}
}

// Checks whether the ray intersects the circle with the specified origin, normal, and radius
bool intersectCircle( const Ray3d& ray, 
					 double& t, 
					 const vl::Vec3d& origin, 
					 const vl::Vec3d& normal, 
					 double radius )
{
	double local_t = dot( origin - ray.getPosition(), normal ) 
		/ dot( ray.getDirection(), normal );
	vl::Vec3d pos = ray.at( local_t );
	if( vl::len( pos - origin ) < radius )
	{
		t = local_t;
		return true;
	}

	return false;
}

struct BoundsForVec3d
{
	BoundingBox3d operator()( const vl::Vec3d& value ) const
	{
		return BoundingBox3d(value, value);
	}
};

class PointIntersectionHandler
{
public:
	PointIntersectionHandler( const BoundingBox2d& box, const ScreenSpaceConverter& converter )
		: box_(box), converter_(converter)
	{
	}

	bool operator()( const vl::Vec3d& pos )
	{
		return contains( box_, strip(converter_.toScreenSpace( pos )) );
	}
	
private:
	const BoundingBox2d box_;
	const ScreenSpaceConverter& converter_;
};

bool intersectSphere( const BoundingBox2d& box, const ScreenSpaceConverter& converter, 
					 const vl::Mat4d& transform )
{
	// warning not threadsafe
	static boost::shared_ptr< AABBTree<vl::Vec3d> > tree;
	if( !tree )
	{
		std::vector<vl::Vec3d> points;
		for( int i = 0; i <= sphereSubdivisions; ++i )
		{
			for( int j = 0; j <= sphereSubdivisions; ++j )
			{
				double theta = 2.0 * M_PI * 
					(boost::numeric_cast<double>(i) / boost::numeric_cast<double>(sphereSubdivisions));
				double phi = M_PI * 
					(boost::numeric_cast<double>(j) / boost::numeric_cast<double>(sphereSubdivisions));

				points.push_back( 
					vl::Vec3d( cos(theta) * sin(phi),
						sin(theta) * sin(phi),
						cos(phi) ) );
			}
		}

		tree.reset( new AABBTree<vl::Vec3d>(points, BoundsForVec3d()) );
	}

	WrappedScreenSpaceConverter wrappedConverter( converter, transform );

	PointIntersectionHandler handler( box, wrappedConverter );
	return tree->intersect( box, wrappedConverter, handler );
}

bool SphereObject::intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const
{
	return intersectSphere( ray, t );
}

bool SphereObject::intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const
{
	return intersectSphere( box, converter, vl::vl_1 );
}

SphereObject::~SphereObject()
{
}

void SphereObject::renderLocalWire( const vl::Vec4f& wireColor, const char* statePtr ) const
{
#ifdef GUI
	GLUQuadricWrapper wireQuadric;
	gluQuadricDrawStyle( wireQuadric, GLU_LINE );

	gluSphere( wireQuadric, 1.0, sphereSubdivisions, sphereSubdivisions );
#endif
}

void SphereObject::renderLocalFill( float alpha, const char* statePtr ) const
{
#ifdef GUI
    glMatrixMode(GL_TEXTURE);
    glPushMatrix();
    glScaled( 1.0, 0.5, 1.0 );

	GLUQuadricWrapper fillQuadric;
	gluQuadricDrawStyle( fillQuadric, GLU_FILL );
    gluQuadricTexture( fillQuadric, GL_TRUE );

    const GLTexture& texture = checkerboardTextureMap();
    GLBindTextureHandler bindTexture(texture);

	gluSphere( fillQuadric, 1.0, sphereSubdivisions, sphereSubdivisions );
    glPopMatrix();
#endif
}

ConeObject::ConeObject( const std::string& name, boost::shared_ptr<Material> material )
	: PhysicsObject( name, material ), height_(1.0f)
{
}

ConeObject::~ConeObject()
{
}

bool ConeObject::equalsLocal( const SceneGraphElement& other ) const
{
	const ConeObject* cone = dynamic_cast<const ConeObject*>( &other );
	assert( cone != 0 );
	if( !withinEpsilon( this->height_, cone->height_ ) )
		return false;

	return PhysicsObject::equalsLocal( other );
}

const char* HEIGHT_STRING = "Height";

std::vector<std::string> ConeObject::floatAttributes() const
{
	std::vector<std::string> result = PhysicsObject::floatAttributes();
	result.reserve( result.size()+1 );
	result.push_back( HEIGHT_STRING );
	return result;
}

void ConeObject::setAttribute( const std::string& attributeName, float value )
{
	if( attributeName == HEIGHT_STRING )
		this->height_ = value;
	else
		PhysicsObject::setAttribute(attributeName, value);
}

float ConeObject::getAttribute( const std::string& attributeName ) const
{
	if( attributeName == HEIGHT_STRING )
		return this->height_;
	else
		return PhysicsObject::getAttribute(attributeName);
}

#ifdef USE_NOVODEX
// novodex doesn't have a cone primitive
std::vector< boost::shared_ptr<NxShapeDesc> > ConeObject::getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const
{
	std::vector< boost::shared_ptr<NxShapeDesc> > result;
	return result;
}
#endif

#ifdef USE_BULLET
btCollisionShape* ConeObject::bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const
{
	vl::Vec3f scale = this->fullScale();

	btConeShape* cone = new btConeShape( scale[0], height_*scale[0] );

	vl::Mat4f transform( vl::vl_1 );
	transform[1][3] = 0.5 * scale[0] * height_;
	btCompoundShape* result = new btCompoundShape();
	result->addChildShape( toBtTransform( transform ), cone );

	return result;
}
#endif

ConeObject* ConeObject::cloneLocal() const
{
	return new ConeObject(*this);
}

BoundingBox3f ConeObject::bounds( const vl::Mat4f& transform ) const
{
	// @todo make this a tighter bound
	BoundingBox3f result;
	result.expand( vl::xform(transform, vl::Vec3f(0.0f, 1.0f, 0.0f) ) );
	vl::Vec3f fullScale = this->fullScale();
	result.expand( std::max( std::max(fullScale[0], fullScale[1]), fullScale[2] ) );
	return result;
}

RibObjectPtr ConeObject::ribObject() const
{
	return ribObjectWithAttribute("ObjectInstance", 4);
}

boost::shared_ptr<MechModelDeclaredObject> ConeObject::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "cone", this->name() ) );
	result->add( mechModelValuePair("height", this->height_) );

	return result;
}

void ConeObject::setHeight( float height )
{
	height_ = height;
}

float ConeObject::height() const
{
	return height_;
}


#ifdef GUI
void drawCone()
{
	// warning: not threadsafe
	{
		static std::vector<vl::Vec3f> points;
		static std::vector<vl::Vec3f> normals;

		if( points.empty() )
		{
			points.reserve( 2*(sphereSubdivisions+1) );
			normals.reserve( 2*(sphereSubdivisions+1) );

			for( size_t i = 0; i <= sphereSubdivisions; ++i )
			{
				float theta = 2.0f*M_PI*(float(i) / float(sphereSubdivisions));
				float midTheta = 2.0f*M_PI*( (float(i) + 0.5f) / float(sphereSubdivisions));
				normals.push_back( 
					vl::norm( vl::Vec3f(cos( theta ), 1.0f, -sin( theta ) ) ) );
				normals.push_back( 
					vl::norm( vl::Vec3f(cos( midTheta ), 1.0f, -sin( midTheta ) ) ) );
				points.push_back( vl::Vec3f( 0.0f, 1.0f, 0.0f ) );
				points.push_back( 
					vl::Vec3f( sin( theta ), 0.0f, cos( theta ) ) );
			}
		}
		

		// draw the cone part
		{
			GLEnableClientStateHandler vertexArray(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, 3*sizeof(float), &points[0]);

			GLEnableClientStateHandler normalArray(GL_NORMAL_ARRAY);
			glNormalPointer(GL_FLOAT, 3*sizeof(float), &normals[0]);

			glDrawArrays( GL_TRIANGLE_STRIP, 0, points.size() );
		}
	}

	// now draw the base
	{
		static std::vector<vl::Vec3f> points;
		static std::vector<vl::Vec3f> normals;

		if( points.empty() )
		{
			points.reserve( 2*(sphereSubdivisions+1) + 1 );
			normals.reserve( 2*(sphereSubdivisions+1) + 1 );

			points.push_back( vl::Vec3f(0.0f, 0.0f, 0.0f) );
			normals.push_back( vl::Vec3f(0.0f, -1.0f, 0.0f) );

			for( size_t i = 0; i <= sphereSubdivisions; ++i )
			{
				float theta = 2.0f*M_PI*(float(sphereSubdivisions - i) / float(sphereSubdivisions));
				normals.push_back( vl::Vec3f(0.0f, -1.0f, 0.0f) );
				points.push_back( 
					vl::Vec3f( sin( theta ), 0.0f, cos( theta ) ) );
			}
		}

		{
			GLEnableClientStateHandler vertexArray(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, 3*sizeof(float), &points[0]);

			GLEnableClientStateHandler normalArray(GL_NORMAL_ARRAY);
			glNormalPointer(GL_FLOAT, 3*sizeof(float), &normals[0]);

			glDrawArrays( GL_TRIANGLE_FAN, 0, points.size() );
		}
	}
}
#endif

void ConeObject::renderLocalFill(float alpha, const char* statePtr) const
{
	GLMatrixStackHandler pushMatrix;
	glScalef( 1.0f, height_, 1.0f );

	glEnable( GL_LIGHTING );
	glPolygonMode( GL_FRONT, GL_FILL );
	glShadeModel( GL_SMOOTH );

	drawCone();
}

void ConeObject::renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const
{
	GLMatrixStackHandler pushMatrix;
	glScalef( 1.0f, height_, 1.0f );

	glDisable( GL_LIGHTING );
	glPolygonMode( GL_FRONT, GL_LINE );
	glColor4fv( wireframeColor.Ref() );
	glLineWidth( 1.0 );

	drawCone();
}


// check whether the box intersects the cylinder of unit length and radius along the
//   z axis (plus an optional transform)
bool intersectCone( const BoundingBox2d& box, const ScreenSpaceConverter& converter, 
					   const vl::Mat4d& transform = vl::Mat4d(vl::vl_1) )
{
	// warning not threadsafe
	static boost::shared_ptr< AABBTree<vl::Vec3d> > tree;
	if( !tree )
	{
		std::vector<vl::Vec3d> points;
		points.reserve( sphereSubdivisions * (sphereSubdivisions-1) + 1 );
		for( int i = 0; i <= sphereSubdivisions; ++i )
		{
			for( int j = 0; j < sphereSubdivisions; ++j )
			{
				double theta = 2.0 * M_PI * 
					(boost::numeric_cast<double>(i) / boost::numeric_cast<double>(sphereSubdivisions));
				double y = 0.0f + 
					(boost::numeric_cast<double>(j) / boost::numeric_cast<double>(sphereSubdivisions));

				points.push_back( 
					vl::Vec3d( cos(theta)*(1.0 - y),
						y,
						sin(theta)*(1.0 - y) ) );
			}
		}

		points.push_back( vl::Vec3d(0.0, 1.0, 0.0) );

		tree.reset( new AABBTree<vl::Vec3d>(points, BoundsForVec3d()) );
	}

	WrappedScreenSpaceConverter wrappedConverter( converter, transform );

	PointIntersectionHandler handler( box, wrappedConverter );
	return tree->intersect( box, wrappedConverter, handler );
}

// intersect with the cylinder with specified length along the z axis
bool intersectCone( const Ray3d& ray, double& t, double height )
{
	double ts[2];
	int num_t = 0;

	{
		double cone_ts[2];
		int num_cone_ts = 0;

		const double h = height;
		const double h2 = h*h;

		const vl::Vec3d p = ray.getPosition();
		const vl::Vec3d d = ray.getDirection();

		const double A = d[0]*d[0] + d[2]*d[2] - d[1]*d[1]/h2;
		const double B = 2.0*p[0]*d[0] + 2.0*p[2]*d[2] - 2.0*p[1]*d[1]/h2 + 2.0*d[1]/h;
		const double C = p[0]*p[0] + p[2]*p[2] - p[1]*p[1]/h2 + 2.0*p[1]/h - 1.0;
		double det = B*B - 4*A*C;
		if( det >= 0.0 )
		{
			double t0 = (-B - sqrt(det))/(2.0*A);
			double t1 = (-B + sqrt(det))/(2.0*A);

			if( t0 >= 0.0 )
				cone_ts[ num_cone_ts++ ] = t0;
			if( t1 >= 0.0 )
				cone_ts[ num_cone_ts++ ] = t1;

			for( int i = 0; i < num_cone_ts; ++i )
			{
				double t = cone_ts[i];
				vl::Vec3d pos = ray.at( t );
				if( pos[1] > 0.0 && pos[1] < h )
					ts[ num_t++ ] = t;
			}
		}
	}

	if( num_t == 0 )
		return false;

	t = ts[0];
	for( int i = 1; i < num_t; ++i )
		t = std::min( t, ts[i] );

	return true;
}

bool ConeObject::intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const
{
	double ts[2];
	int num_t = 0;

	{
		double t0;
		if( intersectCone(ray, t0, height_) )
			ts[ num_t++ ] = t0;
	}

	/*
	{
		double t0;
		if( intersectCircle(ray, t0, vl::Vec3d(0, 0, 0), vl::Vec3d(0, 0, -1), 1.0) )
			ts[ num_t++ ] = t0;
	}
	*/

	if( num_t == 0 )
		return false;

	t = ts[0];
	for( int i = 1; i < num_t; ++i )
		t = std::min( t, ts[i] );

	return true;
}



bool ConeObject::intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const
{
	vl::Mat4d transform( vl::vl_1 );
	transform[1][1] = height_;
	return intersectCone(box, converter, transform);
}

ODEGeomResult ConeObject::getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const
{
	// ODE doesn't support a cone primitive
	ODEGeomResult result;
	return result;
}

dMass ConeObject::getOdeMassProps() const
{
	// need to implement
	// @otod fix this to be correct
	dMass result;
	vl::Vec3f scale = this->fullScale();
	dMassSetSphere( &result, this->material()->density(), scale[0] );
	return result;
}

CappedCylinderObject::CappedCylinderObject( const std::string& name, boost::shared_ptr<Material> material )
	: PhysicsObject( name, material ), lengthRatio_(1.0f)
{
}

CappedCylinderObject::~CappedCylinderObject()
{
}

CappedCylinderObject* CappedCylinderObject::cloneLocal() const
{
	return new CappedCylinderObject( *this );
}

bool CappedCylinderObject::equalsLocal( const SceneGraphElement& other ) const
{
	const CappedCylinderObject* cyl = dynamic_cast<const CappedCylinderObject*>( &other );
	assert( cyl != 0 );
	if( !withinEpsilon( this->lengthRatio_, cyl->lengthRatio_ ) )
		return false;

	return PhysicsObject::equalsLocal( other );
}

void CappedCylinderObject::toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const
{
	std::string name = prefix + this->name();
	{
		std::string cylinderName = name + "Cylinder";
		ofs << "createNode transform -n \"" << cylinderName << "\" -p \"" << name << "\";\n";
		std::string cylinderSurfaceName = cylinderName + "Shape";
		ofs << "createNode nurbsSurface -n \"" 
			<< cylinderSurfaceName << "\" -p \"" << cylinderName << "\";\n";

		std::string makeName = "make" + cylinderName;
		ofs << "createNode makeNurbCylinder -n \"" << makeName << "\";\n";
		setAttr( ofs, ".ax", vl::Vec3f(0, 0, 1) );
		ofs << "connectAttr \"" << makeName << ".os\" \"" << cylinderSurfaceName << ".cr\";\n";
		setAttr( ofs, ".heightRatio", this->lengthRatio_ );
		ofs << "connectAttr \"" << cylinderSurfaceName << ".iog\" \"" << materialName << ".dsm\" -na;\n";
	}
	{
		std::string capName = name + "Cap1";
		ofs << "createNode transform -n \"" << capName << "\" -p \"" << name << "\";\n";
		setAttr( ofs, ".t", vl::Vec3f(0, 0, 0.5*lengthRatio_) );
		std::string capSurfaceName = capName + "Shape";
		ofs << "createNode nurbsSurface -n \"" 
			<< capSurfaceName << "\" -p \"" << capName << "\";\n";

		std::string makeName = "make" + capName;
		ofs << "createNode makeNurbSphere -n \"" << makeName << "\";\n";
		setAttr( ofs, ".ax", vl::Vec3f(0, 1, 0) );
		setAttr( ofs, ".startSweep", 0.0f );
		setAttr( ofs, ".endSweep", 180.0f );
		ofs << "connectAttr \"" << makeName << ".os\" \"" << capSurfaceName << ".cr\";\n";
		ofs << "connectAttr \"" << capSurfaceName << ".iog\" \"" << materialName << ".dsm\" -na;\n";
	}
	{
		std::string capName = name + "Cap2";
		ofs << "createNode transform -n \"" << capName << "\" -p \"" << name << "\";\n";
		setAttr( ofs, ".t", vl::Vec3f(0, 0, -0.5*lengthRatio_) );
		std::string capSurfaceName = capName + "Shape";
		ofs << "createNode nurbsSurface -n \"" 
			<< capSurfaceName << "\" -p \"" << capName << "\";\n";

		std::string makeName = "make" + capName;
		ofs << "createNode makeNurbSphere -n \"" << makeName << "\";\n";
		setAttr( ofs, ".ax", vl::Vec3f(0, -1, 0) );
		setAttr( ofs, ".startSweep", 0.0f );
		setAttr( ofs, ".endSweep", 180.0f );
		ofs << "connectAttr \"" << makeName << ".os\" \"" << capSurfaceName << ".cr\";\n";
		ofs << "connectAttr \"" << capSurfaceName << ".iog\" \"" << materialName << ".dsm\" -na;\n";
	}
}

BoundingBox3f CappedCylinderObject::bounds(const vl::Mat4f& transform) const
{
	vl::Vec3f top = vl::xform( transform,
		vl::Vec3f(1, 1, 1.0f + 0.5f*lengthRatio_) );
	vl::Vec3f bottom = vl::xform( transform, 
		vl::Vec3f(-1.0f, -1.0f, -1.0f - 0.5f*lengthRatio_) );

	float radius = 
		vl::len( vl::Vec3f( transform[0][0], transform[1][0], transform[2][0] ) );
	BoundingBox3f result;
	result.expand( top );
	result.expand( bottom );
	result.expand( radius );

	return result;
}

boost::shared_ptr<MechModelDeclaredObject> CappedCylinderObject::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "capsule", this->name() ) );
	result->add( mechModelValuePair("lengthRatio", this->lengthRatio_) );

	return result;
}

const char* LENGTH_RATIO_STRING = "Length ratio";

std::vector<std::string> CappedCylinderObject::floatAttributes() const
{
	std::vector<std::string> result = PhysicsObject::floatAttributes();
	result.reserve( result.size()+1 );
	result.push_back( LENGTH_RATIO_STRING );
	return result;
}

void CappedCylinderObject::setAttribute( const std::string& attributeName, float value )
{
	if( attributeName == LENGTH_RATIO_STRING )
		this->lengthRatio_ = value;
	else
		PhysicsObject::setAttribute(attributeName, value);
}

float CappedCylinderObject::getAttribute( const std::string& attributeName ) const
{
	if( attributeName == LENGTH_RATIO_STRING )
		return this->lengthRatio_;
	else
		return PhysicsObject::getAttribute(attributeName);
}


#ifdef GUI
void drawCappedCylinder( GLUQuadricWrapper& quadric, float lengthRatio, size_t numSubdivisions )
{
	glTranslatef( 0.0f, 0.0f, -0.5*lengthRatio );
	gluCylinder( quadric, 1.0, 1.0, lengthRatio, numSubdivisions, 1 );
	glTranslatef( 0.0f, 0.0f, 0.5*lengthRatio );

	glTranslatef( 0.0f, 0.0f, -0.5*lengthRatio );
	gluSphere( quadric, 1.0, numSubdivisions, numSubdivisions );
	glTranslatef( 0.0f, 0.0f, 0.5*lengthRatio );

	glTranslatef( 0.0f, 0.0f, 0.5*lengthRatio );
	gluSphere( quadric, 1.0, numSubdivisions, numSubdivisions );
	glTranslatef( 0.0f, 0.0f, -0.5*lengthRatio );
}
#endif


void CappedCylinderObject::renderLocalFill(float alpha, const char* statePtr) const
{
#ifdef GUI
	GLUQuadricWrapper quadric;
	gluQuadricDrawStyle( quadric, GLU_FILL );
	drawCappedCylinder( quadric, this->lengthRatio_, sphereSubdivisions );
#endif
}

void CappedCylinderObject::renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const
{
#ifdef GUI
	GLUQuadricWrapper quadric;
	gluQuadricDrawStyle( quadric, GLU_LINE );
	drawCappedCylinder( quadric, this->lengthRatio_, sphereSubdivisions/2 );
#endif
}

// intersect with the cylinder with specified length along the z axis
bool intersectCylinder( const Ray3d& ray, double& t, const double length )
{
	double ts[2];
	int num_t = 0;

	{
		double cyl_ts[2];
		int num_cyl_ts = 0;

		const vl::Vec3d p = ray.getPosition();
		const vl::Vec3d d = ray.getDirection();

		double B = 2*(p[0]*d[0] + p[1]*d[1]);
		double A = d[0]*d[0] + d[1]*d[1];
		double C = p[0]*p[0] + p[1]*p[1] - 1.0;
		double det = B*B - 4*A*C;
		if( det >= 0.0 )
		{
			double t0 = (-B - sqrt(det))/(2.0*A);
			double t1 = (-B + sqrt(det))/(2.0*A);

			if( t0 >= 0.0 )
				cyl_ts[ num_cyl_ts++ ] = t0;
			if( t1 >= 0.0 )
				cyl_ts[ num_cyl_ts++ ] = t1;

			for( int i = 0; i < num_cyl_ts; ++i )
			{
				double t = cyl_ts[i];
				vl::Vec3d pos = ray.at( t );
				if( pos[2] > -0.5*length && pos[2] < 0.5*length )
					ts[ num_t++ ] = t;
			}
		}
	}

	if( num_t == 0 )
		return false;

	t = ts[0];
	for( int i = 1; i < num_t; ++i )
		t = std::min( t, ts[i] );

	return true;
}

bool CappedCylinderObject::intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const
{
	double ts[4];
	int num_t = 0;

	// first the caps
	{
		double t0;
		if( intersectSphere( ray, t0, vl::Vec3d(0.0, 0.0, 0.5*this->lengthRatio_) ) )
			ts[ num_t++ ] = t0;
	}

	{
		double t0;
		if( intersectSphere( ray, t0, vl::Vec3d(0.0, 0.0, -0.5*this->lengthRatio_) ) )
			ts[ num_t++ ] = t0;
	}

	// now the cylinder
	{
		double t0;
		if( intersectCylinder( ray, t0, this->lengthRatio_ ) )
			ts[ num_t++ ] = t0;
	}

	if( num_t == 0 )
		return false;

	t = ts[0];
	for( int i = 1; i < num_t; ++i )
		t = std::min( t, ts[i] );

	return true;
}

// check whether the box intersects the cylinder of unit length and radius along the
//   z axis (plus an optional transform)
bool intersectCylinder( const BoundingBox2d& box, const ScreenSpaceConverter& converter, 
					   const vl::Mat4d& transform = vl::Mat4d(vl::vl_1) )
{
	// warning not threadsafe
	static boost::shared_ptr< AABBTree<vl::Vec3d> > tree;
	if( !tree )
	{
		std::vector<vl::Vec3d> points;
		for( int i = 0; i <= sphereSubdivisions; ++i )
		{
			for( int j = 0; j <= sphereSubdivisions; ++j )
			{
				double theta = 2.0 * M_PI * 
					(boost::numeric_cast<double>(i) / boost::numeric_cast<double>(sphereSubdivisions));
				double z = -0.5 + 
					(boost::numeric_cast<double>(j) / boost::numeric_cast<double>(sphereSubdivisions));

				points.push_back( 
					vl::Vec3d( cos(theta),
						sin(theta),
						z ) );
			}
		}

		tree.reset( new AABBTree<vl::Vec3d>(points, BoundsForVec3d()) );
	}

	WrappedScreenSpaceConverter wrappedConverter( converter, transform );

	PointIntersectionHandler handler( box, wrappedConverter );
	return tree->intersect( box, wrappedConverter, handler );
}

bool CappedCylinderObject::intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const
{
	{
		vl::Mat4d cap1_xform( vl::vl_1 );
		cap1_xform[2][3] = -this->lengthRatio_;
		if( intersectSphere( box, converter, cap1_xform ) )
			return true;
	}

	{
		vl::Mat4d cap2_xform( vl::vl_1 );
		cap2_xform[2][3] = this->lengthRatio_;
		if( intersectSphere( box, converter, cap2_xform ) )
			return true;
	}

	{
		vl::Mat4d transform( vl::vl_1 );
		transform[2][2] = this->lengthRatio_;
		return intersectCylinder( box, converter, transform );
	}
}

ODEGeomResult CappedCylinderObject::getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const
{
	assert( centerOfMass == vl::vl_0 );
	vl::Vec3f scale = this->fullScale();
	return ODEGeomResult( dCreateCCylinder( space, scale[0] + padding, scale[0]*this->lengthRatio_ + padding ) );
}

#ifdef USE_NOVODEX
std::vector< boost::shared_ptr<NxShapeDesc> > CappedCylinderObject::getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const
{
	vl::Mat4f transform = this->transform();
	vl::Mat3f upperLeftBlock = toMat3f( transform );

	vl::Vec3f scale = this->fullScale();
	boost::shared_ptr<NxCapsuleShapeDesc> shapeDesc( new NxCapsuleShapeDesc );
	shapeDesc->radius = 1.0f*scale[0] + padding + defaultSkinWidth();
	shapeDesc->height = 1.0f*scale[0]*lengthRatio_ + padding + defaultSkinWidth();

	NxMat33 rotMat;
	rotMat.rotX(0.5*M_PI);

	shapeDesc->localPose = NxMat34( rotMat, NxVec3(0.0f, 0.0f, 0.0f) );

	return std::vector< boost::shared_ptr<NxShapeDesc> > (1, shapeDesc);
}
#endif

#ifdef USE_NEWTON
NewtonCollision* CappedCylinderObject::getNewtonCollision( 
	const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const
{
	vl::Vec3f scale = this->fullScale();
	dFloatMat4 rotateMat;
	rotateMat.MakeHRot( dFloatVec3(0, 1, 0), M_PI/2.0 );
	dFloatMat4 localMat = trans(rotateMat) * offsetMatrix;
	return NewtonCreateCapsule( newtonWorld, scale[0] + padding, scale[0]*(this->lengthRatio_ + 2.0) + padding, localMat.Ref() );
}
#endif

#ifdef USE_BULLET
btCollisionShape* CappedCylinderObject::bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const
{
	vl::Vec3f scale = this->fullScale();

	// These don't really matter since we're using ODE to compute inertia:
	btVector3 inertiaHalfExtents( 
		scale[0] + padding, scale[1] + padding, scale[2]*(this->lengthRatio_ + 2.0) + padding );

	vl::Mat4f scaleMatrix = vl::vl_1;
	for( vl::Int i = 0; i < 3; ++i )
		scaleMatrix[i][i] = scale[0];
	vl::Mat4f fullTransform = offsetMatrix * scaleMatrix;

	vl::Vec3f pos1 = vl::Vec3f(0, 0, -0.5*this->lengthRatio_);
	vl::Vec3f pos2 = vl::Vec3f(0, 0, 0.5*this->lengthRatio_);

	btVector3 positions[2] = { toBtVector3(pos1), toBtVector3(pos2) };
	btScalar radii[2] = { scale[0] + padding, scale[0] + padding };
	btMultiSphereShape* result = new btMultiSphereShape( 
		inertiaHalfExtents,
		&positions[0],
		&radii[0], 
		2 );
	return result;
}
#endif

dMass CappedCylinderObject::getOdeMassProps() const
{
	dMass result;
	vl::Vec3f scale = this->fullScale();
	dMassSetCappedCylinder( &result, this->material()->density(), 3, scale[0], lengthRatio_*scale[0] );
	return result;
}

RibObjectPtr CappedCylinderObject::ribObject() const
{
	return ribObjectWithAttribute("ObjectInstance", 3);
}


CylinderObject::CylinderObject( const std::string& name, boost::shared_ptr<Material> material )
	: PhysicsObject( name, material ), lengthRatio_(1.0f)
{
}

CylinderObject::~CylinderObject()
{
}

bool CylinderObject::equalsLocal( const SceneGraphElement& other ) const
{
	const CylinderObject* cyl = dynamic_cast<const CylinderObject*>( &other );
	assert( cyl != 0 );
	if( !withinEpsilon( this->lengthRatio_, cyl->lengthRatio_ ) )
		return false;

	return PhysicsObject::equalsLocal( other );
}


CylinderObject* CylinderObject::cloneLocal() const
{
	return new CylinderObject( *this );
}


void CylinderObject::toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const
{
	// @todo put caps on it
	std::string name = prefix + this->name();
	{
		std::string cylinderName = name + "Cylinder";
		ofs << "createNode transform -n \"" << cylinderName << "\" -p \"" << name << "\";\n";
		std::string cylinderSurfaceName = cylinderName + "Shape";
		ofs << "createNode nurbsSurface -n \"" 
			<< cylinderSurfaceName << "\" -p \"" << cylinderName << "\";\n";

		std::string makeName = "make" + cylinderName;
		ofs << "createNode makeNurbCylinder -n \"" << makeName << "\";\n";
		setAttr( ofs, ".ax", vl::Vec3f(0, 0, 1) );
		ofs << "connectAttr \"" << makeName << ".os\" \"" << cylinderSurfaceName << ".cr\";\n";
		setAttr( ofs, ".heightRatio", this->lengthRatio_ );
		ofs << "connectAttr \"" << cylinderSurfaceName << ".iog\" \"" << materialName << ".dsm\" -na;\n";
	}
}

BoundingBox3f CylinderObject::bounds(const vl::Mat4f& transform) const
{
	vl::Vec3f top = vl::xform( transform,
		vl::Vec3f(1, 1, 1.0f + 0.5f*lengthRatio_) );
	vl::Vec3f bottom = vl::xform( transform, 
		vl::Vec3f(-1.0f, -1.0f, -1.0f - 0.5f*lengthRatio_) );

	float radius = 
		vl::len( vl::Vec3f( transform[0][0], transform[1][0], transform[2][0] ) );
	BoundingBox3f result;
	result.expand( top );
	result.expand( bottom );
	result.expand( radius );

	return result;
}

boost::shared_ptr<MechModelDeclaredObject> CylinderObject::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "cylinder", this->name() ) );
	if( lengthRatio_ != 1.0f )
		result->add( mechModelValuePair("lengthRatio", this->lengthRatio_) );

	return result;
}

std::vector<std::string> CylinderObject::floatAttributes() const
{
	std::vector<std::string> result = PhysicsObject::floatAttributes();
	result.reserve( result.size()+1 );
	result.push_back( LENGTH_RATIO_STRING );
	return result;
}

void CylinderObject::setAttribute( const std::string& attributeName, float value )
{
	if( attributeName == LENGTH_RATIO_STRING )
		this->lengthRatio_ = value;
	else
		PhysicsObject::setAttribute(attributeName, value);
}

float CylinderObject::getAttribute( const std::string& attributeName ) const
{
	if( attributeName == LENGTH_RATIO_STRING )
		return this->lengthRatio_;
	else
		return PhysicsObject::getAttribute(attributeName);
}


#ifdef GUI
void drawCylinder( GLUQuadricWrapper& quadric, float lengthRatio, size_t numSubdivisions )
{
	GLMatrixStackHandler pushMatrix;
	glTranslatef( 0.0f, 0.0f, -lengthRatio );
	gluCylinder( quadric, 1.0, 1.0, 2.0*lengthRatio, numSubdivisions, 1 );
	glTranslatef( 0.0f, 0.0f, lengthRatio );

	glTranslatef( 0.0f, 0.0f, lengthRatio );
	gluDisk( quadric, 0.0, 1.0, numSubdivisions, 1 );
	glTranslatef( 0.0f, 0.0f, -lengthRatio );

	glRotatef( 180, 0, 1, 0 );
	glTranslatef( 0.0f, 0.0f, lengthRatio );
	gluDisk( quadric, 0.0, 1.0, numSubdivisions, 1 );
}
#endif


void CylinderObject::renderLocalFill(float alpha, const char* statePtr) const
{
#ifdef GUI
	GLUQuadricWrapper quadric;
	gluQuadricDrawStyle( quadric, GLU_FILL );
	drawCylinder( quadric, this->lengthRatio_, sphereSubdivisions );
#endif
}

void CylinderObject::renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const
{
#ifdef GUI
	GLUQuadricWrapper quadric;
	gluQuadricDrawStyle( quadric, GLU_LINE );
	drawCylinder( quadric, this->lengthRatio_, sphereSubdivisions/2 );
#endif
}

bool CylinderObject::intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const
{
	double ts[3];
	int num_t = 0;

	{
		double t0;
		if( intersectCylinder(ray, t0, 2.0*this->lengthRatio_) )
			ts[ num_t++ ] = t0;
	}

	{
		double t0;
		if( intersectCircle(ray, t0, vl::Vec3d(0, 0, -1), vl::Vec3d(0, 0, -1), 1.0) )
			ts[ num_t++ ] = t0;
	}

	{
		double t0;
		if( intersectCircle(ray, t0, vl::Vec3d(0, 0, 1), vl::Vec3d(0, 0, 1), 1.0) )
			ts[ num_t++ ] = t0;
	}

	if( num_t == 0 )
		return false;

	t = ts[0];
	for( int i = 1; i < num_t; ++i )
		t = std::min( t, ts[i] );

	return true;
}

bool CylinderObject::intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const
{
	vl::Mat4d transform( vl::vl_1 );
	transform[2][2] = 2.0 * this->lengthRatio_;
	return intersectCylinder( box, converter, transform );
}

ODEGeomResult CylinderObject::getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const
{
	assert( centerOfMass == vl::vl_0 );
	vl::Vec3f scale = this->fullScale();
#if ODE_VERSION > 5
	return ODEGeomResult( dCreateCylinder( space, scale[0] + padding, scale[0]*this->lengthRatio_ + padding ) );
#else
	return ODEGeomResult( dCreateCCylinder( space, scale[0] + padding, scale[0]*this->lengthRatio_ + padding ) );
#endif
}

#ifdef USE_NEWTON
NewtonCollision* CylinderObject::getNewtonCollision( 
	const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const
{
	vl::Vec3f scale = this->fullScale();

	dFloatMat4 rotateMat;
	rotateMat.MakeHRot( dFloatVec3(0, 1, 0), M_PI/2.0 );
	dFloatMat4 localMat = trans(rotateMat) * offsetMatrix;

	return NewtonCreateCylinder( newtonWorld, 
		scale[0] + padding, 
		2.0*scale[0]*this->lengthRatio_ + padding, 
		localMat.Ref() );
}
#endif

#ifdef USE_NOVODEX
// novodex doesn't have a cylinder primitive
std::vector< boost::shared_ptr<NxShapeDesc> > CylinderObject::getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const
{
	std::vector< boost::shared_ptr<NxShapeDesc> > result;
	return result;
}
#endif

#ifdef USE_BULLET
btCollisionShape* CylinderObject::bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const
{
	vl::Vec3f scale = this->fullScale();
	vl::Mat3f rotMat = toMat3f( offsetMatrix );

	const float epsilon = 1e-4;
	{
		// dont know how to apply a translation offset
		vl::Vec3f translate( offsetMatrix[0][3], offsetMatrix[1][3], offsetMatrix[2][3] );
		for( vl::Int i = 0; i < 3; ++i )
			assert( fabs(translate[i]) < epsilon );
	}

	{
		// needs to be a permutation matrix of some sort
		for( vl::Int i = 0; i < 3; ++i )
			for( vl::Int j = 0; j < 3; ++j )
				assert( fabs(rotMat[j][i]) > (1.0f - epsilon) 
					|| fabs(rotMat[i][j]) < epsilon );
	}

	// these don't matter
	btVector3 inertiaExtents(
		scale[0] + padding, 
		scale[0] + padding, 
		scale[0]*this->lengthRatio_ + padding);

	btVector3 localScale( scale[0] + padding, scale[0] + padding, 
		scale[0] * this->lengthRatio_ + padding );

	btCollisionShape* result = new btCylinderShapeZ( inertiaExtents );
	return result;
}
#endif

dMass CylinderObject::getOdeMassProps() const
{
	dMass result;
	vl::Vec3f scale = this->fullScale();
	dMassSetCylinder( &result, this->material()->density(), 3, scale[0], 2.0*lengthRatio_*scale[0] );
	return result;
}

RibObjectPtr CylinderObject::ribObject() const
{
	return ribObjectWithAttribute("ObjectInstance", 4);
}

PlaneObject::PlaneObject( const std::string& name, boost::shared_ptr<Material> material )
	: PhysicsObject( name, material )
{
	setStatic(true);
}

PlaneObject::~PlaneObject()
{
}

BoundingBox3f PlaneObject::bounds(const vl::Mat4f& transform) const
{
	BoundingBox3f result;
	vl::Vec3f scale;
	for( vl::Int i = 0; i < 3; ++i )
		scale[i] = vl::len( vl::Vec3f( transform[0][i], transform[1][i], transform[2][i] ) );
		
	result.expand( vl::xform( transform, vl::Vec3f(-scale[0], 0.0f, -scale[2]) ) );
	result.expand( vl::xform( transform, vl::Vec3f(-scale[0], 0.0f, scale[2]) ) );
	result.expand( vl::xform( transform, vl::Vec3f(scale[0], 0.0f, scale[2]) ) );
	result.expand( vl::xform( transform, vl::Vec3f(scale[0], 0.0f, -scale[2]) ) );
	return result;
}

PlaneObject* PlaneObject::cloneLocal() const
{
	return new PlaneObject( *this );
}

RibObjectPtr PlaneObject::ribObject() const
{
	return ribObjectWithAttribute("ObjectInstance", 2);
}

boost::shared_ptr<MechModelDeclaredObject> PlaneObject::mechModelObject() const
{
	return boost::shared_ptr<MechModelDeclaredObject>(
		new MechModelDeclaredObject( "plane", this->name() ) );
}

#ifdef USE_NOVODEX
std::vector< boost::shared_ptr<NxShapeDesc> > PlaneObject::getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const
{
	vl::Mat4f transform = this->transform();
	vl::Mat3f upperLeftBlock = toMat3f( transform );

	boost::shared_ptr<NxPlaneShapeDesc> shapeDesc( new NxPlaneShapeDesc );

	vl::Vec3f normal = trans( inv( upperLeftBlock ) ) * vl::Vec3f( 0.0, 1.0, 0.0 );
	normal = vl::norm( normal );
	shapeDesc->normal.set( normal.Ref() );

	vl::Vec3f position = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	float d = dot( normal, position );
	shapeDesc->d = d + padding + defaultSkinWidth();

	return std::vector< boost::shared_ptr<NxShapeDesc> > (1, shapeDesc);
}
#endif

#ifdef USE_BULLET
btCollisionShape* PlaneObject::bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const
{
	/*
	vl::Mat4f transform = this->transform();
	vl::Mat3f upperLeftBlock = toMat3f( transform );

	vl::Vec3f normal = trans( inv( upperLeftBlock ) ) * vl::Vec3f( 0.0, 1.0, 0.0 );
	normal = vl::norm( normal );

	vl::Vec3f position = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	float d = dot( normal, position );

	return new btStaticPlaneShape( toBtVector3(normal), d );
	*/
	return new btStaticPlaneShape( btVector3(0.0, 1.0, 0.0), padding );
}
#endif


ODEGeomResult PlaneObject::getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const
{
	assert( centerOfMass == vl::vl_0 );

	vl::Mat4f transform = this->transform();
	vl::Mat3f upperLeftBlock = toMat3f( transform );

	vl::Vec3f normal = trans( inv( upperLeftBlock ) ) * vl::Vec3f( 0.0, 1.0, 0.0 );
	normal = vl::norm( normal );

	vl::Vec3f position = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	float d = dot( normal, position ) + padding;

	return ODEGeomResult( dCreatePlane( space, 
		normal[0], normal[1], normal[2], 
		d ) );
}

#ifdef USE_NEWTON
NewtonCollision* PlaneObject::getNewtonCollision( 
	const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const
{
	return 0;
}
#endif

dMass PlaneObject::getOdeMassProps() const
{
	// should never reach here!!
	assert( false );

	dMass result;
	dMassSetSphere( &result, this->material()->density(), 100.0f );
	return result;
}


void PlaneObject::renderLocalFill(float alpha, const char* statePtr) const
{
#ifdef GUI
	double textureRepeat = 10.0;
	double extent = 500.0;

    const GLTexture& texture = checkerboardTextureMap();
    GLBindTextureHandler bindTexture(texture);

	GLActionHandler quads( GL_QUADS );
	glNormal3f( 0.0, textureRepeat, 0.0 );
	glTexCoord2f( 0.0, 0.0 );
	glVertex3f( -extent, 0.0, -extent );
	glTexCoord2f( 0.0, textureRepeat );
	glVertex3f( -extent, 0.0,  extent );
	glTexCoord2f( textureRepeat, textureRepeat );
	glVertex3f(  extent, 0.0,  extent );
	glTexCoord2f( textureRepeat, 0.0 );
	glVertex3f(  extent, 0.0, -extent );
#endif
}

void PlaneObject::renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const
{
#ifdef GUI
	double extent = 500.0;

	glDisable( GL_LIGHTING );
	glColor4fv( wireframeColor.Ref() );
	glLineWidth( 1.0 );

	unsigned int x = 0;
	unsigned int y = 2;

	for( int i = 0; i < 2; ++i )
	{
		vl::Vec3 dir_x = vl::vl_0;
		dir_x[x] = 1.0;

		vl::Vec3 dir_y = vl::vl_0;
		dir_y[y] = 1.0;

		{
			GLActionHandler lines( GL_LINES );

			int numLines = 10;

			for( int i = -numLines; i <= numLines; ++i )
			{
				vl::Vec3 begin = extent*static_cast<double>(i) / static_cast<double>(numLines) * dir_x - extent*dir_y;
				vl::Vec3 end   = extent*static_cast<double>(i) / static_cast<double>(numLines) * dir_x + extent*dir_y;

				glVertex3dv( begin.Ref() );
				glVertex3dv( end.Ref() );
			}
		}

		std::swap(x, y);
	}
#endif
}

bool PlaneObject::intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const
{
	vl::Vec3d N( 0.0, 1.0, 0.0);
	double Vd = dot( ray.getDirection(), N );

	const double epsilon = 1e-6;
	if( fabs(Vd) < epsilon )
		return false;			// parallel to plane

	if( Vd > 0.0 )
		return false;			// underside of plane

	double V0 = -dot( ray.getPosition(), N );
	t = V0/Vd;

	return (t >= 0.0);
}

bool PlaneObject::intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const
{
	return false;
}


class TriangleIntersectionHandler
{
public:
	TriangleIntersectionHandler( const vl::Vec3f* positions )
		: positions_(positions), t_( std::numeric_limits<double>::max() ) {}

	bool operator()( const Triangle& tri, const Ray3d& ray )
	{
		std::pair< bool, TriangleIntersection3d > result =
			intersectTriangle( ray, 
				toVec3d(positions_[tri[0]]), 
				toVec3d(positions_[tri[1]]), 
				toVec3d(positions_[tri[2]]) );
		if( result.first && result.second.t() >= 0.0 && result.second.t() < t_ )
		{
			t_ = result.second.t();
			intersectedTri_ = tri;
			baryCoords_ = result.second.barycentricCoordinates();
			return true;
		}
		else
			return false;
	}

	double t() const
	{
		return t_;
	}

	Triangle triangle() const
	{
		return intersectedTri_;
	}

	vl::Vec3d barycentricCoordinates() const
	{
		return baryCoords_;
	}

private:
	const vl::Vec3f* positions_;
	vl::Vec3d baryCoords_;
	double t_;
	Triangle intersectedTri_;
};


MeshObject::MeshObject( const std::string& name,
	TriangleMeshPtr mesh,
	boost::shared_ptr<Material> material )
	:	PhysicsObject( name, material )
{
	setMesh( mesh );
}

MeshObject::MeshObject( const std::string& name,
	boost::shared_ptr<Material> material )
	:	PhysicsObject( name, material )
{
}

void MeshObject::setConvexHulls(const ConvexHullList& hulls)
{
	this->convexHulls_ = hulls;
}

ConvexHullList MeshObject::getConvexHulls() const
{
	return this->convexHulls_;
}

ConvexHullPtr MeshObject::convexHullForPoints( const std::deque<unsigned int>& vertices )
{
	std::vector<vl::Vec3d> verts;
	verts.reserve( vertices.size() );

	const std::vector<vl::Vec3f>& positions = this->mesh_->positions();
	for( std::deque<unsigned int>::const_iterator vIter = vertices.begin();
		vIter != vertices.end(); ++vIter )
	{
		verts.push_back( toVec3d(positions[*vIter]) );
	}

	return this->convexHullForPoints( verts );
}

ConvexHullPtr MeshObject::convexHullForPoints( const std::vector<vl::Vec3d>& verts )
{
	ConvexDecomposition::HullDesc desc;

	desc.mVcount = verts.size();
	desc.mVertices = reinterpret_cast<const double*>( &verts[0] );
	desc.mVertexStride = sizeof( vl::Vec3d );
	desc.mSkinWidth = 0.0;
	desc.mFlags = ConvexDecomposition::QF_TRIANGLES;

	ConvexDecomposition::HullResult hullRes;
	ConvexDecomposition::HullLibrary lib;
	lib.CreateConvexHull( desc, hullRes );

	boost::shared_ptr<ConvexHull> hullPtr( new ConvexHull );
	ConvexHull& hull = *hullPtr;

	hull.vertices.resize( hullRes.mNumOutputVertices );
	std::copy( hullRes.mOutputVertices, hullRes.mOutputVertices + 3*hullRes.mNumOutputVertices,
		reinterpret_cast<float*>(&hull.vertices[0]) );

	assert( !hullRes.mPolygons );
	hull.triangles.reserve( hullRes.mNumFaces );
	for( size_t iTriangle = 0; iTriangle < hullRes.mNumFaces; ++iTriangle )
	{
		const unsigned int* startVert = hullRes.mIndices + 3*iTriangle;
		hull.triangles.push_back( ConvexHull::Triangle(startVert[0], startVert[1], startVert[2]) );
	}

	// need to compute face normals
	hull.normals.resize( hull.triangles.size() );
	for( size_t iTriangle = 0; iTriangle < hull.triangles.size(); ++iTriangle )
	{
		const TinyVec<int, 3>& tri = hull.triangles[iTriangle];
		hull.normals[ iTriangle ] = vl::norm( vl::cross( 
			hull.vertices[ tri[1] ] - hull.vertices[ tri[0] ],
			hull.vertices[ tri[2] ] - hull.vertices[ tri[0] ] ) );
	}

	return hullPtr;
}


vl::Vec3f MeshObject::furthestPoint( const vl::Vec3f& point ) const
{
	vl::Vec3f furthest = point;

	vl::Vec3f scale = this->fullScale();
	const std::vector<vl::Vec3f>& positions = mesh_->positions();
	for( std::vector<vl::Vec3f>::const_iterator posItr = positions.begin();
		posItr != positions.end(); ++posItr )
	{
		vl::Vec3f newPos = scale * (*posItr);

		if( vl::len(point - newPos) > vl::len(point - furthest) )
			furthest = newPos;
	}

	return furthest;
}

bool MeshObject::equalsLocal( const SceneGraphElement& other ) const
{
	// unimplemented
	return PhysicsObject::equalsLocal( other );
}


TriangleMeshPtr MeshObject::getMesh() const
{
	return this->mesh_;
}

void MeshObject::setMesh( TriangleMeshPtr mesh )
{
	this->mesh_ = mesh;
}

MeshObject* MeshObject::cloneLocal() const
{
	return new MeshObject( *this );
}

BoundingBox3f MeshObject::bounds(const vl::Mat4f& transform) const
{
	return transformBounds( this->mesh_->bounds(), transform );
}

void MeshObject::updateScaleLocal()
{
#if dTRIMESH_ENABLED
	this->triMeshData_.reset();
#endif
}

#ifdef USE_BULLET
std::vector<btVector3> expandHull( const ConvexHull& hull, const vl::Vec3f& scale, float amount )
{
	const std::vector<vl::Vec3f>& hullVertices = hull.vertices;

	btAlignedObjectArray<btVector3> vertices;
	vertices.reserve( hullVertices.size() );
	for( size_t i = 0; i < hullVertices.size(); ++i )
		vertices.push_back( toBtVector3( scale * hullVertices[i] ) );

	btAlignedObjectArray<btVector3> planeEquations;
	planeEquations.reserve( 4*vertices.size() );
	btGeometryUtil::getPlaneEquationsFromVertices(vertices,planeEquations);

	// Bullet gets unhappy if we don't have a margin, so we'll pull it out of the convex hull.
	btAlignedObjectArray<btVector3> shiftedPlaneEquations;
	shiftedPlaneEquations.reserve( planeEquations.size() );
	for (int p = 0; p < planeEquations.size(); p++)
	{
		btVector3 plane = planeEquations[p];
		plane[3] -= amount;
		shiftedPlaneEquations.push_back(plane);
	}
	btAlignedObjectArray<btVector3> shiftedVertices;
	btGeometryUtil::getVerticesFromPlaneEquations(shiftedPlaneEquations, shiftedVertices);

	// It's not safe to copy a btAlignedObjectArray:
	std::vector<btVector3> result( &shiftedVertices[0], &shiftedVertices[0] + shiftedVertices.size() );
	return result;
}

btCollisionShape* MeshObject::bulletCollisionShape( const vl::Mat4f& offsetMatrix, float padding ) const
{
	vl::Vec3f scale = this->fullScale();
	std::deque<btCollisionShape*> meshShapes;
	btScalar collisionMargin = 0.04;
	for( ConvexHullList::const_iterator itr = this->convexHulls_.begin();
		itr != this->convexHulls_.end(); ++itr )
	{
		std::vector<btVector3> shiftedVertices = expandHull( *(*itr), scale, padding - collisionMargin );
        /*
	    const std::vector<vl::Vec3f>& hullVertices = (*itr)->vertices;

	    std::vector<btVector3> shiftedVertices;
	    shiftedVertices.reserve( hullVertices.size() );
	    for( size_t i = 0; i < hullVertices.size(); ++i )
		    shiftedVertices.push_back( toBtVector3( scale * hullVertices[i] ) );
            */

		if( shiftedVertices.empty() )
			continue;

		btCollisionShape* convexShape = new btConvexHullShape(&(shiftedVertices[0].getX()),shiftedVertices.size());
		assert( convexShape->getMargin() == collisionMargin );
		convexShape->setMargin( collisionMargin );
		meshShapes.push_back( convexShape );
	}

	if( meshShapes.empty() )
		return 0;

	vl::Mat4f invOffsetMatrix = rigidInverse( offsetMatrix );
	btCompoundShape* result = new btCompoundShape();
	for( std::deque<btCollisionShape*>::const_iterator itr = meshShapes.begin();
		itr != meshShapes.end(); ++itr )
	{
		result->addChildShape( toBtTransform( invOffsetMatrix ), *itr );
	}

	return result;
}
#endif

ODEGeomResult MeshObject::getOdeGeom(dSpaceID space, const vl::Vec3f& centerOfMass, float padding) const
{
	/*
	BoundingBox3f bounds = this->mesh_->bounds();
	assert( !bounds.empty() );
	vl::Vec3f mid = 0.5*(bounds.minimum() + bounds.maximum());

	std::vector<dGeomID> result;
	std::vector<SphereWrapperTree::Sphere> spheres = 
		this->mesh_->sphereWrapperTree().level(3);
	result.reserve( spheres.size() );

	vl::Vec3f s = this->scale();
	assert( s[0] == s[1] && s[1] == s[2] );
	float scale = s[0];

	for( std::vector<SphereWrapperTree::Sphere>::const_iterator iter = spheres.begin();
		iter != spheres.end(); ++iter )
	{
		assert( iter->radius > 0.0f );

		dGeomID sphereId = dCreateSphere(0, scale*iter->radius);
		dGeomID transformId = dCreateGeomTransform(space);
		dGeomTransformSetGeom( transformId, sphereId );
		dGeomTransformSetCleanup( transformId, 1 );

		vl::Vec3f offset = scale*(mid + iter->center);

		dGeomSetPosition( sphereId, offset[0], offset[1], offset[2] );
		result.push_back( transformId );
	}
	*/

	ODEGeomResult result;
#if dTRIMESH_ENABLED
	if( !triMeshData_ )
	{
		std::vector<vl::Vec3f> positions; positions.reserve( this->mesh_->positions().size() );
		const vl::Vec3f scale = this->fullScale();
		for( std::vector<vl::Vec3f>::const_iterator itr = this->mesh_->positions().begin();
			itr != this->mesh_->positions().end(); ++itr )
		{
			positions.push_back( (*itr) * scale - centerOfMass );
		}

		triMeshData_.reset( 
			new ODETriMeshDataWrapper( positions, this->mesh_->triangles() ) );
	}

	boost::shared_ptr<ODETriMeshWrapper> triMeshWrapper(
		new ODETriMeshWrapper(triMeshData_) );
	result.triMeshes.push_back( triMeshWrapper );
	// the 0s are for callback functions:
	result.geoms.push_back( dCreateTriMesh(space, triMeshWrapper->get(), 0, 0, 0) );
#endif

	return result;
}

#ifdef USE_NEWTON
NewtonCollision* MeshObject::getNewtonCollision( 
	const NewtonWorld* newtonWorld, const dFloatMat4& offsetMatrix, float padding ) const
{
	return 0;
}
#endif

dMass MeshObject::getOdeMassProps() const
{
	MassProperties massProps = computeMassProperties( 
		*this->mesh_->objFile(), 
		this->material()->density(),
		toVec3d(this->fullScale()) );

	dMass result;
	dMassSetParameters(
		&result,
		massProps.mass,
		0.0f, 0.0f, 0.0f,
		massProps.inertiaTensor[0][0], massProps.inertiaTensor[1][1], massProps.inertiaTensor[2][2], 
		massProps.inertiaTensor[0][1], massProps.inertiaTensor[1][2], massProps.inertiaTensor[2][0] );

	dMassTranslate( &result, 
		massProps.centerOfMass[0], massProps.centerOfMass[1], massProps.centerOfMass[2] );

	return result;	
}

#ifdef USE_NOVODEX
class InitCooking
{
public:
	InitCooking()
	{
		// Recursive lock should be ok
		boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
		gCooking = NxGetCookingLib(NX_PHYSICS_SDK_VERSION);
		gCooking->NxInitCooking();
	}

	~InitCooking()
	{
		boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
		gCooking->NxCloseCooking();
	}

	NxCookingInterface* get()
	{
		return this->gCooking;
	}

private:
	NxCookingInterface *gCooking;
};


std::vector< boost::shared_ptr<NxShapeDesc> > MeshObject::getNxShapeDesc(NxPhysicsSDKPtr sdk, float padding) const
{
	vl::Vec3f scale = this->fullScale();
	for( vl::Int i = 0; i < 3; ++i )
		scale[i] += padding;

	std::vector< boost::shared_ptr<NxShapeDesc> > result;
	{
		boost::recursive_mutex::scoped_lock lock( novodexCreateDeleteMutex );
		static InitCooking cooking;

		for( ConvexHullList::const_iterator itr = this->convexHulls_.begin();
			itr != this->convexHulls_.end(); ++itr )
		{
			const ConvexHull& hull = *(*itr);

			std::vector<btVector3> expandedVertices = 
				expandHull( *(*itr), scale, padding /* + defaultSkinWidth() */ );
			if( expandedVertices.empty() )
				continue;

			std::vector<NxVec3> vertices( expandedVertices.size() );
			for( size_t i = 0; i < vertices.size(); ++i )
			{
				for( size_t j = 0; j < 3; ++j )
					vertices[i][j] = expandedVertices[i][j];
			}

			NxConvexMeshDesc convexDesc;
			convexDesc.numVertices = vertices.size();
			convexDesc.pointStrideBytes = sizeof( NxVec3 );
			convexDesc.points = &vertices[0];
			convexDesc.flags = NX_CF_COMPUTE_CONVEX;
			assert( convexDesc.isValid() );

			MemoryWriteBuffer buf;
			bool status = cooking.get()->NxCookConvexMesh(convexDesc, buf);
			if( !status )
				continue;

			boost::shared_ptr<NxConvexShapeDesc> shapeDesc( new NxConvexShapeDesc );
			shapeDesc->meshData = sdk->createConvexMesh(MemoryReadBuffer(buf.data));
			result.push_back( shapeDesc );
		}
	}

	return result;
}
#endif

MeshObject::~MeshObject()
{
}

void MeshObject::toMayaAsciiLocal( std::ofstream& ofs, const std::string& materialName, const std::string& prefix ) const
{
	std::string name = prefix + this->name();
	std::string shapeName = name + "Shape";
	ofs << "createNode mesh -n \"" << shapeName << "\" -p \"" << name << "\";\n";
//	ofs << "\tsetAttr -k off \".v\";\n";
//	ofs << "\tsetAttr \".io\" yes;\n";

	const std::vector<vl::Vec3f>& positions = this->mesh_->positions();
	const std::vector<Triangle>& triangles = this->mesh_->triangles();
	setAttr( ofs, ".vt", positions );

	// gather up all edges
	typedef std::pair<size_t, size_t> Edge;
	std::deque< Edge > edges;
	for( std::vector<Triangle>::const_iterator itr = triangles.begin();
		itr != triangles.end(); ++itr )
	{
		for( size_t j = 0; j < 3; ++j )
		{
			Edge edge( (*itr)[j], (*itr)[(j+1)%3] );
			if( edge.second < edge.first )
				std::swap( edge.first, edge.second );
			edges.push_back( edge );
		}
	}

	// dump edges
	std::sort( edges.begin(), edges.end() );
	edges.erase( std::unique( edges.begin(), edges.end() ), edges.end() );
	{
		std::vector< TinyVec<size_t, 3> > edgeTemp;
		for( size_t iEdge = 0; iEdge < edges.size(); ++iEdge )
		{
			edgeTemp.push_back( 
				TinyVec<size_t, 3>( 
					edges[iEdge].first, edges[iEdge].second, 0 ) );
		}
		setAttr( ofs, ".ed", edgeTemp );
	}

	// dump normals
	std::vector<vl::Vec3f> normals( 2*edges.size(), 
		vl::Vec3f( 1e20, 1e20, 1e20 ) );
	setAttr( ofs, ".n", normals );

	ofs << "\tsetAttr -s " << triangles.size() << " \".fc\";\n";
	for( size_t start = 0; start < triangles.size(); start += 500 )
	{
		const size_t end = std::min( start+500, triangles.size() );
		ofs << "\tsetAttr \".fc[" << start << ":" << (end-1) << "]\" -type \"polyFaces\"";
		// now dump faces
		for( size_t iTriangle = start; iTriangle < end; ++iTriangle )
		{
			ofs << "\n\t\tf 3 ";

			Triangle tri = triangles.at( iTriangle );
			for( size_t j = 0; j < 3; ++j )
			{
				Edge edge( tri[j], tri[(j+1)%3] );
				size_t offset = 0;
				if( edge.second < edge.first )
				{
					std::swap( edge.first, edge.second );
					ofs << "-";
					offset = 1;
				}

				typedef std::deque<Edge>::iterator Iter;
				typedef std::pair<Iter, Iter> IterPr;
				IterPr pr = std::equal_range( edges.begin(), edges.end(), edge );
				assert( pr.second != pr.first );

				size_t iEdge = std::distance( edges.begin(), pr.first ) + offset;
				ofs << iEdge << " ";
			}
		}
		ofs << ";\n";
	}

}

RibObjectPtr MeshObject::ribObject() const
{
	boost::shared_ptr<RibNamedBlock> result( new RibNamedBlock("Transform") );
	boost::shared_ptr<RibObjectWithAttributes> translateObject(
		new RibObjectWithAttributes( "Translate" ) );

	boost::shared_ptr<RibMetaComment> objComment(
		new RibMetaComment( "METARIB-OBJ " + this->mesh_->objFile()->filename() ) );
	result->addChild( objComment );
	return result;
}

boost::shared_ptr<MechModelDeclaredObject> MeshObject::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "mesh", this->name() ) );

	result->add( mechModelValuePair("objFile", this->mesh_->objFile()->filename()) );
	for( ConvexHullList::const_iterator hullItr = this->convexHulls_.begin();
		hullItr != this->convexHulls_.end(); ++hullItr )
	{
		const std::vector<vl::Vec3f>& verts = (*hullItr)->vertices;
		result->add( mechModelValuePair("hull", verts) );
	}

	return result;
}

void MeshObject::renderLocalFill(float alpha, const char* statePtr) const
{
#ifdef GUI
	mesh_->displayList()->render();
#endif
}

void MeshObject::renderPoints( const std::deque<unsigned int>& selected, 
	const vl::Mat4f& transform, 
	const char* state ) const
{
	if( !this->visible() )
		return;

	glDisable( GL_LIGHTING );
	glPointSize( 5.0 );
	std::vector<vl::Vec3f> colors( this->mesh_->positions().size(), 
		vl::Vec3f(253.0/255.0, 255.0/255.0, 58.0/255.0) );
	for( std::deque<unsigned int>::const_iterator selectedItr = selected.begin();
		selectedItr != selected.end(); ++selectedItr )
	{
		colors[ *selectedItr ] = vl::Vec3f(253.0/255.0, 58.0/255.0, 58.0/255.0);
	}

	vl::Mat4f transpose = trans(transform);
	GLMatrixStackHandler handler;
	glMultMatrixf( transpose.Ref() );

	GLEnableClientStateHandler vertexArray(GL_VERTEX_ARRAY);
	GLEnableClientStateHandler colorArray(GL_COLOR_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &mesh_->positions()[0]);
	glColorPointer(3, GL_FLOAT, 0, &colors[0]);
	glDrawArrays( GL_POINTS, 0, mesh_->positions().size() );
}

std::pair<bool, unsigned int> MeshObject::selectPoint( const vl::Vec2d& clicked, 
	const ScreenSpaceConverter& converter, 
	const vl::Mat4f& transform,
	const char* state ) const
{
	if( !visible() )
		return std::make_pair( false, 0 );

	double closest = boost::numeric::bounds<double>::highest();
	unsigned int closestPoint;

	const std::vector<vl::Vec3f> positions = this->mesh_->positions();
	for( size_t iVertex = 0; iVertex < positions.size(); ++iVertex )
	{
		vl::Vec2d screenPos = vl::strip( converter.toScreenSpace( toVec3d(vl::xform( transform, positions[iVertex] )) ) );
		if( vl::len( screenPos - clicked ) < closest )
		{
			closest = vl::len( screenPos - clicked );
			closestPoint = iVertex;
		}
	}

	if( closest < 10.0 )
		return std::make_pair(true, closestPoint);
	else
		return std::make_pair(false, 0);
}

void MeshObject::renderLocalHulls(
	const float alpha, 
	const char* statePtr) const
{
#ifdef GUI
	glShadeModel( GL_FLAT );
	glPolygonMode( GL_FRONT, GL_FILL );

	for( size_t iHull = 0; iHull < this->convexHulls_.size(); ++iHull )
	{
		const ConvexHull& hull = *convexHulls_[iHull];
		const std::vector<Color>& colors = mapColors();
		Color color = colors.at( iHull % colors.size() );
		vl::Vec4f c(vl::vl_1);
		for( vl::Int i = 0; i < 3; ++i )
			c[i] = static_cast<float>(color[i]) / 255.0f;
		glMaterialfv( GL_FRONT, GL_DIFFUSE, c.Ref() );

		GLActionHandler handler( GL_TRIANGLES );
		for( unsigned int iTriangle = 0; iTriangle < hull.triangles.size(); ++iTriangle )
		{
			assert( !isBad( hull.normals[iTriangle] ) );
			glNormal3fv( hull.normals[iTriangle].Ref() );

			const ConvexHull::Triangle& tri = hull.triangles[iTriangle];
			for( unsigned int i = 0; i < 3; ++i )
				glVertex3fv( hull.vertices[ tri[i] ].Ref() );
		}
	}
	glShadeModel( GL_SMOOTH );
#endif
}

void MeshObject::renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const
{
#ifdef GUI
	glDisable( GL_LIGHTING );
	glPolygonMode( GL_FRONT, GL_LINE );
	glColor4fv( wireframeColor.Ref() );
	glLineWidth( 1.0 );

	GLEnableClientStateHandler vertexArray(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &mesh_->positions()[0]);
	/*
	will reenable eventually, once I figure out what is causing the weird bug in the triangle stripping

	for( std::deque<PrimGroup>::const_iterator stripItr = mesh_->strips.begin();
		stripItr != mesh_->strips.end(); ++stripItr )
	{
		glDrawElements(stripItr->glPrimitiveType(), stripItr->indices.size(), GL_UNSIGNED_SHORT, 
			&(stripItr->indices)[0]);
	}
	*/

	glDrawElements( GL_TRIANGLES, 3*mesh_->triangles().size(), GL_UNSIGNED_INT, &(mesh_->triangles())[0] );
#endif
}

bool MeshObject::intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const
{
	TriangleIntersectionHandler handler( &(this->mesh_->positions())[0] );
	if( mesh_->triangleTree().intersectRay( ray, handler ) )
	{
		t = handler.t();
		return true;
	}

	return false;
}

class BoxTriangleIntersectionHandler
{
public:
	BoxTriangleIntersectionHandler( const BoundingBox2d& box, const ScreenSpaceConverter& converter, 
		const vl::Vec3f* positions )
		: box_(box), converter_(converter), positions_(positions)
	{
	}

	bool operator()( const Triangle& tri )
	{
		for( unsigned int i = 0; i < 3; ++i )
		{
			vl::Vec3f pos = positions_[ tri[i] ];
			if( contains( box_, strip(converter_.toScreenSpace( toVec3d(pos) ) ) ) )
				return true;
		}

		return false;
	}
	
private:
	const BoundingBox2d box_;
	const ScreenSpaceConverter& converter_;
	const vl::Vec3f* positions_;
};

bool MeshObject::intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const
{
	BoxTriangleIntersectionHandler handler( box, converter, &(this->mesh_->positions())[0] );
	return mesh_->triangleTree().intersect( box, converter, handler );
}

// defaults
Material::Material( const std::string& name )
	:	Named(name), 
		density_(200.0),
		dynamicFriction_(0.05),
		staticFriction_(0.1), 
		restitution_(0.2),
		color_(0.8, 0.8, 0.8),
		youngsModulus_(1.9E11), // steel
		poissonRatio_(0.3),
		raleighAlpha_(0.001),
		raleighBeta_(0.001),
		coulombFriction_(true)
{
}

Material::Material( const std::string& name, float density, 
				   float dynamicFriction, float staticFriction, float restitution )
	:	Named(name), 
		density_(density), 
		dynamicFriction_(dynamicFriction), 
		staticFriction_(staticFriction), 
		restitution_(restitution),
		youngsModulus_(1.9E11), // steel
		poissonRatio_(0.3),
		color_( 0.8, 0.8, 0.8 ),
		raleighAlpha_(0.001),
		raleighBeta_(0.001),
		coulombFriction_(true)
{
}

bool Material::equals( const Material& other ) const
{
	return withinEpsilon( this->density_, other.density_ )
		&& withinEpsilon( this->dynamicFriction_, other.dynamicFriction_ )
		&& withinEpsilon( this->staticFriction_, other.staticFriction_ )
		&& withinEpsilon( this->restitution_, other.restitution_ )
		&& withinEpsilon( this->youngsModulus_, other.youngsModulus_ )
		&& withinEpsilon( this->poissonRatio_, other.poissonRatio_ )
		&& withinEpsilon( this->restitution_, other.restitution_ )
		&& withinEpsilon( this->raleighAlpha_, other.raleighAlpha_ )
		&& withinEpsilon( this->raleighBeta_, other.raleighBeta_ )
		&& this->coulombFriction_ == other.coulombFriction_;
}

Material::~Material()
{
}

MechModelObjectPtr Material::toMechModel(MechModelType type) const
{
	boost::shared_ptr<MechModelDeclaredObject> object(
		new MechModelDeclaredObject( "material", this->name() ) );

	object->add( mechModelValuePair("density", density_) );

	if( type != MECH_MODEL_AUDIO && type != MECH_MODEL_AUDIO_ROOT )
	{
		object->add( mechModelValuePair("dynamicFriction", dynamicFriction_) );
		object->add( mechModelValuePair("staticFriction", staticFriction_) );
		object->add( mechModelValuePair("restitution", restitution_) );
		object->add( mechModelValuePair("color", color_) );
		if( !coulombFriction_ )
			object->add( mechModelValuePair("coulomb", coulombFriction_) );

		object->add( mechModelValuePair("raleighAlpha", raleighAlpha_) );
		object->add( mechModelValuePair("raleighBeta", raleighBeta_) );
	}

	// don't need damping parameters when computing modes
	object->add( mechModelValuePair("youngsModulus", youngsModulus_) );
	object->add( mechModelValuePair("poissonRatio", poissonRatio_) );

	return object;
}

float Material::youngsModulus() const
{
	return this->youngsModulus_;
}

float Material::poissonRatio() const
{
	return this->poissonRatio_;
}

float Material::raleighAlpha() const
{
	return this->raleighAlpha_;
}

float Material::raleighBeta() const
{
	return this->raleighBeta_;
}

float Material::density() const
{
	return density_;
}

float Material::dynamicFriction() const
{
	return dynamicFriction_;
}

float Material::staticFriction() const
{
	return staticFriction_;
}

float Material::restitution() const
{
	return restitution_;
}

vl::Vec3f Material::color() const
{
	return color_;
}

bool Material::coulombFriction() const
{
	return this->coulombFriction_;
}

void Material::setDensity( float density )
{
	density_ = density;
}

void Material::setDynamicFriction( float friction )
{
	dynamicFriction_ = friction;
}

void Material::setStaticFriction( float friction )
{
	staticFriction_ = friction;
}

void Material::setRestitution( float restitution )
{
	restitution_ = restitution;
}

void Material::setCoulombFriction( bool value )
{
	this->coulombFriction_ = value;
}

void Material::setColor( const vl::Vec3f& color )
{
	color_ = color;
}

void Material::setYoungsModulus(float value)
{
	this->youngsModulus_ = value;
}

void Material::setPoissonRatio(float value)
{
	this->poissonRatio_ = value;
}

void Material::setRaleighAlpha(float value)
{
	this->raleighAlpha_ = value;
}

void Material::setRaleighBeta(float value)
{
	this->raleighBeta_ = value;
}


class BoundsForTriangle
{
public:
	BoundsForTriangle( const vl::Vec3f* positions )
		: positions_(positions) {}

	BoundingBox3d operator()( const Triangle& tri ) const
	{
		BoundingBox3d result;
		for( unsigned int iVertex = 0; iVertex < 3; ++iVertex )
			result.expand( toVec3d(positions_[ tri[iVertex] ]) );
		return result;
	}

private:
	const vl::Vec3f* positions_;
};

TriangleMesh::TriangleMesh( ObjFilePtr objFile )
{
	init( objFile );
}

void TriangleMesh::init( ObjFilePtr objFile )
{
	objFile_ = objFile;

	positions_.reserve( objFile->vertexCount() );
	for( unsigned int iVertex = 0; iVertex < objFile->vertexCount(); ++iVertex )
	{
		vl::Vec3f pos = toVec3f(objFile->vertexPosition(iVertex + 1));
		positions_.push_back( pos );
		bounds_.expand( pos );
	}

	normals_.resize( objFile->vertexCount(), vl::vl_0 );

	triangles_.reserve( objFile->faceCount() );
	for( unsigned int iGroup = 0; iGroup < objFile->numGroups(); ++iGroup )
	{
		ObjFile::Group group = objFile->group( iGroup );
		for( unsigned int iFace = 0; iFace < group.faceCount(); ++iFace )
		{
			ObjFile::Face<unsigned int> face = group.face( iFace );

			vl::Vec3d faceNormal = cross(
				objFile->vertexPosition( face.vertex(1).position() ) - 
					objFile->vertexPosition( face.vertex(0).position() ),
				objFile->vertexPosition( face.vertex(2).position() ) - 
					objFile->vertexPosition( face.vertex(0).position() ) );
			faceNormal = vl::norm( faceNormal );

			// triangulate it:
			for( size_t i = 1; (i+1) < face.vertexCount(); ++i )
			{
				triangles_.push_back( 
					Triangle( 
						face.vertex(0).position() - 1,
						face.vertex(i).position() - 1,
						face.vertex(i+1).position() - 1 ) );
			}

			for( unsigned int iVertex = 0; iVertex < face.vertexCount(); ++iVertex )
			{
				vl::Vec3d normal = faceNormal;
				if( face.vertex(iVertex).hasNormal() )
					normal = objFile->normal( face.vertex(iVertex).normal() );
				normals_[ face.vertex(iVertex).position() - 1 ] += toVec3f(normal);
			}
		}
	}

	triangleTree_.reset( 
		new AABBTree<Triangle>( triangles_, BoundsForTriangle(&positions_[0]) ) );

	for( size_t iNormal = 0; iNormal < normals_.size(); ++iNormal )
		normals_[iNormal] = vl::norm( normals_[iNormal] );


	/*
	just linking against the triangle stripping library 
	causes a really bizarre memory error, so we'll leave it out.

	{
		using namespace triangle_stripper;
		std::vector<unsigned int> triangleIndices;
		triangleIndices.reserve( 3*mesh->triangles.size() );
		for( std::vector<Triangle>::const_iterator triangleItr = mesh->triangles.begin();
			triangleItr != mesh->triangles.end(); ++triangleItr )
		{
			for( unsigned int i = 0; i < 3; ++i )
				triangleIndices.push_back( (*triangleItr)[i] );
		}

		tri_stripper triStripper(triangleIndices);
		triStripper.SetMinStripSize(2);
		triStripper.SetCacheSize(16);
		triStripper.SetBackwardSearch(false);

		primitive_vector PrimitivesVector;
		triStripper.Strip(&PrimitivesVector);

		for (size_t i = 0; i < PrimitivesVector.size(); ++i)
		{
			PrimGroup primGroup;
			if (PrimitivesVector[i].Type == TRIANGLE_STRIP)
				primGroup.type = PrimGroup::STRIP;
			else
				primGroup.type = PrimGroup::LIST;

			primGroup.indices.reserve( PrimitivesVector[i].Indices.size() );
			std::copy( PrimitivesVector[i].Indices.begin(), 
				PrimitivesVector[i].Indices.end(), 
				std::back_inserter(primGroup.indices) );

			mesh->strips.push_back( primGroup );
		}
	}
	*/
}

ConvexHullList TriangleMesh::computeConvexHulls( int depth, int maxVertsPerHull ) const
{
	ConvexHullList result;

	// perform convex decomposition
	{
		ConvexDecomposition::DecompDesc desc;

		desc.mDepth = depth;
		desc.mMaxVertices = maxVertsPerHull;

		desc.mVcount = this->positions_.size();

		const float* vertPositions = reinterpret_cast<const float*>(&this->positions_[0]);
		std::vector<double> vertices( vertPositions, vertPositions + 3*desc.mVcount );
		desc.mVertices = &vertices[0];

		desc.mTcount = this->triangles_.size();
		const unsigned int* trianglePtr = 
			reinterpret_cast<const unsigned int*>(&this->triangles_[0]);
		std::vector<unsigned int> triangles( trianglePtr, trianglePtr + 3*triangles_.size() );
		desc.mIndices = &triangles[0];

		class Callback
			: public ConvexDecomposition::ConvexDecompInterface
		{
		public:
			Callback( ConvexHullList& hulls )
				: hulls_(hulls) {}

			void ConvexDecompResult(ConvexDecomposition::ConvexResult &result)
			{
				boost::shared_ptr<ConvexHull> hullPtr( new ConvexHull );
				ConvexHull& hull = *hullPtr;

				hull.vertices.resize( result.mHullVcount );
				std::copy( result.mHullVertices, result.mHullVertices + 3*result.mHullVcount,
					reinterpret_cast<float*>(&hull.vertices[0]) );

				hull.triangles.reserve( result.mHullTcount );
				// winding is backwards
				for( size_t iTriangle = 0; iTriangle < result.mHullTcount; ++iTriangle )
				{
					const unsigned int* startVert = result.mHullIndices + 3*iTriangle;
					hull.triangles.push_back( ConvexHull::Triangle(startVert[0], startVert[1], startVert[2]) );
				}

				// need to compute face normals
				hull.normals.resize( hull.triangles.size() );
				for( size_t iTriangle = 0; iTriangle < hull.triangles.size(); ++iTriangle )
				{
					const TinyVec<int, 3>& tri = hull.triangles[iTriangle];
					hull.normals[ iTriangle ] = vl::norm( vl::cross( 
						hull.vertices[ tri[1] ] - hull.vertices[ tri[0] ],
						hull.vertices[ tri[2] ] - hull.vertices[ tri[0] ] ) );
				}

				hulls_.push_back( hullPtr );
			}

		private:
			ConvexHullList& hulls_;
		};

		Callback callback( result );
		desc.mCallback = &callback;

		performConvexDecomposition( desc );
	}

	return result;
}

BoundingBox3f TriangleMesh::bounds() const
{
	return this->bounds_;
}

#ifdef GUI
boost::shared_ptr<GLDisplayList> TriangleMesh::displayList() const
{
	if( !displayList_ )
	{
		displayList_.reset( new GLDisplayList() );
		GLDisplayList::CreateList createList( *displayList_ );

		for( size_t iGroup = 0; iGroup < objFile_->numGroups(); ++iGroup )
			objFile_->render( objFile_->group(iGroup).name() );
	}

	return displayList_;
}
#endif

const std::vector<vl::Vec3f>& TriangleMesh::positions() const
{
	return positions_;
}

const std::vector<vl::Vec3f>& TriangleMesh::normals() const
{
	return normals_;
}

const std::vector<Triangle>& TriangleMesh::triangles() const
{
	return triangles_;
}

const AABBTree<Triangle>& TriangleMesh::triangleTree() const
{
	return *triangleTree_;
}

ObjFilePtr TriangleMesh::objFile() const
{
	return this->objFile_;
}

StaticTriangleMesh::StaticTriangleMesh( 
	const float* positions, 
	size_t numVertices,
	const float* normals, 
	const unsigned int* triangles,
	size_t numTriangles,
	const unsigned int* boundarySegments,
	size_t numSegments )
	:	positions_( reinterpret_cast<const vl::Vec3f*>( positions ) ),
		normals_( reinterpret_cast<const vl::Vec3f*>( normals ) ),
		triangles_( reinterpret_cast<const Triangle*>( triangles ) ),
		segments_( boundarySegments ),
		numVertices_( numVertices ),
		numTriangles_( numTriangles ),
		numSegments_( numSegments )
{
	triangleTree_.reset( 
		new AABBTree<Triangle>( 
			std::vector<Triangle>(this->triangles_, this->triangles_+numTriangles_), 
			BoundsForTriangle(positions_) ) );

	for( size_t i = 0; i < this->numVertices_; ++i )
		this->bounds_.expand( this->positions_[i] );
}

const AABBTree<Triangle>& StaticTriangleMesh::triangleTree() const
{
	return *this->triangleTree_;
}

BoundingBox3f StaticTriangleMesh::bounds() const
{
	return this->bounds_;
}

#ifdef GUI
boost::shared_ptr<GLDisplayList> StaticTriangleMesh::displayList() const
{
	if( !displayList_ )
	{
		displayList_.reset( new GLDisplayList() );
		GLDisplayList::CreateList createList( *displayList_ );

		GLEnableClientStateHandler vertexArray(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, this->positions_);

		GLEnableClientStateHandler normalArray(GL_NORMAL_ARRAY);
		glNormalPointer(GL_FLOAT, 0, this->normals_);

		glDrawElements( GL_TRIANGLES, 3*this->numTriangles_, GL_UNSIGNED_INT, 
			this->triangles_ );
	}

	return displayList_;
}
#endif

const unsigned int* StaticTriangleMesh::segments() const
{
	return this->segments_;
}

size_t StaticTriangleMesh::numSegments() const
{
	return this->numSegments_;
}

const vl::Vec3f* StaticTriangleMesh::positions() const
{
	return this->positions_;
}

size_t StaticTriangleMesh::numVertices() const
{
	return this->numVertices_;
}


int dMyTriMeshClass;

/*
int collideMeshBox(dGeomID o1, dGeomID o2, int flags,
	dContactGeom *contact, int skip)
{
}

int collideMeshPlane(dGeomID o1, dGeomID o2, int flags,
	dContactGeom *contact, int skip)
{
}

int collideMeshMesh(dGeomID o1, dGeomID o2, int flags,
	dContactGeom *contact, int skip)
{
}
*/

/*
void dGetAABBFn(dGeomID g, dReal aabb[6])
{
}
*/

struct ODETriMeshStructure
{
	ODETriMeshStructure( boost::shared_ptr<TriangleMesh> m, const vl::Vec3f& s )
		: mesh(m), scale(s) {}

	boost::shared_ptr<TriangleMesh> mesh;
	vl::Vec3f scale;
};

class PlaneSphereHandler
{
public:
	PlaneSphereHandler( const vl::Vec3f& n, float d )
		: normal_(n), d_(d) {}


	typedef std::pair<vl::Vec3f, float> Collision;

	bool operator()( const vl::Vec3f& point, float radius ) const
	{
		float value = dot(normal_, point) - d_;
		return( value < radius );

		return true;
	}

	bool operator()( const vl::Vec3f& point )
	{
		float value = dot(normal_, point) - d_;
		if( value < 0.0f )
			collisions_.push_back( Collision(point, -value) );

		// never want it to exit early:
		return false;
	}

	typedef std::deque<Collision>::const_iterator collision_iter;
	collision_iter collisions_begin() const { return collisions_.begin(); }
	collision_iter collisions_end() const   { return collisions_.end();   }

private:
	std::deque<Collision> collisions_;

	vl::Vec3f normal_;
	float d_;
};

int dCollideMeshPlane(dGeomID o1, dGeomID o2, int flags,
              dContactGeom *contact, int skip)
{
	assert(skip >= (int)sizeof(dContactGeom));
	assert(dGeomGetClass(o1) == dMyTriMeshClass);
	assert(dGeomGetClass(o2) == dPlaneClass);

	dVector4 planeParams;
	dGeomPlaneGetParams( o2, planeParams );
	vl::Vec3f normal( planeParams[0], planeParams[1], planeParams[2] );
	float d = planeParams[3];
	{
		float nLen = len( normal );
		normal /= nLen;
		d /= nLen;
	}

	/*
	// first, transform the plane according to the appropriate
	//   translation:

	const dReal* bodyPosV = dGeomGetPosition(o1);
	vl::Vec3f bodyPos( bodyPosV[0], bodyPosV[1], bodyPosV[2] );
	d -= dot( bodyPos, normal );

	const dReal* bodyRot = dGeomGetRotation(o1);
	vl::Mat3f rotMat( 
		bodyRot[0], bodyRot[1], bodyRot[2], 
		bodyRot[4], bodyRot[5], bodyRot[6], 
		bodyRot[8], bodyRot[9], bodyRot[10] );

	ODETriMeshStructure* data = 
		reinterpret_cast<ODETriMeshStructure*>(dGeomGetClassData(o1));
	vl::Mat3f triMeshScale( data->scale[0], 0.0, 0.0,
		0.0, data->scale[1], 0.0,
		0.0, 0.0, data->scale[2] );
	vl::Mat3f normalTransform = triMeshScale * trans(rotMat);

	vl::Vec3f newNormal = normalTransform * normal;

	// now, need to re-normalize
    float nLen = len(newNormal);
	newNormal /= nLen;
	d /= nLen;

	// Okay, we have our modified plane equation.  Now, test against
	//   all the points in the tree.
	PlaneSphereHandler handler( newNormal, d );
	data->mesh->sphereTree().query( handler );

	size_t iCollision = 0;
	for( PlaneSphereHandler::collision_iter iter = handler.collisions_begin();
		iter != handler.collisions_end(); ++iter )
	{
		if( iCollision == flags )
			return iCollision;

		contact[iCollision].depth = iter->second;

		vl::Vec3f pos = iter->first;
		vl::Vec3f globalPos = rotMat*(triMeshScale*pos);
		globalPos += bodyPos;
		for( size_t i = 0; i < 3; ++i )
			contact[iCollision].pos[i] = globalPos[i];

		for( size_t i = 0; i < 3; ++i )
			contact[iCollision].normal[i] = normal[i];

		contact[iCollision].g1 = o1;
		contact[iCollision].g2 = o2;

		++iCollision;
	}

	return iCollision;
	*/

	return 0;
}

/*
int dCollideSpherePlane (dxGeom *o1, dxGeom *o2, int flags,
			 dContactGeom *contact, int skip)
{
  assert(skip >= (int)sizeof(dContactGeom));
  assert(o1->type == dMyTriMeshClass);
  assert(o2->type == dPlaneClass);
  dxSphere *sphere = (dxSphere*) o1;
  dxPlane *plane = (dxPlane*) o2;

  contact->g1 = o1;
  contact->g2 = o2;
  dReal k = dDOT (o1->pos,plane->p);
  dReal depth = plane->p[3] - k + sphere->radius;
  if (depth >= 0) {
    contact->normal[0] = plane->p[0];
    contact->normal[1] = plane->p[1];
    contact->normal[2] = plane->p[2];
    contact->pos[0] = o1->pos[0] - plane->p[0] * sphere->radius;
    contact->pos[1] = o1->pos[1] - plane->p[1] * sphere->radius;
    contact->pos[2] = o1->pos[2] - plane->p[2] * sphere->radius;
    contact->depth = depth;
    return 1;
  }
  else return 0;
}
*/

dColliderFn* getMeshColliderFn(int num)
{
	switch( num )
	{
	case dPlaneClass:
		return &dCollideMeshPlane;
	default:
		return 0;
	}
}

void dMyTriMeshDestructor(dGeomID id)
{
	ODETriMeshStructure* data = 
		reinterpret_cast<ODETriMeshStructure*>(dGeomGetClassData(id));
	data->~ODETriMeshStructure();
}

void initMyTriMeshClass()
{
	dGeomClass myTriMeshClass;
	myTriMeshClass.aabb = dInfiniteAABB;
	myTriMeshClass.aabb_test = 0;
	myTriMeshClass.bytes = sizeof(ODETriMeshStructure);
	myTriMeshClass.collider = &getMeshColliderFn;
	myTriMeshClass.dtor = &dMyTriMeshDestructor;

	dMyTriMeshClass = dCreateGeomClass( &myTriMeshClass );
}

boost::once_flag initTriMeshClass_onceFlag = BOOST_ONCE_INIT;

dGeomID dCreateMyTriMesh( dSpaceID space, boost::shared_ptr<TriangleMesh> m, const vl::Vec3f& s )
{
	boost::call_once( &initMyTriMeshClass, initTriMeshClass_onceFlag );

	dGeomID id = dCreateGeom( dMyTriMeshClass );
	void* data = dGeomGetClassData(id);
	ODETriMeshStructure* triMesh = new (data) ODETriMeshStructure(m, s);

	dSpaceAdd(space, id);

	return id;
}

const float jointScale = 0.1;

Joint::Joint( const std::string& name, 
			 const StaticTriangleMesh* leftMeshObj, 
			 const StaticTriangleMesh* rightMeshObj )
	:	TransformedSceneGraphElement(name),
		disableAfterTime_( boost::numeric::bounds<float>::highest() )
{
	this->meshes_.reserve(2);
    if( leftMeshObj )
    	this->meshes_.push_back( leftMeshObj );

    if( rightMeshObj )
    	this->meshes_.push_back( rightMeshObj );
}

bool Joint::equalsLocal( const SceneGraphElement& other ) const
{
	const Joint* joint = dynamic_cast<const Joint*>( &other );
	assert( joint != 0 );

	assert( this->numControllableAxes() == joint->numControllableAxes() );
	for( size_t i = 0; i < this->numControllableAxes(); ++i )
	{
		if( !withinEpsilon(this->minAngle_[i], joint->minAngle_[i]) ||
			!withinEpsilon(this->maxAngle_[i], joint->maxAngle_[i]) ||
			!withinEpsilon(this->cfm_[i], joint->cfm_[i]) ||
			!withinEpsilon(this->jointFriction_[i], joint->jointFriction_[i]) ||
			!withinEpsilon(this->goalAngle_[i], joint->goalAngle_[i]) ||
			!withinEpsilon(this->maxForce_[i], joint->maxForce_[i]) ||
			!withinEpsilon(this->disableAfterTime_, joint->disableAfterTime_) )
			return false;
	}

	// We no longer check the linked objects for equality since this would be duplicating effort
	//   (and, if we checked joints for equality in the objects, would cause infinite loops). 
	//   Therefore it is up to the calling routine to take care of fixing up connections.

	return TransformedSceneGraphElement::equalsLocal( other );
}

void Joint::setObjects( PhysicsObjectPtr object )
{
    assert( this->numObjects() == 1 );

    this->clearObjects();
    objects_.push_back( object );
    object->addTransformChangedTarget( this );

	this->updateTransformLocal();
}

void Joint::setObjects( PhysicsObjectPair objects )
{
    assert( this->numObjects() == 2 );

    this->clearObjects();

	objects_.resize( 2 );
	objects_[0] = objects.first;
	objects_[1] = objects.second;

	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
		objects_[iObject].lock()->addTransformChangedTarget( this );

	this->updateTransformLocal();
}

void Joint::setObjects( const std::vector<PhysicsObjectPtr>& objects )
{
    assert( objects.size() == this->numObjects() );
    this->clearObjects();

    objects_.reserve( objects.size() );
    std::copy( objects.begin(), objects.end(), 
        std::back_inserter(objects_) );

	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
		objects_[iObject].lock()->addTransformChangedTarget( this );

	this->updateTransformLocal();
}

void Joint::clearObjects()
{
	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
	{
		PhysicsObjectPtr object = objects_[iObject].lock();
		if( object )
			object->removeTransformChangedTarget( this );
	}

    objects_.clear();
}

void Joint::init()
{
	this->minAngle_.resize( this->numControllableAxes() );
	this->maxAngle_.resize( this->numControllableAxes() );
	this->cfm_.resize( this->numControllableAxes() );
	this->jointFriction_.resize( this->numControllableAxes() );
	this->goalAngle_.resize( this->numControllableAxes() );
	this->maxForce_.resize( this->numControllableAxes() );

	for( size_t iAxis = 0; iAxis < numControllableAxes(); ++iAxis )
	{
		this->minAngle_[iAxis] = -dInfinity;
		this->maxAngle_[iAxis] = dInfinity;
		this->cfm_[iAxis] = 1e-5;
		this->jointFriction_[iAxis] = 0.0f;
		this->goalAngle_[iAxis] = 0.0f;
		this->maxForce_[iAxis] = 0.0f;
	}

	this->setScale( vl::Vec3f(jointScale, jointScale, jointScale) );
}

std::vector<ConstPhysicsObjectPtr> Joint::objects() const
{
	std::vector<ConstPhysicsObjectPtr> result;
	for( std::vector< boost::weak_ptr<PhysicsObject> >::const_iterator itr = this->objects_.begin();
		itr != this->objects_.end(); ++itr )
	{
		result.push_back( itr->lock() );
	}

	return result;
}

const char* JOINT_FRICTION_STRING = "jointFriction";
const char* CFM_STRING = "cfm";
const char* DISABLE_AFTER_TIME_STRING = "disableAfterTime";
const char* HIGH_STOP_STRING = "highStop";
const char* LOW_STOP_STRING = "lowStop";
const char* GOAL_ANGLE_STRING = "goalAngle";
const char* MAX_FORCE_STRING = "maxForce";

boost::shared_ptr<MechModelDeclaredObject> Joint::toMechModelLocal(MechModelType type) const
{
	boost::shared_ptr<MechModelDeclaredObject> object = TransformedSceneGraphElement::toMechModelLocal(type);
	// joints have nothing to do with audio
	if( type == MECH_MODEL_AUDIO || type == MECH_MODEL_AUDIO_ROOT )
		return object;

    assert( !this->objects_.empty() );
	object->add( mechModelValuePair("first", this->objects_[0].lock()->name()) );
    if( this->objects_.size() > 1 )
    	object->add( mechModelValuePair("second", this->objects_[1].lock()->name()) );

	for( size_t iAxis = 0; iAxis < numControllableAxes(); ++iAxis )
	{
		std::string number = boost::lexical_cast<std::string>(iAxis+1);
		if( this->cfm_[iAxis] != 1e-5 )
			object->add( mechModelValuePair(CFM_STRING + number, this->cfm_[iAxis]) );
		if( this->jointFriction_[iAxis] != 0.0f )
			object->add( mechModelValuePair(JOINT_FRICTION_STRING + number, this->jointFriction_[iAxis]) );
		if( this->minAngle_[iAxis] > -M_PI )
			object->add( mechModelValuePair(LOW_STOP_STRING + number, this->minAngle_[iAxis]) );
		if( this->maxAngle_[iAxis] < M_PI )
			object->add( mechModelValuePair(HIGH_STOP_STRING + number, this->maxAngle_[iAxis]) );
		if( this->goalAngle_[iAxis] != 0.0f || this->maxForce_[iAxis] != 0.0f )
			object->add( mechModelValuePair(GOAL_ANGLE_STRING + number, this->goalAngle_[iAxis]) );
		if( this->maxForce_[iAxis] != 0.0f )
			object->add( mechModelValuePair(MAX_FORCE_STRING + number, this->maxForce_[iAxis]) );
	}
	if( this->disableAfterTime_ < boost::numeric::bounds<float>::highest() )
		object->add( mechModelValuePair(DISABLE_AFTER_TIME_STRING, this->disableAfterTime_) );
	return object;
}

std::vector<std::string> Joint::jointAttributes() const
{
	std::vector<std::string> result;
	result.reserve( numControllableAxes() * 6 );
	for( size_t iAxis = 0; iAxis < numControllableAxes(); ++iAxis )
	{
		std::string number = boost::lexical_cast<std::string>(iAxis+1);
		result.push_back( CFM_STRING + number );
		result.push_back( JOINT_FRICTION_STRING + number );
		result.push_back( HIGH_STOP_STRING + number );
		result.push_back( LOW_STOP_STRING + number );
		result.push_back( GOAL_ANGLE_STRING + number );
		result.push_back( MAX_FORCE_STRING + number );
	}
	result.push_back( DISABLE_AFTER_TIME_STRING );
	return result;
}


std::vector<std::string> Joint::floatAttributes() const
{
	std::vector<std::string> result = TransformedSceneGraphElement::floatAttributes();
	std::vector<std::string> localAtts = this->jointAttributes();
	result.reserve( result.size() + localAtts.size() );
	std::copy( localAtts.begin(), localAtts.end(), std::back_inserter(result) );
	return result;
}

void Joint::setAttribute( const std::string& attributeName, float value )
{
	assert( !attributeName.empty() );
	char lastChar = attributeName[ attributeName.size() - 1 ];
	if( attributeName == DISABLE_AFTER_TIME_STRING )
		this->disableAfterTime_ = value;
	else if( lastChar >= '1' && lastChar < ('1' + numControllableAxes()) )
	{
		int iAxis = lastChar - '1';
		std::string prefix = attributeName.substr( 0, attributeName.size() - 1 );
		if( prefix == CFM_STRING )
			this->cfm_[iAxis] = value;
		else if( prefix == JOINT_FRICTION_STRING )
			this->jointFriction_[iAxis] = value;
		else if( prefix == HIGH_STOP_STRING )
			this->maxAngle_[iAxis] = value;
		else if( prefix == LOW_STOP_STRING )
			this->minAngle_[iAxis] = value;
		else if( prefix == GOAL_ANGLE_STRING )
			this->goalAngle_[iAxis] = value;
		else if( prefix == MAX_FORCE_STRING )
			this->maxForce_[iAxis] = value;
		else
			TransformedSceneGraphElement::setAttribute( attributeName, value );
	}
	else
		TransformedSceneGraphElement::setAttribute( attributeName, value );
}

float Joint::getAttribute( const std::string& attributeName ) const
{
	assert( !attributeName.empty() );
	char lastChar = attributeName[ attributeName.size() - 1 ];
	if( attributeName == DISABLE_AFTER_TIME_STRING )
		return this->disableAfterTime_;
	else if( lastChar >= '1' && lastChar < ('1' + numControllableAxes()) )
	{
		int iAxis = lastChar - '1';
		std::string prefix = attributeName.substr( 0, attributeName.size() - 1 );
		if( prefix == CFM_STRING )
			return this->cfm_[iAxis];
		else if( prefix == JOINT_FRICTION_STRING )
			return this->jointFriction_[iAxis];
		else if( prefix == HIGH_STOP_STRING )
			return this->maxAngle_[iAxis];
		else if( prefix == LOW_STOP_STRING )
			return this->minAngle_[iAxis];
		else if( prefix == GOAL_ANGLE_STRING )
			return this->goalAngle_[iAxis];
		else if( prefix == MAX_FORCE_STRING )
			return this->maxForce_[iAxis];
		else
			return TransformedSceneGraphElement::getAttribute( attributeName );
	}
	else
		return TransformedSceneGraphElement::getAttribute( attributeName );
}

void Joint::updateTransformLocal()
{
	if( objects_.empty() )
		return;

	this->objectFrameToLocal_.resize( objects_.size() );
	vl::Mat4f myTransform = this->transform();
	for( size_t iObject = 0; iObject < objects_.size(); ++iObject )
	{
		vl::Mat4f objectTransform = this->objects_[iObject].lock()->transform();
		vl::Mat4f localJointTransform = this->localJointTransform(iObject, myTransform, 0);
		vl::Mat4f xform = vl::inv(objectTransform) * (myTransform * localJointTransform);
		this->objectFrameToLocal_[iObject] = xform;
	}

	renderedMeshes_ = meshes_;
	assert( meshes_.size() == this->numObjects() && 
        objects_.size() == this->numObjects() );

	/*
	vl::Vec3f pos1 = this->objectPositionInLocalFrame(0, myTransform, 0);
	vl::Vec3f pos2 = this->objectPositionInLocalFrame(1, myTransform, 0);

	if( pos1[2] < pos2[2] )
		std::swap( renderedMeshes_.at(0), renderedMeshes_.at(1) );
		*/
}

vl::Vec3f Joint::objectPositionInLocalFrame( size_t iObject, 
	const vl::Mat4f& currentTransform, const char* statePtr ) const
{
	if( statePtr != 0 )
	{
		typedef const char* StatePtr;
		const char* currentPtr = statePtr + iObject*sizeof(const char*);
		const StatePtr* currentState = 
			reinterpret_cast<const StatePtr*>( currentPtr );
		vl::Mat4f objectTransform = objects_[iObject].lock()->transform( *currentState );
		vl::Vec3f objectPosGlobal = vl::xform( objectTransform, vl::Vec3f(vl::vl_0) );
		return vl::xform( inv( currentTransform ), objectPosGlobal );
	}
	else
	{
		vl::Mat4f objectTransform = this->objects_[iObject].lock()->transform();
		vl::Vec3f objectGlobalPos = vl::xform( objectTransform, vl::Vec3f(vl::vl_0) );
		vl::Vec3f objectLocalPos = vl::xform( vl::inv( currentTransform ), objectGlobalPos );

		return objectLocalPos;
	}
}

vl::Mat4f Joint::localJointTransform( size_t iObject,
	const vl::Mat4f& myTransform, const char* statePtr ) const
{
	if( statePtr )
	{
		typedef const char* StatePtr;
		const char* currentPtr = statePtr + iObject*sizeof(const char*);
		const StatePtr* currentState = 
			reinterpret_cast<const StatePtr*>( currentPtr );
		vl::Mat4f objectTransform = this->objects_[iObject].lock()->transform( *currentState );
		return objectTransform * this->objectFrameToLocal_[iObject];
	}

	/*
	vl::Vec3f pos = this->objectPositionInLocalFrame( iObject, myTransform, statePtr );
	if( iObject % 2 )
		pos *= -1.0f;

	vl::Mat4f transform = this->jointTransform( pos );
	return transform;
	*/

	return vl::vl_1;
}

Joint::~Joint()
{
	for( std::vector< boost::weak_ptr<PhysicsObject> >::iterator itr = this->objects_.begin();
		itr != this->objects_.end(); ++itr )
	{
		boost::shared_ptr<PhysicsObject> object = itr->lock();
		if( object )
			object->removeTransformChangedTarget(this);
	}
}

dJointID Joint::odeJoint(dWorldID worldId, dJointGroupID jointGroupId, const ODESimulation& sim) const
{
	dJointID id = this->createOdeJoint(worldId, jointGroupId);

    assert( !this->objects_.empty() && this->objects_.size() <= 2 );
    if( this->objects_.size() == 1 )
    {
	    boost::shared_ptr<PhysicsObject> first = this->objects_[0].lock();

	    dJointAttach( id, 
                    first->isStatic() ? 0 : sim.odeBodyID(first), 
                    0 );
    }
    else
    {
	    boost::shared_ptr<PhysicsObject> first = this->objects_[0].lock();
	    boost::shared_ptr<PhysicsObject> second = this->objects_[1].lock();

	    dJointAttach( id, 
                    first->isStatic() ? 0 : sim.odeBodyID(first), 
                    second->isStatic() ? 0 : sim.odeBodyID(second) );
    }

	this->setJointParameters(id);

	// set standard params
	for( size_t iAxis = 0; iAxis < numControllableAxes(); ++iAxis )
	{
		int offset = dParamGroup * iAxis;
		ODEParamFunction dJointSetParam = this->odeParamFunction();
		assert( this->minAngle_[iAxis] < this->maxAngle_[iAxis] );
		dJointSetParam( id, dParamLoStop + offset, this->minAngle_[iAxis] );
		dJointSetParam( id, dParamHiStop + offset, this->maxAngle_[iAxis] );
		assert( this->cfm_[iAxis] > 0.0f );
		dJointSetParam( id, dParamCFM + offset, this->cfm_[iAxis] );
		dJointSetParam( id, dParamVel + offset, 0.0f );
		dJointSetParam( id, dParamFMax, this->jointFriction_[iAxis] );
	}

	return id;
}

void Joint::setCFM( float value )
{
	for( size_t i = 0; i < this->cfm_.size(); ++i )
		this->cfm_[i] = value;
}

void Joint::setJointFriction( float value )
{
	for( size_t i = 0; i < this->jointFriction_.size(); ++i )
		this->jointFriction_[i] = value;
}

void Joint::updateControllers( dJointID odeJoint, float time, float timestep ) const
{
	/*
	const dReal maxVel = 4.0*M_PI;
	for( size_t iAxis = 0; iAxis < numControllableAxes(); ++iAxis )
	{
		if( this->maxForce_[iAxis] <= 0.0f )
			continue;

		ODEGetAngleFunction dJointGetAngle = this->odeGetAngleFunction(iAxis);
		ODEGetAngleRateFunction dJointGetAngleRate = this->odeGetAngleRateFunction(iAxis);

		dReal angle = dJointGetAngle(odeJoint);
		dReal rate = dJointGetAngleRate(odeJoint);

		dReal diff = this->goalAngle_[iAxis] - angle;
		dReal goalVel = 0.1 * diff / timestep;
		goalVel = std::min( std::max( -maxVel, goalVel ), maxVel );
		dReal velDiff = goalVel - rate;

		// the 0.9 is to relax it a bit
		int offset = dParamGroup * iAxis;
		ODEParamFunction dJointSetParam = this->odeParamFunction();
		dJointSetParam( odeJoint, dParamFMax + offset, this->maxForce_[iAxis] );
		dJointSetParam( odeJoint, dParamVel + offset, velDiff );
	}
	*/
}

#ifdef USE_NOVODEX
void Joint::initNxJointDesc(NxJointDesc& joint, const NxSimulation& sim) const
{
	ConstPhysicsObjectPtr object0 = this->objects_[0].lock();
	ConstPhysicsObjectPtr object1 = this->objects_[1].lock();

	if( object0->isStatic() )
		joint.actor[0] = 0;
	else
		joint.actor[0] = sim.nxActor( object0 ).get();

	if( object1->isStatic() )
		joint.actor[1] = 0;
	else
		joint.actor[1] = sim.nxActor( object1 ).get();
}
#endif

RibObjectPtr Joint::ribObject() const
{
	return RibObjectPtr();
}

float Joint::disableAfterTime() const
{
	return this->disableAfterTime_;
}

void Joint::setDisableAfterTime( float value )
{
	this->disableAfterTime_ = value;
}

BoundingBox3f Joint::bounds(const vl::Mat4f& transform) const
{
	BoundingBox3f result;
	// @todo this isn't strictly accurate, because there might be a local
	//   transform due to the object in question
	for( size_t i = 0; i < objects_.size(); ++i )
	{
		vl::Mat4f transform = this->localJointTransform( i, transform, 0 );
		result.expand( this->boundsMesh(this->renderedMeshes_[i], transform) );
	}

	return result;
}

vl::Mat4f Joint::localTransform( const char* state ) const
{
	if( state == 0 )
		return TransformedSceneGraphElement::localTransform(0);

	return vl::Mat4f( vl::vl_1 );
}

void Joint::renderLocalFill(float alpha, const char* statePtr) const
{
#ifdef GUI
	vl::Mat4f transform = this->transform( statePtr );

	for( size_t i = 0; i < objects_.size(); ++i )
	{
		vl::Vec4f color( 0.9, 0.9, 0.9, 1.0 );
		if( (i % 2) )
			color[2] = 0.0;

		glMaterialfv( GL_FRONT, GL_DIFFUSE, color.Ref() );

		GLMatrixStackHandler handler;
		vl::Mat4f jointTransform = vl::trans(  this->localJointTransform( i, transform, statePtr ) );
		glMultMatrixf( jointTransform.Ref() );
		this->renderedMeshes_[i]->displayList()->render();
	}
#endif
}

void Joint::renderLocalWire(const vl::Vec4f& wireframeColor, const char* statePtr) const
{
#ifdef GUI
	vl::Mat4f transform = this->transform( statePtr );

	/*
	// first draw the lines
	if( statePtr == 0 )
	{
		vl::Vec4f white( 0.9f, 0.9f, 0.9f, 1.0f );
		vl::Vec4f yellow( 0.9f, 0.9, 0.0f, 1.0f );

		{
			GLDisableHandler lighting( GL_LIGHTING );
			glLineWidth( 4.0f );

			GLActionHandler lines(GL_LINES);

			for( size_t i = 0; i < objects_.size(); ++i )
			{
				vl::Vec3f pos = this->objectPositionInLocalFrame(i, transform, statePtr);

				if( i == 0 )
					glColor3fv( white.Ref() );
				else
					glColor3fv( yellow.Ref() );

				glVertex3f(0.0f, 0.0f, 0.0f);
				glVertex3fv( pos.Ref() );
			}
		}
	}
	*/

	glDisable( GL_LIGHTING );
	glPolygonMode( GL_FRONT, GL_LINE );
	glColor4fv( wireframeColor.Ref() );
	glLineWidth( 1.0 );

	GLEnableClientStateHandler vertexArray(GL_VERTEX_ARRAY);

	for( size_t i = 0; i < objects_.size(); ++i )
	{
		GLMatrixStackHandler handler;
		glVertexPointer(3, GL_FLOAT, 0, renderedMeshes_[i]->positions());
		vl::Mat4f jointTransform = vl::trans(  this->localJointTransform( i, transform, statePtr ) );
		glMultMatrixf( jointTransform.Ref() );
		glDrawElements( GL_LINES, 2*this->renderedMeshes_[i]->numSegments(), GL_UNSIGNED_INT, 
			this->renderedMeshes_[i]->segments() );
	}

#endif
}

bool Joint::intersectLocal(const Ray3d& ray, double& t, const char* statePtr) const
{
	vl::Mat4f transform = this->transform( statePtr );

	double tMin = boost::numeric::bounds<double>::highest();
	bool intersected = false;
	for( size_t i = 0; i < objects_.size(); ++i )
	{
		double t;
		vl::Mat4f jointTransform = this->localJointTransform( i, transform, statePtr );

		if( this->intersectMesh(ray, t, this->renderedMeshes_[i], jointTransform ) )
		{
			tMin = std::min(tMin, t);
			intersected = true;
		}
	}

	if( intersected )
	{
		t = tMin;
		return true;
	}
	else
	{
		return false;
	}
}

bool Joint::intersectLocal(const BoundingBox2d& box, const ScreenSpaceConverter& converter, const char* statePtr) const
{
	vl::Mat4f transform = this->transform( statePtr );
	for( size_t i = 0; i < objects_.size(); ++i )
	{
		vl::Mat4f jointTransform = this->localJointTransform( i, transform, statePtr );
		if( this->intersectMesh(box, converter, this->renderedMeshes_[i], jointTransform) )
			return true;
	}

	return false;
}

bool Joint::intersectMesh(const Ray3d& ray, double& t, 
	const StaticTriangleMesh* mesh, vl::Mat4f transform ) const
{
	vl::Mat4d localToGlobal = toMat4d( transform );
	vl::Mat4d globalToLocal = inv( localToGlobal );

	vl::Vec3d pos = xform(globalToLocal, ray.getPosition());
	vl::Vec3d dir = xform(globalToLocal, ray.getPosition() + ray.getDirection()) - pos;
	float length = len( dir );
	dir /= length;

    Ray3d localRay( pos, dir );
	TriangleIntersectionHandler handler( mesh->positions() );
	if( mesh->triangleTree().intersectRay( localRay, handler ) )
	{
		t = handler.t();
		t /= length;
		return true;
	}
	else
	{
		return false;
	}
}

bool Joint::intersectMesh(const BoundingBox2d& box, const ScreenSpaceConverter& converter, 
	const StaticTriangleMesh* mesh, vl::Mat4f transform ) const
{
	vl::Mat4d localToGlobal = toMat4d( transform );
	WrappedScreenSpaceConverter wrappedConverter( converter, localToGlobal );

	BoxTriangleIntersectionHandler handler( box, wrappedConverter, mesh->positions() );
	return mesh->triangleTree().intersect( box, wrappedConverter, handler );
}

BoundingBox3f Joint::boundsMesh(const StaticTriangleMesh* mesh, vl::Mat4f transform) const
{
	return transformBounds( mesh->bounds(), transform );
}

BallAndSocketJoint::BallAndSocketJoint( const std::string& name )
	: Joint(name, &ballJointLeft_mesh, &ballJointRight_mesh)
{
	this->init();
}

BallAndSocketJoint* BallAndSocketJoint::cloneLocal() const
{
	return (new BallAndSocketJoint(*this) );
}

boost::shared_ptr<MechModelDeclaredObject> BallAndSocketJoint::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "ballJoint", this->name() ) );
	return result;
}

vl::Mat4f BallAndSocketJoint::jointTransform( const vl::Vec3f& direction ) const
{
	vl::Vec3f crossProd( vl::vl_0 );
	for( size_t i = 0; i < 3; ++i )
	{
		vl::Vec3f axisDir( vl::vl_0 );
		axisDir[i] = 1.0f;
		vl::Vec3f newCrossProd = vl::cross( direction, axisDir );
		if( vl::sqrlen(newCrossProd) > vl::sqrlen(crossProd) )
			crossProd = newCrossProd;
	}

	// @todo check this
	vl::Mat3f rotMat;
	rotMat[2] = vl::norm(direction);
	rotMat[0] = vl::norm(crossProd);
	rotMat[1] = vl::norm( vl::cross( rotMat[2], rotMat[0] ) );

	vl::Mat4f result( vl::vl_1 );
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			result[i][j] = rotMat[i][j];

	return trans(result);
}

#ifdef USE_NOVODEX
std::vector< boost::shared_ptr<NxJointDesc> > BallAndSocketJoint::getNxJointDesc(const NxSimulation& sim) const
{
	boost::shared_ptr<NxSphericalJointDesc> result( new NxSphericalJointDesc );
	initNxJointDesc( *result, sim );

	vl::Vec3f center = vl::xform( this->transform(), vl::Vec3f(vl::vl_0) );
	result->setGlobalAnchor( toNxVec(center) );

	return std::vector< boost::shared_ptr<NxJointDesc> >(1, result);
}
#endif

dJointID BallAndSocketJoint::createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const
{
	dJointID id = dJointCreateBall( worldId, jointGroupId );
	return id;
}

void BallAndSocketJoint::setJointParameters(dJointID id) const
{
	vl::Vec3f center = vl::xform( this->transform(), vl::Vec3f(vl::vl_0) );
	dJointSetBallAnchor( id, center[0], center[1], center[2] );

	// not allowed for ball joints

	//dJointSetBallParam( id, dParamCFM, this->CFM() );
	//dJointSetBallParam( id, dParamLoStop, this->getMinAngle() );
	//dJointSetBallParam( id, dParamHiStop, this->getMaxAngle() );
	
}

HingeJoint::HingeJoint( const std::string& name )
	: Joint( name, &hingeJointLeft_mesh, &hingeJointRight_mesh )
{
	this->init();
}

HingeJoint* HingeJoint::cloneLocal() const
{
	return (new HingeJoint(*this) );
}

vl::Mat4f HingeJoint::jointTransform( const vl::Vec3f& direction ) const
{
	vl::Vec3f crossProd( vl::vl_0 );
	for( size_t i = 0; i < 3; ++i )
	{
		vl::Vec3f axisDir( vl::vl_0 );
		axisDir[i] = 1.0f;
		vl::Vec3f newCrossProd = vl::cross( direction, axisDir );
		if( vl::sqrlen(newCrossProd) > vl::sqrlen(crossProd) )
			crossProd = newCrossProd;
	}

	// @todo check this
	vl::Mat3f rotMat;
	rotMat[1] = vl::Vec3f(0.0f, 1.0f, 0.0f);
	rotMat[2] = vl::norm(direction);
	rotMat[0] = vl::cross( rotMat[1], rotMat[2] );
	if( vl::sqrlen(rotMat[0]) < 1e-4 )
	{
		rotMat[1] = vl::Vec3f(1.0f, 0.0f, 0.0f);
		rotMat[0] = vl::cross( rotMat[1], rotMat[2] );
	}
	rotMat[0] = vl::norm( rotMat[0] );
	rotMat[2] = vl::norm( vl::cross( rotMat[0], rotMat[1] ) );

	vl::Mat4f result( vl::vl_1 );
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			result[i][j] = rotMat[i][j];

	return trans(result);
}

boost::shared_ptr<MechModelDeclaredObject> HingeJoint::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "hingeJoint", this->name() ) );
	return result;
}

#ifdef USE_NOVODEX
std::vector< boost::shared_ptr<NxJointDesc> > HingeJoint::getNxJointDesc(const NxSimulation& sim) const
{
	boost::shared_ptr<NxRevoluteJointDesc> result( new NxRevoluteJointDesc );
	initNxJointDesc( *result, sim );

	vl::Mat4f transform = this->transform();
	vl::Vec3f center = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	vl::Vec3f axis = vl::strip( transform * vl::Vec4f( 0.0f, 1.0f, 0.0f, 0.0f ) );

	result->setGlobalAxis( toNxVec(axis) );
	result->setGlobalAnchor( toNxVec(center) );

	return std::vector< boost::shared_ptr<NxJointDesc> >(1, result);
}
#endif

dJointID HingeJoint::createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const
{
	dJointID id = dJointCreateHinge( worldId, jointGroupId );
	return id;
}

void HingeJoint::setJointParameters(dJointID id) const
{
	vl::Mat4f transform = this->transform();
	vl::Vec3f center = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	vl::Vec3f axis = vl::strip( transform * vl::Vec4f( 0.0f, 1.0f, 0.0f, 0.0f ) );

	dJointSetHingeAnchor( id, center[0], center[1], center[2] );
	dJointSetHingeAxis( id, axis[0], axis[1], axis[2] );
}


UniversalJoint::UniversalJoint( const std::string& name )
	: Joint( name, &univJointLeft_mesh, &univJointRight_mesh )
{
	this->init();
}

UniversalJoint* UniversalJoint::cloneLocal() const
{
	return (new UniversalJoint(*this) );
}

boost::shared_ptr<MechModelDeclaredObject> UniversalJoint::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "universalJoint", this->name() ) );
	return result;
}

#ifdef USE_NOVODEX
std::vector< boost::shared_ptr<NxJointDesc> > UniversalJoint::getNxJointDesc(const NxSimulation& sim) const
{
	boost::shared_ptr<NxD6JointDesc> result( new NxD6JointDesc  );
	initNxJointDesc( *result, sim );

	vl::Mat4f transform = this->transform();
	vl::Vec3f center = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	vl::Vec3f axis1 = vl::norm( vl::strip( transform * vl::Vec4f(0.0f, 1.0f, 0.0f, 0.0f) ) );
	vl::Vec3f axis2 = vl::norm( vl::strip( transform * vl::Vec4f(1.0f, 0.0f, 0.0f, 0.0f) ) );
	vl::Vec3f rotateAxis = vl::norm( vl::cross( axis2, axis1 ) );

	result->twistMotion = NX_D6JOINT_MOTION_LOCKED;
	result->swing1Motion = NX_D6JOINT_MOTION_FREE;
	result->swing2Motion = NX_D6JOINT_MOTION_FREE;

	result->xMotion = NX_D6JOINT_MOTION_LOCKED;
	result->yMotion = NX_D6JOINT_MOTION_LOCKED;
	result->zMotion = NX_D6JOINT_MOTION_LOCKED;

	result->projectionMode = NX_JPM_NONE;

	result->setGlobalAxis( toNxVec(rotateAxis) );
	result->setGlobalAnchor( toNxVec(center) );

	return std::vector< boost::shared_ptr<NxJointDesc> >(1, result);
}
#endif

dJointID UniversalJoint::createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const
{
	dJointID id = dJointCreateUniversal( worldId, jointGroupId );
	return id;
}

void UniversalJoint::setJointParameters(dJointID id) const
{
	vl::Mat4f transform = this->transform();
	vl::Vec3f center = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	vl::Vec3f axis1 = vl::norm( vl::strip( transform * vl::Vec4f(0.0f, 1.0f, 0.0f, 0.0f) ) );
	vl::Vec3f axis2 = vl::norm( vl::strip( transform * vl::Vec4f(1.0f, 0.0f, 0.0f, 0.0f) ) );

	dJointSetUniversalAnchor( id, center[0], center[1], center[2] );

	// we can't set axis1 to this value or we will cause an error in ODE;
	// see the Wiki
	dReal defaultAxis1_tmp[4];
	dJointGetUniversalAxis1( id, defaultAxis1_tmp );
	vl::Vec3f defaultAxis1( 
		defaultAxis1_tmp[0], defaultAxis1_tmp[1], defaultAxis1_tmp[2] );
	// find a vector that is perpendicular to both the default axis2 and
	// our final axis2
	vl::Vec3f tempAxis2 = vl::cross( axis1, defaultAxis1 );
	float tempAxis2Len = vl::len( tempAxis2 );
	const float epsilon = 1e-3;
	if( tempAxis2Len < epsilon ) // in this case, axis1 and defaultAxis1 are in
	                             // the same direction, so defaultAxis2 must
	                             // be ok
		tempAxis2 = axis2;
	else
		tempAxis2 /= tempAxis2Len;

	dJointSetUniversalAxis2( id, tempAxis2[0], tempAxis2[1], tempAxis2[2] );
	dJointSetUniversalAxis1( id, axis1[0], axis1[1], axis1[2] );
	dJointSetUniversalAxis2( id, axis2[0], axis2[1], axis2[2] );
}

vl::Mat4f UniversalJoint::jointTransform( const vl::Vec3f& direction ) const
{
	return vl::Mat4f(vl::vl_1);
}

Hinge2Joint::Hinge2Joint( const std::string& name )
	: Joint( name, &hinge2JointLeft_mesh, &hinge2JointRight_mesh )
{
	this->init();
}

Hinge2Joint* Hinge2Joint::cloneLocal() const
{
	return (new Hinge2Joint(*this) );
}

boost::shared_ptr<MechModelDeclaredObject> Hinge2Joint::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "hinge2Joint", this->name() ) );
	return result;
}

#ifdef USE_NOVODEX
// @todo implement this
std::vector< boost::shared_ptr<NxJointDesc> > Hinge2Joint::getNxJointDesc(const NxSimulation& sim) const
{
	std::vector< boost::shared_ptr<NxJointDesc> > result;
	return result;
}
#endif


dJointID Hinge2Joint::createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const
{
	dJointID id = dJointCreateHinge2( worldId, jointGroupId );
	return id;
}

void Hinge2Joint::setJointParameters(dJointID id) const
{
	vl::Mat4f transform = this->transform();
	vl::Vec3f center = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	vl::Vec3f axis1 = vl::norm( vl::strip( transform * vl::Vec4f(0.0f, 1.0f, 0.0f, 0.0f) ) );
	vl::Vec3f axis2 = vl::norm( vl::strip( transform * vl::Vec4f(0.0f, 0.0f, 1.0f, 0.0f) ) );
	dJointSetHinge2Anchor( id, center[0], center[1], center[2] );

	// we can't set axis1 to this value or we will cause an error in ODE;
	// see the Wiki
	dReal defaultAxis1_tmp[4];
	dJointGetHinge2Axis1( id, defaultAxis1_tmp );
	vl::Vec3f defaultAxis1( 
		defaultAxis1_tmp[0], defaultAxis1_tmp[1], defaultAxis1_tmp[2] );
	// find a vector that is perpendicular to both the default axis2 and
	// our final axis2
	vl::Vec3f tempAxis2 = vl::cross( axis1, defaultAxis1 );
	float tempAxis2Len = vl::len( tempAxis2 );
	const float epsilon = 1e-3;
	if( tempAxis2Len < epsilon )
		tempAxis2 = axis2;
	else
		tempAxis2 /= tempAxis2Len;

	dJointSetHinge2Axis2( id, tempAxis2[0], tempAxis2[1], tempAxis2[2] );
	dJointSetHinge2Axis1( id, axis1[0], axis1[1], axis1[2] );
	dJointSetHinge2Axis2( id, axis2[0], axis2[1], axis2[2] );
}


vl::Mat4f Hinge2Joint::jointTransform( const vl::Vec3f& direction ) const
{
	return vl::Mat4f(vl::vl_1);
}

vl::Mat4f Hinge2Joint::middleTransform(const char* statePtr) const
{
	vl::Mat4f transform = this->transform( statePtr );
	vl::Mat4f firstTransform = 
		this->localJointTransform( 0, transform, statePtr );
	vl::Mat4f secondTransform = 
		this->localJointTransform( 1, transform, statePtr );

	// need to rotate the second y vector up so that it matches the first
	vl::Vec3f yFirst = vl::norm( strip(firstTransform * vl::Vec4f(0.0f, 1.0f, 0.0f, 0.0f)) );
	vl::Vec3f ySecond = vl::norm( strip(secondTransform * vl::Vec4f(0.0f, 1.0f, 0.0f, 0.0f)) );

	vl::Mat3f upperLeft;
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			upperLeft[i][j] = secondTransform[i][j];

	vl::Mat3f rotMat = fromToRotation( ySecond, yFirst ) * upperLeft;
	for( vl::Int i = 0; i < 3; ++i )
		for( vl::Int j = 0; j < 3; ++j )
			secondTransform[i][j] = rotMat[i][j];

	return secondTransform;
}

void Hinge2Joint::renderLocalFill(float alpha, const char* statePtr) const
{
#ifdef GUI
	Joint::renderLocalFill( alpha, statePtr );

	vl::Vec4f color( 0.9, 0.9, 0.9, 1.0 );
	glMaterialfv( GL_FRONT, GL_DIFFUSE, color.Ref() );


	vl::Mat4f transform = this->middleTransform(statePtr);

	GLMatrixStackHandler handler;
	vl::Mat4f glTransform = vl::trans( transform );
	glMultMatrixf( glTransform.Ref() );
	hinge2JointMid_mesh.displayList()->render();
#endif
}

void Hinge2Joint::renderLocalWire(
	const vl::Vec4f& wireframeColor, 
	const char* statePtr) const
{
#ifdef GUI
	Joint::renderLocalWire( wireframeColor, statePtr );

	vl::Mat4f transform = this->middleTransform(statePtr);

	GLMatrixStackHandler handler;
	vl::Mat4f glTransform = vl::trans( transform );
	glMultMatrixf( glTransform.Ref() );
	hinge2JointMid_mesh.displayList()->render();

	glDisable( GL_LIGHTING );
	glColor4fv( wireframeColor.Ref() );
	glLineWidth( 1.0 );

	GLEnableClientStateHandler vertexArray(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, hinge2JointMid_mesh.positions());
	glDrawElements( GL_LINES, 2*hinge2JointMid_mesh.numSegments(), GL_UNSIGNED_INT, 
		hinge2JointMid_mesh.segments() );
#endif
}

SliderJoint::SliderJoint( const std::string& name )
	: Joint(name, &sliderJointLeft_mesh, &sliderJointRight_mesh)
{
	this->init();
}

SliderJoint* SliderJoint::cloneLocal() const
{
	return (new SliderJoint(*this) );
}

boost::shared_ptr<MechModelDeclaredObject> SliderJoint::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "sliderJoint", this->name() ) );
	return result;
}

#ifdef USE_NOVODEX
std::vector< boost::shared_ptr<NxJointDesc> > SliderJoint::getNxJointDesc(const NxSimulation& sim) const
{
	boost::shared_ptr<NxPrismaticJointDesc> result( new NxPrismaticJointDesc );
	initNxJointDesc( *result, sim );

	vl::Mat4f transform = this->transform();
	vl::Vec3f center = vl::xform( transform, vl::Vec3f(vl::vl_0) );
	result->setGlobalAnchor( toNxVec(center) );

	vl::Vec3f axis = vl::strip( transform * vl::Vec4f( 0.0f, 0.0f, 1.0f, 0.0f ) );
	result->setGlobalAxis( toNxVec(axis) );

	return std::vector< boost::shared_ptr<NxJointDesc> >(1, result);
}
#endif

dJointID SliderJoint::createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const
{
	dJointID id = dJointCreateSlider( worldId, jointGroupId );
	return id;
}

void SliderJoint::setJointParameters(dJointID id) const
{
	vl::Mat4f transform = this->transform();
	vl::Vec3f axis = vl::strip( transform * vl::Vec4f( 0.0f, 0.0f, 1.0f, 0.0f ) );

	dJointSetSliderAxis( id, axis[0], axis[1], axis[2] );
}

vl::Mat4f SliderJoint::jointTransform( const vl::Vec3f& direction ) const
{
	return vl::Mat4f(vl::vl_1);
}

PlaneJoint::PlaneJoint( const std::string& name )
    : Joint( name, &planeJoint_mesh  )
{
}

PlaneJoint* PlaneJoint::cloneLocal() const
{
	return (new PlaneJoint(*this) );
}

dJointID PlaneJoint::createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const
{
	dJointID id = dJointCreatePlane2D( worldId, jointGroupId );
	return id;
}

void PlaneJoint::setJointParameters(dJointID id) const
{
    // Nothing to do
}

boost::shared_ptr<MechModelDeclaredObject> PlaneJoint::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "planeJoint", this->name() ) );
	return result;
}

/*
OneAxisRotationalMotor::OneAxisRotationalMotor( const std::string& name )
	: Joint(name, &rotMotorLeft_mesh, &rotMotorRight_mesh)
{
}

vl::Mat4f OneAxisRotationalMotor::jointTransform( const vl::Vec3f& direction ) const
{
	return vl::vl_1;
}

OneAxisRotationalMotor* OneAxisRotationalMotor::cloneLocal() const
{
	return new OneAxisRotationalMotor( *this );
}

dJointID OneAxisRotationalMotor::createOdeJoint(dWorldID worldId, dJointGroupID jointGroupId) const
{
	dJointID id = dJointCreateAMotor( worldId, jointGroupId );
	return id;
}

void OneAxisRotationalMotor::setJointParameters(dJointID id) const
{
	dJointSetAMotorMode( id, dAMotorUser );
	dJointSetAMotorNumAxes( id, 1 );

	vl::Mat4f transform = this->transform();
	vl::Vec3f axis = vl::strip( transform * vl::Vec4f( 0.0f, 1.0f, 0.0f, 0.0f ) );

	dJointSetAMotorAxis( id, 0, 1, axis[0], axis[1], axis[2] );
	dJointSetAMotorAngle( id, 0, 0 );

	dJointSetAMotorParam( id, dParamVel, 0 );
}

boost::shared_ptr<MechModelDeclaredObject> OneAxisRotationalMotor::mechModelObject() const
{
	boost::shared_ptr<MechModelDeclaredObject> result(
		new MechModelDeclaredObject( "rotMotor1", this->name() ) );
	return result;
}
*/

} // namespace planning
