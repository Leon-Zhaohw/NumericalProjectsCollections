#include "stdafx.h"

#ifdef __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include <iostream>
#include <cassert>

#include "twigg/gslWrappers.h"
#include "twigg/camera.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

using namespace std;
using namespace vl;

const float kMouseRotationSensitivity		= 1.0f/90.0f;
const float kMouseTranslationXSensitivity	= 0.01f;
const float kMouseTranslationYSensitivity	= 0.01f;
const float kMouseZoomSensitivity			= 0.06f;

/*
const float kMouseTranslationXSensitivity	= 0.003f;
const float kMouseTranslationYSensitivity	= 0.003f;
const float kMouseZoomSensitivity			= 0.008f;
*/

Camera::Camera()
	: zMin_( 0.1 ), zMax_( 3000.0 ), ratio_(0.0), sensitivity_(1.0), basis_(vl::vl_1), scale_(1.0)
{
}

Camera::~Camera()
{
}

void Camera::setTranslationalSensitivity(float value)
{
	sensitivity_ = value;
}

float Camera::sensitivity() const
{
	return sensitivity_;
}

void Camera::setZRange( double zMin, double zMax )
{
	assert( zMax > zMin );
	assert( zMin > 0.0 );

	zMin_ = zMin;
	zMax_ = zMax;
}

double Camera::scale() const
{
	return this->scale_;
}

void Camera::setScale( double scale )
{
	this->scale_ = scale;
}

double Camera::zMin() const
{
	return zMin_;
}

double Camera::zMax() const
{
	return zMax_;
}

double Camera::ratio() const
{
	if( ratio_ == 0.0 )
	{
		double viewport[4];
		glGetDoublev(GL_VIEWPORT, viewport);
		const double viewportWidth = viewport[2];
		const double viewportHeight = viewport[3];
		return viewportWidth / viewportHeight;
	}
	else
		return ratio_;
}

void Camera::setUpVector(const Vec3f& up)
{
	Vec3f y = up;

	// We need a basis, first of all.
	basis_[1] = y;
	basis_[1] /= len(basis_[1]);

	float maxLen = 0.0;
	for( int i=0; i < 3; ++i )
	{
		Vec3f b = vl_0;
		b[i] = 1.0;
		Vec3f z = cross(b, y);
		if( len(z) > maxLen )
			basis_[2] = z / len(z);
	}

	basis_[0] = cross( basis_[1], basis_[2] );
	basis_[0] /= len(basis_[0]);
}

vl::Vec3f Camera::upVector() const
{
	return basis_[1];
}

vl::Mat3f Camera::basis() const
{
	return basis_;
}

void Camera::setBasis(const vl::Mat3f& basis)
{
	basis_ = basis;
}

vl::Vec3f ProjectiveFixedYCamera::getPosition() const
{
	vl::Mat4f dollyXform;
	vl::Mat4f azimXform;
	vl::Mat4f elevXform;
	vl::Mat4f twistXform;
	vl::Mat4f originXform;

	vl::Mat3f basis = Camera::basis();
	dollyXform.MakeHTrans(mDolly * basis[2]);
	azimXform.MakeHRot(basis[1], mAzimuth);
	elevXform.MakeHRot(basis[0], mElevation);
	twistXform = vl_1;
	originXform.MakeHTrans(mLookAt);

	// grouped for (mat4 * vec3) ops instead of (mat4 * mat4) ops
	Vec4f pos(0.0, 0.0, 0.0, 1.0);
	pos = originXform * (azimXform * (elevXform * (dollyXform * pos)));
	return Vec3f(pos[0], pos[1], pos[2]);
}

vl::Vec3f ProjectiveFixedYCamera::effectiveUpVector() const
{
	vl::Mat3f basis = Camera::basis();
	if ( fmod(mElevation, static_cast<float>(2.0*vl_pi)) < 3*vl_halfPi && 
		fmod(mElevation, static_cast<float>(2.0*vl_pi)) > vl_halfPi )
		return -basis[1];
	else
		return basis[1];
}

ProjectiveFixedYCamera::ProjectiveFixedYCamera(bool fixedCenter)
	: fixedCenter_( fixedCenter )
{
	mElevation = mAzimuth = mTwist = 0.0f;
	mDolly = -5.0f;
	mElevation = 0.2f;
	mAzimuth = (float)vl_pi;

	mLookAt = Vec3f( 0, 0, 0 );
	mCurrentMouseAction = kActionNone;

	fov_ = 40.0;
}

void ProjectiveFixedYCamera::frame( const BoundingBox3f& box )
{
	if( box.empty() )
		return;


	vl::Vec3f boxCenter = 0.5*(box.minimum() + box.maximum());
	vl::Vec3f camPos = this->getPosition();

	if( this->fixedCenter_ )
		this->setLookAt( boxCenter );

	vl::Vec3f dir = boxCenter - camPos;
	float dist = vl::len(dir);
	dir /= dist;

	// make sure each corner of the box is in view
	float minDist = 0.0f;
	{

		const vl::Vec3f mn = box.minimum();
		const vl::Vec3f mx = box.maximum();
		for( int i = 0; i < 2; ++i )
			for( int j = 0; j < 2; ++j )
				for( int k = 0; k < 2; ++k )
				{
					vl::Vec3f point( 
						i == 0 ? mn[0] : mx[0],
						j == 0 ? mn[1] : mx[1],
						k == 0 ? mn[2] : mx[2] );

					vl::Vec3f vec = point - camPos;
					vl::Vec3f perpendicularVec = vec - dot(vec, dir)*dir;
					// shouldn't need sqrt(2) here, compensating for something
					minDist = std::max<float>( minDist, 
						sqrt(2.0)*vl::len(perpendicularVec) / tan(0.5*(M_PI/180.0)*this->FOV()) );
				}


	}

	if( this->fixedCenter_ )
	{
		this->setDolly( signum(this->getDolly()) * minDist );
	}
	else
	{
		float actualDist = minDist + this->getDolly();
		this->setLookAt( boxCenter - actualDist * dir );
	}
}


void ProjectiveFixedYCamera::dumpParametersBinary( std::ostream& os ) const
{
	vl::Mat3f basis = Camera::basis();

	os.write( reinterpret_cast<const char*>(&mElevation), sizeof(float) );
	os.write( reinterpret_cast<const char*>(&mAzimuth), sizeof(float) );
	os.write( reinterpret_cast<const char*>(&mTwist), sizeof(float) );
	os.write( reinterpret_cast<const char*>(&mDolly), sizeof(float) );
	os.write( reinterpret_cast<const char*>(&fov_), sizeof(float) );
	os.write( reinterpret_cast<const char*>(&mLookAt), sizeof(vl::Vec3f) );
	os.write( reinterpret_cast<const char*>(&basis), sizeof(vl::Mat3f) );
}

void ProjectiveFixedYCamera::readParametersBinary( std::istream& is )
{
	vl::Mat3f basis;

	is.read( reinterpret_cast<char*>(&mElevation), sizeof(float) );
	is.read( reinterpret_cast<char*>(&mAzimuth), sizeof(float) );
	is.read( reinterpret_cast<char*>(&mTwist), sizeof(float) );
	is.read( reinterpret_cast<char*>(&mDolly), sizeof(float) );
	is.read( reinterpret_cast<char*>(&fov_), sizeof(float) );
	is.read( reinterpret_cast<char*>(&mLookAt), sizeof(vl::Vec3f) );
	is.read( reinterpret_cast<char*>(&basis), sizeof(vl::Mat3f) );

	this->setBasis( basis );
}

void ProjectiveFixedYCamera::dumpParameters( std::ostream& os ) const
{
	vl::Mat3f basis = Camera::basis();

	os << mElevation << " " << mAzimuth << " " << mTwist << " " << mDolly << " " << fov_ << std::endl;
	os << mLookAt << std::endl;
	os << basis;
}

void ProjectiveFixedYCamera::readParameters( std::istream& is )
{
	vl::Mat3f basis;

	is >> mElevation >> mAzimuth >> mTwist >> mDolly >> fov_;
	is >> mLookAt;
	is >> basis;

	this->setBasis( basis );
}

void ProjectiveFixedYCamera::clickMouse( MouseAction_t action, int x, int y )
{
	mCurrentMouseAction = action;
	mLastMousePosition[0] = x;
	mLastMousePosition[1] = y;
}

void ProjectiveFixedYCamera::dragMouse( int x, int y )
{
	Vec3f mouseDelta   = Vec3f(x,y,0.0f) - mLastMousePosition;
	mLastMousePosition = Vec3f(x,y,0.0f);

	switch(mCurrentMouseAction)
	{
	case kActionTranslate:
		{
			float xTrack =  -mouseDelta[0] * sensitivity() * kMouseTranslationXSensitivity;
			float yTrack =  mouseDelta[1] * sensitivity() * kMouseTranslationYSensitivity;

			Vec3f transXAxis = cross(effectiveUpVector(), (getPosition() - mLookAt));
			transXAxis /= sqrt(dot(transXAxis,transXAxis));
			Vec3f transYAxis = cross((getPosition() - mLookAt), transXAxis);
			transYAxis /= sqrt(dot(transYAxis,transYAxis));

			setLookAt(getLookAt() + transXAxis*xTrack + transYAxis*yTrack);

			break;
		}
	case kActionRotate:
		{
			float dAzimuth		=   -mouseDelta[0] * kMouseRotationSensitivity;
			float dElevation	=   mouseDelta[1] * kMouseRotationSensitivity;

			setAzimuth(getAzimuth() + dAzimuth);
			setElevation(getElevation() + dElevation);

			break;
		}
	case kActionZoom:
		{
			float dDolly = -mouseDelta[1] * sensitivity() * kMouseZoomSensitivity;
			if( fixedCenter_ )
				setDolly(getDolly() + dDolly);
			else
				setLookAt(getLookAt() - dDolly*norm(getPosition() - mLookAt));
			break;
		}
	case kActionTwist:
		// Not implemented
	default:
		break;
	}
}

HasKeyable<float>::AttributeList keyableAttributes( vl::Vec2f& v, const std::string& name )
{
	HasKeyable<float>::AttributeList result;

	for( unsigned int i = 0; i < 2; ++i )
	{
		std::string curName = name;
		curName.push_back('.');
		curName.push_back( 'x' + i );

		result.push_back( HasKeyable<float>::KeyableAttributePtr(
			new SimpleKeyable<float>(curName, v[i]) ) );
	}

	return result;
}

HasKeyable<float>::AttributeList keyableAttributes( vl::Vec3f& v, const std::string& name )
{
	HasKeyable<float>::AttributeList result;

	for( unsigned int i = 0; i < 3; ++i )
	{
		std::string curName = name;
		curName.push_back('.');
		curName.push_back( 'x' + i );

		result.push_back( HasKeyable<float>::KeyableAttributePtr(
			new SimpleKeyable<float>(curName, v[i]) ) );
	}

	return result;
}

HasKeyable<float>::AttributeList keyableAttributes( vl::Mat3f& m, const std::string& name )
{
	HasKeyable<float>::AttributeList result;

	for( unsigned int i = 0; i < 3; ++i )
	{
		for( unsigned int j = 0; j < 3; ++j )
		{
			std::string curName = name;
			curName.push_back('.');
			curName.push_back( 'x' + i );
			curName.push_back( 'x' + j );

			result.push_back( HasKeyable<float>::KeyableAttributePtr(
				new SimpleKeyable<float>(curName, m[i][j]) ) );
		}
	}

	return result;
}


HasKeyable<float>::AttributeList ProjectiveFixedYCamera::keyable()
{
	HasKeyable<float>::AttributeList result;
	result.push_back( HasKeyable<float>::KeyableAttributePtr(
		new SimpleKeyable<float>("azimuth", this->mAzimuth) ) );
	result.push_back( HasKeyable<float>::KeyableAttributePtr(
		new SimpleKeyable<float>("elevation", this->mElevation) ) );
	result.push_back( HasKeyable<float>::KeyableAttributePtr(
		new SimpleKeyable<float>("dolly", this->mDolly) ) );
	result.push_back( HasKeyable<float>::KeyableAttributePtr(
		new SimpleKeyable<float>("fov", this->fov_) ) );

	{
		HasKeyable<float>::AttributeList lookAtAtt = 
			keyableAttributes( this->mLookAt, "lookAt" );
		std::copy( lookAtAtt.begin(), lookAtAtt.end(), std::back_inserter(result) );
	}

	return result;
}

void ProjectiveFixedYCamera::releaseMouse( int x, int y )
{
	mCurrentMouseAction = kActionNone;
}

void ProjectiveFixedYCamera::applyViewingTransform() const
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	/*
	vl::Vec3f mPosition = getPosition();
	vl::Vec3f mUpVector = upVector();
	assert( !isBad(mPosition) );
	assert( !isBad(mLookAt) );
	assert( !isBad(mUpVector) );

	// Place the camera at mPosition, aim the camera at
	// mLookAt, and twist the camera such that mUpVector is up
	gluLookAt(	mPosition[0], mPosition[1], mPosition[2],
				mLookAt[0],   mLookAt[1],   mLookAt[2],
				mUpVector[0], mUpVector[1], mUpVector[2]);
				*/
	vl::Mat4f modelView = trans(this->modelViewTransform());
	glMultMatrixf( modelView.Ref() );
}

vl::Mat4f ProjectiveFixedYCamera::modelViewTransform() const
{
	vl::Mat4f scaleMatrix( vl::vl_1 );
	scaleMatrix *= scale();

	vl::Vec3f mPosition = getPosition();
	vl::Vec3f mUpVector = effectiveUpVector();
	assert( !isBad(mPosition) );
	assert( !isBad(mLookAt) );
	assert( !isBad(mUpVector) );

	vl::Vec3f f = vl::norm(mLookAt - mPosition);
	vl::Vec3f up = vl::norm( mUpVector );
	vl::Vec3f s = vl::norm( cross( f, up ) );
	vl::Vec3f u = cross( s, f );

	vl::Mat4f M( 
		s[0], s[1], s[2], 0.0, 
		u[0], u[1], u[2], 0.0,
		-f[0], -f[1], -f[2], 0.0,
		0.0, 0.0, 0.0, 1.0 );
    
	vl::Mat4f translation;
	translation.MakeHTrans(-mPosition);
	return scaleMatrix * M * translation;
}

vl::Mat4f gldPerspective(GLdouble fovx, GLdouble aspect, GLdouble zNear, GLdouble zFar)
{
	// borrowed from http://glprogramming.com/dgs.php?dg=1 
	GLdouble xmin, xmax, ymin, ymax;

	vl::Mat4f M( vl::vl_0 );

	xmax = zNear * tan(fovx * M_PI / 360.0);
	xmin = -xmax;

	ymin = xmin / aspect;
	ymax = xmax / aspect;

	// Set up the projection matrix
	M[0][0] = (2.0 * zNear) / (xmax - xmin);
	M[1][1] = (2.0 * zNear) / (ymax - ymin);
	M[2][2] = -(zFar + zNear) / (zFar - zNear);

	M[0][2] = (xmax + xmin) / (xmax - xmin);
	M[1][2] = (ymax + ymin) / (ymax - ymin);
	M[3][2] = -1.0;

	M[2][3] = -(2.0 * zFar * zNear) / (zFar - zNear);

	return M;
}

vl::Mat4f ProjectiveFixedYCamera::projectiveTransform() const
{
	return gldPerspective(FOV(), ratio(), this->zMin(), this->zMax());
}

void ProjectiveFixedYCamera::applyProjectiveTransform(TRcontext* tr) const
{
	if( tr )
	{
		trPerspective( tr, FOV(), ratio(), this->zMin(), this->zMax() );
	}
	else
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		vl::Mat4f projectionMatrix = trans(this->projectiveTransform());
		glMultMatrixf( projectionMatrix.Ref() );

//		gluPerspective(FOV(), ratio(),
//				this->zMin(), this->zMax());
	}
}

OrthoCamera::OrthoCamera()
	:	center_(vl::vl_0),
		port_(5.0),
		azimRot_(0.0),
		mCurrentMouseAction(kActionNone)
{
}

void OrthoCamera::setPort( float port )
{
	port_ = port;
}

void OrthoCamera::setCenter( vl::Vec2f center )
{
	center_ = center;
}

void OrthoCamera::setAzimRot( float azimRot )
{
	this->azimRot_ = azimRot;
}


void OrthoCamera::clickMouse( MouseAction_t action, int x, int y )
{
	mCurrentMouseAction = action;
	mLastMousePosition[0] = x;
	mLastMousePosition[1] = y;
}

void OrthoCamera::dragMouse( int x, int y )
{
	Vec3f mouseDelta   = Vec3f(x,y,0.0f) - mLastMousePosition;
	mLastMousePosition = Vec3f(x,y,0.0f);

	switch(mCurrentMouseAction)
	{
	case kActionTranslate:
		{
			float xTrack =  -mouseDelta[0] * port_ * sensitivity() * kMouseTranslationXSensitivity;
			float yTrack =  mouseDelta[1] * port_ * sensitivity() * kMouseTranslationYSensitivity;

			center_[0] += xTrack;
			center_[1] += yTrack;
			break;
		}
	case kActionRotate:
		{
			float dAzimuth		=   -mouseDelta[0] * kMouseRotationSensitivity;
			float dElevation	=   mouseDelta[1] * kMouseRotationSensitivity;

			azimRot_ += dAzimuth;
			break;
		}
	case kActionZoom:
		{
			float dDolly = -mouseDelta[1] * sensitivity() * kMouseZoomSensitivity;
			double temp = log(port_);
			temp -= dDolly;
			port_ = exp(temp);
			//port_ -= dDolly;
			break;
		}
	case kActionTwist:
		// Not implemented
	default:
		break;
	}
}

HasKeyable<float>::AttributeList OrthoCamera::keyable()
{
	HasKeyable<float>::AttributeList result;
	result.push_back( HasKeyable<float>::KeyableAttributePtr(
		new SimpleKeyable<float>("azimuth", this->azimRot_) ) );
	result.push_back( HasKeyable<float>::KeyableAttributePtr(
		new SimpleKeyable<float>("port", this->port_) ) );

	{
		HasKeyable<float>::AttributeList centerAtt = 
			keyableAttributes( this->center_, "center" );
		std::copy( centerAtt.begin(), centerAtt.end(), std::back_inserter(result) );
	}

	return result;
}

void OrthoCamera::releaseMouse( int x, int y )
{
	mCurrentMouseAction = kActionNone;
}

bool OrthoCamera::moving() const
{
	return mCurrentMouseAction != kActionNone;
}

void OrthoCamera::applyProjectiveTransform(TRcontext* tr) const
{
	if( tr )
	{
		trOrtho( tr, -ratio()*port_, ratio()*port_, -port_, port_, this->zMin(), this->zMax() );
	}
	else
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

//		glOrtho(-ratio()*port_, ratio()*port_, -port_, port_, this->zMin(), this->zMax());

		vl::Mat4f orthoMatrix = trans(this->projectiveTransform());
		glMultMatrixf( orthoMatrix.Ref() );
	}
}

vl::Mat4f orthoMatrix(double left, double right, double bottom, double top, double near_, double far_)
{
	double t_x = -(right + left) / (right - left);
	double t_y = -(top + bottom) / (top - bottom);
	double t_z = -(far_ + near_) / (far_ - near_);

	vl::Mat4f result( 
		2.0/(right - left), 0.0, 0.0, t_x,
		0.0, 2.0/(top - bottom), 0.0, t_y,
		0.0, 0.0, -2.0/(far_ - near_), t_z,
		0.0, 0.0, 0.0, 1.0 );

	return result;
}

vl::Mat4f OrthoCamera::projectiveTransform() const
{
	return orthoMatrix(-ratio()*port_, ratio()*port_, -port_, port_, this->zMin(), this->zMax());
}

void OrthoCamera::applyViewingTransform() const
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	/*
	glTranslatef( -center_[0], -center_[1], -20.0 );

	vl::Mat3f basis = this->basis();

	{
		// Quantize azimRot:
		const float quantum = M_PI_2;
		float quantAzim = floor( azimRot_ / quantum ) * quantum;

		vl::Mat3f rotXForm;
		rotXForm.MakeRot(basis[1], quantAzim);
		vl::Mat3f rot = basis*rotXForm;

		Vec3f uAxis = rot[0];
		Vec3f vAxis = rot[1];
		Vec3f wAxis = rot[2];

		GLfloat rotMat[16];
		rotMat[0] = uAxis[0];
		rotMat[1] = vAxis[0];
		rotMat[2] = wAxis[0];
		rotMat[3] = 0.0;
		rotMat[4] = uAxis[1];
		rotMat[5] = vAxis[1];
		rotMat[6] = wAxis[1];
		rotMat[7] = 0.0;
		rotMat[8] = uAxis[2];
		rotMat[9] = vAxis[2];
		rotMat[10] = wAxis[2];
		rotMat[11] = 0.0;
		rotMat[12] = 0.0;
		rotMat[13] = 0.0;
		rotMat[14] = 0.0;
		rotMat[15] = 1.0;

		glMultMatrixf( rotMat );
	}
	*/

	vl::Mat4f transform = trans( this->modelViewTransform() );
	glMultMatrixf( transform.Ref() );
}

vl::Mat4f OrthoCamera::modelViewTransform() const
{
	vl::Mat4f translation;
	translation.MakeHTrans( vl::Vec3f(-center_[0], -center_[1], -20.0) );

	vl::Mat3f basis = this->basis();

	// Quantize azimRot:
	const float quantum = M_PI_2;
	float quantAzim = floor( azimRot_ / quantum ) * quantum;

	vl::Mat3f rotXForm;
	rotXForm.MakeRot(basis[1], quantAzim);
	vl::Mat3f rot = basis*rotXForm;

	vl::Mat4f rotMat( vl::vl_1 );
	for( vl::Int i = 0; i < 3; ++i )
		rotMat[i] = vl::Vec4f( rot[i], 0.0 );

	return translation*rotMat;
}


void OrthoCamera::dumpParametersBinary( std::ostream& os ) const
{
	vl::Mat3f basis = this->basis();

	os.write( reinterpret_cast<const char*>(&azimRot_), sizeof(float) );
	os.write( reinterpret_cast<const char*>(&port_), sizeof(float) );
	os.write( reinterpret_cast<const char*>(&center_), sizeof(vl::Vec2f) );
	os.write( reinterpret_cast<const char*>(&basis), sizeof(vl::Mat3f) );
}

void OrthoCamera::readParametersBinary( std::istream& is )
{
	vl::Mat3f basis;
	is.read( reinterpret_cast<char*>(&azimRot_), sizeof(float) );
	is.read( reinterpret_cast<char*>(&port_), sizeof(float) );
	is.read( reinterpret_cast<char*>(&center_), sizeof(vl::Vec2f) );
	is.read( reinterpret_cast<char*>(&basis), sizeof(vl::Mat3f) );
	this->setBasis( basis );

}

void OrthoCamera::dumpParameters( std::ostream& os ) const
{
	vl::Mat3f basis = this->basis();
	os << azimRot_ << " " << port_ << std::endl;
	os << center_ << std::endl;
	os << basis;
}

void OrthoCamera::readParameters( std::istream& is )
{
	vl::Mat3f basis;
	is >> azimRot_ >> port_;
	is >> center_;
	is >> basis;
	this->setBasis( basis );
}

vl::Vec3f OrthoCamera::getPosition() const
{
	vl::Mat3f basis = this->basis();

	const float quantum = M_PI_2;
	float quantAzim = floor( azimRot_ / quantum ) * quantum;

	vl::Vec3f translation(-center_[0], -center_[1], -20.0);

	vl::Mat3f rotXForm;
	rotXForm.MakeRot(basis[1], quantAzim);
	vl::Mat3f rot = basis*rotXForm;

	return rot*translation;
}

void OrthoCamera::frame( const BoundingBox3f& box )
{
	if( box.empty() )
		return;

	// unimplemented
	assert( false );
}


