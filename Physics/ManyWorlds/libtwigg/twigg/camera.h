// This camera stuff mostly by ehsu

#ifndef CAMERA_H
#define CAMERA_H

#include <iosfwd>
#include <algorithm>
#include "twigg/keyable.h"
#include "twigg/tr.h"
#include "twigg/boundingbox.h"

//==========[ class Camera ]===================================================

typedef enum { kActionNone, kActionTranslate, kActionRotate, kActionZoom, kActionTwist,} MouseAction_t;

class Camera
	: public HasKeyable<float>
{
public:
	Camera();
	virtual ~Camera();

	//---[ Interactive Adjustment ]------------------------
	virtual void clickMouse( MouseAction_t action, int x, int y ) = 0;
	virtual void dragMouse( int x, int y ) = 0;
	virtual void releaseMouse( int x, int y ) = 0;
	virtual bool moving() const = 0;

	//---[ Viewing Transform ]--------------------------------
	virtual vl::Mat4f projectiveTransform() const = 0;
	virtual vl::Mat4f modelViewTransform() const = 0;

	virtual void applyProjectiveTransform(TRcontext* tr) const = 0;
	virtual void applyViewingTransform() const = 0;

	virtual void dumpParameters( std::ostream& os ) const = 0;
	virtual void readParameters( std::istream& is ) = 0;

	virtual void dumpParametersBinary( std::ostream& os ) const = 0;
	virtual void readParametersBinary( std::istream& os ) = 0;

	virtual void setLookAt( const vl::Vec3f &lookAt ) = 0;
	virtual vl::Vec3f getLookAt() const = 0;
	virtual vl::Vec3f getPosition() const = 0;

	// frame the given object in the camera window
	virtual void frame( const BoundingBox3f& box ) = 0;

	void setZRange( double zMin, double zMax );
	double zMin() const;
	double zMax() const;

	vl::Vec3f upVector() const;
	void setUpVector(const vl::Vec3f& up);

	void setRatio(double ratio)		{ ratio_ = ratio; }
	double ratio() const;

	void setTranslationalSensitivity(float value);
	float sensitivity() const;

	double scale() const;
	void setScale( double scale );

protected:
	vl::Mat3f basis() const;
	void setBasis( const vl::Mat3f& basis );

private:
	double zMin_;
	double zMax_;
	double ratio_;
	float sensitivity_;
	double scale_;

	vl::Mat3f basis_;
};

class ProjectiveFixedYCamera
	: public Camera
{
private:
	vl::Mat3f	basis;
	float		fov_;

	float		mElevation;
    float		mAzimuth;
    float		mDolly;
    float		mTwist; // Not implemented yet

    vl::Vec3f	mLookAt;

    vl::Vec3f			mLastMousePosition;
    MouseAction_t		mCurrentMouseAction;

	bool		fixedCenter_;

public:

    //---[ Constructors ]----------------------------------
    // defaults to (0,0,0) facing down negative z axis
    ProjectiveFixedYCamera(bool fixedCenter = true);

    void setElevation( float elevation )
    {
        // don't want elevation to be negative
        if (elevation<0) elevation+=6.28318530717f;
        mElevation = elevation;
    }
    float getElevation() const
    { return mElevation; }

    void setAzimuth( float azimuth )
    { mAzimuth = azimuth; }
    float getAzimuth() const
    { return mAzimuth; }

    void setDolly( float dolly )
    { mDolly = dolly; }
    float getDolly() const
    { return mDolly; }

    void setTwist( float twist )
    { mTwist = twist; }
    float getTwist() const
    { return mTwist; }

    void setLookAt( const vl::Vec3f &lookAt )
    { mLookAt = lookAt; }
    vl::Vec3f getLookAt() const
    { return mLookAt; }

	vl::Vec3f getPosition() const;

	void frame( const BoundingBox3f& box );

	void setFOV( double fov )	{ fov_ = fov;	}
	double FOV() const			{ return fov_; }

    //---[ Interactive Adjustment ]------------------------
    // these should be used from a mouse event handling routine that calls
    // the startX method on a mouse down, updateX on mouse move and finally
    // endX on mouse up.
    //-----------------------------------------------------
    void clickMouse( MouseAction_t action, int x, int y );
    void dragMouse( int x, int y );
    void releaseMouse( int x, int y );

    //---[ Viewing Transform ]--------------------------------
	void applyProjectiveTransform(TRcontext* tr = 0) const;
    void applyViewingTransform() const;

	vl::Mat4f projectiveTransform() const;
	vl::Mat4f modelViewTransform() const;

	void dumpParameters( std::ostream& os ) const;
	void readParameters( std::istream& is );
	void dumpParametersBinary( std::ostream& os ) const;
	void readParametersBinary( std::istream& os );

	bool moving() const	{ return mCurrentMouseAction != kActionNone; }

	HasKeyable<float>::AttributeList keyable();

private:
	vl::Vec3f effectiveUpVector() const;
};

class OrthoCamera
	: public Camera
{
public:
	OrthoCamera();

	//---[ Interactive Adjustment ]------------------------
	void clickMouse( MouseAction_t action, int x, int y );
	void dragMouse( int x, int y );
	void releaseMouse( int x, int y );
	bool moving() const;

	//---[ Viewing Transform ]--------------------------------
	void applyProjectiveTransform(TRcontext* tr) const;
	void applyViewingTransform() const;

	vl::Mat4f projectiveTransform() const;
	vl::Mat4f modelViewTransform() const;

	void dumpParameters( std::ostream& os ) const;
	void readParameters( std::istream& is );
	void dumpParametersBinary( std::ostream& os ) const;
	void readParametersBinary( std::istream& is );

	void frame( const BoundingBox3f& box );
	void setLookAt( const vl::Vec3f &lookAt ) {}
	vl::Vec3f getLookAt() const			{ return vl::Vec3f(vl::vl_0); }
	vl::Vec3f getPosition() const;

	HasKeyable<float>::AttributeList keyable();

	void setPort( float port );
	void setCenter( vl::Vec2f center );
	void setAzimRot( float azimRot );

	vl::Vec2f getCenter() const            { return center_; }
	float getPort() const                  { return port_; }
	float getAzimRot() const               { return azimRot_; }

private:
	vl::Vec2f	center_;

	float		azimRot_;
	float		port_;

    vl::Vec3f		mLastMousePosition;
    MouseAction_t	mCurrentMouseAction;
};


#endif
