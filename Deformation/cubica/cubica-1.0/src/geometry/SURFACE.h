/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
// SURFACE.h: interface for the SURFACE class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SURFACE_H
#define SURFACE_H

#include <cmath>
#include <VEC3.h>
#include <MATRIX.h>
#include <TENSOR3.h>
#include <TIMER.h>

#ifndef M_PI
#define M_PI 3.1415926535897931f
#endif

class SURFACE  
{
public:
	SURFACE();
	virtual ~SURFACE();
  
  virtual bool inside(float* point) {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " inside UNIMPLEMENTED" << endl;
    return false;
  };
  virtual float distance(float* point) = 0;

  virtual bool inside(VEC3F& point) {
    float pointFloat[3];
    pointFloat[0] = point[0];
    pointFloat[1] = point[1];
    pointFloat[2] = point[2];
    return inside(pointFloat);
  };
  virtual float distance(VEC3F& point) {
    float pointFloat[3];
    pointFloat[0] = point[0];
    pointFloat[1] = point[1];
    pointFloat[2] = point[2];
    return distance(pointFloat);
  };

  virtual bool isOccupied(Real* point) {
    float fPoint[] = {point[0], point[1], point[2]};
    return inside(fPoint);
  };

  void setRotateY(float rotateY) {_rotateY = rotateY / 360.0f * 2.0f * 3.14159f;};
  void setRotateX(float rotateX) {_rotateX = rotateX / 360.0f * 2.0f * 3.14159f;};
  void setRotateZ(float rotateZ) {_rotateZ = rotateZ / 360.0f * 2.0f * 3.14159f;};

  // collision response functions
  virtual VEC3F force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity); 
  virtual MATRIX springJacobian(const VEC3F& collisionPoint);
  virtual MATRIX dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity);

  // multibody version -- adds the necessary quantities to 3x3 block matrix 'blockJacobian'
  virtual void springJacobian(const VEC3F& collisionPoint, 
                              const VEC3F& collisionVelocity,
                              const VEC3F& collisionForce,
                              const VEC3F& localPoint,
                              const MATRIX& U, 
                              BLOCK_MATRIX& blockJacobian);

  virtual void springJacobianDebug(const VEC3F& collisionPoint, 
                              const VEC3F& collisionVelocity,
                              const VEC3F& collisionForce,
                              const VEC3F& localPoint,
                              const MATRIX& U, 
                              BLOCK_MATRIX& systemMatrix,
                              map<string,double>& timingBreakdown);

  // debug function -- verify that the spring Jacobian is correct
  virtual void verifySpringJacobian(VEC3F& collisionPoint); 

  // draw the surface to OpenGL
  virtual void draw();

  // retrieve type info
  const string& type() { return _type; };

  Real& collisionStiffness() { return _collisionStiffness; };
  Real& collisionDamping() { return _collisionDamping; };
  Real& stickiness() { return _stickiness; };

  // set collision variables
  Real*& collisionAlpha()             { return _collisionAlpha; };
  MATRIX3& collisionRotation()        { return _collisionRotation; };
  TENSOR3& collisionRotationPartial() { return _collisionRotationPartial; };
  VEC3F& collisionTranslation()       { return _collisionTranslation; };

  // return the bounding box dimensions
  virtual void boundingBox(VEC3F& mins, VEC3F& maxs);

protected:
  float _rotateX;
  float _rotateY;
  float _rotateZ;

  void rotationX(float* rotated);
  void rotationY(float* rotated);
  void rotationZ(float* rotated);
  void rotate(float* rotated);
  void rotate(VEC3F& rotated);

  // collision response constants
  Real _collisionStiffness;
  Real _collisionDamping;

  // string holding the type name -- not pretty, but (hopefully) better 
  // than using RTTI
  string _type;

  // cached collision variables -- passing them in every time appears to incur
  // significant overhead
  Real* _collisionAlpha;
  MATRIX3 _collisionRotation;
  TENSOR3 _collisionRotationPartial;
  VEC3F   _collisionTranslation;

  Real _stickiness;
};

#endif
