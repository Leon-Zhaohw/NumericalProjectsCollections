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
// PLANE.h: interface for the PLANE class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PLANE_H
#define PLANE_H

#include "SURFACE.h"
#include <MATRIX3.h>
#include <QUATERNION.h>
#include <TENSOR3.h>

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

class PLANE : public SURFACE  
{
public:
	PLANE(VEC3F& point, VEC3F& normal);
	virtual ~PLANE();

  VEC3F& point() { return _point; };
  VEC3F& normal() { return _normal; };
  MATRIX3& rotation() { return _rotation; };

  virtual bool inside(float* point);
  virtual bool inside(VEC3F& point);
  virtual float potential();
  virtual float distance(float* point);

  virtual void draw();

  virtual VEC3F force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity);
  virtual MATRIX springJacobian(const VEC3F& collisionPoint);
  virtual MATRIX dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity);

  MATRIX defoSpringJacobian(const VEC3F& collisionPoint, const MATRIX& U, const MATRIX3& rotation, const VEC3F& translation);
  MATRIX defoDampingJacobian(const VEC3F& collisionPoint, const MATRIX& U, const MATRIX3& rotation);

  MATRIX angularSpringJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity, const VEC3F& localPoint, const MATRIX3& rotation, const TENSOR3& rotationPartial);
  MATRIX translationDefoJacobian(const MATRIX& U, const MATRIX3& rotation);
  MATRIX translationDefoDampingJacobian(const MATRIX& U, const MATRIX3& rotation);
  MATRIX translationAngularJacobian(const VEC3F& localPoint, const TENSOR3& partialRotation);
  MATRIX translationAngularDampingJacobian(const VEC3F& localVelocity, const TENSOR3& partialRotation);
  MATRIX angularTranslationJacobian(const VEC3F& localPoint, const MATRIX3& rotation);
  MATRIX defoAngularJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity, const VEC3F& localPoint, const MATRIX3& rotation, const TENSOR3& rotationPartial, const MATRIX& U);
  MATRIX angularDefoJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity, const VEC3F& localPoint, const MATRIX& U, const MATRIX3& rotation);
  MATRIX angularDefoJacobianDamping(const VEC3F& collisionPoint, const VEC3F& collisionVelocity, const VEC3F& localPoint, const MATRIX& U, const MATRIX3& rotation);

  virtual void springJacobian(const VEC3F& collisionPoint, 
                              const VEC3F& collisionVelocity,
                              const VEC3F& collisionForce,
                              const VEC3F& localPoint,
                              const MATRIX& U, 
                              BLOCK_MATRIX& systemMatrix);

  /*
  virtual void springJacobianDebug(const VEC3F& collisionPoint, 
                              const VEC3F& collisionVelocity,
                              const VEC3F& localPoint,
                              const MATRIX& U, 
                              const MATRIX3& rotation,
                              const TENSOR3& rotationPartial,
                              const VEC3F& translation,
                              const Real* accelerationAlpha,
                              BLOCK_MATRIX& systemMatrix,
                              map<string,double>& timingBreakdown);
                              */
  virtual void springJacobianDebug(const VEC3F& collisionPoint, 
                              const VEC3F& collisionVelocity,
                              const VEC3F& collisionForce,
                              const VEC3F& localPoint,
                              const MATRIX& U, 
                              BLOCK_MATRIX& systemMatrix,
                              map<string,double>& timingBreakdown);

  // debug function -- verify that the spring Jacobian is correct
  virtual void verifySpringJacobian(VEC3F& collisionPoint);

private:
  VEC3F _point;
  VEC3F _normal;

  // when drawing the plane, apply this rotation about _point
  MATRIX3 _rotation;

  VEC3F _axis;
  Real _angle;

  MATRIX _springJacobian;
  MATRIX _dampingJacobian;
};

#endif
