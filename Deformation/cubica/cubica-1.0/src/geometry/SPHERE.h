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
// SPHERE.h: interface for the SPHERE class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SPHERE_H
#define SPHERE_H

#include "SURFACE.h"

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

class SPHERE : public SURFACE
{
public:
	SPHERE(float radius = 0.5f, float potential = 0.0f);
	SPHERE(float radius, float* center);
	SPHERE(float radius, VEC3F& center);
  SPHERE(const SPHERE& sphere);
	virtual ~SPHERE();

  virtual bool inside(float* point);
  virtual float potential() { return _potential; };
  virtual float distance(float* point);
  VEC3F& center() { return _center; };
  Real& radius() { return _radius; };
  bool intersect(SPHERE& rightSphere);

  virtual VEC3F force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity);
  virtual MATRIX springJacobian(const VEC3F& collisionPoint);
  virtual MATRIX dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity);

  void draw();

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
  
  // return the bounding box dimensions
  virtual void boundingBox(VEC3F& mins, VEC3F& maxs);

private:
  Real _radius;
  float _potential;
  VEC3F _center;
};

#endif
