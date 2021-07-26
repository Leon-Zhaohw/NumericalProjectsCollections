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
// BOX.h: interface for the BOX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef BOX_H
#define BOX_H

#include "SURFACE.h"
#include "TRIANGLE.h"
#include "PLANE.h"

class BOX : public SURFACE  
{
public:
	BOX(float xPlus = 0.55f, float xMinus = 0.45f, 
      float yPlus = 0.55f, float yMinus = 0.45f, 
      float zPlus = 0.55f, float zMinus = 0.45f);
  BOX(const VEC3F& center, const VEC3F& sizes);
	virtual ~BOX();

  virtual bool inside(float* point);
  virtual float potential();
  virtual float distance(float* point);

  virtual void draw();

  void translateX(float delta = 0.01);
  void translateY(float delta = 0.01);
  void translateZ(float delta = 0.01);

  void scaleX(float delta = 0.01);
  void scaleY(float delta = 0.01);
  void scaleZ(float delta = 0.01);

  VEC3F boxMins() { return VEC3F(_xMinus, _yMinus, _zMinus); };
  VEC3F boxMaxs() { return VEC3F(_xPlus, _yPlus, _zPlus); };

  // explicitly load the overlap test here, since the code is ugly and we
  // want to encapsulate it
  bool overlap(TRIANGLE& triangle);

  virtual VEC3F force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity);
  virtual MATRIX springJacobian(const VEC3F& collisionPoint);
  virtual MATRIX dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity);

  virtual void springJacobian(const VEC3F& collisionPoint, 
                              const VEC3F& collisionVelocity,
                              const VEC3F& collisionForce,
                              const VEC3F& localPoint,
                              const MATRIX& U, 
                              BLOCK_MATRIX& systemMatrix);

private:
  float _xPlus, _xMinus;
  float _yPlus, _yMinus;
  float _zPlus, _zMinus;

  inline float mymin(float i, float j) { return (i < j) ? i : j; };

  // return the plane corresponding to the closest face, assuming
  // the point is inside the box
  PLANE closestFace(const VEC3F& collisionPoint);
};

#endif
