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
// CYLINDER.h: interface for the CYLINDER class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CYLINDER_H
#define CYLINDER_H

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

#include "QUATERNION.h"

class CYLINDER : public SURFACE  
{
public:
  CYLINDER();
	CYLINDER(float* center, float radius, float* caps);
	virtual ~CYLINDER();

  virtual bool inside(float* p);
  virtual float potential() { return 0.0f; };
  virtual float distance(float* point);
  virtual SURFACE* copy();

  virtual void draw();

  // draw caps on the ends
  virtual void drawCapsule();

  QUATERNION& cylinderRotation() { return _cylinderRotation; };
  VEC3F& cylinderTranslation()   { return _cylinderTranslation; };
  const float& radius()          { return _radius; };
  void setRadius(float radius) { _radius = radius; _radiusSq = _radius * _radius; };
  float getLength()             { return _caps[1] - _caps[0]; };
  void setLength(float length);
  float* caps()   { return _caps; };
  
  // return the bounding box dimensions
  virtual void boundingBox(VEC3F& mins, VEC3F& maxs);

  // get the (ODE VERSION) of the cylinder end caps, i.e.
  // aligned along the y axis
  void endPoints(VEC3F& leftCap, VEC3F& rightCap);

private:
  // cylinder rotation
  QUATERNION _cylinderRotation;

  // cylinder translation
  VEC3F _cylinderTranslation;

  // x,z coordinates of center axis
  float _center[2];
  float _radius;
  // bottom and top y coordinate of caps
  float _caps[2];
  float _radiusSq;
};

#endif
