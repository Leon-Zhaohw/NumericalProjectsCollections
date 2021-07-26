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
// CYLINDER.cpp: implementation of the CYLINDER class.
//
//////////////////////////////////////////////////////////////////////

#include "CYLINDER.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CYLINDER::CYLINDER()
{
  _center[0] = 0.5f;
  _center[1] = 0.5f;

  _caps[0] = 0.2f;
  _caps[1] = 0.8f;

  _radius = 0.1f;

  _radiusSq = _radius * _radius;
  _type.assign("CYLINDER");

  _cylinderTranslation = VEC3F(0,0,0);
  _cylinderRotation = QUATERNION(0,0,0,1);
}

CYLINDER::CYLINDER(float* center, float radius, float* caps)
{
  _center[0] = center[0];
  _center[1] = center[1];

  _caps[0] = caps[0];
  _caps[1] = caps[1];

  _radius = radius;
  _radiusSq = _radius * _radius;
  _type.assign("CYLINDER");
  _cylinderTranslation = VEC3F(0,0,0);
  _cylinderRotation = QUATERNION(0,0,0,1);
}

CYLINDER::~CYLINDER()
{

}

//////////////////////////////////////////////////////////////////////
// deep copy
//////////////////////////////////////////////////////////////////////
SURFACE* CYLINDER::copy() {
  CYLINDER* final = new CYLINDER(_center, _radius, _caps);
  final->cylinderTranslation() = _cylinderTranslation;
  final->cylinderRotation() = _cylinderRotation;

  return final;
}

//////////////////////////////////////////////////////////////////////
// inside outside function
//////////////////////////////////////////////////////////////////////
bool CYLINDER::inside(float* p) {
  VEC3F point;
  point[0] = p[0];
  point[1] = p[1];
  point[2] = p[2];

  // transform into local frame
  point -= _cylinderTranslation;
  MATRIX3 RT = _cylinderRotation.toExplicitMatrix3x3().transpose();
  point = RT * point;

  if (point[1] > _caps[1] || point[1] < _caps[0]) return false;

  float diff[2];
  diff[0] = point[0] - _center[0];
  diff[1] = point[2] - _center[1];

  float magnitude = diff[0] * diff[0] + diff[1] * diff[1];

  return magnitude <= _radiusSq;
}

//////////////////////////////////////////////////////////////////////
// distance field function
//////////////////////////////////////////////////////////////////////
float CYLINDER::distance(float* point2) {
  float point[] = {point2[0], point2[1], point2[2]};
  
  // distance from the cylinder center
  float diff[2];
  diff[0] = point[0] - _center[0];
  diff[1] = point[2] - _center[1];
  float distance = sqrtf(diff[0] * diff[0] + diff[1] * diff[1]);

  // do the inside case
  if (inside(point))
  {
    // distance from caps
    float bottom = point[1] - _caps[0];
    float top    = _caps[1] - point[1];

    // distance from the cylinder wall
    distance = _radius - distance;
    
    // return the smallest
    if (bottom < distance) distance = bottom;
    return (top < distance) ? fabs(top) : fabs(distance);
  }

  // if it is outside and projects in the y direction 
  // onto the top or bottom cap
  if (distance < _radius)
  {
    if (point[1] > _caps[1])
      return fabs(point[1] - _caps[1]);

    return fabs(_caps[0] - point[1]);
  }

  // if the nearest point is along the cylinder wall
  if (point[1] < _caps[1] && point[1] > _caps[0])
    return fabs(distance - _radius);

  // else it must be closest to a cap edge
  if (point[1] > _caps[1])
    diff[1] = point[1] - _caps[1];
  else
    diff[1] = _caps[0] - point[1];

  diff[0] = distance - _radius;

  return sqrtf(diff[0] * diff[0] + diff[1] * diff[1]);
}

//////////////////////////////////////////////////////////////////////
// draw the surface to OpenGL
//////////////////////////////////////////////////////////////////////
void CYLINDER::drawCapsule()
{
  float length = fabs(_caps[0] - _caps[1]);
  float radius = _radius;

  glPushMatrix();
    VEC3F axis;
    Real angle;
    _cylinderRotation.axisAngle(axis, angle);
    glTranslatef(_cylinderTranslation[0], _cylinderTranslation[1], _cylinderTranslation[2]);
    glRotatef(angle, axis[0], axis[1], axis[2]);

    Real cylHalfHeight = length / 2.0;

    glBegin(GL_QUAD_STRIP);
    for (int i = 0; i < 17; i++)
    {
      Real angle = i / 16.0 * 2.0 * M_PI;
      Real ca = cos(angle);
      Real sa = sin(angle);

      VEC3F normal(ca, sa, 0);
      normal = _cylinderRotation.toExplicitMatrix3x3() * normal;
      glNormal3f(normal[0], normal[1], normal[2]);
 
      glVertex3f(radius * ca, cylHalfHeight, radius * sa);
      glVertex3f(radius * ca, -cylHalfHeight, radius * sa);
    }

    glEnd();

    // draw the caps
    glTranslated(0, cylHalfHeight,0);
    glutSolidSphere(radius, 16, 12);
    glTranslated(0, -2.0 * cylHalfHeight,0);
    glutSolidSphere(radius, 16, 12);

    // draw a line towards one of the end caps
    glDisable(GL_DEPTH_TEST);
    glColor4f(0,10,0,1);
    glLineWidth(10);
    glBegin(GL_LINES);
      glVertex3f(0, 0, 0);
      glVertex3f(0, cylHalfHeight, 0);
    glEnd();
    glEnable(GL_DEPTH_TEST);
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// draw the surface to OpenGL
//////////////////////////////////////////////////////////////////////
void CYLINDER::draw()
{
  float length = fabs(_caps[0] - _caps[1]);
  float radius = _radius;

  glPushMatrix();
    VEC3F axis;
    Real angle;
    _cylinderRotation.axisAngle(axis, angle);
    glTranslatef(_cylinderTranslation[0], _cylinderTranslation[1], _cylinderTranslation[2]);
    glRotatef(angle, axis[0], axis[1], axis[2]);

    Real cylHalfHeight = length / 2.0;

    glBegin(GL_QUAD_STRIP);
    for (int i = 0; i < 17; i++)
    {
      Real angle = i / 16.0 * 2.0 * M_PI;
      Real ca = cos(angle);
      Real sa = sin(angle);

      VEC3F normal(ca, sa, 0);
      normal = _cylinderRotation.toExplicitMatrix3x3() * normal;
      glNormal3f(normal[0], normal[1], normal[2]);
 
      glVertex3f(radius * ca, cylHalfHeight, radius * sa);
      glVertex3f(radius * ca, -cylHalfHeight, radius * sa);
    }
    glEnd();

    // draw a line towards one of the end caps
    glDisable(GL_DEPTH_TEST);
    glColor4f(0,10,0,1);
    glLineWidth(10);
    glBegin(GL_LINES);
      glVertex3f(0, 0, 0);
      glVertex3f(0, cylHalfHeight, 0);
    glEnd();
    glEnable(GL_DEPTH_TEST);
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// return the bounding box dimensions
//////////////////////////////////////////////////////////////////////
void CYLINDER::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  // get the eight corner of the original bounding box
  VEC3F corners[8];

  Real distance = sqrtf(_radius * _radius + _radius * _radius);

  Real length = fabs(_caps[0] - _caps[1]);

  corners[0] = VEC3F(distance, length * 0.5, distance); 
  corners[1] = VEC3F(-distance, length * 0.5, distance); 
  corners[2] = VEC3F(distance, length * 0.5, -distance); 
  corners[3] = VEC3F(-distance, length * 0.5, -distance);

  corners[4] = VEC3F(distance, -length * 0.5, distance); 
  corners[5] = VEC3F(-distance, -length * 0.5, distance); 
  corners[6] = VEC3F(distance, -length * 0.5, -distance); 
  corners[7] = VEC3F(-distance, -length * 0.5, -distance);

  // transform the corners
  MATRIX3 rotation = _cylinderRotation.toExplicitMatrix3x3();
  for (int x = 0; x < 8; x++)
    corners[x] =  rotation * corners[x] + _cylinderTranslation;

  // find the bounds
  mins = corners[0];
  maxs = corners[0];

  for (int x = 1; x < 8; x++)
    for (int y = 0; y < 3; y++)
    {
      if (corners[x][y] < mins[y])
        mins[y] = corners[x][y];
      if (corners[x][y] > maxs[y])
        maxs[y] = corners[x][y];
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void CYLINDER::setLength(float length)
{
  Real middle = (_caps[1] + _caps[0]) * 0.5;

  _caps[0] = middle - length * 0.5;
  _caps[1] = middle + length * 0.5;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void CYLINDER::endPoints(VEC3F& leftCap, VEC3F& rightCap)
{
  VEC3F& translation = _cylinderTranslation;
  MATRIX3 rotation = _cylinderRotation.toExplicitMatrix3x3();
  Real length = getLength();
  VEC3F beginVertex(0,0, length * 0.5 + _radius);
  VEC3F endVertex(0,0, -length * 0.5 - _radius);

  // apply the ODE-style rotation first
  MATRIX3 first = MATRIX3::rotation(VEC3F(1,0,0), M_PI / 2);
  rotation = rotation * first;

  leftCap  = rotation * beginVertex + translation;
  rightCap = rotation * endVertex + translation;
}
