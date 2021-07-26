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
// BOX.cpp: implementation of the BOX class.
//
//////////////////////////////////////////////////////////////////////

#include "BOX.h"

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BOX::BOX(float xPlus, float xMinus, 
         float yPlus, float yMinus, 
         float zPlus, float zMinus) :
  _xPlus(xPlus), _xMinus(xMinus), _yPlus(yPlus),
  _yMinus(yMinus), _zPlus(zPlus), _zMinus(zMinus)
{
  _type.assign("BOX");
}


BOX::BOX(const VEC3F& center, const VEC3F& sizes)
{
  _type.assign("BOX");
  _xPlus  = center[0] + sizes[0] * 0.5;
  _xMinus = center[0] - sizes[0] * 0.5;
  _yPlus  = center[1] + sizes[1] * 0.5;
  _yMinus = center[1] - sizes[1] * 0.5;
  _zPlus  = center[2] + sizes[2] * 0.5;
  _zMinus = center[2] - sizes[2] * 0.5;
}

BOX::~BOX()
{
}

bool BOX::inside(float* point) {
  float rotated[] = {point[0], point[1], point[2]};
  rotate(rotated);

  if (rotated[0] <= _xPlus && rotated[0] >= _xMinus &&
      rotated[1] <= _yPlus && rotated[1] >= _yMinus &&
      rotated[2] <= _zPlus && rotated[2] >= _zMinus)
    return true;

  return false;
}

float BOX::potential() { 
  return 0.0f; 
}

float BOX::distance(float* point) {
  float rotated[] = {point[0], point[1], point[2]};
  rotate(rotated);

  if (inside(rotated))
  {
    float xMin = mymin((_xPlus - rotated[0]), (rotated[0] - _xMinus));
    float yMin = mymin((_yPlus - rotated[1]), (rotated[1] - _yMinus));
    float zMin = mymin((_zPlus - rotated[2]), (rotated[2] - _zMinus));

    return fabs(mymin(mymin(xMin, yMin), zMin));
  }
 
  float diff[] = {0.0f, 0.0f, 0.0f};

  // handle edges and points
  if (rotated[2] > _zPlus) 
    diff[2] = rotated[2] - _zPlus;
  else if (rotated[2] < _zMinus)
    diff[2] = _zMinus - rotated[2];
  
  if (rotated[1] > _yPlus)
    diff[1] = rotated[1] - _yPlus;
  else if (rotated[1] < _yMinus)
  diff[1] = _yMinus - rotated[1];
  
  if (rotated[0] > _xPlus)
    diff[0] = rotated[0] - _xPlus;
  else if (rotated[0] < _xMinus)
    diff[0] = _xMinus - rotated[0];
      
  return sqrtf(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
}

void BOX::draw()
{
  VEC3F v000(_xMinus, _yMinus, _zMinus); 
  VEC3F v100(_xPlus, _yMinus, _zMinus); 
  VEC3F v010(_xMinus, _yPlus, _zMinus); 
  VEC3F v110(_xPlus, _yPlus, _zMinus); 
  VEC3F v001(_xMinus, _yMinus, _zPlus); 
  VEC3F v101(_xPlus, _yMinus, _zPlus); 
  VEC3F v011(_xMinus, _yPlus, _zPlus); 
  VEC3F v111(_xPlus, _yPlus, _zPlus); 

  rotate(v000);
  rotate(v100);
  rotate(v010);
  rotate(v110);
  rotate(v001);
  rotate(v101);
  rotate(v011);
  rotate(v111);

  glBegin(GL_QUADS);
    // x plus
    VEC3F normal = cross(v000 - v100, v000 - v110);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal);
    glVertex3dv(v010);
    glVertex3dv(v110);
    glVertex3dv(v100);
    glVertex3dv(v000);

    // x minus
    normal = cross(v001 - v101, v001 - v111);
    normal.normalize();
    glNormal3dv(normal);
    glVertex3dv(v001);
    glVertex3dv(v101);
    glVertex3dv(v111);
    glVertex3dv(v011);

    // y minus
    normal = cross(v000 - v100, v000 - v101);
    normal.normalize();
    glNormal3dv(normal);
    glVertex3dv(v000);
    glVertex3dv(v100);
    glVertex3dv(v101);
    glVertex3dv(v001);

    // y plus
    normal = cross(v010 - v110, v010 - v111);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal);
    glVertex3dv(v011);
    glVertex3dv(v111);
    glVertex3dv(v110);
    glVertex3dv(v010);

    // z plus
    normal = cross(v000 - v010, v000 - v011);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal);
    glVertex3dv(v001);
    glVertex3dv(v011);
    glVertex3dv(v010);
    glVertex3dv(v000);

    // z minus
    normal = cross(v100 - v110, v100 - v111);
    normal.normalize();
    glNormal3dv(normal);
    glVertex3dv(v100);
    glVertex3dv(v110);
    glVertex3dv(v111);
    glVertex3dv(v101);
  glEnd();
}

void BOX::translateX(float delta)
{
  _xPlus += delta;
  _xMinus += delta;
}

void BOX::translateY(float delta)
{
  _yPlus += delta;
  _yMinus += delta;
}

void BOX::translateZ(float delta)
{
  _zPlus += delta;
  _zMinus += delta;
}

void BOX::scaleX(float delta)
{
  _xPlus += delta;
  _xMinus -= delta;
}

void BOX::scaleY(float delta)
{
  _yPlus += delta;
  _yMinus -= delta;
}

void BOX::scaleZ(float delta)
{
  _zPlus += delta;
  _zMinus -= delta;
}

#include "BOX_TRIANGLE_INTERSECTION.cpp"

//////////////////////////////////////////////////////////////////////
// explicitly load the overlap test here, since the code is ugly and we
// want to encapsulate it
//////////////////////////////////////////////////////////////////////
bool BOX::overlap(TRIANGLE& triangle)
{
  float boxCenter[3];
  float boxHalfSize[3];
  float triVertex[3][3];

  boxCenter[0] = (_xMinus + _xPlus) * 0.5;
  boxCenter[1] = (_yMinus + _yPlus) * 0.5;
  boxCenter[2] = (_zMinus + _zPlus) * 0.5;

  boxHalfSize[0] = (_xPlus - _xMinus) * 0.5;
  boxHalfSize[1] = (_yPlus - _yMinus) * 0.5;
  boxHalfSize[2] = (_zPlus - _zMinus) * 0.5;

  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
      triVertex[x][y] = (*(triangle.vertex(x)))[y];

  return triBoxOverlap(boxCenter, boxHalfSize, triVertex);
}

//////////////////////////////////////////////////////////////////////
// return the plane corresponding to the closest face, assuming
// the point is inside the box
//////////////////////////////////////////////////////////////////////
PLANE BOX::closestFace(const VEC3F& collisionPoint)
{
  float diffs[6];

  diffs[0] = collisionPoint[0] - _xMinus;
  diffs[1] = _xPlus - collisionPoint[0];

  diffs[2] = collisionPoint[1] - _yMinus;
  diffs[3] = _yPlus - collisionPoint[1];

  diffs[4] = collisionPoint[2] - _zMinus;
  diffs[5] = _zPlus - collisionPoint[2];

  Real minFound = diffs[0];
  int minIndex = 0;

  for (int x = 1; x < 6; x++)
  {
    if (diffs[x] < minFound)
    {
      minFound = diffs[x];
      minIndex = x;
    }
  }

  VEC3F point(_xMinus, 0, 0);
  VEC3F normal(-1, 0, 0);

  switch (minIndex)
  {
    case 1:
      point = VEC3F(_xPlus, 0,0);
      normal = VEC3F(1,0, 0);
      break;
    case 2:
      point = VEC3F(0, _yMinus,0);
      normal = VEC3F(0, -1,0);
      break;
    case 3:
      point = VEC3F(0, _yPlus,0);
      normal = VEC3F(0, 1,0);
      break;
    case 4:
      point = VEC3F(0, 0, _zMinus);
      normal = VEC3F(0, 0, -1);
      break;
    case 5:
      point = VEC3F(0, 0, _zPlus);
      normal = VEC3F(0, 0, 1);
      break;
  }

  PLANE final(point, normal);
  final.collisionStiffness() = _collisionStiffness;
  final.collisionDamping() = _collisionDamping;

  return final;
}

VEC3F BOX::force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity)
{
  PLANE plane = closestFace(collisionPoint);
  return plane.force(collisionPoint, collisionVelocity);
}

MATRIX BOX::springJacobian(const VEC3F& collisionPoint)
{
  PLANE plane = closestFace(collisionPoint);
  return plane.springJacobian(collisionPoint);
}

MATRIX BOX::dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity)
{
  PLANE plane = closestFace(collisionPoint);
  return plane.dampingJacobian(collisionPoint, collisionVelocity);
}

//////////////////////////////////////////////////////////////////////
// multibody version -- adds the necessary quantities to 3x3 block 
// matrix 'blockJacobian'
//////////////////////////////////////////////////////////////////////
void BOX::springJacobian(const VEC3F& collisionPoint, 
                         const VEC3F& collisionVelocity,
                         const VEC3F& collisionForce,
                         const VEC3F& localPoint,
                         const MATRIX& U, 
                         BLOCK_MATRIX& systemMatrix)
{
  PLANE plane = closestFace(collisionPoint);
  plane.collisionAlpha() = _collisionAlpha;
  plane.collisionTranslation() = _collisionTranslation;
  plane.collisionRotation() = _collisionRotation;
  plane.collisionRotationPartial() = _collisionRotationPartial;
  return plane.springJacobian(collisionPoint, collisionVelocity, collisionForce, localPoint, U, systemMatrix);
}
