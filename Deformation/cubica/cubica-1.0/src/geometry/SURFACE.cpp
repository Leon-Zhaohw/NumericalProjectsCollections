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
// SURFACE.cpp: implementation of the SURFACE class.
//
//////////////////////////////////////////////////////////////////////

#include "SURFACE.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SURFACE::SURFACE() :
  _rotateX(0.0f),
  _rotateY(0.0f),
  _rotateZ(0.0f),
  _type("SURFACE")
{

}

SURFACE::~SURFACE()
{

}

void SURFACE::rotationX(float* rotated)
{
  float point[] = {rotated[0], rotated[1], rotated[2]};
  rotated[0] = point[0];
  rotated[1] = point[1] * cos(_rotateX) + point[2] * sin(_rotateX);
  rotated[2] = -point[1] * sin(_rotateX) + point[2] * cos(_rotateX);
}

void SURFACE::rotationY(float* rotated)
{
  float point[] = {rotated[0], rotated[1], rotated[2]};
  rotated[0] = point[0] * cos(_rotateY) - point[2] * sin(_rotateY);
  rotated[1] = point[1];
  rotated[2] = point[0] * sin(_rotateY) + point[2] * cos(_rotateY);
}

void SURFACE::rotationZ(float* rotated)
{
  float point[] = {rotated[0], rotated[1], rotated[2]};
  rotated[0] = point[0] * cos(_rotateZ) + point[1] * sin(_rotateZ);
  rotated[1] = -point[0] * sin(_rotateZ) + point[1] * cos(_rotateZ);
  rotated[2] = point[2];
}

void SURFACE::rotate(float* rotated)
{
  rotationX(rotated);
  rotationY(rotated);
  rotationZ(rotated);
}

void SURFACE::rotate(VEC3F& rotated)
{
  float data[3];
  data[0] = rotated[0];
  data[1] = rotated[1];
  data[2] = rotated[2];

  rotationX(data);
  rotationY(data);
  rotationZ(data);

  rotated[0] = data[0];
  rotated[1] = data[1];
  rotated[2] = data[2];
}

//////////////////////////////////////////////////////////////////////
// compute the force at the current point and velocity
//////////////////////////////////////////////////////////////////////
VEC3F SURFACE::force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " FORCE NOT IMPLEMENTED FOR THIS SURFACE " << endl;
  return VEC3F();
}

//////////////////////////////////////////////////////////////////////
// compute the force jacobian at the current point
//////////////////////////////////////////////////////////////////////
MATRIX SURFACE::springJacobian(const VEC3F& collisionPoint)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " SPRING JACOBIAN NOT IMPLEMENTED FOR THIS SURFACE " << endl;
  return MATRIX(3,3);
}

//////////////////////////////////////////////////////////////////////
// compute the damping jacobian at the current point
//////////////////////////////////////////////////////////////////////
MATRIX SURFACE::dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " DAMPING JACOBIAN NOT IMPLEMENTED FOR THIS SURFACE " << endl;
  return MATRIX(3,3);
}
  
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SURFACE::verifySpringJacobian(VEC3F& collisionPoint)
{
  VEC3F velocity;
  VEC3F originalForce = force(collisionPoint, velocity);
  VEC3F originalPoint = collisionPoint;

  cout << " original force: " << originalForce << endl;

  Real delta = 1e-6;
  vector<VECTOR> columns;
  for (int x = 0; x < 3; x++)
  {
    VEC3F perturb = originalPoint;
    perturb[x] += delta;

    VEC3F newForce = force(perturb, velocity);
    VEC3F diff = newForce - originalForce;
    diff *= 1.0 / delta;
    columns.push_back(diff.toVector());
  }

  MATRIX finiteDiff(columns);
  cout << " finite diff: " << finiteDiff << endl;
  cout << " computed: " << springJacobian(collisionPoint) << endl;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SURFACE::draw()
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " OPENGL RENDERING NOT ENABLED FOR THIS SURFACE " << endl;
}

//////////////////////////////////////////////////////////////////////
// multibody version -- adds the necessary quantities to 3x3 block 
// matrix 'blockJacobian'
//////////////////////////////////////////////////////////////////////
void SURFACE::springJacobian(const VEC3F& collisionPoint, 
                             const VEC3F& collisionVelocity,
                             const VEC3F& collisionForce,
                             const VEC3F& localPoint,
                             const MATRIX& U, 
                             BLOCK_MATRIX& systemMatrix)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " MULTIBODY JACOBIANS NOT IMPLEMENTED FOR THIS SURFACE" << endl;
}

//////////////////////////////////////////////////////////////////////
// multibody version -- adds the necessary quantities to 3x3 block 
// matrix 'blockJacobian'
//////////////////////////////////////////////////////////////////////
void SURFACE::springJacobianDebug(const VEC3F& collisionPoint, 
                             const VEC3F& collisionVelocity,
                             const VEC3F& collisionForce,
                             const VEC3F& localPoint,
                             const MATRIX& U, 
                             BLOCK_MATRIX& systemMatrix,
                             map<string, double>& timingBreakdown)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " MULTIBODY DEBUG JACOBIANS NOT IMPLEMENTED FOR THIS SURFACE" << endl;
}

//////////////////////////////////////////////////////////////////////
// return the bounding box dimensions
//////////////////////////////////////////////////////////////////////
void SURFACE::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " NOT IMPLEMENTED FOR THIS SURFACE " << endl;
}
