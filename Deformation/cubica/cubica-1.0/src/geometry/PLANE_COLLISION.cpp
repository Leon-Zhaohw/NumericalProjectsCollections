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
// PLANE_COLLISION.cpp: implementation of the PLANE_COLLISION class.
//
//////////////////////////////////////////////////////////////////////

#include "PLANE_COLLISION.h"
#include <iostream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor for PLANE_COLLISION
//////////////////////////////////////////////////////////////////////
PLANE_COLLISION::PLANE_COLLISION(VEC3F down, Real planePosition, Real stiffness, Real damping, Real slidingDamping) :
  _down(down),
  _up(-down),
  _planePosition(planePosition),
  _stiffness(stiffness),
  _damping(damping),
  _slidingDamping(slidingDamping)
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MATRIX3 PLANE_COLLISION::forceJacobian(VEC3F& vertex, VEC3F& velocity)
{
  // do the hacky thing for the time being -- assume the
  // plane is axis-aligned and just find out which axis it is
  int axis = 0;
  if (fabs(_down[1]) > 0.5) axis = 1;
  if (fabs(_down[2]) > 0.5) axis = 2;
  
  MATRIX3 jacobian;
  if (vertex[axis] < _planePosition)
  {
    float diff = _planePosition - vertex[axis];
    VEC3F response = _up * _stiffness;

    float sign = velocity[axis] > 0.0 ? 1.0 : -1.0;
    response -= _up * _damping;

    for (int x = 0; x < 3; x++)
      jacobian(x,x) = response[x];

    /*
    for (int x = 0; x < 3; x++)
      if (x != axis)
      {
        sign = velocity[x] > 0.0 ? 1.0 : -1.0;
        jacobian(x,x) -= _newmarkAlpha * sign * _slidingDamping;
      }
      */

    /*
    cout << " newmark: " << _newmarkAlpha << endl;
    cout << " velocity: " << velocity << endl;
    cout << " jacobian diagonal: " << response << endl << endl;
    */
  }

  return jacobian;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VEC3F PLANE_COLLISION::force(VEC3F& position, VEC3F& velocity)
{
  // do the hacky thing for the time being -- assume the
  // plane is axis-aligned and just find out which axis it is
  int axis = 0;
  if (fabs(_down[1]) > 0.5) axis = 1;
  if (fabs(_down[2]) > 0.5) axis = 2;
  
  VEC3F response;
  if (position[axis] < _planePosition)
  {
    float diff = _planePosition - position[axis];
    response = _up * (diff * _stiffness);

    // The absolute value is because the dashpot should always act in the
    // direction opposite to the spring force, which is always up.
    response -= _up * (fabs(velocity[axis]) * _damping);

    /*
    for (int x = 0; x < 3; x++)
      if (x != axis)
        response[x] -= velocity[x] * _slidingDamping;
    //cout << " force response: " << response << endl;
    */
  }

  return response;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VEC3F PLANE_COLLISION::friction(VEC3F& position, VEC3F& velocity, VEC3F& force)
{
  // do the hacky thing for the time being -- assume the
  // plane is axis-aligned and just find out which axis it is
  int axis = 0;
  if (fabs(_down[1]) > 0.5) axis = 1;
  if (fabs(_down[2]) > 0.5) axis = 2;

  VEC3F response;
  if (position[axis] < _planePosition)
  {
    for (int x = 0; x < 3; x++)
      if (x != axis)
        response[x] -= velocity[x] * _slidingDamping;
  }

  //return response * norm(force);
  return response;
}
