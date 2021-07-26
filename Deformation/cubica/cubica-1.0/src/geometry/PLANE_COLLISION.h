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
// PLANE_COLLISION.h: interface for the PLANE_COLLISION class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PLANE_COLLISION_H
#define PLANE_COLLISION_H

#include <VEC3.h>
#include <MATRIX3.h>
#include "COLLISION_RESPONSE.h"

//////////////////////////////////////////////////////////////////////
// Collision forces with respect to a plane (ie the floor)
//////////////////////////////////////////////////////////////////////
class PLANE_COLLISION : public COLLISION_RESPONSE {

public:
  PLANE_COLLISION(VEC3F down, Real planePosition, Real stiffness = 100, Real damping = 0.1, Real slidingDamping = 0.1);
  ~PLANE_COLLISION() {};

  MATRIX3 forceJacobian(VEC3F& vertex, VEC3F& velocity);
  VEC3F force(VEC3F& vertex, VEC3F& velocity);
  VEC3F friction(VEC3F& vertex, VEC3F& velocity, VEC3F& force);

  Real& planePosition() { return _planePosition; };
  Real& stiffness()    { return _stiffness; };

private:
  VEC3F _down;
  VEC3F _up;
  Real _planePosition;

  Real _stiffness;
  Real _damping;
  Real _slidingDamping;
};

#endif
