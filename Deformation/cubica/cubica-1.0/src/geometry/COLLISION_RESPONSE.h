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
// COLLISION_RESPONSE.h: interface for the COLLISION_RESPONSE class.
//
//////////////////////////////////////////////////////////////////////

#ifndef COLLISION_RESPONSE_H
#define COLLISION_RESPONSE_H

#include <VEC3.h>
#include <MATRIX3.h>

//////////////////////////////////////////////////////////////////////
// abstract class for returning the collision response of a
// vertex
//////////////////////////////////////////////////////////////////////
class COLLISION_RESPONSE {

public:
  COLLISION_RESPONSE() {};
  virtual ~COLLISION_RESPONSE() {};

  virtual MATRIX3 forceJacobian(VEC3F& position, VEC3F& velocity) = 0;
  virtual VEC3F force(VEC3F& position, VEC3F& velocity) = 0;
  virtual VEC3F friction(VEC3F& vertex, VEC3F& velocity, VEC3F& force) = 0;

  Real& newmarkAlpha() { return _newmarkAlpha; };

protected:
  Real _newmarkAlpha;
};

#endif

