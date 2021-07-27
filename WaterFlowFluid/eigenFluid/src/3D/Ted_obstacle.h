/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.

 Copyright 2013 Theodore Kim
*/
#ifndef TED_OBSTACLE_H
#define TED_OBSTACLE_H

#include "Eigen"
#include "Alg/VEC3.h"
#include "Alg/MATRIX3.h"
#include "3D/FIELD_3D.h"
#include <iostream>

using namespace std;

enum OBSTACLE_FLAGS {
	EMPTY = 0, 
	MARCHED = 2, 
	RETIRED = 4 
};  

class OBSTACLE  
{
public:
	OBSTACLE() {};
	virtual ~OBSTACLE() {};

  virtual bool inside(float x, float y, float z) = 0;
  virtual void draw() = 0;
  virtual bool inside(const VEC3F& point) { return inside(point[0], point[1], point[2]); };
  virtual Real distance(const VEC3F& point) const { return 0; };

  VEC3F& translation()     { return _translation; };
  MATRIX3& rotation()      { return _rotation; };
  VEC3F& velocity()        { return _velocity; };
  VEC3F& angularVelocity() { return _angularVelocity; };

  virtual void scale(const Real alpha) = 0;
  virtual void boundingBox(VEC3F& mins, VEC3F& maxs) const = 0;

  void rotate(const MATRIX3& rotation) {
    _translation = rotation * _translation;
    _rotation = rotation * _rotation;
    _velocity = rotation * _velocity;
    _angularVelocity = rotation * _angularVelocity;
  };

  VEC3F gradient(const VEC3F& point) const {
    VEC3F copy = point;
    const Real dx = 1e-4;
    const Real invDx = 1.0 / dx;
    Real center = distance(point);
    copy[0] += dx; 
    Real right = distance(copy);
    copy = point;
    copy[1] += dx;
    Real up = distance(copy);
    copy = point;
    copy[2] += dx;
    Real out = distance(copy);

    return VEC3F((right - center) * invDx,
                 (up - center) * invDx,
                 (out - center) * invDx);
  };

  FIELD_3D distanceField(const FIELD_3D& example) {
    FIELD_3D final(example);
    final.clear();

    int xRes = final.xRes();
    int yRes = final.yRes();
    int zRes = final.zRes();

    for (int z = 0; z < zRes; z++)
      for (int y = 0; y < yRes; y++)
        for (int x = 0; x < xRes; x++)
        {
          VEC3F point = final.cellCenter(x,y,z);
          final(x,y,z) = this->distance(point);
        }
    return final;
  };

protected:
  VEC3F _translation;
  MATRIX3 _rotation;

  VEC3F _velocity;
  VEC3F _angularVelocity;
};

#endif

