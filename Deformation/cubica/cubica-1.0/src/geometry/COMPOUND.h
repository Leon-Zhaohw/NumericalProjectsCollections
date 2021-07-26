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
// COMPOUND.h: interface for the COMPOUND class.
//
//////////////////////////////////////////////////////////////////////

#ifndef COMPOUND_H
#define COMPOUND_H

#include "SURFACE.h"
#include <vector>

using namespace std;

class COMPOUND : public SURFACE  
{
public:
	COMPOUND();
	virtual ~COMPOUND();

  virtual bool inside(float* point);
  virtual float potential();
  virtual float distance(float* point);

  bool insideHead(float* point);
  bool insideBones(float* point);
  float distanceHead(float* point);

  // does its own cleanup, so call with something like
  // compound->addSurface(new SPHERE(0.5f));
  void addSurface(SURFACE* surface);

  void subtractSurface(SURFACE* surface);

  int& bccRes() { return _bccRes; };
  bool& meshingHead() { return _meshingHead; };
  bool& meshingBones() { return _meshingBones; };
  vector<SURFACE*>& surfaces() { return _surfaces; };

  // normalize all the objects in exactly the same way
  void normalize(int res);

  // call all the constituent draws
  virtual void draw();
  
  // return the bounding box dimensions
  void boundingBox(VEC3F& mins, VEC3F& maxs);

private:
  vector<SURFACE*> _surfaces;
  vector<SURFACE*> _subtractions;
  inline float mymin(float i, float j) { return (i < j) ? i : j; };

  // head-specific meshing hacks
  int _bccRes;
  bool _meshingHead;
  bool _meshingBones;
};

#endif

