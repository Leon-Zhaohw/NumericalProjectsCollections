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
// COMPOUND.cpp: implementation of the COMPOUND class.
//
//////////////////////////////////////////////////////////////////////

#include "COMPOUND.h"
#include "OBJ.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

COMPOUND::COMPOUND() :
  _bccRes(-1),
  _meshingHead(false),
  _meshingBones(false)
{
}

COMPOUND::~COMPOUND() 
{
  for (unsigned int x = 0; x < _surfaces.size(); x++)
    delete _surfaces[x];
}

bool COMPOUND::inside(float* point) {
  if (_meshingHead)
    return insideHead(point);
  if (_meshingBones)
    return insideBones(point);

  bool result = false;
  for (unsigned int x = 0; x < _surfaces.size(); x++)
    if (_surfaces[x]->inside(point))
      result = true;

  return result;
}

float COMPOUND::potential() { 
  return 0.0f; 
}

float COMPOUND::distance(float* point) {
  if (_surfaces.size() == 0)
    return 123456;

  if (_meshingHead)
    return distanceHead(point);

  float minDistance = fabs(_surfaces[0]->distance(point));
  for (unsigned int x = 1; x < _surfaces.size(); x++)
  {
    float currentDist = fabs(_surfaces[x]->distance(point));
    minDistance = (minDistance < currentDist) ? minDistance : currentDist;
  }

  for (unsigned int x = 0; x < _subtractions.size(); x++)
  {
    float currentDist = fabs(_subtractions[x]->distance(point));
    minDistance = (minDistance < currentDist) ? minDistance : currentDist;
  }

  return minDistance;
}

void COMPOUND::addSurface(SURFACE* surface) {
  _surfaces.push_back(surface);
}

void COMPOUND::subtractSurface(SURFACE* surface) {
  _subtractions.push_back(surface);
}

//////////////////////////////////////////////////////////////////////
// Normalize all the components in exactly the same way
//
// This assumes everything in a COMPOUND is an OBJ - not pretty
//////////////////////////////////////////////////////////////////////
void COMPOUND::normalize(int res)
{
  // dump everything into a big vector
  vector<vector<VEC3>*> surfaces;
  for (unsigned int x = 0; x < _subtractions.size(); x++)
    surfaces.push_back(&((OBJ*)_subtractions[x])->vertices);
  for (unsigned int x = 0; x < _surfaces.size(); x++)
    surfaces.push_back(&((OBJ*)_surfaces[x])->vertices);

  // first get the center of mass
  VEC3F centerOfMass;
  int totalVertices = 0;
  int size = surfaces.size();

  for (int x = 0; x < size; x++)
  {
    vector<VEC3>& vertices = *surfaces[x];
    if (_meshingBones && x == 1) continue;
    if (_meshingBones && x == 2) continue;
    for (unsigned int y = 0; y < vertices.size(); y++)
      centerOfMass += vertices[y];
    totalVertices += vertices.size();
  }
  centerOfMass *= 1.0 / totalVertices;
  cout << " center of mass: " << centerOfMass << endl;

  if (_meshingBones)
  {
    centerOfMass[0] = 0.0197512;
    centerOfMass[1] = 0.427389;
    centerOfMass[2] = -0.0937424;
    cout << " after clamping: " << centerOfMass << endl;
  }

  // translate everything to the center of mass
  for (unsigned int x = 0; x < surfaces.size(); x++)
  {
    vector<VEC3>& vertices = *surfaces[x];
    for (unsigned int y = 0; y < vertices.size(); y++)
      vertices[y] -= centerOfMass;
  }

  // find the maximum magnitude
  double maxVal = 0.0f;
  size = surfaces.size();
  for (int x = 0; x < size; x++)
  {
    vector<VEC3>& vertices = *surfaces[x];
    for (unsigned int y = 0; y < vertices.size(); y++)
    {
      maxVal = (fabs(vertices[y][0]) > maxVal) ? fabs(vertices[y][0]) : maxVal;
      maxVal = (fabs(vertices[y][1]) > maxVal) ? fabs(vertices[y][1]) : maxVal;
      maxVal = (fabs(vertices[y][2]) > maxVal) ? fabs(vertices[y][2]) : maxVal;
    }
  }

  // scale everything
  double scale = 0.5 - 4.0 / res;
  cout << " scale factor: " << scale / maxVal << endl;
  for (unsigned int x = 0; x < surfaces.size(); x++)
  {
    vector<VEC3>& vertices = *surfaces[x];
    for (unsigned int y = 0; y < vertices.size(); y++)

      vertices[y] *= scale / maxVal;
  }

  cout << " scale: " << scale / maxVal << endl;

  // translate everything to 0.5, 0.5, 0.5
  VEC3F half(0.5, 0.5, 0.5);
  for (unsigned int x = 0; x < surfaces.size(); x++)
  {
    vector<VEC3>& vertices = *surfaces[x];
    for (unsigned int y = 0; y < vertices.size(); y++)
      vertices[y] += half;
  }
}
//////////////////////////////////////////////////////////////////////
// Head-specific hacks
//////////////////////////////////////////////////////////////////////
bool COMPOUND::insideHead(float* point) {
  float distanceHead = _surfaces[0]->distance(point);
  //float tolerance = 1.0 / (float)_bccRes;
  float tolerance = 4.0 / (float)_bccRes;

  bool result = false;
  for (unsigned int x = 0; x < _surfaces.size(); x++)
    if (_surfaces[x]->inside(point))
      result = true;

  if (distanceHead <= 2.0 * tolerance && result == true)
    return true;

  if (result == false)
    return false;

  for (unsigned int x = 0; x < _subtractions.size(); x++)
  {
    if (_subtractions[x]->inside(point) && distanceHead > 2.0 * tolerance)
      result = false;
  }

  return result;
}

float COMPOUND::distanceHead(float* point) {
  float tolerance = 4.0 / (float)_bccRes;

  float minDistance = _surfaces[0]->distance(point);
  for (unsigned int x = 1; x < _surfaces.size(); x++)
  {
    float currentDist = _surfaces[x]->distance(point);
    minDistance = (minDistance < currentDist) ? minDistance : currentDist;
  }
  if (minDistance <= 1.0 * tolerance)
    return minDistance;

  for (unsigned int x = 0; x < _subtractions.size(); x++)
  {
    float currentDist = fabs(_subtractions[x]->distance(point));
    minDistance = (minDistance < currentDist) ? minDistance : currentDist;
  }
  return minDistance;
}

bool COMPOUND::insideBones(float* point) {

  bool insideEye   = _surfaces[1]->inside(point);
  bool insideTeeth = _surfaces[2]->inside(point);
  if (insideEye || insideTeeth) return true;

  float distanceHead = _subtractions[0]->distance(point);
  float tolerance = 4.0 / (float)_bccRes;
  if (distanceHead <= 1.0 * tolerance)
    return false;

  bool insideSkull = _surfaces[3]->inside(point);
  return insideSkull;
}

//////////////////////////////////////////////////////////////////////
// call all the constituent draws
//////////////////////////////////////////////////////////////////////
void COMPOUND::draw()
{
  for (unsigned int x = 0; x < _surfaces.size(); x++)
  {
    glColor4f(1,1,1,1);
    _surfaces[x]->draw();
  }
}

//////////////////////////////////////////////////////////////////////
// return the bounding box dimensions
//////////////////////////////////////////////////////////////////////
void COMPOUND::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  for (unsigned int x = 0; x < _surfaces.size(); x++)
  {
    VEC3F thisMin;
    VEC3F thisMax;

    _surfaces[x]->boundingBox(thisMin, thisMax);

    if (x == 0)
    {
      mins = thisMin;
      maxs = thisMax;
      continue;
    }

    for (int y = 0; y < 3; y++)
    {
      if (thisMin[y] < mins[y]) mins[y] = thisMin[y];
      if (thisMax[y] > maxs[y]) maxs[y] = thisMax[y];
    }
  }
}
