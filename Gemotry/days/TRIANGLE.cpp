/*
QUIJIBO: Source code for the paper Symposium on Geometry Processing
         2015 paper "Quaternion Julia Set Shape Optimization"
Copyright (C) 2015  Theodore Kim

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "TRIANGLE.h"
#include <algorithm>

//////////////////////////////////////////////////////////////////////
// Constructor for triangle
//////////////////////////////////////////////////////////////////////
TRIANGLE::TRIANGLE(VEC3F* v0, VEC3F* v1, VEC3F* v2)
{
  _vertices.push_back(v0);
  _vertices.push_back(v1);
  _vertices.push_back(v2);

  _normal = cross(*v1 - *v0, *v2 - *v0);
  _area   = 0.5f * norm(_normal);
  _color = VEC3F(1, 1, 1);
  _normal.normalize();
}

//////////////////////////////////////////////////////////////////////
// Constructor for triangle
//////////////////////////////////////////////////////////////////////
TRIANGLE::TRIANGLE(TRIANGLE* triangle)
{
  _vertices.push_back(triangle->_vertices[0]);
  _vertices.push_back(triangle->_vertices[1]);
  _vertices.push_back(triangle->_vertices[2]);

  _normal = cross(*(_vertices[1]) - *(_vertices[0]), 
                  *(_vertices[2]) - *(_vertices[0]));
  _area   = 0.5f * norm(_normal);

  _color = triangle->color();
  _normal.normalize();
}

void TRIANGLE::recomputeNormal()
{
	_normal = cross(*(_vertices[1]) - *(_vertices[0]), 
                  *(_vertices[2]) - *(_vertices[0]));
  _area   = 0.5f * norm(_normal);
	_normal.normalize();
}

//////////////////////////////////////////////////////////////////////
// Constructor for triangle
//////////////////////////////////////////////////////////////////////
TRIANGLE::TRIANGLE()
{
  _vertices.resize(3);
  for (int x = 0; x < 3; x++)
    _vertices[x] = NULL;
  _color = VEC3F(1, 1, 1);
}

//////////////////////////////////////////////////////////////////////
// for equality, sort the indices first and then determine
// if they are all the same
//////////////////////////////////////////////////////////////////////
bool TRIANGLE::operator==(const TRIANGLE& RHS) const {
  vector<VEC3F*> copyLHS = _vertices;
  vector<VEC3F*> copyRHS = RHS._vertices;
  sort(copyLHS.begin(), copyLHS.end());
  sort(copyRHS.begin(), copyRHS.end());

  bool equal = true;
  for (unsigned int x =0; x < copyLHS.size(); x++)
    if (copyLHS[x] != copyRHS[x])
      equal = false;

  return equal;
}

//////////////////////////////////////////////////////////////////////
// set the current triangle equal to the RHS
//////////////////////////////////////////////////////////////////////
TRIANGLE& TRIANGLE::operator=(const TRIANGLE& RHS) {
  _vertices.clear();
  _vertices.resize(0);

  _vertices.push_back(RHS._vertices[0]);
  _vertices.push_back(RHS._vertices[1]);
  _vertices.push_back(RHS._vertices[2]);

  _normal = RHS._normal;
  _area = RHS._area;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// centroid of the triangle
//////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE::centroid()
{
  VEC3F final;
  final += *(_vertices[0]);
  final += *(_vertices[1]);
  final += *(_vertices[2]);

  final *= 1.0 / 3.0;

  return final;
}

//////////////////////////////////////////////////////////////////////
// dump triangle vertices to iostream
//////////////////////////////////////////////////////////////////////
ostream &operator<<(ostream &out, TRIANGLE& triangle)
{
  assert(triangle.vertices().size() == 3);

  for (int x = 0; x < 3; x++)
  {
    assert(triangle.vertex(x) != NULL);
    out << *(triangle.vertex(x)) << endl;
  }
  out << endl;
  return out;
}

//////////////////////////////////////////////////////////////////////
// See if the actual vertex positions are the same
//////////////////////////////////////////////////////////////////////
bool TRIANGLE::positionsEqual(TRIANGLE& RHS)
{
  int match[] = {-1, -1, -1};
  for (int x = 0; x < 3; x++)
  {
    Real original = (*_vertices[x]) * (*_vertices[x]);
    for (int y = 0; y < 3; y++)
    {
      VEC3F diff = (*_vertices[x] - *(RHS.vertex(y)));
      if ((diff * diff) / original < 1e-4)
      {
        match[x] = y;
        break;
      }
    }
  }
  for (int x = 0; x < 3; x++)
    if (match[x] < 0)
      return false;

  if (match[0] == match[1])
    return false;
  if (match[0] == match[2])
    return false;
  if (match[1] == match[2])
    return false;

  return true;
}

//////////////////////////////////////////////////////////////////////
// Find the maximum edge length of all the edges
//////////////////////////////////////////////////////////////////////
Real TRIANGLE::maxEdgeLength()
{
  VEC3F diff0 = (*_vertices[0]) - (*_vertices[1]);
  VEC3F diff1 = (*_vertices[0]) - (*_vertices[2]);
  VEC3F diff2 = (*_vertices[1]) - (*_vertices[2]);

  Real norms[] = { norm(diff0), norm(diff1), norm(diff2) };

  Real maxLength = norms[0];
  maxLength = (norms[1] > maxLength) ? norms[1] : maxLength;
  maxLength = (norms[2] > maxLength) ? norms[2] : maxLength;

  return maxLength;
}
//////////////////////////////////////////////////////////////////////
// get the bounding box
//////////////////////////////////////////////////////////////////////
void TRIANGLE::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  mins = *_vertices[0];
  maxs = *_vertices[0];

  for (int x = 1; x < 3; x++)
    for (int y = 0; y < 3; y++)
    {
      if ((*_vertices[x])[y] < mins[y])
        mins[y] = (*_vertices[x])[y];
      if ((*_vertices[x])[y] > maxs[y])
        maxs[y] = (*_vertices[x])[y];
    }
}
