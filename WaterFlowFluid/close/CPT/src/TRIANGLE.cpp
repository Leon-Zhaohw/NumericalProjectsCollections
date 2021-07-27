//////////////////////////////////////////////////////////////////////
// This file is part of Closest Point Turbulence.
// 
// Closest Point Turbulence is free software: you can redistribute it 
// and/or modify it under the terms of the GNU General Public License 
// as published by the Free Software Foundation, either version 3 of 
// the License, or (at your option) any later version.
// 
// Closest Point Turbulence is distributed in the hope that it will 
// be useful, but WITHOUT ANY WARRANTY; without even the implied 
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Closest Point Turbulence. 
// If not, see <http://www.gnu.org/licenses/>.
// 
// Copyright 2013 Theodore Kim and Nils Thuerey
//////////////////////////////////////////////////////////////////////

#include "TRIANGLE.h"

//////////////////////////////////////////////////////////////////////
// Constructor for triangle
//////////////////////////////////////////////////////////////////////
TRIANGLE::TRIANGLE()
{
  _vertices.resize(3);
  for (int x = 0; x < 3; x++)
    _vertices[x] = NULL;
}

TRIANGLE::TRIANGLE(VEC3F* v0, VEC3F* v1, VEC3F* v2)
{
  _vertices.push_back(v0);
  _vertices.push_back(v1);
  _vertices.push_back(v2);

  _normal = cross(*v1 - *v0, *v2 - *v0);
  _area   = 0.5f * norm(_normal);
  _normal.normalize();
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
