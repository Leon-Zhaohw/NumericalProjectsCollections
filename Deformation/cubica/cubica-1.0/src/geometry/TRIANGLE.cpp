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
// TRIANGLE.h: interface for the TRIANGLE class.
//
//////////////////////////////////////////////////////////////////////

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
// draw this triangle to GL
//////////////////////////////////////////////////////////////////////
void TRIANGLE::draw() {
  // recompute normals every time
  _normal = cross(*_vertices[1] - *_vertices[0], 
                  *_vertices[2] - *_vertices[0]);
  _normal.normalize();

  glBegin(GL_TRIANGLES);
#ifdef SINGLE_PRECISION
    // PetSc doesn't like glNormal3fv
    glNormal3f(_normal[0], _normal[1], _normal[2]);
    glVertex3fv(*_vertices[0]);
    glVertex3fv(*_vertices[1]);
    glVertex3fv(*_vertices[2]);
#else
    // PetSc doesn't like glNormal3fv
    glNormal3d(_normal[0], _normal[1], _normal[2]);
    glVertex3dv(*_vertices[0]);
    glVertex3dv(*_vertices[1]);
    glVertex3dv(*_vertices[2]);
#endif
  glEnd();
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

#include "TRIANGLE_TRIANGLE_INTERSECTION.cpp"

//////////////////////////////////////////////////////////////////////
// See if two triangles overlap
//////////////////////////////////////////////////////////////////////
bool TRIANGLE::intersects(TRIANGLE& RHS)
{
  // Call Moller's algorithm
  assert(_vertices.size() == 3);
  assert(RHS.vertices().size() == 3);
  vector<VEC3F*> verticesRHS = RHS.vertices();
  return NoDivTriTriIsect(*_vertices[0], *_vertices[1], *_vertices[2],
                          *verticesRHS[0], *verticesRHS[1], *verticesRHS[2]);
}

//////////////////////////////////////////////////////////////////////
// intersect with line segment
//////////////////////////////////////////////////////////////////////
bool TRIANGLE::intersects(VEC3F& start, VEC3F& end)
{
  VEC3F& a = *_vertices[0];
  VEC3F& b = *_vertices[1];
  VEC3F& c = *_vertices[2];

  VEC3F geometricNormal = (b - a) ^ (c - a);
  geometricNormal.normalize();

  VEC3F direction = end - start;
  Real length = norm(direction);
  direction.normalize();
  VEC3F diff = a - start;
  double denom = direction * geometricNormal;

  // catch divide by zero
  if (fabs(denom) <= 0.0) return false;

  double t = (diff * geometricNormal) / denom;
  if (t < 0) return false;

  VEC3F h = start + direction * t;

  VEC3F test = (b - a) ^ (h - a);
  if (test * geometricNormal < 0) return false; 
  test = (c - b) ^ (h - b);
  if (test * geometricNormal < 0) return false; 
  test = (a - c) ^ (h - c);
  if (test * geometricNormal < 0) return false; 

  return (t <= length);
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
