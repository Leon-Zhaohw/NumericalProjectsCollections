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
//////////////////////////////////////////////////////////////////////
// Triangle - usually used for a face of a tetrahedron
//
// One of the key functions is the overloaded ==, which can tell if 
// two triangles with permuted indices are in fact the same.
//////////////////////////////////////////////////////////////////////

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <SETTINGS.h>
#include <VEC3.h>
#include <vector>
#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace std;

class TRIANGLE
{
public:
  TRIANGLE(VEC3F* v0, VEC3F* v1, VEC3F* v2);
  TRIANGLE(TRIANGLE* triangle);
  TRIANGLE();
  void draw();

  // for equality, sort the indices first and then determine
  // if they are all the same
  bool operator==(const TRIANGLE &RHS) const;
  TRIANGLE& operator=(const TRIANGLE &RHS);

  VEC3F& normal() { return _normal; };
  Real area()    { return _area; };
  Real maxEdgeLength();
  VEC3F* vertex(int x) { return _vertices[x]; };
  void setVertex(int x, VEC3F* vertex) { _vertices[x] = vertex; };
  VEC3F centroid();
  vector<VEC3F*>& vertices() { return _vertices; };

  bool positionsEqual(TRIANGLE& RHS);

  bool intersects(TRIANGLE& RHS);

  // intersect with line segment
  bool intersects(VEC3F& start, VEC3F& end);

private:
  vector<VEC3F*> _vertices;

  VEC3F _normal;
  Real _area;
};

//////////////////////////////////////////////////////////////////////
// dump triangle vertices to iostream
//////////////////////////////////////////////////////////////////////
ostream &operator<<(ostream &out, TRIANGLE& triangle);
#endif
