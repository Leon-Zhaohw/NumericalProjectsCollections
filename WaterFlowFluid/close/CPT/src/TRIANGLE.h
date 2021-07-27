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

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <SETTINGS.h>
#include <VEC3.h>
#include <vector>
#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
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
  TRIANGLE();

  VEC3F*& vertex(int x)                 { return _vertices[x]; };
  vector<VEC3F*>& vertices()            { return _vertices; };

  // get the bounding box
  void boundingBox(VEC3F& mins, VEC3F& maxs);

private:
  vector<VEC3F*> _vertices;

  VEC3F _normal;
  Real _area;
};

#endif
