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

#include "BOX.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
BOX::BOX(Real xPlus, Real xMinus, 
         Real yPlus, Real yMinus, 
         Real zPlus, Real zMinus) :
  _xPlus(xPlus), _xMinus(xMinus), _yPlus(yPlus),
  _yMinus(yMinus), _zPlus(zPlus), _zMinus(zMinus)
{
  _rotation = MATRIX3::I();
}

BOX::~BOX()
{
}

///////////////////////////////////////////////////////////////////////
// draw in GL
///////////////////////////////////////////////////////////////////////
void BOX::draw()
{
  VEC3F v000(_xMinus, _yMinus, _zMinus); 
  VEC3F v100(_xPlus, _yMinus, _zMinus); 
  VEC3F v010(_xMinus, _yPlus, _zMinus); 
  VEC3F v110(_xPlus, _yPlus, _zMinus); 
  VEC3F v001(_xMinus, _yMinus, _zPlus); 
  VEC3F v101(_xPlus, _yMinus, _zPlus); 
  VEC3F v011(_xMinus, _yPlus, _zPlus); 
  VEC3F v111(_xPlus, _yPlus, _zPlus); 


  glBegin(GL_QUADS);
    // x plus
    VEC3F normal = cross(v000 - v100, v000 - v110);
    normal.normalize();
    normal *= -1.0;
    my_glNormal3fv(normal);
    my_glVertex3fv(v010);
    my_glVertex3fv(v110);
    my_glVertex3fv(v100);
    my_glVertex3fv(v000);

    // x minus
    normal = cross(v001 - v101, v001 - v111);
    normal.normalize();
    my_glNormal3fv(normal);
    my_glVertex3fv(v001);
    my_glVertex3fv(v101);
    my_glVertex3fv(v111);
    my_glVertex3fv(v011);

    // y minus
    normal = cross(v000 - v100, v000 - v101);
    normal.normalize();
    my_glNormal3fv(normal);
    my_glVertex3fv(v000);
    my_glVertex3fv(v100);
    my_glVertex3fv(v101);
    my_glVertex3fv(v001);

    // y plus
    normal = cross(v010 - v110, v010 - v111);
    normal.normalize();
    normal *= -1.0;
    my_glNormal3fv(normal);
    my_glVertex3fv(v011);
    my_glVertex3fv(v111);
    my_glVertex3fv(v110);
    my_glVertex3fv(v010);

    // z plus
    normal = cross(v000 - v010, v000 - v011);
    normal.normalize();
    normal *= -1.0;
    my_glNormal3fv(normal);
    my_glVertex3fv(v001);
    my_glVertex3fv(v011);
    my_glVertex3fv(v010);
    my_glVertex3fv(v000);

    // z minus
    normal = cross(v100 - v110, v100 - v111);
    normal.normalize();
    my_glNormal3fv(normal);
    my_glVertex3fv(v100);
    my_glVertex3fv(v110);
    my_glVertex3fv(v111);
    my_glVertex3fv(v101);
  glEnd();
}
