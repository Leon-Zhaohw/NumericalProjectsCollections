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

#ifndef SURFACE_H
#define SURFACE_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <cmath>
#include "SETTINGS.h"
#include "VEC3.h"
#include "FIELD_3D.h"
#include "VECTOR3_FIELD_3D.h"

#ifndef M_PI
#define M_PI 3.1415926535897931f
#endif

class SURFACE  
{
public:
	SURFACE() {};
	virtual ~SURFACE() {};
};

#endif
