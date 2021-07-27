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

#ifndef BOX_H
#define BOX_H

#include "SURFACE.h"
#include <MATRIX3.h>

class BOX : public SURFACE
{
public:
	BOX(Real xPlus = 0.5f, Real xMinus = -0.5f, 
      Real yPlus = 0.5f, Real yMinus = -0.5f, 
      Real zPlus = 0.5f, Real zMinus = -0.5f);
	virtual ~BOX();

  MATRIX3& rotation() { return _rotation; };

  // yes, this is lazy. deadline is looming and propagating through the inheritence tree is error-prone though.
  bool insideConst(const VEC3F& point) const;

  void stompInterior(FIELD_3D& field);
  void draw();

private:
  Real _xPlus, _xMinus;
  Real _yPlus, _yMinus;
  Real _zPlus, _zMinus;

  MATRIX3 _rotation;

  inline Real mymin(Real i, Real j) { return (i < j) ? i : j; };
};

#endif
