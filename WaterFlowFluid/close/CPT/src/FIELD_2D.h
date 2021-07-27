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

#ifndef FIELD_2D_H
#define FIELD_2D_H

#include <cmath>
#include <string>
#include <iostream>
#include "SETTINGS.h"
#include <VEC3.h>
#include <MATRIX3.h>

#ifndef VARNAME
#define VARNAME(x) #x
#endif
#ifndef FIELDVIEW2D
#define FIELDVIEW2D(x) FIELD_2D::fieldViewer(x, VARNAME(x)); sleep(1);
#endif

using namespace std;

class FIELD_2D {
public:
  FIELD_2D();
  FIELD_2D(const int& rows, const int& cols);
  FIELD_2D(const FIELD_2D& m);
  ~FIELD_2D();

  // accessors
  inline Real& operator()(int x, int y) { return _data[y * _xRes + x]; };
  const Real operator()(int x, int y) const { return _data[y * _xRes + x]; };
  inline Real& operator[](int x) { return _data[x]; };
  const Real operator[](int x) const { return _data[x]; };
  Real* data() { return _data; };
  const int xRes() const { return _xRes; };
  const int yRes() const { return _yRes; };
  const int totalCells() const { return _totalCells; };

  // common field operations
  void clear();

  // field minimum cell index
  VEC3F minIndex();

  void resizeAndWipe(int xRes, int yRes);

  // overloaded operators
  FIELD_2D& operator=(const FIELD_2D& A);

private:
  int _xRes;
  int _yRes;
  int _totalCells;
  Real* _data;
};

#endif
