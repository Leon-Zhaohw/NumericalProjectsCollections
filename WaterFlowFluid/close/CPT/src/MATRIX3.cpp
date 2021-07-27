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

#include "MATRIX3.h"

MATRIX3 MATRIX3::I() { return MATRIX3(VEC3F(1,0,0), VEC3F(0,1,0), VEC3F(0,0,1)); }

MATRIX3 &MATRIX3::diag(Real d)
{
  *this = 0.0;
  row[0][0] = row[1][1] = row[2][2] = d;
  return *this;
}

MATRIX3 MATRIX3::outer_product(const VEC3F& u, const VEC3F& v)
{
  MATRIX3 A;

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      A(i, j) = u[i]*v[j];

  return A;
}

MATRIX3 operator*(const MATRIX3& n, const MATRIX3& m)
{
  MATRIX3 A;

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      A(i,j) = n[i]*m.col(j);
  return A;
}
