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
#include <MATRIX3.h>
#include <iostream>

using namespace std;

MATRIX3::MATRIX3(Real* data)
{
  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
      (*this)(x,y) = data[x + y * 3];
}

MATRIX3 MATRIX3::I() { return MATRIX3(VEC3F(1,0,0), VEC3F(0,1,0), VEC3F(0,0,1)); }

MATRIX3 &MATRIX3::diag(Real d)
{
  *this = 0.0;
  row[0][0] = row[1][1] = row[2][2] = d;
  return *this;
}

MATRIX3 diag(const VEC3F& v)
{
  return MATRIX3(VEC3F(v[0],0,0),  VEC3F(0,v[1],0),  VEC3F(0,0,v[2]));
}

MATRIX3 MATRIX3::outer_product(const VEC3F& v)
{
  MATRIX3 A;
  Real x=v[0], y=v[1], z=v[2];

  A(0,0) = x*x;  A(0,1) = x*y;  A(0,2) = x*z;
  A(1,0)=A(0,1); A(1,1) = y*y;  A(1,2) = y*z;
  A(2,0)=A(0,2); A(2,1)=A(1,2); A(2,2) = z*z;

  return A;
}

MATRIX3 MATRIX3::outer_product(const VEC3F& u, const VEC3F& v)
{
  MATRIX3 A;

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      A(i, j) = u[i]*v[j];

  return A;
}

MATRIX3 MATRIX3::cross(const VEC3F& vec)
{
  MATRIX3 final;
  final(0,0) = 0.0;
  final(1,0) = vec[2];
  final(2,0) = -vec[1];

  final(0,1) = -vec[2];
  final(1,1) = 0.0;
  final(2,1) = vec[0];

  final(0,2) = vec[1];
  final(1,2) = -vec[0];
  final(2,2) = 0.0;

  return final;
}

MATRIX3 operator*(const MATRIX3& n, const MATRIX3& m)
{
  MATRIX3 A;

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      A(i,j) = n[i]*m.col(j);
  return A;
}

MATRIX3 MATRIX3::exp(VEC3F omega)
{
  Real theta = norm(omega);
  if (fabs(theta) <= 1e-8) return MATRIX3::I();

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " NOT IMPLEMENTED " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  return MATRIX3::I();
}

MATRIX3 MATRIX3::dexp(VEC3F omega)
{
  Real theta = norm(omega);
  if (fabs(theta) <= 1e-8) return MATRIX3::I();

  MATRIX3 hat = MATRIX3::cross(omega);

  MATRIX3 final = MATRIX3::I();
  final += (1.0 - cos(theta)) * hat * 1.0 / (theta * theta);
  final += (theta - sin(theta)) * hat * hat * (1.0 / (theta * theta * theta));

  return final;
}

MATRIX3 MATRIX3::cayley(VEC3F omega)
{
  Real theta = norm(omega);
  if (fabs(theta) <= 1e-4) return MATRIX3::I();

  MATRIX3 hat = MATRIX3::cross(omega);

  MATRIX3 final = MATRIX3::I();
  final += (4.0 / (4.0 + theta)) * hat;
  final += (2.0 / (4.0 + theta)) * hat * hat;

  return final;
}

void MATRIX3::read(FILE* file)
{
  row[0].read(file); 
  row[1].read(file); 
  row[2].read(file); 
}
