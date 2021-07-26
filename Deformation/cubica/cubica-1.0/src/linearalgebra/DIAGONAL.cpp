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
// DIAGONAL.h: interface for the DIAGONAL class.
//
//////////////////////////////////////////////////////////////////////

#include "DIAGONAL.h"

//////////////////////////////////////////////////////////////////////
// Constructor for the sparse matrix
//////////////////////////////////////////////////////////////////////
DIAGONAL::DIAGONAL(SPARSE_MATRIX& matrix) :
  _matrix(matrix),
  _diagonal(matrix.cols())
{
}

DIAGONAL::~DIAGONAL()
{
}

//////////////////////////////////////////////////////////////////////
// compute diagonal factorization
//////////////////////////////////////////////////////////////////////
void DIAGONAL::init()
{
  for (int x = 0; x < _matrix.cols(); x++)
    _diagonal[x] = 1.0 / _matrix(x,x);
}

//////////////////////////////////////////////////////////////////////
// apply the Incomplete Cholesky factorization
//////////////////////////////////////////////////////////////////////
void DIAGONAL::solve(VECTOR& x, VECTOR& b)
{
  for (int i = 0; i < b.size(); i++)
    x[i] = _diagonal[i] * b[i];
}
