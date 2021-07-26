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

#ifndef DIAGONAL_H
#define DIAGONAL_H

#include <PRECONDITIONER.h>

//////////////////////////////////////////////////////////////////////
// Diagonal (Jacobi) preconditioner
//////////////////////////////////////////////////////////////////////
class DIAGONAL : public PRECONDITIONER {

public:
  DIAGONAL(SPARSE_MATRIX& matrix);
  ~DIAGONAL();

  // compute Incomplete Cholesky factorization
  void init();
  void solve(VECTOR& x, VECTOR& b);

protected:
  VECTOR _diagonal;
  SPARSE_MATRIX& _matrix;
};

#endif
