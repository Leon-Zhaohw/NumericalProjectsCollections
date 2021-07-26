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
// INCOMPLETE_CHOLESKY.h: interface for the INCOMPLETE_CHOLESKY class.
//
//////////////////////////////////////////////////////////////////////

#ifndef INCOMPLETE_CHOLESKY_H
#define INCOMPLETE_CHOLESKY_H

#include <PRECONDITIONER.h>

//////////////////////////////////////////////////////////////////////
// Incomplete Cholesky preconditioner
//////////////////////////////////////////////////////////////////////
class INCOMPLETE_CHOLESKY : public PRECONDITIONER {

public:
  INCOMPLETE_CHOLESKY(SPARSE_MATRIX& matrix);
  ~INCOMPLETE_CHOLESKY();

  // compute Incomplete Cholesky factorization
  void init();
  void solve(VECTOR& x, VECTOR& b);

protected:
  SPARSE_MATRIX _IC;
  SPARSE_MATRIX& _matrix;
  vector<int>* _columnIndices;  // sparse indices per column, after the diagonal
  vector<int>* _rowIndices;     // sparse indices per row, after the diagonal

  bool _firstSolve;
};

#endif

