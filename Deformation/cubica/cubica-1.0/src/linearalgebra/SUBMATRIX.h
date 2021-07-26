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
// SUBMATRIX.h: interface for the SUBMATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#include <MATRIX.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// SUBMATRIX is the same as MATRIX with the key difference
// that it does NOT own the memory pointed to by _matrix, 
// and thus does not try to destroy it.
//
// The sole purpose is to make submatrix manipulations more readable.
// It is assumed that SUBMATRIX points to a contiguous subset of 
// some MATRIX object, so it unfortunately cannot be used for
// subbasis matrices, because those rows are generally scattered
// throughout _Ubasis.
//////////////////////////////////////////////////////////////////////
class SUBMATRIX : public MATRIX {

public:
  SUBMATRIX(MATRIX& matrix, int startingRow, int totalRows) 
  { 
    _rows = totalRows;
    _cols = matrix.cols();
    _matrix = &(matrix.data()[startingRow * _cols]);
  };

  // set the pointer to a dummy so the parent
  // class destructor does not destroy it
  ~SUBMATRIX() { 
    _matrix = NULL;
  };

};

#endif
