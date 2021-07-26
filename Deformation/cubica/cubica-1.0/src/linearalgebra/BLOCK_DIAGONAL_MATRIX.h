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
// BLOCK_DIAGONAL_MATRIX.h: interface for the BLOCK_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef BLOCK_DIAGONAL_MATRIX_H
#define BLOCK_DIAGONAL_MATRIX_H

#include <MATRIX.h>
#include <BLOCK_VECTOR.h>
#include <MATRIX3.h>

//////////////////////////////////////////////////////////////////////
// Block matrix -- essentially a 2D storage class for a bunch of 
// matrices of the same size.
//
// The purpose is to allow matrices to be added and subtracted to
// blocks without having to explicitly check if the blocks exist
// already.
//
// There is no dimension checking, so you better know what you're
// doing.
//////////////////////////////////////////////////////////////////////
class BLOCK_DIAGONAL_MATRIX {

public:
  BLOCK_DIAGONAL_MATRIX(MATRIX& diagonal, int total);
  BLOCK_DIAGONAL_MATRIX(const BLOCK_DIAGONAL_MATRIX& A);
  virtual ~BLOCK_DIAGONAL_MATRIX();

  int blockRows() const { return _blockRows; };
  int blockCols() const { return _blockCols; };

  int subRows(int x) const { return _subRows[x]; };
  int subCols(int x) const { return _subCols[x]; };

  bool exists(int row, int col) const {
    if (row != col) return false;
    return true;
  };

  MATRIX& entry(int row, int col) const {
    assert(row == col);
    return _blocks[row];
  };
  
  BLOCK_DIAGONAL_MATRIX transpose() const;

protected:
  MATRIX* _blocks;

  // dimensions of the high level 2D block array
  int _blockRows;
  int _blockCols;

  // dimensions of each of the subblocks
  int* _subRows;
  int* _subCols;

  // is zero -- nothing has been added to this matrix
  bool _isEmpty;
};

MATRIX operator*(const MATRIX& A, const BLOCK_DIAGONAL_MATRIX& B);
MATRIX operator*(const BLOCK_DIAGONAL_MATRIX& A, const MATRIX& B);

#endif
