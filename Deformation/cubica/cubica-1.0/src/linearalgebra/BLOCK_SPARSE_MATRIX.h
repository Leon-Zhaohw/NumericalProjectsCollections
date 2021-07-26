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
// BLOCK_SPARSE_MATRIX.h: interface for the BLOCK_SPARSE_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef BLOCK_SPARSE_MATRIX_H
#define BLOCK_SPARSE_MATRIX_H

#include <SPARSE_MATRIX.h>
#include <BLOCK_VECTOR.h>

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
class BLOCK_SPARSE_MATRIX {

public:
  BLOCK_SPARSE_MATRIX();
  BLOCK_SPARSE_MATRIX(const BLOCK_SPARSE_MATRIX& A);

  // create a [blockRows x blockCols] block matrix
  BLOCK_SPARSE_MATRIX(int blockRows, int blockCols);
  
  virtual ~BLOCK_SPARSE_MATRIX();

  // set dimensions
  void resizeAndWipe(int blockRows, int blockCols);
  void resizeAndWipeBlock(int blockRow, int blockCol, int rows, int cols);
  
  // check if a block entry exists
  bool exists(int row, int col) const { return (_blocks[row * _blockCols + col] != NULL); };

  // get the matrix at this block entry
  // if the matrix does not exist, it will return NULL
  SPARSE_MATRIX* entry(int row, int col) { return _blocks[row * _blockCols + col]; };
  const SPARSE_MATRIX* constEntry(int row, int col) const { return _blocks[row * _blockCols + col]; };

  // naked accessor -- doesn't do any checks
  inline SPARSE_MATRIX& operator()(int row, int col) {
    return *_blocks[row * _blockCols + col];
  };
  BLOCK_SPARSE_MATRIX& operator+=(BLOCK_SPARSE_MATRIX& A);
  BLOCK_SPARSE_MATRIX& operator-=(BLOCK_SPARSE_MATRIX& A);
  BLOCK_SPARSE_MATRIX& operator*=(const Real alpha);
  BLOCK_SPARSE_MATRIX& operator=(BLOCK_SPARSE_MATRIX A);
  
  // add this matrix to the block at (row,col) 
  // 
  // If this is the first time this block is called, the matrix
  // will be set to the size of the matrix you pass in.
  // 
  // It is your responsibility to ensure that the matrix dimensions
  // match for subsequent calls.
  void add(const SPARSE_MATRIX& matrix, int row, int col);
  void add(const MATRIX& matrix, int row, int col);

  // subtract this matrix from the block at (row,col)
  // 
  // If this is the first time this block is called, the matrix
  // will be set to the size of the matrix you pass in.
  //
  // It is your responsibility to ensure that the matrix dimensions
  // match for subsequent calls.
  void subtract(const SPARSE_MATRIX& matrix, int row, int col);

  // set the block at (row,col) to this matrix
  // 
  // If this is the first time this block is called, the matrix
  // will be set to the size of the matrix you pass in.
  //
  // It is your responsibility to ensure that the matrix dimensions
  // match for subsequent calls.
  void equals(SPARSE_MATRIX& matrix, int row, int col);

  // accessors
  int blockRows() const { return _blockRows; };
  int blockCols() const { return _blockCols; };
  int subRows(int blockRow) const { return _subRows[blockRow]; };
  int subCols(int blockCol) const { return _subCols[blockCol]; };
  bool isEmpty() const { return _isEmpty; };
  SPARSE_MATRIX** data() { return _blocks; };

  int rows();
  int cols();

  // wipe all entries to zero
  void clear();
  
  // summed dims across all blocks
  int totalRows();
  int totalCols();

  // sum squared of all entries
  Real sum2();

  // convert to a sparse matrix
  SPARSE_MATRIX toSparseMatrix();

  // convert to a dense matrix
  MATRIX full();
  
  // populate a full matrix with the entries, and if the dimensions of a block are
  // ambiguous, force it to 1x1 with a zero entry
  MATRIX forceFull();

protected:
  SPARSE_MATRIX** _blocks;

  // dimensions of the high level 2D block array
  int _blockRows;
  int _blockCols;

  // dimensions of each of the subblocks
  int* _subRows;
  int* _subCols;

  // is zero -- nothing has been added to this matrix
  bool _isEmpty;
};

BLOCK_VECTOR operator*(BLOCK_SPARSE_MATRIX& A, BLOCK_VECTOR& v);
VECTOR operator*(BLOCK_SPARSE_MATRIX& A, VECTOR& v);
BLOCK_SPARSE_MATRIX operator+(const BLOCK_SPARSE_MATRIX& A, const BLOCK_SPARSE_MATRIX& B);

#endif
