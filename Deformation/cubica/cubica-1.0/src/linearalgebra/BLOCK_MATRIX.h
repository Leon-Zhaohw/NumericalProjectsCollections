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
// BLOCK_MATRIX.h: interface for the BLOCK_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef BLOCK_MATRIX_H
#define BLOCK_MATRIX_H

#include <MATRIX.h>
#include <BLOCK_VECTOR.h>
#include <MATRIX3.h>
#include <SPARSE_MATRIX.h>

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
class BLOCK_MATRIX {

public:
  BLOCK_MATRIX();

  // create a [blockRows x blockCols] block matrix
  BLOCK_MATRIX(int blockRows, int blockCols);
  
  // create a block matrix with "diagonal" along the diagonal
  BLOCK_MATRIX(MATRIX** diagonal, int total);
  
  // create a block matrix with "diagonal" along the diagonal
  BLOCK_MATRIX(vector<MATRIX>& diagonal);
  
  // create a block matrix with "diagonal" *repeated* along the diagonal
  BLOCK_MATRIX(MATRIX& diagonal, int total);
  BLOCK_MATRIX(MATRIX3& diagonal3x3, int total);

  // copy constructor
  BLOCK_MATRIX(const BLOCK_MATRIX& A);

  // create a matrix with the entries in the vector stacked up on
  // top of each other
  BLOCK_MATRIX(vector<BLOCK_MATRIX>& rows);

  virtual ~BLOCK_MATRIX();

  // set dimensions
  void resizeAndWipe(int blockRows, int blockCols);
  void resize(MATRIX** diagonal, int total);
  void resizeAndWipeBlock(int x, int y, int rows, int cols);
  
  // check if a block entry exists
  bool exists(int row, int col) const { 
    assert(row >= 0);
    assert(row < _blockRows);
    assert(col >= 0);
    assert(col < _blockCols);
    return (_blocks[row * _blockCols + col] != NULL); 
  };

  // get the matrix at this block entry
  // if the matrix does not exist, it will return NULL
  MATRIX* entry(int row, int col) { return _blocks[row * _blockCols + col]; };
  MATRIX& entry(int row, int col) const { return *_blocks[row * _blockCols + col]; };
  const MATRIX& constEntry(int row, int col) const { return *_blocks[row * _blockCols + col]; };

  // get the scalar entry, as if this is a full matrix
  Real scalarEntry(int row, int col) const;

  // naked accessor
  inline MATRIX* operator()(int row, int col) {
    assert(row >= 0); assert(col >= 0);
    assert(row < _blockRows); assert(col < _blockCols);
    return _blocks[row * _blockCols + col];
  };
  BLOCK_MATRIX& operator-=(BLOCK_MATRIX& A);
  BLOCK_MATRIX& operator+=(BLOCK_MATRIX& A);
  BLOCK_MATRIX& operator*=(const Real alpha);
  BLOCK_MATRIX& operator=(const BLOCK_MATRIX& m);
  
  // add this matrix to the block at (row,col) 
  // 
  // If this is the first time this block is called, the matrix
  // will be set to the size of the matrix you pass in.
  // 
  // It is your responsibility to ensure that the matrix dimensions
  // match for subsequent calls.
  void add(const MATRIX& matrix, int row, int col);

  // subtract this matrix from the block at (row,col)
  // 
  // If this is the first time this block is called, the matrix
  // will be set to the size of the matrix you pass in.
  //
  // It is your responsibility to ensure that the matrix dimensions
  // match for subsequent calls.
  void subtract(const MATRIX& matrix, int row, int col);

  // set the block at (row,col) to this matrix
  // 
  // If this is the first time this block is called, the matrix
  // will be set to the size of the matrix you pass in.
  //
  // It is your responsibility to ensure that the matrix dimensions
  // match for subsequent calls.
  void equals(const MATRIX& matrix, int row, int col);

  // this matrix times the transpose of B
  BLOCK_MATRIX timesTranspose(BLOCK_MATRIX& B);

  // accessors
  int blockRows() const { return _blockRows; };
  int blockCols() const { return _blockCols; };
  int subRows(int blockRow) const { return _subRows[blockRow]; };
  int subCols(int blockCol) const { return _subCols[blockCol]; };
  bool isEmpty() const { return _isEmpty; };
  MATRIX** data() { return _blocks; };
  VECTOR dims(int blockRow, int blockCol);
  VECTOR rowDims();
  VECTOR colDims();
 
  // populate a full matrix with the entries
  MATRIX full() const;

  // populate a full matrix with the entries, and if the dimensions of a block are
  // ambiguous, force it to 1x1 with a zero entry
  MATRIX forceFull();

  // populate a full matrix with entries over a range of block rows and block columns
  MATRIX full(int startRow, int endRow, int startCol, int endCol);
  VECTOR fullColumn(int column) const;

  // wipe all entries to zero
  void clear();
  
  // BLAS axpy operation: B += alpha * A, where B is this matrix
  //
  // Note that axpy actually applies to vectors, but in this
  // case we can just treat the matrix as a vector and multiply
  // all its elements by alpha
  void axpy(Real alpha, BLOCK_MATRIX& A);

  Real sum2();
  Real offDiagonalSum2();

  // summed dims across all blocks
  int totalRows() const;
  int totalCols() const;

  // get the transpose matrix
  BLOCK_MATRIX transpose();

  static void multiply(const BLOCK_MATRIX& A, const BLOCK_MATRIX& B, BLOCK_MATRIX& C);

  // copy upper triangle into lower
  void copyUpperToLower();

  // multiply vector times a specific block column of the block matrix
  BLOCK_VECTOR timesBlockColumn(int column, VECTOR& x);

  // return a sparse UMFPACK and TAUCS-ready representation of the matrix
  void entries(vector<int>& rows, vector<int>& cols, vector<Real>& values);

  // return a sparse SuperLU-ready representation of the matrix
  void entriesColMajor(vector<int>& rows, vector<int>& cols, vector<Real>& values);
  void entriesColMajor(double* values);

  // just repopulate the values vector here -- assume it is the correct size from
  // a previous call to entries(rows, cols, values)
  void entries(vector<Real>& values);

  // convert to a sparse matrix
  SPARSE_MATRIX toSparseMatrix();

  // clamp the entries of the matrix smaller than a threshold to zero
  void clampToZero(const Real threshold);

  // find the max absolute entry in the whole matrix
  Real maxAbsEntry();

protected:
  MATRIX** _blocks;

  // dimensions of the high level 2D block array
  int _blockRows;
  int _blockCols;

  // dimensions of each of the subblocks
  int* _subRows;
  int* _subCols;

  // is zero -- nothing has been added to this matrix
  bool _isEmpty;
};

ostream& operator<<(ostream &out, BLOCK_MATRIX& matrix);
VECTOR operator*(const BLOCK_MATRIX& A, const VECTOR& x);
BLOCK_VECTOR operator*(const BLOCK_MATRIX& A, const BLOCK_VECTOR& v);
VECTOR operator^(BLOCK_MATRIX& A, VECTOR& x);
BLOCK_VECTOR operator^(BLOCK_MATRIX& A, BLOCK_VECTOR& x);
BLOCK_MATRIX operator+(BLOCK_MATRIX& A, BLOCK_MATRIX& B);
BLOCK_MATRIX operator*(const BLOCK_MATRIX& A, const BLOCK_MATRIX& B);
BLOCK_MATRIX operator^(const BLOCK_MATRIX& A, const BLOCK_MATRIX& B);
MATRIX operator*(const BLOCK_MATRIX& A, const MATRIX& B);

// UNOPTIMIZED!!!!
MATRIX operator*(const MATRIX& A, const BLOCK_MATRIX& B);

MATRIX operator^(BLOCK_MATRIX& A, MATRIX& B);

#endif
