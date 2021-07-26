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
// TENSOR3.h: interface for the TENSOR3 class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TENSOR3_H
#define TENSOR3_H

#include <SETTINGS.h>
#include <map>
#include <iostream>
#include <cstdio>
#include <VECTOR.h>
#include <MATRIX.h>
#include <BLOCK_MATRIX.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// A rank 3 tensor
//////////////////////////////////////////////////////////////////////
class TENSOR3 {

public:
  TENSOR3();
  TENSOR3(int rows, int cols, int slabs);
  TENSOR3(const TENSOR3& tensor);
  TENSOR3(const vector<MATRIX>& slabs);

  // a tensor repeated 'repeat' times along the diagonal
  TENSOR3(const TENSOR3& tensor, const int repeats);

  // the vector exponential derivative of a matrix
  TENSOR3(const MATRIX& A, const vector<MATRIX>& dAs);

  virtual ~TENSOR3();

  inline Real& operator()(int row, int col, int slab) {
    return _tensor[slab](row, col);
  };
  Real operator()(int row, int col, int slab) const {
    return _tensor[slab](row, col);
  };

  VEC3F dims() const { return VEC3F(_rows, _cols, _slabs); }; 
  int& rows() { return _rows; };
  int& cols() { return _cols; };
  int& slabs() { return _slabs; };
  int rows() const { return _rows; };
  int cols() const { return _cols; };
  int slabs() const { return _slabs; };
  MATRIX* data() { return _tensor; };
  const MATRIX* dataConst() const { return _tensor; };
  MATRIX& slab(int x) { return _tensor[x]; };
  const MATRIX& constSlab(int x) const { return _tensor[x]; };
  TENSOR3& operator+=(const TENSOR3& m);
  TENSOR3& operator-=(const TENSOR3& m);
  TENSOR3& operator*=(const Real& scalar);
  TENSOR3& operator=(const TENSOR3& m);

  // set the mode 3 entries according to the data
  // this isn't that safe, but it avoids having to add a new function to MATRIX
  void set(int row, int col, const Real* data);

  // scale along mode 3
  void scale(int row, int col, Real alpha);

  // wipe the whole matrix
  void clear();

  // resize the matrix and wipe to zero
  void resizeAndWipe(int rows, int cols, int slabs);

  // stomp the current tensor with the given tensor
  // starting at row number "row". It is your responsibility
  // to ensure that you don't fall off the end of this tensor.
  void setSubtensor(TENSOR3& tensor, int row);

  // squared sum of all the entries
  Real sum2();

  // transpose of the tensor -- just transpose each of the matrices (correct?)
  TENSOR3 transpose() const;

  // flip the slabs and columns of the tensor
  TENSOR3 flipSlabsAndColumns();

  // flip the slabs and rows of the tensor
  TENSOR3 flipSlabsAndRows();

  // take the product w.r.t. mode one, i.e. matrix-vector multiply with each slab
  MATRIX modeOneProduct(const VECTOR& x) const;
  void modeOneProductInplace(const VECTOR& x, MATRIX& final);
  void modeOneAxpy(const VECTOR& x, const Real scalar, MATRIX& final);
  
  // take the product w.r.t. mode one, i.e. matrix-matrix multiply with each slab,
  // but with the tensor's slabs transposed
  TENSOR3 modeOneTransposeProduct(const MATRIX& m) const;
  TENSOR3 modeOneProduct(const MATRIX& m) const;

  // take the product w.r.t. mode three, i.e. scale each slab by vector entries, and sum them
  MATRIX modeThreeProduct(const VECTOR& x);
  void modeThreeAxpy(const VECTOR& x, const Real scalar, MATRIX& final);
  void modeThreeProductInplace(const VECTOR& x, MATRIX& final);

  // get the matrix corresponding to column number 'column'
  MATRIX getColumnSlab(int column) const;

  // set the matrix corresponding to column number 'column'
  void setColumnSlab(const MATRIX& slab, int column);
  
  // add the matrix to corresponding to column number 'column'
  void addColumnSlab(const MATRIX& slab, int column);

  // compute the cross product tensor
  static TENSOR3 cross(const MATRIX& Uij);

  // compute the transpose of the cross product tensor
  static TENSOR3 crossTranspose(const MATRIX& Uij);

  void read(FILE* file);
  void write(FILE* file);

  // high performance routines
  inline void axpy(const Real alpha, const TENSOR3& RHS);

  // do an in-place subtract of the first slab from a matrix
  void subtractCollapsedSlabs(MATRIX& matrix);
  
  // do an in-place add of the first slab from a matrix
  void addCollapsedSlabs(MATRIX& matrix);
  
  // do an in-place add of the first slab from a matrix
  void axpyCollapsedSlabs(const Real& scalar, MATRIX& matrix);

  // in-place tensor-matrix multiply
  void multiplyInplace(const TENSOR3& A, const MATRIX& B);

protected:
  int _rows;
  int _cols;
  int _slabs;

  MATRIX* _tensor;

  VECTOR _workspaceColumn;
};


// These operators are ambiguous -- better to just use 
// modeOneProduct and modeThreeProduct
//
// mode 3 product (ie each slab is scaled, and they are summed)
//MATRIX operator*(TENSOR3& A, VECTOR& x);
// mode 1 product (ie matrix-vector multiply with each slab)
//MATRIX operator^(TENSOR3& A, VECTOR& x);

// unsafe version
//TENSOR3 operator*(TENSOR3& A, MATRIX& B);

// pointwise add/subtract of two tensors
TENSOR3 operator-(const TENSOR3& A, const TENSOR3& B);
TENSOR3 operator+(const TENSOR3& A, const TENSOR3& B);

// multiply each slab in A by B
TENSOR3 operator*(const TENSOR3& A, const MATRIX& B);

// multiply A by B, but cutting slabs in the other direction
TENSOR3 operator^(const TENSOR3& A, const MATRIX& B);
ostream& operator<<(ostream &out, const TENSOR3& matrix);

// mode 1 product with a block matrix
TENSOR3 operator*(const BLOCK_MATRIX& A, const TENSOR3& B);

// mode 1 transpose product with a block matrix
TENSOR3 operator^(const TENSOR3& A, const BLOCK_MATRIX& B);

// mode 1 product with a matrix
TENSOR3 operator*(const MATRIX& A, const TENSOR3& B);

// mode 1 transpose product with a matrix
TENSOR3 operator^(const MATRIX& A, const TENSOR3& B);

// scale the entire tensor
TENSOR3 operator*(const Real& a, const TENSOR3& B);

// high performance routines
inline void TENSOR3::axpy(const Real alpha, const TENSOR3& RHS)
{
  for (int x = 0; x < _slabs; x++)
    _tensor[x].axpy(alpha, RHS._tensor[x]);
}

#endif
