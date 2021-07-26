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
// SPARSE_MATRIX.h: interface for the SPARSE_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <SETTINGS.h>
#include <MATRIX.h>
#include <vector>
#include <map>
#include <iostream>
#include <cstdio>

using namespace std;

//////////////////////////////////////////////////////////////////////
// A sparse matrix class based on maps
//////////////////////////////////////////////////////////////////////
class SPARSE_MATRIX {

public:
  SPARSE_MATRIX();
  SPARSE_MATRIX(const MATRIX& matrix);
  SPARSE_MATRIX(MATRIX& matrix, Real threshold);
  SPARSE_MATRIX(int rows, int cols);

  virtual ~SPARSE_MATRIX() {};

  // see if an entry exists in the matrix
  bool exists(int row, int col) const;

  // get the reference to an entry -
  // note that if an entry doesn't exist, it will be created and
  // set to zero.
  // 
  // to check if an entry already exists, use the exists() function
  Real& operator()(int row, int col);

  // const version of above
  Real constEntry(int row, int col) const;

  const int& rows() const { return _rows; };
  const int& cols() const { return _cols; };
  void resize(int rows, int cols) { _rows = rows; _cols = cols; };
  SPARSE_MATRIX& operator*=(const Real& alpha);
  SPARSE_MATRIX& operator+=(const SPARSE_MATRIX& A);
  SPARSE_MATRIX& operator-=(const SPARSE_MATRIX& A);

  // Set the matrix to zero. Note this will *NOT* stomp the underlying map!
  // It will instead set all current entries to zero so that we are not
  // forced to reallocate the sparsity structure again
  void clear();

  // write to a simple binary format
  void writeToBinary(string filename);

  // write the sparse matrix to Matlab format
  void writeToMatlab(string filename = string("mysparse.m"),
                     string varName = string("A"));

  // direct access to the matrix
  const map<pair<int,int>, Real>& matrix() const { return _matrix; };

  // All the matrix entries repackaged into vectors,
  // ie _matrix(rows[x], cols[x]) = values[x];
  void entries(vector<int>& rows, vector<int>& cols, vector<Real>& values);
  void entries(vector<int>& rows, vector<int>& cols, vector<const Real*>& values);
  
  // All the matrix entries repackaged into vectors,
  // ie _matrix(rows[x], cols[x]) = values[x];
  // lower or upper triangular entries only
  void entriesLower(vector<int>& rows, vector<int>& cols, vector<Real>& values);
  void entriesUpper(vector<int>& rows, vector<int>& cols, vector<Real>& values);

  // count of non-zero elements per row -- PetSc needs this to
  // efficiently preallocate
  vector<int> nonZerosPerRow();
  
  // BLAS axpy operation: B += alpha * A, where B is this matrix
  //
  // Note that axpy actually applies to vectors, but in this
  // case we can just treat the matrix as a vector and multiply
  // all its elements by alpha
  void axpy(Real alpha, SPARSE_MATRIX& A);

  // copy A into the current object
  void copies(SPARSE_MATRIX& A);
  void equals(SPARSE_MATRIX& A) { copies(A); };

  // copy A into the current object, where A(0,0) starts copying 
  // into (*this)(row,col)
  void subcopies(MATRIX& A, int row, int col);

  // add A into the current object, where A(0,0) starts copying 
  // into (*this)(row,col)
  void add(MATRIX& A, int row, int col);
  void add(Real& entry, int row, int col);

  // subtract A from the current object, where A(0,0) starts copying 
  // into (*this)(row,col)
  void subtract(MATRIX& A, int row, int col);
  void subtract(Real& entry, int row, int col);
  
  // print a specific row
  void printRow(int row);

  // convert to a full matrix
  MATRIX full() const;

  // erase a row/column
  void eraseRowColumn(int rowCol);
  
  virtual int size() { return _matrix.size(); };

  Real norm1();
  Real normInf();
  Real sum2();
  Real sum();

  bool isSymmetric();

  // return the transpose
  SPARSE_MATRIX transpose();
  
protected:
  int _rows;
  int _cols;

  // a dud Real to pass back if the index is out of bounds
  Real _dud;

  // pair is <row,col>
  map<pair<int,int>, Real> _matrix;
};

ostream& operator<<(ostream &out, SPARSE_MATRIX& matrix);
VECTOR operator*(const SPARSE_MATRIX& A, const VECTOR& x);
SPARSE_MATRIX operator*(SPARSE_MATRIX& A, Real& alpha);
MATRIX operator*(const SPARSE_MATRIX& A, const MATRIX& B);
SPARSE_MATRIX operator*(const Real& alpha, const SPARSE_MATRIX& A);
SPARSE_MATRIX operator*(const SPARSE_MATRIX& A, const SPARSE_MATRIX& B);
MATRIX operator^(MATRIX& A, SPARSE_MATRIX& B);
MATRIX operator*(const MATRIX& A, const SPARSE_MATRIX& B);
VECTOR operator^(const SPARSE_MATRIX& A, const VECTOR& x);

#endif
