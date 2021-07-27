/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2013 Theodore Kim
 */
// SPARSE_MATRIX.h: interface for the SPARSE_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#include <map>
#include "SPARSE_MATRIX.h"

#include <fstream>

//////////////////////////////////////////////////////////////////////
// Constructor for the sparse matrix
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX::SPARSE_MATRIX(int rows, int cols) :
  _rows(rows), _cols(cols), _dud(0.0f)
{
}

SPARSE_MATRIX::SPARSE_MATRIX() :
  _rows(0), _cols(0), _dud(0.0f)
{
}

//////////////////////////////////////////////////////////////////////
// return a reference to an entry
//////////////////////////////////////////////////////////////////////
Real& SPARSE_MATRIX::operator()(int row, int col)
{
  // bounds check
  assert(col >= 0);
  assert(row >= 0); 
  assert(col < _cols);
  assert(row < _rows);

  // lookup the entry
  pair<int,int> index(row, col);
  map<pair<int,int>, Real>::iterator i = _matrix.find(index);
  if (i != _matrix.end())
    return i->second;

  // if it doesn't exist, create it
  _matrix[index] = 0.0f;
  return _matrix[index];
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, SPARSE_MATRIX& matrix)
{
  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator i;
  const map<pair<int,int>, Real>& data = matrix.matrix();
  for (i = data.begin(); i != data.end(); i++)
  {
    const pair<int,int> index = i->first;
    const Real value = i->second;
    out << "(" << index.first << "," << index.second << ") = " <<  value << endl;
  }
  return out;
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
VectorXd operator*(const SPARSE_MATRIX& A, const VectorXd& x) 
{
  assert(A.cols() == x.size());

  VectorXd y(A.rows());
  y.setZero();

  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator iter;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    const pair<int,int> index = iter->first;
    const Real value = iter->second;
    int i = index.first;
    int j = index.second;
    assert(j >= 0);
    assert(j < x.size());

    y[i] += x[j] * value;
  }
  return y;
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator*(const SPARSE_MATRIX& A, const VECTOR& x) 
{
  assert(A.cols() == x.size());

  VECTOR y(A.rows());

  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator iter;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    const pair<int,int> index = iter->first;
    const Real value = iter->second;
    int i = index.first;
    int j = index.second;
    y(i) += x[j] * value;
  }
  return y;
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX operator*(const Real& alpha, const SPARSE_MATRIX& A) 
{
  SPARSE_MATRIX final(A);
  final *= alpha;
  return final;
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator^(const SPARSE_MATRIX& A, const VECTOR& x) 
{
  assert(A.rows() == x.size());

  VECTOR y(A.cols());

  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator iter;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    const pair<int,int> index = iter->first;
    const Real value = iter->second;
    int i = index.first;
    int j = index.second;
    y(j) += x[i] * value;
  }
  return y;
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX operator*(SPARSE_MATRIX& A, Real& alpha) 
{
  SPARSE_MATRIX B(A);
  B *= alpha;
  return B;
}

//////////////////////////////////////////////////////////////////////
// sparse-sparse matrix subtract
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX operator-(const SPARSE_MATRIX& A, const SPARSE_MATRIX& B)
{
  assert(A.rows() == B.rows() && A.cols() == B.cols());

  SPARSE_MATRIX final(A);
  final -= B;

  return final;
}

//////////////////////////////////////////////////////////////////////
// sparse-full matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(const SPARSE_MATRIX& A, const MATRIX& B)
{
  TIMER functionTimer(__FUNCTION__);
  // create an array of row beginnings
  vector<map<pair<int, int>, Real>::const_iterator> rowBegins;

  int currentRow = 0;
  map<pair<int,int>, Real>::const_iterator iter;
  for (iter = A.matrix().begin(); iter != A.matrix().end(); iter++)
  {
    pair<int,int> index = iter->first;
    int row = index.first;

    if (row != currentRow || iter == A.matrix().begin())
    {
      rowBegins.push_back(iter);
      currentRow = row;
    }
  }

  MATRIX C(A.rows(), B.cols());
  const int rowBeginsSize = rowBegins.size();
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < rowBeginsSize; x++)
  {
    // get the beginning of the row
    map<pair<int,int>, Real>::const_iterator iter = rowBegins[x];
    int row = iter->first.first;
    int currentRow = row;

    // loop until we hit the next row
    while (row == currentRow)
    {
      pair<int,int> index = iter->first;
      int row = index.first;
      int col = index.second;
      Real entry = iter->second;
      for (int i = 0; i < B.cols(); i++)
        C(row, i) += B(col, i) * entry;

      iter++;
    }
  }
  return C;
}

//////////////////////////////////////////////////////////////////////
// sparse-full matrix multiply
//////////////////////////////////////////////////////////////////////
MatrixXd operator*(const SPARSE_MATRIX& A, const MatrixXd& B)
{
  MatrixXd C(A.rows(), B.cols());
  C.setZero();

  const map<pair<int,int>, Real>& matrix = A.matrix(); 

  map<pair<int,int>, Real>::const_iterator iter;
  for (iter = matrix.begin(); iter != matrix.end(); iter++)
  {
    pair<int,int> index = iter->first;
    int row = index.first;
    int col = index.second;
    Real entry = iter->second;

    for (int x = 0; x < B.cols(); x++)
      C(row, x) += B(col, x) * entry;
  }

  return C;
}

//////////////////////////////////////////////////////////////////////
// full-sparse matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(const MATRIX& A, const SPARSE_MATRIX& B)
{
  MATRIX C(A.rows(), B.cols());
  const map<pair<int,int>, Real>& matrix = B.matrix(); 

  map<pair<int,int>, Real>::const_iterator iter;
  for (iter = matrix.begin(); iter != matrix.end(); iter++)
  {
    pair<int,int> index = iter->first;
    int row = index.first;
    int col = index.second;
    Real entry = iter->second;

    for (int x = 0; x < A.rows(); x++)
      C(x, col) += A(x, row) * entry;
  }

  return C;
}

//////////////////////////////////////////////////////////////////////
// full^T-sparse matrix multiply
//////////////////////////////////////////////////////////////////////
MatrixXd operator^(const MatrixXd& A, const SPARSE_MATRIX& B)
{
  assert(B.rows() == A.rows());

  MatrixXd C(A.cols(), B.cols());
  C.setZero();
  const map<pair<int,int>, Real>& matrix = B.matrix(); 

  map<pair<int,int>, Real>::const_iterator iter;
  for (iter = matrix.begin(); iter != matrix.end(); iter++)
  {
    pair<int,int> index = iter->first;
    int row = index.first;
    int col = index.second;
    Real entry = iter->second;

    for (int x = 0; x < A.cols(); x++)
      C(x, col) += A(row, x) * entry;
  }

  return C;
}

//////////////////////////////////////////////////////////////////////
// full-sparse matrix multiply
//////////////////////////////////////////////////////////////////////
MatrixXd operator*(const MatrixXd& A, const SPARSE_MATRIX& B)
{
  MatrixXd C(A.rows(), B.cols());
  C.setZero();
  const map<pair<int,int>, Real>& matrix = B.matrix(); 

  map<pair<int,int>, Real>::const_iterator iter;
  for (iter = matrix.begin(); iter != matrix.end(); iter++)
  {
    pair<int,int> index = iter->first;
    int row = index.first;
    int col = index.second;
    Real entry = iter->second;

    for (int x = 0; x < A.rows(); x++)
      C(x, col) += A(x, row) * entry;
  }

  return C;
}

//////////////////////////////////////////////////////////////////////
// full^T-sparse matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator^(const MATRIX& A, const SPARSE_MATRIX& B)
{
  MATRIX C(A.cols(), B.cols());
  const map<pair<int,int>, Real>& matrix = B.matrix(); 

  map<pair<int,int>, Real>::const_iterator iter;
  for (iter = matrix.begin(); iter != matrix.end(); iter++)
  {
    pair<int,int> index = iter->first;
    int row = index.first;
    int col = index.second;
    Real entry = iter->second;

    for (int x = 0; x < A.cols(); x++)
      C(x, col) += A(row, x) * entry;
  }

  return C;
}

//////////////////////////////////////////////////////////////////////
// scale matrix by a scalar
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& SPARSE_MATRIX::operator*=(const Real& alpha) 
{
  // iterate through all the entries
  map<pair<int,int>, Real>::iterator iter;
  for (iter = _matrix.begin(); iter != _matrix.end(); iter++)
    iter->second *= alpha;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// add two sparse matrices together
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& SPARSE_MATRIX::operator+=(const SPARSE_MATRIX& A) 
{
  assert(A.rows() == _rows && A.cols() == _cols);

  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator iter;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    const pair<int,int> index = iter->first;
    const Real value = iter->second;
    int i = index.first;
    int j = index.second;
    (*this)(i,j) += value;
  }

  return *this;
}

//////////////////////////////////////////////////////////////////////
// subtact two sparse matrices
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& SPARSE_MATRIX::operator-=(const SPARSE_MATRIX& A) 
{
  assert(A.rows() == _rows && A.cols() == _cols);

  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator iter;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    const pair<int,int> index = iter->first;
    const Real value = iter->second;
    int i = index.first;
    int j = index.second;
    (*this)(i,j) -= value;
  }

  return *this;
}

//////////////////////////////////////////////////////////////////////
// Set the matrix to zero. Note this will *NOT* stomp the underlying map!
// It will instead set all current entries to zero so that we are not
// forced to reallocate the sparsity structure again.
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::clear()
{
  map<pair<int,int>, Real>::iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
    i->second = 0.0;
}

//////////////////////////////////////////////////////////////////////
// set to the identity matrix
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::setToIdentity()
{
  _matrix.clear();

  int smaller = (_rows < _cols) ? _rows : _cols;

  for (int x = 0; x < smaller; x++)
  {
    //(*this)(x,x) = 1;
    pair<int,int> index(x,x);
    _matrix[index] = 1.0;
  }
}

//////////////////////////////////////////////////////////////////////
// project by the given left and right matrices, where we assume the 
// left one needs to be transposed
//////////////////////////////////////////////////////////////////////
MatrixXd SPARSE_MATRIX::project(const MatrixXd& left, const MatrixXd& right) const
{
  MatrixXd final(left.cols(), right.cols());

  // cache for reset later
  int totalThreads = omp_get_max_threads();

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " OMP PROJECTION LIMITED " << endl;
  omp_set_num_threads(4);

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int r = 0; r < right.cols(); r++)
  {
    VectorXd column = right.col(r);
    VectorXd product = (*this) * column;

    for (int l = 0; l < left.cols(); l++)
      final(l,r) = left.col(l).dot(product);
  }

  omp_set_num_threads(totalThreads);

  return final;
}

//////////////////////////////////////////////////////////////////////
// read from a file stream
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::readGz(gzFile& file)
{
  _matrix.clear();

  // write dimensions
  gzread(file, (void*)&_rows, sizeof(int));
  gzread(file, (void*)&_cols, sizeof(int));

  // write out how many sparse entries there are
  int totalEntries;
  gzread(file, (void*)&totalEntries, sizeof(int));

  if (_rows == 0 && _cols == 0) return;
  map<pair<int,int>, Real>::const_iterator i;
  for (int i = 0; i < totalEntries; i++)
  {
    int row;
    int col;
    Real entry;

    gzread(file, (void*)&row, sizeof(int));
    gzread(file, (void*)&col, sizeof(int));
    gzread(file, (void*)&entry, sizeof(Real));

    (*this)(row,col) = entry;
  }
}

//////////////////////////////////////////////////////////////////////
// write to a file stream
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::write(FILE* file) const
{
  // write dimensions
  fwrite((void*)&_rows, sizeof(int), 1, file);
  fwrite((void*)&_cols, sizeof(int), 1, file);

  // write out how many sparse entries there are
  int totalEntries = _matrix.size();
  fwrite((void*)&totalEntries, sizeof(int), 1, file);

  if (_rows == 0 && _cols == 0) return;

  map<pair<int,int>, Real>::const_iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    const pair<int, int> index = i->first;
    const int row = index.first;
    const int col = index.second;
    const Real entry = i->second;

    fwrite((void*)&row, sizeof(int), 1, file);
    fwrite((void*)&col, sizeof(int), 1, file);
    fwrite((void*)&entry, sizeof(Real), 1, file);
  }
}

//////////////////////////////////////////////////////////////////////
// write to a file stream
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::writeGz(gzFile& file) const
{
  // write dimensions
  gzwrite(file, (void*)&_rows, sizeof(int));
  gzwrite(file, (void*)&_cols, sizeof(int));

  // write out how many sparse entries there are
  int totalEntries = _matrix.size();
  gzwrite(file, (void*)&totalEntries, sizeof(int));

  if (_rows == 0 && _cols == 0) return;

  map<pair<int,int>, Real>::const_iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    const pair<int, int> index = i->first;
    const int row = index.first;
    const int col = index.second;
    const Real entry = i->second;

    gzwrite(file, (void*)&row, sizeof(int));
    gzwrite(file, (void*)&col, sizeof(int));
    gzwrite(file, (void*)&entry, sizeof(Real));
  }
}
