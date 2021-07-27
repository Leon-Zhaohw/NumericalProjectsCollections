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
// SPARSE_MATRIX_ARRAY.h: interface for the SPARSE_MATRIX_ARRAY class.
//
//////////////////////////////////////////////////////////////////////

#include <map>
#include "SPARSE_MATRIX_ARRAY.h"

#include <fstream>

//////////////////////////////////////////////////////////////////////
// Constructor for the sparse matrix
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX_ARRAY::SPARSE_MATRIX_ARRAY(int rows, int cols) :
  _rows(rows), _cols(cols)
{
  //_matrix.resize(_rows);
  _matrix = new vector<pair<int, Real> >[_rows];
}

//////////////////////////////////////////////////////////////////////
// return a reference to an entry
//////////////////////////////////////////////////////////////////////
Real& SPARSE_MATRIX_ARRAY::operator()(int row, int col)
{
  // bounds check
  assert(col >= 0);
  assert(row >= 0); 
  assert(col < _cols);
  assert(row < _rows);

  // get the row
  vector<pair<int, Real> >& currentRow = _matrix[row];

  for (unsigned int x = 0; x < currentRow.size(); x++)
    if (currentRow[x].first == col)
      return currentRow[x].second;

  currentRow.push_back(pair<int, Real>(col, 0.0));

  return currentRow.back().second;
}

//////////////////////////////////////////////////////////////////////
// Set the matrix to zero. Note this will *NOT* stomp the underlying map!
// It will instead set all current entries to zero so that we are not
// forced to reallocate the sparsity structure again.
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX_ARRAY::clear()
{
  for (int x = 0; x < _rows; x++)
    _matrix[x].clear();
}

//////////////////////////////////////////////////////////////////////
// project by the given left and right matrices, where we assume the 
// left one needs to be transposed
//////////////////////////////////////////////////////////////////////
MatrixXd SPARSE_MATRIX_ARRAY::projectVerySparse(const MatrixXd& left, const MatrixXd& right) const
{
  assert(_rows == left.rows());
  assert(_cols == right.rows());

  TIMER functionTimer(__FUNCTION__);
  const SPARSE_MATRIX_ARRAY& A = *this;
  MatrixXd final(left.cols(), right.cols());
  final.setZero();

  vector<MatrixXd> finals;
  int totalThreads = 1; 
  for (int x = 0; x < totalThreads; x++)
    finals.push_back(final);
  for (int x = 0; x < totalThreads; x++)
    finals[x].setZero();

  cout << " Projecting ... " << flush;

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int row = 0; row < _rows; row++)
  {
    int threadID = 0;
    for (unsigned int i = 0; i < _matrix[row].size(); i++)
    {
      int col = _matrix[row][i].first;
      Real entry = _matrix[row][i].second;
      
      // skip the zeros
      //
      //if (fabs(entry) < 1e-8) continue;

      for (int x = 0; x < right.cols(); x++)
      {
        Real outer = right(row, x) * entry;
        for (int y = 0; y < left.cols(); y++)
          finals[threadID](y,x) += left(col, y) * outer;
      }
    }
  }
  cout << " done. " << endl;
  
  for (int x = 0; x < totalThreads; x++)
    final += finals[x];

  return final;
}
