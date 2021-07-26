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

#include "INCOMPLETE_CHOLESKY.h"

//////////////////////////////////////////////////////////////////////
// Constructor for the sparse matrix
//////////////////////////////////////////////////////////////////////
INCOMPLETE_CHOLESKY::INCOMPLETE_CHOLESKY(SPARSE_MATRIX& matrix) :
  _matrix(matrix), _columnIndices(NULL), _rowIndices(NULL), _firstSolve(true)
{
  // sparse indices per column, after the diagonal, not including
  // the diagonal
  _columnIndices = new vector<int>[_matrix.cols()];

  // construct the index lists
  map<pair<int,int>, Real>& data = _matrix.matrix();
  map<pair<int,int>, Real>::iterator iter;
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    int row = iter->first.first;
    int col = iter->first.second;

    if (row > col)
      _columnIndices[col].push_back(row);
  }
}

INCOMPLETE_CHOLESKY::~INCOMPLETE_CHOLESKY()
{
  delete[] _columnIndices;
}

//////////////////////////////////////////////////////////////////////
// compute Incomplete Cholesky factorization
//////////////////////////////////////////////////////////////////////
void INCOMPLETE_CHOLESKY::init()
{
  if (!_firstSolve) return;
  _firstSolve = false;

  int currentRow = -1;
  Real diagonal = 0.0;
  if (_IC.rows() != _matrix.rows() || _IC.cols() != _matrix.cols())
    _IC.resize(_matrix.rows(), _matrix.cols());
  _IC.equals(_matrix);

  // See Golub and Van Loan, 10.3.2
  int n = _matrix.rows();
  for (int k = 0; k < n; k++)
  {
    _IC(k,k) = sqrt(_IC(k,k));
    Real diagonal = 1.0 / _IC(k,k);

    for (int i = 0; i < _columnIndices[k].size(); i++)
    {
      int row = _columnIndices[k][i];
      _IC(row, k) *= diagonal;
    }

    for (int i = 0; i < _columnIndices[k].size(); i++)
    {
      int j = _columnIndices[k][i];
      Real jk = _IC(j,k);
      for (int l = 0; l < _columnIndices[j].size(); l++)
      {
        int row = _columnIndices[j][l];
        if (_IC.exists(row,k))
          _IC(row,j) -= _IC(row,k) * jk;
      }
    }
  }

  // reflect L into L^T
  map<pair<int, int>, Real>::iterator iter;
  for (iter = _IC.matrix().begin(); iter != _IC.matrix().end(); iter++)
  {
    int row = iter->first.first;
    int col = iter->first.second;

    if (col < row)
      _IC(col, row) = iter->second;
  }
}

//////////////////////////////////////////////////////////////////////
// apply the Incomplete Cholesky factorization
//////////////////////////////////////////////////////////////////////
void INCOMPLETE_CHOLESKY::solve(VECTOR& x, VECTOR& b)
{
  map<pair<int, int>, Real>::iterator iter;
  map<pair<int, int>, Real>::reverse_iterator riter;
  map<pair<int,int>, Real>& matrix = _IC.matrix();

  // forward substitute
  int currentRow = -1;
  Real sum = 0.0;
  for (iter = matrix.begin(); iter != matrix.end(); iter++)
  {
    int row = iter->first.first;
    int col = iter->first.second;

    if (row != currentRow)
    {
      // finalize the x if we've moved on to a new row
      if (currentRow != -1)
        x(currentRow) = sum / _IC(currentRow, currentRow);

      // start a new sum for this row
      currentRow = row;
      sum = b(row);
    }
    if (col < row)
      sum -= iter->second * x(col);
  }
  // finalize the last row
  if (currentRow != -1)
    x(currentRow) = sum / _IC(currentRow, currentRow);
  
  // backward substitute
  currentRow = _matrix.rows();
  for (riter = matrix.rbegin(); riter != matrix.rend(); riter++)
  {
    int row = riter->first.first;
    int col = riter->first.second;
    if (row != currentRow)
    {
      // finalize the x if we've moved on to a new row
      if (currentRow != _matrix.rows())
        x(currentRow) = sum / _IC(currentRow, currentRow);

      // start a new sum for this row
      currentRow = row;
      sum = x(row);
    }
    if (col > row)
      sum -= riter->second * x(col);
  }
  // finalize the last row
  if (currentRow != _matrix.rows())
    x(currentRow) = sum / _IC(currentRow, currentRow);
}
