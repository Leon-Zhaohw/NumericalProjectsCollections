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
// BLOCK_DIAGONAL_MATRIX.h: interface for the BLOCK_DIAGONAL_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#include "BLOCK_DIAGONAL_MATRIX.h"

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
BLOCK_DIAGONAL_MATRIX::BLOCK_DIAGONAL_MATRIX(MATRIX& diagonal, int total) :
  _blockRows(total), _blockCols(total), _isEmpty(false)
{
  _blocks = new MATRIX[total];
  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = diagonal.rows();
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = diagonal.cols();

  // put these matrices along the diagonal
  for (int x = 0; x < total; x++)
    _blocks[x] = diagonal;
}

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
BLOCK_DIAGONAL_MATRIX::BLOCK_DIAGONAL_MATRIX(const BLOCK_DIAGONAL_MATRIX& A)
{
  _blockRows = A.blockRows();
  _blockCols = A.blockCols();
  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = A.subRows(x);
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = A.subCols(x);
  
  for (int x = 0; x < _blockRows; x++)
    _blocks[x] = A.entry(x,x);
}

BLOCK_DIAGONAL_MATRIX::~BLOCK_DIAGONAL_MATRIX()
{
  if (_blocks)
    delete[] _blocks;
  delete[] _subRows;
  delete[] _subCols;
}

//////////////////////////////////////////////////////////////////////
// full matrix - block matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(const MATRIX& A, const BLOCK_DIAGONAL_MATRIX& B)
{
  MATRIX trans = B.transpose() * A.transpose();

  return trans.transpose();

 // test this if the above starts running out of memory 
 /* 
  vector<VECTOR> finalColumns;
  for (int x = 0; x < B.totalCols(); x++)
  {
    VECTOR product = A * B.fullColumn(x);
    finalColumns.push_back(product);
  }

  return MATRIX(finalColumns);
  */
}

//////////////////////////////////////////////////////////////////////
// block matrix - full matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(const BLOCK_DIAGONAL_MATRIX& A, const MATRIX& B)
{
  int blockCols = A.blockRows();
  int blockRows = A.blockCols();

  int totalCols = 0;
  int totalRows = 0;
  for (int x = 0; x < blockCols; x++)
    totalCols += A.subCols(x);
  for (int x = 0; x < blockRows; x++)
    totalRows += A.subRows(x);

  MATRIX final(totalRows, B.cols());

  assert(totalCols == B.rows());
 
  for (int iBlock = 0; iBlock < A.blockRows(); iBlock++) 
    for (int kBlock = 0; kBlock < A.blockCols(); kBlock++)
    {
      if (A.exists(iBlock, kBlock))
      {
        int iStart = 0;
        int kStart = 0;
        for (int x = 0; x < iBlock; x++)
          iStart += A.subRows(x);
        for (int x = 0; x < kBlock; x++)
          kStart += A.subCols(x);
        int iEnd = iStart + A.subRows(iBlock);
        int kEnd = kStart + A.subCols(kBlock);

        // make a subblock of B
        int totalSubrows= kEnd - kStart;
        MATRIX subB(totalSubrows, B.cols());
        int subrow = 0;
        for (int k = kStart; k < kEnd; k++, subrow++)
          for (int j = 0; j < B.cols(); j++)
            subB(subrow, j) = B(k,j);

        // do the fast LAPACK multiply
        MATRIX subFinal = (A.entry(iBlock, kBlock)) * subB;

        // copy the result into the final
        subrow = 0;
        for (int i = iStart; i < iEnd; i++, subrow++)
          for (int j = 0; j < B.cols(); j++)
            final(i,j) += subFinal(subrow, j);
      }
    }

  return final;
}

//////////////////////////////////////////////////////////////////////
// get the transpose
//////////////////////////////////////////////////////////////////////
BLOCK_DIAGONAL_MATRIX BLOCK_DIAGONAL_MATRIX::transpose() const
{
  BLOCK_DIAGONAL_MATRIX final(*this);
  for (int x = 0; x < final._blockRows; x++)
    final._blocks[x] = final._blocks[x].transpose();
  return final;
}
