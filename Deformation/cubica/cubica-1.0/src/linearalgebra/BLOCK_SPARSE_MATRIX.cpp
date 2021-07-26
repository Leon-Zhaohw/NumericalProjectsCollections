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

#include "BLOCK_SPARSE_MATRIX.h"
#include <BLOCK_MATRIX.h>

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX::BLOCK_SPARSE_MATRIX(int blockRows, int blockCols) :
  _blockRows(blockRows), _blockCols(blockCols), _isEmpty(true)
{
  _blocks = new SPARSE_MATRIX*[blockRows * blockCols];
  for (int x = 0; x < blockRows * blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
}

BLOCK_SPARSE_MATRIX::BLOCK_SPARSE_MATRIX() :
  _blocks(NULL), _blockRows(0), _blockCols(0),
  _subRows(NULL), _subCols(NULL), _isEmpty(true)
{
}

BLOCK_SPARSE_MATRIX::BLOCK_SPARSE_MATRIX(const BLOCK_SPARSE_MATRIX& A) :
  _blockRows(A._blockRows), _blockCols(A._blockCols), _isEmpty(false)
{
  _blocks = new SPARSE_MATRIX*[_blockRows * _blockCols];
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];

  // copy A
  for (int x = 0; x < _blockCols; x++)
    for (int y = 0; y < _blockRows; y++)
    {
      if (A.exists(y,x))
        add(*(A.constEntry(y,x)), y, x);
    }
}

BLOCK_SPARSE_MATRIX::~BLOCK_SPARSE_MATRIX()
{
  if (_blocks)
    for (int x = 0; x < _blockRows * _blockCols; x++)
      if (_blocks[x])
        delete _blocks[x];
  delete[] _blocks;
  delete[] _subRows;
  delete[] _subCols;
}

//////////////////////////////////////////////////////////////////////
// add this matrix to the block at (row,col)
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::add(const SPARSE_MATRIX& matrix, int row, int col)
{
  SPARSE_MATRIX* block = entry(row, col);
  
  if (block == NULL)
  {
    _blocks[row * _blockCols + col] = new SPARSE_MATRIX(matrix);
    _subRows[row] = matrix.rows();
    _subCols[col] = matrix.cols();
    _isEmpty = false;
    return;
  }

  (*block) += matrix;
}

//////////////////////////////////////////////////////////////////////
// add this matrix to the block at (row,col)
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::add(const MATRIX& matrix, int row, int col)
{
  SPARSE_MATRIX* block = entry(row, col);
  
  if (block == NULL)
  {
    _blocks[row * _blockCols + col] = new SPARSE_MATRIX(matrix);
    _subRows[row] = matrix.rows();
    _subCols[col] = matrix.cols();
    _isEmpty = false;
    return;
  }

  (*block) += matrix;
}

//////////////////////////////////////////////////////////////////////
// subtract this matrix from the block at (row,col)
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::subtract(const SPARSE_MATRIX& matrix, int row, int col)
{
  SPARSE_MATRIX* block = entry(row, col);
  
  if (block == NULL)
  {
    _blocks[row * _blockCols + col] = new SPARSE_MATRIX(matrix);
    (*_blocks[row * _blockCols + col]) *= -1.0;
    _subRows[row] = matrix.rows();
    _subCols[col] = matrix.cols();
    _isEmpty = false;
    return;
  }

  (*block) -= matrix;
}

//////////////////////////////////////////////////////////////////////
// set the block at (row,col) to this matrix
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::equals(SPARSE_MATRIX& matrix, int row, int col)
{
  SPARSE_MATRIX* block = entry(row, col);
  
  if (block == NULL)
  {
    _blocks[row * _blockCols + col] = new SPARSE_MATRIX(matrix);
    _subRows[row] = matrix.rows();
    _subCols[col] = matrix.cols();
    _isEmpty = false;
    return;
  }

  (*block) = matrix;
}

//////////////////////////////////////////////////////////////////////
// block matrix block vector multiply
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR operator*(BLOCK_SPARSE_MATRIX& A, BLOCK_VECTOR& x)
{
  assert(x.totalBlocks() == A.blockCols());

  int blockCols = A.blockCols();
  int blockRows = A.blockRows();
  BLOCK_VECTOR y(blockRows);

  for (int i = 0; i < blockRows; i++)
  {
    VECTOR rowSum(A.subRows(i));
    for (int j = 0; j < blockCols; j++)
    {
      // make sure the block matrix has an entry here
      if (A.exists(i,j) && x.exists(j))
      {
        rowSum += A(i,j) * (*(x(j)));
      }
    }

    y.add(rowSum, i);
  }
  return y;
}

//////////////////////////////////////////////////////////////////////
// block sparse-block sparse ass
//////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX operator+(const BLOCK_SPARSE_MATRIX& A, const BLOCK_SPARSE_MATRIX& B)
{
  assert(A.blockRows() == B.blockRows());
  assert(A.blockCols() == B.blockCols());

  BLOCK_SPARSE_MATRIX C(A);

  for (int x = 0; x < B.blockCols(); x++)
    for (int y = 0; y < B.blockRows(); y++)
      if (B.exists(y,x))
        C.add(*(B.constEntry(y,x)), y, x);

  return C;
}

//////////////////////////////////////////////////////////////////////
// block matrix vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator*(BLOCK_SPARSE_MATRIX& A, VECTOR& x)
{
  assert(x.size() == A.cols());

  int blockCols = A.blockCols();
  int blockRows = A.blockRows();
  VECTOR y(A.cols());

  int rowOffset = 0;
  for (int i = 0; i < blockRows; i++)
  {
    int colOffset = 0;
    for (int j = 0; j < blockCols; j++)
    {
      // make sure the block matrix has an entry here
      if (A.exists(i,j))
      {
        // create a subvector for multiplication
        VECTOR subVector(A.subCols(j));
        for (int k = 0; k < A.subCols(j); k++)
          subVector[k] = x[colOffset + k];

        // do the multiply
        VECTOR product =  A(i,j) * subVector;
        assert(A.subRows(i) == product.size());

        // patch it back into the final vector
        for (int k = 0; k < product.size(); k++)
          y[rowOffset + k] += product[k];
      }
      colOffset += A.subCols(j);
    }
    rowOffset += A.subRows(i);
  }
  return y;
}


//////////////////////////////////////////////////////////////////////
// Wipe all entries to zero
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::clear()
{
  int index = 0;
  for (int i = 0; i < _blockRows; i++)
    for (int j = 0; j < _blockCols; j++, index++)
      if (_blocks[index]) 
        _blocks[index]->clear();
}

//////////////////////////////////////////////////////////////////////
// Summed rows across all blocks
//////////////////////////////////////////////////////////////////////
int BLOCK_SPARSE_MATRIX::totalRows()
{
  int total = 0;
  for (int x = 0; x < _blockRows; x++)
    total += _subRows[x];
  return total;
}

//////////////////////////////////////////////////////////////////////
// Summed cols across all blocks
//////////////////////////////////////////////////////////////////////
int BLOCK_SPARSE_MATRIX::totalCols()
{
  int total = 0;
  for (int x = 0; x < _blockCols; x++)
    total += _subCols[x];
  return total;
}

//////////////////////////////////////////////////////////////////////
// block matrix scalar multiply
//////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX& BLOCK_SPARSE_MATRIX::operator*=(const Real alpha)
{
  for (int x = 0; x < _blockRows * _blockCols; x++)
    if (_blocks[x])
      (*_blocks[x]) *= alpha;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// block matrix scalar multiply
//////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX& BLOCK_SPARSE_MATRIX::operator=(BLOCK_SPARSE_MATRIX A)
{
  if (_blocks)
    for (int x = 0; x < _blockRows * _blockCols; x++)
      delete _blocks[x];
  delete[] _blocks;
  delete[] _subRows;
  delete[] _subCols;

  _blocks = new SPARSE_MATRIX*[_blockRows * _blockCols];
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];

  // copy A
  for (int x = 0; x < _blockCols; x++)
    for (int y = 0; y < _blockRows; y++)
    {
      if (A.exists(y,x))
        add(*(A.entry(y,x)), y, x);
    }

  return *this;
}

//////////////////////////////////////////////////////////////////////
// sum squared of all entries
//////////////////////////////////////////////////////////////////////
Real BLOCK_SPARSE_MATRIX::sum2()
{
  Real final = 0.0;
  for (int x = 0; x < _blockRows * _blockCols; x++)
    if (_blocks[x])
      final += _blocks[x]->sum2();

  return final;
}

//////////////////////////////////////////////////////////////////////
// get the total rows
//////////////////////////////////////////////////////////////////////
int BLOCK_SPARSE_MATRIX::rows()
{
  int totalRows = 0;
  for (int x = 0; x < _blockRows; x++)
    totalRows += _subRows[x];

  return totalRows;
}

//////////////////////////////////////////////////////////////////////
// get the total columns
//////////////////////////////////////////////////////////////////////
int BLOCK_SPARSE_MATRIX::cols()
{
  int totalCols = 0;
  for (int x = 0; x < _blockCols; x++)
    totalCols += _subCols[x];

  return totalCols;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::resizeAndWipe(int blockRows, int blockCols)
{
  // clean up old versions
  if (_blocks)
    for (int x = 0; x < _blockRows * _blockCols; x++)
      if (_blocks[x])
        delete _blocks[x];
  delete[] _blocks;
  delete[] _subRows;
  delete[] _subCols;
  
  // resize and init
  _blockRows = blockRows;
  _blockCols = blockCols;

  _blocks = new SPARSE_MATRIX*[blockRows * blockCols];
  for (int x = 0; x < blockRows * blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
}

//////////////////////////////////////////////////////////////////////
// Convert to just a sparse matrix
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX BLOCK_SPARSE_MATRIX::toSparseMatrix()
{
  SPARSE_MATRIX sparseMatrix(totalRows(), totalCols());

  int rowOffset = 0;
  for (int x = 0; x < _blockRows; x++)
  {
    int colOffset = 0;
    for (int y = 0; y < _blockCols; y++)
    {
      if (exists(x,y))
      {
        // get the sub-block values
        vector<int> blockRows;
        vector<int> blockCols;
        vector<Real> blockValues;
        entry(x,y)->entries(blockRows, blockCols, blockValues);

        // doctor them with the row and column offsets before committing them
        for (unsigned int z = 0; z < blockRows.size(); z++)
          sparseMatrix(blockRows[z] + rowOffset, blockCols[z] + colOffset) = blockValues[z];
      }
      colOffset += _subCols[y];
    }
    rowOffset += _subRows[x];
  }

  return sparseMatrix;
}

//////////////////////////////////////////////////////////////////////
// Convert to a dense matrix
//////////////////////////////////////////////////////////////////////
MATRIX BLOCK_SPARSE_MATRIX::full()
{
  BLOCK_MATRIX blockMatrix(_blockRows, _blockCols);
  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
      if (exists(x,y))
      {
        MATRIX fullEntry = entry(x,y)->full();
        blockMatrix.add(fullEntry, x, y);
      }

  return blockMatrix.full();
}


//////////////////////////////////////////////////////////////////////
// Convert to a dense matrix
//////////////////////////////////////////////////////////////////////
MATRIX BLOCK_SPARSE_MATRIX::forceFull()
{
  BLOCK_MATRIX blockMatrix(_blockRows, _blockCols);
  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
      if (exists(x,y))
      {
        MATRIX fullEntry = entry(x,y)->full();
        blockMatrix.add(fullEntry, x, y);
      }

  return blockMatrix.forceFull();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX& BLOCK_SPARSE_MATRIX::operator+=(BLOCK_SPARSE_MATRIX& A)
{
  assert(_blockRows == A.blockRows());
  assert(_blockCols == A.blockCols());

  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
      if (A.exists(x,y))
        add(A(x,y), x,y);

  return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
BLOCK_SPARSE_MATRIX& BLOCK_SPARSE_MATRIX::operator-=(BLOCK_SPARSE_MATRIX& A)
{
  assert(_blockRows == A.blockRows());
  assert(_blockCols == A.blockCols());

  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
      if (A.exists(x,y))
        subtract(A(x,y), x,y);

  return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void BLOCK_SPARSE_MATRIX::resizeAndWipeBlock(int blockRow, int blockCol, int rows, int cols)
{
  if (exists(blockRow, blockCol))
  {
    (*this)(blockRow, blockCol).clear();
    (*this)(blockRow, blockCol).resize(rows, cols);
    return;
  }

  _blocks[blockRow * _blockCols + blockCol] = new SPARSE_MATRIX(rows, cols);
  _subRows[blockRow] = rows;
  _subCols[blockCol] = cols;
  _isEmpty = false;
}
