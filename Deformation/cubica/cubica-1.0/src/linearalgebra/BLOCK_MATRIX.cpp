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

#include "BLOCK_MATRIX.h"

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX::BLOCK_MATRIX(int blockRows, int blockCols) :
  _blockRows(blockRows), _blockCols(blockCols), _isEmpty(true)
{
  _blocks = new MATRIX*[blockRows * blockCols];
  for (int x = 0; x < blockRows * blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = -1;
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = -1;
}

BLOCK_MATRIX::BLOCK_MATRIX() :
  _blocks(NULL), _blockRows(0), _blockCols(0),
  _subRows(NULL), _subCols(NULL), _isEmpty(true)
{
}

BLOCK_MATRIX::BLOCK_MATRIX(MATRIX** diagonal, int total) :
  _blockRows(total), _blockCols(total), _isEmpty(false)
{
  _blocks = new MATRIX*[_blockRows * _blockCols];
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = -1;
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = -1;

  // put these matrices along the diagonal
  for (int x = 0; x < total; x++)
    add(*(diagonal[x]), x, x);
}

BLOCK_MATRIX::BLOCK_MATRIX(vector<MATRIX>& diagonal)
{
  _blockRows = diagonal.size();
  _blockCols = diagonal.size();

  _blocks = new MATRIX*[_blockRows * _blockCols];
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = -1;
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = -1;

  // put these matrices along the diagonal
  for (unsigned int x = 0; x < diagonal.size(); x++)
    add(diagonal[x], x, x);
}


BLOCK_MATRIX::BLOCK_MATRIX(MATRIX& diagonal, int total) :
  _blockRows(total), _blockCols(total), _isEmpty(false)
{
  _blocks = new MATRIX*[_blockRows * _blockCols];
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = -1;
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = -1;

  // put these matrices along the diagonal
  for (int x = 0; x < total; x++)
    add(diagonal, x, x);
}

BLOCK_MATRIX::BLOCK_MATRIX(MATRIX3& diagonal3x3, int total) :
  _blockRows(total), _blockCols(total), _isEmpty(false)
{
  _blocks = new MATRIX*[_blockRows * _blockCols];
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = -1;
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = -1;

  // put these matrices along the diagonal
  MATRIX diagonal(diagonal3x3);
  for (int x = 0; x < total; x++)
    add(diagonal, x, x);
}

BLOCK_MATRIX::BLOCK_MATRIX(const BLOCK_MATRIX& A) :
  _blocks(NULL), _blockRows(0), _blockCols(0),
  _subRows(NULL), _subCols(NULL), _isEmpty(true)
{
  *this = A;
}

BLOCK_MATRIX::BLOCK_MATRIX(vector<BLOCK_MATRIX>& rows)
{
  assert(rows.size() > 0);

  // check all the column dims match
  _blockCols = rows[0].blockCols();
  for (unsigned int x = 0; x < rows.size(); x++)
    assert(_blockCols == rows[x].blockCols());

  // count total block rows
  _blockRows = rows[0].blockRows();
  for (unsigned int x = 1; x < rows.size(); x++)
    _blockRows += rows[x].blockRows();

  // allocate and initialize
  _blocks = new MATRIX*[_blockRows * _blockCols];
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = -1;
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = -1;

  int currentBlockRow = 0;
  for (unsigned int x = 0; x < rows.size(); x++)
  {
    BLOCK_MATRIX& row = rows[x];
    for (int i = 0; i < row.blockCols(); i++)
      for (int j = 0; j < row.blockRows(); j++)
        if (row.exists(j,i))
          this->add(*row(j,i), currentBlockRow + j, i);

    currentBlockRow += row.blockRows();
  }
}

BLOCK_MATRIX::~BLOCK_MATRIX()
{
  if (_blocks)
  {
    for (int x = 0; x < _blockRows * _blockCols; x++)
      if (_blocks[x]) delete _blocks[x];
    delete[] _blocks;
  }
  delete[] _subRows;
  delete[] _subCols;
}

BLOCK_MATRIX& BLOCK_MATRIX::operator=(const BLOCK_MATRIX& A)
{
  if (_blocks)
    for (int x = 0; x < _blockRows * _blockCols; x++)
      delete _blocks[x];
  delete[] _blocks;
  delete[] _subRows;
  delete[] _subCols;

  _blockRows = A.blockRows();
  _blockCols = A.blockCols();

  _blocks = new MATRIX*[_blockRows * _blockCols];
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = -1;
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = -1;

  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
      if (A.exists(x,y))
        this->add(A.entry(x,y), x, y);

  return *this;
}

//////////////////////////////////////////////////////////////////////
// set dimensions
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::resizeAndWipe(int blockRows, int blockCols)
{
  // stomp old contents
  if (_blocks)
  {
    for (int x = 0; x < _blockRows * _blockCols; x++)
      delete _blocks[x];
    delete[] _blocks;
  }
  delete[] _subRows;
  delete[] _subCols;

  // allocate new blocks
  _blockRows = blockRows;
  _blockCols = blockCols;
  _blocks = new MATRIX*[blockRows * blockCols];
  for (int x = 0; x < blockRows * blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = -1;
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = -1;
}

//////////////////////////////////////////////////////////////////////
// set dimensions
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::resize(MATRIX** diagonal, int total)
{ 
  // stomp old contents
  if (_blocks)
    for (int x = 0; x < _blockRows * _blockCols; x++)
      delete _blocks[x];
  delete[] _blocks;
  delete[] _subRows;
  delete[] _subCols;

  _blockRows = total;
  _blockCols = total;

  _blocks = new MATRIX*[_blockRows * _blockCols];
  for (int x = 0; x < _blockRows * _blockCols; x++)
    _blocks[x] = NULL;

  _subRows = new int[_blockRows];
  for (int x = 0; x < _blockRows; x++)
    _subRows[x] = -1;
  _subCols = new int[_blockCols];
  for (int x = 0; x < _blockCols; x++)
    _subCols[x] = -1;

  // put these matrices along the diagonal
  for (int x = 0; x < total; x++)
    add(*(diagonal[x]), x, x);
}

//////////////////////////////////////////////////////////////////////
// resize block (x,y) to size (rows, cols)
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::resizeAndWipeBlock(int x, int y, int rows, int cols)
{
  assert(x >= 0);
  assert(y >= 0);
  assert(x < _blockRows);
  assert(y < _blockCols);

  MATRIX* block = entry(x, y);
  
  if (block == NULL)
  {
    // if some other block in this row or column has already been
    // set, make sure that the current matrix dimensions match
    assert(_subRows[x] == -1 || _subRows[x] == rows);
    assert(_subCols[y] == -1 || _subCols[y] == cols);

    _blocks[x * _blockCols + y] = new MATRIX(rows, cols);
    _subRows[x] = rows;
    _subCols[y] = cols;
    _isEmpty = false;
    return;
  }

  assert(_subRows[x] == rows);
  assert(_subCols[y] == cols);
  block->resizeAndWipe(rows, cols);
}

//////////////////////////////////////////////////////////////////////
// add this matrix to the block at (row,col)
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::add(const MATRIX& matrix, int row, int col)
{
  assert(row >= 0);
  assert(col >= 0);
  assert(row < _blockRows);
  assert(col < _blockCols);

  MATRIX* block = entry(row, col);
  
  if (block == NULL)
  {
    // if some other block in this row or column has already been
    // set, make sure that the current matrix dimensions match
    assert(_subRows[row] == -1 || _subRows[row] == matrix.rows());
    assert(_subCols[col] == -1 || _subCols[col] == matrix.cols());

    _blocks[row * _blockCols + col] = new MATRIX(matrix);
    _subRows[row] = matrix.rows();
    _subCols[col] = matrix.cols();
    _isEmpty = false;
    return;
  }

  assert(matrix.rows() == _subRows[row]);
  assert(matrix.cols() == _subCols[col]);

  (*block) += matrix;
}

//////////////////////////////////////////////////////////////////////
// subtract this matrix from the block at (row,col)
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::subtract(const MATRIX& matrix, int row, int col)
{
  MATRIX* block = entry(row, col);
  
  if (block == NULL)
  {
    _blocks[row * _blockCols + col] = new MATRIX(matrix);
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
void BLOCK_MATRIX::equals(const MATRIX& matrix, int row, int col)
{
  MATRIX* block = entry(row, col);
  
  if (block == NULL)
  {
    _blocks[row * _blockCols + col] = new MATRIX(matrix);
    _subRows[row] = matrix.rows();
    _subCols[col] = matrix.cols();
    _isEmpty = false;
    return;
  }

  (*block) = matrix;
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, BLOCK_MATRIX& matrix)
{
  for (int x = 0; x < matrix.blockRows(); x++)
    for (int y = 0; y < matrix.blockCols(); y++)
    {
      out << "Block (" << x << ", " << y << ")" << endl;
      if (matrix.exists(x,y))
        out << *(matrix(x,y)) << endl;
      else
        out << "Entry is all zeros." << endl;
    }
  
  return out;
}

//////////////////////////////////////////////////////////////////////
// block matrix vector multiply
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR operator*(const BLOCK_MATRIX& A, const BLOCK_VECTOR& x)
{
  assert(x.totalBlocks() == A.blockCols());

  int blockCols = A.blockCols();
  int blockRows = A.blockRows();
  BLOCK_VECTOR y(blockRows);

  for (int i = 0; i < blockRows; i++)
  {
    VECTOR rowSum(A.subRows(i));
    for (int j = 0; j < blockCols; j++)
      // make sure the block matrix has an entry here
      if (A.exists(i,j) && x.exists(j))
      {
        //rowSum += *(A(i,j)) * (*(x(j)));
        rowSum += A.entry(i,j) * x.entry(j);
      }
    y.set(rowSum, i);
  }
  return y;
}

//////////////////////////////////////////////////////////////////////
// block matrix^T block vector multiply
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR operator^(BLOCK_MATRIX& A, BLOCK_VECTOR& v)
{
  assert(v.totalBlocks() == A.blockRows());

  BLOCK_VECTOR y(A.blockCols());

  for (int i = 0; i < A.blockCols(); i++)
    for (int j = 0; j < A.blockRows(); j++)
      if (A.exists(j,i))
      {
        VECTOR product = (*A.entry(j, i)) ^ (*v.entry(j));
        y.add(product, i);
      }
  return y;
}

//////////////////////////////////////////////////////////////////////
// block matrix^T vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator^(BLOCK_MATRIX& A, VECTOR& v)
{
  assert(A.totalRows() == v.size());
  
  // create a block vector
  BLOCK_VECTOR vBlock(A.blockRows());
  int i = 0;
  for (int x = 0; x < A.blockRows(); x++)
  {
    vBlock.resizeAndWipeBlock(x, A.subRows(x));
    VECTOR& block = *vBlock.entry(x);
    for (int y = 0; y < A.subRows(x); y++, i++)
      block[y] = v[i];
  }

  // just call the block version
  BLOCK_VECTOR result = A ^ vBlock;
  return result.full();
}

//////////////////////////////////////////////////////////////////////
// block matrix vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator*(const BLOCK_MATRIX& A, const VECTOR& v)
{
  assert(A.totalCols() == v.size());
  
  // create a block vector
  BLOCK_VECTOR vBlock(A.blockCols());
  int i = 0;
  for (int x = 0; x < A.blockCols(); x++)
  {
    vBlock.resizeAndWipeBlock(x, A.subCols(x));
    VECTOR& block = *vBlock.entry(x);
    for (int y = 0; y < A.subCols(x); y++, i++)
      block[y] = v[i];
  }

  // just call the block version
  BLOCK_VECTOR result = A * vBlock;
  return result.full();
}

//////////////////////////////////////////////////////////////////////
// block matrix-matrix multiply
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::multiply(const BLOCK_MATRIX& A, const BLOCK_MATRIX& B, BLOCK_MATRIX& C)
{
  for (int i = 0; i < A.blockRows(); i++)
    for (int j = 0; j < B.blockCols(); j++)
      for (int k = 0; k < A.blockCols(); k++)
        if (A.exists(i,k) && B.exists(k,j))
        {
          MATRIX product = (A.entry(i,k)) * (B.entry(k,j));
          C.add(product, i,j);
        }
}

//////////////////////////////////////////////////////////////////////
// block matrix-matrix multiply
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX operator*(const BLOCK_MATRIX& A, const BLOCK_MATRIX& B)
{
  BLOCK_MATRIX C(A.blockRows(), B.blockCols());

  for (int i = 0; i < A.blockRows(); i++)
    for (int j = 0; j < B.blockCols(); j++)
      for (int k = 0; k < A.blockCols(); k++)
        if (A.exists(i,k) && B.exists(k,j))
        {
          MATRIX product = (A.entry(i,k)) * (B.entry(k,j));
          C.add(product, i,j);
        }
  return C;
}

//////////////////////////////////////////////////////////////////////
// block matrix-transpose-matrix multiply
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX operator^(const BLOCK_MATRIX& A, const BLOCK_MATRIX& B)
{
  assert(A.blockRows() == B.blockRows());
  assert(A.totalRows() == B.totalRows());

  BLOCK_MATRIX C(A.blockCols(), B.blockCols());

  for (int i = 0; i < A.blockCols(); i++)
    for (int j = 0; j < B.blockCols(); j++)
      for (int k = 0; k < A.blockRows(); k++)
        if (A.exists(k,i) && B.exists(k,j))
        {
          MATRIX product = (A.entry(k,i)) ^ (B.entry(k,j));
          C.add(product, i,j);
        }
  return C;
}

//////////////////////////////////////////////////////////////////////
// block matrix-matrix multiply
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX BLOCK_MATRIX::timesTranspose(BLOCK_MATRIX& B)
{
  BLOCK_MATRIX& A = *this;
  BLOCK_MATRIX C(A.blockRows(), B.blockRows());

  for (int i = 0; i < A.blockRows(); i++)
    for (int j = 0; j < B.blockRows(); j++)
      for (int k = 0; k < A.blockCols(); k++)
        if (A.exists(i,k) && B.exists(j,k))
        {
          //MATRIX transpose = (*B.entry(j,k)).transpose();
          //MATRIX product = (*A.entry(i,k)) * transpose;
          MATRIX product = A.entry(i,k)->timesTranspose(*B.entry(j,k));
          C.add(product, i,j);
        }
  return C;
}

//////////////////////////////////////////////////////////////////////
// block matrix transpose - full matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator^(BLOCK_MATRIX& A, MATRIX& B)
{
  int totalCols = A.totalCols();
  int totalRows = A.totalRows();

  MATRIX final(totalCols, B.cols());

  assert(totalRows == B.rows());
 
  for (int i = 0; i < totalCols; i++)
    for (int j = 0; j < B.cols(); j++)
      for (int k = 0; k < totalRows; k++)
        final(i,j) += A.scalarEntry(k, i) * B(k, j);

  return final;
}

//////////////////////////////////////////////////////////////////////
// full matrix - block matrix multiply
//
// UNOPTIMIZED!!!!
//////////////////////////////////////////////////////////////////////
MATRIX operator*(const MATRIX& A, const BLOCK_MATRIX& B)
{
  return A * B.full();

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
MATRIX operator*(const BLOCK_MATRIX& A, const MATRIX& B)
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
// block matrix add
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX operator+(BLOCK_MATRIX& A, BLOCK_MATRIX& B)
{
  int blockCols = A.blockRows();
  int blockRows = A.blockCols();
  BLOCK_MATRIX C(blockRows, blockCols);

  for (int i = 0; i < blockRows; i++)
    for (int j = 0; j < blockCols; j++)
    {
      MATRIX* Amatrix = A.entry(i,j);
      MATRIX* Bmatrix = B.entry(i,j);
      if (Amatrix == NULL && Bmatrix == NULL)
        continue;
      if (Amatrix == NULL)
      {
        C.add(*Bmatrix, i, j);
        continue;
      }
      if (Bmatrix == NULL)
      {
        C.add(*Amatrix, i, j);
        continue;
      }
      C.add(*Amatrix, i,j);
      C.add(*Bmatrix, i,j);
    }
  return C;
}

//////////////////////////////////////////////////////////////////////
// block matrix add
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX& BLOCK_MATRIX::operator+=(BLOCK_MATRIX& A)
{
  assert(A._blockRows == _blockRows);
  assert(A._blockCols == _blockCols);
  
  for (int i = 0; i < _blockRows; i++)
    for (int j = 0; j < _blockCols; j++)
    {
      MATRIX* Amatrix = A.entry(i,j);
      if (Amatrix == NULL)
        continue;
      this->add(*Amatrix, i,j);
    }
  return *this;
}

//////////////////////////////////////////////////////////////////////
// block matrix subtract
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX& BLOCK_MATRIX::operator-=(BLOCK_MATRIX& A)
{
  assert(A._blockRows == _blockRows);
  assert(A._blockCols == _blockCols);
  
  for (int i = 0; i < _blockRows; i++)
    for (int j = 0; j < _blockCols; j++)
    {
      MATRIX* Amatrix = A.entry(i,j);
      if (Amatrix == NULL)
        continue;
      this->subtract(*Amatrix, i,j);
    }
  return *this;
}

//////////////////////////////////////////////////////////////////////
// block matrix scalar multiply
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX& BLOCK_MATRIX::operator*=(const Real alpha)
{
  for (int i = 0; i < _blockRows; i++)
    for (int j = 0; j < _blockCols; j++)
    {
      MATRIX* Amatrix = this->entry(i,j);
      if (Amatrix == NULL)
        continue;
      *Amatrix *= alpha;
    }
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Wipe all entries to zero
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::clear()
{
  int index = 0;
  for (int i = 0; i < _blockRows; i++)
    for (int j = 0; j < _blockCols; j++, index++)
      if (_blocks[index]) 
        _blocks[index]->clear();
}

//////////////////////////////////////////////////////////////////////
// Populate a full matrix with this block matrix
//////////////////////////////////////////////////////////////////////
MATRIX BLOCK_MATRIX::full() const
{
  /*
  int rows = 0;
  int cols = 0;
  for (int x = 0; x < _blockRows; x++)
    rows += _subRows[x];
  for (int x = 0; x < _blockCols; x++)
    cols += _subCols[x];

  MATRIX fullMatrix(rows, cols);
  for (int x = 0; x < rows; x++)
    for (int y = 0; y < cols; y++)
      fullMatrix(x,y) = scalarEntry(x,y);

  return fullMatrix;
  */
  MATRIX fullMatrix(totalRows(), totalCols());
  int currentRow = 0;
  for (int x = 0; x < _blockRows; x++)
  {
    int currentCol = 0;
    for (int y = 0; y < _blockCols; y++)
    {
      if (exists(x,y))
      {
        const MATRIX& toCopy = constEntry(x,y);
        for (int i = 0; i < toCopy.rows(); i++)
          for (int j = 0; j < toCopy.cols(); j++)
            fullMatrix(currentRow + i, currentCol + j) = toCopy(i,j);
      }
      currentCol += _subCols[y];
    }
    currentRow += _subRows[x];
  }
  return fullMatrix;
}

//////////////////////////////////////////////////////////////////////
// populate a full matrix with the entries, and if the dimensions of a block are
// ambiguous, force it to 1x1 with a zero entry
//////////////////////////////////////////////////////////////////////
MATRIX BLOCK_MATRIX::forceFull()
{
  int rows = 0;
  for (int x = 0; x < _blockRows; x++)
    rows += (_subRows[x] > 0) ? _subRows[x] : 1;
  int cols = 0;
  for (int x = 0; x < _blockCols; x++)
    cols += (_subCols[x] > 0) ? _subCols[x] : 1;

  MATRIX fullMatrix(rows, cols);
  int currentRow = 0;
  for (int x = 0; x < _blockRows; x++)
  {
    int currentCol = 0;
    for (int y = 0; y < _blockCols; y++)
    {
      if (exists(x,y))
      {
        MATRIX& toCopy = *entry(x,y);
        for (int i = 0; i < toCopy.rows(); i++)
          for (int j = 0; j < toCopy.cols(); j++)
            fullMatrix(currentRow + i, currentCol + j) = toCopy(i,j);
      }
      currentCol += _subCols[y] > 0 ? _subCols[y] : 1;
    }
    currentRow += _subRows[x] > 0 ? _subRows[x] : 1;
  }
  return fullMatrix;
}

//////////////////////////////////////////////////////////////////////
// Get the scalar entry like this is a full matrix
//////////////////////////////////////////////////////////////////////
Real BLOCK_MATRIX::scalarEntry(int row, int col) const
{
  int totalRows = 0;
  int totalCols = 0;
  int blockRow = 0;
  int blockCol = 0;

  // find which row we're in
  for (int x = 0; x < _blockRows; x++)
    if (totalRows + _subRows[x] <= row)
      totalRows += _subRows[x];
    else
    {
      blockRow = x;
      break;
    }

  // find which column we're in
  for (int x = 0; x < _blockCols; x++)
    if (totalCols + _subCols[x] <= col)
      totalCols += _subCols[x];
    else
    {
      blockCol = x;
      break;
    }

  // get the matrix
  if (!exists(blockRow, blockCol)) return 0.0;
  MATRIX& matrix = entry(blockRow, blockCol);

  // return the offset into that matrix
  return matrix(row - totalRows, col - totalCols);
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B += alpha * A, where B is this matrix
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::axpy(Real alpha, BLOCK_MATRIX& A)
{
  assert(A._blockRows == _blockRows);
  assert(A._blockCols == _blockCols);
  
  for (int i = 0; i < _blockRows; i++)
    for (int j = 0; j < _blockCols; j++)
    {
      MATRIX* Amatrix = A.entry(i,j);
      if (Amatrix == NULL)
        continue;
      MATRIX scaled = alpha * (*Amatrix);
      this->add(scaled, i,j);
    }
}

//////////////////////////////////////////////////////////////////////
// Sum squared of all entries
//////////////////////////////////////////////////////////////////////
Real BLOCK_MATRIX::sum2()
{
  Real totalSum = 0;
  for (int x = 0; x < _blockRows * _blockCols; x++)
    if (_blocks[x] != NULL)
      totalSum += _blocks[x]->sum2();

  return totalSum;
}

//////////////////////////////////////////////////////////////////////
// Summed rows across all blocks
//////////////////////////////////////////////////////////////////////
int BLOCK_MATRIX::totalRows() const
{
  int total = 0;
  for (int x = 0; x < _blockRows; x++)
  {
    // if this is tripped, one of the rows has not been initialized
    assert(_subRows[x] >= 0);
    total += _subRows[x];
  }
  return total;
}

//////////////////////////////////////////////////////////////////////
// Summed cols across all blocks
//////////////////////////////////////////////////////////////////////
int BLOCK_MATRIX::totalCols() const
{
  int total = 0;
  for (int x = 0; x < _blockCols; x++)
  {
    // if this is tripped, one of the cols has not been initialized
    assert(_subCols[x] >= 0);
    total += _subCols[x];
  }
  return total;
}

//////////////////////////////////////////////////////////////////////
// get the transpose matrix
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX BLOCK_MATRIX::transpose()
{
  BLOCK_MATRIX A(_blockCols, _blockRows);

  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
    {
      if (exists(x,y))
      {
        MATRIX trans = entry(x,y)->transpose();
        A.add(trans, y,x);
      }
    }
  return A;
}

//////////////////////////////////////////////////////////////////////
// copy upper triangle into lower
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::copyUpperToLower()
{
  assert(_blockRows == _blockCols);

  for (int x = 0; x < _blockRows; x++)
    for (int y = x + 1; y < _blockCols; y++)
    {
      if (entry(y,x))
        entry(y,x)->clear();
     
      if (entry(x,y))
      {
        MATRIX A = entry(x,y)->transpose();
        add(A, y, x);
      }
    }
}

//////////////////////////////////////////////////////////////////////
// sum squared of off-diagonal blocks
//////////////////////////////////////////////////////////////////////
Real BLOCK_MATRIX::offDiagonalSum2()
{
  Real final = 0.0;
  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
    {
      if (x == y) continue;
      if (entry(x,y)) final += entry(x,y)->sum2();
    }
  return final;
}

//////////////////////////////////////////////////////////////////////
// get the dimensions of a specific block
//////////////////////////////////////////////////////////////////////
VECTOR BLOCK_MATRIX::dims(int blockRow, int blockCol)
{
  VECTOR final(2);
  final[0] = _subRows[blockRow];
  final[1] = _subCols[blockCol];

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a vector of all the subrow dimensions
//////////////////////////////////////////////////////////////////////
VECTOR BLOCK_MATRIX::rowDims()
{
  VECTOR final(_blockRows);
  for (int x = 0; x < _blockRows; x++)
    final[x] = _subRows[x];
  return final;
}

//////////////////////////////////////////////////////////////////////
// return a vector of all the sub column dimensions
//////////////////////////////////////////////////////////////////////
VECTOR BLOCK_MATRIX::colDims()
{
  VECTOR final(_blockCols);
  for (int x = 0; x < _blockCols; x++)
    final[x] = _subCols[x];
  return final;
}

//////////////////////////////////////////////////////////////////////
// multiply vector times a specific block column of the block matrix
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR BLOCK_MATRIX::timesBlockColumn(int column, VECTOR& x)
{
  assert(column >= 0);
  assert(column < _blockCols);
  assert(x.size() == _subCols[column]);

  BLOCK_VECTOR final(_blockRows);
  for (int row  = 0; row < _blockRows; row++)
  {
    if (!exists(row, column)) continue;
    final.set(constEntry(row, column) * x, row);
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// populate a full matrix with entries over a range of block rows 
// and block columns
//////////////////////////////////////////////////////////////////////
MATRIX BLOCK_MATRIX::full(int startRow, int endRow, int startCol, int endCol)
{
  assert(startRow < endRow);
  assert(startCol < endCol);
  assert(startRow >= 0);
  assert(startCol >= 0);
  assert(startRow < _blockRows);
  assert(startCol < _blockCols);
  assert(endRow <= _blockRows);
  assert(endCol <= _blockCols);

  int totalRows = 0;
  for (int x = startRow; x < endRow; x++)
    totalRows += _subRows[x];

  int totalCols = 0;
  for (int x = startCol; x < endCol; x++)
    totalCols += _subCols[x];

  MATRIX final(totalRows, totalCols);

  int currentRow = 0;
  int currentCol = 0;
  for (int x = startRow; x < endRow; x++)
  {
    for (int y = startCol; y < endCol; y++)
    {
      if (exists(x, y))
      {
        const MATRIX& block = constEntry(x, y);
        for (int i = 0; i < block.rows(); i++)
          for (int j = 0; j < block.cols(); j++)
            final(currentRow + i, currentCol + j) = block(i,j);
      }
      currentCol += _subCols[y];
    }
    currentRow += _subRows[x];
    currentCol = 0;
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// populate a vector with the full column
//////////////////////////////////////////////////////////////////////
VECTOR BLOCK_MATRIX::fullColumn(int column) const
{
  assert(column >= 0);

  int totalRows = 0;
  for (int x = 0; x < _blockRows; x++)
    totalRows += _subRows[x];

  int blockCol = 0;
  int totalCols = 0;

  // find which block column we're in
  for (int x = 0; x < _blockCols; x++)
    if (totalCols + _subCols[x] <= column)
      totalCols += _subCols[x];
    else
    {
      blockCol = x;
      break;
    }

  VECTOR final(totalRows);

  //for (int x = 0; x < totalRows; x++)
  //  final[x] = scalarEntry(x, column);
  int i = 0;
  for (int x = 0; x < _blockRows; x++)
  {
    if (exists(x, blockCol))
      for (int y = 0; y < _subRows[x]; y++)
        final[i + y] = scalarEntry(i + y, column);
    i += _subRows[x];
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a sparse UMFPACK and TAUCS-ready representation of 
// the matrix
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::entries(vector<int>& rows, vector<int>& cols, vector<Real>& values)
{
  rows.clear();
  cols.clear();
  values.clear();

  int currentRow = 0;
  for (int x = 0; x < _blockRows; x++)
  {
    // for UMFPACK, *all* entries in a row must be enumerated
    // contiguously
    int subRows = _subRows[x];
    for (int i = 0; i < subRows; i++)
    {
      int currentCol = 0;
      for (int y = 0; y < _blockCols; y++)
      {
        if (exists(x,y))
        {
          const MATRIX& toCopy = constEntry(x,y);
          for (int j = 0; j < toCopy.cols(); j++)
          {
            rows.push_back(currentRow + i);
            cols.push_back(currentCol + j);
            values.push_back(toCopy(i,j));
          }
        }
        currentCol += _subCols[y];
      }
    }
    currentRow += _subRows[x];
  }
}

//////////////////////////////////////////////////////////////////////
// return a sparse SuperLU-ready, column major representation of 
// the matrix
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::entriesColMajor(vector<int>& rows, vector<int>& cols, vector<Real>& values)
{
  rows.clear();
  cols.clear();
  values.clear();

  int currentCol = 0;
  for (int x = 0; x < _blockCols; x++)
  {
    int subCols = _subCols[x];
    for (int i = 0; i < subCols; i++)
    {
      int currentRow = 0;
      for (int y = 0; y < _blockRows; y++)
      {
        if (exists(y,x))
        {
          const MATRIX& toCopy = constEntry(y,x);
          for (int j = 0; j < toCopy.rows(); j++)
          {
            rows.push_back(currentRow + j);
            cols.push_back(currentCol + i);
            values.push_back(toCopy(j,i));
          }
        }
        currentRow += _subRows[y];
      }
    }
    currentCol += _subCols[x];
  }
}

//////////////////////////////////////////////////////////////////////
// return a sparse SuperLU-ready, column major representation of 
// the matrix
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::entriesColMajor(double* values)
{
  int currentCol = 0;
  int entry = 0;
  for (int x = 0; x < _blockCols; x++)
  {
    int subCols = _subCols[x];
    for (int i = 0; i < subCols; i++)
    {
      int currentRow = 0;
      for (int y = 0; y < _blockRows; y++)
      {
        if (exists(y,x))
        {
          const MATRIX& toCopy = constEntry(y,x);
          for (int j = 0; j < toCopy.rows(); j++, entry++)
            values[entry] = toCopy(j,i);
        }
        currentRow += _subRows[y];
      }
    }
    currentCol += _subCols[x];
  }
}

//////////////////////////////////////////////////////////////////////
// just repopulate the values vector here -- assume it is the correct 
// size from a previous call to entries(rows, cols, values)
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::entries(vector<Real>& values)
{
  int entry = 0;
  int currentRow = 0;
  for (int x = 0; x < _blockRows; x++)
  {
    // for UMFPACK, *all* entries in a row must be enumerated
    // contiguously
    int subRows = _subRows[x];
    for (int i = 0; i < subRows; i++)
    {
      int currentCol = 0;
      for (int y = 0; y < _blockCols; y++)
      {
        if (exists(x,y))
        {
          const MATRIX& toCopy = constEntry(x,y);
          for (int j = 0; j < toCopy.cols(); j++)
          {
            values[entry] = (toCopy(i,j));
            entry++;
          }
        }
        currentCol += _subCols[y];
      }
    }
    currentRow += _subRows[x];
  }
}

//////////////////////////////////////////////////////////////////////
// convert to a sparse matrix
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX BLOCK_MATRIX::toSparseMatrix()
{
  int rows = totalRows();
  int cols = totalCols();
  SPARSE_MATRIX final(rows, cols);

  int currentRow = 0;
  for (int x = 0; x < _blockRows; x++)
  {
    int currentCol = 0;
    for (int y = 0; y < _blockCols; y++)
    {
      if (exists(x,y))
      {
        const MATRIX& toCopy = constEntry(x,y);
        for (int i = 0; i < toCopy.rows(); i++)
          for (int j = 0; j < toCopy.cols(); j++)
            final(currentRow + i, currentCol + j) = toCopy(i,j);
      }
      currentCol += _subCols[y];
    }
    currentRow += _subRows[x];
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// clamp the entries of the matrix smaller than a threshold to zero
//////////////////////////////////////////////////////////////////////
void BLOCK_MATRIX::clampToZero(const Real threshold)
{
  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
      if (exists(x,y))
        entry(x,y)->clampToZero(threshold);
}

//////////////////////////////////////////////////////////////////////
// find the max absolute entry in the whole matrix
//////////////////////////////////////////////////////////////////////
Real BLOCK_MATRIX::maxAbsEntry()
{
  Real maxFound = 0.0;
  for (int x = 0; x < _blockRows; x++)
    for (int y = 0; y < _blockCols; y++)
      if (exists(x,y))
      {
        Real currentMax = entry(x,y)->maxAbsEntry();
        if (currentMax > maxFound)
          maxFound = currentMax;        
      }
  return maxFound;
}
