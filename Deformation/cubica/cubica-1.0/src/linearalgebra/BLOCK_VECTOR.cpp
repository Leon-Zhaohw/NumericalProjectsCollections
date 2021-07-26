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
// BLOCK_VECTOR.h: interface for the BLOCK_VECTOR class.
//
//////////////////////////////////////////////////////////////////////

#include "BLOCK_VECTOR.h"

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR::BLOCK_VECTOR(int blocks) :
  _totalBlocks(blocks)
{
  _blocks = new VECTOR*[_totalBlocks];
  for (int x = 0; x < _totalBlocks; x++)
    _blocks[x] = NULL;
}

BLOCK_VECTOR::BLOCK_VECTOR(const BLOCK_VECTOR& vec)
{
  _totalBlocks = vec._totalBlocks;
  _blocks = new VECTOR*[_totalBlocks];
  for (int x = 0; x < _totalBlocks; x++)
    _blocks[x] = new VECTOR(*(vec._blocks[x]));
}

BLOCK_VECTOR::BLOCK_VECTOR(const vector<int>& blockSizes) :
  _totalBlocks(blockSizes.size())
{
  _blocks = new VECTOR*[_totalBlocks];
  for (int x = 0; x < _totalBlocks; x++)
  {
    assert(blockSizes[x] > 0);
    _blocks[x] = new VECTOR(blockSizes[x]);
  }
}

BLOCK_VECTOR::BLOCK_VECTOR(const vector<VECTOR>& blocks) :
  _totalBlocks(blocks.size())
{
  _blocks = new VECTOR*[_totalBlocks];
  for (int x = 0; x < _totalBlocks; x++)
    _blocks[x] = NULL;
  for (int x = 0; x < _totalBlocks; x++)
  {
    set(blocks[x], x);
  }
}

BLOCK_VECTOR::BLOCK_VECTOR() :
  _totalBlocks(0), _blocks(NULL)
{
}

BLOCK_VECTOR::~BLOCK_VECTOR()
{
  if (_blocks)
    for (int x = 0; x < _totalBlocks; x++)
      delete _blocks[x];
  delete[] _blocks;
}

//////////////////////////////////////////////////////////////////////
// set dimensions
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::resizeAndWipe(int blocks)
{
  // stomp old contents
  if (_blocks)
    for (int x = 0; x < _totalBlocks; x++)
      delete _blocks[x];
  delete[] _blocks;

  // allocate new blocks
  _totalBlocks = blocks;
  _blocks = new VECTOR*[_totalBlocks];
  for (int x = 0; x < _totalBlocks; x++)
    _blocks[x] = NULL;
}

//////////////////////////////////////////////////////////////////////
// set dimensions
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::resizeAndWipeBlock(int block, int size)
{
  assert(block < _totalBlocks);
  if (_blocks[block])
    delete _blocks[block];

  _blocks[block] = new VECTOR(size);
}

//////////////////////////////////////////////////////////////////////
// set dimensions
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::resize(BLOCK_VECTOR& vec)
{
  // stomp old contents
  if (_blocks)
  {
    for (int x = 0; x < _totalBlocks; x++)
      delete _blocks[x];
    delete[] _blocks;
  }

  _totalBlocks = vec._totalBlocks;
  _blocks = new VECTOR*[_totalBlocks];
  for (int x = 0; x < _totalBlocks; x++)
    _blocks[x] = new VECTOR(*(vec._blocks[x]));
}

//////////////////////////////////////////////////////////////////////
// set dimensions
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::resizeAndWipe(BLOCK_VECTOR& vec)
{
  // stomp old contents
  if (_blocks)
  {
    for (int x = 0; x < _totalBlocks; x++)
      delete _blocks[x];
    delete[] _blocks;
  }

  _totalBlocks = vec._totalBlocks;
  _blocks = new VECTOR*[_totalBlocks];
  for (int x = 0; x < _totalBlocks; x++)
    _blocks[x] = new VECTOR(vec.entry(x)->size());
}

//////////////////////////////////////////////////////////////////////
// add this matrix to the block at (whichBlock)
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::add(const VECTOR& vec, int whichBlock)
{
  assert(whichBlock >= 0);
  assert(whichBlock < _totalBlocks);

  VECTOR* block = entry(whichBlock);
  
  if (block == NULL)
  {
    _blocks[whichBlock] = new VECTOR(vec);
    return;
  }

  (*block) += vec;
}

//////////////////////////////////////////////////////////////////////
// subtract this matrix to the block at (whichBlock)
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::subtract(const VECTOR& vec, int whichBlock)
{
  assert(whichBlock >= 0);
  assert(whichBlock < _totalBlocks);

  VECTOR* block = entry(whichBlock);
  
  if (block == NULL)
  {
    _blocks[whichBlock] = new VECTOR(vec);
    *(_blocks[whichBlock]) *= -1.0;
    return;
  }

  (*block) -= vec;
}

//////////////////////////////////////////////////////////////////////
// set this block to the current vector
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::set(const VECTOR& vec, int whichBlock)
{
  assert(whichBlock >= 0);
  assert(whichBlock < _totalBlocks);

  VECTOR* block = _blocks[whichBlock];
  
  if (block == NULL)
  {
    _blocks[whichBlock] = new VECTOR(vec);
    return;
  }

  /*
  if (block->size() != vec.size())
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "Block vector sizes do not match!" << endl;
    cout << " current block size is: " << block->size() << endl;
    cout << " input block size is: " << vec.size() << endl;
  }
  */
  if (block->size() != vec.size())
    block->resizeAndWipe(vec.size());

  *(_blocks[whichBlock]) = vec;
}

//////////////////////////////////////////////////////////////////////
// set this block to the current vector
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::set(const VEC3& vec, int whichBlock)
{
  assert(whichBlock >= 0);
  assert(whichBlock < _totalBlocks);

  VECTOR* block = entry(whichBlock);
  
  if (block == NULL)
  {
    _blocks[whichBlock] = new VECTOR(3);
    block = _blocks[whichBlock];
  }

  if (block->size() != 3) 
    block->resizeAndWipe(3);

  (*block)[0] = vec[0];
  (*block)[1] = vec[1];
  (*block)[2] = vec[2];
}

//////////////////////////////////////////////////////////////////////
// compute 2 norm of this vector
//////////////////////////////////////////////////////////////////////
Real BLOCK_VECTOR::norm2()
{
  Real total = 0.0;
  for (int x = 0; x < _totalBlocks; x++)
  {
    if (_blocks[x])
      total += _blocks[x]->sum2();
  }
  return sqrt(total);
}

//////////////////////////////////////////////////////////////////////
// compute sum squared of this vector
//////////////////////////////////////////////////////////////////////
Real BLOCK_VECTOR::sum2()
{
  Real total = 0.0;
  for (int x = 0; x < _totalBlocks; x++)
  {
    if (_blocks[x])
      total += _blocks[x]->sum2();
  }
  return total;
}

//////////////////////////////////////////////////////////////////////
// block matrix vector multiply -- assumes A points to a 1D array
// of matrices that correspond to the diagonal
//
// This is pretty dangerous, so use sparingly
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR operator*(MATRIX** A, BLOCK_VECTOR& x)
{
  int totalBlocks = x.totalBlocks();
  BLOCK_VECTOR y(totalBlocks);

  for (int i = 0; i < totalBlocks; i++)
    if (x.exists(i))
    {
      VECTOR product = *(A[i]) * (*(x(i)));
      y.set(product, i);
    }
  return y;
}

//////////////////////////////////////////////////////////////////////
// Self-add
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR& BLOCK_VECTOR::operator+=(const BLOCK_VECTOR& m)
{
  assert(m.totalBlocks() == _totalBlocks);

  int totalBlocks = m.totalBlocks();
  for (int x = 0; x < totalBlocks; x++)
    if (m.exists(x))
      //add(*(m.entry(x)), x);
      add(*(m._blocks[x]), x);

  return *this;
}

//////////////////////////////////////////////////////////////////////
// Self-add
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR& BLOCK_VECTOR::operator+=(const VECTOR& m)
{
  int i = 0;
  for (int x = 0; x < _totalBlocks; x++)
    for (int y = 0; y < _blocks[x]->size(); y++, i++)
      (*_blocks[x])[y] += m[i];

  // check after the fact that the vectors were in fact the same size
  assert(i == m.size());
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Self-subtract
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR& BLOCK_VECTOR::operator-=(const VECTOR& m)
{
  int i = 0;
  for (int x = 0; x < _totalBlocks; x++)
    for (int y = 0; y < _blocks[x]->size(); y++, i++)
      (*_blocks[x])[y] -= m[i];

  // check after the fact that the vectors were in fact the same size
  assert(i == m.size());
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Self-multiply
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR& BLOCK_VECTOR::operator*=(const Real scalar)
{
  for (int x = 0; x < _totalBlocks; x++)
    if (exists(x))
      *(entry(x)) *= scalar;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// Scale block vector
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR operator*(const Real& scalar, const BLOCK_VECTOR& x)
{
  BLOCK_VECTOR final(x);

  final *= scalar;

  return final;
}

//////////////////////////////////////////////////////////////////////
// Self-subtract
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR& BLOCK_VECTOR::operator-=(const BLOCK_VECTOR& m)
{
  assert(_totalBlocks == m._totalBlocks);

  int totalBlocks = m._totalBlocks;
  for (int x = 0; x < totalBlocks; x++)
    if (m.exists(x))
      //subtract(*(m.entry(x)), x);
      subtract(*(m._blocks[x]), x);

  return *this;
}

//////////////////////////////////////////////////////////////////////
// subtract
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR operator-(const BLOCK_VECTOR& x, const BLOCK_VECTOR& y)
{
  assert(x.totalBlocks() == y.totalBlocks());

  int totalBlocks = x.totalBlocks();
  BLOCK_VECTOR final(totalBlocks);
  for (int i = 0; i < totalBlocks; i++)
  {
    const VECTOR* left = x.constEntry(i);
    const VECTOR* right = y.constEntry(i);

    if (left && right)
    {
      VECTOR diff = (*left) - (*right);
      final.add(diff, i);
    }
    else if (left)
      final.add(*left, i);
    else if (right)
      final.subtract(*right, i);
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// add
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR operator+(const BLOCK_VECTOR& x, const BLOCK_VECTOR& y)
{
  assert(x.totalBlocks() == y.totalBlocks());

  int totalBlocks = x.totalBlocks();
  BLOCK_VECTOR final(totalBlocks);
  for (int i = 0; i < totalBlocks; i++)
  {
    const VECTOR* left = x.constEntry(i);
    const VECTOR* right = y.constEntry(i);

    if (left && right)
    {
      VECTOR sum = (*left) + (*right);
      final.add(sum, i);
    }
    else if (left)
      final.add(*left, i);
    else if (right)
      final.add(*right, i);
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// Equality
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR& BLOCK_VECTOR::operator=(const BLOCK_VECTOR& m)
{
  if (_totalBlocks != m._totalBlocks)
  {
    delete[] _blocks;
    _totalBlocks = m._totalBlocks;
    _blocks = new VECTOR*[_totalBlocks];

    for (int x = 0; x < _totalBlocks; x++)
      _blocks[x] = NULL;
  }
  
  int totalBlocks = m._totalBlocks;
  for (int x = 0; x < totalBlocks; x++)
    set(*(m._blocks[x]), x);

  return *this;
}

//////////////////////////////////////////////////////////////////////
// Equality
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR& BLOCK_VECTOR::operator=(VECTOR& m)
{
  int totalEntries = 0;
  for (int x = 0; x < _totalBlocks; x++)
    totalEntries += _blocks[x]->size();

  assert(totalEntries == m.size());
 
  int i = 0;
  for (int x = 0; x < _totalBlocks; x++)
    for (int y = 0; y < _blocks[x]->size(); y++, i++)
      (*(_blocks[x]))[y] = m[i];

  return *this;
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, BLOCK_VECTOR& vec)
{
  out << endl;
  for (int x = 0; x < vec.totalBlocks(); x++)
  {
    out << "Block (" << x << ")" << endl;
    if (vec.exists(x))
      out << *(vec(x)) << endl;
    else
      out << "Entry is all zeros." << endl;
  }
  
  return out;
}

//////////////////////////////////////////////////////////////////////
// wipe all the sub-blocks
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::clear()
{
  for (int x = 0; x < _totalBlocks; x++)
    if (_blocks[x])
      _blocks[x]->clear();
}

//////////////////////////////////////////////////////////////////////
// Dump everything to a big vector
//////////////////////////////////////////////////////////////////////
VECTOR BLOCK_VECTOR::full()
{
  int size = 0;
  for (int x = 0; x < _totalBlocks; x++)
  {
    assert(_blocks[x] != NULL);
    size += _blocks[x]->size();
  }

  VECTOR fullVector(size);
  int i = 0;
  for (int x = 0; x < _totalBlocks; x++)
    for (int y = 0; y < _blocks[x]->size(); y++, i++)
      fullVector(i) = (*_blocks[x])(y);

  return fullVector;
}

//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::axpy(Real scalar, BLOCK_VECTOR& m)
{
  for (int x = 0; x < _totalBlocks; x++)
    _blocks[x]->axpy(scalar, *m(x));
}

//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::clearingAxpy(Real scalar, BLOCK_VECTOR& m)
{
  for (int x = 0; x < _totalBlocks; x++)
    _blocks[x]->clearingAxpy(scalar, *m(x));
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
Real operator^(BLOCK_VECTOR& m, BLOCK_VECTOR& n)
{
  assert(m.totalBlocks() == n.totalBlocks());
  
  Real alpha = 0.0;
  for (int x = 0; x < m.totalBlocks(); x++)
    alpha += (*m(x)) ^ (*n(x));

  return alpha;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
Real BLOCK_VECTOR::maxValue()
{
  Real maxFound = _blocks[0]->maxValue();
  for (int x = 1; x < _totalBlocks; x++)
  {
    Real currentMax = _blocks[x]->maxValue();
    maxFound = (currentMax > maxFound) ? currentMax : maxFound;
  }
  return maxFound;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int BLOCK_VECTOR::totalEntries() const
{
  int size = 0;
  for (int x = 0; x < _totalBlocks; x++)
    if (_blocks[x] != NULL)
      size += _blocks[x]->size();
  return size;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR BLOCK_VECTOR::blockSizes()
{
  VECTOR sizes(_totalBlocks);
  for (int x = 0; x < _totalBlocks; x++)
    sizes[x] = blockSize(x);
  return sizes;
}

//////////////////////////////////////////////////////////////////////
// populate a vector with a subset of the block vectors
//////////////////////////////////////////////////////////////////////
VECTOR BLOCK_VECTOR::full(int startBlock, int endBlock)
{
  assert(startBlock >= 0);
  assert(startBlock < endBlock);
  assert(startBlock <= _totalBlocks);
  assert(endBlock <= _totalBlocks);

  int totalEntries = 0;
  for (int x = startBlock; x < endBlock; x++)
  {
    assert(_blocks[x] != NULL);
    totalEntries += _blocks[x]->size();
  }

  VECTOR final(totalEntries);
  int i = 0;
  for (int x = startBlock; x < endBlock; x++)
  {
    for (int y = 0; y < _blocks[x]->size(); y++)
      final[i + y] = (*_blocks[x])[y];
    i += _blocks[x]->size();
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// clamp entries smaller than a threshold to zero
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::clampToZero(const Real threshold)
{
  for (int x = 0; x < _totalBlocks; x++)
    if (exists(x))
      _blocks[x]->clampToZero(threshold);
}

//////////////////////////////////////////////////////////////////////
// write block vector to a stream
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::write(FILE* file)
{
  // write dimensions
  fwrite((void*)&_totalBlocks, sizeof(int), 1, file);

  // cycle through the block
  int nullEntry = 0;
  int entryExists = 1;
  for (int x = 0; x < _totalBlocks; x++)
  {
    // write out a null flag
    if (_blocks[x] == NULL)
      fwrite((void*)&nullEntry, sizeof(int), 1, file);
    else
    {
      fwrite((void*)&entryExists, sizeof(int), 1, file);
      _blocks[x]->write(file);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// read block vector from a stream
//////////////////////////////////////////////////////////////////////
void BLOCK_VECTOR::read(FILE* file)
{
  // stomp any previously created block vector
  if (_blocks)
    for (int x = 0; x < _totalBlocks; x++)
      delete _blocks[x];
  delete[] _blocks;

  // read in dimensions
  fread((void*)&_totalBlocks, sizeof(int), 1, file);

  _blocks = new VECTOR*[_totalBlocks];
  for (int x = 0; x < _totalBlocks; x++)
    _blocks[x] = NULL;

  int existsFlag;
  for (int x = 0; x < _totalBlocks; x++)
  {
    fread((void*)&existsFlag, sizeof(int), 1, file);

    if (existsFlag == 1)
    {
      _blocks[x] = new VECTOR();
      _blocks[x]->read(file);
    }
  }
}
