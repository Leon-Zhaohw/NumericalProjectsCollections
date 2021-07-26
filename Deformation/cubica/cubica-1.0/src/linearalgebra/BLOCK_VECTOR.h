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

#ifndef BLOCK_VECTOR_H
#define BLOCK_VECTOR_H

#include <MATRIX.h>
#include <VEC3.h>

//////////////////////////////////////////////////////////////////////
// Block vector -- essentially a 1D storage class for a bunch of 
// vectors
//
// This was created so that the logic of checking all the indices
// when multiplying with a BLOCK_MATRIX would only have to be written
// once
//
// There is no dimension checking, so you better know what you're
// doing.
//////////////////////////////////////////////////////////////////////
class BLOCK_VECTOR {

public:
  BLOCK_VECTOR();
  BLOCK_VECTOR(int blocks);
  BLOCK_VECTOR(const vector<int>& blockSizes);
  BLOCK_VECTOR(const vector<VECTOR>& blocks);
  BLOCK_VECTOR(const BLOCK_VECTOR& vec);
  virtual ~BLOCK_VECTOR();

  // set dimensions
  void resizeAndWipe(int blocks);
  void resizeAndWipeBlock(int block, int size);
  void resize(BLOCK_VECTOR& vec);
  void resizeAndWipe(BLOCK_VECTOR& vec);
  
  // check if a block entry exists
  bool exists(int block) const { return (_blocks[block] != NULL); };

  // get the matrix at this block entry
  // if the matrix does not exist, it will return NULL
  VECTOR* entry(int block) { return _blocks[block]; };
  VECTOR& entry(int block) const { return *_blocks[block]; };
  const VECTOR* constEntry(int block) const { return _blocks[block]; };

  // naked accessor -- doesn't do any checks
  inline VECTOR* operator()(int block) {
    assert(block >= 0);
    assert(block < _totalBlocks);
    return _blocks[block];
  };
  
  inline VECTOR& operator[] (int block) {
    assert(block >= 0);
    assert(block < _totalBlocks);
    return *_blocks[block];
  };

  // add this vector to a block
  // 
  // If this is the first time this block is called, the vector
  // will be set to the size of the vector you pass in.
  // 
  // It is your responsibility to ensure that the vector dimensions
  // match for subsequent calls.
  void add(const VECTOR& vec, int whichBlock);

  // subtract this vector from a block
  // 
  // If this is the first time this block is called, the vector
  // will be set to the size of the vector you pass in.
  // 
  // It is your responsibility to ensure that the vector dimensions
  // match for subsequent calls.
  void subtract(const VECTOR& vec, int whichBlock);

  // set this block to the current vector
  void set(const VECTOR& vec, int whichBlock);
  void set(const VEC3& vec, int whichBlock);
  void equals(const VECTOR& vec, int whichBlock) { set(vec, whichBlock); };

  // 2 norm of this vector
  Real norm2();

  // max norm
  Real maxValue();
 
  // wipe all the entries
  void clear();
 
  // sum squared of this vector
  Real sum2();

  // accessors
  int totalBlocks() const { return _totalBlocks; };
  VECTOR** data() { return _blocks; };
  int totalEntries() const;
  int blockSize(int x) const { if (_blocks[x] == NULL) return -1; return _blocks[x]->size(); };
  VECTOR blockSizes();

  // overload operators
  BLOCK_VECTOR& operator+=(const BLOCK_VECTOR& m);
  BLOCK_VECTOR& operator+=(const VECTOR& m);
  BLOCK_VECTOR& operator-=(const VECTOR& m);
  BLOCK_VECTOR& operator*=(const Real scalar);
  BLOCK_VECTOR& operator-=(const BLOCK_VECTOR& m);
  BLOCK_VECTOR& operator=(const BLOCK_VECTOR& m);
  BLOCK_VECTOR& operator=(VECTOR& m);

  // dump to a full vector
  VECTOR full();
  VECTOR full(int startBlock, int endBlock);

  void axpy(Real scalar, BLOCK_VECTOR& m);
  void clearingAxpy(Real scalar, BLOCK_VECTOR& m);

  // clamp entries smaller than a threshold to zero
  void clampToZero(const Real threshold);

  // file IO
  void write(FILE* file);
  void read(FILE* file);

private:
  int _totalBlocks;
  VECTOR** _blocks;
};

// block matrix vector multiply -- assumes A points to a 1D array
// of matrices that correspond to the diagonal
//
// This is pretty dangerous, so use sparingly
BLOCK_VECTOR operator*(MATRIX** Adiagonal, BLOCK_VECTOR& x);

BLOCK_VECTOR operator*(const Real& scalar, const BLOCK_VECTOR& x);
BLOCK_VECTOR operator-(const BLOCK_VECTOR& x, const BLOCK_VECTOR& y);
BLOCK_VECTOR operator+(const BLOCK_VECTOR& x, const BLOCK_VECTOR& y);
Real operator^(BLOCK_VECTOR& m, BLOCK_VECTOR& n);
ostream& operator<<(ostream &out, BLOCK_VECTOR& vec);


#endif
