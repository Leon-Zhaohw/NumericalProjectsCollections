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
// BIG_MATRIX.h: interface for the BIG_MATRIX class, which tries to
// get around the 2 GB limit on a malloc, and does some out-of-core
// things
//
//////////////////////////////////////////////////////////////////////

#ifndef BIG_MATRIX_H
#define BIG_MATRIX_H

#include "EIGEN.h"

#include <SETTINGS.h>
#include <map>
#include <iostream>
#include <cstdio>
#include <vector>
#include <VECTOR.h>
#include <MATRIX.h>

#include <memory.h>
#include <zlib.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class BIG_MATRIX 
{

public:
  BIG_MATRIX();
  BIG_MATRIX(int rows, int cols);
  BIG_MATRIX(const BIG_MATRIX& m);
  BIG_MATRIX(const MATRIX& m);
  BIG_MATRIX(const string& filename);

  virtual ~BIG_MATRIX() {};

  inline Real& operator()(int row, int col) {
    assert(row < _rows);
    assert(col < _cols);
    return _columns[col][row];
  };
 
  // const version
  inline Real operator()(int row, int col) const {
    assert(row < _rows);
    assert(col < _cols);
    return _columns[col][row];
  };

  // access a column directly
  inline VECTOR& operator[](int col) { return _columns[col]; };
  inline const VECTOR& operator[](int col) const { return _columns[col]; };

  int& rows() { return _rows; };
  int& cols() { return _cols; };
  int rows() const { return _rows; };
  int cols() const { return _cols; };
  static int& blockSize() { return _blockSize; };
  static string& scratchPath() { return _scratchPath; };

  // read in columns from a stream
  bool readColumns(const int rows, const int totalColumns, FILE* file);
  
  // write the matrix to a binary file
  void write(const string& filename);
  void read(const string& filename);

  // overload operators
  BIG_MATRIX& operator=(const BIG_MATRIX m);

  // do an out-of-core Gram-Schmidt QR factorization
  static MATRIX outOfCoreQR(const string& filenamePrefix, int& qRows, int& qCols);

  // do an out-of-core SVD using the ooc QR
  static void outOfCoreSVD(const string& filenamePrefix, const string& reducedPath, const Real& discardThreshold);

  // read/write a dimensions file to the scratch space
  static void writeDimensions(const string filename, const int& rows, const int& cols);
  static void readDimensions(const string filename, int& rows, int& cols);

  // compute a feasbile block size
  static int computeBlockSize(int rows, int cols);

  // write out the final big matrix
  static void writeFinalU(const string& filename);

protected:
  int _rows;
  int _cols;

  vector<VECTOR> _columns;

  // out-of-core variables
  static int _blockSize;
  static string _scratchPath;
};

ostream& operator<<(ostream &out, const BIG_MATRIX& matrix);

#endif
