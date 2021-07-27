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

#ifndef SPARSE_MATRIX_ARRAY_H
#define SPARSE_MATRIX_ARRAY_H

#include "EIGEN.h"

#include <SETTINGS.h>
#include <MATRIX.h>
#include <vector>
#include <map>
#include <iostream>
#include <cstdio>
#include "SPARSE_MATRIX.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// A sparse matrix class based on maps
//////////////////////////////////////////////////////////////////////
class SPARSE_MATRIX_ARRAY {

public:
  SPARSE_MATRIX_ARRAY(int rows, int cols);

  virtual ~SPARSE_MATRIX_ARRAY() { if (_matrix) delete[] _matrix; };

  // get the reference to an entry -
  // note that if an entry doesn't exist, it will be created and
  // set to zero.
  // 
  // to check if an entry already exists, use the exists() function
  Real& operator()(int row, int col);
  const int rows() const { return _rows; };
  const int cols() const { return _cols; };

  const vector<pair<int, Real> >& row(int index) const { return _matrix[index]; };

  // return the dimensions packed into a vector (for output, usually)
  VECTOR dims() const { VECTOR final(2); final[0] = _rows; final[1] = _cols; return final; };

  // Set the matrix to zero. Note this will *NOT* stomp the underlying map!
  // It will instead set all current entries to zero so that we are not
  // forced to reallocate the sparsity structure again
  void clear();

  // project by the given left and right matrices, where we assume the left one needs
  // to be transposed
  MatrixXd projectVerySparse(const MatrixXd& left, const MatrixXd& right) const;

protected:
  int _rows;
  int _cols;

  // row, col, entry
  vector<pair<int, Real> >* _matrix;
};

#endif
