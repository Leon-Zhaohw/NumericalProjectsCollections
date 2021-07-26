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
// SANDWICH_TRANSFORM.h: interface for the SANDWICH_TRANSFORM class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SANDWICH_TRANSFORM_H
#define SANDWICH_TRANSFORM_H

#include <SETTINGS.h>
#include <MATRIX.h>
#include <MATRIX3.h>
#include <TENSOR3.h>
#include <BLOCK_MATRIX.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Sandwich transform class -- assumes the middle matrix is 3x3
//////////////////////////////////////////////////////////////////////
class SANDWICH_TRANSFORM {

public:
  SANDWICH_TRANSFORM();

  // note that left matrix gets transposed!
  SANDWICH_TRANSFORM(MATRIX& leftMatrix, MATRIX& rightMatrix);
  SANDWICH_TRANSFORM(MATRIX& leftMatrix, MATRIX& rightMatrix, int middleRows, int middleCols);
  SANDWICH_TRANSFORM(MATRIX& leftMatrix, VECTOR& rightVector);
  void init(MATRIX& leftMatrix, MATRIX& rightMatrix, int middleRows, int middleCols);
  void init3x3(MATRIX& leftMatrix, MATRIX& rightMatrix);
  void init3x3(MATRIX& leftMatrix, VECTOR& rightVector);
  void init3x3(VECTOR& leftVector, MATRIX& rightMatrix);
  void init3x3(VECTOR& leftVector, VECTOR& rightVector);

  SANDWICH_TRANSFORM(const SANDWICH_TRANSFORM& sandwich);
  virtual ~SANDWICH_TRANSFORM();
  
  MATRIX transform(MATRIX3& middleMatrix);
  MATRIX transform(const MATRIX& middleMatrix);
  VECTOR vectorTransform(MATRIX3& middleMatrix);
  Real scalarTransform(MATRIX3& middleMatrix);
  MATRIX transform(TENSOR3& middleTensor);
  TENSOR3 tensorTransform(const TENSOR3& middleTensor);
  bool initialized() { return _initialized; };
  int rows() const { return _leftMatrix.rows(); };
  int cols() const { return _rightMatrix.cols(); };
  int middleRows() const { return _middleRows; };
  int middleCols() const { return _middleCols; };

  // in-place versions of above
  void transformInplace(const MATRIX& middleMatrix, MATRIX& final);
  void transformInplace(TENSOR3& middleTensor, MATRIX& final);
  void tensorTransformInplace(const TENSOR3& middleTensor, TENSOR3& final);
  void vectorTransformAxpy(MATRIX3& middleMatrix, VECTOR& final);
  void vectorTransformAxpy(const Real scalar, MATRIX3& middleMatrix, VECTOR& final);

  int sandwichRows() const { return _slices[0].rows(); };
  int sandwichCols() const { return _slices[0].cols(); };
  MATRIX* slices() { return _slices; };
  MATRIX slice(int x) { return _slices[x]; };
  const MATRIX& slice(int x) const { return _slices[x]; };

  // file IO
  void read(FILE* file);
  void write(FILE* file);
  
  SANDWICH_TRANSFORM& operator=(const SANDWICH_TRANSFORM& sandwich);

  SANDWICH_TRANSFORM& operator+=(const SANDWICH_TRANSFORM& sandwich);

  VECTOR dims() { return _slices[0].dims(); };

protected:
  MATRIX* _slices;
  bool _initialized;

  MATRIX fullTransform3x3(MATRIX& leftMatrix, MATRIX3& middleMatrix, MATRIX& rightMatrix);
  MATRIX fullTransform(MATRIX& leftMatrix, MATRIX& middleMatrix, MATRIX& rightMatrix);

  // for a non-3x3 matrix, the matrix dimensions
  int _middleRows;
  int _middleCols;

  // DEBUG
  MATRIX _leftMatrix;
  MATRIX _rightMatrix;

  MATRIX _vectorWorkspace;
};

SANDWICH_TRANSFORM operator+(const SANDWICH_TRANSFORM& A, const SANDWICH_TRANSFORM& B);
SANDWICH_TRANSFORM operator*(const Real& scalar, const SANDWICH_TRANSFORM& B);

#endif
