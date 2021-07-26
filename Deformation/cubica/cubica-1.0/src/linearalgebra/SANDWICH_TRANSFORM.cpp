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

#include "SANDWICH_TRANSFORM.h"

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
SANDWICH_TRANSFORM::SANDWICH_TRANSFORM(MATRIX& leftMatrix, MATRIX& rightMatrix):
  _slices(NULL),
  _initialized(false),
  _middleRows(3),
  _middleCols(3)
{
  init3x3(leftMatrix, rightMatrix);
}

SANDWICH_TRANSFORM::SANDWICH_TRANSFORM(MATRIX& leftMatrix, VECTOR& rightVector):
  _slices(NULL),
  _initialized(false),
  _middleRows(3),
  _middleCols(3)
{
  MATRIX rightMatrix(rightVector.size(), 1, rightVector.data());
  init3x3(leftMatrix, rightMatrix);
}

SANDWICH_TRANSFORM::SANDWICH_TRANSFORM(MATRIX& leftMatrix, MATRIX& rightMatrix, int middleRows, int middleCols):
  _slices(NULL),
  _initialized(false),
  _middleRows(middleRows),
  _middleCols(middleCols)
{
  init(leftMatrix, rightMatrix, _middleRows, _middleCols);
}

SANDWICH_TRANSFORM::SANDWICH_TRANSFORM():
  _slices(NULL),
  _initialized(false),
  _middleRows(0),
  _middleCols(0)
{
}

SANDWICH_TRANSFORM::SANDWICH_TRANSFORM(const SANDWICH_TRANSFORM& sandwich) :
  _initialized(sandwich._initialized),
  _middleRows(sandwich._middleRows),
  _middleCols(sandwich._middleCols)
{
  if (sandwich._slices == NULL)
  {
    _slices = NULL;
    return;
  }

  int entries = _middleRows * _middleCols;
  _slices = new MATRIX[entries];
  for (int x = 0; x < entries; x++)
    _slices[x] = sandwich._slices[x];
}

SANDWICH_TRANSFORM::~SANDWICH_TRANSFORM()
{
  delete[] _slices;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::init3x3(MATRIX& leftMatrix, VECTOR& rightVector)
{
  assert(leftMatrix.rows() == rightVector.size());

  MATRIX rightMatrix(rightVector.size(), 1, rightVector.data());
  init3x3(leftMatrix, rightMatrix);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::init3x3(VECTOR& leftVector, MATRIX& rightMatrix)
{
  assert(leftVector.size() == rightMatrix.rows());

  MATRIX leftMatrix(leftVector.size(), 1, leftVector.data());
  init3x3(leftMatrix, rightMatrix);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::init3x3(VECTOR& leftVector, VECTOR& rightVector)
{
  assert(leftVector.size() == rightVector.size());

  MATRIX rightMatrix(rightVector.size(), 1, rightVector.data());
  MATRIX leftMatrix(leftVector.size(), 1, leftVector.data());
  init3x3(leftMatrix, rightMatrix);
}

//////////////////////////////////////////////////////////////////////
// initialize the sandwich
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::init3x3(MATRIX& leftMatrix, MATRIX& rightMatrix)
{
  // DEBUG
  //_leftMatrix = leftMatrix;
  //_rightMatrix = rightMatrix;
  
  _initialized = true;

  int rows = rightMatrix.rows();

  // sandwich currently only works for 3x3 submatrices!
  assert(rows % 3 == 0);

  _middleRows = 3;
  _middleCols = 3;

  // if necessary, stomp previous results
  if (_slices)
    delete[] _slices;

  // allocate the slices
  _slices = new MATRIX[9];

  int i = 0;
  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++, i++)
    {
      MATRIX3 middleMatrix;
      middleMatrix(x,y) = 1.0;
      _slices[i] = fullTransform3x3(leftMatrix, middleMatrix, rightMatrix);
    }
}

//////////////////////////////////////////////////////////////////////
// initialize the arbitrarily sized sandwich
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::init(MATRIX& leftMatrix, MATRIX& rightMatrix, int middleRows, int middleCols)
{
  _middleRows = middleRows;
  _middleCols = middleCols;

  assert(_middleRows > 0);
  assert(_middleCols > 0);
  _initialized = true;

  _slices = new MATRIX[_middleCols * _middleRows];

  int i = 0;
  for (int y = 0; y < _middleCols; y++)
    for (int x = 0; x < _middleRows; x++, i++)
    {
      MATRIX middleMatrix(_middleRows, _middleCols);
      middleMatrix(x,y) = 1.0;
      _slices[i] = fullTransform(leftMatrix, middleMatrix, rightMatrix);
    }
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich the inefficient way -- this is needed to compute
// the slices, and is also useful for debugging
//////////////////////////////////////////////////////////////////////
MATRIX SANDWICH_TRANSFORM::fullTransform(MATRIX& leftMatrix, MATRIX& middleMatrix, MATRIX& rightMatrix)
{
  /*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " HAS NOT BEEN DEBUGGED! " << endl;
  MATRIX middleRight(rightMatrix);
  for (int y = 0; y < rightMatrix.cols(); y++)
    for (int x = 0; x < rightMatrix.rows() / _middleRows; x++)
    {
      VECTOR rightVec(_middleRows);
      for (int z = 0; z < _middleRows; z++)
        rightVec[z] = rightMatrix(_middleRows * x + z, y);

      VECTOR final = middleMatrix * rightVec;
      for (int z = 0; z < _middleRows; z++)
        middleRight(_middleRows * x + z, y) = final[z];
    }
    */
  int blockRows = rightMatrix.rows() / middleMatrix.cols();
  BLOCK_MATRIX blockedMiddleMatrix(middleMatrix, blockRows);
  MATRIX middleRight = blockedMiddleMatrix * rightMatrix;
 
  return leftMatrix ^ middleRight;
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich the inefficient way -- this is needed to compute
// the slices, and is also useful for debugging
//////////////////////////////////////////////////////////////////////
MATRIX SANDWICH_TRANSFORM::fullTransform3x3(MATRIX& leftMatrix, MATRIX3& middleMatrix, MATRIX& rightMatrix)
{
  MATRIX middleRight(rightMatrix);
  for (int y = 0; y < rightMatrix.cols(); y++)
    for (int x = 0; x < rightMatrix.rows() / 3; x++)
    {
      VEC3F rightVec;
      rightVec[0] = rightMatrix(3 * x, y);
      rightVec[1] = rightMatrix(3 * x + 1, y);
      rightVec[2] = rightMatrix(3 * x + 2, y);

      VEC3F final = middleMatrix * rightVec;
      middleRight(3 * x, y) = final[0];
      middleRight(3 * x + 1, y) = final[1];
      middleRight(3 * x + 2, y) = final[2];
    }
  
  return leftMatrix ^ middleRight;
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::transformInplace(const MATRIX& middleMatrix, MATRIX& final)
{
  assert(final.rows() == _slices[0].rows());
  assert(final.cols() == _slices[0].cols());

  int i = 0;
  final.clear();
  for (int y = 0; y < _middleCols; y++)
    for (int x = 0; x < _middleRows; x++, i++)
      final.axpy(middleMatrix(x,y), _slices[i]);
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
MATRIX SANDWICH_TRANSFORM::transform(const MATRIX& middleMatrix)
{
  // debug version
  //return fullTransform(_leftMatrix, middleMatrix, _rightMatrix);
  
  MATRIX final(_slices[0].rows(), _slices[0].cols());

  int i = 0;
  for (int y = 0; y < _middleCols; y++)
    for (int x = 0; x < _middleRows; x++, i++)
      //final += middleMatrix(x,y) * _slices[i];
      final.axpy(middleMatrix(x,y), _slices[i]);

  return final;
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
MATRIX SANDWICH_TRANSFORM::transform(TENSOR3& middleTensor)
{
  assert(_initialized);

  MATRIX final(_slices[0].rows(), middleTensor.slabs());

  /*
  for (int z = 0; z < middleTensor.slabs(); z++)
  {
    MATRIX slab = middleTensor.slab(z);
    int slice = 0;
    for (int y = 0; y < slab.cols(); y++)
      for (int x = 0; x < slab.rows(); x++, slice++)
      {
        MATRIX product = slab(x,y) * _slices[slice];

        assert(product.cols() == 1);

        final(0, z) += product(0,0);
        final(1, z) += product(1,0);
        final(2, z) += product(2,0);
      }
  }
  */
  for (int z = 0; z < middleTensor.slabs(); z++)
  {
    MATRIX& slab = middleTensor.slab(z);
    int slice = 0;
    for (int y = 0; y < slab.cols(); y++)
      for (int x = 0; x < slab.rows(); x++, slice++)
      {
        assert(_slices[slice].cols() == 1);
        final(0, z) += slab(x,y) * _slices[slice](0,0);
        final(1, z) += slab(x,y) * _slices[slice](1,0);
        final(2, z) += slab(x,y) * _slices[slice](2,0);
      }
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::transformInplace(TENSOR3& middleTensor, MATRIX& final)
{
  assert(_initialized);

  final.clear();
  for (int z = 0; z < middleTensor.slabs(); z++)
  {
    MATRIX& slab = middleTensor.slab(z);
    int slice = 0;
    for (int y = 0; y < slab.cols(); y++)
      for (int x = 0; x < slab.rows(); x++, slice++)
      {
        final(0, z) += slab(x,y) * _slices[slice](0,0);
        final(1, z) += slab(x,y) * _slices[slice](1,0);
        final(2, z) += slab(x,y) * _slices[slice](2,0);
      }
  }
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
TENSOR3 SANDWICH_TRANSFORM::tensorTransform(const TENSOR3& middleTensor)
{
  TENSOR3 final(_slices[0].rows(), _slices[0].cols(), middleTensor.slabs());

  for (int z = 0; z < middleTensor.slabs(); z++)
    final.slab(z) = transform(middleTensor.constSlab(z));
  return final;
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::tensorTransformInplace(const TENSOR3& middleTensor, TENSOR3& final)
{
  final.clear();
  for (int z = 0; z < middleTensor.slabs(); z++)
    transformInplace(middleTensor.constSlab(z), final.slab(z));
}

/*
//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
MATRIX SANDWICH_TRANSFORM::transform(MATRIX3& middleMatrix)
{
  // debug version
  //return fullTransform(_leftMatrix, middleMatrix, _rightMatrix);
  
  MATRIX final(_slices[0].rows(), _slices[0].cols());

  final += middleMatrix(0,0) * _slices[0];
  final += middleMatrix(1,0) * _slices[1];
  final += middleMatrix(2,0) * _slices[2];
  final += middleMatrix(0,1) * _slices[3];
  final += middleMatrix(1,1) * _slices[4];
  final += middleMatrix(2,1) * _slices[5];
  final += middleMatrix(0,2) * _slices[6];
  final += middleMatrix(1,2) * _slices[7];
  final += middleMatrix(2,2) * _slices[8];

  return final;
}
*/

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
Real SANDWICH_TRANSFORM::scalarTransform(MATRIX3& middleMatrix)
{
  Real final = 
    middleMatrix(0,0) * _slices[0](0,0) +
    middleMatrix(1,0) * _slices[1](0,0) +
    middleMatrix(2,0) * _slices[2](0,0) +
    middleMatrix(0,1) * _slices[3](0,0) +
    middleMatrix(1,1) * _slices[4](0,0) +
    middleMatrix(2,1) * _slices[5](0,0) +
    middleMatrix(0,2) * _slices[6](0,0) +
    middleMatrix(1,2) * _slices[7](0,0) +
    middleMatrix(2,2) * _slices[8](0,0);
  return final;
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::vectorTransformAxpy(MATRIX3& middleMatrix, VECTOR& final)
{
  assert(_initialized == true);

  // should only wipe if the size is already correct
  _vectorWorkspace.resizeAndWipe(final.size(), 1); 

  transformInplace(middleMatrix, _vectorWorkspace);

  for (int x = 0; x < final.size(); x++)
    final[x] += _vectorWorkspace(x, 0);
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::vectorTransformAxpy(const Real scalar, MATRIX3& middleMatrix, VECTOR& final)
{
  assert(_initialized == true);

  // should only wipe if the size is already correct
  _vectorWorkspace.resizeAndWipe(final.size(), 1); 

  transformInplace(middleMatrix, _vectorWorkspace);

  for (int x = 0; x < final.size(); x++)
    final[x] += scalar * _vectorWorkspace(x, 0);
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
VECTOR SANDWICH_TRANSFORM::vectorTransform(MATRIX3& middleMatrix)
{
  assert(_initialized == true);

  MATRIX final = transform(middleMatrix);
  assert(final.cols() == 1 || final.rows() == 1);

  if (final.cols() == 1)
  {
    VECTOR finalVector(final.rows());
    for (int x = 0; x < final.rows(); x++)
      finalVector[x] = final(x,0);

    return finalVector;
  }

  VECTOR finalVector(final.cols());
  for (int x = 0; x < final.cols(); x++)
    finalVector[x] = final(0,x);

  return finalVector;
}

//////////////////////////////////////////////////////////////////////
// read from a file stream
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::read(FILE* file)
{
  fread((void*)&_middleRows, sizeof(int), 1, file);
  fread((void*)&_middleCols, sizeof(int), 1, file);
  fread((void*)&_initialized, sizeof(bool), 1, file);

  if (_slices) delete[] _slices;

  _slices = new MATRIX[_middleRows * _middleCols];
  for (int x = 0; x < _middleRows * _middleCols; x++)
    _slices[x].read(file);
}

//////////////////////////////////////////////////////////////////////
// write to a file stream
//////////////////////////////////////////////////////////////////////
void SANDWICH_TRANSFORM::write(FILE* file)
{
  assert(_middleRows > 0);
  assert(_middleCols > 0);

  fwrite((void*)&_middleRows, sizeof(int), 1, file);
  fwrite((void*)&_middleCols, sizeof(int), 1, file);
  fwrite((void*)&_initialized, sizeof(bool), 1, file);
  for (int x = 0; x < _middleRows * _middleCols; x++)
    _slices[x].write(file);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
SANDWICH_TRANSFORM& SANDWICH_TRANSFORM::operator=(const SANDWICH_TRANSFORM& sandwich)
{
  if (_slices)
    delete[] _slices;

  _initialized = sandwich._initialized;
  _middleRows = sandwich._middleRows;
  _middleCols = sandwich._middleCols;

  if (sandwich._slices == NULL)
    return *this;

  int entries = _middleRows * _middleCols;
  _slices = new MATRIX[entries];
  for (int x = 0; x < entries; x++)
    _slices[x] = sandwich._slices[x];

  return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
SANDWICH_TRANSFORM& SANDWICH_TRANSFORM::operator+=(const SANDWICH_TRANSFORM& sandwich)
{
  assert(_middleRows == sandwich.middleRows());
  assert(_middleCols == sandwich.middleCols());

  int totalSlices = _middleRows * _middleCols;
  for (int x = 0; x < totalSlices; x++)
    _slices[x] += sandwich.slice(x);

  return *this;
}

//////////////////////////////////////////////////////////////////////
// Apply the sandwich
//////////////////////////////////////////////////////////////////////
MATRIX SANDWICH_TRANSFORM::transform(MATRIX3& middleMatrix)
{
  assert(_middleCols != 0);
  assert(_middleRows != 0);
  assert(_middleCols == 3);
  assert(_middleRows == 3);

  MATRIX final(_slices[0].rows(), _slices[0].cols());

  int i = 0;
  for (int y = 0; y < _middleCols; y++)
    for (int x = 0; x < _middleRows; x++, i++)
      final.axpy(middleMatrix(x,y), _slices[i]);

  return final;
}
/*
{
  // debug version
  //return fullTransform(_leftMatrix, middleMatrix, _rightMatrix);
  
  MATRIX final(_slices[0].rows(), _slices[0].cols());

  final.axpy(middleMatrix(0,0), _slices[0]);
  final.axpy(middleMatrix(1,0), _slices[1]);
  final.axpy(middleMatrix(2,0), _slices[2]);
  final.axpy(middleMatrix(0,1), _slices[3]);
  final.axpy(middleMatrix(1,1), _slices[4]);
  final.axpy(middleMatrix(2,1), _slices[5]);
  final.axpy(middleMatrix(0,2), _slices[6]);
  final.axpy(middleMatrix(1,2), _slices[7]);
  final.axpy(middleMatrix(2,2), _slices[8]);

  return final;
}
*/

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
SANDWICH_TRANSFORM operator+(const SANDWICH_TRANSFORM& A, const SANDWICH_TRANSFORM& B)
{
  assert(A.middleRows() == B.middleRows());
  assert(A.middleCols() == B.middleCols());
  SANDWICH_TRANSFORM final(A);

  int totalSlices = A.middleRows() * A.middleCols();
  for (int x = 0; x < totalSlices; x++)
    final.slice(x) += B.slice(x);
  return final; 
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
SANDWICH_TRANSFORM operator*(const Real& scalar, const SANDWICH_TRANSFORM& A)
{
  SANDWICH_TRANSFORM final(A);

  int totalSlices = A.middleRows() * A.middleCols();
  for (int x = 0; x < totalSlices; x++)
    final.slices()[x] *= scalar;
  return final; 
}
