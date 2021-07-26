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
// TENSOR3.h: interface for the TENSOR3 class.
//
//////////////////////////////////////////////////////////////////////

#include "TENSOR3.h"

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
TENSOR3::TENSOR3() :
  _rows(0), _cols(0), _slabs(0)
{
  _tensor = NULL;
}

TENSOR3::TENSOR3(int rows, int cols, int slabs) :
  _rows(rows), _cols(cols), _slabs(slabs)
{
  _tensor = new MATRIX[slabs];

  for (int x = 0; x < slabs; x++)
    _tensor[x].resizeAndWipe(rows, cols);

  clear();
}

TENSOR3::TENSOR3(const TENSOR3& tensor) :
  _rows(tensor.rows()), _cols(tensor.cols()), _slabs(tensor.slabs())
{
  _tensor = new MATRIX[_slabs];

  for (int x = 0; x < _slabs; x++)
    _tensor[x] = tensor._tensor[x];
}


TENSOR3::TENSOR3(const TENSOR3& tensor, const int repeats) :
  _rows(tensor.rows() * repeats), 
  _cols(tensor.cols() * repeats), 
  _slabs(tensor.slabs())
{
  _tensor = new MATRIX[_slabs];

  for (int x = 0; x < _slabs; x++)
    _tensor[x].resizeAndWipe(_rows, _cols);

  for (int r = 0; r < repeats; r++)
    for (int z = 0; z < tensor.slabs(); z++)
      for (int y = 0; y < tensor.cols(); y++)
        for (int x = 0; x < tensor.rows(); x++)
        {
          int startRow = r * tensor.rows();
          int startCol = r * tensor.cols();
          _tensor[z](startRow + x, startCol + y) = tensor(x,y,z);
        }
}

TENSOR3::TENSOR3(const vector<MATRIX>& slabs)
{
  assert(slabs.size() > 0);

  _rows = slabs[0].rows();
  _cols = slabs[0].cols();
  _slabs = slabs.size();

  _tensor = new MATRIX[slabs.size()];

  for (int x = 0; x < _slabs; x++)
  {
    assert(slabs[x].rows() == _rows);
    assert(slabs[x].cols() == _cols);
    _tensor[x] = slabs[x];
  }
}

// the vector exponential derivative of a matrix
TENSOR3::TENSOR3(const MATRIX& A, const vector<MATRIX>& dAs) :
  _rows(A.rows()), _cols(A.cols()), _slabs(dAs.size())
{
  _tensor = new MATRIX[_slabs];
  for (int x = 0; x < _slabs; x++)
    _tensor[x].resizeAndWipe(_rows, _cols);

  // set to exponential derivative
  for (int x = 0; x < _slabs; x++)
    _tensor[x] = A.dexp(dAs[x]);
}

TENSOR3::~TENSOR3()
{
  delete[] _tensor;
}

void TENSOR3::clear()
{
  for (int x = 0; x < _slabs; x++)
    _tensor[x].clear();
}

void TENSOR3::resizeAndWipe(int rows, int cols, int slabs)
{
  delete[] _tensor;
  _tensor = new MATRIX[slabs];

  for (int x = 0; x < slabs; x++)
    _tensor[x].resizeAndWipe(rows, cols);

  _rows = rows;
  _cols = cols;
  _slabs = slabs;

  clear();
}

// These operators are ambiguous -- better to just use 
// modeOneProduct and modeThreeProduct
#if 0
//////////////////////////////////////////////////////////////////////
// multiply by vector by tensor -- directly multiply each
// slab by the vector and sum
//////////////////////////////////////////////////////////////////////
MATRIX operator*(TENSOR3& A, VECTOR& x)
{
  assert(A.slabs() == x.size());
  assert(A.slabs() > 0);

  int rows = A.rows();
  int cols = A.cols();
  int slabs = A.slabs();

  MATRIX result(rows, cols);
  result.clearingAxpy(x[0], A.data()[0]);
  for (int slab = 1; slab < slabs; slab++)
    result.axpy(x[slab], A.data()[slab]);

  /*
  MATRIX result(rows, slabs);
  for (int i = 0; i < slabs; i++)
  {
    VECTOR column = A.data()[i] * x;
    for (int j = 0; j < rows; j++)
      result(j, i) = column[j];
  }
  */

  return result;
}

//////////////////////////////////////////////////////////////////////
// Multiply in the other direction -- each slab is multiplied by
// the vector, and the resulting vectors are then packed into the
// columns for a new matrix
//////////////////////////////////////////////////////////////////////
MATRIX operator^(TENSOR3& A, VECTOR& x)
{
  assert(A.cols() == x.size());
  assert(A.slabs() > 0);

  int rows = A.rows();
  int cols = A.cols();
  int slabs = A.slabs();

  MATRIX result(rows, cols);
  for (int slab = 0; slab < slabs; slab++)
  {
    VECTOR column = A.data()[slab] * x;
    result.setColumn(column, slab);
  }
  return result;
}
#endif

//////////////////////////////////////////////////////////////////////
// take the product w.r.t. mode one, i.e. matrix-vector multiply 
// with each slab
//////////////////////////////////////////////////////////////////////
MATRIX TENSOR3::modeOneProduct(const VECTOR& x) const
{
  assert(_cols == x.size());
  assert(_slabs > 0);

  MATRIX result(_rows, _slabs);
  for (int slab = 0; slab < _slabs; slab++)
  {
    VECTOR column = _tensor[slab] * x;
    result.setColumn(column, slab);
  }
  return result;
}

//////////////////////////////////////////////////////////////////////
// take the product w.r.t. mode one, i.e. matrix-vector multiply 
// with each slab
//////////////////////////////////////////////////////////////////////
void TENSOR3::modeOneProductInplace(const VECTOR& x, MATRIX& final)
{
  assert(_cols == x.size());
  assert(_slabs > 0);
  assert(final.cols() == _slabs);
  assert(final.rows() == _rows);

  /*
  VECTOR column(_tensor[0].rows());
  for (int slab = 0; slab < _slabs; slab++)
  {
    column.clear();
    _tensor[slab].gemvInplace(1.0, x, column);
    final.setColumn(column, slab);
  }
  */
  final.clear();
  if (_workspaceColumn.size() == 0)
    _workspaceColumn.resizeAndWipe(_tensor[0].rows());
  for (int slab = 0; slab < _slabs; slab++)
  {
    _workspaceColumn.clear();
    _tensor[slab].gemvInplace(1.0, x, _workspaceColumn);
    final.addColumn(_workspaceColumn, slab);
  }
}

//////////////////////////////////////////////////////////////////////
// take the product w.r.t. mode one, i.e. matrix-vector multiply 
// with each slab
//////////////////////////////////////////////////////////////////////
void TENSOR3::modeOneAxpy(const VECTOR& x, const Real scalar, MATRIX& final)
{
  assert(_cols == x.size());
  assert(_slabs > 0);
  assert(final.cols() == _slabs);
  assert(final.rows() == _rows);

  /*
  VECTOR column(_tensor[0].rows());
  for (int slab = 0; slab < _slabs; slab++)
  {
    VECTOR column = _tensor[slab] * x;
    column.clear();
    _tensor[slab].gemvInplace(scalar, x, column);
    final.addColumn(column, slab);
  }
  */

  if (_workspaceColumn.size() == 0)
    _workspaceColumn.resizeAndWipe(_tensor[0].rows());
  for (int slab = 0; slab < _slabs; slab++)
  {
    _workspaceColumn.clear();
    _tensor[slab].gemvInplace(scalar, x, _workspaceColumn);
    final.addColumn(_workspaceColumn, slab);
  }
}

//////////////////////////////////////////////////////////////////////
// take the product w.r.t. mode three, i.e. scale each slab by 
// vector entries, and sum them
//////////////////////////////////////////////////////////////////////
MATRIX TENSOR3::modeThreeProduct(const VECTOR& x)
{
  assert(_slabs == x.size());
  assert(_slabs > 0);

  MATRIX result(_rows, _cols);
  /*
  result.clearingAxpy(x[0], _tensor[0]);
  for (int slab = 1; slab < _slabs; slab++)
    result.axpy(x[slab], _tensor[slab]);
    */
  for (int slab = 0; slab < _slabs; slab++)
    result.axpy(x[slab], _tensor[slab]);

  return result;
}

//////////////////////////////////////////////////////////////////////
// take the product w.r.t. mode three, i.e. scale each slab by 
// vector entries, and sum them
//////////////////////////////////////////////////////////////////////
void TENSOR3::modeThreeAxpy(const VECTOR& x, const Real scalar, MATRIX& final)
{
  assert(_slabs == x.size());
  assert(_slabs > 0);

  /*
  final.clearingAxpy(scalar * x[0], _tensor[0]);
  //final.clearingAxpy(x[0], _tensor[0]);
  for (int slab = 1; slab < _slabs; slab++)
    final.axpy(scalar * x[slab], _tensor[slab]);
    //final.axpy(x[slab], _tensor[slab]);
    */
  for (int slab = 0; slab < _slabs; slab++)
    final.axpy(scalar * x[slab], _tensor[slab]);
}

//////////////////////////////////////////////////////////////////////
// take the product w.r.t. mode three, i.e. scale each slab by 
// vector entries, and sum them
//////////////////////////////////////////////////////////////////////
void TENSOR3::modeThreeProductInplace(const VECTOR& x, MATRIX& final)
{
  assert(_slabs == x.size());
  assert(_slabs > 0);

  assert(final.rows() == _rows);
  assert(final.cols() == _cols);

  final.clear();
  for (int slab = 0; slab < _slabs; slab++)
    final.axpy(x[slab], _tensor[slab]);
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, const TENSOR3& tensor)
{
  for (int x = 0; x < tensor.slabs(); x++)
  {
    out << " Slab " << x << ":";
    out << tensor.dataConst()[x];
  }
  return out;
}

//////////////////////////////////////////////////////////////////////
// Sum of all the entries
//////////////////////////////////////////////////////////////////////
Real TENSOR3::sum2()
{
  Real final = 0.0;
  for (int x = 0; x < _slabs; x++)
    final += _tensor[x].sum2();

  return final;
}

//////////////////////////////////////////////////////////////////////
// pointwise add of two tensors
//////////////////////////////////////////////////////////////////////
TENSOR3 operator+(const TENSOR3& A, const TENSOR3& B)
{
  assert(A.rows() == B.rows());
  assert(A.cols() == B.cols());
  assert(A.slabs() == B.slabs());

  TENSOR3 result(A.rows(), A.cols(), A.slabs());
  MATRIX* dataResult = result.data();
  const MATRIX* dataA = A.dataConst();
  const MATRIX* dataB = B.dataConst();
  for (int x = 0; x < A.slabs(); x++)
    dataResult[x] = dataA[x] + dataB[x];

  return result;
}

//////////////////////////////////////////////////////////////////////
// pointwise subtract of two tensors
//////////////////////////////////////////////////////////////////////
TENSOR3 operator-(const TENSOR3& A, const TENSOR3& B)
{
  assert(A.rows() == B.rows());
  assert(A.cols() == B.cols());
  assert(A.slabs() == B.slabs());

  TENSOR3 result(A.rows(), A.cols(), A.slabs());
  MATRIX* dataResult = result.data();
  const MATRIX* dataA = A.dataConst();
  const MATRIX* dataB = B.dataConst();
  for (int x = 0; x < A.slabs(); x++)
    dataResult[x] = dataA[x] - dataB[x];

  return result;
}

//////////////////////////////////////////////////////////////////////
// multiply each slab in A by B
//////////////////////////////////////////////////////////////////////
TENSOR3 operator*(const TENSOR3& A, const MATRIX& B)
{
  assert(A.cols() == B.rows());

  int rows = A.rows();
  int cols = B.cols();
  int slabs = A.slabs();

  TENSOR3 result(rows, cols, slabs);
  MATRIX* dataResult = result.data();
  const MATRIX* dataA = A.dataConst();
  for (int x = 0; x < slabs; x++)
    dataResult[x] = dataA[x] * B;

  return result;
}

//////////////////////////////////////////////////////////////////////
// multiply each slab in A by B
//////////////////////////////////////////////////////////////////////
void TENSOR3::multiplyInplace(const TENSOR3& A, const MATRIX& B)
{
  assert(A.cols() == B.rows());

  int rows = A.rows();
  int cols = B.cols();
  int slabs = A.slabs();

  assert(rows == _rows);
  assert(cols = _cols);
  assert(slabs = _slabs);

  clear();
  const MATRIX* dataA = A.dataConst();
  for (int x = 0; x < slabs; x++)
    //_tensor[x] = dataA[x] * B;
    _tensor[x].gemm(dataA[x], B);
}

//////////////////////////////////////////////////////////////////////
// multiply A by B, but cutting slabs in the other direction
//////////////////////////////////////////////////////////////////////
TENSOR3 operator^(const TENSOR3& A, const MATRIX& B)
{
  assert(A.slabs() == B.rows());

  int rows = A.rows();
  int cols = A.cols();
  int slabs = B.cols();

  TENSOR3 result(rows, cols, slabs);

  MATRIX slice(A.rows(), A.slabs());
  for (int x = 0; x < cols; x++)
  {
    // copy entries into slice
    for (int y = 0; y < A.rows(); y++)
      for (int z = 0; z < A.slabs(); z++)
        slice(y,z) = A(y, x, z);

    // do the multiply
    MATRIX product = slice * B;

    // copy the results into the final result
    for (int y = 0; y < rows; y++)
      for (int z = 0; z < slabs; z++)
        result(y, x, z) = product(y,z);
  }
  return result;
}

/*
//////////////////////////////////////////////////////////////////////
// mode 1 product (ie matrix-matrix multiply with each slab)
//////////////////////////////////////////////////////////////////////
TENSOR3 operator*(TENSOR3& A, MATRIX& B)
{
  assert(A.cols() == B.rows());

  int rows = A.rows();
  int cols = B.cols();
  int slabs = A.slabs();

  TENSOR3 result(rows, cols, slabs);

  for (int x = 0; x < slabs; x++)
    result.slab(x) = A.slab(x) * B;

  return result;
}
*/

//////////////////////////////////////////////////////////////////////
// transpose of the tensor -- just transpose each of the matrices
//////////////////////////////////////////////////////////////////////
TENSOR3 TENSOR3::transpose() const
{
  /*
  vector<MATRIX> slabs;

  for (int x = 0; x < _slabs; x++)
    slabs.push_back(_tensor[x].transpose());
    */
  TENSOR3 final(_cols, _rows, _slabs);
  for (int x = 0; x < _slabs; x++)
    _tensor[x].transpose(final._tensor[x]);

  return final;
}

//////////////////////////////////////////////////////////////////////
// flip the slabs and columns of the tensor
//////////////////////////////////////////////////////////////////////
TENSOR3 TENSOR3::flipSlabsAndColumns()
{
  TENSOR3 final(_rows, _slabs, _cols);
  for (int x = 0; x < _rows; x++)
    for (int y = 0; y < _cols; y++)
      for (int z = 0; z < _slabs; z++)
        final(x,z,y) = (*this)(x,y,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// flip the slabs and rows of the tensor
//////////////////////////////////////////////////////////////////////
TENSOR3 TENSOR3::flipSlabsAndRows()
{
  TENSOR3 final(_slabs, _cols, _rows);
  for (int x = 0; x < _rows; x++)
    for (int y = 0; y < _cols; y++)
      for (int z = 0; z < _slabs; z++)
        final(z,y,x) = (*this)(x,y,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// set the mode 3 entries according to the data
// this isn't that safe, but it avoids having to add a new function to MATRIX
//////////////////////////////////////////////////////////////////////
void TENSOR3::set(int row, int col, const Real* data)
{
  for (int z = 0; z < _slabs; z++)
    (*this)(row, col, z) = data[z];
}

//////////////////////////////////////////////////////////////////////
// scale along mode 3
//////////////////////////////////////////////////////////////////////
void TENSOR3::scale(int row, int col, Real alpha)
{
  for (int z = 0; z < _slabs; z++)
    (*this)(row, col, z) *= alpha;
}

//////////////////////////////////////////////////////////////////////
// self plus
//////////////////////////////////////////////////////////////////////
TENSOR3& TENSOR3::operator+=(const TENSOR3& m)
{
  assert(_rows == m.rows());
  assert(_cols == m.cols());
  assert(_slabs == m.slabs());

  for (int z = 0; z < _slabs; z++)
    _tensor[z] += m._tensor[z];

  return *this;
}

//////////////////////////////////////////////////////////////////////
// self minus
//////////////////////////////////////////////////////////////////////
TENSOR3& TENSOR3::operator-=(const TENSOR3& m)
{
  assert(_rows == m.rows());
  assert(_cols == m.cols());
  assert(_slabs == m.slabs());

  for (int z = 0; z < _slabs; z++)
    _tensor[z] -= m._tensor[z];

  return *this;
}

//////////////////////////////////////////////////////////////////////
// self scale
//////////////////////////////////////////////////////////////////////
TENSOR3& TENSOR3::operator*=(const Real& scalar)
{
  for (int z = 0; z < _slabs; z++)
    _tensor[z] *= scalar;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// compute the cross product tensor
//////////////////////////////////////////////////////////////////////
TENSOR3 TENSOR3::cross(const MATRIX& Uij)
{
  assert(Uij.rows() == 3);

  TENSOR3 final(3,3,Uij.cols());
  final.set(0,1, Uij.row(2));
  final.set(0,2, Uij.row(1));

  final.set(1,0, Uij.row(2));
  final.set(1,2, Uij.row(0));

  final.set(2,0, Uij.row(1));
  final.set(2,1, Uij.row(0));

  final.scale(0,1, -1.0);
  final.scale(1,2, -1.0);
  final.scale(2,0, -1.0);

  return final;
}

//////////////////////////////////////////////////////////////////////
// compute the cross product tensor
//////////////////////////////////////////////////////////////////////
TENSOR3 TENSOR3::crossTranspose(const MATRIX& Uij)
{
  assert(Uij.rows() == 3);

  TENSOR3 final(3,3,Uij.cols());
  final.set(1,0, Uij.row(2));
  final.set(2,0, Uij.row(1));

  final.set(0,1, Uij.row(2));
  final.set(2,1, Uij.row(0));

  final.set(0,2, Uij.row(1));
  final.set(1,2, Uij.row(0));

  final.scale(1,0, -1.0);
  final.scale(2,1, -1.0);
  final.scale(0,2, -1.0);

  return final;
}

//////////////////////////////////////////////////////////////////////
// mode 1 product with a block matrix
//////////////////////////////////////////////////////////////////////
TENSOR3 operator*(const BLOCK_MATRIX& A, const TENSOR3& B)
{
  TENSOR3 final(A.totalRows(), B.cols(), B.slabs());

  for (int x = 0; x < B.slabs(); x++)
    final.slab(x) = A * B.constSlab(x);

  return final;
}

//////////////////////////////////////////////////////////////////////
// mode 1 transpose product with a block matrix
//////////////////////////////////////////////////////////////////////
TENSOR3 operator^(const TENSOR3& A, const BLOCK_MATRIX& B)
{
  TENSOR3 final(A.cols(), B.totalCols(), A.slabs());

  for (int x = 0; x < A.slabs(); x++)
    final.slab(x) = A.constSlab(x).transpose() * B;

  return final;
}

//////////////////////////////////////////////////////////////////////
// mode 1 transpose product with a matrix
//////////////////////////////////////////////////////////////////////
TENSOR3 operator^(const MATRIX& A, const TENSOR3& B)
{
  TENSOR3 final(A.cols(), B.cols(), B.slabs());

  for (int x = 0; x < B.slabs(); x++)
    final.slab(x) = A ^ B.constSlab(x);

  return final;
}

//////////////////////////////////////////////////////////////////////
// mode 1 product with a matrix
//////////////////////////////////////////////////////////////////////
TENSOR3 operator*(const MATRIX& A, const TENSOR3& B)
{
  TENSOR3 final(A.rows(), B.cols(), B.slabs());

  for (int x = 0; x < B.slabs(); x++)
    final.slab(x) = A * B.constSlab(x);

  return final;
}

//////////////////////////////////////////////////////////////////////
// scale the entire tensor
//////////////////////////////////////////////////////////////////////
TENSOR3 operator*(const Real& a, const TENSOR3& B)
{
  TENSOR3 final(B.rows(), B.cols(), B.slabs());

  for (int x = 0; x < B.slabs(); x++)
    final.slab(x) = a * B.constSlab(x);

  return final;
}

//////////////////////////////////////////////////////////////////////
// copy the tensor
//////////////////////////////////////////////////////////////////////
TENSOR3& TENSOR3::operator=(const TENSOR3& m)
{
  if (_rows != m.rows() || _cols != m.cols() || _slabs != m.slabs())
  {
    _rows = m.rows();
    _cols = m.cols();
    _slabs = m.slabs();

    if (_tensor) delete[] _tensor;
    _tensor = new MATRIX[_slabs];
  }

  for (int x = 0; x < _slabs; x++)
    _tensor[x] = m._tensor[x];

  return *this;
}

//////////////////////////////////////////////////////////////////////
// stomp the current tensor with the given tensor
// starting at row number "row". It is your responsibility
// to ensure that you don't fall off the end of this tensor.
//////////////////////////////////////////////////////////////////////
void TENSOR3::setSubtensor(TENSOR3& tensor, int row)
{
  assert(tensor.slabs() == _slabs);
  assert(tensor.cols() == _cols);

  for (int x = 0; x < _slabs; x++)
  {
    MATRIX& slab = tensor.slab(x);
    _tensor[x].setSubmatrix(slab, row);
  }
}

//////////////////////////////////////////////////////////////////////
// take the product w.r.t. mode one, i.e. matrix-matrix multiply with each slab,
// but with the tensor's slabs transposed
//////////////////////////////////////////////////////////////////////
TENSOR3 TENSOR3::modeOneTransposeProduct(const MATRIX& m) const
{
  vector<MATRIX> slabs;
  for (int x = 0; x < _slabs; x++)
  {
    assert(_tensor[x].rows() == m.rows());
    slabs.push_back(_tensor[x] ^ m);
  }

  return TENSOR3(slabs);
}

//////////////////////////////////////////////////////////////////////
// take the product w.r.t. mode one, i.e. matrix-matrix multiply with each slab,
// but with the tensor's slabs transposed
//////////////////////////////////////////////////////////////////////
TENSOR3 TENSOR3::modeOneProduct(const MATRIX& m) const
{
  vector<MATRIX> slabs;
  for (int x = 0; x < _slabs; x++)
  {
    assert(_tensor[x].cols() == m.rows());
    slabs.push_back(_tensor[x] * m);
  }

  return TENSOR3(slabs);
}

//////////////////////////////////////////////////////////////////////
// get the matrix corresponding to column number 'column'
//////////////////////////////////////////////////////////////////////
MATRIX TENSOR3::getColumnSlab(int column) const
{
  MATRIX final(_rows, _slabs);

  for (int z = 0; z < _slabs; z++)
    for (int x = 0; x < _rows; x++)
      final(x,z) = (*this)(x, column, z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// set the matrix corresponding to column number 'column'
//////////////////////////////////////////////////////////////////////
void TENSOR3::setColumnSlab(const MATRIX& slab, int column)
{
  for (int z = 0; z < _slabs; z++)
    for (int x = 0; x < _rows; x++)
      (*this)(x, column, z) = slab(x,z);
}

//////////////////////////////////////////////////////////////////////
// add the matrix to corresponding to column number 'column'
//////////////////////////////////////////////////////////////////////
void TENSOR3::addColumnSlab(const MATRIX& slab, int column)
{
  for (int z = 0; z < _slabs; z++)
    for (int x = 0; x < _rows; x++)
      (*this)(x, column, z) += slab(x,z);
}

//////////////////////////////////////////////////////////////////////
// write to a file stream
//////////////////////////////////////////////////////////////////////
void TENSOR3::write(FILE* file)
{
  fwrite((void*)&_rows, sizeof(int), 1, file);
  fwrite((void*)&_cols, sizeof(int), 1, file);
  fwrite((void*)&_slabs, sizeof(int), 1, file);

  for (int x = 0; x < _slabs; x++)
    _tensor[x].write(file);
}

//////////////////////////////////////////////////////////////////////
// read from a file stream
//////////////////////////////////////////////////////////////////////
void TENSOR3::read(FILE* file)
{
  fread((void*)&_rows, sizeof(int), 1, file);
  fread((void*)&_cols, sizeof(int), 1, file);
  fread((void*)&_slabs, sizeof(int), 1, file);

  if (_tensor) delete[] _tensor;
  _tensor = new MATRIX[_slabs];

  for (int x = 0; x < _slabs; x++)
    _tensor[x].read(file);
}

//////////////////////////////////////////////////////////////////////
// do an in-place subtract of the first slab from a matrix
//////////////////////////////////////////////////////////////////////
void TENSOR3::subtractCollapsedSlabs(MATRIX& matrix)
{
  assert(matrix.rows() == _rows);
  assert(matrix.cols() == _slabs);

  for (int z = 0; z < _slabs; z++)
    for (int x = 0; x < _rows; x++)
      matrix(x,z) -= (*this)(x, 0, z);
}

//////////////////////////////////////////////////////////////////////
// do an in-place add of the first slab from a matrix
//////////////////////////////////////////////////////////////////////
void TENSOR3::addCollapsedSlabs(MATRIX& matrix)
{
  assert(matrix.rows() == _rows);
  assert(matrix.cols() == _slabs);

  for (int z = 0; z < _slabs; z++)
    for (int x = 0; x < _rows; x++)
      matrix(x,z) += (*this)(x, 0, z);
}

//////////////////////////////////////////////////////////////////////
// do an in-place axpy of the first slab from a matrix
//////////////////////////////////////////////////////////////////////
void TENSOR3::axpyCollapsedSlabs(const Real& scalar, MATRIX& matrix)
{
  assert(matrix.rows() == _rows);
  assert(matrix.cols() == _slabs);

  for (int z = 0; z < _slabs; z++)
    for (int x = 0; x < _rows; x++)
      matrix(x,z) += scalar * (*this)(x, 0, z);
}
