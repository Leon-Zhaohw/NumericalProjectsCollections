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
// MATRIX.h: interface for the MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#include "MATRIX.h"

#include <memory.h>

#ifdef USING_MKL
#include <mkl_lapack.h>
#include <mkl_types.h>
#include <mkl_cblas.h>
#endif
#ifdef USING_ATLAS
extern "C" {
#include <clapack.h>
#include <cblas.h>
}
#endif

#ifdef SINGLE_PRECISION
#define REAL_EPSILON FLT_EPSILON
#else
#define DBL_EPSILON 1e-14
#define REAL_EPSILON DBL_EPSILON
#endif

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
MATRIX::MATRIX() :
  _rows(0), _cols(0), _pivots(NULL)
{
  _matrix = NULL;
}

MATRIX::MATRIX(int rows, int cols) :
  _rows(rows), _cols(cols), _pivots(NULL)
{
  _matrix = new Real[_rows * _cols];
  clear();
}

MATRIX::MATRIX(int rows, int cols, Real* data) :
  _rows(rows), _cols(cols), _pivots(NULL)
{
  _matrix = new Real[_rows * _cols];
	memcpy(_matrix, data, _rows * _cols * sizeof(Real));
}

MATRIX::MATRIX(const char* filename) :
  _pivots(NULL)
{
  read(filename);
}

MATRIX::MATRIX(const MATRIX& m) :
  _pivots(NULL)
{
  _cols = m._cols;
  _rows = m._rows;

  _matrix = new Real[_rows * _cols];
	memcpy(_matrix, m._matrix, _rows * _cols * sizeof(Real));
}

MATRIX::MATRIX(VECTOR& vec) :
  _pivots(NULL)
{
  _cols = vec.size();
  _rows = vec.size();

  _matrix = new Real[_rows * _cols];
  clear();

  for (int x = 0; x < vec.size(); x++)
    (*this)(x,x) = vec(x);
}

MATRIX::MATRIX(vector<VECTOR>& cols) :
  _pivots(NULL)
{
  _cols = cols.size();
  _rows = cols[0].size();
  _matrix = new Real[_rows * _cols];

  for (unsigned int x = 0; x < cols.size(); x++)
  {
    assert(cols[x].size() == _rows);
    for (int y = 0; y < _rows; y++)
      (*this)(y,x) = cols[x][y];
  }
}

MATRIX::MATRIX(MATRIX3 matrix3) :
  _pivots(NULL)
{
  _cols = 3;
  _rows = 3;

  _matrix = new Real[9];
  clear();

  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
      (*this)(x,y) = matrix3(x,y);
}

MATRIX::MATRIX(const MATRIX3& rotation, const VEC3F& translation) :
  _pivots(NULL)
{
  _cols = 4;
  _rows = 4;

  _matrix = new Real[16];
  clear();

  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
      (*this)(x,y) = rotation(x,y);

  (*this)(0,3) = translation[0];
  (*this)(1,3) = translation[1];
  (*this)(2,3) = translation[2];
  (*this)(3,3) = 1.0;
}

MATRIX::~MATRIX()
{
  delete[] _matrix;
  delete[] _pivots;
}

//////////////////////////////////////////////////////////////////////
// wipe the whole matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::clear()
{
  memset(_matrix, 0, _rows * _cols * sizeof(Real));
}

//////////////////////////////////////////////////////////////////////
// resize and wipe the matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::resizeAndWipe(int rows, int cols)
{
  if (_rows != rows || _cols != cols)
  {
    delete[] _matrix;
    _rows = rows;
    _cols = cols;
    _matrix = new Real[_rows * _cols];
  }
  clear();
}

//////////////////////////////////////////////////////////////////////
// write the matrix to a file
//////////////////////////////////////////////////////////////////////
void MATRIX::write(const char* filename)
{
  FILE* file;
  file = fopen(filename, "wb");

  // write dimensions
  fwrite((void*)&_rows, sizeof(int), 1, file);
  fwrite((void*)&_cols, sizeof(int), 1, file);

  // always write out as a double
  if (sizeof(Real) != sizeof(double))
  {
    double* matrixDouble = new double[_rows * _cols];
    for (int x = 0; x < _rows * _cols; x++)
      matrixDouble[x] = _matrix[x];

    fwrite((void*)matrixDouble, sizeof(double), _rows * _cols, file);
    delete[] matrixDouble;
    fclose(file);
  }
  else
    fwrite((void*)_matrix, sizeof(Real), _rows * _cols, file);
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read matrix from a file
//////////////////////////////////////////////////////////////////////
bool MATRIX::read(const char* filename)
{
  FILE* file;
  file = fopen(filename, "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << filename << " not found! " << endl;
    return false;
  }

  // read dimensions
  fread((void*)&_rows, sizeof(int), 1, file);
  fread((void*)&_cols, sizeof(int), 1, file);

  // read data
  if (_matrix) delete[] _matrix;
  _matrix = new Real[_rows * _cols];

  // always read in a double
  if (sizeof(Real) != sizeof(double))
  {
    double* matrixDouble = new double[_rows * _cols];
    fread((void*)matrixDouble, sizeof(double), _rows * _cols, file);
    for (int x = 0; x < _rows * _cols; x++)
      _matrix[x] = matrixDouble[x];
    delete[] matrixDouble;
  }
  else 
    fread((void*)_matrix, sizeof(Real), _rows * _cols, file);
  fclose(file);

  return true;
}

//////////////////////////////////////////////////////////////////////
// read matrix from a file
//////////////////////////////////////////////////////////////////////
bool MATRIX::readSubcols(const char* filename, int maxCols)
{
  FILE* file;
  file = fopen(filename, "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << filename << " not found! " << endl;
    return false;
  }

  // read dimensions
  fread((void*)&_rows, sizeof(int), 1, file);
  fread((void*)&_cols, sizeof(int), 1, file);

  int totalCols = _cols;
  if (_cols > maxCols) _cols = maxCols;
  cout << " Reading " << _cols << " of " << totalCols << " columns" << endl;

  // read data
  if (_matrix) delete[] _matrix;
  _matrix = new Real[_rows * _cols];

  // always read in a double
  for (int y = 0; y < _rows; y++)
    for (int x = 0; x < totalCols; x++)
    {
      double entry = 0.0;
      fread((void*)&entry, sizeof(double), 1, file);
      if (x < _cols)
        _matrix[x + y * _cols] = entry;
    } 
  fclose(file);

  return true;
}

//////////////////////////////////////////////////////////////////////
// read matrix from a file
//////////////////////////////////////////////////////////////////////
void MATRIX::read(FILE* file)
{
  int oldRows = _rows;
  int oldCols = _cols;

  // read dimensions
  fread((void*)&_rows, sizeof(int), 1, file);
  fread((void*)&_cols, sizeof(int), 1, file);

  // read data -- if the size is the same, try to preserve
  // the pointer address
  if (oldRows != _rows || oldCols != _cols)
  {
    delete[] _matrix;
    _matrix = new Real[_rows * _cols];
  }

  // always read in a double
  if (sizeof(Real) != sizeof(double))
  {
    double* matrixDouble = new double[_rows * _cols];
    fread((void*)matrixDouble, sizeof(double), _rows * _cols, file);
    for (int x = 0; x < _rows * _cols; x++)
      _matrix[x] = matrixDouble[x];
    delete[] matrixDouble;
  }
  else 
    fread((void*)_matrix, sizeof(Real), _rows * _cols, file);
}

//////////////////////////////////////////////////////////////////////
// write matrix to a file
//////////////////////////////////////////////////////////////////////
void MATRIX::write(FILE* file)
{
  // write dimensions
  fwrite((void*)&_rows, sizeof(int), 1, file);
  fwrite((void*)&_cols, sizeof(int), 1, file);

  // always write in a double
  if (sizeof(Real) != sizeof(double))
  {
    double* matrixDouble = new double[_rows * _cols];
    for (int x = 0; x < _rows * _cols; x++)
      matrixDouble[x] = _matrix[x];
    fwrite((void*)matrixDouble, sizeof(double), _rows * _cols, file);
    delete[] matrixDouble;
  }
  else 
    fwrite((void*)_matrix, sizeof(Real), _rows * _cols, file);
}

//////////////////////////////////////////////////////////////////////
// return transpose of current matrix
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::transpose() const
{
  MATRIX toReturn(_cols, _rows);

  for (int y = 0; y < _cols; y++)
    for (int x = 0; x < _rows; x++)
      toReturn(y,x) = (*this)(x,y);

  return toReturn;
}

//////////////////////////////////////////////////////////////////////
// Put the transpose of this matrix in output
//////////////////////////////////////////////////////////////////////
void MATRIX::transpose( MATRIX &output )
{
  for (int y = 0; y < _cols; y++)
    for (int x = 0; x < _rows; x++)
      output(y,x) = (*this)(x,y);
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator*(const MATRIX& A, const VECTOR& x) 
{
  assert(A.cols() == x.size());
  VECTOR y(A.rows());

#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasNoTrans, 
		A.rows(), A.cols(),
		1.0f, A.dataConst(), A.cols(), x.dataConst(), 1, 
		0.0f, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasNoTrans, 
		A.rows(), A.cols(),
		1.0, A.dataConst(), A.cols(), x.dataConst(), 1, 
		0.0, y.data(), 1);
#endif

  return y;
}

//////////////////////////////////////////////////////////////////////
// Vector-matrix multiply
//
// This is really just A^T * x -- add cblas call later
//////////////////////////////////////////////////////////////////////
VECTOR operator*(VECTOR& x, MATRIX& A)
{
  assert(A.rows() == x.size());

  VECTOR y(A.cols());
  for (int i = 0; i < A.cols(); i++)
    for (int j = 0; j < A.rows(); j++)
      y[i] += A(j, i) * x(j);
  return y;
}

//////////////////////////////////////////////////////////////////////
// Scale matrix
//////////////////////////////////////////////////////////////////////
MATRIX operator*(const MATRIX& A, Real alpha) 
{
  MATRIX y(A.rows(), A.cols());

  /*
  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      y(i,j) = A(i, j) * alpha;
      */
  y.axpy(alpha, A);

  return y;
}

//////////////////////////////////////////////////////////////////////
// Scale matrix
//////////////////////////////////////////////////////////////////////
MATRIX operator*(Real alpha, const MATRIX& A) 
{
  MATRIX y(A.rows(), A.cols());

  /*
  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      y(i,j) = A(i, j) * alpha;
      */
  y.axpy(alpha, A);

  return y;
}

//////////////////////////////////////////////////////////////////////
// scale matrix by a constant
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator*=(const Real& alpha)
{
#ifdef SINGLE_PRECISION
  cblas_sscal(_cols * _rows, alpha, _matrix, 1);
#else
  cblas_dscal(_cols * _rows, alpha, _matrix, 1);
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Matrix-Matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(const MATRIX& A, const MATRIX& B) 
{
  assert(A.cols() == B.rows());

  MATRIX y(A.rows(), B.cols());

#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			A.rows(), 
      B.cols(),
      A.cols(), 
			1.0f,
			A.dataConst(), A.cols(),
			B.dataConst(), B.cols(),
			0.0f,
			y.data(), y.cols());
#else
	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			A.rows(), 
      B.cols(),
      A.cols(), 
			1.0,
			A.dataConst(), A.cols(),
			B.dataConst(), B.cols(),
			0.0,
			y.data(), y.cols());
#endif
  return y;
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply where A is transposed
//////////////////////////////////////////////////////////////////////
VECTOR operator^(const MATRIX& A, const VECTOR& x) 
{
  assert(A.rows() == x.size());

  VECTOR y(A.cols());

#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasTrans, 
		A.rows(), A.cols(), 
		1.0f, A.dataConst(), A.cols(), x.dataConst(), 1, 
		0.0f, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasTrans, 
		A.rows(), A.cols(), 
		1.0, A.dataConst(), A.cols(), x.dataConst(), 1, 
		0.0, y.data(), 1);
#endif
  return y;
}

//////////////////////////////////////////////////////////////////////
// Matrix^T -Matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator^(const MATRIX& A, const MATRIX& B) 
{
  assert(A.rows() == B.rows());

  MATRIX y(A.cols(), B.cols());

#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			CblasTrans,
			CblasNoTrans,
			A.cols(), 
      B.cols(), 
      A.rows(),
			1.0f,
			A.dataConst(), A.cols(),
			B.dataConst(), B.cols(),
			0.0f,
			y.data(), y.cols());
#else
	cblas_dgemm(
			CblasRowMajor,
			CblasTrans,
			CblasNoTrans,
			A.cols(), 
      B.cols(), 
      A.rows(),
			1.0,
			A.dataConst(), A.cols(),
			B.dataConst(), B.cols(),
			0.0,
			y.data(), y.cols());
#endif
  return y;
}

/*
// moved to MATRIX.cpp to support variable precision cout
//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, const MATRIX& matrix)
{
  out << "[" << endl; 
  for (int row = 0; row < matrix.rows(); row++)
  {
    for (int col = 0; col < matrix.cols(); col++)
      out << matrix(row, col) << " ";
    out << endl;
  }
  out << "]" << endl;
  return out;
}
*/

//////////////////////////////////////////////////////////////////////
// Deep copy equality operator
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator=(const MATRIX3 m)
{
  if (3 != _cols || 3 != _rows)
  {
    delete[] _matrix;
    _cols = 3;
    _rows = 3;

    _matrix = new Real[9];
  }
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
      (*this)(x,y) = m(x,y);
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Deep copy equality operator
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator=(const MATRIX m)
{
  if (m._cols != _cols || m._rows != _rows)
  {
    delete[] _matrix;
    _cols = m._cols;
    _rows = m._rows;

    _matrix = new Real[_rows * _cols];
  }
	memcpy(_matrix, m._matrix, _rows * _cols * sizeof(Real));
  return *this;
}

//////////////////////////////////////////////////////////////////////
// self minus
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator-=(const MATRIX& m)
{
  assert(m._cols == _cols && m._rows == _rows);

#ifdef SINGLE_PRECISION
	cblas_saxpy(_cols * _rows, -1.0f, m._matrix, 1, _matrix, 1);
#else
	cblas_daxpy(_cols * _rows, -1.0, m._matrix, 1, _matrix, 1);
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// self plus
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator+=(const MATRIX& m)
{
  assert(m._cols == _cols && m._rows == _rows);

#ifdef SINGLE_PRECISION
	cblas_saxpy(_cols * _rows, 1.0f, m._matrix, 1, _matrix, 1);
#else
	cblas_daxpy(_cols * _rows, 1.0, m._matrix, 1, _matrix, 1);
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Return the matrix diagonal
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::diagonal()
{
  int minDim = (_rows > _cols) ? _cols : _rows;
  VECTOR diag(minDim);
  for (int x = 0; x < minDim; x++)
    diag(x) = (*this)(x,x);

  return diag;
}

//////////////////////////////////////////////////////////////////////
// stomp the current matrix with the given matrix starting at "row". 
// It is your responsibility to ensure that you don't fall off the 
// end of this matrix.
//////////////////////////////////////////////////////////////////////
void MATRIX::setSubmatrix(MATRIX& matrix, int row)
{
  assert(row >= 0);
  assert(row < _rows);

  int totalSize = matrix.rows() * matrix.cols();
  int index = row * _cols;

  for (int x = 0; x < totalSize; x++, index++)
    _matrix[index] = matrix._matrix[x];
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B += alpha * A, where B is this matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::axpy(const Real alpha, const MATRIX& A)
{
  assert(_rows == A.rows() && _cols == A.cols());

#ifdef SINGLE_PRECISION
	cblas_saxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#else
	cblas_daxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B = alpha * A, where B is this matrix, and 
// current contents of B are stomped
//////////////////////////////////////////////////////////////////////
void MATRIX::clearingAxpy(Real alpha, MATRIX& A)
{
  assert(_rows == A.rows() && _cols == A.cols());

  memset(_matrix, 0, _rows * _cols * sizeof(Real));
#ifdef SINGLE_PRECISION
	cblas_saxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#else
	cblas_daxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// BLAS gemm operation: C += alpha * A * B where C is this matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::gemm(Real alpha, const MATRIX& A, const MATRIX& B, Real beta,
                  bool transposeA, bool transposeB)
{
  int requiredRows = transposeA ? A.cols() : A.rows();
  int requiredCols = transposeB ? B.rows() : B.cols();

  assert(_rows == requiredRows && _cols == requiredCols);

#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			transposeA ? CblasTrans : CblasNoTrans,
			transposeB ? CblasTrans : CblasNoTrans,
			A.rows(), 
      B.cols(),
      A.cols(), 
			alpha,
			A.dataConst(), A.cols(),
			B.dataConst(), B.cols(),
			beta,
			_matrix, _cols);
#else
	cblas_dgemm(
			CblasRowMajor,
			transposeA ? CblasTrans : CblasNoTrans,
			transposeB ? CblasTrans : CblasNoTrans,
			A.rows(), 
      B.cols(),
      A.cols(), 
			alpha,
			A.dataConst(), A.cols(),
			B.dataConst(), B.cols(),
			beta,
			_matrix, _cols);
#endif
}

/*
 * Untested -- don't uncomment until a test case comes up
 * 
//////////////////////////////////////////////////////////////////////
// BLAS gemv operation: y += A * x where A is this matrix
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::gemv(VECTOR& x)
{
  if (x.size() != _rows)
    cout << __FILE__ << " " << __LINE__ << " : Matrix-Vector dimensions do not match!" << endl;

  VECTOR y(x.size());
  for (int j = 0; j < _cols; j++)
    for (int i = 0; i < _rows; i++)
      y(j) += (*this)(i,j) * x(j);

  return y;
}
*/

//////////////////////////////////////////////////////////////////////
// BLAS gemv operation: y += A * x where A is this matrix
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::gemv(VEC3F& x)
{
  VECTOR y(_rows);
#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		1.0f, _matrix, _cols, x, 1, 
		0.0f, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		1.0, _matrix, _cols, x, 1, 
		0.0, y.data(), 1);
#endif
  return y;
}

//////////////////////////////////////////////////////////////////////
// BLAS gemv operation: y += alpha * A * x where A is this matrix
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::gemv(Real alpha, VEC3F& x)
{
  VECTOR y(_rows);
#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		alpha, _matrix, _cols, x, 1, 
		0.0f, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		alpha, _matrix, _cols, x, 1, 
		0.0, y.data(), 1);
#endif
  return y;
}

//////////////////////////////////////////////////////////////////////
// BLAS gemv operation:
//  y = alpha * A * x + beta * y
//////////////////////////////////////////////////////////////////////
void MATRIX::gemvInplace(Real alpha, const VECTOR &x, VECTOR &y,
                         Real beta, bool transpose)
{
#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		transpose ? CblasTrans : CblasNoTrans,
    _rows, _cols,
		alpha, _matrix, _cols, x.dataConst(), 1, 
		beta, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		transpose ? CblasTrans : CblasNoTrans,
    _rows, _cols,
		alpha, _matrix, _cols, x.dataConst(), 1, 
		beta, y.data(), 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// A bit of a hack
//////////////////////////////////////////////////////////////////////
void MATRIX::gemvInplace(Real alpha, VEC3F &x, VECTOR &y,
                         Real beta, bool transpose)
{
#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		transpose ? CblasTrans : CblasNoTrans,
    _rows, _cols,
		alpha, _matrix, _cols, x, 1, 
		beta, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		transpose ? CblasTrans : CblasNoTrans,
    _rows, _cols,
		alpha, _matrix, _cols, x, 1, 
		beta, y.data(), 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-Matrix^T multiply
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::timesTranspose(MATRIX& B) 
{
  MATRIX& A = *this;
  assert(A.cols() == B.cols());

  MATRIX y(A.rows(), B.rows());

#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasTrans,
			A.rows(), 
      B.rows(), 
      A.cols(),
			1.0f,
			A.data(), A.cols(),
			B.data(), B.cols(),
			0.0f,
			y.data(), y.cols());
#else
	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasTrans,
			A.rows(), 
      B.rows(), 
      A.cols(),
			1.0,
			A.data(), A.cols(),
			B.data(), B.cols(),
			0.0,
			y.data(), y.cols());
#endif
  return y;
}

//////////////////////////////////////////////////////////////////////
// BLAS gemm operation: C = alpha * A * B where C is this matrix and
// current contents of C are stomped
//////////////////////////////////////////////////////////////////////
void MATRIX::clearingGemm(Real alpha, MATRIX& A, MATRIX& B)
{
  assert(_rows == A.rows() && _cols == B.cols());

#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			A.rows(), 
      B.cols(),
      A.cols(), 
			alpha,
			A.data(), A.cols(),
			B.data(), B.cols(),
			0.0f,
			_matrix, _cols);
#else
	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			A.rows(), 
      B.cols(),
      A.cols(), 
			alpha,
			A.data(), A.cols(),
			B.data(), B.cols(),
			0.0,
			_matrix, _cols);
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//
// Do *NOT* let x == y!
//////////////////////////////////////////////////////////////////////
void MATRIX::multiplyInplace(VECTOR& x, VECTOR& y) 
{
  assert(x.data() != y.data());

#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		1.0f, _matrix, _cols, x.data(), 1, 
		0.0f, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		1.0, _matrix, _cols, x.data(), 1, 
		0.0, y.data(), 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// solve the linear system Ax = b, return x in the passed in b
//////////////////////////////////////////////////////////////////////
void MATRIX::solve(VECTOR& b)
{
  assert(_rows == b.size());

  char uplo = 'U';
  int nrhs = 1;
  int info = 0;
  int R = _rows;

#ifdef USING_MKL  
  // call the general LU solver from MKL
#ifdef SINGLE_PRECISION
  sposv(&uplo, &R, &nrhs, _matrix, &R, b.data(), &R, &info);
#else
  dposv(&uplo, &R, &nrhs, _matrix, &R, b.data(), &R, &info);
#endif
  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "Solve bombed: " << info << endl;
  }
#endif
#ifdef USING_ATLAS
  // call the general LU solver from ATLAS
#ifdef SINGLE_PRECISION
  clapack_sposv(CblasRowMajor, CblasUpper, R, nrhs, _matrix, R, b.data(), R);
#else
  clapack_dposv(CblasRowMajor, CblasUpper, R, nrhs, _matrix, R, b.data(), R);
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// solve the linear system Ax = b, return x in the passed in b
//////////////////////////////////////////////////////////////////////
void MATRIX::solveSPD(VECTOR& b)
{
  assert(_rows == b.size());

  char uplo = 'U';
  int nrhs = 1;
  int info = 0;
  int R = _rows;

#ifdef USING_MKL  
  // call the general LU solver from MKL
#ifdef SINGLE_PRECISION
  sposv(&uplo, &R, &nrhs, _matrix, &R, b.data(), &R, &info);
#else
  dposv(&uplo, &R, &nrhs, _matrix, &R, b.data(), &R, &info);
#endif
#endif
#ifdef USING_ATLAS
  // call the general LU solver from ATLAS
#ifdef SINGLE_PRECISION
  clapack_sposv(CblasRowMajor, CblasUpper, R, nrhs, _matrix, R, b.data(), R);
#else
  clapack_dposv(CblasRowMajor, CblasUpper, R, nrhs, _matrix, R, b.data(), R);
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// solve for the eigensystem of the matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem(VECTOR& eigenvalues, MATRIX& eigenvectors)
{
  // basic error checking
  if (_rows != _cols)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Matrix must be square to get eigenvalues! " << endl;
    return;
  }

  // resize result space
  // FIXME: Just check for the right size
  eigenvalues.resizeAndWipe(_rows);
  eigenvectors.resizeAndWipe(_rows, _rows);

  int rowsize = _rows;
  int worksize = 5 * _rows;
  Real* work = new Real[worksize];
  Real* valuesReal = new Real[2 * _rows];
  Real* valuesImag = valuesReal + _rows;
  Real* vectors = new Real[_rows * _rows];

  // the actual LAPACK call
  int error;
#ifdef SINGLE_PRECISION
  sgeev("V","N", &rowsize, _matrix, &rowsize, 
        valuesReal, valuesImag, 
        vectors, &rowsize, NULL, &rowsize,
        work, &worksize, &error);
#else
  dgeev("V","N", &rowsize, _matrix, &rowsize, 
        valuesReal, valuesImag, 
        vectors, &rowsize, NULL, &rowsize,
        work, &worksize, &error);
#endif

  if (error != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " eigenvalue solver bombed!" << endl;
  }

  // copy out results
  for (int x = 0; x < _rows; x++)
    eigenvalues(x) = valuesReal[x];

  for (int x = 0; x < _rows; x++)
    for (int y = 0; y < _rows; y++)
      eigenvectors(x,y) = vectors[x + y * _cols];
 
  // cleanup
  delete[] work;
  delete[] valuesReal;
  delete[] vectors;
}

//////////////////////////////////////////////////////////////////////
// perform PCA on this matrix
// Note this destroys the matrix!
//////////////////////////////////////////////////////////////////////
void MATRIX::verbosePCA(MATRIX& components, VECTOR& values, MATRIX& U, VECTOR& S, MATRIX& VT)
{
  subtractRowMeans();

  MATRIX& data = *this;
  int rows = data.rows();
  int cols = data.cols();

  // do a full SVD on the whole shebang
  data.SVD(U,S,VT);

  // store everything
  components.resizeAndWipe(rows, cols);
  values.resizeAndWipe(cols);
  for (int x = 0; x < cols; x++)
  {
    for (int y = 0; y < rows; y++)
      components(y,x) = U(y,x);
    values[x] = S[x] * S[x];
  }
}

//////////////////////////////////////////////////////////////////////
// perform PCA on this matrix
// Note this destroys the matrix!
//////////////////////////////////////////////////////////////////////
void MATRIX::PCA(MATRIX& components, VECTOR& values)
{
  subtractRowMeans();
  meanNormalizedPCA(components, values);
}

//////////////////////////////////////////////////////////////////////
// perform PCA on this matrix
// Note this destroys the matrix!
//////////////////////////////////////////////////////////////////////
void MATRIX::meanNormalizedPCA(MATRIX& components, VECTOR& values)
{
  MATRIX& data = *this;
  int rows = data.rows();
  int cols = data.cols();

  MATRIX U;
  MATRIX VT;
  VECTOR S;

  // do a full SVD on the whole shebang
  data.SVD(U,S,VT);

  // store everything
  components.resizeAndWipe(rows, cols);
  values.resizeAndWipe(cols);
  for (int x = 0; x < cols; x++)
  {
    for (int y = 0; y < rows; y++)
      components(y,x) = U(y,x);
    values[x] = S[x] * S[x];
  }
}

//////////////////////////////////////////////////////////////////////
// subtract the mean from each row in prep for PCA
//////////////////////////////////////////////////////////////////////
void MATRIX::subtractRowMeans()
{
  MATRIX& data = *this;

  // subtract out the means
  int rows = data.rows();
  int cols = data.cols();
  for (int x = 0; x < rows; x++)
  {
    Real mean = 0.0;
    for (int y = 0; y < cols; y++)
      mean += data(x,y);
    mean /= cols;

    for (int y = 0; y < cols; y++)
      data(x,y) -= mean;
  }

  // divide through by sqrt(N - 1)
  data *= 1.0 / sqrt(cols - 1);
}

//////////////////////////////////////////////////////////////////////
// copy this matrix to MATRIX3 type
//////////////////////////////////////////////////////////////////////
void MATRIX::copiesInto(MATRIX3& matrix3)
{
  assert(_rows == 3 && _cols == 3);

  /*
  if (_rows != 3 || _cols != 3)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "Trying to copy a MATRIX that is not 3x3 to a MATRIX3!" << endl;
    return;
  }
  */
  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
      matrix3(x,y) = (*this)(x,y);
}

//////////////////////////////////////////////////////////////////////
// Copy the given matrix (or matrix data) directly,
// only reallocating if necessary.
//////////////////////////////////////////////////////////////////////
void MATRIX::copyInplace( const Real *data, int rows, int cols )
{
  if ( rows != _rows || cols != _cols )
  {
    resizeAndWipe( rows, cols );
  }

	memcpy(_matrix, data, _rows * _cols * sizeof(Real));
}

//////////////////////////////////////////////////////////////////////
// stomp the contents and convert this to an identity matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::setToIdentity()
{
  this->clear();
  int least = (_rows > _cols) ? _cols : _rows;
  for (int x = 0; x < least; x++)
    (*this)(x,x) = 1.0;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void MATRIX::eigenMinMax(Real& minEig, Real& maxEig)
{
  assert(_rows == _cols);
  
  VECTOR values(_rows);
  MATRIX vectors(_rows, _cols);

  eigensystem(values, vectors);

  minEig = values[0];
  maxEig = values[0];
  for (int x = 1; x < values.size(); x++)
  {
    if (values[x] < minEig) minEig = values[x];
    if (values[x] > maxEig) maxEig = values[x];
  }
}

//////////////////////////////////////////////////////////////////////
// Solve for eigenvalues in the symmetric case.
// Behaviour is undefined for an unsymmetric matrix.
// Puts all eigenvalues in the specified range (vl, vu]
// in the provided vector and returns the number of
// eigenvalues computed.
// 
// Destroys the matrix.
//
// This call expects a workspace.
//////////////////////////////////////////////////////////////////////
int MATRIX::eigensystemSymmetricRange(Real *eigenvalues,
                                      Real *work, int lwork,
                                      int *iwork, int liwork,
                                      Real vl, Real vu,
                                      Real tolerance )
{
  assert( _rows == _cols );

  int n = _rows;
  int zero = 0;
  int one = 1;
  int numFound = -1;
  Real *w = eigenvalues;
  int info = -1;

#ifdef SINGLE_PRECISION
  // FIXME: Fill this in once the double precision version is working...
#else
#if 0
  dsyevr( "N",           // Eigenvalues only
          "V",           // Only within a certain interval
          "U",           // Upper triangular part stored ("L" would work too)
          &n,            // Dimension of matrix
          _matrix,       // Actual data
          &n,            // Leading dimension of _matrix
          &vl, &vu,      // Range to find eigenvalues in
          &zero, &zero,  // Index range for eigenvalues (not used)
          &tolerance,    // Solution tolerance
          &numFound,     // Number of eigenvalues found
          w,             // Stores eigenvalues on return
          NULL,          // Eigenvector storage (not needed)
          &one,          // Leading dimension of vector array (not needed)
          NULL,          // Eigenvector support (not needed)
          work,          // Work array
          &lwork,        // Work array dimension
          iwork,         // Integer work array
          &liwork,       // Integer work array dimension
          &info          // Return status
        );
#endif
  dsyevx( "N",           // Eigenvalues only
          "V",           // Only within a certain interval
          "U",           // Upper triangular part stored ("L" would work too)
          &n,            // Dimension of matrix
          _matrix,       // Actual data
          &n,            // Leading dimension of _matrix
          &vl, &vu,      // Range to find eigenvalues in
          &zero, &zero,  // Index range for eigenvalues (not used)
          &tolerance,    // Solution tolerance
          &numFound,     // Number of eigenvalues found
          w,             // Stores eigenvalues on return
          NULL,          // Eigenvector storage (not needed)
          &one,          // Leading dimension of vector array (not needed)
          work,          // Work array
          &lwork,        // Work array dimension
          iwork,         // Integer work array
          NULL,          // Failure array (not needed)
          &info          // Return status
        );

  if ( info != 0 )
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " symmetric eigenvalue solver bombed!" << endl;
    numFound = -1;

    if ( info < 0 )
    {
      cout << "Parameter " << (-info) << " had an illegal value!" << endl;
    }
  }
#endif

  return numFound;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
Real MATRIX::conditionNumber()
{
  Real minEig;
  Real maxEig;
  eigenMinMax(minEig, maxEig);

  return fabs(maxEig / minEig);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void MATRIX::normalizeColumns()
{
  for (int x = 0; x < _cols; x++)
  {
    Real sumSq = 0.0;
    for (int y = 0; y < _rows; y++)
      sumSq += (*this)(x,y) * (*this)(x,y);

    Real invSumSq = 1.0 / sqrt(sumSq);
    for (int y = 0; y < _rows; y++)
      (*this)(x,y) *= invSumSq;
  }
}

//////////////////////////////////////////////////////////////////////
// quickly solve for the eigensystem of a 2x2 matrix
//
// Adapted from: 
// http://www.cds.caltech.edu/~andrea/research/sw/2x2_eigen.html
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem2x2(VECTOR& eigenvalues, MATRIX& eigenvectors)
{
  Real A = (*this)(0,0);
  Real B = (*this)(0,1);
  Real C = (*this)(1,0);
  Real D = (*this)(1,1);

  eigenvalues.resizeAndWipe(2);
  eigenvectors.resizeAndWipe(2,2);
  
  if (B * C <= 0.1e-20) {
    eigenvalues(0) = A; 
    eigenvectors(0,0) = 1; 
    eigenvectors(1,0) = 0;
    eigenvalues(1) = D; 
    eigenvectors(0,1) = 0; 
    eigenvectors(1,1) = 1;
    return;
  }

  Real tr = A + D;
  Real det = A * D - B * C;
  Real S = sqrt(tr * tr * 0.25 - det );
  eigenvalues(0) = tr * 0.5 + S;
  eigenvalues(1) = tr * 0.5 - S;

  Real temp = (A - D) * (A - D) * 0.25 + B * C;
  Real SS = (temp < 0.0) ? 0.0 : sqrt(temp);
  if (A - D < 0.0) {
    eigenvectors(0,0) = C;
    eigenvectors(1,0) = -(A - D) * 0.5 + SS;
    eigenvectors(0,1) = (A - D) * 0.5 - SS;
    eigenvectors(1,1) = B;
  } 
  else {
    eigenvectors(0,1) = C;
    eigenvectors(1,1) = -(A - D) * 0.5 - SS;
    eigenvectors(0,0) = (A - D) * 0.5 + SS;
    eigenvectors(1,0) = B;
  }

  Real n1 = sqrt(eigenvectors(0,0) * eigenvectors(0,0) +
                 eigenvectors(1,0) * eigenvectors(1,0));
  Real inv = 1.0 / n1;
  eigenvectors(0,0) *= inv; 
  eigenvectors(1,0) *= inv;
  Real n2 = sqrt(eigenvectors(0,1) * eigenvectors(0,1) +
                 eigenvectors(1,1) * eigenvectors(1,1));
  inv = 1.0 / n2;
  eigenvectors(0,1) *= inv; 
  eigenvectors(1,1) *= inv;
}

//////////////////////////////////////////////////////////////////////
// squared sum of the whole matrix
//////////////////////////////////////////////////////////////////////
Real MATRIX::sum2()
{
  Real sum = 0.0;
  for (int y = 0; y < _cols; y++)
    for (int x = 0; x < _rows; x++)
      sum += (*this)(x,y) * (*this)(x,y);
  return sum;
}

//////////////////////////////////////////////////////////////////////
// 1 norm of the matrix
// - maximum absolute column sum
//////////////////////////////////////////////////////////////////////
Real MATRIX::norm1()
{
  Real maxSum = 0.0;
  for (int x = 0; x < _cols; x++)
  {
    Real colSum = 0.0;
    for (int y = 0; y < _rows; y++)
      colSum += fabs((*this)(y,x));
    if (colSum > maxSum)
      maxSum = colSum;
  }
  return maxSum;
}

//////////////////////////////////////////////////////////////////////
// inf norm of the matrix
// - maximum absolute row sum
//////////////////////////////////////////////////////////////////////
Real MATRIX::normInf()
{
  Real maxSum = 0.0;
  for (int x = 0; x < _rows; x++)
  {
    Real rowSum = 0.0;
    for (int y = 0; y < _cols; y++)
      rowSum += fabs((*this)(x,y));
    if (rowSum > maxSum)
      maxSum = rowSum;
  }
  return maxSum;
}

//////////////////////////////////////////////////////////////////////
// sum of two matrices
//////////////////////////////////////////////////////////////////////
MATRIX operator+(const MATRIX& A, const MATRIX& B)
{
  MATRIX result(A.rows(), A.cols());
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++)
      result(x,y) = A(x,y) + B(x,y);
  return result;
}

//////////////////////////////////////////////////////////////////////
// difference of two matrices
//////////////////////////////////////////////////////////////////////
MATRIX operator-(const MATRIX& A, const MATRIX& B)
{
  MATRIX result(A.rows(), A.cols());
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++)
      result(x,y) = A(x,y) - B(x,y);
  return result;
}

//////////////////////////////////////////////////////////////////////
// Get LU factorization of the matrix
//////////////////////////////////////////////////////////////////////
bool MATRIX::factorLU()
{
  if (_pivots) delete[] _pivots;
  int smaller = (_rows < _cols) ? _rows : _cols;
  int info;
  _pivots = new int[smaller]; 
  
#ifdef SINGLE_PRECISION
  sgetrf(&_rows, &_cols, _matrix, &_rows, _pivots, &info);
#else
  dgetrf(&_rows, &_cols, _matrix, &_rows, _pivots, &info);
#endif

  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " LU factorization failed!" << endl;
    if (info < 0)
      cout << " Argument " << info << " is bad" << endl;
    else
      cout << " Matrix is not positive definite " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////
// Use LU factorization to solve for rhs b
//////////////////////////////////////////////////////////////////////
void MATRIX::solveLU(VECTOR& b)
{
  assert(_rows == b.size());

  if (_pivots == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " You need to call solveLU before you can call factorLU!" << endl;
    return;
  }

  // do we want to solve the transposed problem? (No)
  char transposed = 'N';

  // how many B's are we solving for? (one)
  int howManyBs = 1;

  int info;
  
#ifdef SINGLE_PRECISION
  sgetrs(&transposed, &_rows, &howManyBs, _matrix, &_rows, _pivots, b.data(), &_rows, &info);
#else
  dgetrs(&transposed, &_rows, &howManyBs, _matrix, &_rows, _pivots, b.data(), &_rows, &info);
#endif

  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " LU solve failed!" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// Use LU factorization to solve for rhs b
//////////////////////////////////////////////////////////////////////
void MATRIX::solveLU(MATRIX& B)
{
  assert(_cols == B.rows());

  if (_pivots == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " You need to call factorLU before you can call solveLU!" << endl;
    return;
  }

  // do we want to solve the transposed problem? (No)
  char transposed = 'N';

  // how many B's are we solving for? (the total number of columns in B)
  int howManyBs = B.cols();

  int info;
 
  // repack for the LAPACK call 
  B.formatFortran();
#ifdef SINGLE_PRECISION
  sgetrs(&transposed, &_rows, &howManyBs, _matrix, &_rows, _pivots, B.data(), &_rows, &info);
#else
  dgetrs(&transposed, &_rows, &howManyBs, _matrix, &_rows, _pivots, B.data(), &_rows, &info);
#endif
  // unpack for the LAPACK call
  B.unformatFortran();

  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " LU solve failed!" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// explicitly compute the inverse, not just the LU decomposition
//////////////////////////////////////////////////////////////////////
bool MATRIX::invert()
{
  // first do the LU decomposition
  if (_pivots) delete[] _pivots;
  int smaller = (_rows < _cols) ? _rows : _cols;
  int info;
  _pivots = new int[smaller]; 

#ifdef SINGLE_PRECISION
  sgetrf(&_rows, &_cols, _matrix, &_rows, _pivots, &info);
#else
  dgetrf(&_rows, &_cols, _matrix, &_rows, _pivots, &info);
#endif
  
  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " LU factorization failed!" << endl;

    if (info < 0)
      cout << " The " << -1 * info << "th argument is invalid! " << endl;
    else
      cout << " The " << info << "th diagonal entry in U is zero! " << endl;
    return false;
  }

  // compute the final inv(A)
  int lwork = _cols;
  Real* work = new Real[lwork];
#ifdef SINGLE_PRECISION
  sgetri(&_rows, _matrix, &_rows, 
          _pivots, work, &lwork, &info);
#else
  dgetri(&_rows, _matrix, &_rows, 
          _pivots, work, &lwork, &info);
#endif

  delete[] work;

  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " LU substitution failed!" << endl;
    return false;
  }
  this->unformatFortran();
  return true;
}

//////////////////////////////////////////////////////////////////////
// Get Cholesky factorization of the matrix
//////////////////////////////////////////////////////////////////////
bool MATRIX::factorCholesky()
{
  // do the upper or lower factorization? (lower)
  char upperOrLower = 'L';
  int info;
  
#ifdef SINGLE_PRECISION
  spotrf(&upperOrLower, &_rows, _matrix, &_rows, &info); 
#else
  dpotrf(&upperOrLower, &_rows, _matrix, &_rows, &info); 
#endif

  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Cholesky factorization failed!" << endl;
    if (info < 0)
      cout << " Argument " << info << " is bad" << endl;
    else
      cout << " Matrix is not symmetric positive definite " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////
// Solve matrix using precomputed Cholesky factorization
//////////////////////////////////////////////////////////////////////
void MATRIX::solveCholesky(VECTOR& b)
{
  assert(_rows == b.size());

  // do the upper or lower factorization? (lower)
  char upperOrLower = 'L';

  // how many B's are we solving for? (one)
  int howManyBs = 1;

  int info;
 
#ifdef SINGLE_PRECISION
  spotrs(&upperOrLower, &_rows, &howManyBs, _matrix, &_rows, b.data(), &_rows, &info);
#else
  dpotrs(&upperOrLower, &_rows, &howManyBs, _matrix, &_rows, b.data(), &_rows, &info);
#endif

  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Cholesky solve failed!" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// Solve for least squares
//////////////////////////////////////////////////////////////////////
void MATRIX::solveLeastSquares(VECTOR& b)
{
  char trans = 'N';
  int info = 0;
  
  // how many B's are we solving for? (one)
  int howManyBs = 1;

  // calculate work size (assuming 1 b)
  int workSize = (_rows < _cols) ? _rows : _cols;
  workSize += (_rows > _cols) ? _rows : _cols;
  Real* work = new Real[workSize];

#ifdef SINGLE_PRECISION
  sgels(&trans, &_rows, &_cols, &howManyBs, _matrix, &_rows, b.data(), &_rows, work, &workSize, &info);
#else
  dgels(&trans, &_rows, &_cols, &howManyBs, _matrix, &_rows, b.data(), &_rows, work, &workSize, &info);
#endif

  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " element " << info << " was singular, can't compute full rank least squares!" << endl;
  }

  delete[] work;
}

//////////////////////////////////////////////////////////////////////
// Get the matrix SVD
// Note this destroys the matrix!
//////////////////////////////////////////////////////////////////////
void MATRIX::SVD(MATRIX& U, VECTOR& S, MATRIX& VT)
{
  MATRIX M(_cols, _rows);
  for (int y = 0; y < _cols; y++)
    for (int x = 0; x < _rows; x++)
      M(y,x) = (*this)(x,y);

  int bigger = (_rows > _cols) ? _rows : _cols;
  int smaller = (_rows > _cols) ? _cols: _rows;

  // resize everything to the proper size
  U.resizeAndWipe(_cols, _rows);
  VT.resizeAndWipe(_cols, _cols);
  S.resizeAndWipe(bigger);

  // compute only the singular vectors
  char jobu = 'S';
  char jobvt = 'S';
  int ldu = _rows;
  int ldvt = smaller;

  // allocate the work array
  int lwork = 3 * smaller + bigger;
  lwork = (lwork < 5 * smaller) ? 5 * smaller : lwork;
  Real* work = new Real[lwork];
 
  int info;
#ifdef SINGLE_PRECISION
  sgesvd(&jobu, &jobvt, &_rows, &_cols, M.data(), &_rows, S.data(), U.data(), &ldu, VT.data(), &ldvt, work, &lwork, &info);
#else
  dgesvd(&jobu, &jobvt, &_rows, &_cols, M.data(), &_rows, S.data(), U.data(), &ldu, VT.data(), &ldvt, work, &lwork, &info);
#endif

  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "SVD DID NOT CONVERGE" << endl;
  }
  U = U.transpose();
  VT = VT.transpose();

  delete[] work;
}

//////////////////////////////////////////////////////////////////////
// Get the trace of the matrix
//////////////////////////////////////////////////////////////////////
Real MATRIX::trace()
{
  int size = (_rows > _cols) ? _cols : _rows;
  Real result = 0.0;
  for (int x = 0; x < size; x++)
    result += (*this)(x,x);
  return result;
}

//////////////////////////////////////////////////////////////////////
// Get the absolute max entry in the matrix
//////////////////////////////////////////////////////////////////////
Real MATRIX::maxAbsEntry()
{
  Real maxFound = 0.0;

  for (int x = 0; x < _rows * _cols; x++)
    if (fabs(_matrix[x]) > maxFound)
      maxFound = fabs(_matrix[x]);

  return maxFound;
}

//////////////////////////////////////////////////////////////////////
// Get the absolute max entry in the matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::orthogonalize()
{
  // going to do this in an easy to debug but not
  // terribly efficient way
  //
  // Better to do it in place, or call MKL QR?
  for (int x = 0; x < _cols; x++)
  {
    VECTOR currentColumn = getColumn(x);

    // project out other components
    for (int y = 0; y < x; y++)
    {
      VECTOR previousColumn = getColumn(y);
      Real dot = previousColumn * currentColumn;
      currentColumn.axpy(-dot, previousColumn);
    }

    Real norm2 = currentColumn.norm2();
    currentColumn *= 1.0 / norm2;

    // save out the new orthogonalized column
    setColumn(currentColumn, x);
  }
}

//////////////////////////////////////////////////////////////////////
// Get a vector with the entires of the desired column
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::getColumn(int index)
{
  VECTOR column(_rows);
  for (int x = 0; x < _rows; x++)
    column[x] = (*this)(x, index);
  return column;
}

//////////////////////////////////////////////////////////////////////
// Get a vector with the entires of the desired column
//////////////////////////////////////////////////////////////////////
void MATRIX::setColumn(const VECTOR& column, int index)
{
  assert(column.size() == _rows);

  for (int x = 0; x < _rows; x++)
    (*this)(x, index) = column[x];
}

//////////////////////////////////////////////////////////////////////
// add the current vector to the desired column
//////////////////////////////////////////////////////////////////////
void MATRIX::addColumn(const VECTOR& column, int index)
{
  assert(column.size() == _rows);

  for (int x = 0; x < _rows; x++)
    (*this)(x, index) += column[x];
}

//////////////////////////////////////////////////////////////////////
// write the vector to a Matlab file
//////////////////////////////////////////////////////////////////////
void MATRIX::writeMatlab(string filename, string varname)
{
  FILE* file;
  file = fopen(filename.c_str(), "w");
  fprintf(file, "%s = [", varname.c_str());
  for (int x = 0; x < _rows; x++)
  {
    for (int y = 0; y < _cols; y++)
      fprintf(file, "%f ", (*this)(x,y));
    fprintf(file, "; ");
  }
  fprintf(file, "];\n");

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Get a submatrix of the current matrix
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::getSubmatrix(int rowBegin, int totalRows, int colBegin, int totalCols)
{
  assert(rowBegin < _rows);
  assert(colBegin < _cols);
  assert(rowBegin >= 0);
  assert(colBegin >= 0);
  assert(totalRows > 0);
  assert(totalCols > 0);
  assert(rowBegin + totalRows <= _rows);
  assert(colBegin + totalCols <= _rows);

  MATRIX toReturn(totalRows, totalCols);

  for (int row = 0; row < totalRows; row++)
    for (int col = 0; col < totalCols; col++)
    {
      int thisRow = row + rowBegin;
      int thisCol = col + colBegin;

      toReturn(row, col) = (*this)(thisRow, thisCol);
    }

  return toReturn;
}

//////////////////////////////////////////////////////////////////////
// Compute the QR factorization of the current matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::qr()
{
  // repack into fortran format for QR call
  formatFortran();

  int tauSize = (_rows < _cols) ? _rows : _cols;
  int lwork = _cols;
  Real* work = new Real[lwork];
  Real* tau = new Real[tauSize];

  // compute the QR
  int info;
#ifdef SINGLE_PRECISION
  sgeqrf(&_rows, &_cols, _matrix, &_rows, tau, work, &lwork, &info);
#else
  dgeqrf(&_rows, &_cols, _matrix, &_rows, tau, work, &lwork, &info);
#endif
  assert(info >= 0);

  // form Q
#ifdef SINGLE_PRECISION
  sorgqr(&_rows, &_cols, &_cols, _matrix, &_rows, tau, work, &lwork, &info);
#else
  dorgqr(&_rows, &_cols, &_cols, _matrix, &_rows, tau, work, &lwork, &info);
#endif
  assert(info >= 0);

  delete[] work;
  delete[] tau;

  // unpack fortran format results
  unformatFortran();
}

//////////////////////////////////////////////////////////////////////
// repack matrix into fortran ordering in preparation for a LAPACK call
// that does not support transposes
//////////////////////////////////////////////////////////////////////
void MATRIX::formatFortran()
{
  Real* repacked = new Real[_rows * _cols];

  int i = 0;
  for (int x = 0; x < _cols; x++)
    for (int y = 0; y < _rows; y++, i++)
      repacked[i] = (*this)(y,x);
  
  delete[] _matrix;
  _matrix = repacked;
}

//////////////////////////////////////////////////////////////////////
// unpack matrix from a LAPACK call that only supports fortran ordering
//////////////////////////////////////////////////////////////////////
void MATRIX::unformatFortran()
{
  Real* unpacked = new Real[_rows * _cols];

  int i = 0;
  for (int x = 0; x < _cols; x++)
    for (int y = 0; y < _rows; y++, i++)
      //repacked[i] = (*this)(y,x);
      unpacked[y * _cols + x] = _matrix[i];
  
  delete[] _matrix;
  _matrix = unpacked;
}

//////////////////////////////////////////////////////////////////////
// Are any of the entries a nan?
//////////////////////////////////////////////////////////////////////
bool MATRIX::anyNans()
{
  for (int x = 0; x < _rows * _cols; x++)
#ifdef _WIN32
    if (_isnan(_matrix[x]))
#else
    if (isnan(_matrix[x]))
#endif
      return true;
  return false;
}

//////////////////////////////////////////////////////////////////////
// scale all the columns by the entries in the vector
//////////////////////////////////////////////////////////////////////
void MATRIX::scaleColumns(VECTOR scale)
{
  for (int x = 0; x < _cols; x++)
    for (int y = 0; y < _rows; y++)
      (*this)(y,x) *= scale[x];
}

//////////////////////////////////////////////////////////////////////
// copy upper triangle into lower
//////////////////////////////////////////////////////////////////////
void MATRIX::copyUpperToLower()
{
  assert(_cols == _rows);

  for (int x = 0; x < _rows; x++)
    for (int y = x + 1; y < _cols; y++)
      (*this)(y,x) = (*this)(x,y);
}

//////////////////////////////////////////////////////////////////////
// compute the Moore-Penrose pseudo-inverse
//////////////////////////////////////////////////////////////////////
void MATRIX::pseudoInverse(MATRIX& inverse)
{
  MATRIX M(_cols, _rows);
  for (int y = 0; y < _cols; y++)
    for (int x = 0; x < _rows; x++)
      M(y,x) = (*this)(x,y);

  int bigger = (_rows > _cols) ? _rows : _cols;
  int smaller = (_rows > _cols) ? _cols: _rows;

  // resize everything to the proper size
  MATRIX U(_rows, _rows);     U.clear();
  MATRIX VT(_cols, _cols);    VT.clear();
  VECTOR S(bigger);           S.clear();
  inverse.resizeAndWipe(_cols, _rows);

  // The pseudo inverse matrix of S diagonal matrix
  MATRIX SM(_cols, _rows);    SM.clear();

  // compute all vectors
  char jobu = 'S';
  char jobvt = 'S';
  int ldu = _rows;
  int ldvt = _cols;
  
  // allocate the work array
  int lwork = 3 * smaller + bigger;
  lwork = (lwork < 5 * smaller) ? 5 * smaller : lwork;
  Real* work = new Real[lwork];
 
  int info;
#ifdef SINGLE_PRECISION
  //sgesvd_(&jobu, &jobvt, (__CLPK_integer*)&_rows, (__CLPK_integer*)&_cols, M.data(), (__CLPK_integer*)&_rows, S.data(), U.data(), (__CLPK_integer*)&ldu, VT.data(), (__CLPK_integer*)&ldvt, work, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
  sgesvd_(&jobu, &jobvt, &_rows, &_cols, M.data(), &_rows, S.data(), U.data(), &ldu, VT.data(), &ldvt, work, &lwork, &info);
#else
  //dgesvd_(&jobu, &jobvt, (__CLPK_integer*)&_rows, (__CLPK_integer*)&_cols, M.data(), (__CLPK_integer*)&_rows, S.data(), U.data(), (__CLPK_integer*)&ldu, VT.data(), (__CLPK_integer*)&ldvt, work, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
  dgesvd_(&jobu, &jobvt, &_rows, &_cols, M.data(), &_rows, S.data(), U.data(), &ldu, VT.data(), &ldvt, work, &lwork, &info);
#endif

  if (info != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "SVD DID NOT CONVERGE WHILE COMPUTING PSEUDO-INVERSE" << endl;
  }
  
  // Don't transpose these, so U is really UT, VT is really V.
  // U = U.transpose();
  // VT = VT.transpose();
  
  delete[] work;

  Real threshold = 0.5 * sqrt((_cols + _rows) + 1.0) * S(0) * REAL_EPSILON;
  for (int x = 0; x < S.size(); x++)
    SM(x,x) = (S(x) > threshold) ? 1.0 / S(x) : 0.0;

  inverse = (VT * SM) * U;
}

//////////////////////////////////////////////////////////////////////
// get the vector of sums of the rows in the matrix
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::rowSums()
{
  VECTOR sums(_rows);
  for (int x = 0; x < _rows; x++)
    for (int y = 0; y < _cols; y++)
      sums[x] += (*this)(x,y);

  return sums;
}

//////////////////////////////////////////////////////////////////////
// get the cross product matrix of a vector
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::cross(const VEC3F& vec)
{
  MATRIX final(3,3);
  final(0,0) = 0.0;
  final(1,0) = vec[2];
  final(2,0) = -vec[1];

  final(0,1) = -vec[2];
  final(1,1) = 0.0;
  final(2,1) = vec[0];

  final(0,2) = vec[1];
  final(1,2) = -vec[0];
  final(2,2) = 0.0;

  return final;
}

//////////////////////////////////////////////////////////////////////
// get the cross product matrix of a vector
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::cross(const VECTOR& vec)
{
  MATRIX final(3,3);
  final(0,0) = 0.0;
  final(1,0) = vec[2];
  final(2,0) = -vec[1];

  final(0,1) = -vec[2];
  final(1,1) = 0.0;
  final(2,1) = vec[0];

  final(0,2) = vec[1];
  final(1,2) = -vec[0];
  final(2,2) = 0.0;

  return final;
}

//////////////////////////////////////////////////////////////////////
// Get the Euclidean norm of a matrix
//////////////////////////////////////////////////////////////////////
Real MATRIX::norm2() const
{
  MATRIX copy(*this);
  MATRIX U;
  VECTOR S;
  MATRIX VT;
  copy.SVD(U,S,VT);

  return S.maxValue();
}

//////////////////////////////////////////////////////////////////////
// Static routines for handling array-based matrices.
// Use with case - nothing here does any bound checking
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Copy one matrix in to another (A <- B)
//////////////////////////////////////////////////////////////////////
void MATRIX::copy( Real *A, const Real *B, int rows, int cols )
{
  memcpy( A, B, rows * cols * sizeof( Real ) );
}

//////////////////////////////////////////////////////////////////////
// Copy a 3x3 matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::copy( Real *A, const MATRIX3 &B )
{
  for ( int i = 0; i < 3; i++ )
  for ( int j = 0; j < 3; j++ )
  {
    access( A, 3, 3, i, j ) = B( i, j );
  }
}

//////////////////////////////////////////////////////////////////////
// Add one matrix to another (A = beta * A + alpha * B
//////////////////////////////////////////////////////////////////////
void MATRIX::axpy( Real *A, const Real *B, int rows, int cols,
                  Real alpha, Real beta )
{
#ifdef SINGLE_PRECISION
	cblas_saxpy(rows * cols, alpha, B, beta, A, 1);
#else
	cblas_daxpy(rows * cols, alpha, B, beta, A, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-matrix multiplication (C = beta * C + alpha * A * B)
//////////////////////////////////////////////////////////////////////
void MATRIX::gemm( const Real *A, const Real *B, Real *C,
                  int rowsA, int colsA, int rowsB, int colsB,
                  bool transposeA, bool transposeB,
                  Real alpha, Real beta )
{
  int requiredCols = transposeB ? rowsB : colsB;

#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			transposeA ? CblasTrans : CblasNoTrans,
			transposeB ? CblasTrans : CblasNoTrans,
			rowsA, 
      colsB,
      colsA, 
			alpha,
			A, colsA,
			B, colsB,
			beta,
			C, requiredCols);
#else
	cblas_dgemm(
			CblasRowMajor,
			transposeA ? CblasTrans : CblasNoTrans,
			transposeB ? CblasTrans : CblasNoTrans,
			rowsA, 
      colsB,
      colsA, 
			alpha,
			A, colsA,
			B, colsB,
			beta,
			C, requiredCols);
#endif
}

//////////////////////////////////////////////////////////////////////
// Get transpose (B = A^T)
//////////////////////////////////////////////////////////////////////
void MATRIX::transpose( Real *A, const Real *B, int rows, int cols )
{
  for (int y = 0; y < cols; y++)
    for (int x = 0; x < rows; x++)
      access( A, cols, rows, y, x ) = access( B, rows, cols, x, y );
}

//////////////////////////////////////////////////////////////////////
// In place transpose for square matrices
//////////////////////////////////////////////////////////////////////
void MATRIX::transpose( Real *A, int rows )
{
  for ( int i = 1; i < rows; i++ )
  {
    for ( int j = 0; j < i; j++ )
    {
      Real temp = access( A, rows, rows, i, j );
      access( A, rows, rows, i, j ) = access( A, rows, rows, j, i );
      access( A, rows, rows, j, i ) = temp;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Compute eigenvalues/vectors.  We require everything here,
// including workspaces, etc.
// workspace should have size (7 * rows)
// vectorWorkspace should have size (rows * rows)
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem( Real *A, int rows,
                          Real *eigenvalues, Real *eigenvectors,
                          Real *workspace, Real *vectorWorkspace )
{
  int rowsize = rows;
  int worksize = 5 * rows;

  Real *work = workspace;
  Real *valuesReal = workspace + worksize;
  Real *valuesImag = valuesReal + rows;
  Real *vectors = vectorWorkspace;

  // the actual LAPACK call
  int error;
#ifdef SINGLE_PRECISION
  sgeev("V","N", &rowsize, A, &rowsize, 
        valuesReal, valuesImag, 
        vectors, &rowsize, NULL, &rowsize,
        work, &worksize, &error);
#else
  dgeev("V","N", &rowsize, A, &rowsize, 
        valuesReal, valuesImag, 
        vectors, &rowsize, NULL, &rowsize,
        work, &worksize, &error);
#endif

  if (error != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " eigenvalue solver bombed!" << endl;
  }

  // copy out results
  for (int x = 0; x < rows; x++)
    eigenvalues[x] = valuesReal[x];

  for (int x = 0; x < rows; x++)
    for (int y = 0; y < rows; y++)
      MATRIX::access( eigenvectors, rows, rows, x, y ) = vectors[x + y * rows];
}

//////////////////////////////////////////////////////////////////////
// 2x2 eigensolver
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem2x2( Real *matrix, Real *eigenvalues, Real *eigenvectors )
{
  Real A = access( matrix, 2, 2, 0, 0);
  Real B = access( matrix, 2, 2, 0, 1);
  Real C = access( matrix, 2, 2, 1, 0);
  Real D = access( matrix, 2, 2, 1, 1);

  if (B * C <= 0.1e-20) {
    eigenvalues[0] = A; 
    access( eigenvectors, 2, 2, 0,0) = 1; 
    access( eigenvectors, 2, 2, 1,0) = 0;
    eigenvalues[1] = D; 
    access( eigenvectors, 2, 2, 0,1) = 0; 
    access( eigenvectors, 2, 2, 1,1) = 1;
    return;
  }

  Real tr = A + D;
  Real det = A * D - B * C;
  Real S = sqrt(tr * tr * 0.25 - det );
  eigenvalues[0] = tr * 0.5 + S;
  eigenvalues[1] = tr * 0.5 - S;

  Real temp = (A - D) * (A - D) * 0.25 + B * C;
  Real SS = (temp < 0.0) ? 0.0 : sqrt(temp);
  if (A - D < 0.0) {
    access( eigenvectors, 2, 2, 0,0) = C;
    access( eigenvectors, 2, 2, 1,0) = -(A - D) * 0.5 + SS;
    access( eigenvectors, 2, 2, 0,1) = (A - D) * 0.5 - SS;
    access( eigenvectors, 2, 2, 1,1) = B;
  } 
  else {
    access( eigenvectors, 2, 2, 0,1) = C;
    access( eigenvectors, 2, 2, 1,1) = -(A - D) * 0.5 - SS;
    access( eigenvectors, 2, 2, 0,0) = (A - D) * 0.5 + SS;
    access( eigenvectors, 2, 2, 1,0) = B;
  }

  Real n1 = sqrt(access( eigenvectors, 2, 2, 0,0) * access( eigenvectors, 2, 2, 0,0) +
                 access( eigenvectors, 2, 2, 1,0) * access( eigenvectors, 2, 2, 1,0));
  Real inv = 1.0 / n1;
  access( eigenvectors, 2, 2, 0,0) *= inv; 
  access( eigenvectors, 2, 2, 1,0) *= inv;
  Real n2 = sqrt(access( eigenvectors, 2, 2, 0,1) * access( eigenvectors, 2, 2, 0,1) +
                 access( eigenvectors, 2, 2, 1,1) * access( eigenvectors, 2, 2, 1,1));
  inv = 1.0 / n2;
  access( eigenvectors, 2, 2, 0,1) *= inv; 
  access( eigenvectors, 2, 2, 1,1) *= inv;
}

/*
//////////////////////////////////////////////////////////////////////
// 3x3 low precision eigensolver
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem3x3( Real *A, Real *eigenvalues, Real *eigenvectors )
{
	register int i, j, k, l;

  //float TOL = 1e-8f;
  //int MAX_SWEEPS = 500;
  float TOL = 1e-3f;
  int MAX_SWEEPS = 5;
  unsigned int n = 3;

  float a[9];
  float d[3];
  float v[9];
  i = 0;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++, i++)
    {
      a[i] = access( A, 3, 3, x, y);
      v[i] = (x == y) ? 1.0 : 0.0;
    }
  
  
	float onorm, dnorm;
	float b, dma, q, t, c, s;
	float atemp, vtemp, dtemp;

	// Set v to the identity matrix, set the vector d to contain the
	// diagonal elements of the matrix a
	d[0] = a[0];
	d[1] = a[4];
	d[2] = a[8];

	for (l = 1; l <= MAX_SWEEPS; l++)
	{
		// Set dnorm to be the maximum norm of the diagonal elements, set
		// onorm to the maximum norm of the off-diagonal elements
		
		dnorm = (float)fabs(d[0]) + (float)fabs(d[1]) + (float)fabs(d[2]);
		onorm = (float)fabs(a[1]) + (float)fabs(a[2]) + (float)fabs(a[5]);
		// Normal end point of this algorithm.
		if((onorm/dnorm) <= TOL)
			goto Exit_now;

		for (j = 1; j < static_cast<int>(n); j++)
		{
			for (i = 0; i <= j - 1; i++)
			{

				b = a[n*i+j];
				if(fabs(b) > 0.0f)
				{
					dma = d[j] - d[i];
					if((fabs(dma) + fabs(b)) <= fabs(dma))
						t = b / dma;
					else
					{
						q = 0.5f * dma / b;
						t = 1.0f/((float)fabs(q) + (float)sqrt(1.0f+q*q));
						if (q < 0.0)
							t = -t;
					}

					c = 1.0f/(float)sqrt(t*t + 1.0f);
					s = t * c;
					a[n*i+j] = 0.0f;

					for (k = 0; k <= i-1; k++)
					{
						atemp = c * a[n*k+i] - s * a[n*k+j];
						a[n*k+j] = s * a[n*k+i] + c * a[n*k+j];
						a[n*k+i] = atemp;
					}

					for (k = i+1; k <= j-1; k++)
					{
						atemp = c * a[n*i+k] - s * a[n*k+j];
						a[n*k+j] = s * a[n*i+k] + c * a[n*k+j];
						a[n*i+k] = atemp;
					}

					for (k = j+1; k < static_cast<int>(n); k++)
					{
						atemp = c * a[n*i+k] - s * a[n*j+k];
						a[n*j+k] = s * a[n*i+k] + c * a[n*j+k];
						a[n*i+k] = atemp;
					}

					for (k = 0; k < static_cast<int>(n); k++)
					{
						vtemp = c * v[n*k+i] - s * v[n*k+j];
						v[n*k+j] = s * v[n*k+i] + c * v[n*k+j];
						v[n*k+i] = vtemp;
					}

					dtemp = c*c*d[i] + s*s*d[j] - 2.0f*c*s*b;
					d[j] = s*s*d[i] + c*c*d[j] + 2.0f*c*s*b;
					d[i] = dtemp;
				}
			}
		}
	}

Exit_now:
  for (int x = 0; x < 3; x++)
    eigenvalues[x] = d[x];

  i = 0;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++, i++)
      access( eigenvectors, 3, 3, x, y) = v[i];

	return;
}
*/

//////////////////////////////////////////////////////////////////////
// Scales a matrix.  A *= alpha
//////////////////////////////////////////////////////////////////////
void MATRIX::scale( Real *A, int rows, int cols, Real alpha )
{
#ifdef SINGLE_PRECISION
  cblas_sscal(cols * rows, alpha, A, 1);
#else
  cblas_dscal(cols * rows, alpha, A, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// take the absolute value of all entries
//////////////////////////////////////////////////////////////////////
void MATRIX::absoluteValue()
{
  // what's the mac version of cblas_dabs1?
  for (int x = 0; x < _rows * _cols; x++)
    _matrix[x] = fabs(_matrix[x]);
}
