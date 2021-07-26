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
// MATRIX_DEBUG.cpp: 
//
// This is a *slow* implementation of the MATRIX.h header. It does not
// use ATLAS, LAPACK, BLAS, or MKL to facilitate debugging. As it also
// does not need to link to those libraries, it is also self-contained
// and much easier to compile.
//
// Note that the "solve" function is *not* implemented however, since
// I can't seem to find slow, simple, reliable LU code anywhere.
//
// Support for tet inversion makes the inclusion of PetSc's eigenvalue
// solver necessary unfortunately, which affects the self-containment
// of the code. Suggestions for alternatives are welcome.
//
//////////////////////////////////////////////////////////////////////

#include "MATRIX.h"
//#include <petscblaslapack.h>
//#include <slepceps.h>

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
  for (int x = 0; x < _rows * _cols; x++)
    _matrix[x] = data[x];
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
  for (int x = 0; x < _rows * _cols; x++)
    _matrix[x] = m._matrix[x];
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

MATRIX::MATRIX(MATRIX3& matrix3) :
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

//////////////////////////////////////////////////////////////////////
// construct a homogeneous transformation matrix
//////////////////////////////////////////////////////////////////////
MATRIX::MATRIX(MATRIX3& rotation, VEC3F& translation)
{
  cout << " UNIMPLEMENTED" << endl;
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
    delete[] _matrix;
  _rows = rows;
  _cols = cols;

  _matrix = new Real[_rows * _cols];
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
void MATRIX::read(const char* filename)
{
  FILE* file;
  file = fopen(filename, "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << filename << " not found! " << endl;
    return;
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
}

//////////////////////////////////////////////////////////////////////
// read matrix from a file
//////////////////////////////////////////////////////////////////////
void MATRIX::read(FILE* file)
{
  // read dimensions
  fread((void*)&_rows, sizeof(int), 1, file);
  fread((void*)&_cols, sizeof(int), 1, file);

  // read data
  delete[] _matrix;
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
MATRIX MATRIX::transpose()
{
  MATRIX toReturn(_cols, _rows);

  for (int y = 0; y < _cols; y++)
    for (int x = 0; x < _rows; x++)
      toReturn(y,x) = (*this)(x,y);

  return toReturn;
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator*(MATRIX& A, VECTOR& x) 
{
  VECTOR y(A.rows());

#if BOUNDS_CHECKING_ENABLED
  if (A.cols() != x.size())
    cout << __FILE__ << " " << __LINE__ << " : Matrix-Vector dim mismatch! " << endl;
#endif

  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      y(i) += x(j) * A(i, j);

  return y;
}

//////////////////////////////////////////////////////////////////////
// Vector-matrix multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator*(VECTOR& x, MATRIX& A)
{
  VECTOR y(A.cols());
  for (int i = 0; i < A.cols(); i++)
    for (int j = 0; j < A.rows(); j++)
      y[i] += A(j, i) * x(j);
  return y;
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(MATRIX& A, Real alpha) 
{
  MATRIX y(A.rows(), A.cols());

  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      y(i,j) = A(i, j) * alpha;

  return y;
}

//////////////////////////////////////////////////////////////////////
// Scale matrix
//////////////////////////////////////////////////////////////////////
MATRIX operator*(Real alpha, MATRIX& A) 
{
  MATRIX y(A.rows(), A.cols());

  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      y(i,j) = A(i, j) * alpha;

  return y;
}

//////////////////////////////////////////////////////////////////////
// Matrix-Matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(MATRIX& A, MATRIX& B) 
{
  MATRIX y(A.rows(), B.cols());

#if BOUNDS_CHECKING_ENABLED
  if (A.cols() != B.rows())
  {
    cout << __FILE__ << " " << __LINE__ << " : Matrix-Matrix dimensions do not match! " << endl;
    return y;
  }
#endif

  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < B.cols(); j++)
      for (int k = 0; k < A.cols(); k++)
        y(i,j) += A(i, k) * B(k, j);

  return y;
}

//////////////////////////////////////////////////////////////////////
// scale vector by a constant
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator*=(const Real& alpha)
{
  for (int x = 0; x < _cols * _rows; x++)
    _matrix[x] *= alpha;
  
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply where A is transposed
//////////////////////////////////////////////////////////////////////
VECTOR operator^(MATRIX& A, VECTOR& x) 
{
  VECTOR y(A.cols());

#if BOUNDS_CHECKING_ENABLED
  if (A.rows() != x.size())
  {
    cout << __FILE__ << " " << __LINE__ << " : Transposed Matrix-Vector dim mismatch! " << endl;
    cout << "Matrix: " << A.rows() << " " << A.cols() << endl;
    cout << "Vector: " << x.size() << endl;
  }
#endif

  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      y(j) += x(i) * A(i, j);

  return y;
}

//////////////////////////////////////////////////////////////////////
// Matrix^T -Matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator^(MATRIX& A, MATRIX& B) 
{
  MATRIX y(A.cols(), B.cols());

#if BOUNDS_CHECKING_ENABLED
  if (A.rows() != B.rows())
  {
    cout << __FILE__ << " " << __LINE__ << " : Transposed Matrix-Matrix dimensions do not match! " << endl;
    return y;
  }
#endif

  for (int i = 0; i < A.cols(); i++)
    for (int j = 0; j < B.cols(); j++)
      for (int k = 0; k < A.rows(); k++)
        y(i,j) += A(k, i) * B(k, j);

  return y;
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, MATRIX& matrix)
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
  for (int x = 0; x < _rows * _cols; x++)
    _matrix[x] = m._matrix[x];

  return *this;
}

//////////////////////////////////////////////////////////////////////
// self minus
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator-=(const MATRIX& m)
{
  if (m._cols != _cols || m._rows != _rows)
  {
    delete[] _matrix;
    _cols = m._cols;
    _rows = m._rows;

    _matrix = new Real[_rows * _cols];
  }
  for (int x = 0; x < _rows * _cols; x++)
    _matrix[x] -= m._matrix[x];

  return *this;
}

//////////////////////////////////////////////////////////////////////
// self plus
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator+=(const MATRIX& m)
{
  if (m._cols != _cols || m._rows != _rows)
  {
    delete[] _matrix;
    _cols = m._cols;
    _rows = m._rows;

    _matrix = new Real[_rows * _cols];
  }
  for (int x = 0; x < _rows * _cols; x++)
    _matrix[x] += m._matrix[x];

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
  int totalSize = matrix.rows() * matrix.cols();
  int index = row * _cols;

  for (int x = 0; x < totalSize; x++, index++)
    _matrix[index] = matrix._matrix[x];
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B += alpha * A, where B is this matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::axpy(Real alpha, MATRIX& A)
{
#if BOUNDS_CHECKING_ENABLED
  if (A.rows() != _rows || A.cols() != _cols)
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : Matrix sizes do not match! " << endl;
#endif

  int total = _rows * _cols;
  for (int x = 0; x < total; x++)
    _matrix[x] += alpha * A._matrix[x];
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B = alpha * A, where B is this matrix, and 
// current contents of B are stomped
//////////////////////////////////////////////////////////////////////
void MATRIX::clearingAxpy(Real alpha, MATRIX& A)
{
#if BOUNDS_CHECKING_ENABLED
  if (A.rows() != _rows || A.cols() != _cols)
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : Matrix sizes do not match! " << endl;
#endif

  int total = _rows * _cols;
  for (int x = 0; x < total; x++)
    _matrix[x] = alpha * A._matrix[x];
}

//////////////////////////////////////////////////////////////////////
// BLAS gemm operation: C += alpha * A * B where C is this matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::gemm(Real alpha, MATRIX& A, MATRIX& B)
{
#if BOUNDS_CHECKING_ENABLED
  if (A.cols() != B.rows() || A.rows() != _rows || B.cols() != _cols)
  {
    cout << __FILE__ << " " << __LINE__ << " : Matrix-Matrix dimensions do not match! " << endl;
    cout << "Matrix A: " << A.cols() << " " << A.cols() << endl;
    cout << "Matrix B: " << B.cols() << " " << B.cols() << endl;
    return;
  }
#endif

  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < B.cols(); j++)
    {
      Real sum = 0.0;
      for (int k = 0; k < A.cols(); k++)
        sum += A(i, k) * B(k, j);

      (*this)(i,j) += alpha * sum;
    }
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
#if BOUNDS_CHECKING_ENABLED
  if (_cols != 3)
    cout << __FILE__ << " " << __LINE__ << " : Matrix-Vector dimensions do not match!" << endl;
#endif

  VECTOR y(_rows);
  for (int j = 0; j < _cols; j++)
    for (int i = 0; i < _rows; i++)
      y(i) += (*this)(i,j) * x[i];

  return y;
}

//////////////////////////////////////////////////////////////////////
// BLAS gemv operation: y += alpha * A * x where A is this matrix
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::gemv(Real alpha, VEC3F& x)
{
#if BOUNDS_CHECKING_ENABLED
  if (_cols != 3)
    cout << __FILE__ << " " << __LINE__ << " : Matrix-Vector dimensions do not match!" << endl;
#endif

  VECTOR y(_rows);
  for (int j = 0; j < _cols; j++)
  {
    for (int i = 0; i < _rows; i++)
      y(i) += (*this)(i,j) * x[i];
    y(j) *= alpha;
  }

  return y;
}

//////////////////////////////////////////////////////////////////////
// BLAS gemm operation: C = alpha * A * B where C is this matrix and
// current contents of C are stomped
//////////////////////////////////////////////////////////////////////
void MATRIX::clearingGemm(Real alpha, MATRIX& A, MATRIX& B)
{
#if BOUNDS_CHECKING_ENABLED
  if (A.cols() != B.rows() || A.rows() != _rows || B.cols() != _cols)
  {
    cout << __FILE__ << " " << __LINE__ << " : Matrix-Matrix dimensions do not match! " << endl;
    cout << "Matrix A: " << A.rows() << " " << A.cols() << endl;
    cout << "Matrix B: " << B.rows() << " " << B.cols() << endl;
    cout << "Matrix C: " << this->rows() << " " << this->cols() << endl;
    return;
  }
#endif

  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < B.cols(); j++)
    {
      Real sum = 0.0;
      for (int k = 0; k < A.cols(); k++)
        sum += A(i, k) * B(k, j);

      (*this)(i,j) = alpha * sum;
    }
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//////////////////////////////////////////////////////////////////////
void MATRIX::multiplyInplace(VECTOR& x, VECTOR& y) 
{
#if BOUNDS_CHECKING_ENABLED
  if (_cols != x.size())
    cout << __FILE__ << " " << __LINE__ << " : Matrix-Vector dim mismatch! " << endl;
#endif

  // do the product into another vector in case x is also y
  VECTOR z(y.size());
  z.clear();
  for (int i = 0; i < _rows; i++)
    for (int j = 0; j < _cols; j++)
      z(i) += x(j) * (*this)(i, j);
  y = z;
}

//////////////////////////////////////////////////////////////////////
// solve the linear system Ax = b, return x in the passed in b
//////////////////////////////////////////////////////////////////////
void MATRIX::solve(VECTOR& b)
{
  cout << __FILE__ << " " << __LINE__ << " : MATRIX_DEBUG.cpp is being used, so" 
       << " there is no LU solver. Use MATRIX_FAST.cpp." << endl;
}

//////////////////////////////////////////////////////////////////////
// solve for the eigensystem of the matrix
//////////////////////////////////////////////////////////////////////
//#include "petsc.h"
//#include "petscblaslapack.h"
void MATRIX::eigensystem(VECTOR& eigenvalues, MATRIX& eigenvectors)
{
  /*
  // This is commented out to make the code more easily compiled in Cygwin.
  // If you really want to 
  // basic error checking
  if (_rows != _cols)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Matrix must be square to get eigenvalues! " << endl;
    return;
  }

  // resize result space
  eigenvalues.resizeAndWipe(_rows);
  eigenvectors.resizeAndWipe(_rows, _rows);

  PetscBLASInt rowsize = _rows;
  PetscBLASInt worksize = 5 * _rows;
  PetscScalar* work;
  PetscReal*   valuesReal;
  PetscReal*   valuesImag;
  PetscScalar* input;
  PetscScalar* vectors;

  // allocate space
  PetscMalloc(2 * _rows * sizeof(PetscReal),&valuesReal);
  PetscMalloc(worksize * sizeof(PetscReal),&work);
  PetscMalloc(_rows * _rows * sizeof(PetscReal),&input);
  PetscMalloc(_rows * _rows * sizeof(PetscScalar),&vectors);
  valuesImag = valuesReal + _rows;

  // copy matrix into PetSc array
  for (int x = 0; x < _rows * _rows; x++)
    input[x] = _matrix[x];
 
  // the actual LAPACK call
  PetscBLASInt error;
  LAPACKgeev_("V","N", &rowsize, input, &rowsize, 
              valuesReal, valuesImag, 
              vectors, &rowsize, NULL, &rowsize,
              work, &worksize, &error);

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
  PetscFree(input);
  PetscFree(valuesReal);
  PetscFree(work);
  PetscFree(vectors);
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " Not supported in this implementation! Install Intel MKL!" << endl;
  */
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
// copy this matrix to MATRIX3 type
//////////////////////////////////////////////////////////////////////
void MATRIX::copiesInto(MATRIX3& matrix3)
{
  if (_rows != 3 || _cols != 3)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "Trying to copy a MATRIX that is not 3x3 to a MATRIX3!" << endl;
    return;
  }
  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
      matrix3(x,y) = (*this)(x,y);
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
// 2-norm of the whole matrix
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
// difference of two matrices
//////////////////////////////////////////////////////////////////////
MATRIX operator-(MATRIX& A, MATRIX& B)
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
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Not supported in this implementation! Install Intel MKL!" << endl;

  return false;
}

//////////////////////////////////////////////////////////////////////
// Use LU factorization to solve for rhs b
//////////////////////////////////////////////////////////////////////
void MATRIX::solveLU(VECTOR& b)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Not supported in this implementation! Install Intel MKL!" << endl;
}
//////////////////////////////////////////////////////////////////////
// Get Cholesky factorization of the matrix
//////////////////////////////////////////////////////////////////////
bool MATRIX::factorCholesky()
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Not supported in this implementation! Install Intel MKL!" << endl;
  return false;
}
//////////////////////////////////////////////////////////////////////
// Solve matrix using precomputed Cholesky factorization
//////////////////////////////////////////////////////////////////////
void MATRIX::solveCholesky(VECTOR& b)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Not supported in this implementation! Install Intel MKL!" << endl;
}
//////////////////////////////////////////////////////////////////////
// Get the matrix SVD
//////////////////////////////////////////////////////////////////////
void MATRIX::SVD(MATRIX& U, VECTOR& S, MATRIX& VT)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " Not supported in this implementation! Install Intel MKL!" << endl;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
Real MATRIX::norm1()
{
  cout << " UNIMPLEMENTED. " << endl;
  return 0.0;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
Real MATRIX::normInf()
{
  cout << " UNIMPLEMENTED. " << endl;
  return 0.0;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
Real MATRIX::maxAbsEntry()
{
  cout << " UNIMPLEMENTED. " << endl;
  return 0.0;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void MATRIX::orthogonalize()
{
  cout << " UNIMPLEMENTED. " << endl;
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

/*
//////////////////////////////////////////////////////////////////////
// quickly solve for the eigensystem of a 3x3 matrix
// This code courtesy of Alec Rivers -- where did he get it?
//
// The Irving et al. diagonalization needs more digits than this in
// practice for the Newton solver to converge though, so I am not
// using it. For a fixed number of iterations (3) it is about twice
// as fast as MKL LAPACK.
// 
// This static function computes the eigenvalues and
// eigenvectors of a SYMMETRIC nxn matrix. This method is used
// internally by OpenBabel, but may be useful as a general
// eigenvalue finder.
// 
// The algorithm uses Jacobi transformations. It is described
// e.g. in Wilkinson, Reinsch "Handbook for automatic computation,
// Volume II: Linear Algebra", part II, contribution II/1. The
// implementation is also similar to the implementation in this
// book. This method is adequate to solve the eigenproblem for
// small matrices, of size perhaps up to 10x10. For bigger
// problems, you might want to resort to the sophisticated routines
// of LAPACK.
// 
// \note If you plan to find the eigenvalues of a symmetric 3x3
// matrix, you will probably prefer to use the more convenient
// method findEigenvectorsIfSymmetric()
// 
// @param n the size of the matrix that should be diagonalized
// 
// @param a array of size n^2 which holds the symmetric matrix
// whose eigenvectors are to be computed. The convention is that
// the entry in row r and column c is addressed as a[n*r+c] where,
// of course, 0 <= r < n and 0 <= c < n. There is no check that the
// matrix is actually symmetric. If it is not, the behaviour of
// this function is undefined. On return, the matrix is
// overwritten with junk.
// 
// @param d pointer to a field of at least n doubles which will be
// overwritten. On return of this function, the entries d[0]..d[n-1]
// will contain the eigenvalues of the matrix.
// 
// @param v an array of size n^2 where the eigenvectors will be
// stored. On return, the columns of this matrix will contain the
// eigenvectors. The eigenvectors are normalized and mutually
// orthogonal.
//
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem3x3(VECTOR& eigenvalues, MATRIX& eigenvectors)
{
	register int i, j, k, l;

  float TOL = 1e-8f;
  int MAX_SWEEPS = 500;
  //float TOL = 1e-3f;
  //int MAX_SWEEPS = 50;
  unsigned int n = 3;

  float a[9];
  float d[3];
  float v[9];
  i = 0;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++, i++)
    {
      a[i] = (*this)(x,y);
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
      eigenvectors(x,y) = v[i];

	return;
}
*/

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
// sum of two matrices
//////////////////////////////////////////////////////////////////////
MATRIX operator+(MATRIX& A, MATRIX& B)
{
  MATRIX result(A.rows(), A.cols());
  for (int y = 0; y < A.cols(); y++)
    for (int x = 0; x < A.rows(); x++)
      result(x,y) = A(x,y) + B(x,y);
  return result;
}


