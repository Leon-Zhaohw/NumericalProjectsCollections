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

#ifndef MATRIX_H
#define MATRIX_H

#include <SETTINGS.h>
#include <map>
#include <iostream>
#include <cstdio>
#include <vector>
#include <VECTOR.h>
#include <VEC3.h>
#include <MATRIX3.h>

#include <memory.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// An arbitrary dimension matrix class
//////////////////////////////////////////////////////////////////////
class MATRIX {

public:
  MATRIX();
  MATRIX(int rows, int cols);
  MATRIX(int rows, int cols, Real* data);
  MATRIX(const char* filename);
  MATRIX(const MATRIX& m);
  // matrix with "vec" along the diagonal
  MATRIX(VECTOR& vec);

  // matrix with "cols" as columns
  MATRIX(vector<VECTOR>& cols);

  // make a copy of a 3x3 matrix
  //MATRIX(MATRIX3& matrix3);
  MATRIX(MATRIX3 matrix3);

  // construct a homogeneous transformation matrix
  MATRIX(const MATRIX3& rotation, const VEC3F& translation);
  
  virtual ~MATRIX();

  inline Real& operator()(int row, int col) {
    return _matrix[row * _cols + col];
  };
 
  // const version
  inline Real operator()(int row, int col) const {
    return _matrix[row * _cols + col];
  };

  int& rows() { return _rows; };
  int& cols() { return _cols; };
  int rows() const { return _rows; };
  int cols() const { return _cols; };
  static int& coutPrecision() { return _precision; };
  static int& coutWidth() { return _width; };

  // return a pointer to the beginning of a row
  Real* row(int index) { return &_matrix[index * _cols]; };
  const Real *row(int index) const { return &_matrix[index * _cols]; };
  VECTOR getRow(int index) const;

  // wipe the whole matrix
  void clear();

  // write the matrix to a binary file
  // everything is always written as a double
  void write(const char* filename);
  void write(FILE* file);
  void writeMatlab(string filename, string varname);

  // read from a binary file
  // everything is always read in as a double, then
  // converted if necessary
  bool read(const char* filename);
  void read(FILE* file);
  bool readSubcols(const char* filename, int maxCols);

  // resize the matrix and wipe to zero
  void resizeAndWipe(int rows, int cols);

  // overload operators
  MATRIX& operator=(const MATRIX m);
  MATRIX& operator=(const MATRIX3 m);
  MATRIX& operator-=(const MATRIX& m);
  MATRIX& operator+=(const MATRIX& m);
  MATRIX& operator*=(const Real& alpha);

  // return the matrix diagonal
  VECTOR diagonal();

  // return the transpose of the current matrix
  MATRIX transpose() const;

  void transpose( MATRIX &output );

  // raw data pointer
  Real* data() { return _matrix; };
  const Real* dataConst() const { return _matrix; };

  // stomp the current matrix with the given matrix 
  // starting at row number "row". It is your responsibility
  // to ensure that you don't fall off the end of this matrix.
  void setSubmatrix(MATRIX& matrix, int row);
  
  // stomp the current matrix with the given matrix 
  // starting at row number "row" and "col". It is your responsibility
  // to ensure that you don't fall off the end of this matrix.
  void setSubmatrix(MATRIX& matrix, int row, int col);

  // add the given matrix to the current matrix
  // starting at row number "row". It is your responsibility
  // to ensure that you don't fall off the end of this matrix.
  void addSubmatrix(MATRIX& matrix, int row);

  // subtract the given matrix from the current matrix
  // starting at row number "row". It is your responsibility
  // to ensure that you don't fall off the end of this matrix.
  void subtractSubmatrix(MATRIX& matrix, int row);

  // BLAS axpy operation: B += alpha * A, where B is this matrix
  //
  // Note that axpy actually applies to vectors, but in this
  // case we can just treat the matrix as a vector and multiply
  // all its elements by alpha
  void axpy(const Real alpha, const MATRIX& A);

  // same as axpy above, but matrix contents are stomped as well
  void clearingAxpy(Real alpha, MATRIX& A);

  // BLAS gemm operation: C += alpha * A * B
  // where C is this matrix
  void gemm(Real alpha, const MATRIX& A, const MATRIX& B, Real beta = 1.0,
            bool transposeA = false, bool transposeB = false);
  void gemm(const MATRIX& A, const MATRIX& B) { gemm(1.0f, A, B); };

  // compute A * B^T, where A is the current matrix
  MATRIX timesTranspose(MATRIX& B);

  // same as gemm above, but matrix contents are stomped as well
  void clearingGemm(Real alpha, MATRIX& A, MATRIX& B);
  void clearingGemm(MATRIX& A, MATRIX& B) { clearingGemm(1.0f, A, B); };

  // BLAS gemv operation: y = alpha * A * x
  // where A is this matrix
  VECTOR gemv(VEC3F& x);
  VECTOR gemv(Real alpha, VEC3F& x);
  //Untested -- don't uncomment until a test case comes up
  //VECTOR gemv(VECTOR& x);

  // BLAS gemv operation:
  //  y = alpha * A * x + beta * y
  //  No error checking in here - make sure y is the right size!
  void gemvInplace(Real alpha, const VECTOR &x, VECTOR &y, Real beta = 1.0,
                   bool transpose = false);

  // A bit of a hack
  void gemvInplace(Real alpha, VEC3F &x, VECTOR &y, Real beta = 1.0,
                   bool transpose = false);
  
  // solve the general linear system Ax = b, return x in the passed in b
  void solve(VECTOR& b);

  // solve the SPD linear system Ax = b, return x in the passed in b
  void solveSPD(VECTOR& b);

  // multiply in place
  // this * x = y
  // Do *NOT* let x == y!
  void multiplyInplace(VECTOR& x, VECTOR& y);

  // solve for eigenvalues
  void eigensystem(VECTOR& eigenvalues, MATRIX& eigenvectors);
  void eigensystem2x2(VECTOR& eigenvalues, MATRIX& eigenvectors);
  void eigensystem3x3(VECTOR& eigenvalues, MATRIX& eigenvectors);
  void eigenMinMax(Real& minEig, Real& maxEig);

  // Solve for eigenvalues in the symmetric case.
  // Behaviour is undefined for an unsymmetric matrix.
  // Puts all eigenvalues in the specified range (vl, vu]
  // in the provided vector and returns the number of
  // eigenvalues computed.
  // 
  // Destroys the matrix.
  //
  // This call expects a workspace.
  int eigensystemSymmetricRange(Real *eigenvalues,
                                Real *work, int lwork,
                                int *iwork, int liwork,
                                Real vl, Real vu,
                                Real tolerance = 1e-10 );

  inline int eigensystemSymmetricRange(
                                VECTOR &eigenvalues,
                                Real *work, int lwork,
                                int *iwork, int liwork,
                                Real vl, Real vu,
                                Real tolerance = 1e-10 )
  {
    return eigensystemSymmetricRange( eigenvalues.data(),
                                      work, lwork,
                                      iwork, liwork,
                                      vl, vu, tolerance );
  }

  // NOTE: this one doesn't destroy the matrix
  VECTOR eigenvalues();

  // solve SVD
  // Note this destroys the matrix!
  void SVD(MATRIX& U, VECTOR& S, MATRIX& VT);
 
  // perform PCA on this matrix
  // Note this destroys the matrix!
  void PCA(MATRIX& components, VECTOR& values);
  void meanNormalizedPCA(MATRIX& components, VECTOR& values);
  void verbosePCA(MATRIX& components, VECTOR& values, MATRIX& U, VECTOR& S, MATRIX& VT);

  // compute the Moore-Penrose pseudo-inverse
  void pseudoInverse(MATRIX& inverse);

  // subtract the mean from all the rows (in preparation for PCA)
  void subtractRowMeans();
  VECTOR rowMeans();
  VECTOR rowSums();

  // copy this matrix to MATRIX3 type
  void copiesInto(MATRIX3& matrix3);

  // Copy the given matrix (or matrix data) directly,
  // only reallocating if necessary.
  void copyInplace( const Real *data, int rows, int cols );
  void copyInplace( const MATRIX &M )
  {
    copyInplace( M._matrix, M.rows(), M.cols() );
  }

  // stomp everything and set to identity
  void setToIdentity();

  // squared sum of matrix entries
  Real sum2();

  // norms of matrix
  Real norm1();
  Real normInf();
  Real norm2() const;
 
  // trace of the matrix
  Real trace();

  // max absolute entry
  Real maxAbsEntry();

  // LU factorization routines -
  //
  // Note that this *WILL* stomp your current matrix with the
  // factorization!
  // 
  // The call to factorLU performs the LU factorization and stomps
  // the original matrix.
  bool factorLU();
  // The call to solveLU uses the precomputed factorization to perform
  // the solve. The solve stomps the original b to store the solution.
  void solveLU(VECTOR& b);
 
  // The call to solveLU uses the precomputed factorization to perform
  // the solve on multiple B's packed into an array. 
  // The solve stomps the original matrix B to store the solution.
  void solveLU(MATRIX& B);

  // Cholesky factorization routines 
  //
  // Note that this *WILL* stomp your current matrix with the
  // factorization!
  // 
  // Same as factorLU, solveLU, but it is YOUR responsibility to 
  // ensure the matrix is symmetric
  bool factorCholesky();
  void solveCholesky(VECTOR& b);

  // solve for linear least squares
  void solveLeastSquares(VECTOR& b);

  // orthogonalize the columns
  void orthogonalize();

  // get a vector containing the desired column
  VECTOR getColumn(int index);
  void setColumn(const VECTOR& column, int index);
  void addColumn(const VECTOR& column, int index);

  // compute the QR factorization of this matrix
  void qr();

  // get submatrix
  MATRIX getSubmatrix(int rowBegin, int totalRows, int colBegin, int totalCols);

  // explicitly compute the inverse, not just the LU decomposition
  //
  // Note -- this is a LOT more expensive than just computing the LU
  // and calling solveLU()!
  bool invert();
  
  // explicitly compute the inverse, not just the Cholesky decomposition
  bool invertSymmetric();

  // are any of the entries a nan?
  bool anyNans();

  // compute the condition number using the full eigendecomposition
  Real conditionNumber();

  // normalize the columns so that their magnitudes are all 1
  void normalizeColumns();

  // scale the columns by the entries in this vector
  void scaleColumns(VECTOR scale);

  // copy upper triangle into lower
  void copyUpperToLower();

  // return the dimensions packed into a vector (for output, usually)
  VECTOR dims() const { VECTOR final(2); final[0] = _rows; final[1] = _cols; return final; };

  static MATRIX cross(const VEC3F& vec);
  static MATRIX cross(const VECTOR& vec);

  // get the matrix exponential
  MATRIX exp();

  // get the derivative of the matrix exponential, assuming that the
  // derivative of the matrix itself is already known
  MATRIX dexp(const MATRIX& dA) const;

  // take the absolute value of all entries
  void absoluteValue();

  // concatecate the columns of the passed in matrix to the current matrix
  void concatenateColumns(MATRIX& newColumns);

  // a unit test function for matrix exponentials
  static bool unitTestExp();

  //////////////////////////////////////////////////////////////////////
  // Static routines for handling array-based matrices.
  // Use with care - nothing here does any bound checking.
  // We assume column-major format in all cases.
  //////////////////////////////////////////////////////////////////////

  // Zero out a matrix
  static inline void clear( Real *A, int rows, int cols )
  {
    memset( (void *)A, 0, rows * cols * sizeof( Real ) );
  }

  // Generate a diagonal matrix
  static inline void diagonal( Real *A, Real *D, int rows )
  {
    for ( int i = 0; i < rows; i++ )
    {
      access( A, rows, rows, i, i ) = D[ i ];
    }
  }

  // Copy one matrix in to another (A <- B)
  static void copy( Real *A, const Real *B, int rows, int cols );

  // Copy a 3x3 matrix
  static void copy( Real *A, const MATRIX3 &B );

  // Access element from a matrix
  static inline Real &access( Real *A, int rows, int cols, int i, int j )
  {
    return A[ i * cols + j ];
  }
  static inline Real access( const Real *A, int rows, int cols, int i, int j )
  {
    return A[ i * cols + j ];
  }

  // Add one matrix to another (A = beta * A + alpha * B
  static void axpy( Real *A, const Real *B, int rows, int cols,
                    Real alpha = 1.0, Real beta = 0.0 );

  // Matrix-matrix multiplication (C = beta * C + alpha * A * B)
  static void gemm( const Real *A, const Real *B, Real *C,
                    int rowsA, int colsA, int rowsB, int colsB,
                    bool transposeA, bool transposeB,
                    Real alpha = 1.0, Real beta = 0.0 );

  // Get transpose (A = B^T)
  static void transpose( Real *A, const Real *B, int rows, int cols );

  // In place transpose for square matrices
  static void transpose( Real *A, int rows );

  // Compute eigenvalues/vectors.  We require everything here,
  // including workspaces, etc.
  // workspace should have size (7 * rows)
  // vectorWorkspace should have size (rows * rows)
  static void eigensystem( Real *A, int rows,
                           Real *eigenvalues, Real *eigenvectors,
                           Real *workspace, Real *vectorWorkspace );

  // 2x2 eigensolver
  static void eigensystem2x2( Real *A, Real *eigenvalues, Real *eigenvectors );

  // 3x3 low precision eigensolver
  static void eigensystem3x3( Real *A, Real *eigenvalues, Real *eigenvectors );

  // Scales a matrix.  A *= alpha
  static void scale( Real *A, int rows, int cols, Real alpha );

  // copy the contents of this submatrix into the passed in
  // matrix, starting at the specified row
  void copiesInto(MATRIX& matrix, int startingRow);

  // copy the contents of this submatrix into the passed in
  // matrix, starting at the specified row and column
  void copiesInto(MATRIX& matrix, int startingRow, int startingCol);

  // return a matrix with 'n' 3x3 identity matrices stacked in a row
  static MATRIX rowOfIdentities(int n);
  
  // return a matrix with 'n' 3x3 identity matrices stacked in a column
  static MATRIX columnOfIdentities(int n);

  // return the outer product of two vectors
  static MATRIX outerProduct(const VECTOR& left, const VECTOR& right);

  // compute the adjoint rigid body matrix form Eqn. 15 of
  // "Lie Group Integratos for Animation and Control of Vehicles"
  static MATRIX rigidAdjoint(const VEC3F& angular, const VEC3F& linear);

  // assuming this is a 4x4 homogeneous matrix, unpack it into the rotation/translation
  void unpackHomogeneous(MATRIX3& rotation, VEC3F& translation);

  // clamp the entries of the matrix smaller than a threshold to zero
  void clampToZero(const Real threshold);

  // stack all the matrices in the vector into one big matrix -- all the column
  // dimensions of all the matrices must match!
  static MATRIX columnOfMatrices(const vector<MATRIX>& matrices);

protected:
  int _rows;
  int _cols;

  Real* _matrix;

  // this stores the pivots if LAPACK LU factorization is called. 
  // If unused, the storage is modest -- certainly worth saving the
  // headache of having to track this memory explicitly.
  int* _pivots;

  // operator << precision parameters
  static int _precision;
  static int _width;

  // repack matrix into fortran ordering in preparation for a LAPACK call
  // that does not support transposes
  void formatFortran();

  // unpack matrix from a LAPACK call that only supports fortran ordering
  void unformatFortran();
};

// overloaded operators
VECTOR operator*(const MATRIX& A, const VECTOR& x);
VECTOR operator*(VECTOR& x, MATRIX& A);
MATRIX operator*(const MATRIX& A, const Real alpha);
MATRIX operator*(const Real alpha, const MATRIX& A);
MATRIX operator*(const MATRIX& A, const MATRIX& B);
ostream& operator<<(ostream &out, const MATRIX& matrix);
MATRIX operator-(const MATRIX& A, const MATRIX& B);
MATRIX operator+(const MATRIX& A, const MATRIX& B);

// multiply by the transpose of A
VECTOR operator^(const MATRIX& A, const VECTOR& x);

// multiply by the transpose of A -- NOT B!!!
MATRIX operator^(const MATRIX& A, const MATRIX& B);

#endif
