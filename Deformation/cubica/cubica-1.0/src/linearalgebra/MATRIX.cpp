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
#include "SETTINGS.h"

#if USING_MKL
#include "MATRIX_FAST.cpp"
#elif USING_OSX
#include "MATRIX_OSX.cpp"
#else
#include "MATRIX_DEBUG.cpp"
#endif

#include <iomanip>

int MATRIX::_precision = 4;
int MATRIX::_width = 14;

//////////////////////////////////////////////////////////////////////
// Get the matrix exponential
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::exp()
{
  /*
  // only implementing Taylor series version here. As long as the norm
  // is small enough, this is the right computation method. See:
  // 
  // Fung, Computation of the matrix exponential and its derivatives by
  // scaling and squaring, international journal for numerical methods
  // in engineering, 2004, 59:1273-1286.
  assert(this->norm2() < 0.2);

  MATRIX final(_rows, _cols);
  final.setToIdentity();
  int k = 1;
  MATRIX& A = *this;
  MATRIX Aminus1 = final;
  Real maxEntry = 1.0;

  while (maxEntry > 1e-16)
  {
    // Eqn. 22
    Aminus1 = (1.0 / (Real)k) * (Aminus1 * A);
    final += Aminus1;
    maxEntry = Aminus1.maxAbsEntry();
    k++;
  }

  return final;
  */

  // implementing full scale and square version here. See:
  // 
  // Fung, Computation of the matrix exponential and its derivatives by
  // scaling and squaring, international journal for numerical methods
  // in engineering, 2004, 59:1273-1286.
  //
  int m = 10;
  int N = 8;

  MATRIX Bk(_rows, _cols);
  Bk.setToIdentity();
  MATRIX& A = *this;

  MATRIX BN = Bk;
  BN *= 0;
  Real inv = 1.0 / pow(2.0, N);
  for (int k = 1; k <= m; k++)
  {
    Bk = (1.0 / k) * inv * Bk * A;
    BN += Bk; 
  }

  MATRIX final = BN;
  for (int k = 0; k < N; k++)
    final = 2 * final + final * final;

  MATRIX identity(_rows, _cols);
  identity.setToIdentity();
  final += identity;

  return final;
}

//////////////////////////////////////////////////////////////////////
// get the derivative of the matrix exponential, assuming that the
// derivative of the matrix itself is already known
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::dexp(const MATRIX& dA) const
{
  // only implementing Taylor series version here. As long as the norm
  // is small enough, this is the right computation method. See:
  // 
  // Fung, Computation of the matrix exponential and its derivatives by
  // scaling and squaring, international journal for numerical methods
  // in engineering, 2004, 59:1273-1286.
  if (this->norm2() < 0.2)
  {
    MATRIX final(_rows, _cols);
    int k = 1;
    const MATRIX& A = *this;
    MATRIX Aminus1 = final;
    Aminus1.setToIdentity();

    MATRIX dAminus1 = final;
    Real maxEntry = 1.0;

    while (maxEntry > 1e-16)
    {
      // Eqn. 24a
      dAminus1 = (1.0 / (Real)k) * (Aminus1 * dA + dAminus1 * A);

      // Eqn. 22
      Aminus1 = (1.0 / (Real)k) * (Aminus1 * A);

      final += dAminus1;
      maxEntry = dAminus1.maxAbsEntry();
      k++;
    }

    return final;
  }

  // implementing full scale and square version here. See:
  // 
  // Fung, Computation of the matrix exponential and its derivatives by
  // scaling and squaring, international journal for numerical methods
  // in engineering, 2004, 59:1273-1286.
  //
  int m = 10;
  int N = 8;

  MATRIX Bk(_rows, _cols);
  Bk.setToIdentity();
  const MATRIX& A = *this;

  MATRIX BkAlpha(_rows, _cols);

  MATRIX BN = Bk;
  BN *= 0;
  MATRIX BNAlpha(_rows, _cols);

  Real inv = 1.0 / pow(2.0, N);
  for (int k = 1; k <= m; k++)
  {
    Real invK = (1.0 / k) * inv;
    BkAlpha = invK * (Bk * dA + BkAlpha * A); 
    Bk = invK * Bk * A;
    BNAlpha += BkAlpha;
    BN += Bk; 
  }

  MATRIX Bm1 = BN;
  MATRIX BAlpham1 = BNAlpha;
  for (int k = 0; k < N; k++)
  {
    BAlpham1 = 2 * BAlpham1 + BAlpham1 * Bm1 + Bm1 * BAlpham1;
    Bm1 = 2 * Bm1 + Bm1 * Bm1;
  }

  return BAlpham1;
}

//////////////////////////////////////////////////////////////////////
// a unit test function for matrix exponentials -- this exploits
// the analytical solution from Section 6 in:
//
// Fung, Computation of the matrix exponential and its derivatives by
// scaling and squaring, international journal for numerical methods
// in engineering, 2004, 59:1273-1286.
//////////////////////////////////////////////////////////////////////
bool MATRIX::unitTestExp()
{
  MATRIX A(2,2);

  // you can plug in whatever small numbers you want here
  //Real a = 0.0321;
  Real a = 1e-8;
  //Real a = 0.000123;
  //Real a = 0.00123;
  //Real a = 0.0123;
  //Real a = 0.123;
  //Real a = 1.23;
  //Real a = 10;
  Real a2 = a * a;
  Real a3 = a * a * a;
  Real ea = std::exp(a);
  A(0,0) = a3;
  A(0,1) = a2 * (1 - a);
  A(1,0) = a2 * (1 + a);
  A(1,1) = a * (1 - a2);
  cout << " A norm: " << A.norm2() << endl;

  MATRIX B(2,2);
  B(0,0) = 1 + a2 * (ea - 1);
  B(0,1) = (ea - 1) * (a - a2);
  B(1,0) = (ea - 1) * (a + a2);
  B(1,1) = ea - (ea - 1) * a2;

  MATRIX dA(2,2);
  dA(0,0) = 3 * a2;
  dA(0,1) = 2 * a - 3 * a2;
  dA(1,0) = 3 * a2 + 2 * a;
  dA(1,1) = 1 - 3 * a2;

  MATRIX dB(2,2);
  dB(0,0) = (a2 + 2 * a) * ea - 2 * a;
  dB(0,1) = -(a2 + a - 1) * ea + 2 * a - 1;
  dB(1,0) = (a2 + 3 * a + 1) * ea - 2 * a - 1;
  dB(1,1) = -(a2 + 2 * a - 1) * ea + 2 * a;

  // should be the same
  MATRIX diff = A.exp() - B;
  cout << " exp ground truth: " << B << endl;
  cout << " exp computed: " << A.exp() << endl;
  cout << " diff: " << sqrt(diff.sum2()) << endl;
  cout << " relative error: " << sqrt(diff.sum2()) / sqrt(B.sum2()) << endl;
  if (sqrt(diff.sum2()) > 1e-14) return false;

  // should be the same
  MATRIX test = A.dexp(dA);
  diff = test - dB;
  cout << " dexp ground truth: " << dB << endl;
  cout << " dexp computed: " << A.dexp(dA) << endl;
  cout << " diff: " << sqrt(diff.sum2()) << endl;
  cout << " relative error: " << sqrt(diff.sum2())  / sqrt(dB.sum2()) << endl;
  if (sqrt(diff.sum2()) > 1e-14) return false;

  // peek at the finite difference version
  {
    Real delta = 1e-6 * a;
    a -= delta;
    Real a2 = a * a;
    Real a3 = a * a * a;
    MATRIX newA(2,2);
    newA(0,0) = a3;
    newA(0,1) = a2 * (1 - a);
    newA(1,0) = a2 * (1 + a);
    newA(1,1) = a * (1 - a2);

    MATRIX finiteDiff = A.exp() - newA.exp();
    finiteDiff *= 1.0 / delta;
    cout << endl;
    cout << " finite diff derivative: " << finiteDiff << endl;
  }

  return true;
}

//////////////////////////////////////////////////////////////////////
// copy the contents of this submatrix into the passed in
// matrix, starting at the specified row
//////////////////////////////////////////////////////////////////////
void MATRIX::copiesInto(MATRIX& matrix, int startingRow)
{
  for (int x = 0; x < _rows; x++)
    for (int y = 0; y < _cols; y++)
      matrix(startingRow + x, y) = (*this)(x,y);
}

//////////////////////////////////////////////////////////////////////
// copy the contents of this submatrix into the passed in
// matrix, starting at the specified row and column
//////////////////////////////////////////////////////////////////////
void MATRIX::copiesInto(MATRIX& matrix, int startingRow, int startingCol)
{
  for (int x = 0; x < _rows; x++)
    for (int y = 0; y < _cols; y++)
      matrix(startingRow + x, startingCol + y) = (*this)(x,y);
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, const MATRIX& matrix)
{
  int oldPrecision = out.precision();

  out << "[" << endl;
  for (int row = 0; row < matrix.rows(); row++)
  {
    for (int col = 0; col < matrix.cols(); col++)
    {
      //out << matrix(row, col) << " ";
      out.width(MATRIX::coutWidth());
      out << setprecision(MATRIX::coutPrecision()) << matrix(row, col) << " ";
    }
    out << endl;
  }
  out << "]" << endl;
  out << resetiosflags(ios::floatfield);
  out.precision(oldPrecision);
  return out;
}

//////////////////////////////////////////////////////////////////////
// get the rows means of all the rows
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::rowMeans()
{
  MATRIX& data = *this;
  VECTOR means(_rows);

  // subtract out the means
  int rows = data.rows();
  int cols = data.cols();
  for (int x = 0; x < rows; x++)
  {
    Real mean = 0.0;
    for (int y = 0; y < cols; y++)
      mean += data(x,y);
    mean /= cols;

    means[x] = mean;
  }
  return means;
}

//////////////////////////////////////////////////////////////////////
// get a vector corrsesponding to a matrix row
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::getRow(int index) const
{
  VECTOR final(_cols);
  const Real* start = row(index);
  for (int x = 0; x < _cols; x++)
    final[x] = start[x];

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a matrix with 'n' 3x3 identity matrices stacked in a row
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::rowOfIdentities(int n)
{
  MATRIX final(3, 3 * n);

  for (int x = 0; x < 3 * n; x++)
    final(x % 3, x) = 1.0;

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a matrix with 'n' 3x3 identity matrices stacked in a column
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::columnOfIdentities(int n)
{
  MATRIX final(3 * n, 3);

  for (int x = 0; x < 3 * n; x++)
    final(x, x % 3) = 1.0;

  return final;
}

//////////////////////////////////////////////////////////////////////
// return the outer product of two vectors
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::outerProduct(const VECTOR& left, const VECTOR& right)
{
  MATRIX final(left.size(), right.size());
  for (int x = 0; x < right.size(); x++)
    for (int y = 0; y < left.size(); y++)
      final(y,x) = left[y] * right[x];

  return final;
}

//////////////////////////////////////////////////////////////////////
// compute the adjoint rigid body matrix form Eqn. 15 of
// "Lie Group Integrators for Animation and Control of Vehicles"
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::rigidAdjoint(const VEC3F& angular, const VEC3F& linear)
{
  MATRIX final(6,6);

  MATRIX3 angularHat = MATRIX3::cross(angular);
  MATRIX3 linearHat = MATRIX3::cross(linear);

  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
    {
      final(x,y) = angularHat(x,y);
      final(3 + x, 3 + y) = angularHat(x,y);
      final(3 + x, y) = linearHat(x,y);
    }
  return final;
}

//////////////////////////////////////////////////////////////////////
// assuming this is a 4x4 homogeneous matrix, unpack it into the 
// rotation/translation
//////////////////////////////////////////////////////////////////////
void MATRIX::unpackHomogeneous(MATRIX3& rotation, VEC3F& translation)
{
  assert(_rows == 4);
  assert(_cols == 4);

  for (int x = 0; x < 3; x++)
  {
    for (int y = 0; y < 3; y++)
      rotation(x,y) = (*this)(x,y);
    translation[x] = (*this)(x, 3);
  }
}

//////////////////////////////////////////////////////////////////////
// clamp the entries of the matrix smaller than a threshold to zero
//////////////////////////////////////////////////////////////////////
void MATRIX::clampToZero(const Real threshold)
{
  for (int x = 0; x < _rows; x++)
    for (int y = 0; y < _cols; y++)
    {
      if (fabs((*this)(x,y)) < threshold)
        (*this)(x,y) = 0.0;
    }
}

#include "MATRIX_EIGS3X3.cpp"

//////////////////////////////////////////////////////////////////////
// Call one of the eigensolvers from 
//
// Joachim Kopp
// Efficient numerical diagonalization of hermitian 3x3 matrices
// Int. J. Mod. Phys. C 19 (2008) 523-548
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem3x3(VECTOR& eigenvalues, MATRIX& eigenvectors)
{
  /*
  double A[3][3];
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
      A[x][y] = (*this)(x,y);

  double Q[3][3];
  double w[3];

  hybrid(A, Q, w);
  */

  eigenvectors.resizeAndWipe(3,3);
  eigenvalues.resizeAndWipe(3);
  hybrid(*this, eigenvectors, eigenvalues);

  /*
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
      eigenvectors(x,y) = Q[x][y];

  for (int x = 0; x < 3; x++)
    eigenvalues[x] = w[x];
    */
}

//////////////////////////////////////////////////////////////////////
// 3x3 low precision eigensolver
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem3x3( Real *A, Real *eigenvalues, Real *eigenvectors )
{
  MATRIX matrix(3,3);
  MATRIX vectors(3,3);
  VECTOR values(3);

  int i = 0;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++, i++)
      matrix(x,y) = A[i];

  hybrid(matrix, vectors, values);

  i = 0;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++, i++)
      eigenvectors[i] = vectors(x,y);

  for (int x = 0; x < 3; x++)
    eigenvalues[x] = values[x];
  /*
  double matrix[3][3];
  int i = 0;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++, i++)
      matrix[x][y] = A[i];

  double Q[3][3];
  double w[3];

  hybrid(matrix, Q, w);

  i = 0;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++, i++)
      eigenvectors[i] = Q[x][y];

  for (int x = 0; x < 3; x++)
    eigenvalues[x] = w[x];
    */
}

//////////////////////////////////////////////////////////////////////
// stack all the matrices in the vector into one big matrix -- all the column
// dimensions of all the matrices must match!
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::columnOfMatrices(const vector<MATRIX>& matrices)
{
  assert(matrices.size() > 0);
  int totalCols = matrices[0].cols();
  int totalRows = 0;
  for (unsigned int x = 0; x < matrices.size(); x++)
  {
    assert(matrices[x].cols() == totalCols);
    totalRows += matrices[x].rows();
  }

  MATRIX final(totalRows, totalCols);

  int currentRow = 0;
  for (unsigned int x = 0; x < matrices.size(); x++)
  {
    const MATRIX& current = matrices[x];
    for (int i = 0; i < current.rows(); i++)
      for (int j = 0; j < current.cols(); j++)
        final(currentRow + i, j) = current(i,j);
    currentRow += matrices[x].rows();
  }
  
  return final;
}

//////////////////////////////////////////////////////////////////////
// concatecate the columns of the passed in matrix to the current 
// matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::concatenateColumns(MATRIX& newColumns)
{
  assert(newColumns.rows() == this->rows());

  int rows = this->rows();
  int cols = this->cols() + newColumns.cols();
  MATRIX newMatrix(this->rows(), this->cols() + newColumns.cols());

  // copy in the current contents
  for (int y = 0; y < this->cols(); y++)
    for (int x = 0; x < rows; x++)
      newMatrix(x,y) = (*this)(x,y);

  // copy in the new columns
  for (int y = 0; y < newColumns.cols(); y++)
    for (int x = 0; x < rows; x++)
      newMatrix(x,y + this->cols()) = newColumns(x,y);

  // set the current matrix to the final contents
  *this = newMatrix;
}
