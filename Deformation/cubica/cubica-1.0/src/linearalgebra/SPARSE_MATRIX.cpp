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
// SPARSE_MATRIX.h: interface for the SPARSE_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#include <map>
#include "SPARSE_MATRIX.h"

#include <fstream>

//////////////////////////////////////////////////////////////////////
// Constructor for the sparse matrix
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX::SPARSE_MATRIX(int rows, int cols) :
  _rows(rows), _cols(cols), _dud(0.0f)
{
}

SPARSE_MATRIX::SPARSE_MATRIX() :
  _rows(0), _cols(0), _dud(0.0f)
{
}

SPARSE_MATRIX::SPARSE_MATRIX(const MATRIX& matrix) :
  _dud(0.0f)
{
  _rows = matrix.rows();
  _cols = matrix.cols();

  for (int y = 0; y < _rows; y++)
    for (int x = 0; x < _cols; x++)
      (*this)(y,x) = matrix(y,x);
}

SPARSE_MATRIX::SPARSE_MATRIX(MATRIX& matrix, Real threshold) :
  _dud(0.0f)
{
  _rows = matrix.rows();
  _cols = matrix.cols();

  for (int y = 0; y < _rows; y++)
    for (int x = 0; x < _cols; x++)
      if (fabs(matrix(y,x)) > threshold)
        (*this)(y,x) = matrix(y,x);
}

//////////////////////////////////////////////////////////////////////
// check if an entry already exists
//////////////////////////////////////////////////////////////////////
bool SPARSE_MATRIX::exists(int row, int col) const
{
  pair<int,int> index(row, col);
  map<pair<int,int>, Real>::const_iterator i = _matrix.find(index);
  if (i != _matrix.end())
    return true;
  return false;
}

//////////////////////////////////////////////////////////////////////
// return a reference to an entry
//////////////////////////////////////////////////////////////////////
Real& SPARSE_MATRIX::operator()(int row, int col)
{
  // bounds check
  assert(col >= 0);
  assert(row >= 0); 
  assert(col < _cols);
  assert(row < _rows);

  // lookup the entry
  pair<int,int> index(row, col);
  map<pair<int,int>, Real>::iterator i = _matrix.find(index);
  if (i != _matrix.end())
    return i->second;

  // if it doesn't exist, create it
  _matrix[index] = 0.0f;
  return _matrix[index];
}

//////////////////////////////////////////////////////////////////////
// const accessor to matrix entries
//////////////////////////////////////////////////////////////////////
Real SPARSE_MATRIX::constEntry(int row, int col) const
{
  // bounds check
  assert(col >= 0);
  assert(row >= 0); 
  assert(col < _cols);
  assert(row < _rows);

  // lookup the entry
  pair<int,int> index(row, col);
  map<pair<int,int>, Real>::const_iterator i = _matrix.find(index);
  if (i != _matrix.end())
    return i->second;

  return 0.0;
}

//////////////////////////////////////////////////////////////////////
// dump all the matrix entries to a Matlab file
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::writeToMatlab(string filename, string varName)
{
  FILE* file;
  file = fopen(filename.c_str(), "w");

  // create the matrix and sparsify it
  fprintf(file, "%s = [];\n", varName.c_str());
  fprintf(file, "%s = sparse(%s);\n", varName.c_str(), varName.c_str());

  // iterate through all the entries
  map<pair<int,int>, Real>::iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    pair<int,int> index = i->first;
    Real value = i->second;
    fprintf(file, "%s(%i,%i) = %.16f;\n", varName.c_str(), index.first + 1, index.second + 1, value);
  }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// write to a simple binary format
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::writeToBinary(string filename)
{
	ofstream f( filename.c_str(), ios::binary );
	if( f.fail() )
	{
		cerr << "Couldn't write sparse matrix to binary file: " << filename.c_str() << endl;
	}
	else
	{
		// (stevenan) 2008-11-24 13:28:02
		// We write a -1 as a header, to specify that this is a
		// struct-of-arrays-format binary sparse matrix.
		// This is what our matlab read function expects.
		// This is to maintain backwards compatibility with the old
		// array-of-structs format.
		int header = -1;
		f.write( (char*)&header, sizeof(int) );

		f.write( (char*)&_rows, sizeof(int) );
		f.write( (char*)&_cols, sizeof(int) );
		int nnz = _matrix.size();
		f.write( (char*)&nnz, sizeof(int) );

		// Write row indices
		map<pair<int,int>, Real>::iterator i;
		for( i = _matrix.begin(); i != _matrix.end(); ++i ) {
			pair<int,int> index = i->first;
			// This does NOT add 1 to the index. The matlab read function will do that on its own.
			f.write( (char*)&index.first, sizeof(int) );
		}
		// Write columns
		for( i = _matrix.begin(); i != _matrix.end(); ++i ) {
			pair<int,int> index = i->first;
			// This does NOT add 1 to the index. The matlab read function will do that on its own.
			f.write( (char*)&index.second, sizeof(int) );
		}
		// Write values
		for( i = _matrix.begin(); i != _matrix.end(); ++i ) {
			pair<int,int> index = i->first;

			// Read and write everything as a double
			double valueBin = (double)i->second;

			f.write( (char*)&valueBin, sizeof(double) );
		}
	}
	f.close();
}

//////////////////////////////////////////////////////////////////////
// count the number of non-zeros per row and return them
//////////////////////////////////////////////////////////////////////
vector<int> SPARSE_MATRIX::nonZerosPerRow()
{
  vector<int> count(_rows);

  // iterate through all the entries
  map<pair<int,int>, Real>::iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
    count[i->first.first]++;

  return count;
}

//////////////////////////////////////////////////////////////////////
// Print a specific row
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::printRow(int row)
{
  if (row < 0) return;
  if (row > _rows) return;

  for (int x = 0; x < _cols; x++)
    if (exists(x, row))
      cout << "A(" << row << "," << x << "): " << (*this)(x, row) << endl;
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, SPARSE_MATRIX& matrix)
{
  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator i;
  const map<pair<int,int>, Real>& data = matrix.matrix();
  for (i = data.begin(); i != data.end(); i++)
  {
    const pair<int,int> index = i->first;
    const Real value = i->second;
    out << "(" << index.first << "," << index.second << ") = " <<  value << endl;
  }
  return out;
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B += alpha * A, where B is this matrix
//
// Note that axpy actually applies to vectors, but in this
// case we can just treat the matrix as a vector and multiply
// all its elements by alpha
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::axpy(Real alpha, SPARSE_MATRIX& A)
{
  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator i;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (i = data.begin(); i != data.end(); i++)
  {
    const pair<int,int> index = i->first;
    const Real value = i->second;
    (*this)(index.first, index.second) += alpha * value;
  }
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator*(const SPARSE_MATRIX& A, const VECTOR& x) 
{
  assert(A.cols() == x.size());

  VECTOR y(A.rows());

  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator iter;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    const pair<int,int> index = iter->first;
    const Real value = iter->second;
    int i = index.first;
    int j = index.second;
    //y(i) += x(j) * value;
    y(i) += x[j] * value;
  }
  return y;
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX operator*(const Real& alpha, const SPARSE_MATRIX& A) 
{
  SPARSE_MATRIX final(A);
  final *= alpha;
  return final;
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator^(const SPARSE_MATRIX& A, const VECTOR& x) 
{
  assert(A.rows() == x.size());

  VECTOR y(A.cols());

  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator iter;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    const pair<int,int> index = iter->first;
    const Real value = iter->second;
    int i = index.first;
    int j = index.second;
    y(j) += x[i] * value;
  }
  return y;
}

//////////////////////////////////////////////////////////////////////
// sparse matrix-vector multiply
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX operator*(SPARSE_MATRIX& A, Real& alpha) 
{
  SPARSE_MATRIX B(A);
  B *= alpha;
  return B;
}

//////////////////////////////////////////////////////////////////////
// sparse-sparse matrix multiply
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX operator*(const SPARSE_MATRIX& A, const SPARSE_MATRIX& B)
{
  SPARSE_MATRIX C(A.rows(), B.cols());

  const map<pair<int,int>, Real>& matrix = A.matrix(); 

  map<pair<int,int>, Real>::const_iterator iter;
  for (iter = matrix.begin(); iter != matrix.end(); iter++)
  {
    pair<int,int> index = iter->first;
    int row = index.first;
    int col = index.second;
    Real entry = iter->second;

    for (int x = 0; x < B.cols(); x++)
      if (B.exists(col, x))
        //C(row, x) += B(col, x) * entry;
        C(row, x) += B.constEntry(col, x) * entry;
  }

  return C;
}

//////////////////////////////////////////////////////////////////////
// sparse-full matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(const SPARSE_MATRIX& A, const MATRIX& B)
{
  MATRIX C(A.rows(), B.cols());

  const map<pair<int,int>, Real>& matrix = A.matrix(); 

  map<pair<int,int>, Real>::const_iterator iter;
  for (iter = matrix.begin(); iter != matrix.end(); iter++)
  {
    pair<int,int> index = iter->first;
    int row = index.first;
    int col = index.second;
    Real entry = iter->second;

    for (int x = 0; x < B.cols(); x++)
      C(row, x) += B(col, x) * entry;
  }

  return C;
} 

//////////////////////////////////////////////////////////////////////
// full-sparse matrix multiply
// (Really, really, not optimized)
//////////////////////////////////////////////////////////////////////
MATRIX operator*(const MATRIX& A, const SPARSE_MATRIX& B)
{
  return A * B.full();
}

//////////////////////////////////////////////////////////////////////
// full^T-sparse matrix multiply
// (Not terribly optimized)
//////////////////////////////////////////////////////////////////////
MATRIX operator^(MATRIX& A, SPARSE_MATRIX& B)
{
  MATRIX C(A.cols(), B.cols());

  for (int i = 0; i < A.cols(); i++)
    for (int j = 0; j < B.cols(); j++)
      for (int k = 0; k < A.rows(); k++)
        if (B.exists(k,j))
          C(i,j) += A(k, i) * B(k, j);

  return C;
}

//////////////////////////////////////////////////////////////////////
// scale matrix by a scalar
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& SPARSE_MATRIX::operator*=(const Real& alpha) 
{
  // iterate through all the entries
  map<pair<int,int>, Real>::iterator iter;
  for (iter = _matrix.begin(); iter != _matrix.end(); iter++)
    iter->second *= alpha;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// add two sparse matrices together
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& SPARSE_MATRIX::operator+=(const SPARSE_MATRIX& A) 
{
  assert(A.rows() == _rows && A.cols() == _cols);

  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator iter;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    const pair<int,int> index = iter->first;
    const Real value = iter->second;
    int i = index.first;
    int j = index.second;
    (*this)(i,j) += value;
  }

  return *this;
}

//////////////////////////////////////////////////////////////////////
// subtact two sparse matrices
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& SPARSE_MATRIX::operator-=(const SPARSE_MATRIX& A) 
{
  assert(A.rows() == _rows && A.cols() == _cols);

  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator iter;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    const pair<int,int> index = iter->first;
    const Real value = iter->second;
    int i = index.first;
    int j = index.second;
    (*this)(i,j) -= value;
  }

  return *this;
}

//////////////////////////////////////////////////////////////////////
// copy A into the current object
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::copies(SPARSE_MATRIX& A)
{
  assert(A.cols() == _cols && A.rows() == _rows);

  if (A._cols == _cols && A._rows == _rows)
  {
    // first zero out what's already there, but don't do
    // a total clear, or else it will have to reallocate all the
    // entries
    map<pair<int,int>, Real>::iterator i;
    for (i = _matrix.begin(); i != _matrix.end(); i++)
      i->second = 0.0;

    // next copy in all the entries from A. If an entry already
    // existed at the matrix entry, we will save the memory
    // allocation
    for (i = A._matrix.begin(); i != A._matrix.end(); i++)
    {
      const pair<int,int> index = i->first;
      (*this)(index.first, index.second) = i->second;
    }
    return;
  }
  
  // else do an expensive deep copy
  _matrix = map<pair<int,int>, Real>(A.matrix());
}

//////////////////////////////////////////////////////////////////////
// copy A into the current object
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::subcopies(MATRIX& A, int row, int col)
{
  for (int y = 0; y < A.rows(); y++)
    for (int x = 0; x < A.cols(); x++)
      (*this)(y + row, x + col) = A(y,x);
}

//////////////////////////////////////////////////////////////////////
// add entry into the current object
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::add(Real& entry, int row, int col)
{
  (*this)(row, col) += entry;
}

//////////////////////////////////////////////////////////////////////
// subtract entry into the current object
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::subtract(Real& entry, int row, int col)
{
  (*this)(row, col) -= entry;
}

//////////////////////////////////////////////////////////////////////
// add A into the current object
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::add(MATRIX& A, int row, int col)
{
  for (int y = 0; y < A.rows(); y++)
    for (int x = 0; x < A.cols(); x++)
      (*this)(y + row, x + col) += A(y,x);
}

//////////////////////////////////////////////////////////////////////
// subtract A from the current object
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::subtract(MATRIX& A, int row, int col)
{
  for (int y = 0; y < A.rows(); y++)
    for (int x = 0; x < A.cols(); x++)
      (*this)(y + row, x + col) -= A(y,x);
}

//////////////////////////////////////////////////////////////////////
// All the matrix entries repackaged into vectors,
// ie _matrix(rows[x], cols[x]) = vals[x];
// upper triangular entries only
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::entriesUpper(vector<int>& rows, vector<int>& cols, vector<Real>& values)
{
  // wipe previous entries
  rows.clear();
  cols.clear();
  values.clear();
  
  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    const pair<int,int> index = i->first;
    const Real value = i->second;

    // ensure lower triangle only
    int row = index.first;
    int col = index.second;
    if (row > col) continue;

    rows.push_back(index.first);
    cols.push_back(index.second);
    values.push_back(value);
  }
}

//////////////////////////////////////////////////////////////////////
// All the matrix entries repackaged into vectors,
// ie _matrix(rows[x], cols[x]) = vals[x];
// lower triangular entries only
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::entriesLower(vector<int>& rows, vector<int>& cols, vector<Real>& values)
{
  // wipe previous entries
  rows.clear();
  cols.clear();
  values.clear();
 
  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    const pair<int,int> index = i->first;
    const Real value = i->second;

    // ensure lower triangle only
    int row = index.first;
    int col = index.second;
    if (col > row) continue;

    rows.push_back(index.first);
    cols.push_back(index.second);
    values.push_back(value);
  }
}

//////////////////////////////////////////////////////////////////////
// All the matrix entries repackaged into vectors,
// ie _matrix(rows[x], cols[x]) = vals[x];
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::entries(vector<int>& rows, vector<int>& cols, vector<Real>& values)
{
  // wipe previous entries
  rows.clear();
  cols.clear();
  values.clear();
  
  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    const pair<int,int> index = i->first;
    const Real value = i->second;
    rows.push_back(index.first);
    cols.push_back(index.second);
    values.push_back(value);
  }
}

//////////////////////////////////////////////////////////////////////
// All the matrix entries repackaged into vectors,
// ie _matrix(rows[x], cols[x]) = vals[x];
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::entries(vector<int>& rows, vector<int>& cols, vector<const Real*>& values)
{
  // wipe previous entries
  rows.clear();
  cols.clear();
  values.clear();
  
  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    const pair<int,int> index = i->first;
    //Real value = i->second;
    rows.push_back(index.first);
    cols.push_back(index.second);
    values.push_back(&(i->second));
  }
}

//////////////////////////////////////////////////////////////////////
// Set the matrix to zero. Note this will *NOT* stomp the underlying map!
// It will instead set all current entries to zero so that we are not
// forced to reallocate the sparsity structure again.
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::clear()
{
  map<pair<int,int>, Real>::iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
    i->second = 0.0;
}

//////////////////////////////////////////////////////////////////////
// Convert to a full matrix
//////////////////////////////////////////////////////////////////////
MATRIX SPARSE_MATRIX::full() const
{
  MATRIX matrix(_rows, _cols);
  
  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator i;
  const map<pair<int,int>, Real>& data = this->matrix();
  for (i = data.begin(); i != data.end(); i++)
  {
    const pair<int,int> index = i->first;
    const Real value = i->second;
    matrix(index.first, index.second) = value;
  }
  return matrix;
}

//////////////////////////////////////////////////////////////////////
// Erase a row/column
//////////////////////////////////////////////////////////////////////
void SPARSE_MATRIX::eraseRowColumn(int rowCol)
{
  pair<int,int> entry(rowCol, rowCol);

  map<pair<int,int>, Real>::iterator i = _matrix.find(entry);

  if (i == _matrix.end())
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Row/Col erasure failed!" << endl;
    return;
  }
  while (i->first.first == rowCol)
  {
    i->second = 0.0;

    // zero out the equivalent column entry
    pair<int,int> transpose(i->first.second, i->first.first);
    _matrix[transpose] = 0.0;

    i++;
  }

  map<pair<int,int>, Real>::iterator ri = _matrix.find(entry);
  while (ri->first.first == rowCol)
  {
    ri->second = 0.0;

    // zero out the equivalent column entry
    pair<int,int> transpose(ri->first.second, ri->first.first);
    _matrix[transpose] = 0.0;
    ri--;
  }
}

//////////////////////////////////////////////////////////////////////
// 1 norm of the matrix
// - maximum absolute column sum
//////////////////////////////////////////////////////////////////////
Real SPARSE_MATRIX::norm1()
{
  VECTOR colSums(_cols);
  map<pair<int,int>, Real>::iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    pair<int,int> index = i->first;
    Real entry = i->second;
    colSums[index.second] += fabs(entry);
  }
  return colSums.maxValue();
}

//////////////////////////////////////////////////////////////////////
// Inf norm of the matrix
// - maximum absolute row sum
//////////////////////////////////////////////////////////////////////
Real SPARSE_MATRIX::normInf()
{
  VECTOR rowSums(_rows);
  map<pair<int,int>, Real>::iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    pair<int,int> index = i->first;
    Real entry = i->second;
    rowSums[index.first] += fabs(entry);
  }
  return rowSums.maxValue();
}

//////////////////////////////////////////////////////////////////////
// Squared sum of the matrix entries
//////////////////////////////////////////////////////////////////////
Real SPARSE_MATRIX::sum2()
{
  Real sum = 0.0;
  map<pair<int,int>, Real>::iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    pair<int,int> index = i->first;
    Real entry = i->second;
    sum += entry * entry;
  }
  return sum;
}

//////////////////////////////////////////////////////////////////////
// sum of the matrix entries
//////////////////////////////////////////////////////////////////////
Real SPARSE_MATRIX::sum()
{
  Real sum = 0.0;
  map<pair<int,int>, Real>::iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    pair<int,int> index = i->first;
    Real entry = i->second;
    sum += entry;
  }
  return sum;
}

//////////////////////////////////////////////////////////////////////
// is the matrix symmetric?
//////////////////////////////////////////////////////////////////////
bool SPARSE_MATRIX::isSymmetric()
{
  map<pair<int,int>, Real>::iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    pair<int,int> index = i->first;
    Real value = i->second;

    pair<int,int> mirror;
    mirror.first = index.second;
    mirror.second = index.first;
    map<pair<int,int>, Real>::iterator j = _matrix.find(mirror);
    if (j == _matrix.end())
      return false;

    Real mirrorValue = j->second;

    if (fabs(value - mirrorValue) > 1e-8)
      return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////
// compute and return the transpose
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX SPARSE_MATRIX::transpose()
{
  SPARSE_MATRIX toReturn(_cols, _rows);

  map<pair<int,int>, Real>::iterator i;
  for (i = _matrix.begin(); i != _matrix.end(); i++)
  {
    pair<int,int> index = i->first;
    Real entry = i->second;
    toReturn(index.second, index.first) = entry;
  }

  return toReturn;
}
