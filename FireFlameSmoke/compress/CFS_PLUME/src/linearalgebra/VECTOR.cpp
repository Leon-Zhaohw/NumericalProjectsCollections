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
#include "VECTOR.h"
bool VECTOR::printVertical = true;

#include <algorithm>
using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor for the full vector
//////////////////////////////////////////////////////////////////////
VECTOR::VECTOR() :
  _size(0), _sseSize(0)
{
  _vector = NULL;
}

VECTOR::VECTOR(int size) :
  _size(size)
{
  _sseSize = (size % 2) ? size : size + 1;

  _vector = (double*)malloc(sizeof(double) * _sseSize);
  clear();
}

VECTOR::VECTOR(const char* filename) :
  _size(0),
  _sseSize(0),
  _vector(NULL)
{
  read(filename);
}

VECTOR::VECTOR(const VECTOR& v) 
{
  _size = v._size;
  _sseSize = v._sseSize;

  _vector = (double*)malloc(sizeof(double) * _sseSize);
  for (int x = 0; x < _size; x++)
    _vector[x] = v._vector[x];
}

VECTOR::VECTOR(const vector<Real>& v)
{
  _size = v.size();
  _sseSize = (_size % 2) ? _size : _size + 1;

  _vector = (double*)malloc(sizeof(double) * _sseSize);
  for (int x = 0; x < _size; x++)
    _vector[x] = v[x];
}

VECTOR::~VECTOR()
{
  if (_vector) {
    free(_vector);
    _vector = NULL;
  }
}

//////////////////////////////////////////////////////////////////////
// resize and wipe the vector
//////////////////////////////////////////////////////////////////////
void VECTOR::resizeAndWipe(int size)
{
  if (_size != size)
  {
    if (_vector) free(_vector);
    _size = size;
    _sseSize = (_size % 2) ? _size : _size + 1;

    _vector = (double*)malloc(sizeof(double) * _sseSize);
  }
  clear();
}

//////////////////////////////////////////////////////////////////////
// write the vector to a file
//////////////////////////////////////////////////////////////////////
void VECTOR::write(FILE* file)
{
  // write dimensions
  fwrite((void*)&_size, sizeof(int), 1, file);

  // write data
  if (sizeof(Real) == sizeof(float))
  {
    double* vecDouble = new double[_size];
    for (int x = 0; x < _size; x++)
      vecDouble[x] = _vector[x];
    fwrite((void*)vecDouble, sizeof(double), _size, file);
    delete[] vecDouble;
  } 
  else
    fwrite((void*)_vector, sizeof(double), _size, file);
}

//////////////////////////////////////////////////////////////////////
// write the vector to a file
//////////////////////////////////////////////////////////////////////
void VECTOR::write(const char* filename)
{
  FILE* file;
  file = fopen(filename, "wb");
  if (file == NULL)
  {
    cout << " Could not open file " << filename << "!" << endl;
    return;
  }

  // write dimensions
  fwrite((void*)&_size, sizeof(int), 1, file);

  // write data
  if (sizeof(Real) == sizeof(float))
  {
    double* vecDouble = new double[_size];
    for (int x = 0; x < _size; x++)
      vecDouble[x] = _vector[x];
    fwrite((void*)vecDouble, sizeof(double), _size, file);
    delete[] vecDouble;
  } 
  else
    fwrite((void*)_vector, sizeof(Real), _size, file);
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read the vector to a file
//////////////////////////////////////////////////////////////////////
void VECTOR::read(FILE* file)
{
  if (_vector) free(_vector);

  // read dimensions
  fread((void*)&_size, sizeof(int), 1, file);
  _sseSize = (_size % 2) ? _size : _size + 1;
  _vector = (double*)malloc(sizeof(double) * _sseSize);

  // read data
  if (sizeof(Real) == sizeof(float))
  {
    double* vecDouble = new double[_size];
    fread((void*)vecDouble, sizeof(double), _size, file);
    for (int x = 0; x < _size; x++)
      vecDouble[x] = _vector[x];
    delete[] vecDouble;
  } 
  else
    fread((void*)_vector, sizeof(double), _size, file);
}

//////////////////////////////////////////////////////////////////////
// read vector from a file
//////////////////////////////////////////////////////////////////////
bool VECTOR::read(const char* filename)
{
  FILE* file;
  file = fopen(filename, "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << filename << " not found! " << endl;
    return false;
  }

  // read dimensions
  fread((void*)&_size, sizeof(int), 1, file);
  _sseSize = (_size % 2) ? _size : _size + 1;

  // read data
  free(_vector);
  _vector = (double*)malloc(sizeof(double) * _sseSize);

  if (sizeof(Real) == sizeof(float))
  {
    double* vecDouble = new double[_size];
    fread((void*)vecDouble, sizeof(double), _size, file);
    for (int x = 0; x < _size; x++)
      _vector[x] = vecDouble[x];
    delete[] vecDouble;
  }
  else
    fread((void*)_vector, sizeof(Real), _size, file);
  fclose(file);

  return true;
}

//////////////////////////////////////////////////////////////////////
// Deep copy equality operator
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator=(const VECTOR& m)
{
  if (m.size() != _size)
  { 
    free(_vector);
    _size = m.size();
    _sseSize = (_size % 2) ? _size : _size + 1;
    _vector = (double*)malloc(sizeof(double) * _sseSize);
  }
  __builtin_memcpy (_vector, m._vector, _size * sizeof(Real));
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Deep copy equality operator
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator=(vector<Real> m)
{
  if (m.size() != (unsigned int)_size)
  {
    free(_vector);
    _size = m.size();
    _sseSize = (_size % 2) ? _size : _size + 1;
    _vector = (double*)malloc(sizeof(double) * _sseSize);
  }
  for (unsigned int x = 0; x < m.size(); x++)
    _vector[x] = m[x];

  return *this;
}

//////////////////////////////////////////////////////////////////////
// wipe the whole vector
//////////////////////////////////////////////////////////////////////
void VECTOR::clear()
{
  __builtin_memset(_vector, 0.0, _size * sizeof(Real));
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator=(const float& m)
{
  for (int x = 0; x < _size; x++)
    _vector[x] = m;

  return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator=(const double& m)
{
  for (int x = 0; x < _size; x++)
    _vector[x] = m;

  return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator=(const int& m)
{
  for (int x = 0; x < _size; x++)
    _vector[x] = m;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// add two vectors
//////////////////////////////////////////////////////////////////////
VECTOR operator+(const VECTOR& x, const VECTOR& y) 
{
  VECTOR z(x.size());

  for (int i = 0; i < x.size(); i++)
    //z(i) = x(i) + y(i);
    z(i) = x[i] + y[i];

  return z;
}

//////////////////////////////////////////////////////////////////////
// subtract two vectors
//////////////////////////////////////////////////////////////////////
VECTOR operator-(const VECTOR& x, const VECTOR& y) 
{
  VECTOR z(x.size());

  for (int i = 0; i < x.size(); i++)
    //z(i) = x(i) - y(i);
    z(i) = x[i] - y[i];

  return z;
}

//////////////////////////////////////////////////////////////////////
// scale a vector
//////////////////////////////////////////////////////////////////////
VECTOR operator*(const VECTOR& x, const Real& scalar) 
{
  VECTOR z(x.size());

  for (int i = 0; i < x.size(); i++)
    //z(i) = x(i) * scalar;
    z(i) = x[i] * scalar;

  return z;
}

//////////////////////////////////////////////////////////////////////
// scale a vector
//////////////////////////////////////////////////////////////////////
VECTOR operator*(const Real& scalar, const VECTOR& x) 
{
  VECTOR z(x.size());

  for (int i = 0; i < x.size(); i++)
    //z(i) = x(i) * scalar;
    z(i) = x[i] * scalar;

  return z;
}

//#include "VECTOR_FAST.cpp"

#include <algorithm>
#include <memory.h>
#if __APPLE__
#include <Accelerate/Accelerate.h>
#endif

//////////////////////////////////////////////////////////////////////
// dot product with another vector
//////////////////////////////////////////////////////////////////////
Real VECTOR::operator*(const VECTOR& vector) const
{
  assert(vector._size == _size);

#if __APPLE__
#ifdef SINGLE_PRECISION
	return cblas_sdot(_size, _vector, 1, vector._vector, 1);
#else
	return cblas_ddot(_size, _vector, 1, vector._vector, 1);
#endif
#else
  Real final = 0;
  for (int x = 0; x < _size; x++)
    final += vector[x] * (*this)[x];

  return final;
#endif
}

//////////////////////////////////////////////////////////////////////
// scale vector by a constant
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator*=(const Real& alpha)
{
#if __APPLE__
#ifdef SINGLE_PRECISION
  cblas_sscal(_size, alpha, _vector, 1);
#else
  cblas_dscal(_size, alpha, _vector, 1);
#endif
#else
  for (int x = 0; x < _size; x++)
    (*this)[x] *= alpha;
#endif

  return *this;
}

//////////////////////////////////////////////////////////////////////
// Self-add operator
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator+=(const VECTOR& m)
{
  assert(m._size == _size);

#if __APPLE__
#ifdef SINGLE_PRECISION
	cblas_saxpy (_size, 1.0, m._vector, 1, _vector, 1);
#else
	cblas_daxpy (_size, 1.0, m._vector, 1, _vector, 1);
#endif
#else
  for (int x = 0; x < _size; x++)
    (*this)[x] += m[x];
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Self-subtract operator
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator-=(const VECTOR& m)
{
  assert(m._size == _size);

#if __APPLE__
#ifdef SINGLE_PRECISION
	cblas_saxpy (_size, -1.0, m._vector, 1, _vector, 1);
#else
	cblas_daxpy (_size, -1.0, m._vector, 1, _vector, 1);
#endif
#else
  for (int x = 0; x < _size; x++)
    (*this)[x] -= m[x];
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// compute the 2 norm
//////////////////////////////////////////////////////////////////////
Real VECTOR::norm2() const
{
#if __APPLE__
#ifdef SINGLE_PRECISION
	return sqrt(cblas_sdot(_size, _vector, 1, _vector, 1));
#else
	return sqrt(cblas_ddot(_size, _vector, 1, _vector, 1));
#endif
#else
  return sqrt((*this) * (*this));
#endif
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: y += alpha * x, where y is this vector
//////////////////////////////////////////////////////////////////////
void VECTOR::axpy(Real alpha, const VECTOR& x)
{
  assert(_size == x._size);

#if __APPLE__
#ifdef SINGLE_PRECISION
	cblas_saxpy (_size, alpha, x._vector, 1, _vector, 1);
#else
	cblas_daxpy (_size, alpha, x._vector, 1, _vector, 1);
#endif
#else
  for (int i = 0; i < _size; i++)
    (*this)[i] += alpha * x[i];
#endif
}

//////////////////////////////////////////////////////////////////////
// out of place component-wise absolute value
//////////////////////////////////////////////////////////////////////
VECTOR VECTOR::abs() const
{
  VECTOR return_vector(_size);
  for (int x = 0; x < _size; x++) {
    return_vector[x] = std::abs(_vector[x]);
  }
  return return_vector;
}

//////////////////////////////////////////////////////////////////////
// return the minimum entry
//////////////////////////////////////////////////////////////////////
Real VECTOR::min() const
{
  Real m = 0.0;
  for (int x = 0; x < _size; x++) {
    if ((*this)[x] < m) {
      m = (*this)[x];
    }
  }
  return m;
}

//////////////////////////////////////////////////////////////////////
// return the maximum entry
//////////////////////////////////////////////////////////////////////
Real VECTOR::max() const
{
  Real m = 0.0;
  for (int x = 0; x < _size; x++) {
    if ((*this)[x] > m) {
      m = (*this)[x];
    }
  }
  return m;
}

//////////////////////////////////////////////////////////////////////
// compute the 2 norm
//////////////////////////////////////////////////////////////////////
Real operator^(const VECTOR& x, const VECTOR& y)
{
  assert(x.size() == y.size());
#if __APPLE__
#ifdef SINGLE_PRECISION
	return cblas_sdot (x.size(), x.dataConst(), 1, y.dataConst(), 1);
#else
	return cblas_ddot (x.size(), x.dataConst(), 1, y.dataConst(), 1);
#endif
#else
  return x * y;
#endif
}
