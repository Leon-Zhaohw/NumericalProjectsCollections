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
// VECTOR.h: interface for the VECTOR class.
//
//////////////////////////////////////////////////////////////////////

#ifndef VECTOR_H
#define VECTOR_H

#include <SETTINGS.h>
#include <map>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace std;

//////////////////////////////////////////////////////////////////////
// An arbitrary dimension vector class
//////////////////////////////////////////////////////////////////////
class VECTOR {

public:
  VECTOR();
  VECTOR(int size);
  VECTOR(const char* filename);
  VECTOR(const VECTOR& v);
  VECTOR(FILE* file);
  VECTOR(const vector<Real>& v);
  VECTOR(const vector<int>& v);
  ~VECTOR();

  inline Real& operator()(int index) { 
    assert(index >= 0);
    assert(index < _size);
    return _vector[index];
  };
  inline Real& operator[](int index) { 
    assert(index >= 0);
    assert(index < _size);
    return _vector[index];
  };
  inline Real operator[] (int index) const { 
    assert(index >= 0);
    assert(index < _size);
    return _vector[index];
  };
 
  int size() const { return _size; };

  // wipe the whole vector
  void clear();

  // write the vector to a binary file
  // everything is always written as a double
  void write(const char* filename);
  void write(FILE* file);
  void writeMatlab(string filename, string varname);

  // read the vector to a binary file
  // everything is always read in as a double, then
  // converted if necessary
  bool read(const char* filename);
  void read(FILE* file);

  // resize the vector and wipe to zero
  void resizeAndWipe(int size);

  // resize the vector and copy its current contents
  void resizeAndCopy(int size);

  // overloaded operators
  VECTOR& operator=(VECTOR m);
  VECTOR& operator=(vector<Real> m);
  VECTOR& operator+=(const VECTOR& v);
  VECTOR& operator-=(const VECTOR& v);
  VECTOR& operator*=(const Real& alpha);

  // 2 norm
  Real norm2();
  
  // Infinity norm
  Real normInf();

  Real rms();

  // dot product
  Real operator*(const VECTOR& vector);

  // swap contents with another vector --
  // it's your job to ensure that they are the same size
  void swap(VECTOR& vector);

  // raw data pointer
  Real* data() { return _vector; };
  const Real* dataConst() const { return _vector; };

  // BLAS axpy operation: y += alpha * x, where y is this vector
  void axpy(Real alpha, VECTOR& x);

  // same as axpy above, but vector contents are stomped as well
  void clearingAxpy(Real alpha, VECTOR& x);

  // sum of all the elements
  Real sum();

  // squared sum of all the elements
  Real sum2();

  // in-place copy, since operator= must allocate a new VECTOR
  void copyInplace(VECTOR& vector, bool resizeIfNeeded = false);
  void equals(VECTOR& vector) { copyInplace(vector); };
  //void equals(VECTOR vector) { copyInplace(vector); };

  // mean of all the entries
  Real mean() { return sum() / _size; };

  // max of all the elements
  Real maxValue();

  // min of all the elements
  Real minValue();

  // take the absolute value of all entires
  void fabs();

  static bool printVertical;

  // print a subset of the vector to cout
  void printSubset(int start, int end);

  // get a sorted version of this vector
  VECTOR sorted();

  // clamp entries smaller than a threshold to zero
  void clampToZero(const Real threshold);

private:
  int _size;
  Real* _vector;
};

//////////////////////////////////////////////////////////////////////
// dump vector to iostream
//////////////////////////////////////////////////////////////////////
inline ostream &operator<<(ostream &out, VECTOR vector)
{
  out << "[";
  if (VECTOR::printVertical)
  {
    out << endl;
    for (int x = 0; x < vector.size(); x++)
      out << vector(x) << endl;
  }
  else
    for (int x = 0; x < vector.size(); x++)
      out << vector(x) << " ";
  out << "]";
  if (VECTOR::printVertical)
    out << endl;
  return out;
}

// overloaded operators
VECTOR operator-(const VECTOR& x, const VECTOR& y);
VECTOR operator+(const VECTOR& x, const VECTOR& y);
VECTOR operator*(const VECTOR& x, const Real& scalar);
VECTOR operator*(const Real& scalar, const VECTOR& x);

// x^T . y
Real operator^(const VECTOR& x, const VECTOR& y);

#endif
