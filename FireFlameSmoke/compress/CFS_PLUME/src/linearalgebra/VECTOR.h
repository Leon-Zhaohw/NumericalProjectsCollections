/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2013 Theodore Kim
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
  VECTOR(const vector<Real>& v);
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
  void write(const string& filename) { write(filename.c_str()); };
  void write(const char* filename);
  void write(FILE* file);

  // read the vector to a binary file
  // everything is always read in as a double, then
  // converted if necessary
  bool read(const char* filename);
  void read(FILE* file);

  // resize the vector and wipe to zero
  void resizeAndWipe(int size);

  // overloaded operators
  VECTOR& operator=(const VECTOR& m);
  VECTOR& operator=(vector<Real> m);
  VECTOR& operator=(const double& m);
  VECTOR& operator=(const float& m);
  VECTOR& operator=(const int& m);
  VECTOR& operator+=(const VECTOR& v);
  VECTOR& operator-=(const VECTOR& v);
  VECTOR& operator*=(const Real& alpha);

  // 2 norm
  Real norm2() const;
  
  // dot product
  Real operator*(const VECTOR& vector) const;

  // raw data pointer
  Real* data() { return _vector; };
  const Real* dataConst() const { return _vector; };

  // BLAS axpy operation: y += alpha * x, where y is this vector
  void axpy(Real alpha, const VECTOR& x);
  
  // out of place component-wise absolute value
  VECTOR abs() const;

  // get the minimum entry
  Real min() const;

  // get the maximum entry
  Real max() const;

  static bool printVertical;

private:
  int _size;
  int _sseSize;
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
