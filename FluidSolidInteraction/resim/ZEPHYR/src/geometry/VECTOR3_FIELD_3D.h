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
#ifndef VECTOR3_FIELD_3D_H
#define VECTOR3_FIELD_3D_H

#include "EIGEN.h"

#include <cmath>
#include <string>
#include <map>
#include <iostream>
#include "VECTOR.h"

#include "VEC3.h"
#include "FIELD_3D.h"

using namespace std;

class VECTOR3_FIELD_3D {
public:
  VECTOR3_FIELD_3D();
  VECTOR3_FIELD_3D(const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  ~VECTOR3_FIELD_3D();

  // accessors
  inline VEC3F& operator()(int x, int y, int z) { return _data[z * _slabSize + y * _xRes + x]; };
  const VEC3F operator()(int x, int y, int z) const { return _data[z * _slabSize + y * _xRes + x]; };
  const VEC3F operator()(const VEC3F& position) const;
  inline VEC3F& operator[](int x) { 
    assert(x >= 0);
    assert(x < _totalCells);
    return _data[x]; 
  };
  const VEC3F operator[](int x) const { 
    assert(x >= 0);
    assert(x < _totalCells);
    return _data[x]; 
  };
  VEC3F* data() { return _data; };
  VEC3F*& dataRef() { return _data; };
  const int xRes() const { return _xRes; };
  const int yRes() const { return _yRes; };
  const int zRes() const { return _zRes; };
  VEC3F dims() { return VEC3F(_xRes, _yRes, _zRes); };
  const Real dx() const { return _dx; };
  const Real dy() const { return _dy; };
  const Real dz() const { return _dz; };
  VEC3F dxs() const { return VEC3F(_dx, _dy, _dz); };
  const int slabSize() const { return _slabSize; };
  const VEC3F center() const { return _center; };
  const VEC3F lengths() const { return _lengths; };
  const int totalCells() const { return _totalCells; };
  const bool initialized() const { return _initialized; };
  const VEC3F constEntry(int index) const { return _data[index]; };
  void clear();
 
  void setCenter(const VEC3F& center) { _center = center; };

  // what's the maximum resolution in any direction?
  int maxRes();

  // retrieve the components
  FIELD_3D magnitudeField() const;

  // absolute max single entry
  Real maxAbsScalar();

  // 2 norm of the whole field
  Real twoNorm();

  // overloaded operators
  VECTOR3_FIELD_3D& operator-=(const VECTOR3_FIELD_3D& input);
  VECTOR3_FIELD_3D& operator+=(const VECTOR3_FIELD_3D& input);
  VECTOR3_FIELD_3D& operator+=(const Real& value);
  VECTOR3_FIELD_3D& operator*=(const Real& value);
  VECTOR3_FIELD_3D& operator*=(const VEC3F& value);
  VECTOR3_FIELD_3D& operator*=(const FIELD_3D& input);
  VECTOR3_FIELD_3D& operator=(const Real& value);
  VECTOR3_FIELD_3D& operator=(const VECTOR3_FIELD_3D& input);

  // extend some vector quantity off of a front, given a signed distance function
  void fastExtension(const FIELD_3D& signedDistance);

  // file stream IO
  void write(FILE* file) const;
  void writeGz(gzFile& file) const;
  void readGz(gzFile& file);

  // advect using first order semi-Lagrangian
  static void advect(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField);
  static void advect(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const VECTOR3_FIELD_3D& oldField, VECTOR3_FIELD_3D& newField);

  // normalize all the vectors in the field
  void normalizeToLargest();

  // set various components
  void setZeroX();
  void setZeroY();
  void setZeroZ();
  void setZeroBorder();

  void setNeumannX();  
  void setNeumannY();  
  void setNeumannZ();

  void copyBorderX();
  void copyBorderY();
  void copyBorderZ();

  // BLAS-like interface, output += alpha * input
  void axpy(const Real& alpha, const VECTOR3_FIELD_3D& field);

  // swap the contents with another object
  void swapPointers(VECTOR3_FIELD_3D& field);

  // return a flattened array of all the field contents
  VECTOR flattened();
  VectorXd flattenedEigen();

  // peel off the outer boundary of grid cells in preparation for PCA
  VECTOR3_FIELD_3D peelBoundary() const;

  // set the field innards to a peeled version
  void setWithPeeled(const VectorXd& data);

  // do a projection of the peeled field, using the passed in basis
  VectorXd peeledProject(const MatrixXd& U);

  // unproject the reduced coordinate into the peeled cells in this field
  void peeledUnproject(const MatrixXd& U, const VectorXd& q);

  // take the dot product of the current field with another vector field
  // and return the scalar field
  FIELD_3D dot(const VECTOR3_FIELD_3D& rhs);

private:
  int _xRes;
  int _yRes;
  int _zRes;

  int _slabSize;
  int _totalCells;

  VEC3F* _data;

  // center position of the grid
  VEC3F _center;

  // lengths of the x,y,z dimensions of the grid
  VEC3F _lengths;

  // physical lengths
  Real _dx;
  Real _dy;
  Real _dz;

  // has this field been allocated?
  bool _initialized;
};

// take the field dot product
FIELD_3D operator*(const VECTOR3_FIELD_3D&u, const VECTOR3_FIELD_3D& v);
VECTOR3_FIELD_3D operator*(const Real& a, const VECTOR3_FIELD_3D& v);
VECTOR3_FIELD_3D operator*(const VECTOR3_FIELD_3D& v, const Real& a);
VECTOR3_FIELD_3D operator*(const FIELD_3D& u, const VECTOR3_FIELD_3D& v);
VECTOR3_FIELD_3D operator+(const VECTOR3_FIELD_3D& u, const VECTOR3_FIELD_3D& v);

// diff two vector fields
VECTOR3_FIELD_3D operator-(const VECTOR3_FIELD_3D& u, const VECTOR3_FIELD_3D& v);

#endif
