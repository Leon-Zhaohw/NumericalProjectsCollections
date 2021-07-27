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
#ifndef FIELD_3D_H
#define FIELD_3D_H

#include "EIGEN.h"

#include <cmath>
#include <string>
#include <map>
#include <iostream>

#include "VEC3.h"
#include "MATRIX3.h"
#include "VECTOR.h"

using namespace std;

class FIELD_3D {
public:
  FIELD_3D();
  FIELD_3D(const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  FIELD_3D(const FIELD_3D& m);
  FIELD_3D(const bool* data, const int& xRes, const int& yRes, const int& zRes);
  ~FIELD_3D();

  // accessors
  inline Real& operator()(int x, int y, int z) { 
    assert(z >= 0 && z < _zRes && y >= 0 && y < _yRes && x >= 0 && x < _xRes);
    return _data[z * _slabSize + y * _xRes + x]; 
  };
  const Real operator()(int x, int y, int z) const { 
    assert(z >= 0 && z < _zRes && y >= 0 && y < _yRes && x >= 0 && x < _xRes);
    return _data[z * _slabSize + y * _xRes + x]; 
  };
  const Real operator()(const VEC3F& position) const; 

  inline Real& operator[](int x) { return _data[x]; };
  const Real operator[](int x) const { return _data[x]; };
  Real* data() { return _data; };
  Real*& dataRef() { return _data; };
  const Real* dataConst() const { return _data; };
  const int xRes() const { return _xRes; };
  const int yRes() const { return _yRes; };
  const int zRes() const { return _zRes; };
  VEC3F dxs() const { return VEC3F(_dx, _dy, _dz); };
  VEC3F invDxs() const { return VEC3F(_invDx, _invDy, _invDz); };
  const Real dx() const { return _dx; };
  const Real dy() const { return _dy; };
  const Real dz() const { return _dz; };
  const Real invDx() const { return _invDx; };
  const Real invDy() const { return _invDy; };
  const Real invDz() const { return _invDz; };
  const int slabSize() const { return _slabSize; };
  const VEC3F center() const { return _center; };
  const VEC3F lengths() const { return _lengths; };
  const int totalCells() const { return _totalCells; };
  const int outside() const { return _outside; };

  const int maxRes() const { return (_xRes > _yRes) ? ((_zRes > _xRes) ? _zRes : _xRes) : ((_zRes > _yRes) ? _zRes : _yRes); };

  // reset dimensions
  void setCenter(const VEC3F& center) { _center = center; };
  
  void clear();
  
  // real-valued cell center coordinates
  VEC3F cellCenter(int x, int y, int z) const;
  
  // overloaded operators
  FIELD_3D& operator=(const Real& alpha);
  FIELD_3D& operator=(const FIELD_3D& A);
  FIELD_3D& operator*=(const Real& alpha);
  FIELD_3D& operator/=(const Real& alpha);
  FIELD_3D& operator+=(const Real& alpha);
  FIELD_3D& operator-=(const Real& alpha);
  FIELD_3D& operator-=(const FIELD_3D& input);
  FIELD_3D& operator+=(const FIELD_3D& input);
  FIELD_3D& operator*=(const FIELD_3D& input);

  // BLAS-like interface, output += alpha * input
  void axpy(const Real& alpha, const FIELD_3D& input);

  // IO functions
  void writeGz(gzFile& file) const;
  void readGz(gzFile& file);
	
  void resizeAndWipe(int xRes, int yRes, int zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));

  // what's the maximum resolution in any direction?
  int maxRes();

  // assuming that SURFACE.initializeSignedDistanceField() has been called on this field,
  // do the fast marching method
  void fastMarchingMethod();

  // extend some scalar quantity off of a front, given a signed distance function
  void fastExtension(const FIELD_3D& signedDistance);
  
  // norms
  Real max();

  // sum of the entire field
  Real sum();

  // vector of the field's dimesions
  VEC3F dims() const { return VEC3F(_xRes, _yRes, _zRes); };

  // field minimum
  Real fieldMin();

  // field maximum
  Real fieldMax();
  
  // field maximum cell index
  VEC3F maxIndex();

  // field minimum cell index
  VEC3F minIndex();

  // draw to OpenGL
  void draw() const;
  void drawBoundingBox() const;

  // copy out the boundary
  void copyBorderAll();

  // do a cubic Hermite interpolation, but turn off monotonic clamping
  static Real cubicInterpUnclamped(const Real interp, const Real* points);
 
  // set a border of size 1 to zero
  void setZeroBorder();

  // swap the contents with another object
  void swapPointers(FIELD_3D& field);

  // take the dot product with respect to another field
  Real dot(const FIELD_3D& input) const;

  // peel off the outer boundary of grid cells
  FIELD_3D peelBoundary() const;

  // return a flattened array of all the field contents
  VECTOR flattened();
  VectorXd flattenedEigen();

private:
  int _xRes;
  int _yRes;
  int _zRes;

  int _slabSize;
  int _totalCells;

  Real* _data;

  // what fast marching considers "outside"
  int _outside;

  // center position of the grid
  VEC3F _center;

  // lengths of the x,y,z dimensions of the grid
  VEC3F _lengths;

  // physical lengths
  Real _dx;
  Real _dy;
  Real _dz;

  Real _invDx;
  Real _invDy;
  Real _invDz;
};

FIELD_3D operator^(const FIELD_3D& A, const Real alpha);
FIELD_3D operator*(const FIELD_3D& A, const Real alpha);
FIELD_3D operator/(const FIELD_3D& A, const Real alpha);
FIELD_3D operator/(const FIELD_3D& A, const FIELD_3D& B);
FIELD_3D operator+(const FIELD_3D& A, const Real alpha);
FIELD_3D operator*(const Real alpha, const FIELD_3D& A);
FIELD_3D operator+(const Real alpha, const FIELD_3D& A);
FIELD_3D operator-(const Real alpha, const FIELD_3D& A);
FIELD_3D operator-(const FIELD_3D& A, const FIELD_3D& B);
FIELD_3D operator+(const FIELD_3D& A, const FIELD_3D& B);
FIELD_3D operator*(const FIELD_3D& A, const FIELD_3D& B);

#endif
