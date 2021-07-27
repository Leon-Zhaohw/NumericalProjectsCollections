//////////////////////////////////////////////////////////////////////
// This file is part of Closest Point Turbulence.
// 
// Closest Point Turbulence is free software: you can redistribute it 
// and/or modify it under the terms of the GNU General Public License 
// as published by the Free Software Foundation, either version 3 of 
// the License, or (at your option) any later version.
// 
// Closest Point Turbulence is distributed in the hope that it will 
// be useful, but WITHOUT ANY WARRANTY; without even the implied 
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Closest Point Turbulence. 
// If not, see <http://www.gnu.org/licenses/>.
// 
// Copyright 2013 Theodore Kim and Nils Thuerey
//////////////////////////////////////////////////////////////////////

#ifndef FIELD_3D_H
#define FIELD_3D_H

#include <cmath>
#include <string>
#include <map>
#include <iostream>
#include <vector>

#include "VEC3.h"
#include "FIELD_2D.h"

using namespace std;

class FIELD_3D {
public:
  FIELD_3D();
  FIELD_3D(const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  FIELD_3D(const FIELD_3D& m);
  ~FIELD_3D();

  // accessors
  const int xRes() const { return _xRes; };
  const int yRes() const { return _yRes; };
  const int zRes() const { return _zRes; };
  const VEC3F center() const { return _center; };
  const VEC3F lengths() const { return _lengths; };
  VEC3F dxs() const { return VEC3F(_dx, _dy, _dz); };
  const Real dx() const { return _dx; };
  const Real dy() const { return _dy; };
  const Real dz() const { return _dz; };
  Real* data() { return _data; };
  const int slabSize() const { return _slabSize; };

  // hack to make Houdini work
  void setDx(const Real newDx) { _dx = newDx; };

  // vector of the field's dimesions
  VEC3F dims() const { return VEC3F(_xRes, _yRes, _zRes); };

  // overloaded accessor operators
  const Real operator()(const VEC3F& position) const;
  inline Real& operator()(int x, int y, int z) { 
    assert(z >= 0 && z < _zRes && y >= 0 && y < _yRes && x >= 0 && x < _xRes);
    return _data[z * _slabSize + y * _xRes + x]; 
  };
  const Real operator()(int x, int y, int z) const { 
    assert(z >= 0 && z < _zRes && y >= 0 && y < _yRes && x >= 0 && x < _xRes);
    return _data[z * _slabSize + y * _xRes + x]; 
  };
  inline Real& operator[](int x) { return _data[x]; };
  const Real operator[](int x) const { return _data[x]; };
  
  // various order lookup functions
  Real quarticLookup(const VEC3F& position) const;
  Real cubicLookupUnclamped(const VEC3F& position) const;
  Real nearestNeighborLookup(const VEC3F& position) const;

  // do a cubic Hermite interpolation, but turn off monotonic clamping
  static Real cubicInterpUnclamped(const Real interp, const Real* points);

  void clear();
  const bool initialized() const;
  
  // real-valued cell center coordinates
  VEC3F cellCenter(int x, int y, int z) const;
  
  // overloaded operators
  FIELD_3D& operator=(const Real& alpha);
  FIELD_3D& operator=(const FIELD_3D& A);
  FIELD_3D& operator*=(const Real& alpha);
  FIELD_3D& operator-=(const FIELD_3D& input);
  FIELD_3D& operator+=(const FIELD_3D& input);
  FIELD_3D& operator*=(const FIELD_3D& input);

  // IO functions
  void writeGz(string filename) const;
  void readHoudini12Surf(string filename);
  void writeGz(gzFile& file) const;
  void readPhysBAMGz(const char* filename);

  void resizeAndWipe(int xRes, int yRes, int zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));

  // return the projection of the field in different directions
  FIELD_2D zProjection();
  
  // what's the maximum resolution in any direction?
  int maxRes();

  // convolve this field with a smaller field
  FIELD_3D convolveNarrowBandFast15(const FIELD_3D& filter, const FIELD_3D& distance, const Real maxRadius);

  // compute the narrow band indices for this object, which is assumed to be a SDF
  // the returned vector is (x,y,z) triplets
  vector<int> computeNarrowBand(Real maxCellDistance) const;

  // field minimum
  Real fieldMin();

  // field maximum
  Real fieldMax();
  
  // field maximum cell index
  VEC3F maxIndex();

  // field minimum cell index
  VEC3F minIndex();

  // copy values out into the border, assuming that "borderSize" is the width of the grid padding
  void copyIntoBorder(int borderSize);

  // first order spatial derivatives
  // on the border, difference falls back to first order (not centered) difference
  inline Real Dx(int x, int y, int z) const;
  inline Real Dy(int x, int y, int z) const;
  inline Real Dz(int x, int y, int z) const;
  
  // get the curvature
  void principalCurvatures(FIELD_3D& minCurvature, FIELD_3D& maxCurvature) const;

  // clamp the field to a min and max
  void clamp(const Real minValue, const Real maxValue);

  // get a resampled version
  FIELD_3D resampleCubicUnclamped(int xRes, int yRes, int zRes) const;
  FIELD_3D resampleCubicUnclampedNarrowBand(int upResFactor, const FIELD_3D& distanceField, const int maxRadius) const;

  // get the normal at a point
  VEC3F normal(int x, int y, int z) const;

  // get the derivative of the normal at a point
  MATRIX3 Dnormal(int x, int y, int z) const;

  // do a cubic Hermite interpolation
  static Real cubicInterp(const Real interp, const Real* points);
  
  // do a cubic Hermite that clamps to the immediate neighborhood
  static Real cubicInterpClamped(const Real interp, const Real* points);

  // do a quartic WENO interpolation
  static Real quarticInterp(const Real interp, const Real* points);

  // do a soft bandpass where there's a gradual falloff
  void softBandPass(const Real band, const Real falloff);

  // pass back a field with a new padding of size "paddingSize"
  FIELD_3D withAddedPadding(int paddingSize) const;

  // stomp the border to zero
  void stompBorder(int borderSize);

  // set everything in the specified interval to a given value
  void setInterval(const int xMin, const int xMax, const int yMin, const int yMax, const int zMin, const int zMax, const Real value = 0);

  // stomp everything outside a narrow band to zero
  void stompOutsideNarrowBand(const FIELD_3D& distance, const int maxRadius);

private:
  int _xRes;
  int _yRes;
  int _zRes;

  int _slabSize;
  int _totalCells;

  Real* _data;

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

#endif
