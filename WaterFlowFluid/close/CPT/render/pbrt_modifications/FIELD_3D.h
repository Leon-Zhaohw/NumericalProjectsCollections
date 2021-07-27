#ifndef FIELD_3D_H
#define FIELD_3D_H

#include <cmath>
#include <string>
#include <map>
#include <iostream>
#include <cassert>

#include "VEC3.h"

#define Real float

using namespace std;

class FIELD_3D {
public:
  FIELD_3D();
  FIELD_3D(const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  FIELD_3D(const double* data, const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  FIELD_3D(const FIELD_3D& m);
  FIELD_3D(const char* filename);
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
  const Real* dataConst() const { return _data; };
  const int xRes() const { return _xRes; };
  const int yRes() const { return _yRes; };
  const int zRes() const { return _zRes; };
  VEC3F dxs() const { return VEC3F(_dx, _dy, _dz); };
  const Real dx() const { return _dx; };
  const Real dy() const { return _dy; };
  const Real dz() const { return _dz; };
  const int slabSize() const { return _slabSize; };
  const VEC3F center() const { return _center; };
  const VEC3F lengths() const { return _lengths; };
  const int totalCells() const { return _totalCells; };
  const int outside() const { return _outside; };
  const bool initialized() const;

  // reset dimensions
  void setCenter(const VEC3F& center) { _center = center; };
  void setLengths(const VEC3F& lengths);
  
  // interpolation lookup
  Real quarticLookup(const VEC3F& position) const;

  void clear();
  
  // real-valued cell center coordinates
  VEC3F cellCenter(int x, int y, int z) const;
  
  // cell index of a real-valued position
  int cellIndex(VEC3F& position) const;

  // overloaded operators
  FIELD_3D& operator=(const Real& alpha);
  FIELD_3D& operator=(const FIELD_3D& A);
  FIELD_3D& operator*=(const Real& alpha);
  FIELD_3D& operator/=(const Real& alpha);
  FIELD_3D& operator+=(const Real& alpha);
  FIELD_3D& operator-=(const FIELD_3D& input);
  FIELD_3D& operator+=(const FIELD_3D& input);
  FIELD_3D& operator*=(const FIELD_3D& input);

  // IO functions
  void write(string filename) const;
  void writeGz(string filename) const;
  void read(string filename);
  void readGz(string filename);

  // streaming IO functions
  void write(FILE* file) const;
  void read(FILE* file);
  void writeGz(gzFile& file) const;
  void readGz(gzFile& file);

  void resizeAndWipe(int xRes, int yRes, int zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));

  // what's the maximum resolution in any direction?
  int maxRes();

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

  // copy out the boundary
  void copyBorderAll();
  
  // do a quartic WENO interpolation
  static Real quarticInterp(const Real interp, const Real* points);
};

FIELD_3D operator^(const FIELD_3D& A, const Real alpha);
FIELD_3D operator*(const FIELD_3D& A, const Real alpha);
FIELD_3D operator/(const FIELD_3D& A, const Real alpha);
FIELD_3D operator/(const FIELD_3D& A, const FIELD_3D& B);
FIELD_3D operator+(const FIELD_3D& A, const Real alpha);
FIELD_3D operator*(const Real alpha, const FIELD_3D& A);
FIELD_3D operator+(const Real alpha, const FIELD_3D& A);
FIELD_3D operator-(const FIELD_3D& A, const FIELD_3D& B);
FIELD_3D operator+(const FIELD_3D& A, const FIELD_3D& B);
FIELD_3D operator*(const FIELD_3D& A, const FIELD_3D& B);

#endif
