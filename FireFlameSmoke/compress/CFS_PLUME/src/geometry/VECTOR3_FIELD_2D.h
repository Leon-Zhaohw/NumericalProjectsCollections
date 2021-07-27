#ifndef VECTOR3_FIELD_2D_H
#define VECTOR3_FIELD_2D_H
#include "EIGEN.h"

#include <cmath>
#include <string>
#include <map>
#include <iostream>

#include <SETTINGS.h>
#include "VEC3.h"
#include "FIELD_2D.h"

using namespace std;

class VECTOR3_FIELD_2D {
public:
  VECTOR3_FIELD_2D();
  VECTOR3_FIELD_2D(const int& xRes, const int& yRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  VECTOR3_FIELD_2D(double* data, const int& xRes, const int& yRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  VECTOR3_FIELD_2D(float* xData, float* yData, const int& xRes, const int& yRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  VECTOR3_FIELD_2D(const VECTOR3_FIELD_2D& m);
  VECTOR3_FIELD_2D(const FIELD_2D& m);
  ~VECTOR3_FIELD_2D();

  // accessors
  inline VEC3F& operator()(int x, int y) { return _data[y * _xRes + x]; };
  const VEC3F operator()(int x, int y) const { return _data[y * _xRes + x]; };
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
  const Real dx() const { return _dx; };
  const Real dy() const { return _dy; };
  const VEC3F center() const { return _center; };
  const VEC3F lengths() const { return _lengths; };
  const int totalCells() const { return _totalCells; };
  const VEC3F constEntry(int index) const { return _data[index]; };
  void clear();
 
  // overloaded operators
  VECTOR3_FIELD_2D& operator=(const VECTOR3_FIELD_2D& input);

  // what's the maximum resolution in any direction?
  int maxRes();

  // retrieve the components
  FIELD_2D scalarField(int component) const;
  FIELD_2D magnitudeField() const;

  // check if any entry is a nan
  bool isNan();

  // real-valued cell center coordinates
  VEC3F cellCenter(int x, int y) const;

  // normalize all the vectors in the field
  void normalize();

  // draw to GL
  void draw();

  // write a LIC iamge
  // http://www.zhanpingliu.org/research/flowvis/LIC/MiniLIC.c
  void writeLIC(const char* filename);
  
  void writeLIC(int scaleUp, const char* filename);

private:
  int _xRes;
  int _yRes;

  int _totalCells;

  VEC3F* _data;

  // center position of the grid
  VEC3F _center;

  // lengths of the x,y,z dimensions of the grid
  VEC3F _lengths;

  bool _initialized;

  // physical lengths
  Real _dx;
  Real _dy;
};

#endif
