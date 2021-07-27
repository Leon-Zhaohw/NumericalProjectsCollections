/*
QUIJIBO: Source code for the paper Symposium on Geometry Processing
         2015 paper "Quaternion Julia Set Shape Optimization"
Copyright (C) 2015  Theodore Kim

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef FIELD_2D_H
#define FIELD_2D_H

#include <cmath>
#include <string>
#include <iostream>

#include <SETTINGS.h>
#include <VEC3.h>

#ifndef VARNAME
#define VARNAME(x) #x
#endif
#ifndef FIELDVIEW2D
#define FIELDVIEW2D(x) FIELD_2D::fieldViewer(x, VARNAME(x)); sleep(1);
#endif
#ifndef OVERLAYFIELDVIEW2D
#define OVERLAYFIELDVIEW2D(x,y) FIELD_2D::overlayFieldViewer(x, y, VARNAME(x)); sleep(1);
#endif

using namespace std;

class FIELD_2D {
public:
  FIELD_2D();
  FIELD_2D(const int& rows, const int& cols, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  FIELD_2D(const FIELD_2D& m);
  ~FIELD_2D();

  // accessors
  inline Real& operator()(int x, int y) { return _data[y * _xRes + x]; };
  inline Real operator()(int x, int y) const { return _data[y * _xRes + x]; };
  inline Real& operator[](int x) { return _data[x]; };
  inline Real operator[](int x) const { return _data[x]; };
  inline Real& operator()(const VEC3F& v) { return _data[cellIndex(v)]; };
  Real* data() { return _data; };
  int xRes() const { return _xRes; };
  int yRes() const { return _yRes; };
  const VEC3F& center() const { return _center; };
  const VEC3F& lengths() const { return _lengths; };
  int totalCells() const { return _totalCells; };
  const Real& dx() const { return _dx; };
  const Real& dy() const { return _dy; };
  const Real& invDx() const { return _invDx; };
  const Real& invDy() const { return _invDy; };
  VEC3F dims() const { return VEC3F(_xRes, _yRes, 0); };

  // get the real-space location of a cell
  VEC3F cellCenter(int index) const;
  VEC3F cellCenter(int x, int y) const;
  int cellIndex(const VEC3F& position) const;
  void cellIndex(const VEC3F& position, int& x, int& y) const;

  // a safe, slow, toroidal data accessor -- does all bounds checking for you
  inline Real safe(int x, int y);

  // common field operations
  void clear();
  void normalize();
  FIELD_2D& abs();

  Real min();
  Real max();

  // field maximum cell index
  VEC3F maxIndex();

  // field minimum cell index
  VEC3F minIndex();

  // take the log
  void log(Real base = 2.0);
 
  // IO functions
  void writeJPG(string filename);
  void writePPM(string filename);
  void writeMatlab(string filename, string variableName) const;
  void write(string filename) const;
  void write(FILE* file) const;
  void read(string filename);
  void read(FILE* file);
  void writePNG(string filename);
  void readPNG(string filename);

  // FFT functions
  void FFT(FIELD_2D& real, FIELD_2D& im);
  void shiftFFT();

  // set this field to the result of convolving filter and input
  void convolve(const FIELD_2D& filter, const FIELD_2D& input);

  void resizeAndWipe(int xRes, int yRes, const VEC3F& center, const VEC3F& lengths);

  // overloaded operators
  FIELD_2D& operator=(const Real& alpha);
  FIELD_2D& operator=(const FIELD_2D& A);
  FIELD_2D& operator*=(const Real& alpha);
  FIELD_2D& operator/=(const Real& alpha);
  FIELD_2D& operator+=(const Real& alpha);
  FIELD_2D& operator-=(const Real& alpha);
  FIELD_2D& operator-=(const FIELD_2D& input);
  FIELD_2D& operator+=(const FIELD_2D& input);
  FIELD_2D& operator*=(const FIELD_2D& input);
  FIELD_2D& operator/=(const FIELD_2D& input);

  // sum of all entries
  Real sum();
  
  // Compute the elements of the vertical derivative convolution kernel
  void verticalDerivativeKernel(double kMax = 10, double dk = 0.01, double sigma = 1.0, double L = 0);
  
  // Compute a radial Bessel function
  void radialBessel();
  
  // set to a bessel function
  void setToBessel(float k);

  FIELD_2D nearestNeighborUpsample(int factor);

  // set to a checkboard for debugging
  void setToCheckerboard(int xChecks = 10, int yChecks = 10);
  
  // set to a checkboard for debugging
  void setToRampedCheckerboard(int xChecks = 10, int yChecks = 10);
  
  // set to a ramp for debugging
  void setToRampX();
  
  // set to a ramp for debugging
  void setToRampY();

  // pass a field to fieldViewer2D
  static void fieldViewer(const FIELD_2D& field, string name);
  static void overlayFieldViewer(const FIELD_2D& distanceField, const FIELD_2D& field, string name);

  // get the projection of the field in the x direction
  //VECTOR projectionX();

  // return a field for the Laplacian of this field
  FIELD_2D laplacian();
  
  // return a field for the 4th order accurate Laplacian of this field
  FIELD_2D laplacian4th();
  
  // return a field for the gradient of this field
  FIELD_2D gradient();

  // return the transpose (flip x and y)
  FIELD_2D transpose() const;

  // use wavelet upsampling to double resolution
  FIELD_2D doubleRes() const;

  // get the sum of the absolute value of all the entries
  Real fabsSum();

  // get the mean curvature
  FIELD_2D meanCurvature() const;
  FIELD_2D gaussianCurvature() const;

  // first order spatial derivatives
  // on the border, difference falls back to first order (not centered) difference
  Real Dx(int x, int y) const; 
  Real Dy(int x, int y) const; 

  // second order spatial derivatives
  // on the border, center cell is copied to outside, and centered difference
  // is still taken
  Real DDx(int x, int y) const; 
  Real DDy(int x, int y) const; 

  // mixed derivatives
  // on the border, center cell is copied to outside, and centered difference
  // is still taken
  Real DDxy(int x, int y) const;

  // get the normal at a point
  VEC3F normal(int x, int y) const;

  // count the number of infs in the field
  int totalInfs() const;
  
  // count the number of NaNs in the field
  int totalNans() const;

private:
  int _xRes;
  int _yRes;
  int _totalCells;
  Real* _data;

  VEC3F _center;
  VEC3F _lengths;
  Real _dx, _dy;
  Real _invDx, _invDy;
};

FIELD_2D operator*(const FIELD_2D& A, const Real alpha);
FIELD_2D operator/(const FIELD_2D& A, const Real alpha);
FIELD_2D operator+(const FIELD_2D& A, const Real alpha);
FIELD_2D operator*(const Real alpha, const FIELD_2D& A);
FIELD_2D operator+(const Real alpha, const FIELD_2D& A);
FIELD_2D operator-(const FIELD_2D& A, const FIELD_2D& B);
FIELD_2D operator+(const FIELD_2D& A, const FIELD_2D& B);
FIELD_2D operator*(const FIELD_2D& A, const FIELD_2D& B);

#endif
