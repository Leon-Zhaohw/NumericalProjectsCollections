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
#ifndef FIELD_3D_H
#define FIELD_3D_H

#include <cmath>
#include <string>
#include <map>
#include <iostream>

#include "VEC3.h"
#include "MATRIX3.h"
#include "QUATERNION.h"
#include <vector>

#ifndef VARNAME
#define VARNAME(x) #x
#endif

#ifndef FIELDVIEW3D
#define FIELDVIEW3D(x) FIELD_3D::fieldViewer(x, VARNAME(x)); sleep(1);
#endif
#ifndef OVERLAYFIELDVIEW3D
#define OVERLAYFIELDVIEW3D(x,y) FIELD_3D::overlayFieldViewer(x, y, VARNAME(x)); sleep(1);
#endif
#ifndef FIELDVIEW3DYZ
#define FIELDVIEW3DYZ(x) FIELD_3D::fieldViewerYZ(x, VARNAME(x)); sleep(1);
#endif
#ifndef OVERLAYFIELDVIEW3DYZ
#define OVERLAYFIELDVIEW3DYZ(x,y) FIELD_3D::overlayFieldViewerYZ(x, y, VARNAME(x)); sleep(1);
#endif
#ifndef VOLUMEFIELDRENDER3D
#define VOLUMEFIELDRENDER3D(x) FIELD_3D::volumeRender(x, VARNAME(x)); sleep(1);
#endif

using namespace std;

class FIELD_3D {
public:
  FIELD_3D();
  FIELD_3D(const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  FIELD_3D(const Real* data, const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  FIELD_3D(const FIELD_3D& m);
  FIELD_3D(const char* filename);

  // this is only for debugging if you want to viz a bool field
  FIELD_3D(const bool* data, const int& xRes, const int& yRes, const int& zRes);
  FIELD_3D(const unsigned char* data, const int& xRes, const int& yRes, const int& zRes);
  
  ~FIELD_3D();

  // accessors
  inline Real& operator()(int x, int y, int z) { 
    assert(z >= 0 && z < _zRes && y >= 0 && y < _yRes && x >= 0 && x < _xRes);
    return _data[z * _slabSize + y * _xRes + x]; 
  };
  Real operator()(int x, int y, int z) const { 
    assert(z >= 0 && z < _zRes && y >= 0 && y < _yRes && x >= 0 && x < _xRes);
    return _data[z * _slabSize + y * _xRes + x]; 
  };
  Real operator()(const VEC3F& position) const; 

  inline Real& operator[](int x) { return _data[x]; };
  Real operator[](int x) const { return _data[x]; };
  Real* data() { return _data; };
  Real*& dataRef() { return _data; };
  const Real* dataConst() const { return _data; };
  int xRes() const { return _xRes; };
  int yRes() const { return _yRes; };
  int zRes() const { return _zRes; };
  VEC3I res() const { return VEC3I(_xRes, _yRes, _zRes); };
  VEC3F dxs() const { return VEC3F(_dx, _dy, _dz); };
  VEC3F invDxs() const { return VEC3F(_invDx, _invDy, _invDz); };
  Real dx() const { return _dx; };
  Real dy() const { return _dy; };
  Real dz() const { return _dz; };
  Real invDx() const { return _invDx; };
  Real invDy() const { return _invDy; };
  Real invDz() const { return _invDz; };
  int slabSize() const { return _slabSize; };
  const VEC3F center() const { return _center; };
  const VEC3F lengths() const { return _lengths; };
  VEC3F& center() { return _center; };
  VEC3F& lengths() { return _lengths; };
  int totalCells() const { return _totalCells; };
  int outside() const { return _outside; };
  const int& quinticClamps() const { return _quinticClamps; };
  bool initialized() const;
  const QUATERNION rotation() const { return _rotation; };
  QUATERNION& rotation() { return _rotation; };

  // scale the length of the field, update the _dx and _invDx as well
  void scaleLengths(const Real& scale);

  int maxRes() const { return (_xRes > _yRes) ? ((_zRes > _xRes) ? _zRes : _xRes) : ((_zRes > _yRes) ? _zRes : _yRes); };

  // reset dimensions
  void setCenter(const VEC3F& center) { _center = center; };
  void setLengths(const VEC3F& lengths);

  Real maxLength() const { Real final = (_lengths[0] > _lengths[1]) ? _lengths[0] : _lengths[1]; return (final > _lengths[2]) ? final : _lengths[2]; };
  
  // tricubic interpolation lookup
  Real quinticLookup(const VEC3F& position) const;
  Real quarticLookup(const VEC3F& position) const;
  Real cubicLookup(const VEC3F& position) const;
  Real cubicLookupUnclamped(const VEC3F& position) const;
  Real cubicNewtonLookup(const VEC3F& position) const;

  void clear();
  
  // real-valued cell center coordinates
  VEC3F cellCenter(const VEC3I& index) const { return cellCenter(index[0], index[1], index[2]); };
  VEC3F cellCenter(int x, int y, int z) const;
  VEC3F cellCenter(int index) const;
  
  // cell index of a real-valued position
  int cellIndex(VEC3F& position) const;
  void cellIndex(const VEC3F& position, VEC3I& indices) const;
  void cellNeighborhood(const VEC3F& position, vector<VEC3I>& indices, vector<Real>& weights) const;

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
  static void axpy(const Real& alpha, const FIELD_3D& input, FIELD_3D& output);
  void axpy(const Real& alpha, const FIELD_3D& input);

  // IO functions
  void write(string filename) const;
  void writeGz(string filename) const;
  void writeMatlab(string filename, string variableName) const;
  void read(string filename);
  void readGz(string filename);
  void readHoudini(string filename);
  void readHoudiniSurf(string filename);
  void readHoudiniSurfGz(string filename);
  static void readHoudiniVel(const string filename, FIELD_3D& xVelocity, FIELD_3D& yVelocity, FIELD_3D& zVelocity);

  // streaming IO functions
  void write(FILE* file) const;
  void read(FILE* file);
  void writeGz(gzFile& file) const;
  void readGz(gzFile& file);
	
  // load a PhysBAM level set
  void readPhysBAM(const char* filename);
  void readPhysBAMGz(const char* filename);

  void resizeAndWipe(int xRes, int yRes, int zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));

  // do a union with 'field', assuming both this and field are signed distance fields
  FIELD_3D signedDistanceUnion(const FIELD_3D& field);

  // what's the maximum resolution in any direction?
  int maxRes();

  // extend some scalar quantity off of a front, given a signed distance function
  void fastExtension(const FIELD_3D& signedDistance);
  
  // return the indices of the grid points insides a world-space bounding box
  void boundingBoxIndices(const VEC3F& mins, const VEC3F& maxs, VEC3I& iMins, VEC3I& iMaxs);

  // norms
  Real sumSq();
  Real max();
  Real absMax();

  // build a const field with the given dims
  static FIELD_3D constField(const FIELD_3D& dims, Real value);

  // check if any entry is a nan
  bool isNan();

  // compute the inverse
  FIELD_3D inverse();

  // print the neighborhood of a cell for debugging
  void printNeighborhood(int x, int y, int z) const;
  void printNeighborhood(int index) const;

  // clamp nans to some specified value
  void clampNans(Real value = 0);

  // clamp nansto values in this field
  void clampNans(FIELD_3D& clampField);

  // clamp infinities to some specified value
  void clampInfs(Real value = 0);
  
  // clamp infinities to values in this field
  void clampInfs(FIELD_3D& clampField);
  
  // clamp infinities to mean of the field
  void clampInfsToMean();
  
  // clamp NaNs to mean of the field
  void clampNansToMean();

  // count the number of infs in the field
  int totalInfs() const;
  
  // count the number of NaNs in the field
  int totalNans() const;

  // dump to a viewer
  static void fieldViewer(const FIELD_3D& field, const string name);
  static void fieldViewerYZ(const FIELD_3D& field, const string name);
  
  static void overlayFieldViewer(const FIELD_3D& field, const FIELD_3D& distance, string name);
  static void overlayFieldViewerYZ(const FIELD_3D& field, const FIELD_3D& distance, string name);

  static void volumeRender(const FIELD_3D& field, const string name);

  // get the integer indices of a spatial position
  void indices(const VEC3F& position, int* x);

  // set to a checkerboard solid texture
  void setToSolidCheckboard(int xChecks = 10, int yChecks = 10, int zChecks = 10);
  void setToGrayCheckerboard(int xChecks = 10, int yChecks = 10, int zChecks = 10);

  // set to vertical derivative kernel
  void setToVerticalDerivativeKernel(double kMax = 10, double dk = 0.01, double sigma = 1.0, double L = 1);
  
  // set whole field to a Gaussian
  void setToGaussian(Real amplitude = 1.0, VEC3F sigmas = VEC3F(0.1, 0.1, 0.1));
  
  // set to a Gaussian
  void insertGaussian(const VEC3F& center, const Real amplitude = 1.0, const VEC3F sigmas = VEC3F(0.1, 0.1, 0.1));

  // convolve this field with a smaller field
  FIELD_3D convolve(const FIELD_3D& filter);
  FIELD_3D convolveToroidal(const FIELD_3D& filter);
  FIELD_3D convolveNarrowBand(const FIELD_3D& filter, const FIELD_3D& distance, int maxCells);

  // stomp all the value outside a narrow band to zero in order to boost compression
  void stompOutsideNarrowBand(const FIELD_3D& distance, int maxCells);

  // sum of the entire field
  Real sum() const;
  Real absSum() const;

  // vector of the field's dimesions
  VEC3F dims() const { return VEC3F(_xRes, _yRes, _zRes); };

  // determine how many non-zero entries are in the filter
  int nonZeroEntries();

  // normalize the data
  void normalize();

  // field minimum
  Real fieldMin() const;

  // field maximum
  Real fieldMax() const;
  
  // field maximum cell index
  VEC3I maxIndex() const;

  // field minimum cell index
  VEC3I minIndex() const;
  
  // field maximum real-valued cell center
  VEC3F maxCellCenter() const;

  // flip the z and y coordinates
  FIELD_3D flipZY() const;
  
  // flip the x and y coordinates
  FIELD_3D flipXY() const;

  // flip the x and z coordinates
  FIELD_3D flipXZ() const;
  
  // create a mirror image along different axes
  FIELD_3D mirrorZ() const;
  FIELD_3D mirrorY() const;
  FIELD_3D mirrorX() const;

  // unit tests
  // returns the true value, "ground", and the interpolation "computed",
  // at point "x"
  static void cubicUnitTest(Real& ground, Real& computed, const Real x = 0.5);

  // draw to OpenGL
  void draw() const;
  void drawBoundingBox() const;

  // copy out the boundary
  void copyBorderAll();

  // copy values out into the border, assuming that "borderSize" is the width of the grid padding
  void copyIntoBorder(int borderSize);

  // first order spatial derivatives
  // on the border, difference falls back to first order (not centered) difference
  Real Dx(int x, int y, int z) const; 
  Real Dy(int x, int y, int z) const; 
  Real Dz(int x, int y, int z) const; 
  Real Dx(const VEC3F& point) const; 
  Real Dy(const VEC3F& point) const; 
  Real Dz(const VEC3F& point) const; 
  
  // second order spatial derivatives
  // on the border, center cell is copied to outside, and centered difference
  // is still taken
  Real DDx(int x, int y, int z) const; 
  Real DDy(int x, int y, int z) const; 
  Real DDz(int x, int y, int z) const; 

  // mixed derivatives
  // on the border, center cell is copied to outside, and centered difference
  // is still taken
  Real DDxy(int x, int y, int z) const;
  Real DDxz(int x, int y, int z) const; 
  Real DDyz(int x, int y, int z) const;

  // get a field of the entire derivative
  FIELD_3D Dx() const;
  FIELD_3D Dy() const;
  FIELD_3D Dz() const;
  FIELD_3D DDx() const;
  FIELD_3D DDy() const;
  FIELD_3D DDz() const;
  FIELD_3D DDxy() const;
  FIELD_3D DDxz() const;
  FIELD_3D DDyz() const;

  // get the curvature
  FIELD_3D meanCurvature() const;
  FIELD_3D gaussianCurvature() const;
  FIELD_3D principalCurvature() const;
  Real principalCurvature(const VEC3F& point) const;

  // get the sum of gradients
  FIELD_3D divergence() const;
  
  // get the Laplacian
  FIELD_3D laplacian() const;

  // mask out any values past a certain distance
  void maskByDistance(const FIELD_3D& distanceField, const Real distance);

  // clamp the field to a min and max
  void clamp(const Real minValue, const Real maxValue);

  // get a resampled version
  FIELD_3D resampleCubic(int xRes, int yRes, int zRes) const;
  FIELD_3D resampleQuintic(int xRes, int yRes, int zRes) const;

  // take the square root of the field
  FIELD_3D squareRoot();

  // get the normal at a point
  VEC3F normal(int x, int y, int z) const;
  VEC3F normal(const VEC3F& point) const;

  // get the derivative of the normal at a point
  MATRIX3 Dnormal(int x, int y, int z) const;

  // do a cubic Hermite interpolation
  static Real cubicInterp(const Real interp, const Real* points);
  
  // do a cubic Hermite that clamps to the immediate neighborhood
  static Real cubicInterpClamped(const Real interp, const Real* points);

  // do a cubic Hermite interpolation, but turn off monotonic clamping
  static Real cubicInterpUnclamped(const Real interp, const Real* points);
 
  // do a quartic WENO interpolation
  static Real quarticInterp(const Real interp, const Real* points);

  // do a quintic Hermite interpolation
  static Real quinticInterp(const Real interp, const Real* points);
  
  // stomp to zero anything in the field that is not between the min and max
  void bandPass(const Real minValue, const Real maxValue);
  
  // isolate values near the current value, within a certain width
  void isolateBand(const Real target, const Real width);

  // single explicit diffusion step
  void blur(Real dt = 1.0);

  // get band-limited curvature, with some diffusion
  FIELD_3D bandLimitedCurvature(const Real target, const Real width);

  // set to the absolute value
  void absoluteValue();

  // do a soft bandpass where there's a gradual falloff
  void softBandPass(const Real band, const Real falloff);

  // pass back a field with a new padding of size "paddingSize"
  FIELD_3D withAddedPadding(int paddingSize) const;

  // stomp the border to zero
  void stompBorder(int borderSize);

  // set a border of size 1 to zero
  void setZeroBorder();

  // swap the contents with another object
  void swapPointers(FIELD_3D& field);

  // take the dot product with respect to another field
  Real dot(const FIELD_3D& input) const;

  // get a subfield of this field
  FIELD_3D subfield(const int xBegin, const int xEnd, 
                    const int yBegin, const int yEnd, 
                    const int zBegin, const int zEnd) const;

  // peel off the outer boundary of grid cells
  FIELD_3D peelBoundary() const;

  // is this point outside the bounding box?
  bool inside(const VEC3F& point) const;

  // splat a density particle to this grid
  void splat(const Real weight, const VEC3F& point, const int cellRadius = 1);

  // return a turbulence field
  static FIELD_3D turbulence(const int xRes, const int yRes, const int zRes, const int seed);
  //static FIELD_3D turbulentPlume(const int xRes, const int yRes, const int zRes, const int seed);
  static FIELD_3D yRampField(const FIELD_3D& example, const Real plumeBase = 0.0);

  // linearly downsampled -- halve the size of the field
  FIELD_3D linearHalfSample();

  // downsample multigrid-style so that half the sample still remain on the grid
  FIELD_3D multigridDownsample();

  // take the power of a field
  FIELD_3D power(const Real& exponent) const;

  // assuming that this field is a signed distance field, get the shell
  // of cells that are immediately on the surface
  FIELD_3D outerShell() const;
  FIELD_3D innerShell() const;

  // distance from a point to the surface of the bounding box
  Real distanceToBoundingBox(const VEC3F& v);

  // invert a distance field so that there is 1 at the border, and tapers off afterwards
  FIELD_3D invertedDistance() const;
  
  // invert a distance field so that there is 1 at the border, and tapers off afterwards
  // making sure inside and outside both sum to 0.5
  FIELD_3D invertedDistanceNormalized() const;

  // set all the borders to a constant
  void setBorders(const Real& value);

  // do some uniform transforms
  void translateX(const Real& value) { _center += _rotation.toRotationMatrix() * VEC3F(value,0,0); };
  void translateY(const Real& value) { _center += _rotation.toRotationMatrix() * VEC3F(0,value,0); };
  void translateZ(const Real& value) { _center += _rotation.toRotationMatrix() * VEC3F(0,0,value); };
  //void translateX(const Real& value) { _center += VEC3F(value,0,0); };
  //void translateY(const Real& value) { _center += VEC3F(0,value,0); };
  //void translateZ(const Real& value) { _center += VEC3F(0,0,value); };
  void rotateX(const Real& angle);
  void rotateY(const Real& angle);
  void rotateZ(const Real& angle);

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

  // retirement hash for fast marching
  // note that entry existence is the important things, not whether
  // the value is true or false. If the entry even exists, it was retired
  map<int, bool> _retired;

  // track how many times the quintic is clamped
  static int _quinticClamps;

  // has this field been initialized?
  //bool _initialized;

  // the current field rotation
  QUATERNION _rotation;

  // do a cubic Newton interpolation
  static Real cubicNewtonInterp(Real interp, Real* points);

  // read a Houdini field off from a file stream -- it is assumed that the
  // file is already advanced to the beginning of the field
  void readHoudiniField(FILE* file, bool storeValues);

  // a static version for reading in velocity
  static FIELD_3D readHoudiniField(FILE* file);
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
