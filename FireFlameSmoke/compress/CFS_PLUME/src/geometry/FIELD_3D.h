#ifndef FIELD_3D_H
#define FIELD_3D_H

#include "EIGEN.h"

#include <cmath>
#include <string>
#include <map>
#include <iostream>

#include "VEC3.h"
#include "MATRIX3.h"
#include "FIELD_2D.h"
#include "MIN_HEAP.h"
#include "VECTOR.h"

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

using namespace std;

class FIELD_3D {
public:
  FIELD_3D();
  FIELD_3D(const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  FIELD_3D(const double* data, const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  FIELD_3D(const FIELD_3D& m);
  FIELD_3D(const vector<FIELD_2D>& slices);
  FIELD_3D(const char* filename);

  // this is only for debugging if you want to viz a bool field
  FIELD_3D(const bool* data, const int& xRes, const int& yRes, const int& zRes);
  FIELD_3D(const unsigned char* data, const int& xRes, const int& yRes, const int& zRes);
  
  // strip the dimensions from m, but populate using data
  FIELD_3D(const VECTOR& data, const FIELD_3D& m);
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
  const int& quinticClamps() const { return _quinticClamps; };
  const bool initialized() const;
  static bool& usingFastPow() { return _usingFastPow; };

  const int maxRes() const { return (_xRes > _yRes) ? ((_zRes > _xRes) ? _zRes : _xRes) : ((_zRes > _yRes) ? _zRes : _yRes); };

  // reset dimensions
  void setCenter(const VEC3F& center) { _center = center; };
  void setLengths(const VEC3F& lengths);
  
  // tricubic interpolation lookup
  Real quinticLookup(const VEC3F& position) const;
  Real quarticLookup(const VEC3F& position) const;
  Real cubicLookup(const VEC3F& position) const;
  Real cubicLookupUnclamped(const VEC3F& position) const;
  Real cubicNewtonLookup(const VEC3F& position) const;

  void clear();
  void clear(const vector<int>& nonZeros);
  void clear(const vector<int>& nonZeros, int size);
  
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
  FIELD_3D& operator-=(const Real& alpha);
  FIELD_3D& operator^=(const Real& alpha);
  FIELD_3D& operator-=(const FIELD_3D& input);
  FIELD_3D& operator+=(const FIELD_3D& input);
  FIELD_3D& operator*=(const FIELD_3D& input);
  FIELD_3D& operator/=(const FIELD_3D& input);

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

  // return a slice in the form of a FIELD_2D
  FIELD_2D zSlice(int z) const;
  
  // return the projection of the field in different directions
  FIELD_2D zProjection();
  FIELD_2D yProjection();
  FIELD_2D xProjection();
  
  // do a union with 'field', assuming both this and field are signed distance fields
  FIELD_3D signedDistanceUnion(const FIELD_3D& field);

  // what's the maximum resolution in any direction?
  int maxRes();

  // assuming that SURFACE.initializeSignedDistanceField() has been called on this field,
  // do the fast marching method
  void fastMarchingMethod();

  // extend some scalar quantity off of a front, given a signed distance function
  void fastExtension(const FIELD_3D& signedDistance);
  
  // return the indices of the grid points insides a world-space bounding box
  void boundingBoxIndices(const VEC3F& mins, const VEC3F& maxs, VEC3I& iMins, VEC3I& iMaxs);

  // norms
  Real sumSq() const;
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

  // round each entry to nearest integer
  void roundInt();

  // round each entry to nearest integer out of place
  FIELD_3D castInt() const;

  // clamp nans to some specified value
  void clampNans(Real value = 0);

  // clamp nansto values in this field
  void clampNans(FIELD_3D& clampField);

  // clamp infinities to some specified value
  void clampInfs(Real value = 0);
  
  // clamp infinities to values in this field
  void clampInfs(FIELD_3D& clampField);

  // dump to a viewer
  static void fieldViewer(const FIELD_3D& field, const string name);
  static void fieldViewerYZ(const FIELD_3D& field, const string name);
  
  static void overlayFieldViewer(const FIELD_3D& field, const FIELD_3D& distance, string name);
  static void overlayFieldViewerYZ(const FIELD_3D& field, const FIELD_3D& distance, string name);

  // get the integer indices of a spatial position
  void indices(const VEC3F& position, int* x);
  
  // set to uniform random from 0 to 1
  void setToRandom();

  // take each element to the specified power
  void toPower(double power);

  void toFastPower(double power);

  void toPower(double power, const vector<int>& nonZeros);
  void toPower(double power, const vector<int>& nonZeros, const int size);

  // set to a checkerboard solid texture
  void setToSolidCheckboard(int xChecks = 10, int yChecks = 10, int zChecks = 10);
  void setToGrayCheckerboard(int xChecks = 10, int yChecks = 10, int zChecks = 10);

  // set to vertical derivative kernel
  void setToVerticalDerivativeKernel(double kMax = 10, double dk = 0.01, double sigma = 1.0, double L = 1);
  
  // set to Kolmogorov-Zakharov kernel
  void setToKolmogorovZakharov(double kMin = 1, double kMax = 10, double dk = 0.01);
  
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
  Real sum();

  // vector of the field's dimesions
  VEC3F dims() const { return VEC3F(_xRes, _yRes, _zRes); };

  // determine how many non-zero entries are in the filter
  int nonZeroEntries();

  // set a given z slice to the given 2D field
  void setSliceZ(const int& z, const FIELD_2D& slice);

  // normalize the data
  void normalize();

  // field minimum
  Real fieldMin();

  // field maximum
  Real fieldMax();
  
  // field maximum cell index
  VEC3F maxIndex();

  // field minimum cell index
  VEC3F minIndex();

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
  
  // pass back a field with a new uni-sided padding of size "paddingSize" in the z direction
  FIELD_3D pad_z(int paddingSize) const;
  
  // same as above but in the y direction
  FIELD_3D pad_y(int paddingSize) const;

  // same as above but in the x direction
  FIELD_3D pad_x(int paddingSize) const; 

  // do x, then y, then z
  FIELD_3D pad_xyz(const VEC3I& paddings) const;

  // pad with zeros rather than contiguous values
  FIELD_3D zeroPad_xyz(const VEC3I& paddings) const;

  // pass back a field with a new uni-sided zero padding of size "paddingSize" in the z direction
  FIELD_3D zeroPad_z(int paddingSize) const;
  
  // same as above but in the y direction
  FIELD_3D zeroPad_y(int paddingSize) const;

  // same as above but in the x direction;
  FIELD_3D zeroPad_x(int paddingSize) const;

  // stomp the border to zero
  void stompBorder(int borderSize);

  // set a border of size 1 to zero
  void setZeroBorder();

  // swap the contents with another object
  void swapPointers(FIELD_3D& field);

  // write out to an image file
  void writeImageSliceXY(int slice, string prefix, int picCnt, float scale=1.);
  void writeImageSliceYZ(int slice, string prefix, int picCnt, float scale=1.);
  void writeImageSliceXZ(int slice, string prefix, int picCnt, float scale=1.);

  // take the dot product with respect to another field
  Real dot(const FIELD_3D& input) const;

  // wavelet support
  FIELD_3D downsample() const;
  FIELD_3D upsample() const;

  // get a subfield of this field
  FIELD_3D subfield(const int xBegin, const int xEnd, 
                    const int yBegin, const int yEnd, 
                    const int zBegin, const int zEnd) const;

  // peel off the outer boundary of grid cells
  FIELD_3D peelBoundary() const;

  // return a flattened array of all the field contents
  VECTOR flattened() const;
  VECTOR flattenedRow() const;
  VectorXd flattenedEigen() const;

  // set the field innards to a peeled version
  void setWithPeeled(const VECTOR& data);
  void setWithPeeled(const VectorXd& data);

  // is this point outside the bounding box?
  bool inside(const VEC3F& point) const;

  // splat a density particle to this grid
  void splat(const Real weight, const VEC3F& point);

  // return a turbulence field
  static FIELD_3D turbulence(const int xRes, const int yRes, const int zRes, const int seed);
  //static FIELD_3D turbulentPlume(const int xRes, const int yRes, const int zRes, const int seed);
  static FIELD_3D yRampField(const FIELD_3D& example, const Real plumeBase = 0.0);

  // from here: 
  // http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/

  
  static inline double fastPow(double a, double b) 
  {
    //printf(" calling fast pow. ");
    union {
      double d;
      int x[2];
    } u = { a };
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    return u.d;
  }
  
  // write out to a PBRT file for rendering
  static void exportPbrt(const FIELD_3D& density, const char* out);

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

  static bool _usingFastPow;

  // do fast marching in one direction
  void marchOneway(bool forward, MIN_HEAP& minHeap);
  
  void marchOneway2ndOrder(bool forward, MIN_HEAP& minHeap);

  // insert the front in preparation for reinitialization or extension
  void insertFront(const bool forward, FIELD_3D& distance, MIN_HEAP& minHeap);
  
  // do fast extension in one direction
  void extendOneway(bool forward, FIELD_3D& distance, MIN_HEAP& minHeap);

  // do a cubic Newton interpolation
  static Real cubicNewtonInterp(Real interp, Real* points);

  // read a Houdini field off from a file stream -- it is assumed that the
  // file is already advanced to the beginning of the field
  void readHoudiniField(FILE* file, bool storeValues);

  // a static version for reading in velocity
  static FIELD_3D readHoudiniField(FILE* file);
  
  // helper to write out to an image file
  void writeProjectedIntern(int dir1, int dir2, string prefix, int picCnt, float scale=1.); 


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
