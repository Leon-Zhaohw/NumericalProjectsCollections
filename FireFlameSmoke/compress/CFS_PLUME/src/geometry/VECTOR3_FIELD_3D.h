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
#include "VECTOR3_FIELD_2D.h"

#ifndef VECTORFIELDVIEW3D
#define VECTORFIELDVIEW3D(x,y) VECTOR3_FIELD_3D::fieldViewer(x, y, VARNAME(x)); sleep(1);
#endif

#ifndef GRADIENTVIEW3D
#define GRADIENTVIEW3D(x,y) VECTOR3_FIELD_3D::gradientViewer(x, y, VARNAME(x)); sleep(1);
#endif

#ifndef OVERLAYVECTORFIELDVIEW3D
#define OVERLAYVECTORFIELDVIEW3D(x,y) VECTOR3_FIELD_3D::overlayFieldViewer(x, y, VARNAME(x)); sleep(1);
#endif

#ifndef CLOSESTPOINTVIEW3D
#define CLOSESTPOINTVIEW3D(x,y) VECTOR3_FIELD_3D::closestPointViewer(x, y, VARNAME(x)); sleep(1);
#endif
using namespace std;

class VECTOR3_FIELD_3D {
public:
  VECTOR3_FIELD_3D();
  VECTOR3_FIELD_3D(const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  VECTOR3_FIELD_3D(double* data, const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  VECTOR3_FIELD_3D(const VECTOR& data, const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  VECTOR3_FIELD_3D(const VectorXd& data, const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  // VECTOR3_FIELD_3D(float* xData, float* yData, float* zData, const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  VECTOR3_FIELD_3D(Real* xData, Real* yData, Real* zData, const int& xRes, const int& yRes, const int& zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));
  VECTOR3_FIELD_3D(const VECTOR3_FIELD_3D& m);
  VECTOR3_FIELD_3D(const FIELD_3D& m);

  // strip the dimensions from m, but populate using data
  VECTOR3_FIELD_3D(const VECTOR& data, const VECTOR3_FIELD_3D& m);
  VECTOR3_FIELD_3D(const VectorXd& data, const VECTOR3_FIELD_3D& m);
  ~VECTOR3_FIELD_3D();

  // reallocate everything
  void resizeAndWipe(int xRes, int yRes, int zRes, const VEC3F& center = VEC3F(0,0,0), const VEC3F& lengths = VEC3F(1,1,1));

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
 
  // reset the lengths to something else, and recompute all the dimesions as well
  void setLengths(const VEC3F& lengths);
  
  void setCenter(const VEC3F& center) { _center = center; };

  // what's the maximum resolution in any direction?
  int maxRes();

  // create a field of the grid positions of the passed in grid
  static VECTOR3_FIELD_3D cellCenters(const FIELD_3D& input);
  
  // take the gradient of a scalar field
  static VECTOR3_FIELD_3D gradient(const FIELD_3D& input);

  // return a grid of values at the given spatial positions
  static FIELD_3D compose(const FIELD_3D& values, const VECTOR3_FIELD_3D& positions);
  static VECTOR3_FIELD_3D compose(const VECTOR3_FIELD_3D& values, const VECTOR3_FIELD_3D& positions);

  // retrieve the components
  FIELD_3D scalarField(int component) const;
  FIELD_3D magnitudeField() const;

  // norms
  Real sumMagnitudes();
  Real maxMagnitude();
  Real maxMagnitude(int& index);

  // check if any entry is a nan
  bool isNan();

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

  // set the field values to uniform random doubles from 0 to 1
  void setToRandom();

  // set the values in the field to the values at the closest points
  static FIELD_3D setToClosestPointValues(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints);
  static FIELD_3D setToClosestPointValuesNarrowBand(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints, const FIELD_3D& distance, int maxCells);

  // real-valued cell center coordinates
  VEC3F cellCenter(int x, int y, int z) const;

  // file IO
  void write(const string filename) const;
  void read(const string filename);
  void readPhysBAMGz(const string filename);
  void readPhysBAMGz(const FIELD_3D& distanceField, const string filename);
  void readHoudini(const string filename);
  void readHoudiniGz(const string filename);

  // file stream IO
  void write(FILE* file) const;
  void read(FILE* file);
  void writeGz(gzFile& file) const;
  void readGz(gzFile& file);

  // draw to GL, with the cell centers of 'field' as the vector origins
  void draw(const VECTOR3_FIELD_3D& origins) const;
  void drawZSlice(const VECTOR3_FIELD_3D& origins, const int zSlice, const Real scale = 1, const int stride = 4) const;
  void drawClosestPointZSlice(const VECTOR3_FIELD_3D& origins, const int zSlice, const Real scale = 1, const int stride = 4) const;
  void drawBoundingBox();

  // dump to a viewer
  static void fieldViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name);
  static void gradientViewer(const VECTOR3_FIELD_3D& field, const VECTOR3_FIELD_3D& gradientField, string name);
  static void overlayFieldViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name);
  static void closestPointViewer(const VECTOR3_FIELD_3D& field, const FIELD_3D& distanceField, string name);

  // flip the z and y coordinates
  VECTOR3_FIELD_3D flipZY() const;
  
  // flip the x and y coordinates
  VECTOR3_FIELD_3D flipXY() const;
  
  // flip the x and z coordinates
  VECTOR3_FIELD_3D flipXZ() const;

  // advect using first order semi-Lagrangian
  static void advect(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField);
  static void advect(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const VECTOR3_FIELD_3D& oldField, VECTOR3_FIELD_3D& newField);

  // advect with second order MacCormack
  static void advectMacCormack(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, FIELD_3D& oldField, FIELD_3D& newField, FIELD_3D& temp1, FIELD_3D& temp2);
  static void advectMacCormack(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, VECTOR3_FIELD_3D& oldField, VECTOR3_FIELD_3D& newField, VECTOR3_FIELD_3D& temp1, VECTOR3_FIELD_3D& temp2);

  // version of the () operator for debugging purposes
  VEC3F debugPositionOperator(const VEC3F& position) const;

  // normalize all the vectors in the field
  void normalize();

  // normalize all the vectors in the field
  void normalizeToLargest();

  // copy values out into the border, assuming that "borderSize" is the width of the grid padding
  void copyIntoBorder(int borderSize);

  // pass back a field with a new padding of size "paddingSize"
  VECTOR3_FIELD_3D withAddedPadding(int paddingSize) const;

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
  void copyBorderAll();
  
  void setZeroSphere(const VEC3I& center, double radius);

  // BLAS-like interface, output += alpha * input
  void axpy(const Real& alpha, const VECTOR3_FIELD_3D& field);

  // swap the contents with another object
  void swapPointers(VECTOR3_FIELD_3D& field);

  // return a flattened array of all the field contents
  VECTOR flattened() const;
  VectorXd flattenedEigen() const;
  VectorXd flattenedEigenHomogeneous() const;

  // get gradient field
  VECTOR3_FIELD_3D gradient();

  // peel off the outer boundary of grid cells in preparation for PCA
  VECTOR3_FIELD_3D peelBoundary() const;

  // peel off and return a single boundary
  VECTOR3_FIELD_3D leftBoundary() const;
  VECTOR3_FIELD_3D rightBoundary() const;
  VECTOR3_FIELD_3D topBoundary() const;
  VECTOR3_FIELD_3D bottomBoundary() const;
  VECTOR3_FIELD_3D nearBoundary() const;
  VECTOR3_FIELD_3D farBoundary() const;

  // retrieve the boundaries in an array
  void getPeeledBoundaries(vector<VECTOR3_FIELD_3D>& boundaries) const;

  // reassemble a velocity field from a bunch of peeled boundaries and a center
  void reassemblePeels(const VECTOR3_FIELD_3D& middle, const vector<VECTOR3_FIELD_3D>& boundaries);

  // set the field innards to a peeled version
  void setWithPeeled(const VECTOR& data);
  void setWithPeeled(const VectorXd& data);

  // do a projection of the peeled field, using the passed in basis
  VectorXd peeledProject(const MatrixXd& U);

  // unproject the reduced coordinate into the peeled cells in this field
  void peeledUnproject(const MatrixXd& U, const VectorXd& q);

  // higher order interpolations
  VEC3F cubicLookup(const VEC3F& position) const;
  VEC3F quarticLookup(const VEC3F& position) const;

  // take the dot product of the current field with another vector field
  // and return the scalar field
  FIELD_3D dot(const VECTOR3_FIELD_3D& rhs);

  // take the curl to get a divergence free velocity field
  static VECTOR3_FIELD_3D curl(const FIELD_3D& scalar1, const FIELD_3D& scalar2, const FIELD_3D& scalar3);
  static VECTOR3_FIELD_3D curl(const VECTOR3_FIELD_3D& potential);
  static VECTOR3_FIELD_3D noisePotential(const FIELD_3D& example);

  // get a z slice
  VECTOR3_FIELD_2D zSlice(const int z) const;

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
  
  // retirement hash for fast extension
  // note that entry existence is the important things, not whether
  // the value is true or false. If the entry even exists, it was retired
  map<int, bool> _retired;
  
  // insert the front in preparation for reinitialization or extension
  void insertFront(const bool forward, FIELD_3D& distance, MIN_HEAP& minHeap);
  
  // do fast extension in one direction
  void extendOneway(bool forward, FIELD_3D& distance, MIN_HEAP& minHeap);

  // Clamp the extrema generated by the BFECC error correction
  static void clampExtrema(const Real dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, FIELD_3D& newField);
  static void clampOutsideRays(const Real dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, const FIELD_3D& oldAdvection, FIELD_3D& newField);
  static void clampExtrema(const Real dt, const VECTOR3_FIELD_3D& velocityField, const VECTOR3_FIELD_3D& oldField, VECTOR3_FIELD_3D& newField);
  static void clampOutsideRays(const Real dt, const VECTOR3_FIELD_3D& velocityField, const VECTOR3_FIELD_3D& oldField, const VECTOR3_FIELD_3D& oldAdvection, VECTOR3_FIELD_3D& newField);

  // do a cubic Hermite interpolation
  static Real cubicInterp(const Real interp, const Real* points);
  
  // do a quartic WENO interpolation
  static Real quarticInterp(const Real interp, const Real* points);
};

// take the field dot product
FIELD_3D operator*(const VECTOR3_FIELD_3D&u, const VECTOR3_FIELD_3D& v);
VECTOR3_FIELD_3D operator*(const Real& a, const VECTOR3_FIELD_3D& v);
VECTOR3_FIELD_3D operator*(const VECTOR3_FIELD_3D& v, const Real& a);
VECTOR3_FIELD_3D operator*(const FIELD_3D& u, const VECTOR3_FIELD_3D& v);
//VECTOR3_FIELD_3D operator*(const VECTOR3_FIELD_3D& u, const FIELD_3D& v);
VECTOR3_FIELD_3D operator+(const VECTOR3_FIELD_3D& u, const VECTOR3_FIELD_3D& v);

// diff two vector fields
VECTOR3_FIELD_3D operator-(const VECTOR3_FIELD_3D& u, const VECTOR3_FIELD_3D& v);

#endif
