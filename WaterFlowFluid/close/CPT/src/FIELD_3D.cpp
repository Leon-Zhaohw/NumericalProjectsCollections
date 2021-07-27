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

#include "FIELD_3D.h"
#include <omp.h>
#include <zlib.h>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D::FIELD_3D(const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
 
  try {
    _data = new Real[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(Real);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;
}

FIELD_3D::FIELD_3D(const FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths())
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  try {
    _data = new Real[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(Real);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  for (int x = 0; x < _totalCells; x++)
    _data[x] = m[x];
}

FIELD_3D::FIELD_3D() :
  _xRes(-1), _yRes(-1), _zRes(-1), _totalCells(-1), _data(NULL)
{
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D::~FIELD_3D()
{
  if (_data)
    delete[] _data;
}
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clear()
{
  TIMER functionTimer(__FUNCTION__);

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        _data[x + y * _xRes + z * _slabSize] = 0.0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::writeGz(string filename) const
{
  // make sure it's named gz
  int size = filename.size();
  if (filename[size - 1] != 'z' || filename[size - 2] != 'g')
    filename = filename + string(".gz");

  gzFile file;
  file = gzopen(filename.c_str(), "wb1"); 
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Writing file " << filename.c_str() << " ... "; flush(cout);

  writeGz(file);

  gzclose(file);

  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void FIELD_3D::resizeAndWipe(int xRes, int yRes, int zRes, const VEC3F& center, const VEC3F& lengths)
{
  if (_xRes == xRes && _yRes == yRes && _zRes == zRes)
  {
    _center = center;
    _lengths = lengths;
    clear();

    _dx = _lengths[0] / _xRes;
    _dy = _lengths[1] / _yRes;
    _dz = _lengths[2] / _zRes;
    _invDx = 1.0 / _dx;
    _invDy = 1.0 / _dy;
    _invDz = 1.0 / _dz;
    return;
  }

  if (_data) delete[] _data;

  _xRes = xRes;
  _yRes = yRes;
  _zRes = zRes;
  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  _center = center;
  _lengths = lengths;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  try {
    _data = new Real[_totalCells];
  }
  catch(std::bad_alloc& exc)
  {
    cout << " Failed to allocate " << _xRes << " " << _yRes << " " << _zRes << " FIELD_3D!" << endl;
    double bytes = _xRes * _yRes * _zRes * sizeof(Real);
    cout <<  bytes / pow(2.0,20.0) << " MB needed" << endl;
    exit(0);
  }

  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator=(const Real& alpha)
{
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] = alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator*=(const Real& alpha)
{
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator-=(const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator+=(const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] += input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator*=(const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator=(const FIELD_3D& A)
{
  resizeAndWipe(A.xRes(), A.yRes(), A.zRes(), A.center(), A.lengths());

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] = A[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
// real-valued cell center coordinates
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_3D::cellCenter(int x, int y, int z) const
{
  VEC3F halfLengths = (Real)0.5 * _lengths;

  // set it to the lower corner
  VEC3F final = _center - halfLengths;

  // displace to the NNN corner
  final[0] += x * _dx;
  final[1] += y * _dy;
  final[2] += z * _dz;

  // displace it to the cell center
  final[0] += _dx * 0.5;
  final[1] += _dy * 0.5;
  final[2] += _dz * 0.5;

  return final;
}

///////////////////////////////////////////////////////////////////////
// what's the maximum resolution in any direction?
///////////////////////////////////////////////////////////////////////
int FIELD_3D::maxRes()
{
  int final = _xRes;
  if (_yRes > final) final = _yRes;
  if (_zRes > final) final = _zRes;
  return final;
}

///////////////////////////////////////////////////////////////////////
// lookup value at some real-valued spatial position
///////////////////////////////////////////////////////////////////////
const Real FIELD_3D::operator()(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  VEC3F corner = _center - (Real)0.5 * _lengths;
  VEC3F dxs(_dx, _dy, _dz);
  corner += (Real)0.5 * dxs;

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1    = x0 + 1;
  int y0 = (int)positionCopy[1];
  int y1    = y0 + 1;
  int z0 = (int)positionCopy[2];
  int z1    = z0 + 1;

  // clamp everything
  x0 = (x0 < 0) ? 0 : x0;
  y0 = (y0 < 0) ? 0 : y0;
  z0 = (z0 < 0) ? 0 : z0;
  
  x1 = (x1 < 0) ? 0 : x1;
  y1 = (y1 < 0) ? 0 : y1;
  z1 = (z1 < 0) ? 0 : z1;

  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;

  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  // get interpolation weights
  const float s1 = positionCopy[0]- x0;
  const float s0 = 1.0f - s1;
  const float t1 = positionCopy[1]- y0;
  const float t0 = 1.0f - t1;
  const float u1 = positionCopy[2]- z0;
  const float u0 = 1.0f - u1;

  const int i000 = x0 + y0 * _xRes + z0 * _slabSize;
  const int i010 = x0 + y1 * _xRes + z0 * _slabSize;
  const int i100 = x1 + y0 * _xRes + z0 * _slabSize;
  const int i110 = x1 + y1 * _xRes + z0 * _slabSize;
  const int i001 = x0 + y0 * _xRes + z1 * _slabSize;
  const int i011 = x0 + y1 * _xRes + z1 * _slabSize;
  const int i101 = x1 + y0 * _xRes + z1 * _slabSize;
  const int i111 = x1 + y1 * _xRes + z1 * _slabSize;

  // interpolate
  // (indices could be computed once)
  return u0 * (s0 * (t0 * _data[i000] + t1 * _data[i010]) +
               s1 * (t0 * _data[i100] + t1 * _data[i110])) +
         u1 * (s0 * (t0 * _data[i001] + t1 * _data[i011]) +
               s1 * (t0 * _data[i101] + t1 * _data[i111]));
}

///////////////////////////////////////////////////////////////////////
// load a PhysBAM level set
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readPhysBAMGz(const char* filename)
{
  // level set contains
  //
  // counts (TV_INT) vector<int, 3>
  // domain (RANGE) probably 6 floats
  // mac_offset (T) single floating point, 0 or 0.5
  //
  // then is reads in a scalar array, which contains
  //
  // length2 - an int, which seems to always equal 1
  // domain (RANGE<TV>) not TV_INT, float ranges, again?
  // the entries - scalars, 
  // everything appears to be single precision

  gzFile file;
  file = gzopen(filename, "rb1");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Could not open " << filename << "!" << endl;
    exit(0);
  }
  cout << " Reading file " << filename << " ... "; flush(cout);

  // counts - grid resolutions without padding
  gzread(file, (void*)&_xRes, sizeof(int));
  gzread(file, (void*)&_yRes, sizeof(int));
  gzread(file, (void*)&_zRes, sizeof(int));

  cout << " Resolution: " << _xRes << " " << _yRes << " " << _zRes << " "; flush(cout);

  // domain
  float xMin, yMin, zMin;
  float xMax, yMax, zMax;
  gzread(file, (void*)&xMin, sizeof(float));
  gzread(file, (void*)&yMin, sizeof(float));
  gzread(file, (void*)&zMin, sizeof(float));
  gzread(file, (void*)&xMax, sizeof(float));
  gzread(file, (void*)&yMax, sizeof(float));
  gzread(file, (void*)&zMax, sizeof(float));

  // MAC offset
  float macOffset;
  gzread(file, (void*)&macOffset, sizeof(float));

  // length2
  int length2;
  gzread(file, (void*)&length2, sizeof(int));

  cout << " MacOff: " << macOffset <<", len2 "<< length2 << "\n"; flush(cout);

  // domain (grid resolutions with padding)
  int xPaddedMin, xPaddedMax;
  int yPaddedMin, yPaddedMax;
  int zPaddedMin, zPaddedMax;
  gzread(file, (void*)&xPaddedMin, sizeof(int));
  gzread(file, (void*)&xPaddedMax, sizeof(int));
  gzread(file, (void*)&yPaddedMin, sizeof(int));
  gzread(file, (void*)&yPaddedMax, sizeof(int));
  gzread(file, (void*)&zPaddedMin, sizeof(int));
  gzread(file, (void*)&zPaddedMax, sizeof(int));

  _xRes = (xPaddedMax - xPaddedMin) + 1;
  _yRes = (yPaddedMax - yPaddedMin) + 1;
  _zRes = (zPaddedMax - zPaddedMin) + 1;

  int temp = _xRes;
  _xRes = _zRes;
  _zRes = temp;

  //cout << " Padded dims: " << endl;
  //cout << xPaddedMin << " " << xPaddedMax << endl;
  //cout << yPaddedMin << " " << yPaddedMax << endl;
  //cout << zPaddedMin << " " << zPaddedMax << endl;

  //_lengths[0] = _xRes;
  //_lengths[1] = _yRes;
  //_lengths[2] = _zRes;
  _lengths[0] = 1;
  _lengths[1] = 1;
  _lengths[2] = 1;

  // set longest dimension to 1
  int biggest = (_xRes > _yRes) ? _xRes : _yRes;
  biggest = (_zRes > biggest) ? _zRes : biggest;
  _lengths[0] = (Real)_xRes / biggest;
  _lengths[1] = (Real)_yRes / biggest;
  _lengths[2] = (Real)_zRes / biggest;

  cout << " lengths: " << _lengths << endl;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  cout << " dx: " << _dx << " dy: " << _dy << " dz: " << _dz << endl;

  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;

  if (_data) delete[] _data;
  _data = new Real[_totalCells];
  
  //cout << " Reading in level set data ... "; flush(cout);


  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        float single;
        gzread(file, (void*)&single, sizeof(float));
        (*this)(x,y,z) = single;
      }

  cout << " done." << endl;
  gzclose(file);
}

///////////////////////////////////////////////////////////////////////
// triquartic interpolation lookup
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::quarticLookup(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  const VEC3F corner = _center - (Real)0.5 * _lengths + (Real)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const float xInterp = positionCopy[0] - x1;
  const float yInterp = positionCopy[1] - y1;
  const float zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const Real p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const Real p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const Real p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const Real p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const Real p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const Real p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const Real p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const Real p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const Real p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const Real p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const Real p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const Real p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const Real p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const Real p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const Real p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const Real p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const Real z0Points[] = {quarticInterp(xInterp, p0), quarticInterp(xInterp, p1), quarticInterp(xInterp, p2), quarticInterp(xInterp, p3)};
  const Real z1Points[] = {quarticInterp(xInterp, p4), quarticInterp(xInterp, p5), quarticInterp(xInterp, p6), quarticInterp(xInterp, p7)};
  const Real z2Points[] = {quarticInterp(xInterp, p8), quarticInterp(xInterp, p9), quarticInterp(xInterp, p10), quarticInterp(xInterp, p11)};
  const Real z3Points[] = {quarticInterp(xInterp, p12), quarticInterp(xInterp, p13), quarticInterp(xInterp, p14), quarticInterp(xInterp, p15)};

  const Real finalPoints[] = {quarticInterp(yInterp, z0Points), quarticInterp(yInterp, z1Points), quarticInterp(yInterp, z2Points), quarticInterp(yInterp, z3Points)};

  return quarticInterp(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// clamp to the nearest neighbor
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::nearestNeighborLookup(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  const VEC3F corner = _center - (Real)0.5 * _lengths + (Real)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  int x0 = (int)positionCopy[0];
  int x1 = x0 + 1;

  int y0 = (int)positionCopy[1];
  int y1 = y0 + 1;
  
  int z0 = (int)positionCopy[2];
  int z1 = z0 + 1;

  const float xInterp = positionCopy[0] - x0;
  const float yInterp = positionCopy[1] - y0;
  const float zInterp = positionCopy[2] - z0;

  x0 = (x0 < 0) ? 0 : x0;
  x0 = (x0 > _xRes - 1) ? _xRes - 1 : x0;
  x1 = (x1 < 0) ? 0 : x1;
  x1 = (x1 > _xRes - 1) ? _xRes - 1 : x1;

  y0 = (y0 < 0) ? 0 : y0;
  y0 = (y0 > _yRes - 1) ? _yRes - 1 : y0;
  y1 = (y1 < 0) ? 0 : y1;
  y1 = (y1 > _yRes - 1) ? _yRes - 1 : y1;

  z0 = (z0 < 0) ? 0 : z0;
  z0 = (z0 > _zRes - 1) ? _zRes - 1 : z0;
  z1 = (z1 < 0) ? 0 : z1;
  z1 = (z1 > _zRes - 1) ? _zRes - 1 : z1;

  int x = (xInterp < 0.5) ? x0 : x1;
  int y = (yInterp < 0.5) ? y0 : y1;
  int z = (zInterp < 0.5) ? z0 : z1;

  return (*this)(x,y,z);
}

///////////////////////////////////////////////////////////////////////
// tricubic interpolation lookup
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::cubicLookupUnclamped(const VEC3F& position) const
{
  VEC3F positionCopy = position;

  // get the lower corner position
  const VEC3F corner = _center - (Real)0.5 * _lengths + (Real)0.5 * dxs();

  // recenter position
  positionCopy -= corner;

  positionCopy[0] *= _invDx;
  positionCopy[1] *= _invDy;
  positionCopy[2] *= _invDz;

  const int x1 = (int)positionCopy[0];
  const int x2    = x1 + 1;
  const int x3    = x1 + 2;
  const int x0    = x1 - 1;

  const int y1 = (int)positionCopy[1];
  const int y2    = y1 + 1;
  const int y3    = y1 + 2;
  const int y0    = y1 - 1;
  
  const int z1 = (int)positionCopy[2];
  const int z2    = z1 + 1;
  const int z3    = z1 + 2;
  const int z0    = z1 - 1;

  if (x0 < 0 || y0 < 0 || z0 < 0 ||
      x3 >= _xRes || y3 >= _yRes || z3 >= _zRes)
    return (*this)(position);

  const float xInterp = positionCopy[0] - x1;
  const float yInterp = positionCopy[1] - y1;
  const float zInterp = positionCopy[2] - z1;

  const int z0Slab = z0 * _slabSize;
  const int z1Slab = z1 * _slabSize;
  const int z2Slab = z2 * _slabSize;
  const int z3Slab = z3 * _slabSize;
  
  const int y0x = y0 * _xRes;
  const int y1x = y1 * _xRes;
  const int y2x = y2 * _xRes;
  const int y3x = y3 * _xRes;

  const int y0z0 = y0x + z0Slab;
  const int y1z0 = y1x + z0Slab;
  const int y2z0 = y2x + z0Slab;
  const int y3z0 = y3x + z0Slab;

  const int y0z1 = y0x + z1Slab;
  const int y1z1 = y1x + z1Slab;
  const int y2z1 = y2x + z1Slab;
  const int y3z1 = y3x + z1Slab;

  const int y0z2 = y0x + z2Slab;
  const int y1z2 = y1x + z2Slab;
  const int y2z2 = y2x + z2Slab;
  const int y3z2 = y3x + z2Slab;

  const int y0z3 = y0x + z3Slab;
  const int y1z3 = y1x + z3Slab;
  const int y2z3 = y2x + z3Slab;
  const int y3z3 = y3x + z3Slab;

  // do the z0 slice
  const Real p0[] = {_data[x0 + y0z0], _data[x1 + y0z0], _data[x2 + y0z0], _data[x3 + y0z0]};
  const Real p1[] = {_data[x0 + y1z0], _data[x1 + y1z0], _data[x2 + y1z0], _data[x3 + y1z0]};
  const Real p2[] = {_data[x0 + y2z0], _data[x1 + y2z0], _data[x2 + y2z0], _data[x3 + y2z0]};
  const Real p3[] = {_data[x0 + y3z0], _data[x1 + y3z0], _data[x2 + y3z0], _data[x3 + y3z0]};

  // do the z1 slice
  const Real p4[] = {_data[x0 + y0z1], _data[x1 + y0z1], _data[x2 + y0z1], _data[x3 + y0z1]};
  const Real p5[]= {_data[x0 + y1z1], _data[x1 + y1z1], _data[x2 + y1z1], _data[x3 + y1z1]};
  const Real p6[] = {_data[x0 + y2z1], _data[x1 + y2z1], _data[x2 + y2z1], _data[x3 + y2z1]};
  const Real p7[] = {_data[x0 + y3z1], _data[x1 + y3z1], _data[x2 + y3z1], _data[x3 + y3z1]};

  // do the z2 slice
  const Real p8[] = {_data[x0 + y0z2], _data[x1 + y0z2], _data[x2 + y0z2], _data[x3 + y0z2]};
  const Real p9[] = {_data[x0 + y1z2], _data[x1 + y1z2], _data[x2 + y1z2], _data[x3 + y1z2]};
  const Real p10[] = {_data[x0 + y2z2], _data[x1 + y2z2], _data[x2 + y2z2], _data[x3 + y2z2]};
  const Real p11[] = {_data[x0 + y3z2], _data[x1 + y3z2], _data[x2 + y3z2], _data[x3 + y3z2]};

  // do the z3 slice
  const Real p12[] = {_data[x0 + y0z3], _data[x1 + y0z3], _data[x2 + y0z3], _data[x3 + y0z3]};
  const Real p13[] = {_data[x0 + y1z3], _data[x1 + y1z3], _data[x2 + y1z3], _data[x3 + y1z3]};
  const Real p14[] = {_data[x0 + y2z3], _data[x1 + y2z3], _data[x2 + y2z3], _data[x3 + y2z3]};
  const Real p15[] = {_data[x0 + y3z3], _data[x1 + y3z3], _data[x2 + y3z3], _data[x3 + y3z3]};
  
  const Real z0Points[] = {cubicInterpUnclamped(xInterp, p0), cubicInterp(xInterp, p1), cubicInterp(xInterp, p2), cubicInterp(xInterp, p3)};
  const Real z1Points[] = {cubicInterpUnclamped(xInterp, p4), cubicInterp(xInterp, p5), cubicInterp(xInterp, p6), cubicInterp(xInterp, p7)};
  const Real z2Points[] = {cubicInterpUnclamped(xInterp, p8), cubicInterp(xInterp, p9), cubicInterp(xInterp, p10), cubicInterp(xInterp, p11)};
  const Real z3Points[] = {cubicInterpUnclamped(xInterp, p12), cubicInterp(xInterp, p13), cubicInterp(xInterp, p14), cubicInterp(xInterp, p15)};

  const Real finalPoints[] = {cubicInterpUnclamped(yInterp, z0Points), cubicInterp(yInterp, z1Points), cubicInterp(yInterp, z2Points), cubicInterp(yInterp, z3Points)};

  return cubicInterpUnclamped(zInterp, finalPoints);
}

///////////////////////////////////////////////////////////////////////
// do a cubic Hermite interpolation
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::cubicInterp(const Real interp, const Real* points)
{
  Real d0 = (points[2] - points[0]) * 0.5;
  Real d1 = (points[3] - points[1]) * 0.5;

  Real deltak = (points[2] - points[1]);

  // do monotonic interpolation
  if (deltak * d0 < 0.0)
    d0 = 0;
  if (deltak * d1 < 0.0)
    d1 = 0;

  Real a0 = points[1];
  Real a1 = d0;
  Real a2 = 3.0 * deltak - 2.0 * d0 - d1;
  Real a3 = -2.0 * deltak + d0 + d1;

  Real squared = interp * interp;
  Real cubed = squared * interp;
  return a3 * cubed + a2 * squared + a1 * interp + a0;
}

///////////////////////////////////////////////////////////////////////
// do a cubic Hermite interpolation
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::cubicInterpUnclamped(const Real interp, const Real* points)
{
  Real d0 = (points[2] - points[0]) * 0.5;
  Real d1 = (points[3] - points[1]) * 0.5;

  Real deltak = (points[2] - points[1]);

  Real a0 = points[1];
  Real a1 = d0;
  Real a2 = 3.0 * deltak - 2.0 * d0 - d1;
  Real a3 = -2.0 * deltak + d0 + d1;

  Real squared = interp * interp;
  Real cubed = squared * interp;
  return a3 * cubed + a2 * squared + a1 * interp + a0;
}

///////////////////////////////////////////////////////////////////////
// do a quartic WENO interpolation
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::quarticInterp(const Real interp, const Real* points)
{
  const Real& fim1 = points[0];
  const Real& fi   = points[1];
  const Real& fip1 = points[2];
  const Real& fip2 = points[3];
  const Real& x = interp;
  const Real xHalf = interp * 0.5;

  const Real p1 = fi + ((fip1 - fim1) + (fip1 - 2.0f * fi + fim1) * x) * xHalf;
  const Real p2 = fi + ((-fip2 + 4.0f * fip1 - 3.0f * fi) + (fip2 - 2.0f * fip1 + fi) * x) * xHalf;

  const Real third = 1.0 / 3.0;
  const Real C1 = (2.0f - x) * third;
  const Real C2 = (x + 1.0f) * third;

  const Real middle = -76.0f * fip1 * fi;
  const Real fip1Sq = fip1 * fip1;
  const Real fiSq = fi * fi;

  const Real twelfth = 1.0 / 12.0;
  const Real IS1 = (26.0f * fip1 * fim1 - 52.0f * fi * fim1 + middle + 25.0f * fip1Sq + 64.0f * fiSq + 13.0f * fim1 * fim1) * twelfth + 1e-6f;
  const Real IS2 = (26.0f * fip2 * fi - 52.0f * fip2 * fip1 + middle + 25.0f * fiSq + 64.0f * fip1Sq + 13.0f * fip2 * fip2) * twelfth + 1e-6f;

  const Real alpha1 = C1 / (IS1 * IS1);
  const Real alpha2 = C2 / (IS2 * IS2);

  return (alpha1 * p1 + alpha2 * p2) / (alpha1 + alpha2);
}

#define FILTERENTRY(X) const int filterIndex##X = filterPartial + X; \
                       const int dataIndex##X   = dataPartial + X; \
                       finalEntry += filter[filterIndex##X] * _data[dataIndex##X];

///////////////////////////////////////////////////////////////////////
// convolve this field with a smaller field
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::convolveNarrowBandFast15(const FIELD_3D& filter, const FIELD_3D& distance, const Real maxRadius)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D final(*this);
  final = 0;
  const Real invDx = 1.0 / distance.dx();

  assert(filter.xRes() < _xRes);
  assert(filter.yRes() < _yRes);
  assert(filter.zRes() < _zRes);

  assert(filter.xRes() % 2);
  assert(filter.yRes() % 2);
  assert(filter.zRes() % 2);

  const int xHalf = filter.xRes() / 2;
  const int yHalf = filter.yRes() / 2;
  const int zHalf = filter.zRes() / 2;

  const int xRes = _xRes;
  const int yRes = _yRes;
  const int zRes = _zRes;

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        Real currentDistance = fabs(distance(x,y,z) * invDx);
        if (currentDistance >= maxRadius) continue;
        if (x < xHalf || y < yHalf || z < zHalf || 
            x > _xRes - xHalf - 1 || y > _yRes - yHalf - 1 || z > _zRes - zHalf - 1) continue;

        float finalEntry = 0;
        const int dataPartialTop = (z - zHalf) * _slabSize - xHalf + (y - yHalf) * _xRes + x;

        for (int fz = 0; fz < filter.zRes(); fz++)
        {
          const int filterPartialZ = fz * 225;
          const int dataPartialZ = fz * _slabSize  + dataPartialTop;
          for (int fy = 0; fy < filter.yRes(); fy++)
          {
            const int filterPartial = fy * 15 + filterPartialZ;
            const int dataPartial = fy * _xRes + dataPartialZ;

            FILTERENTRY(0);
            FILTERENTRY(1);
            FILTERENTRY(2);
            FILTERENTRY(3);
            FILTERENTRY(4);
            FILTERENTRY(5);
            FILTERENTRY(6);
            FILTERENTRY(7);
            FILTERENTRY(8);
            FILTERENTRY(9);
            FILTERENTRY(10);
            FILTERENTRY(11);
            FILTERENTRY(12);
            FILTERENTRY(13);
            FILTERENTRY(14);
          }
        }

        const int finalIndex = x + y * _xRes + z * _slabSize;
        final[finalIndex] = finalEntry;
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// stomp all the value outside a narrow band to zero in order to 
// boost compression
///////////////////////////////////////////////////////////////////////
void FIELD_3D::stompOutsideNarrowBand(const FIELD_3D& distance, const int maxCells)
{
  TIMER functionTimer(__FUNCTION__);
  const Real invDx = 1.0 / distance.dx();

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        Real currentDistance = fabs(distance(x,y,z) * invDx);

        if (currentDistance > maxCells)
          (*this)(x,y,z) = 0;

      }
}

///////////////////////////////////////////////////////////////////////
// return the projection of the filter in Z direction
///////////////////////////////////////////////////////////////////////
FIELD_2D FIELD_3D::zProjection()
{
  FIELD_2D final(_xRes, _yRes);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x,y) += (*this)(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::fieldMax()
{
  assert(_totalCells > 0);

  Real final = _data[0];

  for (int x = 0; x < _totalCells; x++)
    if (_data[x] > final)
      final = _data[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_3D::maxIndex()
{
  Real maxFound = _data[0];

  VEC3F maxFoundIndex;
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        if (_data[index] > maxFound)
        {
          maxFound = _data[index];

          maxFoundIndex[0] = x;
          maxFoundIndex[1] = y;
          maxFoundIndex[2] = z;
        }

  return maxFoundIndex;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::fieldMin()
{
  assert(_totalCells > 0);

  Real final = _data[0];

  for (int x = 0; x < _totalCells; x++)
    if (_data[x] < final)
      final = _data[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_3D::minIndex()
{
  Real minFound = _data[0];

  VEC3F minFoundIndex;
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        if (_data[index] < minFound)
        {
          minFound = _data[index];

          minFoundIndex[0] = x;
          minFoundIndex[1] = y;
          minFoundIndex[2] = z;
        }

  return minFoundIndex;
}

///////////////////////////////////////////////////////////////////////
// This is hard-coded for the paddle example
///////////////////////////////////////////////////////////////////////
void FIELD_3D::readHoudini12Surf(string filename)
{
  cout << " Reading Houdini 12 distance file " << filename.c_str() << " ..."; flush(cout);
  int size = filename.size();
  if (filename[size - 1] == 'z' && filename[size - 2] == 'g')
  {
    cout << " File is Gzipped! " << endl;
    exit(0);
  }
  
  // MAGIC NUMBER - just for the paddle example
	_lengths[0] = 2.5;
	_lengths[1] = 2.5;
	_lengths[2] = 2.5;
	//_lengths[0] = 1;
	//_lengths[1] = 1;
	//_lengths[2] = 1;

  if (_data) delete[] _data;
#if 1
  // MAGIC NUMBER - just for the paddle example
  _xRes = 100; 
  _yRes = 100; 
  _zRes = 100;
#else
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " USING HOUDINI 200 SETTINGS! " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

  _xRes = 200; 
  _yRes = 200; 
  _zRes = 200; 
#endif
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new Real[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  _invDx = 1.0 / _dx;
  _invDy = 1.0 / _dy;
  _invDz = 1.0 / _dz;

  FILE* file;
  file = fopen(filename.c_str(), "r");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Could not open filename " << filename.c_str() << "!!!" << endl;
    exit(0);
  }

  char buffer[512];
  bool tileFound = false;

  while (!tileFound)
  {
    fscanf(file, "%s", buffer);
    string search(buffer); 

    // look for the one that says "tiles"
    size_t findVolume = search.find("\"tiles\"");
    if (findVolume != string::npos)
      tileFound = true;
  }

  int xTotalTiles = (_xRes / 16);
  int yTotalTiles = (_yRes / 16);
  int zTotalTiles = (_zRes / 16);

  if (xTotalTiles * 16 != _xRes) xTotalTiles++;
  if (yTotalTiles * 16 != _yRes) yTotalTiles++;
  if (zTotalTiles * 16 != _zRes) zTotalTiles++;

  for (int zTile = 0; zTile < zTotalTiles; zTile++)
    for (int yTile = 0; yTile < yTotalTiles; yTile++)
      for (int xTile = 0; xTile < xTotalTiles; xTile++)
      {
        // peel off '['
        fscanf(file, "%s", buffer);

        if (buffer[0] != '[')
          break;

        // peel off 'compression'
        fscanf(file, "%s", buffer);

        // peel off whitespace and "
        char single;
        fscanf(file, " %c", &single);
        double data;
        fscanf(file, "data\",[");

        bool closingBracketFound = false;

        // read in a block
        vector<double> blockData;
        while (!closingBracketFound)
        {
          fscanf(file, "%lf", &data);
          blockData.push_back(data);

          // peek at the next token
          fscanf(file, "%c", &single);

          if (single == ']')
            closingBracketFound = true;
        }

        int xTileSize = 16;
        int yTileSize = 16;
        int zTileSize = 16;

        if (xTile == xTotalTiles - 1)
          xTileSize = _xRes % 16;
        if (yTile == yTotalTiles - 1)
          yTileSize = _yRes % 16;
        if (zTile == zTotalTiles - 1)
          zTileSize = _zRes % 16;

        int xTileOffset = xTile * 16;
        int yTileOffset = yTile * 16;
        int zTileOffset = zTile * 16;

        int index = 0;
        for (int z = 0; z < zTileSize; z++)
          for (int y = 0; y < yTileSize; y++)
            for (int x = 0; x < xTileSize; x++, index++)
              (*this)(xTileOffset + x, yTileOffset + y,zTileOffset + z) = blockData[index];

        fscanf(file, "%s", buffer);
      }
  fclose(file);

}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::Dx(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const Real right = (x < _xRes - 1) ? _data[index + 1] : _data[index];
  const Real left  = (x > 0)         ? _data[index - 1] : _data[index];
  const Real denom = (x > 0 && x < _xRes -1) ? 1.0 / (2.0 * _dx) : 1.0 / _dx;
  return (right - left) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::Dy(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const Real up   = (y < _yRes - 1) ? _data[index + _xRes] : _data[index];
  const Real down = (y > 0)         ? _data[index - _xRes] : _data[index];
  const Real denom = (y > 0 && y < _yRes -1) ? 1.0 / (2.0 * _dy) : 1.0 / _dy;
  return (up - down) * denom;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::Dz(int x, int y, int z) const
{
  assert(x >= 0);
  assert(x < _xRes);
  assert(y >= 0);
  assert(y < _yRes);
  assert(z >= 0);
  assert(z < _zRes);

  int index = x + y * _xRes + z * _slabSize;

  const Real in  = (z < _zRes - 1) ? _data[index + _slabSize] : _data[index];
  const Real out = (z > 0)         ? _data[index - _slabSize] : _data[index];
  const Real denom = (z > 0 && z < _zRes -1) ? 1.0 / (2.0 * _dz) : 1.0 / _dz;
  return (in - out) * denom;
}

///////////////////////////////////////////////////////////////////////
// clamp the field to a min and max
///////////////////////////////////////////////////////////////////////
void FIELD_3D::clamp(const Real minValue, const Real maxValue)
{
  TIMER functionTimer(__FUNCTION__);
  assert(minValue <= maxValue);

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        bool greater = (_data[index] > maxValue);
        bool lesser = (_data[index] < minValue);

        if (greater)
        {
          _data[index] = maxValue;
          continue;
        }

        if (lesser)
          _data[index] = minValue;
      }
}

///////////////////////////////////////////////////////////////////////
// get a resampled version
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::resampleCubicUnclamped(int xRes, int yRes, int zRes) const
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D final(xRes, yRes, zRes, _center, _lengths);
  const int slabSize = xRes * yRes;

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        int index = x + y * xRes + z * slabSize;
        VEC3F center = final.cellCenter(x,y,z);
        final[index] = this->cubicLookupUnclamped(center);
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// get a resampled version
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::resampleCubicUnclampedNarrowBand(int upResFactor, const FIELD_3D& distanceField, const int maxRadius) const
{
  TIMER functionTimer(__FUNCTION__);

  int xRes = upResFactor * _xRes;
  int yRes = upResFactor * _yRes;
  int zRes = upResFactor * _zRes;

  int maxRes = (xRes > yRes) ? xRes : yRes;
  maxRes = (zRes > maxRes) ? zRes : maxRes;

  FIELD_3D final(xRes, yRes, zRes, _center, _lengths);
  //const int slabSize = xRes * yRes;
  const Real invDx = 1.0 / distanceField.dx();
  const float outside = sqrt(2.0 * maxRes * maxRes); 

  final = outside;

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        const Real currentDistance = fabs(distanceField(x,y,z) * invDx);
        if (currentDistance >= maxRadius) 
          continue;

        const int xStart = x * upResFactor;
        const int yStart = y * upResFactor;
        const int zStart = z * upResFactor;
        for (int k = 0; k < upResFactor; k++)
        {
          int zUpres = zStart + k;
          for (int j = 0; j < upResFactor; j++)
          {
            int yUpres = yStart + j;
            for (int i = 0; i < upResFactor; i++)
            {
              int xUpres = xStart + i;

              VEC3F center = final.cellCenter(xUpres,yUpres,zUpres);
              final(xUpres, yUpres, zUpres) = this->cubicLookupUnclamped(center);
            }
          }
        }
      }

  return final;
}

///////////////////////////////////////////////////////////////////////
// from "Level Set Surface Editing Operators", Museth et al. 2002
// "Geometric Surface Processing via Normal Maps", Tasdizen 2003
///////////////////////////////////////////////////////////////////////
void FIELD_3D::principalCurvatures(FIELD_3D& minCurvature, FIELD_3D& maxCurvature) const
{
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        int index = x + y * _xRes + z * _slabSize;

        MATRIX3 N = Dnormal(x,y,z);
        VEC3F n = normal(x,y,z);
        MATRIX3 outer = MATRIX3::outer_product(n, n);

        MATRIX3 B = N * (MATRIX3::I() - outer);

        Real D = sqrt(B.squaredSum());
        Real H = trace(B) * 0.5;

        Real discrim = D * D * 0.5 - H * H;

        if (discrim < 0.0)
        {
          minCurvature[index] = 0;
          maxCurvature[index] = 0;
          continue;
        }

        Real root = sqrt(discrim);
        Real k1 = H + root;
        Real k2 = H - root;

        maxCurvature[index] = (fabs(k1) > fabs(k2)) ? k1 : k2;
        minCurvature[index] = (fabs(k1) > fabs(k2)) ? k2 : k1;
      }
}

///////////////////////////////////////////////////////////////////////
// get the normal at a point
///////////////////////////////////////////////////////////////////////
VEC3F FIELD_3D::normal(int x, int y, int z) const
{
  assert (x >= 0 && x < _xRes);
  assert (y >= 0 && y < _yRes);
  assert (z >= 0 && z < _zRes);

  VEC3F final;
  final[0] = Dx(x,y,z);
  final[1] = Dy(x,y,z);
  final[2] = Dz(x,y,z);

  final.normalize();

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the normal at a point
///////////////////////////////////////////////////////////////////////
MATRIX3 FIELD_3D::Dnormal(int x, int y, int z) const
{
  const VEC3F left  = (x == 0)         ? normal(x,y,z) : normal(x-1,y,z);
  const VEC3F right = (x == _xRes - 1) ? normal(x,y,z) : normal(x+1,y,z);
  const Real dx     = (x == 0 || x == _xRes - 1) ? 1.0 / _dx : 0.5 / _dx;

  const VEC3F down = (y == 0)         ? normal(x,y,z) : normal(x,y-1,z);
  const VEC3F up   = (y == _yRes - 1) ? normal(x,y,z) : normal(x,y+1,z);
  const Real dy    = (y == 0 || y == _yRes - 1) ? 1.0 / _dy : 0.5 / _dy;

  const VEC3F out = (z == 0)         ? normal(x,y,z) : normal(x,y,z-1);
  const VEC3F in  = (z == _zRes - 1) ? normal(x,y,z) : normal(x,y,z+1);
  const Real dz   = (z == 0 || z == _zRes - 1) ? 1.0 / _dz : 0.5 / _dz;

  MATRIX3 final((right - left) * dx, (up - down) * dy, (in - out) * dz);

  return final;
}

///////////////////////////////////////////////////////////////////////
// do a soft bandpass where there's a gradual cubic falloff
///////////////////////////////////////////////////////////////////////
void FIELD_3D::softBandPass(const Real band, const Real falloff)
{
  for (int x = 0; x < _totalCells; x++)
    if (fabs(fabs(_data[x]) - band) > falloff)
      _data[x] = 0;
    else
    {
      Real interp = fabs(fabs(_data[x]) - band) / falloff;
      Real squared = interp * interp;
      Real scaling = 2 * squared * interp - 3 * squared + 1;
      _data[x] *= scaling;
    }
}

///////////////////////////////////////////////////////////////////////
// see if the field has any data in it yet
///////////////////////////////////////////////////////////////////////
const bool FIELD_3D::initialized() const
{
  if (_xRes < 0 || _yRes < 0 || _zRes < 0 || _totalCells < 0 || _data == NULL)
    return false;

  return true;
}

///////////////////////////////////////////////////////////////////////
// write out a field to a file stream
///////////////////////////////////////////////////////////////////////
void FIELD_3D::writeGz(gzFile& file) const
{
  // write dimensions
  gzwrite(file, (void*)&_xRes, sizeof(int));
  gzwrite(file, (void*)&_yRes, sizeof(int));
  gzwrite(file, (void*)&_zRes, sizeof(int));
  _center.writeGz(file);
  _lengths.writeGz(file);

  // always write out as a double
  if (sizeof(Real) != sizeof(double))
  {
    double* dataDouble = new double[_totalCells];
    for (int x = 0; x < _totalCells; x++)
      dataDouble[x] = _data[x];

    for (int z = 0; z < _zRes; z++)
      gzwrite(file, (void*)(&dataDouble[z * _slabSize]), _slabSize * sizeof(double));
    
    //gzwrite(file, (void*)dataDouble, _totalCells * sizeof(double));
    delete[] dataDouble;
  }
  else
  {
    for (int z = 0; z < _zRes; z++)
      gzwrite(file, (void*)(&_data[z * _slabSize]), _slabSize * sizeof(double));
    //gzwrite(file, (void*)_data, _totalCells * sizeof(Real));
  }
}

///////////////////////////////////////////////////////////////////////
// copy values out into the border, assuming that "borderSize" is the 
// width of the grid padding
///////////////////////////////////////////////////////////////////////
void FIELD_3D::copyIntoBorder(int borderSize)
{
  TIMER functionTimer(__FUNCTION__);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (x == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x - i,y,z) = value;
        }
        if (x == _xRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x + i,y,z) = value;
        }              
        if (y == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y - i,z) = value;
        }
        if (y == _yRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y+i,z) = value;
        }              
        if (z == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z - i) = value;
        }
        if (z == _zRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z+i) = value;
        }

        if (x == borderSize && y == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x - i,y - j,z) = value;
        }
        if (x == _xRes - 1 - borderSize && y == _yRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x + i,y + j,z) = value;
        }

        // handle the corners
        if (x == borderSize && z == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x - i,y,z-j) = value;
        }
        if (x == _xRes - 1 - borderSize && z == _zRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x + i,y,z+j) = value;
        }

        if (z == borderSize && y == borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y - j,z -i) = value;
        }
        if (z == _xRes - 1 - borderSize && y == _yRes - 1 - borderSize)
        {
          Real value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y + j,z+i) = value;
        }

      }
}

///////////////////////////////////////////////////////////////////////
// pass back a field with a new padding of size "paddingSize"
///////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::withAddedPadding(int paddingSize) const
{
  // new length, with padding
  VEC3F newLengths = _lengths;
  newLengths[0] += paddingSize * 2 * _dx;
  newLengths[1] += paddingSize * 2 * _dx;
  newLengths[2] += paddingSize * 2 * _dx;

  FIELD_3D final(_xRes + 2 * paddingSize, 
                 _yRes + 2 * paddingSize, 
                 _zRes + 2 * paddingSize, _center, newLengths);

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        final(x + paddingSize,
              y + paddingSize,
              z + paddingSize) = (*this)(x,y,z);

  final.copyIntoBorder(paddingSize);

  return final;
}

///////////////////////////////////////////////////////////////////////
// stomp the border to zero
///////////////////////////////////////////////////////////////////////
void FIELD_3D::stompBorder(int borderSize)
{
  TIMER functionTimer(__FUNCTION__);
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (x < borderSize || y < borderSize || z < borderSize ||
            x >= _xRes - borderSize || y >= _yRes - borderSize || z >= _zRes - borderSize)
          (*this)(x,y,z) = 0;
      }
}

///////////////////////////////////////////////////////////////////////
// compute the narrow band indices for this object, which is assumed 
// to be a SDF
///////////////////////////////////////////////////////////////////////
vector<int> FIELD_3D::computeNarrowBand(Real maxCellDistance) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<int> triplets;
  Real invDx = 1.0 / _dx;

  // find which cells are inside the band
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        Real currentDistance = fabs( (*this)(x,y,z) * invDx);
        if (currentDistance < maxCellDistance)
        {
          triplets.push_back(x);
          triplets.push_back(y);
          triplets.push_back(z);
        }
      }

  return triplets;
}

///////////////////////////////////////////////////////////////////////
// set everything in the specified interval to a given value
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setInterval(const int xMin, const int xMax, const int yMin, const int yMax, const int zMin, const int zMax, const Real value)
{
  assert(xMin >= 0);
  assert(xMin <= _xRes);
  assert(xMax >= 0);
  assert(xMax <= _xRes);

  assert(yMin >= 0);
  assert(yMin <= _yRes);
  assert(yMax >= 0);
  assert(yMax <= _yRes);

  assert(zMin >= 0);
  assert(zMin <= _zRes);
  assert(zMax >= 0);
  assert(zMax <= _zRes);

  for (int x = xMin; x < xMax; x++)
    for (int y = yMin; y < yMax; y++)
      for (int z = zMin; z < zMax; z++)
        (*this)(x,y,z) = value;
}
