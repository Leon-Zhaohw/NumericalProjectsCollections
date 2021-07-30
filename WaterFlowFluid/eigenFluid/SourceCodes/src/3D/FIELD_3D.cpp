
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
#include "3D/FIELD_3D.h"
#include "Eigen"
#include <omp.h>
// #include <zlib.h>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "setting.h"

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

  _outside = maxRes() * maxRes();
}

FIELD_3D::FIELD_3D(const bool* data, const int& xRes, const int& yRes, const int& zRes) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(0,0,0), _lengths(1,1,1)
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

  _outside = maxRes() * maxRes();

  for (int x = 0; x < _totalCells; x++)
    _data[x] = data[x];
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

  _outside = maxRes() * maxRes();
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
  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0.0;
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
  _outside = maxRes() * maxRes();
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
  for (int x = 0; x < _totalCells; x++)
    _data[x] = alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator*=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator/=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] /= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator+=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] += alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator-=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator-=(const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);
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

  for (int x = 0; x < _totalCells; x++)
    _data[x] *= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator^(const FIELD_3D& A, const Real alpha)
{
  FIELD_3D final(A);

  for (int x = 0; x < final.totalCells(); x++)
    final[x] = pow(final[x], alpha);

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const FIELD_3D& A, const Real alpha)
{
  FIELD_3D final(A);
  final *= alpha;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator/(const FIELD_3D& A, const FIELD_3D& B)
{
  FIELD_3D final(A);

  for (int x = 0; x < A.totalCells(); x++)
    final[x] *= 1.0 / B[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator/(const FIELD_3D& A, const Real alpha)
{
  FIELD_3D final(A);
  final /= alpha;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator-(const FIELD_3D& A, const FIELD_3D& B)
{
  assert(A.xRes() == B.xRes());
  assert(A.yRes() == B.yRes());
  assert(A.zRes() == B.zRes());

  FIELD_3D final(A);
  final -= B;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator+(const FIELD_3D& A, const FIELD_3D& B)
{
  FIELD_3D final(A);
  final += B;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const FIELD_3D& A, const FIELD_3D& B)
{
  FIELD_3D final(A);
  final *= B;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator+(const FIELD_3D& A, const Real alpha)
{
  FIELD_3D final(A);
  final += alpha;
  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const Real alpha, const FIELD_3D& A)
{
  return A * alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator+(const Real alpha, const FIELD_3D& A)
{
  return A + alpha;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D operator-(const Real alpha, const FIELD_3D& A)
{
  FIELD_3D final(A);

  for (int x = 0; x < final.totalCells(); x++)
    final[x] = alpha - final[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D& FIELD_3D::operator=(const FIELD_3D& A)
{
  resizeAndWipe(A.xRes(), A.yRes(), A.zRes(), A.center(), A.lengths());

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
// copy out the boundary
///////////////////////////////////////////////////////////////////////
void FIELD_3D::copyBorderAll()
{
  int index;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      // front slab
      index = x + y * _xRes;
      _data[index] = _data[index + _slabSize];

      // back slab
      index += _totalCells - _slabSize;
      _data[index] = _data[index - _slabSize];
    }

  for (int z = 0; z < _zRes; z++)
    for (int x = 0; x < _xRes; x++)
    {
      // bottom slab
      index = x + z * _slabSize;
      _data[index] = _data[index + _xRes];

      // top slab
      index += _slabSize - _xRes;
      _data[index] = _data[index - _xRes];
    }

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
    {
      // left slab
      index = y * _xRes + z * _slabSize;
      _data[index] = _data[index + 1];

      // right slab
      index += _xRes - 1;
      _data[index] = _data[index - 1];
    }
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
// maximum entry
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::max()
{
  Real final = _data[0];

  for (int i = 0; i < _totalCells; i++)
    if (_data[i] > final)
      final = _data[i];

  return final;
}

///////////////////////////////////////////////////////////////////////
// get the sum of the field
///////////////////////////////////////////////////////////////////////
Real FIELD_3D::sum()
{
  Real final = 0;
  for (int x = 0; x < _totalCells; x++)
    final += _data[x];

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
// draw to OpenGL
///////////////////////////////////////////////////////////////////////
void FIELD_3D::draw(Real amp) const
{
  Real strideDx = _dx * 0.5;
  Real strideDy = _dy * 0.5;
  Real strideDz = _dz * 0.5;

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        VEC3F center = cellCenter(x,y,z);

        glPushMatrix();
          VEC3F v000 = center + VEC3F(-strideDx, -strideDy, -strideDz);
          VEC3F v001 = center + VEC3F(-strideDx, -strideDy, strideDz);
          VEC3F v010 = center + VEC3F(-strideDx, strideDy, -strideDz);
          VEC3F v011 = center + VEC3F(-strideDx, strideDy, strideDz);
          VEC3F v100 = center + VEC3F(strideDx, -strideDy, -strideDz);
          VEC3F v101 = center + VEC3F(strideDx, -strideDy, strideDz);
          VEC3F v110 = center + VEC3F(strideDx, strideDy, -strideDz);
          VEC3F v111 = center + VEC3F(strideDx, strideDy, strideDz);

          //if (fabs(_data[index]) < 1e-5)
          if (fabs(_data[index])*amp < 1e-3)
          {
            glPopMatrix();
            continue;
          }
          
          glColor4f(0,0,0, fabs(_data[index] * amp));

          glBegin(GL_QUADS);
            // face 1
            glNormal3f(0,0,-1);
#if SINGLE_PRECISION
            glVertex3fv(v000);
            glVertex3fv(v001);
            glVertex3fv(v011);
            glVertex3fv(v010);
#else
            glVertex3dv(v000);
            glVertex3dv(v001);
            glVertex3dv(v011);
            glVertex3dv(v010);
#endif

            // face 2
            glNormal3f(0,0,1);
#if SINGLE_PRECISION
            glVertex3fv(v110);
            glVertex3fv(v111);
            glVertex3fv(v101);
            glVertex3fv(v100);
#else
            glVertex3dv(v110);
            glVertex3dv(v111);
            glVertex3dv(v101);
            glVertex3dv(v100);
#endif

            // face 3
            glNormal3f(-1,0,0);
#if SINGLE_PRECISION
            glVertex3fv(v010);
            glVertex3fv(v110);
            glVertex3fv(v100);
            glVertex3fv(v000);
#else
            glVertex3dv(v010);
            glVertex3dv(v110);
            glVertex3dv(v100);
            glVertex3dv(v000);
#endif

            // face 4
            glNormal3f(1,0,0);
#if SINGLE_PRECISION
            glVertex3fv(v111);
            glVertex3fv(v011);
            glVertex3fv(v001);
            glVertex3fv(v101);
#else
            glVertex3dv(v111);
            glVertex3dv(v011);
            glVertex3dv(v001);
            glVertex3dv(v101);
#endif

            // face 5
            glNormal3f(0,1,0);
#if SINGLE_PRECISION
            glVertex3fv(v111);
            glVertex3fv(v110);
            glVertex3fv(v010);
            glVertex3fv(v011);
#else
            glVertex3dv(v111);
            glVertex3dv(v110);
            glVertex3dv(v010);
            glVertex3dv(v011);
#endif

            // face 6
            glNormal3f(0,-1,0);
#if SINGLE_PRECISION
            glVertex3fv(v001);
            glVertex3fv(v000);
            glVertex3fv(v100);
            glVertex3fv(v101);
#else
            glVertex3dv(v001);
            glVertex3dv(v000);
            glVertex3dv(v100);
            glVertex3dv(v101);
#endif
          glEnd();
        glPopMatrix();
      }
}

///////////////////////////////////////////////////////////////////////
// draw bounding box to OpenGL
///////////////////////////////////////////////////////////////////////
void FIELD_3D::drawBoundingBox() const
{
  glLineWidth(5.0);

  glPushMatrix();
    VEC3F v000 = cellCenter(0,0,0);
    VEC3F v001 = cellCenter(0,0,_zRes);
    VEC3F v010 = cellCenter(0,_yRes,0);
    VEC3F v011 = cellCenter(0,_yRes,_zRes);
    VEC3F v100 = cellCenter(_xRes,0,0);
    VEC3F v101 = cellCenter(_xRes,0,_zRes);
    VEC3F v110 = cellCenter(_xRes,_yRes,0);
    VEC3F v111 = cellCenter(_xRes, _yRes, _zRes);

    glColor4f(10,10,10,1);

    glBegin(GL_LINES);
#ifdef SINGLE_PRECISION
      glVertex3fv(v000);
      glVertex3fv(v001);

      glVertex3fv(v001);
      glVertex3fv(v011);
      
      glVertex3fv(v011);
      glVertex3fv(v010);

      glVertex3fv(v010);
      glVertex3fv(v000);

      glVertex3fv(v110);
      glVertex3fv(v111);
      
      glVertex3fv(v111);
      glVertex3fv(v101);

      glVertex3fv(v101);
      glVertex3fv(v100);
      
      glVertex3fv(v100);
      glVertex3fv(v110);

      glVertex3fv(v010);
      glVertex3fv(v110);

      glVertex3fv(v000);
      glVertex3fv(v100);
      
      glVertex3fv(v111);
      glVertex3fv(v011);
      
      glVertex3fv(v101);
      glVertex3fv(v001);
#else
      glVertex3dv(v000);
      glVertex3dv(v001);

      glVertex3dv(v001);
      glVertex3dv(v011);
      
      glVertex3dv(v011);
      glVertex3dv(v010);

      glVertex3dv(v010);
      glVertex3dv(v000);

      glVertex3dv(v110);
      glVertex3dv(v111);
      
      glVertex3dv(v111);
      glVertex3dv(v101);

      glVertex3dv(v101);
      glVertex3dv(v100);
      
      glVertex3dv(v100);
      glVertex3dv(v110);

      glVertex3dv(v010);
      glVertex3dv(v110);

      glVertex3dv(v000);
      glVertex3dv(v100);
      
      glVertex3dv(v111);
      glVertex3dv(v011);
      
      glVertex3dv(v101);
      glVertex3dv(v001);
#endif
    glEnd();
  glPopMatrix();
}

///////////////////////////////////////////////////////////////////////
// BLAS-like interface, output += alpha * input
///////////////////////////////////////////////////////////////////////
void FIELD_3D::axpy(const Real& alpha, const FIELD_3D& input)
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);

  int totalCells = input.totalCells();
  for (int x = 0; x < totalCells; x++)
    _data[x] += alpha * input[x];
}

///////////////////////////////////////////////////////////////////////
// set a border of size 1 to zero
///////////////////////////////////////////////////////////////////////
void FIELD_3D::setZeroBorder()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
		{
			// left slab
			index = y * _xRes + z * _slabSize;
			_data[index] = 0.0f;

			// right slab
			index += _xRes - 1;
			_data[index] = 0.0f;
		}
	for (int z = 0; z < _zRes; z++)
		for (int x = 0; x < _xRes; x++)
		{
			// bottom slab
			index = x + z * _slabSize;
			_data[index] = 0.0f;

			// top slab
			index += _slabSize - _xRes;
			_data[index] = 0.0f;
		}
	for (int y = 0; y < _yRes; y++)
		for (int x = 0; x < _xRes; x++)
		{
			// front slab
			index = x + y * _xRes;
			_data[index] = 0.0f;

			// back slab
			index += _totalCells - _slabSize;
			_data[index] = 0.0f;
		}
}

//////////////////////////////////////////////////////////////////////
// swap the contents with another object
//////////////////////////////////////////////////////////////////////
void FIELD_3D::swapPointers(FIELD_3D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);
  assert(field.zRes() == _zRes);
  
  Real* temp = _data;
  _data = field._data;
  field._data = temp;
}

//////////////////////////////////////////////////////////////////////
// take the dot product with respect to another field
//////////////////////////////////////////////////////////////////////
Real FIELD_3D::dot(const FIELD_3D& input) const
{
  assert(input.xRes() == _xRes);
  assert(input.yRes() == _yRes);
  assert(input.zRes() == _zRes);
  Real final = 0;

  for (int x = 0; x < _totalCells; x++)
    final += input[x] * _data[x];

  return final;
}

//////////////////////////////////////////////////////////////////////
// peel off the outer boundary of grid cells
//////////////////////////////////////////////////////////////////////
FIELD_3D FIELD_3D::peelBoundary() const
{
  FIELD_3D final(_xRes - 2, _yRes - 2, _zRes - 2);

  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
        final(x - 1, y - 1, z - 1) = (*this)(x,y,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a flattened array of all the field contents
//////////////////////////////////////////////////////////////////////
VECTOR FIELD_3D::flattened()
{
  VECTOR final(_totalCells);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        final[index] = (*this)(x,y,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a flattened array of all the field contents
//////////////////////////////////////////////////////////////////////
Eigen::VectorXd FIELD_3D::flattenedEigen()
{
  Eigen::VectorXd final(_totalCells);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
        final[index] = (*this)(x,y,z);

  return final;
}
