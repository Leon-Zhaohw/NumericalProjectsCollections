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
#include "VECTOR3_FIELD_3D.h"
#include <omp.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h> // OpenGL itself.
//#include <GL/glu.h> // GLU support library.
#include <GL/glut.h> // GLUT support library.
#endif


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const int& xRes, const int& yRes, const int& zRes,
    const VEC3F& center, const VEC3F& lengths) :
  _xRes(xRes), _yRes(yRes), _zRes(zRes), _center(center), _lengths(lengths), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D() :
  _xRes(0), _yRes(0), _zRes(0), _totalCells(0), _data(NULL), _initialized(false)
{
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D::~VECTOR3_FIELD_3D()
{
  delete[] _data;
}
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clear()
{
  TIMER functionTimer(__FUNCTION__);
  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0.0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::magnitudeField() const
{
  FIELD_3D final(_xRes, _yRes, _zRes, _center, _lengths);

  for (int x = 0; x < _totalCells; x++)
    final[x] = norm(_data[x]);

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the field dot product
///////////////////////////////////////////////////////////////////////
FIELD_3D operator*(const VECTOR3_FIELD_3D&u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  FIELD_3D final(u.xRes(), u.yRes(), u.zRes(), u.center(), u.lengths());

  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] * v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator*(const VECTOR3_FIELD_3D& v, const Real& a)
{
  return a * v;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator*(const Real& a, const VECTOR3_FIELD_3D& v)
{
  VECTOR3_FIELD_3D final(v);
  
  int totalCells = v.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = a * v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the field dot product
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator*(const FIELD_3D& u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  VECTOR3_FIELD_3D final(v);

  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] * v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// sum two vector fields
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator+(const VECTOR3_FIELD_3D&u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  VECTOR3_FIELD_3D final(u.xRes(), u.yRes(), u.zRes(), u.center(), u.lengths());
  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] + v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// diff two vector fields
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D operator-(const VECTOR3_FIELD_3D&u, const VECTOR3_FIELD_3D& v)
{
  assert(u.xRes() == v.xRes());
  assert(u.yRes() == v.yRes());
  assert(u.zRes() == v.zRes());

  VECTOR3_FIELD_3D final(u.xRes(), u.yRes(), u.zRes(), u.center(), u.lengths());
  int totalCells = u.totalCells();
  for (int x = 0; x < totalCells; x++)
    final[x] = u[x] - v[x];

  return final;
}

///////////////////////////////////////////////////////////////////////
// lookup value at some real-valued spatial position
///////////////////////////////////////////////////////////////////////
const VEC3F VECTOR3_FIELD_3D::operator()(const VEC3F& position) const
{
  VEC3F positionCopy = position - _center + (Real)0.5 * _lengths - (Real)0.5 * dxs();

  positionCopy[0] *= 1.0 / _dx;
  positionCopy[1] *= 1.0 / _dy;
  positionCopy[2] *= 1.0 / _dz;

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
  const Real s1 = positionCopy[0]- x0;
  const Real s0 = 1.0f - s1;
  const Real t1 = positionCopy[1]- y0;
  const Real t0 = 1.0f - t1;
  const Real u1 = positionCopy[2]- z0;
  const Real u0 = 1.0f - u1;

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
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator-=(const VECTOR3_FIELD_3D& input)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] -= input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator+=(const VECTOR3_FIELD_3D& input)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] += input[x];

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator*=(const FIELD_3D& input)
{
  assert(_xRes == input.xRes());
  assert(_yRes == input.yRes());
  assert(_zRes == input.zRes());

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] *= input[x];
    _data[x][1] *= input[x];
    _data[x][2] *= input[x];
  }

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator*=(const VEC3F& alpha)
{
  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0] *= alpha[0];
    _data[x][1] *= alpha[1];
    _data[x][2] *= alpha[2];
  }

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator+=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] += alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator*=(const Real& alpha)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= alpha;

  return *this;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator=(const Real& value)
{
  for (int x = 0; x < _totalCells; x++)
    _data[x] = value;

  return *this;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D& VECTOR3_FIELD_3D::operator=(const VECTOR3_FIELD_3D& input)
{
  if (input.xRes() != _xRes ||
      input.yRes() != _yRes ||
      input.zRes() != _zRes)
  {
    delete[] _data;

    _xRes = input.xRes();
    _yRes = input.yRes();
    _zRes = input.zRes();

    _totalCells = _xRes * _yRes * _zRes;
    _slabSize = _xRes * _yRes;
    _data = new VEC3F[_totalCells];
  }

  _center = input.center();
  _lengths = input.lengths();

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
  
  for (int x = 0; x < _totalCells; x++)
    _data[x] = input[x];

  _initialized = input._initialized;

  return *this;
}

///////////////////////////////////////////////////////////////////////
// write out a field to a file stream
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::writeGz(gzFile& file) const
{
  // write dimensions
  gzwrite(file, (void*)&_xRes, sizeof(int));
  gzwrite(file, (void*)&_yRes, sizeof(int));
  gzwrite(file, (void*)&_zRes, sizeof(int));
  _center.writeGz(file);
  _lengths.writeGz(file);

  double* dataDouble = new double[3 * _totalCells];
  for (int x = 0; x < _totalCells; x++)
  {
    dataDouble[3 * x] = _data[x][0];
    dataDouble[3 * x + 1] = _data[x][1];
    dataDouble[3 * x + 2] = _data[x][2];
  }

  gzwrite(file, (void*)dataDouble, sizeof(double) * 3 * _totalCells);
  delete[] dataDouble;
}

///////////////////////////////////////////////////////////////////////
// read in a field from a file stream
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::readGz(gzFile& file)
{
  // read dimensions
  gzread(file, (void*)&_xRes, sizeof(int));
  gzread(file, (void*)&_yRes, sizeof(int));
  gzread(file, (void*)&_zRes, sizeof(int));
  _center.readGz(file);
  _lengths.readGz(file);
  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  if (_data) delete[] _data;
  _data = new VEC3F[_totalCells];

  double* dataDouble = new double[3 * _totalCells];

  int totalBytes = sizeof(double) * 3 * _totalCells;
  int readBytes = gzread(file, (void*)dataDouble, totalBytes);

  assert(readBytes == totalBytes);

  for (int x = 0; x < _totalCells; x++)
  {
    _data[x][0]= dataDouble[3 * x];
    _data[x][1]= dataDouble[3 * x + 1];
    _data[x][2]= dataDouble[3 * x + 2];
  }

  delete[] dataDouble;
}

///////////////////////////////////////////////////////////////////////
// advect using first order semi-Lagrangian
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advect(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const VECTOR3_FIELD_3D& oldField, VECTOR3_FIELD_3D& newField)
{
  const int xRes = velocityGrid.xRes();
  const int yRes = velocityGrid.yRes();
  const int zRes = velocityGrid.zRes();
  const int slabSize = velocityGrid.slabSize();

  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        const VEC3F velocity = velocityGrid[index];
        Real xTrace = x - dt * velocity[0];
        Real yTrace = y - dt * velocity[1];
        Real zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
        xTrace = (xTrace > xRes - 2.5) ? xRes - 2.5 : xTrace;
        yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
        yTrace = (yTrace > yRes - 2.5) ? yRes - 2.5 : yTrace;
        zTrace = (zTrace < 1.5) ? 1.5 : zTrace;
        zTrace = (zTrace > zRes - 2.5) ? zRes - 2.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        // get interpolation weights
        const Real s1 = xTrace - x0;
        const Real s0 = 1.0f - s1;
        const Real t1 = yTrace - y0;
        const Real t0 = 1.0f - t1;
        const Real u1 = zTrace - z0;
        const Real u0 = 1.0f - u1;

        const int i000 = x0 + y0 * xRes + z0 * slabSize;
        const int i010 = x0 + y1 * xRes + z0 * slabSize;
        const int i100 = x1 + y0 * xRes + z0 * slabSize;
        const int i110 = x1 + y1 * xRes + z0 * slabSize;
        const int i001 = x0 + y0 * xRes + z1 * slabSize;
        const int i011 = x0 + y1 * xRes + z1 * slabSize;
        const int i101 = x1 + y0 * xRes + z1 * slabSize;
        const int i111 = x1 + y1 * xRes + z1 * slabSize;

        // interpolate
        // (indices could be computed once)
        newField[index] = u0 * (s0 * (t0 * oldField[i000] +
                                      t1 * oldField[i010]) +
                                s1 * (t0 * oldField[i100] +
                                      t1 * oldField[i110])) +
                          u1 * (s0 * (t0 * oldField[i001] +
                                      t1 * oldField[i011]) +
                                s1 * (t0 * oldField[i101] +
                                      t1 * oldField[i111]));
      }
}

///////////////////////////////////////////////////////////////////////
// advect using first order semi-Lagrangian
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advect(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField)
{
  const int xRes = velocityGrid.xRes();
  const int yRes = velocityGrid.yRes();
  const int zRes = velocityGrid.zRes();
  const int slabSize = velocityGrid.slabSize();

  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        const VEC3F velocity = velocityGrid[index];
        Real xTrace = x - dt * velocity[0];
        Real yTrace = y - dt * velocity[1];
        Real zTrace = z - dt * velocity[2];

        // clamp backtrace to grid boundaries
        xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
        xTrace = (xTrace > xRes - 2.5) ? xRes - 2.5 : xTrace;
        yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
        yTrace = (yTrace > yRes - 2.5) ? yRes - 2.5 : yTrace;
        zTrace = (zTrace < 1.5) ? 1.5 : zTrace;
        zTrace = (zTrace > zRes - 2.5) ? zRes - 2.5 : zTrace;

        // locate neighbors to interpolate
        const int x0 = (int)xTrace;
        const int x1 = x0 + 1;
        const int y0 = (int)yTrace;
        const int y1 = y0 + 1;
        const int z0 = (int)zTrace;
        const int z1 = z0 + 1;

        // get interpolation weights
        const Real s1 = xTrace - x0;
        const Real s0 = 1.0f - s1;
        const Real t1 = yTrace - y0;
        const Real t0 = 1.0f - t1;
        const Real u1 = zTrace - z0;
        const Real u0 = 1.0f - u1;

        const int i000 = x0 + y0 * xRes + z0 * slabSize;
        const int i010 = x0 + y1 * xRes + z0 * slabSize;
        const int i100 = x1 + y0 * xRes + z0 * slabSize;
        const int i110 = x1 + y1 * xRes + z0 * slabSize;
        const int i001 = x0 + y0 * xRes + z1 * slabSize;
        const int i011 = x0 + y1 * xRes + z1 * slabSize;
        const int i101 = x1 + y0 * xRes + z1 * slabSize;
        const int i111 = x1 + y1 * xRes + z1 * slabSize;

        // interpolate
        // (indices could be computed once)
        newField[index] = u0 * (s0 * (t0 * oldField[i000] +
                                      t1 * oldField[i010]) +
                                s1 * (t0 * oldField[i100] +
                                      t1 * oldField[i110])) +
                          u1 * (s0 * (t0 * oldField[i001] +
                                      t1 * oldField[i011]) +
                                s1 * (t0 * oldField[i101] +
                                      t1 * oldField[i111]));
      }
}

//////////////////////////////////////////////////////////////////////
// set everything on the border to zero
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setZeroBorder()
{
  TIMER functionTimer(__FUNCTION__);
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
// set x direction to zero
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setZeroX()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
		{
			// left slab
			index = y * _xRes + z * _slabSize;
			_data[index][0] = 0.0f;

			// right slab
			index += _xRes - 1;
			_data[index][0] = 0.0f;
		}
}

//////////////////////////////////////////////////////////////////////
// set y direction to zero
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setZeroY()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int x = 0; x < _xRes; x++)
		{
			// bottom slab
			index = x + z * _slabSize;
			_data[index][1] = 0.0f;

			// top slab
			index += _slabSize - _xRes;
			_data[index][1] = 0.0f;
		}
}

//////////////////////////////////////////////////////////////////////
// set z direction to zero
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setZeroZ()
{
	int index;
	for (int y = 0; y < _yRes; y++)
		for (int x = 0; x < _xRes; x++)
		{
			// front slab
			index = x + y * _xRes;
			_data[index][2] = 0.0f;

			// back slab
			index += _totalCells - _slabSize;
			_data[index][2] = 0.0f;
		}
}

//////////////////////////////////////////////////////////////////////
// set x direction to Neumann boundary conditions
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setNeumannX()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
		{
			// left slab
			index = y * _xRes + z * _slabSize;
			_data[index][0] = _data[index + 2][0];

			// right slab
			index += _xRes - 1;
			_data[index][0] = _data[index - 2][0];
		}
}


//////////////////////////////////////////////////////////////////////
// set y direction to Neumann boundary conditions
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setNeumannY()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int x = 0; x < _xRes; x++)
		{
			// bottom slab
			index = x + z * _slabSize;
			_data[index][1] = _data[index + 2 * _xRes][1];

			// top slab
			index += _slabSize - _xRes;
			_data[index][1] = _data[index - 2 * _xRes][1];
		}
}

//////////////////////////////////////////////////////////////////////
// set z direction to Neumann boundary conditions
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setNeumannZ()
{
	int index;
	for (int y = 0; y < _yRes; y++)
		for (int x = 0; x < _xRes; x++)
		{
			// front slab
			index = x + y * _xRes;
			_data[index][2] = _data[index + 2 * _slabSize][2];

			// back slab
			index += _totalCells - _slabSize;
			_data[index][2] = _data[index - 2 * _slabSize][2];
		}
}

//////////////////////////////////////////////////////////////////////
// copy grid boundary
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::copyBorderX()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
		{
			// left slab
			index = y * _xRes + z * _slabSize;
			_data[index][0] = _data[index + 1][0];

			// right slab
			index += _xRes - 1;
			_data[index][0] = _data[index - 1][0];
		}
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::copyBorderY()
{
	int index;
	for (int z = 0; z < _zRes; z++)
		for (int x = 0; x < _xRes; x++)
		{
			// bottom slab
			index = x + z * _slabSize;
			_data[index][1] = _data[index + _xRes][1]; 
			// top slab
			index += _slabSize - _xRes;
			_data[index][1] = _data[index - _xRes][1];
		}
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::copyBorderZ()
{
	int index;
	for (int y = 0; y < _yRes; y++)
		for (int x = 0; x < _xRes; x++)
		{
			// front slab
			index = x + y * _xRes;
			_data[index][2] = _data[index + _slabSize][2]; 
			// back slab
			index += _totalCells - _slabSize;
			_data[index][2] = _data[index - _slabSize][2];
		}
}

//////////////////////////////////////////////////////////////////////
// BLAS-like interface, output += alpha * input
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::axpy(const Real& alpha, const VECTOR3_FIELD_3D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);
  assert(field.zRes() == _zRes);
  for (int x = 0; x < _totalCells; x++)
    _data[x] += alpha * field[x];
}

//////////////////////////////////////////////////////////////////////
// swap the contents with another object
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::swapPointers(VECTOR3_FIELD_3D& field)
{
  assert(field.xRes() == _xRes);
  assert(field.yRes() == _yRes);
  assert(field.zRes() == _zRes);
  
  VEC3F* temp = _data;
  _data = field._data;
  field._data = temp;
}

//////////////////////////////////////////////////////////////////////
// return a flattened array of all the field contents
//////////////////////////////////////////////////////////////////////
VECTOR VECTOR3_FIELD_3D::flattened()
{
  VECTOR final(_totalCells * 3);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        for (int i = 0; i < 3; i++, index++)
          final[index] = (*this)(x,y,z)[i];

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a flattened array of all the field contents
//////////////////////////////////////////////////////////////////////
VectorXd VECTOR3_FIELD_3D::flattenedEigen()
{
  VectorXd final(_totalCells * 3);

  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        for (int i = 0; i < 3; i++, index++)
          final(index) = (*this)(x,y,z)[i];

  return final;
}

//////////////////////////////////////////////////////////////////////
// peel off the outer boundary of grid cells in preparation for PCA
//////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::peelBoundary() const
{
  VECTOR3_FIELD_3D final(_xRes - 2, _yRes - 2, _zRes - 2);

  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
        final(x - 1, y - 1, z - 1) = (*this)(x,y,z);

  return final;
}

//////////////////////////////////////////////////////////////////////
// do a projection of the peeled field, using the passed in basis
//////////////////////////////////////////////////////////////////////
VectorXd VECTOR3_FIELD_3D::peeledProject(const MatrixXd& U)
{
  TIMER functionTimer(__FUNCTION__);

  const int xPeeled = _xRes - 2;
  const int slabPeeled = (_xRes - 2) * (_yRes - 2);
  const int totalThreads = omp_get_max_threads();

  vector<VectorXd> finals(totalThreads);
  for (unsigned int x = 0; x < finals.size(); x++)
  {
    finals[x].resize(U.cols());
    finals[x].setZero();
  }

  const int totalColumns = U.cols();
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes - 2; z++)
  {
    const int zCached = 3 * z * slabPeeled;
    const int threadID = omp_get_thread_num();
    VectorXd& currentFinal = finals[threadID];

    int stride = 4;
    VectorXd unroll0(finals[threadID].size());
    VectorXd unroll1(finals[threadID].size());
    VectorXd unroll2(finals[threadID].size());
    VectorXd unroll3(finals[threadID].size());

    for (int y = 0; y < _yRes - 2; y++)
    {
      const int yzCached = 3 * y * xPeeled + zCached;

      int xEnd = (_xRes - 2) / stride * stride;

      unroll0.setZero();
      unroll1.setZero();
      unroll2.setZero();
      unroll3.setZero();
      for (int x = 0; x < xEnd; x += stride)
      {
        const int index0 = 3 * x + yzCached;
        const int index1 = 3 * (x + 1) + yzCached;
        const int index2 = 3 * (x + 2) + yzCached;
        const int index3 = 3 * (x + 3) + yzCached;

        const VEC3F v30 = (*this)(x + 1, y + 1, z + 1);
        const VEC3F v31 = (*this)(x + 2, y + 1, z + 1);
        const VEC3F v32 = (*this)(x + 3, y + 1, z + 1);
        const VEC3F v33 = (*this)(x + 4, y + 1, z + 1);

        const Vector3d v0(v30[0], v30[1], v30[2]);
        const Vector3d v1(v31[0], v31[1], v31[2]);
        const Vector3d v2(v32[0], v32[1], v32[2]);
        const Vector3d v3(v33[0], v33[1], v33[2]);

        unroll0.noalias() += U.block(index0, 0, 3, totalColumns).transpose() * v0;
        unroll1.noalias() += U.block(index1, 0, 3, totalColumns).transpose() * v1;
        unroll2.noalias() += U.block(index2, 0, 3, totalColumns).transpose() * v2;
        unroll3.noalias() += U.block(index3, 0, 3, totalColumns).transpose() * v3;
      }
      currentFinal.noalias() += unroll0;
      currentFinal.noalias() += unroll1;
      currentFinal.noalias() += unroll2;
      currentFinal.noalias() += unroll3;

      // unrolling leftovers
      for (int x = xEnd; x < _xRes - 2; x++)
      {
        const int index = 3 * x + yzCached;
        const VEC3F v3 = (*this)(x + 1, y + 1, z + 1);
        const Vector3d v(v3[0], v3[1], v3[2]);
        currentFinal.noalias() += U.block(index, 0, 3, totalColumns).transpose() * v;
      }
    }
  }

  TIMER finalReduce("peeledProject final reduce");

  // do a sum of the thread results
  VectorXd final(U.cols());
  final.setZero();
  for (unsigned int x = 0; x < finals.size(); x++)
    final += finals[x];
  finalReduce.stop();

  return final;
}

//////////////////////////////////////////////////////////////////////
// unproject the reduced coordinate into the peeled cells in this field
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::peeledUnproject(const MatrixXd& U, const VectorXd& q)
{
  TIMER functionTimer(__FUNCTION__);
  const int xPeeled = _xRes - 2;
  const int slabPeeled = (_xRes - 2) * (_yRes - 2);
  const int totalColumns = U.cols();

#pragma omp parallel
#pragma omp for  schedule(static)
  for (int z = 0; z < _zRes - 2; z++)
  {
    const int zCached = 3 * z * slabPeeled;
    for (int y = 0; y < _yRes - 2; y++)
    {
      const int yzCached = 3 * y * xPeeled + zCached;

      int stride = 2;
      int xEnd = (_xRes - 2) / stride * stride;

      // do the unrolling
      for (int x = 0; x < xEnd; x += stride)
      {
        const int index0 = 3 * x + yzCached;
        const int index1 = 3 * (x + 1) + yzCached;

        Real finalX0 = U.row(index0).dot(q);
        Real finalY0 = U.row(index0 + 1).dot(q);
        Real finalZ0 = U.row(index0 + 2).dot(q);
        Real finalX1 = U.row(index1).dot(q);
        Real finalY1 = U.row(index1 + 1).dot(q);
        Real finalZ1 = U.row(index1 + 2).dot(q);

        (*this)(x + 1, y + 1, z + 1) = VEC3F(finalX0, finalY0, finalZ0);
        (*this)(x + 2, y + 1, z + 1) = VEC3F(finalX1, finalY1, finalZ1);
      }

      // catch the leftovers
      for (int x = xEnd; x < _xRes - 2; x++)
      {
        const int index = 3 * x + yzCached;
        Real finalX = U.row(index).dot(q);
        Real finalY = U.row(index + 1).dot(q);
        Real finalZ = U.row(index + 2).dot(q);
        (*this)(x + 1, y + 1, z + 1) = VEC3F(finalX, finalY, finalZ);
      }

    }
  }
}

//////////////////////////////////////////////////////////////////////
// take the dot product of the current field with another vector field
// and return the scalar field
//////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::dot(const VECTOR3_FIELD_3D& rhs)
{
  FIELD_3D final(_xRes, _yRes, _zRes);

  assert(rhs.xRes() == _xRes);
  assert(rhs.yRes() == _yRes);
  assert(rhs.zRes() == _zRes);

  //for (int x = 0; x < _totalCells; x++)
  //  final[x] = rhs[x] * (*this)[x];
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        final(x,y,z) = rhs(x,y,z) * (*this)(x,y,z);
        if (x == 0 && y == 0 && z == 0)
        {
          cout << "RHS: " << rhs(x,y,z) << endl;
          cout << "this: " << (*this)(x,y,z) << endl;
          cout << "final: " << final(x,y,z) << endl;
        }

      }

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::normalizeToLargest()
{
  Real maxMagnitude = 0.0;
  for (int x = 0; x < _totalCells; x++)
  {
    Real candidate = norm(_data[x]);
    if (candidate > maxMagnitude)
      maxMagnitude = candidate;
  }

  Real inverse = 1.0 / maxMagnitude;
  for (int x = 0; x < _totalCells; x++)
    _data[x] *= inverse;
}

//////////////////////////////////////////////////////////////////////
// absolute max single entry
//////////////////////////////////////////////////////////////////////
Real VECTOR3_FIELD_3D::maxAbsScalar()
{
  Real final = 0;
  for (int x = 0; x < _totalCells; x++)
    for (int y = 0; y < 3; y++)
      final = (fabs(_data[x][y]) > final) ? fabs(_data[x][y]) : final;

  return final;
}

//////////////////////////////////////////////////////////////////////
// 2 norm of the whole field
//////////////////////////////////////////////////////////////////////
Real VECTOR3_FIELD_3D::twoNorm()
{
  Real final = 0;
  for (int x = 0; x < _totalCells; x++)
    for (int y = 0; y < 3; y++)
      final += _data[x][y] * _data[x][y];

  return sqrt(final);
}

//////////////////////////////////////////////////////////////////////
// set the field innards to a peeled version
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setWithPeeled(const VectorXd& data)
{
  TIMER functionTimer(__FUNCTION__);
  assert(data.size() == 3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2));

  int index = 0;
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++, index++)
      {
        (*this)(x,y,z)[0] = data[3 * index];
        (*this)(x,y,z)[1] = data[3 * index + 1];
        (*this)(x,y,z)[2] = data[3 * index + 2];
      }
}
