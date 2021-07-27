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

#include "VECTOR3_FIELD_3D.h"
#include <omp.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h> // OpenGL itself.
#include <GL/glut.h> // GLUT support library.
#endif

#include <set>

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

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D(const FIELD_3D& m) :
  _xRes(m.xRes()), _yRes(m.yRes()), _zRes(m.zRes()),
  _center(m.center()), _lengths(m.lengths()), _initialized(true)
{
  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
}

VECTOR3_FIELD_3D::VECTOR3_FIELD_3D() :
  _xRes(-1), _yRes(-1), _zRes(-1), _totalCells(-1), _data(NULL), _initialized(false)
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
  for (int x = 0; x < _totalCells; x++)
    _data[x] = 0.0;
}

///////////////////////////////////////////////////////////////////////
// reset the lengths to something else, and recompute all the 
// dimesions as well
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::setLengths(const VEC3F& lengths)
{
  _lengths = lengths;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;
}

///////////////////////////////////////////////////////////////////////
// create a field of the grid positions of the passed in grid
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::cellCenters(const FIELD_3D& input)
{
  int xRes = input.xRes();
  int yRes = input.yRes();
  int zRes = input.zRes();
  const VEC3F& center = input.center();
  const VEC3F& lengths = input.lengths();
  VECTOR3_FIELD_3D final(xRes, yRes, zRes, center, lengths);

  int index = 0;
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++, index++)
        final[index] = input.cellCenter(x,y,z);

  return final;
}

///////////////////////////////////////////////////////////////////////
// take the gradient of a scalar field
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::gradient(const FIELD_3D& input)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = input.xRes();
  const int yRes = input.yRes();
  const int zRes = input.zRes();
  const int slabSize = xRes * yRes;
  const int totalCells = xRes * yRes * zRes;
  const VEC3F& center = input.center();
  const VEC3F& lengths = input.lengths();
  VECTOR3_FIELD_3D final(xRes, yRes, zRes, center, lengths);
 
  const Real dxHalfInv = 0.5 / input.dx();
  const Real dyHalfInv = 0.5 / input.dy();
  const Real dzHalfInv = 0.5 / input.dz();

  // do the x middles
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 1; x < xRes - 1; x++)
      {
        int index = x + y * xRes + z * slabSize;
        final[index][0] = (input[index + 1] - input[index - 1]) * dxHalfInv;
      }

  // do the y middles
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 1; y < yRes - 1; y++)
      for (int x = 0; x < xRes; x++)
      {
        int index = x + y * xRes + z * slabSize;
        final[index][1] = (input[index + xRes] - input[index - xRes]) * dyHalfInv;
      }

  // do the z middles
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 1; z < zRes - 1; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        int index = x + y * xRes + z * slabSize;
        final[index][2] = (input[index + slabSize] - input[index - slabSize]) * dzHalfInv;
      }

  // reset dx's to a single cell
  Real dxInv = 1.0 / input.dx();
  Real dyInv = 1.0 / input.dy();
  Real dzInv = 1.0 / input.dz();

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++)
    {
      // front slab
      int index = x + y * xRes;
      final[index][2] = (input[index + slabSize] - input[index]) * dzInv;

      // back slab
      index += totalCells - slabSize;
      final[index][2] = (input[index] - input[index - slabSize]) * dzInv;
    }

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int x = 0; x < xRes; x++)
    {
      // bottom slab
      int index = x + z * slabSize;
      final[index][1] = (input[index + xRes] - input[index]) * dyInv;

      // top slab
      index += slabSize - xRes;
      final[index][1] = (input[index] - input[index - xRes]) * dyInv;
    }

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
    {
      // left slab
      int index = y * xRes + z * slabSize;
      final[index][0] = (input[index + 1] - input[index]) * dxInv;

      // right slab
      index += xRes - 1;
      final[index][0] = (input[index] - input[index - 1]) * dxInv;
    }

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
  
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < _totalCells; x++)
    _data[x] = input[x];

  _initialized = input._initialized;

  return *this;
}

//////////////////////////////////////////////////////////////////////
// set the values in the field to the values at the closest points,
// but exclude a central core of cells, so it's a hose, not a band
//
// Might be better to called it a "cored band"
//////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBandFrozenCore(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints, const FIELD_3D& distance, int maxRadius, int coreRadius)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Setting to closest points, narrow band, frozen core..."; flush(cout);
  FIELD_3D final(input);
  const Real invDx = 1.0 / input.dx();

  const int xRes = closestPoints.xRes();
  const int yRes = closestPoints.yRes();
  const int zRes = closestPoints.zRes();

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        Real currentDistance = fabs(distance(x,y,z) * invDx);
        if (currentDistance >= maxRadius || currentDistance <= coreRadius) continue;
        
        final(x,y,z) = input.quarticLookup(closestPoints(x,y,z));
      }
  cout << " done." << endl;

  return final;
}

//////////////////////////////////////////////////////////////////////
// set the values in the field to the values at the closest points
//////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::setToClosestPointValues(const FIELD_3D& input, const VECTOR3_FIELD_3D& closestPoints)
{
  TIMER functionTimer(__FUNCTION__);
  FIELD_3D final(input);

  cout << " Using quartic ... "; flush(cout);

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < closestPoints.zRes(); z++)
    for (int y = 0; y < closestPoints.yRes(); y++)
      for (int x = 0; x < closestPoints.xRes(); x++)
        final(x,y,z) = input.quarticLookup(closestPoints(x,y,z));

  return final;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::readPhysBAMGz(const FIELD_3D& distanceField, const string filename)
{
  VECTOR3_FIELD_3D unpadded;
  unpadded.readPhysBAMGz(filename);

  *this = VECTOR3_FIELD_3D(distanceField);

  // appears to be, at least for these boundary conditions
  for (int z = 0; z < unpadded.zRes(); z++)
    for (int y = 0; y < unpadded.yRes(); y++)
      for (int x = 0; x < unpadded.xRes(); x++)
        (*this)(x + 2, y + 2, z + 2) = unpadded(x,y,z);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::readPhysBAMGz(const string filename)
{
  gzFile file;
  file = gzopen(filename.c_str(), "rb1");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " VECTOR3_FIELD_3D read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }
  cout << " Reading PhysBAM velocity " << filename << " ... "; flush(cout);

  // read in a RANGE<TV_INT>, i.e. index bounds
  int uBegin;
  int uEnd;
  int vBegin;
  int vEnd;
  int wBegin;
  int wEnd;
    
  gzread(file, (void*)&uBegin, sizeof(int));
  gzread(file, (void*)&uEnd, sizeof(int));
  gzread(file, (void*)&vBegin, sizeof(int));
  gzread(file, (void*)&vEnd, sizeof(int));
  gzread(file, (void*)&wBegin, sizeof(int));
  gzread(file, (void*)&wEnd, sizeof(int));

  int temp = uEnd;
  uEnd = wEnd;
  wEnd = temp;

  int scanTotalX = (uEnd + 1) * vEnd * wEnd;
  int scanTotalY = uEnd * (vEnd + 1) * wEnd;
  int scanTotalZ = uEnd * vEnd * (wEnd + 1);
  int totalCells;
  gzread(file, (void*)&totalCells, sizeof(int));

  // am I computing the total cells right?
  assert(guess == totalCells);

  float* data = new float[totalCells];

  gzread(file, data, sizeof(float) * totalCells);
 
  FIELD_3D zCompMAC(uEnd, vEnd, wEnd + 1);
  FIELD_3D yCompMAC(uEnd, vEnd + 1, wEnd);
  FIELD_3D xCompMAC(uEnd + 1, vEnd, wEnd);

  cout << " uRes: " << uEnd << " vEnd: " << vEnd << " wEnd: " << wEnd << endl;

  // Z first!?
  //
  // TK: this is correct -- PhysBAM does a zyx (reverse) ordering for some reason

  for (int x = 0; x < scanTotalZ; x++)
    zCompMAC[x] = data[x];

  for (int x = 0; x < scanTotalY; x++)
    yCompMAC[x] = data[scanTotalZ + x];
  
  for (int x = 0; x < scanTotalX; x++)
    xCompMAC[x] = data[scanTotalZ + scanTotalY + x];

  FIELD_3D xComp(uEnd, vEnd, wEnd);
  FIELD_3D yComp(uEnd, vEnd, wEnd);
  FIELD_3D zComp(uEnd, vEnd, wEnd);

  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++)
        xComp(x,y,z) = (xCompMAC(x,y,z) + xCompMAC(x+1,y,z)) * 0.5;
  
  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++)
        yComp(x,y,z) = (yCompMAC(x,y,z) + yCompMAC(x,y+1,z)) * 0.5;
  
  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++)
        zComp(x,y,z) = (zCompMAC(x,y,z) + zCompMAC(x,y,z+1)) * 0.5;

  _xRes = uEnd;
  _yRes = vEnd;
  _zRes = wEnd;

  // set longest dimension to 1
  int biggest = (_xRes > _yRes) ? _xRes : _yRes;
  biggest = (_zRes > biggest) ? _zRes : biggest;
  _lengths[0] = (Real)_xRes / biggest;
  _lengths[1] = (Real)_yRes / biggest;
  _lengths[2] = (Real)_zRes / biggest;

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;

  if (_data) delete[] _data;
  _data = new VEC3F[_totalCells];

  int index = 0;
  for (int z = 0; z < wEnd; z++)
    for (int y = 0; y < vEnd; y++)
      for (int x = 0; x < uEnd; x++, index++)
      {
        _data[index][0] = xComp(x,y,z);
        _data[index][1] = yComp(x,y,z);
        _data[index][2] = zComp(x,y,z);
      }

  delete[] data;
  gzclose(file);

  cout << " done. " << endl;
}

///////////////////////////////////////////////////////////////////////
// advect using first order semi-Lagrangian
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::advectNarrowBand(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& oldField, FIELD_3D& newField, const FIELD_3D& distance, const int maxCells)
{
  TIMER functionTimer(__FUNCTION__);

  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();
  const Real invDx = 1.0 / oldField.dx();
 
  const VEC3F corner = oldField.center() - (Real)0.5 * oldField.lengths() + (Real)0.5 * oldField.dxs();
  const VEC3F dxs = oldField.dxs();

  // scale dt up to grid resolution
#pragma omp parallel
#pragma omp for schedule(static)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        Real currentDistance = fabs(distance(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;
        const int index = x + y * xRes + z * slabSize;
        
        // backtrace
        const VEC3F velocity = velocityGrid(oldField.cellCenter(x,y,z));
        VEC3F position(x - dt * velocity[0], y - dt * velocity[1], z - dt * velocity[2]);
        position[0] *= dxs[0];
        position[1] *= dxs[1];
        position[2] *= dxs[2];
        position += corner;
        newField[index] = oldField.quarticLookup(position);
      }

}

///////////////////////////////////////////////////////////////////////
// compute the extension field for a masked set of points
///////////////////////////////////////////////////////////////////////
FIELD_3D VECTOR3_FIELD_3D::computeExtensionFieldMasked(const FIELD_3D& distanceField, const FIELD_3D& toExtend, const FIELD_3D& mask)
{
  TIMER functionTimer(__FUNCTION__);

  // get the cell centers 
  VECTOR3_FIELD_3D cellCenters = VECTOR3_FIELD_3D::cellCenters(toExtend);

  // build the gradient field
  VECTOR3_FIELD_3D targetGradient = VECTOR3_FIELD_3D::gradient(distanceField);

  const int xRes = mask.xRes();
  const int yRes = mask.yRes();
  const int zRes = mask.zRes();

  const Real invDx = 1.0 / distanceField.dx();

  // perform Nacelle on all the tagged extension cells
  FIELD_3D extendedOld(toExtend);
  int maxSteps = 100;
  Real fictionalDt = 1.0;
  cout << " Computing masked extension field ..."; flush(cout);
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        VEC3F& position = cellCenters(x,y,z);
        Real currentDistance = fabs(distanceField(position) * invDx);

        //if (mask(x,y,z) < 0.5) continue;
        if (mask(x,y,z) < 0.5 || currentDistance <= 1.0) continue;

        Real diff = 1;
        int steps = 0;
        VEC3F targetDelta = targetGradient(position);

        while (diff > 1e-6 && steps < maxSteps)
        {
          const Real targetDistance = distanceField(position);
          diff = fabs(targetDistance);

          // if the gradient is zero, stop trying, because it probably
          // wandered off the grid
          if (norm2(targetDelta) < 1e-6) break;

          // go ahead and always do first -- second gives the occasional wacky value
          // that throws off the stability of the simulation
          targetDelta.normalize();
          const VEC3F move = (targetDistance) * targetDelta * fictionalDt;

          position = position - move;

          // update the direction
          targetDelta = targetGradient(position);
          steps++;
        }

        // do the extension here
        //extendedOld(x,y,z) = toExtend.quarticLookup(position);
        extendedOld(x,y,z) = toExtend.nearestNeighborLookup(position);
      }

  return extendedOld;
}

///////////////////////////////////////////////////////////////////////
// Clamp the extrema generated by the BFECC error correction
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clampExtrema(const Real dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, FIELD_3D& newField)
{
  TIMER functionTimer(__FUNCTION__);
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();

#pragma omp parallel
#pragma omp for schedule(dynamic)
	for (int z = 1; z < zRes-1; z++)
		for (int y = 1; y < yRes-1; y++)
			for (int x = 1; x < xRes-1; x++)
			{
				// backtrace
        const VEC3F velocity = velocityField(oldField.cellCenter(x,y,z));
        if (norm2(velocity) < 1e-6) continue;

        Real xTrace = x - dt * velocity[0];
        Real yTrace = y - dt * velocity[1];
        Real zTrace = z - dt * velocity[2];

				// clamp backtrace to grid boundaries
        xTrace = (xTrace < 0.5) ? 0.5 : xTrace;
        xTrace = (xTrace > xRes - 1.5) ? xRes - 1.5 : xTrace;
        yTrace = (yTrace < 0.5) ? 0.5 : yTrace;
        yTrace = (yTrace > yRes - 1.5) ? yRes - 1.5 : yTrace;
        zTrace = (zTrace < 0.5) ? 0.5 : zTrace;
        zTrace = (zTrace > zRes - 1.5) ? zRes - 1.5 : zTrace;

				// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;
				const int z0 = (int)zTrace;
				const int z1 = z0 + 1;

        Real minField = oldField(x0,y0,z0);
        Real maxField = minField;

        const Real x0y1z0 = oldField(x0,y1,z0);
        const Real x1y0z0 = oldField(x1,y0,z0);
        const Real x1y1z0 = oldField(x1,y1,z0);
        const Real x0y0z1 = oldField(x0,y0,z1);
        const Real x0y1z1 = oldField(x0,y1,z1);
        const Real x1y0z1 = oldField(x1,y0,z1);
        const Real x1y1z1 = oldField(x1,y1,z1);

        minField = (x0y1z0 < minField) ? x0y1z0 : minField;
        maxField = (x0y1z0 > maxField) ? x0y1z0 : maxField;

        minField = (x1y0z0 < minField) ? x1y0z0 : minField;
        maxField = (x1y0z0 > maxField) ? x1y0z0 : maxField;

        minField = (x1y1z0 < minField) ? x1y1z0 : minField;
        maxField = (x1y1z0 > maxField) ? x1y1z0 : maxField;

        minField = (x0y0z1 < minField) ? x0y0z1 : minField;
        maxField = (x0y0z1 > maxField) ? x0y0z1 : maxField;

        minField = (x0y1z1 < minField) ? x0y1z1 : minField;
        maxField = (x0y1z1 > maxField) ? x0y1z1 : maxField;

        minField = (x1y0z1 < minField) ? x1y0z1 : minField;
        maxField = (x1y0z1 > maxField) ? x1y0z1 : maxField;

        minField = (x1y1z1 < minField) ? x1y1z1 : minField;
        maxField = (x1y1z1 > maxField) ? x1y1z1 : maxField;

        const Real newValue = newField(x,y,z);

        newField(x,y,z) = (newValue > maxField) ? maxField : newValue;
        newField(x,y,z) = (newValue < minField) ? minField : newValue;
      }
}

//////////////////////////////////////////////////////////////////////
// Reverts any backtraces that go into boundaries back to first 
// order -- in this case the error correction term was totally
// incorrect
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::clampOutsideRays(const Real dt, const VECTOR3_FIELD_3D& velocityField, const FIELD_3D& oldField, const FIELD_3D& oldAdvection, FIELD_3D& newField)
{
  TIMER functionTimer(__FUNCTION__);
  // this actually is correct -- should be using the res of the
  // low res grid
  const int sx = velocityField.xRes();
  const int sy = velocityField.yRes();
  const int sz = velocityField.zRes();

#pragma omp parallel
#pragma omp for schedule(dynamic)
	for (int z = 0; z < oldField.zRes(); z++)
		for (int y = 0; y < oldField.yRes(); y++)
			for (int x = 0; x < oldField.xRes(); x++)
			{
				//const int index = x + y * sx + z * slabSize;
				// backtrace
        VEC3F velocity = velocityField(oldField.cellCenter(x,y,z));

        // this line seems to *significantly* reduce banding
        if (norm2(velocity) < 1e-6) continue;

        velocity *= dt;

				float xBackward = x + velocity[0];
				float yBackward = y + velocity[1];
				float zBackward = z + velocity[2];

				float xTrace    = x - velocity[0];
				float yTrace    = x - velocity[1];
				float zTrace    = x - velocity[2];

				// see if it goes outside the boundaries
				bool hasObstacle = 
					(zTrace < 1.0f)    || (zTrace > sz - 2.0f) ||
					(yTrace < 1.0f)    || (yTrace > sy - 2.0f) ||
					(xTrace < 1.0f)    || (xTrace > sx - 2.0f) ||
					(zBackward < 1.0f) || (zBackward > sz - 2.0f) ||
					(yBackward < 1.0f) || (yBackward > sy - 2.0f) ||
					(xBackward < 1.0f) || (xBackward > sx - 2.0f);

				// reuse old advection instead of doing another one...
				if(hasObstacle) { newField(x,y,z) = oldAdvection(x,y,z); continue; }

			} // xyz
}

///////////////////////////////////////////////////////////////////////
// copy values out into the border, assuming that "borderSize" is the 
// width of the grid padding
///////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::copyIntoBorder(int borderSize)
{
  TIMER functionTimer(__FUNCTION__);
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
      {
        if (x == borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x - i,y,z) = value;
        }
        if (x == _xRes - 1 - borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x + i,y,z) = value;
        }              
        if (y == borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y - i,z) = value;
        }
        if (y == _yRes - 1 - borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y+i,z) = value;
        }              
        if (z == borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z - i) = value;
        }
        if (z == _zRes - 1 - borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            (*this)(x,y,z+i) = value;
        }

        // handle the corners
        if (x == borderSize && z == borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x - i,y,z-j) = value;
        }
        if (x == _xRes - 1 - borderSize && z == _zRes - 1 - borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x + i,y,z+j) = value;
        }

        if (z == borderSize && y == borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y - j,z -i) = value;
        }
        if (z == _xRes - 1 - borderSize && y == _yRes - 1 - borderSize)
        {
          VEC3F value = (*this)(x,y,z);
          for (int i = 1; i <= borderSize; i++)
            for (int j = 1; j <= borderSize; j++)
              (*this)(x,y + j,z+i) = value;
        }

      }
}

///////////////////////////////////////////////////////////////////////
// pass back a field with a new padding of size "paddingSize"
///////////////////////////////////////////////////////////////////////
VECTOR3_FIELD_3D VECTOR3_FIELD_3D::withAddedPadding(int paddingSize) const
{
  // new length, with padding
  VEC3F newLengths = _lengths;
  newLengths[0] += paddingSize * 2 * _dx;
  newLengths[1] += paddingSize * 2 * _dx;
  newLengths[2] += paddingSize * 2 * _dx;

  VECTOR3_FIELD_3D final(_xRes + 2 * paddingSize, 
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

//////////////////////////////////////////////////////////////////////
// compute closest point field
//
// For each cell center in input, compute the position of the closest point
// in the surface field
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::computeClosestPointsNarrowBand(const FIELD_3D& input, const FIELD_3D& surfaceField, const int maxCells, VECTOR3_FIELD_3D& final)
{
  TIMER functionTimer(__FUNCTION__);
  const VECTOR3_FIELD_3D targetGradient = VECTOR3_FIELD_3D::gradient(surfaceField);
  
  const Real invDx = 1.0 / surfaceField.dx();

  const int maxSteps = 100;
  const int xRes = input.xRes();
  const int yRes = input.yRes();
  const int zRes = input.zRes();

  const VEC3F lowerCorner = input.center() - (Real)0.5 * input.lengths(); 
  const VEC3F upperCorner = input.center() + (Real)0.5 * input.lengths(); 

#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        Real currentDistance = fabs(surfaceField(x,y,z) * invDx);
        if (currentDistance >= maxCells) continue;

        Real diff = 1;
        int steps = 0;
        VEC3F& position = final(x,y,z);
        position = input.cellCenter(x,y,z);
        VEC3F targetDelta = targetGradient(position);

        // for simplicity we're just doing first order Nacelle here
        while (diff > 1e-6 && steps < maxSteps)
        {
          Real targetDistance = surfaceField(position);
          diff = fabs(targetDistance);
          targetDelta.normalize();
          VEC3F move = (targetDistance) * targetDelta;
          position = position - move;

          // update the direction
          targetDelta = targetGradient(position);
          steps++;
        }

        // clamp position to something reasonable
        position[0] = (position[0] < lowerCorner[0]) ? lowerCorner[0] : position[0]; 
        position[1] = (position[1] < lowerCorner[1]) ? lowerCorner[1] : position[1]; 
        position[2] = (position[2] < lowerCorner[2]) ? lowerCorner[2] : position[2]; 

        position[0] = (position[0] > upperCorner[0]) ? upperCorner[0] : position[0]; 
        position[1] = (position[1] > upperCorner[1]) ? upperCorner[1] : position[1]; 
        position[2] = (position[2] > upperCorner[2]) ? upperCorner[2] : position[2]; 
      }
}

//////////////////////////////////////////////////////////////////////
// scan through file until keyword "tiles" is found
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::scanToTiles(FILE*& file)
{
  char buffer[512];
  bool tileFound = false;
  while (!tileFound)
  {
    int success = fscanf(file, "%s", buffer);

    if (success == EOF)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Scanned off of the end of the file! Houdini 12 velocity file is not valid! " << endl;
      exit(0);
    }

    string search(buffer); 

    // look for the one that says "tiles"
    size_t findVolume = search.find("\"tiles\"");
    if (findVolume != string::npos)
      tileFound = true;
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::readTile(FILE*& file, const int storeIndex)
{
  char buffer[512];
  int xTotalTiles = (_xRes / 16);
  int yTotalTiles = (_yRes / 16);
  int zTotalTiles = (_zRes / 16);

  if (xTotalTiles * 16 != _xRes) xTotalTiles++;
  if (yTotalTiles * 16 != _yRes) yTotalTiles++;
  if (zTotalTiles * 16 != _zRes) zTotalTiles++;

  int success;
  int tilesSeen = 0;

  for (int zTile = 0; zTile < zTotalTiles; zTile++)
    for (int yTile = 0; yTile < yTotalTiles; yTile++)
      for (int xTile = 0; xTile < xTotalTiles; xTile++, tilesSeen++)
      {
        // peel off '['
        success = fscanf(file, "%s", buffer);
        if (success == EOF || success < 0) { cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl; cout << " Houdini 12 velocity file is not valid! " << endl; exit(0); }

        if (buffer[0] != '[')
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " buffer contains: " << buffer << endl;
          cout << " File was not formatted as expected! " << endl;
          exit(0);
        }

        // peel off 'compression'
        success = fscanf(file, "%s,", buffer);
        if (success == EOF || success < 0) { cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl; cout << " Houdini 12 velocity file is not valid! " << endl; exit(0); }

        // if compression type is two, the entire block is just zeros
        if (buffer[14] == '2')
        {
          success = fscanf(file, "%s,", buffer);
          success = fscanf(file, "%s,", buffer);
          continue;
        }

        // peel off whitespace and "
        char single;
        success = fscanf(file, " %c", &single);
        if (success == EOF || success < 0) { cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl; cout << " Houdini 12 velocity file is not valid! " << endl; exit(0); }
        double data;
        success = fscanf(file, "data\",");
        if (success == EOF || success < 0) { cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl; cout << " Houdini 12 velocity file is not valid! " << endl; exit(0); }

        // see if we hit a run of zeros
        success = fscanf(file, "%c", &single);
        if (single == '0')
        {
          // peel off the trailing ],
          success = fscanf(file, "%s", (char*)&buffer);
          continue;
        }

        bool closingBracketFound = false;

        // read in a block
        int totalTileSize = 16 * 16 * 16;
        vector<double> blockData(totalTileSize);
        int index = 0;
        while (!closingBracketFound)
        {
          success = fscanf(file, "%lf", &data);
          if (success == EOF || success < 0) { cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl; cout << " Houdini 12 velocity file is not valid! " << endl; exit(0); }
          blockData[index] = data;
          index++;

          // peek at the next token
          success = fscanf(file, "%c", &single);
          if (success == EOF || success < 0) { cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl; cout << " Houdini 12 velocity file is not valid! " << endl; exit(0); }

          if (single == ']')
            closingBracketFound = true;
          else if (single != ',' && single != '\n')
          {
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << single << " found!!!" << endl;
            cout << " Houdini file not formatted as expected!!!" << endl;
            exit(0);
          }
        }

        int xTileSize = 16;
        int yTileSize = 16;
        int zTileSize = 16;

        if (xTile == xTotalTiles - 1)
          xTileSize = _xRes % 16 + 1;
        if (yTile == yTotalTiles - 1)
          yTileSize = _yRes % 16 + 1;
        if (zTile == zTotalTiles - 1)
          zTileSize = _zRes % 16 + 1;

        int xTileOffset = xTile * 16;
        int yTileOffset = yTile * 16;
        int zTileOffset = zTile * 16;

        index = 0;
        for (int z = 0; z < zTileSize; z++)
          for (int y = 0; y < yTileSize; y++)
            for (int x = 0; x < xTileSize; x++, index++)
            {
              if (xTileOffset + x >= _xRes) continue;
              if (yTileOffset + y >= _yRes) continue;
              if (zTileOffset + z >= _zRes) continue;

              (*this)(xTileOffset + x, yTileOffset + y,zTileOffset + z)[storeIndex] = blockData[index];
            }

        success = fscanf(file, "%s", buffer);
        if (success == EOF || success < 0) { cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl; cout << " Houdini 12 velocity file is not valid! " << endl; exit(0); }
      }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::readHoudini12(const string filename)
{
  cout << " Reading Houdini 12 velocity file " << filename.c_str() << " ..."; flush(cout);
  int size = filename.size();
  if (filename[size - 1] == 'z' && filename[size - 2] == 'g')
  {
    cout << " File is Gzipped! " << endl;
    exit(0);
  }
  
	_lengths[0] = 2.5;
	_lengths[1] = 2.5;
	_lengths[2] = 2.5;

  if (_data) delete[] _data;
  _xRes = 100; 
  _yRes = 100; 
  _zRes = 100; 

  _totalCells = _xRes * _yRes * _zRes;
  _slabSize = _xRes * _yRes;
  _data = new VEC3F[_totalCells];
  clear();

  _dx = _lengths[0] / _xRes;
  _dy = _lengths[1] / _yRes;
  _dz = _lengths[2] / _zRes;

  FILE* file;
  file = fopen(filename.c_str(), "r");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Could not open filename " << filename.c_str() << "!!!" << endl;
    exit(0);
  }

  // scan forward until a "tiles" is seen
  cout << " scanning for x component ..."; flush(cout);
  scanToTiles(file);

  // read in the x components
  cout << " reading x component ..."; flush(cout);
  readTile(file, 0);
  
  // scan forward until a "tiles" is seen
  cout << " scanning for y component ..."; flush(cout);
  scanToTiles(file);

  // read in the y components
  cout << " reading y component ..."; flush(cout);
  readTile(file, 1);

  // scan forward until a "tiles" is seen
  cout << " scanning for z component ..."; flush(cout);
  scanToTiles(file);

  // read in the y components
  cout << " reading z component ..."; flush(cout);
  readTile(file, 2);

  // fix weirdness with the Houdini files
  fixNanInfs();

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// getting nans and infs from reading in the Houdini file -- stomp them here.
//////////////////////////////////////////////////////////////////////
void VECTOR3_FIELD_3D::fixNanInfs()
{
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++)
        for (int i = 0; i < 3; i++)
          if (isinf((*this)(x,y,z)[i]) || isnan((*this)(x,y,z)[i]))
            (*this)(x,y,z)[i] = 0;
}
