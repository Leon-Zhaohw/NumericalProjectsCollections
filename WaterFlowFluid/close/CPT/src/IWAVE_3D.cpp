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
// 
//////////////////////////////////////////////////////////////////////

#include "IWAVE_3D.h"
#include "TIMER.h"
#include <omp.h>

//////////////////////////////////////////////////////////////////////
// Constructor / Destructor
//////////////////////////////////////////////////////////////////////
IWAVE_3D::IWAVE_3D(int xRes, int yRes, int zRes, int filterWidth) :
  _height(xRes, yRes, zRes),
  _heightOld(xRes, yRes, zRes),
  _verticalDerivative(xRes, yRes, zRes),
  _source(xRes, yRes, zRes),
  _filterWidth(filterWidth),
  _xRes(xRes),
  _yRes(yRes),
  _zRes(zRes),
  _slabSize(xRes * yRes),
  _totalCells(xRes * yRes * zRes)
{
  _dt = 0.03;
  _alpha = 0.3;
  _gravity = 9.8 * _dt * _dt;
  _gravityDirection = VEC3F(0,-1,0);
  _frameNumber = -1;

  // compute the 3D iWave kernel
  //
  // if you want to do the traditional wave equation, just swap in the 
  // 3D Laplacian here. You'll probably need to push up the damping
  // or dial back the timestep size though, as it is less stable.
  initializeKernel();
}

IWAVE_3D::IWAVE_3D() :
  _filterWidth(-1),
  _xRes(-1),
  _yRes(-1),
  _zRes(-1)
{
  _dt = 0; 
  _alpha = 0;
  _gravity = 0;
}

IWAVE_3D::~IWAVE_3D()
{
}

//////////////////////////////////////////////////////////////////////
// create the convolution kernel
//////////////////////////////////////////////////////////////////////
void IWAVE_3D::initializeKernel()
{
  assert(_filterWidth > 0);

  double dk = 0.01;
  double sigma = 1.0;
  double norm = 0;
  double startK = dk;
  double endK = 10;

  for (double freq  = startK; freq < endK; freq += dk)
    // the original iWave kernel
    norm += freq * freq * exp(-sigma * freq * freq);

  int halfWidth = _filterWidth / 2;

  _kernel = FIELD_3D(_filterWidth, _filterWidth, _filterWidth);
  for (int i = -halfWidth; i <= halfWidth; i++)
  {
    for (int j = -halfWidth; j <= halfWidth; j++)
      for (int k = -halfWidth; k <= halfWidth; k++)
      {
        double r = sqrt((float)(i * i + j * j + k * k));
        double kern = 0;

        for (double freq = startK; freq < endK; freq += dk)
        {
          double currentSinc = sin(r * freq) / (r * freq);
          kern += freq * freq * freq * exp(-sigma * freq * freq) * currentSinc;
        }

        double interp = cubic(((r /halfWidth) - 0.9) / 0.1);
        kern *= interp;

        _kernel(i + halfWidth, j + halfWidth, k + halfWidth) = kern;
      }
  }

  _kernel *= 1.0 / M_PI;
  _kernel(halfWidth, halfWidth, halfWidth) = 0;

  _kernel *= 1.0 / norm;
  FIELD_2D projection = _kernel.zProjection();
  _kernel(halfWidth, halfWidth, halfWidth) += 1.0 - projection(halfWidth, halfWidth);
}

///////////////////////////////////////////////////////////////////////
// catmull rom interpolation
///////////////////////////////////////////////////////////////////////
Real IWAVE_3D::cubic(Real interp)
{
  if (interp < 0) return 1.0;
  if (interp > 1) return 0;

  Real squared = interp * interp;
  return 2 * squared * interp - 3 * squared + 1;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
string IWAVE_3D::getTrianglePath(string path, int frame, int upResFactor)
{
	string pathSlash = path;
	if (path[path.length() - 1] != '/')
		pathSlash += string("/");
	char buffer[256];
	sprintf(buffer, "%i", frame);
	pathSlash += string(buffer) + string("/");
	sprintf(buffer, "%i", upResFactor);
	string filename = pathSlash + string("triangles.") + string(buffer) + string(".obj");

	return filename;
}

///////////////////////////////////////////////////////////////////////
// load in a simulation frame from PhysBAM
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::loadPhysBAM(int frame, string path, int upResFactor, const Real scaleDistance)
{
  TIMER functionTimer(__FUNCTION__);

  cout << "======================================" << endl;
  cout << " Loading PhysBAM frame " << frame << endl;
  cout << "======================================" << endl;

  _upResFactor = upResFactor;

  // low-res versions
  TRIANGLE_MESH::readPhysBAMFrame(frame, path, _distanceFieldLowRes, _velocityFieldLowRes, _surfaceMesh, false, scaleDistance);
  cout << " Read in PhysBAM SDF dims: " << _distanceFieldLowRes.dims() << endl;

  if (upResFactor != 1)
  {
    // upRes the distance field
    cout << " Upsampling ..."; flush(cout);
    _distanceField = _distanceFieldLowRes.resampleCubicUnclampedNarrowBand(upResFactor, _distanceFieldLowRes, 10);
	  cout << " done. " << endl;
    cout << " Distance field center: " << _distanceField.center() << endl;

    // look for an existing upRes'd mesh otherwise marching cubes it
	  string filename = getTrianglePath(path,frame,upResFactor);
    
    // try to read in a triangle mesh
    bool success = _surfaceMesh.readOBJ(filename.c_str());
    if (!success)
    {
      _surfaceMesh.computeMarchingCubes(_distanceField, true);
      _surfaceMesh.writeOBJ(filename);
    }
  }
  else
    _distanceField = _distanceFieldLowRes;

  // reset the iWave fields too
  reinitialize(_distanceField.xRes(), _distanceField.yRes(), _distanceField.zRes());

  // get the closest point values
  refreshClosestPoints();

  _height = _distanceField;
  _height.clear();

  _heightOld = _height;
  _borderPadding = 3;
  _newFrameLoaded = true;
  _frameNumber = frame;
}

///////////////////////////////////////////////////////////////////////
// load in a simulation frame from Houdini
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::loadHoudini12(int frame, string path, int upResFactor, const Real scaleDistance)
{
  TIMER functionTimer(__FUNCTION__);

  cout << "======================================" << endl;
  cout << " Loading Houdini frame " << frame << endl;
  cout << "======================================" << endl;

  _upResFactor = upResFactor;

  // low-res versions
  TRIANGLE_MESH::readHoudini12(frame, path, _distanceFieldLowRes, _velocityFieldLowRes, _surfaceMesh);

  // the simulation is of size 2, so scale back down to a unit SDF
  _distanceFieldLowRes *= 1.0 / 5.0;

  if (upResFactor != 1)
  {
    // upRes the distance field
    cout << " Upsampling ..."; flush(cout);
    _distanceField = _distanceFieldLowRes.resampleCubicUnclamped(upResFactor * _distanceFieldLowRes.xRes(), upResFactor * _distanceFieldLowRes.yRes(), upResFactor * _distanceFieldLowRes.zRes());
    cout << " done. " << endl;
    cout << " Distance field center: " << _distanceField.center() << endl;


    // look for an existing upres'd mesh otherwise marching cubes it
    cout << " Reading surface mesh ..."; flush(cout);
    string pathSlash = path;
    if (path[path.length() - 1] != '/')
      pathSlash += string("/");
    char buffer[256];
    sprintf(buffer, "%i_%i", frame, upResFactor);
    string filename = pathSlash + string("triangles.") + string(buffer) + string(".obj");

    // try to read in a triangle mesh
    bool success = _surfaceMesh.readOBJ(filename.c_str());
    if (!success)
    {
      _surfaceMesh.computeMarchingCubes(_distanceField, true);
      _surfaceMesh.writeOBJ(filename);
    }
  }
  else
    _distanceField = _distanceFieldLowRes;

  // HOUDINI SPECIFIC: reset the dx values for closest points later
  _distanceField.setDx(1.0 / (upResFactor * 100));

  // reset the iWave fields too
  reinitialize(_distanceField.xRes(), _distanceField.yRes(), _distanceField.zRes());
  refreshClosestPoints();
  
  _height = _distanceField;
  _height.setDx(1.0 / (upResFactor * 100));
  _height.clear();

  _heightOld = _height;
  _borderPadding = 3;
  _newFrameLoaded = true;
}

///////////////////////////////////////////////////////////////////////
// load in a simulation frame from Houdini
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::loadNextHoudini12(int frame, string path, int upResFactor, const Real scaleDistance)
{
  TIMER functionTimer(__FUNCTION__);
  
  cout << "======================================" << endl;
  cout << " Loading next Houdini 12 frame " << frame << endl;
  cout << "======================================" << endl;

  Real beforeMin = _height.fieldMin();
  Real beforeMax = _height.fieldMax();
  Real beforeMinOld = _heightOld.fieldMin();
  Real beforeMaxOld = _heightOld.fieldMax();

  FIELD_3D distanceFieldLowResOld(_distanceFieldLowRes);

  // low-res versions
  TRIANGLE_MESH::readHoudini12(frame, path, _distanceFieldLowRes, _velocityFieldLowRes, _surfaceMesh);

  // the simulation is of size 2, so scale back down to a unit SDF
  _distanceFieldLowRes *= 1.0 / 5.0;

  if (upResFactor != 1)
  {
    // upres the distance field
    cout << " Upsampling ..."; flush(cout);
    _distanceField = _distanceFieldLowRes.resampleCubicUnclamped(upResFactor * _distanceFieldLowRes.xRes(),
                                                                 upResFactor * _distanceFieldLowRes.yRes(),
                                                                 upResFactor * _distanceFieldLowRes.zRes());
    cout << " done. " << endl;

    // look for an existing upres'd mesh otherwise marching cubes it
    cout << " Reading surface mesh ..."; flush(cout);
    string pathSlash = path;
    if (path[path.length() - 1] != '/')
      pathSlash += string("/");
    char buffer[256];
    sprintf(buffer, "%i_%i", frame, upResFactor);
    string filename = pathSlash + string("triangles.") + string(buffer) + string(".obj");

    // try to read in a triangle mesh
    bool success = _surfaceMesh.readOBJ(filename.c_str());
    if (!success)
    {
      _surfaceMesh.computeMarchingCubes(_distanceField, true);
      _surfaceMesh.writeOBJ(filename);
    }
    cout << " done. " << endl;
  }
  else
    _distanceField = _distanceFieldLowRes;

  _heightCached = _height;
  _heightOldCached = _heightOld;

  VEC3F mins, maxs;
  _surfaceMesh.boundingBox(mins, maxs);
  cout << " Triangles bounding box: " << mins << " " << maxs << endl;

  // reset the iWave fields too
  reinitialize(_distanceField.xRes(), _distanceField.yRes(), _distanceField.zRes());

  // clear the height fields -- they've been cached locally at the top of the function
  _height = _distanceField;
  _height.clear();

  _heightOld = _height;
  _heightOld.clear();

  _temp1 = _height;
  _temp2 = _height;
  _temp1.clear();
  _temp2.clear();

  Real dt = upResFactor;
  int maxCells = 3;

  advectMacCormackNarrowBand(dt, _velocityFieldLowRes, distanceFieldLowResOld, _distanceFieldLowRes, maxCells, _heightCached, _height, _temp1, _temp2);
  advectMacCormackNarrowBand(dt, _velocityFieldLowRes, distanceFieldLowResOld, _distanceFieldLowRes, maxCells, _heightOldCached, _heightOld, _temp1, _temp2);
  cout << " done. " << endl;
 
  // clamping to previous min/max
  _height.clamp(beforeMin, beforeMax);
  _heightOld.clamp(beforeMinOld, beforeMaxOld);

  const int narrowRadius = 10;
	cout << " Computing closest points for radius " << narrowRadius << " ..."; flush(cout);
	VECTOR3_FIELD_3D::computeClosestPointsNarrowBand(_height, _distanceField, narrowRadius, _closestPoints); 
	cout << " done." << endl;

  int maxRadius = 10;
  cout << " Propagating height values after advection ... "; flush(cout);
  cout << " Using a narrow band convolution radius of " << maxRadius; flush(cout);
  _height = VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBandFrozenCore(
                                  _height, _closestPoints,
                                  _distanceField, maxRadius, 1);
  _heightOld = VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBandFrozenCore(
                                  _heightOld, _closestPoints,
                                  _distanceField, maxRadius, 1);
  cout << " done. " << endl;

  Real afterMin = _height.fieldMin();
  Real afterMax = _height.fieldMax();
  cout << " Height before advection: " << beforeMin << " " << beforeMax << endl;
  cout << " Height after advection: " << afterMin << " " << afterMax << endl;
  cout << " Height min index: " << _height.minIndex() << endl;
  cout << " Height max index: " << _height.maxIndex() << endl;

  _newFrameLoaded = true;
}

///////////////////////////////////////////////////////////////////////
// refresh the surface textures stored in _surfaceMesh
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::refreshSurfaceTextures()
{
  _surfaceMesh.textureUsingSolidTexture(_height, 3);
}

///////////////////////////////////////////////////////////////////////
// reinitialize after a resolution change
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::reinitialize(int xRes, int yRes, int zRes)
{
  TIMER functionTimer(__FUNCTION__);
  if (xRes == _xRes && yRes == _yRes && zRes == _zRes)
  {
    _height.clear();
    _heightOld.clear();
    _verticalDerivative.clear();
    _source.clear();

    return;
  }

  _xRes = xRes;
  _yRes = yRes;
  _zRes = zRes;
  _slabSize = xRes * yRes;
  _totalCells = xRes * yRes * zRes;

  _height = FIELD_3D(xRes, yRes, zRes),
  _heightOld = FIELD_3D(xRes, yRes, zRes);
  _verticalDerivative = FIELD_3D(xRes, yRes, zRes);
  _source = FIELD_3D(xRes, yRes, zRes);
}

///////////////////////////////////////////////////////////////////////
// refresh the closest point values
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::refreshClosestPoints()
{
  _closestPoints = VECTOR3_FIELD_3D::cellCenters(_height);
}

///////////////////////////////////////////////////////////////////////
// load in a simulation frame from PhysBAM
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::loadNextPouringFrozenCore(int frame, string path, int upResFactor, const Real scaleDistance)
{
  TIMER functionTimer(__FUNCTION__);
  
  cout << "==========================================" << endl;
  cout << " Loading pouring frozen core frame " << frame << " with scale distance " << scaleDistance << endl;
  cout << "==========================================" << endl;

  Real beforeMin = _height.fieldMin();
  Real beforeMax = _height.fieldMax();
  Real beforeMinOld = _heightOld.fieldMin();
  Real beforeMaxOld = _heightOld.fieldMax();

  FIELD_3D distanceFieldLowResOld(_distanceFieldLowRes);

  // low-res versions
  TRIANGLE_MESH::readPhysBAMFrame(frame, path, _distanceFieldLowRes, _velocityFieldLowRes, _surfaceMesh, false, scaleDistance);

  int maxRadius = 10;
    
  if (upResFactor != 1)
  {
    // upres the distance field
    cout << " Upsampling ..."; flush(cout);
    _distanceField = _distanceFieldLowRes.resampleCubicUnclampedNarrowBand(upResFactor, _distanceFieldLowRes, maxRadius);
    cout << " done. " << endl;

    // look for an existing upres'd mesh otherwise marching cubes it
    cout << " Reading surface mesh ..."; flush(cout);
	  string filename = getTrianglePath(path,frame,upResFactor);

    // try to read in a triangle mesh
    bool success = _surfaceMesh.readOBJ(filename.c_str());
    if (!success)
    {
       FIELD_3D marchingOnly = _distanceFieldLowRes.resampleCubicUnclamped(upResFactor * _distanceFieldLowRes.xRes(), upResFactor * _distanceFieldLowRes.yRes(), upResFactor * _distanceFieldLowRes.zRes());
      _surfaceMesh.computeMarchingCubes(marchingOnly, true);
      _surfaceMesh.writeOBJ(filename);
    }
    cout << " done. " << endl;
  }
  else
    _distanceField = _distanceFieldLowRes;
 

  _heightCached = _height;
  _heightOldCached = _heightOld;

  VEC3F mins, maxs;
  _surfaceMesh.boundingBox(mins, maxs);
  cout << " Triangles bounding box: " << mins << " " << maxs << endl;

  // reset the iWave fields too
  reinitialize(_distanceField.xRes(), _distanceField.yRes(), _distanceField.zRes());

  // clear the height fields -- they've been cached locally at the top of the function
  _height = _distanceField;
  _height.clear();

  _heightOld = _height;
  _heightOld.clear();

  _temp1 = _height;
  _temp2 = _height;
  _temp1.clear();
  _temp2.clear();

  Real dt = upResFactor;

  int maxCells = 3; 

  advectMacCormackNarrowBand(dt, _velocityFieldLowRes, distanceFieldLowResOld, _distanceFieldLowRes, maxCells, _heightCached, _height, _temp1, _temp2);
  advectMacCormackNarrowBand(dt, _velocityFieldLowRes, distanceFieldLowResOld, _distanceFieldLowRes, maxCells, _heightOldCached, _heightOld, _temp1, _temp2);

  // clamping to previous min/max
  _height.clamp(beforeMin, beforeMax);
  _heightOld.clamp(beforeMinOld, beforeMaxOld);
  cout << " done. " << endl;

  const int narrowRadius = 10;
	cout << " Computing closest points for radius " << narrowRadius << " ..."; flush(cout);
	VECTOR3_FIELD_3D::computeClosestPointsNarrowBand(_height, _distanceField, narrowRadius, _closestPoints); 
	cout << " done." << endl;

  cout << " Propagating height values after advection ... "; flush(cout);
  cout << " Using a narrow band convolution radius of " << maxRadius; flush(cout);
  _height = VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBandFrozenCore(
                                  _height, _closestPoints,
                                  _distanceField, maxRadius, 1);
  _heightOld = VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBandFrozenCore(
                                  _heightOld, _closestPoints,
                                  _distanceField, maxRadius, 1);
  cout << " done. " << endl;

  Real afterMin = _height.fieldMin();
  Real afterMax = _height.fieldMax();
  cout << " Height before advection: " << beforeMin << " " << beforeMax << endl;
  cout << " Height after advection: " << afterMin << " " << afterMax << endl;
  cout << " Height min index: " << _height.minIndex() << endl;
  cout << " Height max index: " << _height.maxIndex() << endl;

  _newFrameLoaded = true;
  _frameNumber = frame;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::setSourceToHoudiniCurvature(bool frozenCore)
{
  TIMER functionTimer(__FUNCTION__);

  // there are two ways to go here -- either upres the curvature directly, or the 
  // distance field, and then compute the curavture at the higher res.
  //
  // The first option gives smoother results, as it additionally low-pass filters
  // the quantity of interest, i.e. the curvature. If we want slightly more turbulent-looking
  // results, perhaps the second option can be tried.
 
  FIELD_3D maxCurvatureLowRes(_distanceFieldLowRes);
  FIELD_3D minCurvatureLowRes(_distanceFieldLowRes);
  _distanceFieldLowRes.principalCurvatures(minCurvatureLowRes, maxCurvatureLowRes);

  // It also doesn't hurt that computing the curvature at the lower res is faster.
  cout << " Upsampling curvature ... "; flush(cout);
  FIELD_3D maxCurvature = maxCurvatureLowRes.resampleCubicUnclampedNarrowBand(_upResFactor, _distanceFieldLowRes, 10);
  FIELD_3D minCurvature = minCurvatureLowRes.resampleCubicUnclampedNarrowBand(_upResFactor, _distanceFieldLowRes, 10);
  cout << " done." << endl;
  int res = maxCurvatureLowRes.maxRes();

  // set the curvature target to Nyquist, and cutoff near half of Nyquist
  Real targetBand = res * 0.25;
  Real targetWidth = res * 0.2;
  maxCurvature.softBandPass(targetBand, targetWidth);

  // clamp away spurious curvatures that are not feasible at this
  // grid resolution
  maxCurvature *= 1.0 / targetBand;
  maxCurvature.clamp(-targetBand, targetBand);

  cout << " Using Gaussian curvature for seeding." << endl;
  minCurvature.clamp(-1,1);
  _source = maxCurvature;

  // convert what's left to Gaussian curvature
  _source *= minCurvature;
  _source *= 0.1;

  Real stencilRadius = 2.0 * sqrt(2.0);
  cout << " Sourcing with a narrow band stencil radius of " << stencilRadius << endl;
  vector<int> stencilNarrowBand = _distanceField.computeNarrowBand(stencilRadius);

  // make it an extension field
  if (frozenCore)
  {
    int fullRadius = 10;
    cout << " Sourcing with a narrow band radius of " << fullRadius << endl;
    _source = VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBandFrozenCore(_source, _closestPoints, 
        _distanceField, fullRadius, 1);
  }
  else
  {
    cout << " Sourcing with FULL GRID " << endl;
    _source = VECTOR3_FIELD_3D::setToClosestPointValues(_source, _closestPoints);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::setSourceToPouringCurvature(bool frozenCore)
{
  TIMER functionTimer(__FUNCTION__);

  // there are two ways to go here -- either upres the curvature directly, or the 
  // distance field, and then compute the curavture at the higher res.
  //
  // The first option gives smoother results, as it additionally low-pass filters
  // the quantity of interest, i.e. the curvature. If we want slightly more turbulent-looking
  // results, perhaps the second option can be tried.
 
  FIELD_3D maxCurvatureLowRes(_distanceFieldLowRes);
  FIELD_3D minCurvatureLowRes(_distanceFieldLowRes);
  _distanceFieldLowRes.principalCurvatures(minCurvatureLowRes, maxCurvatureLowRes);

  // It also doesn't hurt that computing the curvature at the lower res is faster.
  cout << " Upsampling curvature ... "; flush(cout);
  FIELD_3D maxCurvature = maxCurvatureLowRes.resampleCubicUnclampedNarrowBand(_upResFactor, _distanceFieldLowRes, 10);
  FIELD_3D minCurvature = minCurvatureLowRes.resampleCubicUnclampedNarrowBand(_upResFactor, _distanceFieldLowRes, 10);
  cout << " done." << endl;
  int res = maxCurvatureLowRes.maxRes();

  // set the curvature target to Nyquist, and cutoff near half of Nyquist
  Real targetBand = res * 0.5;
  Real targetWidth = res * 0.2;
  maxCurvature.softBandPass(targetBand, targetWidth);

  // clamp away spurious curvatures that are not feasible at this
  // grid resolution
  maxCurvature *= 1.0 / targetBand;
  maxCurvature.clamp(-targetBand, targetBand);

  cout << " Using Gaussian curvature for seeding." << endl;
  minCurvature.clamp(-1,1);
  _source = maxCurvature;

  // convert what's left to Gaussian curvature
  _source *= minCurvature;
  _source *= 0.1;

  Real stencilRadius = 2.0 * sqrt(2.0);
  cout << " Sourcing with a narrow band stencil radius of " << stencilRadius << endl;
  vector<int> stencilNarrowBand = _distanceField.computeNarrowBand(stencilRadius);

  // make it an extension field
  if (frozenCore)
  {
    int fullRadius = 10;
    cout << " Sourcing with a narrow band radius of " << fullRadius << endl;
    _source = VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBandFrozenCore(_source, _closestPoints, 
        _distanceField, fullRadius, 1);
  }
  else
  {
    cout << " Sourcing with FULL GRID " << endl;
    _source = VECTOR3_FIELD_3D::setToClosestPointValues(_source, _closestPoints);
  }
}

///////////////////////////////////////////////////////////////////////
// do all timestep finite differencing over the narrow band
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::timestepNarrowBand(const int maxRadius)
{
  TIMER functionTimer(__FUNCTION__);
  const int xRes = _height.xRes();
  const int yRes = _height.yRes();
  const int zRes = _height.zRes();

  const bool sourceInit = _source.initialized();
  const Real adt = _alpha * _dt;
  const Real adt2 = 1.0 / (1.0 + adt);
  const Real twoMinusAdt = 2.0 - adt;
  const Real invDx = 1.0 / _distanceField.dx();

#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        Real currentDistance = fabs(_distanceField(x,y,z) * invDx);
        if (currentDistance > maxRadius)
          continue;

        Real& height = _height(x,y,z);
        Real& heightOld = _heightOld(x,y,z);

        Real temp = _height(x,y,z);
        height *= twoMinusAdt;
        height -= heightOld;
        height -= _gravity * _verticalDerivative(x,y,z);
        height *= adt2;

        // check if the source has any entry
        if (!sourceInit) continue;
        
        Real& source = _source(x,y,z);
        height += source;
        temp += source;
        heightOld = temp;
      }

  if (sourceInit)
    _source.clear();
}

///////////////////////////////////////////////////////////////////////
// step the simulation, with per frame extensions, with a frozen core
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::stepPouringFrozenCore()
{
  TIMER functionTimer(__FUNCTION__);

  int stompWidth = _kernel.xRes() / 2 + 2;

  cout << " Using Dirichlet boundaries" << endl;
  _source.stompBorder(stompWidth);
  _height.stompBorder(stompWidth);
  _heightOld.stompBorder(stompWidth);

  // build the stencil narrow band
  Real stencilRadius = 2.0 * sqrt(2.0);
  _verticalDerivative = _height.convolveNarrowBandFast15(_kernel, _distanceField, stencilRadius);

  cout << " Using damping alpha: " << _alpha << endl;
  timestepNarrowBand(10);

  // do the narrow band extension
  int maxRadius = 10;
  cout << " Using a narrow band convolution radius of " << maxRadius << endl;
  _height = VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBandFrozenCore(
                                  _height, _closestPoints,
                                  _distanceField, maxRadius, 1);
}

///////////////////////////////////////////////////////////////////////
// step the simulation, with per frame extensions, with a frozen core
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::stepHoudiniFrozenCore()
{
  TIMER functionTimer(__FUNCTION__);

  // set the border values to zero
  cout << " Using Dirichlet boundaries" << endl;
  int stompWidth = _upResFactor * 4;
  _source.stompBorder(stompWidth);
  _height.stompBorder(stompWidth);
  _heightOld.stompBorder(stompWidth);

  cout << " Using alpha: " << _alpha << endl;

  // build the stencil narrow band
  Real stencilRadius = 2.0 * sqrt(2.0);

  // convolve the height field with the iWave kernel
  // to obtain the "vertical derivative"
  _verticalDerivative = _height.convolveNarrowBandFast15(_kernel, _distanceField, stencilRadius);

  // timestep using the new vertical derivative
  timestepNarrowBand(10);

  // do the frozen core narrow band extension
  //
  // note that the closest points were already cached at the beginning
  // of the timestep in loadNextHoudini12
  int maxRadius = 10;
  cout << " Using a narrow band convolution radius of " << maxRadius << endl;
  _height = VECTOR3_FIELD_3D::setToClosestPointValuesNarrowBandFrozenCore(
                                  _height, _closestPoints,
                                  _distanceField, maxRadius, 1);

  // make sure values don't get out of control
  cout << " Clamp and stomp ...";flush(cout);
  _height.clamp(-1, 1);
  _heightOld.clamp(-1, 1);

  // stomp values outside the narrow band to zero to get better 
  // file compression. Otherwise the ghost of previous timesteps
  // just stick around and make files bigger for no reason
  _height.stompOutsideNarrowBand(_distanceField, 10); 
  _heightOld.stompOutsideNarrowBand(_distanceField, 10); 
  cout << " done." << endl;
}

///////////////////////////////////////////////////////////////////////
// write the current height field out
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::writeHeight(const string& path, const int frame) const
{
  // build the filename
  string filename = path;

  // check for the slash
  if (filename[filename.length() - 1] != '/')
    filename = filename + string("/");

  // add the frame number
  filename = filename + string("height.");
  char buffer[256];
  sprintf(buffer, "%04i", frame);
  filename = filename + string(buffer) + string(".field3d");

  _height.writeGz(filename);
}

///////////////////////////////////////////////////////////////////////
// write the current surface mesh out to PBRT
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::writeSurfaceMesh(const string& path, const int frame)
{
	// build the filename
	string filename = path;

	// check for the slash
	if (filename[filename.length() - 1] != '/')
		filename = filename + string("/");

	// add the frame number
	filename = filename + string("surface.");
	char buffer[256];
	sprintf(buffer, "%04i", frame);
	filename = filename + string(buffer) + string(".pbrt");

	_surfaceMesh.writePBRT(filename.c_str());
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void IWAVE_3D::advectMacCormackNarrowBand(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& distanceOld, const FIELD_3D& distanceNew, const int maxCells, const FIELD_3D& oldField, FIELD_3D& newField, FIELD_3D& temp1, FIELD_3D& temp2)
{
  TIMER functionTimer(__FUNCTION__);
  const int xRes = oldField.xRes();
  const int yRes = oldField.yRes();
  const int zRes = oldField.zRes();
  const int slabSize = oldField.slabSize();
  const Real invDx = 1.0 / oldField.dx();

  // track what cell to do the closest point transform on
  FIELD_3D& forwardExtend = temp1;
  forwardExtend.clear();
  Real* forwardExtendData = forwardExtend.data();

  FIELD_3D& backwardExtend = temp2;
  backwardExtend.clear();
  Real* backwardExtendData = backwardExtend.data();

  // build a list of cells to populate the extension for
  TIMER forwardTimer("Forward band");
  cout << " Computing forward extension band ..."; flush(cout);
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        Real currentDistance = fabs(distanceNew(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        // backtrace, store those cells for extension later.
        const VEC3F velocity = velocityGrid(oldField.cellCenter(x,y,z));
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

        const int z0slabSize = z0 * slabSize;
        const int z1slabSize = z1 * slabSize;
        const int y0xRes = y0 * xRes;
        const int y1xRes = y1 * xRes;

        const int i000 = x0 + y0xRes + z0slabSize;
        const int i010 = x0 + y1xRes + z0slabSize;
        const int i100 = x1 + y0xRes + z0slabSize;
        const int i110 = x1 + y1xRes + z0slabSize;
        const int i001 = x0 + y0xRes + z1slabSize;
        const int i011 = x0 + y1xRes + z1slabSize;
        const int i101 = x1 + y0xRes + z1slabSize;
        const int i111 = x1 + y1xRes + z1slabSize;

#pragma omp atomic
        forwardExtendData[i000] += 1;
#pragma omp atomic
        forwardExtendData[i001] += 1;
#pragma omp atomic
        forwardExtendData[i010] += 1;
#pragma omp atomic
        forwardExtendData[i011] += 1;
#pragma omp atomic
        forwardExtendData[i100] += 1;
#pragma omp atomic
        forwardExtendData[i101] += 1;
#pragma omp atomic
        forwardExtendData[i110] += 1;
#pragma omp atomic
        forwardExtendData[i111] += 1;
      }

  // see what cells show up in the backward advection
  TIMER backwardTimer("Backward band");
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // is the cell in the narrow band?
        Real currentDistance = fabs(distanceNew(oldField.cellCenter(x,y,z)) * invDx);
        if (currentDistance > maxCells)
          continue;

        // backtrace, store those cells for extension later.
        const VEC3F velocity = velocityGrid(oldField.cellCenter(x,y,z));
        Real xTrace = x + dt * velocity[0];
        Real yTrace = y + dt * velocity[1];
        Real zTrace = z + dt * velocity[2];

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

        const int z0slabSize = z0 * slabSize;
        const int z1slabSize = z1 * slabSize;
        const int y0xRes = y0 * xRes;
        const int y1xRes = y1 * xRes;

        const int i000 = x0 + y0xRes + z0slabSize;
        const int i010 = x0 + y1xRes + z0slabSize;
        const int i100 = x1 + y0xRes + z0slabSize;
        const int i110 = x1 + y1xRes + z0slabSize;
        const int i001 = x0 + y0xRes + z1slabSize;
        const int i011 = x0 + y1xRes + z1slabSize;
        const int i101 = x1 + y0xRes + z1slabSize;
        const int i111 = x1 + y1xRes + z1slabSize;

#pragma omp atomic
        backwardExtendData[i000] += 1;
#pragma omp atomic
        backwardExtendData[i001] += 1;
#pragma omp atomic
        backwardExtendData[i010] += 1;
#pragma omp atomic
        backwardExtendData[i011] += 1;
#pragma omp atomic
        backwardExtendData[i100] += 1;
#pragma omp atomic
        backwardExtendData[i101] += 1;
#pragma omp atomic
        backwardExtendData[i110] += 1;
#pragma omp atomic
        backwardExtendData[i111] += 1;
      }

  // what cells need to be extended so that correct values appear 
  // upon backwards advection?
  TIMER backwardExtensionTimer("Backward extension band");
#pragma omp parallel
#pragma omp for schedule(dynamic)
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++)
      {
        // add it to the list of final
        int index = x + y * xRes + z * slabSize;

        if (backwardExtendData[index] < 0.5) 
          continue;

        // backtrace, store those cells for extension later.
        const VEC3F velocity = velocityGrid(oldField.cellCenter(x,y,z));
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

        const int z0slabSize = z0 * slabSize;
        const int z1slabSize = z1 * slabSize;
        const int y0xRes = y0 * xRes;
        const int y1xRes = y1 * xRes;

        const int i000 = x0 + y0xRes + z0slabSize;
        const int i010 = x0 + y1xRes + z0slabSize;
        const int i100 = x1 + y0xRes + z0slabSize;
        const int i110 = x1 + y1xRes + z0slabSize;
        const int i001 = x0 + y0xRes + z1slabSize;
        const int i011 = x0 + y1xRes + z1slabSize;
        const int i101 = x1 + y0xRes + z1slabSize;
        const int i111 = x1 + y1xRes + z1slabSize;

#pragma omp atomic
        forwardExtendData[i000] += 1;
#pragma omp atomic
        forwardExtendData[i001] += 1;
#pragma omp atomic
        forwardExtendData[i010] += 1;
#pragma omp atomic
        forwardExtendData[i011] += 1;
#pragma omp atomic
        forwardExtendData[i100] += 1;
#pragma omp atomic
        forwardExtendData[i101] += 1;
#pragma omp atomic
        forwardExtendData[i110] += 1;
#pragma omp atomic
        forwardExtendData[i111] += 1;
      }

  // add the backward extends to the forward
  forwardExtend += backwardExtend;

  // perform Nacelle on all the tagged extension cells
  FIELD_3D extendedOld = VECTOR3_FIELD_3D::computeExtensionFieldMasked(distanceOld, oldField, forwardExtend); 
  
	FIELD_3D& phiN    = extendedOld;
	FIELD_3D& phiN1   = newField;
	FIELD_3D& phiHatN  = temp1;
	FIELD_3D& phiHatN1 = temp2;
  temp1.clear();
  temp2.clear();

  cout << " advecting narrow band forward ... ";flush(cout);
  VECTOR3_FIELD_3D::advectNarrowBand(dt, velocityGrid, phiN, phiHatN1, distanceNew, maxCells);
  cout << " advecting narrow band backward ... ";flush(cout);
  VECTOR3_FIELD_3D::advectNarrowBand(-1.0 * dt, velocityGrid, phiHatN1, phiHatN, distanceNew, maxCells);

  TIMER arithmeticTimer("Arithmetic field ops");
  phiN1.clear();
  phiN1 += phiN;
  phiN1 -= phiHatN;
  phiN1 *= 0.5;
  phiN1 += phiHatN1;

  cout << " clamping ... ";flush(cout);

	// clamp any newly created extrema
  VECTOR3_FIELD_3D::clampExtrema(dt, velocityGrid, extendedOld, newField);

	// if the error estimate was bad, revert to first order
  VECTOR3_FIELD_3D::clampOutsideRays(dt, velocityGrid, extendedOld, phiHatN1, newField);

  cout << " done." << endl;
}
