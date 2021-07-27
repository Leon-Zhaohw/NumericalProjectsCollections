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

#ifndef IWAVE_3D_H
#define IWAVE_3D_H

#include <SETTINGS.h>
#include <FIELD_3D.h>
#include <VECTOR3_FIELD_3D.h>
#include <TRIANGLE_MESH.h>
#include <BOX.h>

using namespace std;

class IWAVE_3D
{
public:
  IWAVE_3D();
  IWAVE_3D(int xRes, int yRes, int zRes, int filterWidth);
  ~IWAVE_3D();

  // accessors
  Real& alpha()   { return _alpha; };
  TRIANGLE_MESH& surfaceMesh()  { return _surfaceMesh; };
  const FIELD_3D& height() const { return _height; };

  // create the convolution kernel
  void initializeKernel();

  // load in a simulation frame from PhysBAM
  void loadPhysBAM(int frame, string path, int upResFactor = 1, const Real scaleDistance = 1);
  void loadNextPouringFrozenCore(int frame, string path, int upResFactor = 1, const Real scaleDistance = 1);

  // load in a simulation frame from Houdini
  void loadHoudini12(int frame, string path, int upResFactor = 1, const Real scaleDistance = 1);
  void loadNextHoudini12(int frame, string path, int upResFactor = 1, const Real scaleDistance = 1);

  // refresh the surface textures stored in _surfaceMesh
  void refreshSurfaceTextures();

  // set the source term to curvature, but upsample first
  void setSourceToPouringCurvature(bool frozenCore = false);
  void setSourceToHoudiniCurvature(bool frozenCore = false);

  // build the triangle filename
  string getTrianglePath(string path, int frame, int upResFactor);

  // step the simulation, using frozen core extension
  void stepPouringFrozenCore();
  void stepHoudiniFrozenCore();
  
  // add an obstacle to the list -- heights on the interior are zeroed out when stompObstacles is called.
  void addObstacle(const BOX& box) { _obstacles.push_back(box); };

  // write the simulation state out
  void writeHeight(const string& path, const int frame) const;
  void writeSurfaceMesh(const string& path, const int frame);

private:
  // simulation fields
  FIELD_3D _height;
  FIELD_3D _heightOld;
  FIELD_3D _verticalDerivative;
  FIELD_3D _source;
  VECTOR3_FIELD_3D _closestPoints;

  // the 3D iWave convolution kernel
  FIELD_3D _kernel;

  // simulation consts
  Real _dt;
  Real _alpha;
  Real _gravity;
  VEC3F _gravityDirection;

  // filter dimensions
  int _filterWidth;

  // simulation dimensions
  int _xRes;
  int _yRes;
  int _zRes;
  int _slabSize;
  int _totalCells;

  int _upResFactor;
  bool _newFrameLoaded;
  int _borderPadding;

  // input simulation distance field and velocity
  FIELD_3D _distanceField;
  FIELD_3D _distanceFieldLowRes;
  VECTOR3_FIELD_3D _velocityFieldLowRes;

  // current triangle mesh
  TRIANGLE_MESH _surfaceMesh;

  // obstacles -- height values inside them are stomped when stompObstacles is called.
  vector<BOX> _obstacles;

  // cache some temps so we don't take a memory allocation hit every time
  FIELD_3D _temp1;
  FIELD_3D _temp2;
  FIELD_3D _heightCached;
  FIELD_3D _heightOldCached;

  // currently loaded simulation frame number
  int _frameNumber;

  // catmull rom interpolation
  Real cubic(Real interp);

  // reinitialize after a resolution change
  void reinitialize(int xRes, int yRes, int zRes);
  
  // refresh the closest point values
  void refreshClosestPoints();
  
  // advect with second order MacCormack
  void advectMacCormackNarrowBand(const Real dt, const VECTOR3_FIELD_3D& velocityGrid, const FIELD_3D& distanceOld, const FIELD_3D& distanceNew, const int maxCells, const FIELD_3D& oldField, FIELD_3D& newField, FIELD_3D& temp1, FIELD_3D& temp2);
  
  // do all timestep finite differencing over the narrow band
  void timestepNarrowBand(const int maxRadius);
};

#endif
