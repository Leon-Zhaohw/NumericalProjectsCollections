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
#include "TRIANGLE_MESH.h"
#include <algorithm>
#include <fstream>
#include "TIMER.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_math.h>
#include "NONLINEAR_SLICE_CUDA.h"
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::TRIANGLE_MESH() 
{
}

//////////////////////////////////////////////////////////////////////
// do a non-linear marching cubes on just two slabs at a shot
// write out the vertices first, delete them, and then output the faces
//
// This one can currently handle the largest meshes
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::TRIANGLE_MESH(const VEC3F& center, 
                             const VEC3F& lengths, 
                             const VEC3I& res, 
                             const POLYNOMIAL_4D& top, 
                             const POLYNOMIAL_4D& bottom, 
                             const Real expScaling, 
                             const int maxIterations, 
                             const Real slice, 
                             const Real isosurface, 
                             const QUATERNION& rotation, 
                             const string& cacheFilename,
                             const string& objFilename,
                             const double checkpointFrequency) :
  _res(res),
  _lengths(lengths), 
  _center(center), 
  _dxs(lengths[0] / res[0], lengths[1] / res[1], lengths[2] / res[2]),
  _fieldRotation(rotation),
  _checkpointFrequency(checkpointFrequency)
{
  _top = top;
  _bottom = bottom;
  _expScaling = expScaling;
  _maxIterations = maxIterations;
  _quaternionSlice = slice;
  _isosurface = isosurface;
  _escapeRadius = 2000.0;
  _cacheFilename = cacheFilename;

  computeNonlinearMarchingCubesCudaStreamingCheckpointed64gz(objFilename);
}

//////////////////////////////////////////////////////////////////////
// destructor
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::~TRIANGLE_MESH()
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::initSliceCuda()
{
  TIMER functionTimer(__FUNCTION__);

  // copy over the dimensions
  setConsts(_res[0], _res[1], _res[2],
            _dxs[0], _dxs[1], _dxs[2],
            _center[0], _center[1], _center[2],
            _lengths[0], _lengths[1], _lengths[2],
            _escapeRadius, _maxIterations, _expScaling, _quaternionSlice, _isosurface);

  // copy over the root positions
  double4* topRoots = new double4[_top.totalRoots()];
  for (int x = 0; x < _top.totalRoots(); x++)
  {
    topRoots[x].x = _top.roots()[x].x();
    topRoots[x].y = _top.roots()[x].y();
    topRoots[x].z = _top.roots()[x].z();
    topRoots[x].w = _top.roots()[x].w();
  }
  setTopRoots(topRoots, _top.totalRoots());
  delete[] topRoots;

  double4* bottomRoots = new double4[_bottom.totalRoots()];
  for (int x = 0; x < _bottom.totalRoots(); x++)
  {
    bottomRoots[x].x = _bottom.roots()[x].x();
    bottomRoots[x].y = _bottom.roots()[x].y();
    bottomRoots[x].z = _bottom.roots()[x].z();
    bottomRoots[x].w = _bottom.roots()[x].w();
  }
  setBottomRoots(bottomRoots, _bottom.totalRoots());
  delete[] bottomRoots;

  // copy over the powers
  double* topPowers = new double[_top.totalRoots()];
  for (int x = 0; x < _top.totalRoots(); x++)
      topPowers[x] = _top.powers()[x] * _top.powerScalar();
  setTopPowers(topPowers, _top.totalRoots());
  delete[] topPowers;
  
  double* bottomPowers = new double[_bottom.totalRoots()];
  for (int x = 0; x < _bottom.totalRoots(); x++)
      bottomPowers[x] = _bottom.powers()[x] * _bottom.powerScalar();
  setBottomPowers(bottomPowers, _bottom.totalRoots());
  delete[] bottomPowers;

  // allocate the field memory on the device
  checkCudaErrors(cudaMalloc(&_fieldCurrentCuda, sizeof(double) * _xRes * _yRes));
  getLastCudaError("Malloc failed");
  checkCudaErrors(cudaMalloc(&_fieldOldCuda, sizeof(double) * _xRes * _yRes));
  getLastCudaError("Malloc failed");
 
  checkCudaErrors(cudaMalloc(&_occupanciesCuda, sizeof(int) * _xRes * _yRes));
  getLastCudaError("Malloc failed");
  checkCudaErrors(cudaMalloc(&_compactedFlagsCuda, sizeof(int) * _xRes * _yRes));
  getLastCudaError("Malloc failed");
  checkCudaErrors(cudaMalloc(&_indexTranslationCuda, sizeof(int) * _xRes * _yRes));
  getLastCudaError("Malloc failed");
  
  checkCudaErrors(cudaMalloc(&_flagsCuda, sizeof(unsigned char) * _xRes * _yRes));
  getLastCudaError("Malloc failed");

  checkCudaErrors(cudaHostAlloc(&_flagsPinned, sizeof(unsigned char) * _xRes * _yRes, 0));

  // initialize to zero, so the padding sees it once and for all
  checkCudaErrors(cudaMemset(_occupanciesCuda, 0, sizeof(int) * _xRes * _yRes));
}

//////////////////////////////////////////////////////////////////////
// readback the contents of _fieldCurrentCuda into "field"
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::readbackCudaSlice(FIELD_2D& field)
{
  TIMER functionTimer(__FUNCTION__);

  // pull the data back from the GPU
  field.resizeAndWipe(_res[0], _res[1], VEC3F(0,0,0), VEC3F(1,1,1));
  field = DBL_MAX;
  int widthInBytes = sizeof(double)*_res[0];
  int height = _res[1];
  double* fieldHost = field.data();
  checkCudaErrors(cudaMemcpy(fieldHost, _fieldCurrentCuda, widthInBytes * height, cudaMemcpyDeviceToHost));
  getLastCudaError("Memcpy failed");
}

//////////////////////////////////////////////////////////////////////
// readback the contents of _flagsCuda into "flags"
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::readbackCudaFlags(unsigned char* flags)
{
  TIMER functionTimer(__FUNCTION__);

  // pull the data back from the GPU
  int widthInBytes = sizeof(unsigned char)*_xRes;
  checkCudaErrors(cudaMemcpy(flags, _flagsCuda, widthInBytes * _yRes, cudaMemcpyDeviceToHost));
  getLastCudaError("Memcpy failed");
}

//////////////////////////////////////////////////////////////////////
// swap the two CUDA buffers we're using to ping-pong
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::swapCudaBuffers()
{
  // ping-pong the buffers
  double* temp = _fieldCurrentCuda;
  _fieldCurrentCuda = _fieldOldCuda;
  _fieldOldCuda = temp;
}

//////////////////////////////////////////////////////////////////////
// support function for computeNonlinearMarchingCubesLowMemory(), 
// which computes a single slice of the potential function
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeNonlinearSliceCuda(const int z)
{
  TIMER functionTimer(__FUNCTION__);

  static bool firstCall = true;

  if (firstCall)
  {
      initSliceCuda();
      firstCall = false;
  }

  // do the call
  computeSliceOnGPU(_res[0], _res[1], z, _fieldCurrentCuda);
}

//////////////////////////////////////////////////////////////////////
// try to read in an edge cache
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readEdgeCache(vector<pair<int, int> >& flags)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  file = fopen(_cacheFilename.c_str(), "rb");
  if (file == NULL)
  {
    cout << " No cache file found: " << _cacheFilename.c_str() << " " << endl;
    return false;
  }
  cout << " Cache file found: " << _cacheFilename.c_str() << endl;
 
  // read in: map<pair<int, int>, bool> _vertexPairs;
  int totalPairs;
  fread((void*)&(totalPairs), sizeof(int), 1, file);

  _vertexPairs.clear();
  for (int x = 0; x < totalPairs; x++)
  {
    pair<int,int> vertexPair;
    fread((void*)&(vertexPair.first), sizeof(int), 1, file);
    fread((void*)&(vertexPair.second), sizeof(int), 1, file);
    _vertexPairs[vertexPair] = true;
  }

  // read in: vector<pair<int, int> > flags
  int totalFlags;
  fread((void*)&(totalFlags), sizeof(int), 1, file);

  flags.clear();
  for (int x = 0; x < totalFlags; x++)
  {
    pair<int, int> flagPair;
    fread((void*)&(flagPair.first), sizeof(int), 1, file);
    fread((void*)&(flagPair.second), sizeof(int), 1, file);
    flags.push_back(flagPair);
  }

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// try to read in an edge cache from a partially completed run
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readPartialInterpolationCache(int& totalChunks, int& lastCompleted)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  string partialFilename = _cacheFilename + string(".interp.partial");
  file = fopen(partialFilename.c_str(), "rb");
  if (file == NULL)
  {
    cout << " No cache file found: " << partialFilename.c_str() << " " << endl;
    return false;
  }
  cout << " Cache file found: " << partialFilename.c_str() << endl;
  int temp;
  int verticesCompleted;
  fread((void*)&temp, sizeof(int), 1, file);
  totalChunks = temp;
  fread((void*)&temp, sizeof(int), 1, file);
  lastCompleted = temp;
  fread((void*)&verticesCompleted, sizeof(int), 1, file);
  cout << " Got to chunk " << lastCompleted << " of " << totalChunks << " last time" << endl;
 
  cout << " Reading in " << verticesCompleted << " vertices from last time " << endl;
  for (int x = 0; x < verticesCompleted; x++)
  {
    VEC3F second;
    for (int y = 0; y < 3; y++)
    {
      double temp;
      fread((void*)&temp, sizeof(double), 1, file);
      second[y] = temp;
    }
    _vertices[x] = second;
  }

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// try to read in an edge cache from a partially completed run
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readPartialEdgeCache(int& zCurrent, vector<int>& flags, vector<VEC3I>& indices)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  string partialFilename = _cacheFilename + string(".partial");
  file = fopen(partialFilename.c_str(), "rb");
  if (file == NULL)
  {
    cout << " No cache file found: " << partialFilename.c_str() << " " << endl;
    return false;
  }
  cout << " Cache file found: " << partialFilename.c_str() << endl;
  fread((void*)&zCurrent, sizeof(int), 1, file);
  cout << " Got to slice " << zCurrent << " of " << _zRes << " (" << 100.0 * (float)zCurrent / _zRes << "%) last time" << endl;
 
  // read in: map<pair<int, int>, bool> _vertexPairs;
  int totalTriplets;
  fread((void*)&(totalTriplets), sizeof(int), 1, file);
  cout << " Reading in " << totalTriplets << " triplets " << endl;

  _vertexTriplets.clear();
  for (int x = 0; x < totalTriplets; x++)
  {
    pair<VEC3I, VEC3I> vertexTriplet;
    VEC3I first;
    VEC3I second;
    for (int y = 0; y < 3; y++)
      fread((void*)&(first[y]), sizeof(int), 1, file);
    for (int y = 0; y < 3; y++)
      fread((void*)&(second[y]), sizeof(int), 1, file);

    vertexTriplet.first = first;
    vertexTriplet.second = second;
    _vertexTriplets[vertexTriplet] = true;
    if (x % (int)(totalTriplets / 10) == 0)
      cout << 100 * ((Real)x / totalTriplets) << "% " << flush;
  }
  cout << " done." << endl;
  cout << " Stored " << _vertexTriplets.size() << " triplets" << endl;

  // read in: vector<pair<int, int> > flags
  int totalFlags;
  fread((void*)&(totalFlags), sizeof(int), 1, file);
  cout << " Reading in " << totalFlags << " flags " << endl;

  flags.resize(totalFlags);
  indices.resize(totalFlags);
  for (int x = 0; x < totalFlags; x++)
  {
    pair<int, VEC3I> flagPair;
    fread((void*)&(flagPair.first), sizeof(int), 1, file);
    VEC3I& second = flagPair.second;
    for (int y = 0; y < 3; y++)
      fread((void*)&(second[y]), sizeof(int), 1, file);

    flags[x] = flagPair.first;
    indices[x] = flagPair.second;
    if (x % (int)(totalFlags / 10) == 0)
      cout << 100 * ((Real)x / totalFlags) << "% " << flush;
  }
  cout << " done." << endl;

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// try to read in an edge cache
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readEdgeCacheHuge(vector<int>& flags, vector<VEC3I>& indices)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  file = fopen(_cacheFilename.c_str(), "rb");
  if (file == NULL)
  {
    cout << " No cache file found: " << _cacheFilename.c_str() << " " << endl;
    return false;
  }
  cout << " Cache file found: " << _cacheFilename.c_str() << endl;
 
  // read in: map<pair<int, int>, bool> _vertexPairs;
  int totalTriplets;
  fread((void*)&(totalTriplets), sizeof(int), 1, file);
  cout << " Reading in " << totalTriplets << " triplets " << endl;

  _vertexTriplets.clear();
  for (int x = 0; x < totalTriplets; x++)
  {
    pair<VEC3I, VEC3I> vertexTriplet;
    VEC3I first;
    VEC3I second;
    for (int y = 0; y < 3; y++)
      fread((void*)&(first[y]), sizeof(int), 1, file);
    for (int y = 0; y < 3; y++)
      fread((void*)&(second[y]), sizeof(int), 1, file);

    vertexTriplet.first = first;
    vertexTriplet.second = second;
    _vertexTriplets[vertexTriplet] = true;
  }
  cout << " Stored " << _vertexTriplets.size() << " triplets" << endl;

  // read in: vector<pair<int, int> > flags
  int totalFlags;
  fread((void*)&(totalFlags), sizeof(int), 1, file);
  cout << " Reading in " << totalFlags << " flags " << endl;

  flags.clear();
  indices.clear();
  for (int x = 0; x < totalFlags; x++)
  {
    pair<int, VEC3I> flagPair;
    fread((void*)&(flagPair.first), sizeof(int), 1, file);
    VEC3I& second = flagPair.second;
    for (int y = 0; y < 3; y++)
      fread((void*)&(second[y]), sizeof(int), 1, file);
    flags.push_back(flagPair.first);
    indices.push_back(flagPair.second);
  }

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// try to read in an edge cache
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readEdgeCacheHuge(vector<pair<int, VEC3I> >& flags)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  file = fopen(_cacheFilename.c_str(), "rb");
  if (file == NULL)
  {
    cout << " No cache file found: " << _cacheFilename.c_str() << " " << endl;
    return false;
  }
  cout << " Cache file found: " << _cacheFilename.c_str() << endl;
 
  // read in: map<pair<int, int>, bool> _vertexPairs;
  int totalTriplets;
  fread((void*)&(totalTriplets), sizeof(int), 1, file);
  cout << " Reading in " << totalTriplets << " triplets " << endl;

  _vertexTriplets.clear();
  for (int x = 0; x < totalTriplets; x++)
  {
    pair<VEC3I, VEC3I> vertexTriplet;
    VEC3I first;
    VEC3I second;
    for (int y = 0; y < 3; y++)
      fread((void*)&(first[y]), sizeof(int), 1, file);
    for (int y = 0; y < 3; y++)
      fread((void*)&(second[y]), sizeof(int), 1, file);

    vertexTriplet.first = first;
    vertexTriplet.second = second;
    _vertexTriplets[vertexTriplet] = true;
  }
  cout << " Stored " << _vertexTriplets.size() << " triplets" << endl;

  // read in: vector<pair<int, int> > flags
  int totalFlags;
  fread((void*)&(totalFlags), sizeof(int), 1, file);
  cout << " Reading in " << totalFlags << " flags " << endl;

  flags.clear();
  for (int x = 0; x < totalFlags; x++)
  {
    pair<int, VEC3I> flagPair;
    fread((void*)&(flagPair.first), sizeof(int), 1, file);
    VEC3I& second = flagPair.second;
    for (int y = 0; y < 3; y++)
      fread((void*)&(second[y]), sizeof(int), 1, file);
    flags.push_back(flagPair);
  }

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// write out an edge cache
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writeEdgeCache(const vector<pair<int, int> >& flags)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  file = fopen(_cacheFilename.c_str(), "wb");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FAILED TO CACHE OUT EDGES: " << _cacheFilename.c_str() << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    return;
  }

  cout << " Writing out cache file " << _cacheFilename.c_str() << " ... " << flush;
 
  // write out: map<pair<int, int>, bool> _vertexPairs;
  int totalPairs = _vertexPairs.size(); 
  fwrite((void*)&(totalPairs), sizeof(int), 1, file);

  map<pair<int, int>, bool>::iterator iter;
  for (iter = _vertexPairs.begin(); iter != _vertexPairs.end(); iter++)
  {
    pair<int,int> vertexPair = iter->first;
    fwrite((void*)&(vertexPair.first), sizeof(int), 1, file);
    fwrite((void*)&(vertexPair.second), sizeof(int), 1, file);
  }
  
  // write out: vector<pair<int, int> > flags
  int totalFlags = flags.size(); 
  fwrite((void*)&(totalFlags), sizeof(int), 1, file);
  for (unsigned int x = 0; x < flags.size(); x++)
  {
    pair<int, int> flagPair = flags[x];
    fwrite((void*)&(flagPair.first), sizeof(int), 1, file);
    fwrite((void*)&(flagPair.second), sizeof(int), 1, file);
  }

  fclose(file);
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// write out an edge cache
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writeEdgeCacheHuge(const vector<int>& flags, const vector<VEC3I>& indices)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  file = fopen(_cacheFilename.c_str(), "wb");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FAILED TO CACHE OUT EDGES: " << _cacheFilename.c_str() << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    return;
  }

  cout << " Writing out cache file " << _cacheFilename.c_str() << " ... " << flush;
 
  // write out: map<pair<int, int>, bool> _vertexPairs;
  int totalTriplets = _vertexTriplets.size(); 
  fwrite((void*)&(totalTriplets), sizeof(int), 1, file);
  cout << " total triplets: " << totalTriplets << " ";

  map<pair<VEC3I, VEC3I>, bool>::iterator iter;
  for (iter = _vertexTriplets.begin(); iter != _vertexTriplets.end(); iter++)
  {
    pair<VEC3I,VEC3I> vertexTriplet = iter->first;

    VEC3I& first  = vertexTriplet.first;
    VEC3I& second = vertexTriplet.second;
    for (int x = 0; x < 3; x++)
      fwrite((void*)&(first[x]), sizeof(int), 1, file);
    for (int x = 0; x < 3; x++)
      fwrite((void*)&(second[x]), sizeof(int), 1, file);
  }
  
  // write out: vector<pair<int, int> > flags
  int totalFlags = flags.size(); 
  fwrite((void*)&(totalFlags), sizeof(int), 1, file);
  cout << " total flags: " << totalFlags << " ";
  for (unsigned int x = 0; x < flags.size(); x++)
  {
    fwrite((void*)&(flags[x]), sizeof(int), 1, file);
    VEC3I second = indices[x];
    for (int y = 0; y < 3; y++)
      fwrite((void*)&(second[y]), sizeof(int), 1, file);
  }

  fclose(file);
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// write out a partial edge cache, in case computation is 
// interrupted
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writePartialInterpolationCache(const int totalChunks, const int lastCompleted, const int verticesCompleted)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  string partialFilename = _cacheFilename + string(".interp.partial");
  file = fopen(partialFilename.c_str(), "wb");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FAILED TO CACHE OUT INTERPOLATIONS: " << partialFilename.c_str() << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    return;
  }

  cout << " Writing out partial cache file " << partialFilename.c_str() << " ... " << flush;
  fwrite((void*)&totalChunks, sizeof(int), 1, file);
  fwrite((void*)&lastCompleted, sizeof(int), 1, file);
  fwrite((void*)&verticesCompleted, sizeof(int), 1, file);
 
  const int size = _vertices.size();
  //int chunkSize = size / totalChunks;
  cout << verticesCompleted << " vertices of " << size << " ... " << flush;
  for (int x = 0; x < verticesCompleted; x++)
  {
    VEC3F second = _vertices[x];
    for (int y = 0; y < 3; y++)
    {
      // make sure it's always a double
      double value = second[y];
      fwrite((void*)&value, sizeof(double), 1, file);
    }
  }

  fclose(file);
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// write out a partial edge cache, in case computation is 
// interrupted
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writePartialEdgeCache(const int zCurrent, const vector<int>& flags, const vector<VEC3I>& indices)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  string partialFilename = _cacheFilename + string(".partial");
  file = fopen(partialFilename.c_str(), "wb");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FAILED TO CACHE OUT EDGES: " << partialFilename.c_str() << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    return;
  }

  cout << " Writing out partial cache file " << partialFilename.c_str() << " ... " << flush;
  fwrite((void*)&zCurrent, sizeof(int), 1, file);
 
  // write out: map<pair<int, int>, bool> _vertexPairs;
  int totalTriplets = _vertexTriplets.size(); 
  fwrite((void*)&(totalTriplets), sizeof(int), 1, file);
  cout << " total triplets: " << totalTriplets << " ";

  map<pair<VEC3I, VEC3I>, bool>::iterator iter;
  for (iter = _vertexTriplets.begin(); iter != _vertexTriplets.end(); iter++)
  {
    pair<VEC3I,VEC3I> vertexTriplet = iter->first;

    VEC3I& first  = vertexTriplet.first;
    VEC3I& second = vertexTriplet.second;
    for (int x = 0; x < 3; x++)
      fwrite((void*)&(first[x]), sizeof(int), 1, file);
    for (int x = 0; x < 3; x++)
      fwrite((void*)&(second[x]), sizeof(int), 1, file);
  }
  
  // write out: vector<pair<int, int> > flags
  int totalFlags = flags.size(); 
  fwrite((void*)&(totalFlags), sizeof(int), 1, file);
  cout << " total flags: " << totalFlags << " ";
  for (unsigned int x = 0; x < flags.size(); x++)
  {
    int flag = flags[x]; 
    fwrite((void*)&flag, sizeof(int), 1, file);
    VEC3I second = indices[x];
    for (int y = 0; y < 3; y++)
      fwrite((void*)&(second[y]), sizeof(int), 1, file);
  }

  fclose(file);
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// write out an edge cache
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writeEdgeCacheHuge(const vector<pair<int, VEC3I> >& flags)
{
  TIMER functionTimer(__FUNCTION__);
  FILE* file = NULL;
  file = fopen(_cacheFilename.c_str(), "wb");

  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FAILED TO CACHE OUT EDGES: " << _cacheFilename.c_str() << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    return;
  }

  cout << " Writing out cache file " << _cacheFilename.c_str() << " ... " << flush;
 
  // write out: map<pair<int, int>, bool> _vertexPairs;
  int totalTriplets = _vertexTriplets.size(); 
  fwrite((void*)&(totalTriplets), sizeof(int), 1, file);
  cout << " total triplets: " << totalTriplets << " ";

  map<pair<VEC3I, VEC3I>, bool>::iterator iter;
  for (iter = _vertexTriplets.begin(); iter != _vertexTriplets.end(); iter++)
  {
    pair<VEC3I,VEC3I> vertexTriplet = iter->first;

    VEC3I& first  = vertexTriplet.first;
    VEC3I& second = vertexTriplet.second;
    for (int x = 0; x < 3; x++)
      fwrite((void*)&(first[x]), sizeof(int), 1, file);
    for (int x = 0; x < 3; x++)
      fwrite((void*)&(second[x]), sizeof(int), 1, file);
  }
  
  // write out: vector<pair<int, int> > flags
  int totalFlags = flags.size(); 
  fwrite((void*)&(totalFlags), sizeof(int), 1, file);
  cout << " total flags: " << totalFlags << " ";
  for (unsigned int x = 0; x < flags.size(); x++)
  {
    pair<int, VEC3I> flagPair = flags[x];
    fwrite((void*)&(flagPair.first), sizeof(int), 1, file);
    VEC3I second = flagPair.second;
    for (int y = 0; y < 3; y++)
      fwrite((void*)&(second[y]), sizeof(int), 1, file);
  }

  fclose(file);
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// same as computeNonlinearMarchingCubesLowMemoryHuge(), but using CUDA,
// and incrementally write to the OBJ.
//
// This version streams the vertices to files as it goes along, so only
// a single slice and its vertices are ever in memory at the same time.
//
// uses 64-bit ints for counting, in case there are more than 2 billion
// vertices output in the mesh
  
// Writing to a GZipped file as well, because checkpointing is 
// getting really expensive.
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeNonlinearMarchingCubesCudaStreamingCheckpointed64gz(const string& objFilename)
{
  TIMER functionTimer(__FUNCTION__);

  // make sure that long long int actually specifies the width we want
  if (sizeof(int64) != 8)
  {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      std::cout << " long long int is not 64 bits!!!! " << std::endl;
      exit(0);
  }
  cout << " long long int is: " << sizeof(int64) << " bytes, " << sizeof(int64) * 8 << " bits." << std::endl;

  // clear any previous front
  _vertices.clear();
  _triangles.clear();
  _vertexTripletHash.clear();

  // set "outside" to something a lot bigger than the known grid bounds
  _outside = _center[0] + _lengths[0] * 10000;
 
  _xRes = _res[0];
  _yRes = _res[1];
  _zRes = _res[2];
  _slabSize = _xRes * _yRes;

  // hash the existence of edge pairs
  map<pair<VEC3I, VEC3I>, bool> hash0;
  map<pair<VEC3I, VEC3I>, bool> hash1;
  map<pair<VEC3I, VEC3I>, bool>& uniquePairs0 = hash0;
  map<pair<VEC3I, VEC3I>, bool>& uniquePairs1 = hash1;

  // setup the hashes
  _streamingHashCurrent64 = &_streamHash0_64;
  _streamingHashPrevious64 = &_streamHash1_64;

  // track this just for totals later
  int64 triangleVerticesSeen = 0;

  // marching cube triangle edges seen so far
  int64 pairsSeen = 0;

  FIELD_2D& slab0 = _slab0;
  FIELD_2D& slab1 = _slab1;

  gzFile vertexFile;
  gzFile faceFile;
  int zSlab = 0;
  
  // see if a previous checkpointed run occurred
  readStreamingCheckpoint64gz(objFilename, 
                              zSlab,
                              pairsSeen,
                              triangleVerticesSeen,
                              slab0,
                              slab1,
                              uniquePairs0,
                              uniquePairs1,
                              vertexFile, 
                              faceFile);

  if (zSlab == 0)
  {
    computeNonlinearSliceCuda(zSlab);
    readbackCudaSlice(slab1);
    swapCudaBuffers();
  }
  else
    zSlab++;

  // number of nans and infs
  int64 totalNans = 0;
  int64 totalInfs = 0;
  totalNans += slab1.totalNans();
  totalInfs += slab1.totalInfs();

  // record the last time a checkpoint fired
  double lastCheckpointTime = TIMER::totalTimeSeen();

  for (int z = zSlab; z < _zRes - 1; z++)
  {
    // swap in the old "next" slice as the new current ont
    FIELD_2D& old = slab0;
    slab0 = slab1;
    slab1 = old;

    // compute the next needed slice on the fly
    computeNonlinearSliceCuda(z + 1);
    readbackCudaSlice(slab1);

    swapCudaBuffers();

    // swap the hash tables too
    map<pair<VEC3I, VEC3I>, bool>& temp = uniquePairs0;
    uniquePairs0 = uniquePairs1;
    uniquePairs1 = temp;
    totalNans += slab1.totalNans();
    totalInfs += slab1.totalInfs();

    // store all the meaningful flags <flag, index>
    vector<int> flags;
    vector<VEC3I> indices;

    // build out the flags and indices based on the new slab
    for (int y = 0; y < _yRes - 1; y++)
    {
      for (int x = 0; x < _xRes - 1; x++) 
      {
        VEC3I index(x,y,z);

        CUBE cube;
        cube.NNN = slab0(x,y);
        cube.NNP = slab1(x,y);
        cube.NPN = slab0(x,y + 1);
        cube.NPP = slab1(x,y + 1);
        cube.PNN = slab0(x + 1,y);
        cube.PNP = slab1(x + 1,y);
        cube.PPN = slab0(x + 1,y + 1);
        cube.PPP = slab1(x + 1,y + 1);
                
        // construct the flag - can this be an unsigned char?
        int flag =    ((cube.NNN > 0) + 2 *   (cube.NNP > 0) + 4  * (cube.NPN > 0) +
                   8 * (cube.NPP > 0) + 16 *  (cube.PNN > 0) + 32 * (cube.PNP > 0) +
                   64 *(cube.PPN > 0) + 128 * (cube.PPP > 0));

        if (flag == 0 || flag == 255) continue;

        flags.push_back(flag);
        indices.push_back(index);
        
        // three edges are added to _vertexTriplets here
        //
        // flag calls a switch in:
        //      addVertexTriplets(int i, int j, int k, VEC3I index)
        // where the i,j,k are determined by the marching cubes stencil
        switch (flag)
#include "MARCHING_CUBES_VERTICES.include.streaming" 
      }
    }

    // check to see if we've generated this edge already, 
    // both in the current slab and the previous
    uniquePairs1.clear();
    for (unsigned int x = 0; x < _streamingTriplets.size(); x++)
    {
        const pair<VEC3I, VEC3I> triplet = _streamingTriplets[x];
        // make sure we didn't see it in the previous slab
        if (uniquePairs0.find(triplet) != uniquePairs0.end())
            continue;
        uniquePairs1[triplet] = true;
    }

    // add unique edges to the pairs
    vector<pair<VEC3I, VEC3I> > newPairs;
    for (auto iter = uniquePairs1.begin(); iter != uniquePairs1.end(); iter++)
    {
      pair<VEC3I,VEC3I> vertexPair = iter->first;
      newPairs.push_back(vertexPair);
    }

    // clear out the streaming buffer
    _streamingTriplets.clear();

    // do the non-linear interpolations
    vector<VEC3F> newVertices;
    computeSingleSlabEdgeInterpolationsCuda(newPairs, newVertices);

    // write out all the vertices
    gzprintf(vertexFile, "\n# Outputting slab %i\n", z); 
    for (unsigned int i = 0; i < newVertices.size(); i++)
    {
        double vertex[] = {(double)newVertices[i][0], 
                           (double)newVertices[i][1], 
                           (double)newVertices[i][2]};
        gzprintf(vertexFile, "v %.16f %.16f %.16f\n", vertex[0], vertex[1], vertex[2]);
    }

    // hash where each vertex is so the triangle construction looking for it later can find it
    map<pair<VEC3I, VEC3I>, int64>* tempHash = _streamingHashCurrent64;
    _streamingHashCurrent64 = _streamingHashPrevious64;
    _streamingHashPrevious64 = tempHash;
    _streamingHashCurrent64->clear();
    map<pair<VEC3I, VEC3I>, int64>& currentHash = *_streamingHashCurrent64;
    
    for (unsigned int x = 0; x < newPairs.size(); x++)
    {
        pair<VEC3I,VEC3I> vertexTriplet = newPairs[x];
        VEC3I firstIndex = vertexTriplet.first;
        VEC3I secondIndex = vertexTriplet.second;
        currentHash[pair<VEC3I,VEC3I>(firstIndex, secondIndex)] = (int64)x + pairsSeen;
        currentHash[pair<VEC3I,VEC3I>(secondIndex, firstIndex)] = (int64)x + pairsSeen;
    }

    // cache the counter to start on next time
    pairsSeen += (int64)(newPairs.size());

    // go back over the vertex pairs and emit the triangles
    for (unsigned int x = 0; x < flags.size(); x++)
    {
        int flag = flags[x];
        VEC3I index = indices[x];

        // emits using:
        //      addVertexTripletTriangle(int i, int j, int k, VEC3I index)
        // each vertex is emitted using:
        //      VEC3F p0 = computeVertex(i, index);
        // we check whether the vertex was created before using:
        //      int v0 = storeVertex(p0, index);
        // this one does a lookup into _vertexTripletHash
        //
        // then pushes the final result back in:
        //      _triangleVertices64.push_back(v0); 
        switch (flag)
#include "MARCHING_CUBES_TRIANGLES.include.huge.streaming.64"
    }

    // write out the faces
    gzprintf(faceFile, "\n# Outputting slab %i\n", z); 
    for (unsigned int i = 0; i < _triangleVertices64.size() / 3; i++)
    {
      const int64 one = 1;
      int64 indices[] = { _triangleVertices64[3 * i] + one,
                          _triangleVertices64[3 * i + 1] + one,
                          _triangleVertices64[3 * i + 2] + one};

      gzprintf(faceFile, "f %lli %lli %lli\n", indices[0], indices[1], indices[2]);
    }
    triangleVerticesSeen += _triangleVertices64.size();
    _triangleVertices64.clear();

    // has enough time elapsed that we should checkpoint?
    double currentTime = TIMER::totalTimeSeen();
    double timeElapsed = currentTime - lastCheckpointTime;

    // if an hour has elapsed, build a checkpoint
    if (timeElapsed >= _checkpointFrequency)
    {
      cout << " Checkpointing, because at least " << _checkpointFrequency << " seconds have elapsed" << std::endl;
      TIMER::printTimings();
      cout << 100 * ((Real)z / _zRes) << "% COMPLETE" << endl;
      cout << " Currently on slab " << z << " of " << _zRes << endl;
      cout << " Triangles so far: " << triangleVerticesSeen / 3 << endl;
      cout << " infs so far: " << totalInfs << endl;
      cout << " NaNs so far: " << totalNans << endl;

      writeStreamingCheckpoint64gz(objFilename, 
                                   z,
                                   pairsSeen,
                                   triangleVerticesSeen,
                                   slab0,
                                   slab1,
                                   uniquePairs0,
                                   uniquePairs1,
                                   vertexFile, 
                                   faceFile);

      // set the time seen to the exact time now,
      // just in case writing the checkpoint took forever.
      lastCheckpointTime = TIMER::totalTimeSeen();
    }
  }
  cout << " done." << endl;
  cout << " Total infs: " << totalInfs << endl;
  cout << " Total NaNs: " << totalNans << endl;
  cout << " Vertex pairs: " << pairsSeen << endl;
  cout << " Triangle vertices: " << triangleVerticesSeen << endl;
  cout << " Triangles: " << triangleVerticesSeen / 3 << endl;
  gzclose(vertexFile);
  gzclose(faceFile);

  // cat them together
  string vertexFilenameGz = objFilename + string(".vertices.gz");
  string faceFilenameGz = objFilename + string(".faces.gz");
  string unzipVertex = string("gunzip -f ") + vertexFilenameGz;
  cout << " Calling " << unzipVertex.c_str() << endl;
  system(unzipVertex.c_str());

  string unzipFaces = string("gunzip -f ") + faceFilenameGz;
  cout << " Calling " << unzipFaces.c_str() << endl;
  system(unzipFaces.c_str());

  string vertexFilename = objFilename + string(".vertices");
  string faceFilename = objFilename + string(".faces");
  string catCall = string("cat " ) + vertexFilename + string(" ") + faceFilename + string(" > ") + objFilename;
  system(catCall.c_str());
  cout << " Cat'ed together file " << objFilename.c_str() << endl;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::readStreamingCheckpoint64gz(const string& objFilename,
                                                int& zSlab,
                                                int64& pairsSeen,
                                                int64& triangleVerticesSeen,
                                                FIELD_2D& slab0,
                                                FIELD_2D& slab1,
                                                map<pair<VEC3I, VEC3I>, bool>& uniquePairs0,
                                                map<pair<VEC3I, VEC3I>, bool>& uniquePairs1,
                                                gzFile& vertexFile,
                                                gzFile& faceFile)
{
  TIMER functionTimer(__FUNCTION__);

  // reconstruct the file stream filenames
  string vertexFilename = objFilename + string(".vertices.gz");
  string faceFilename = objFilename + string(".faces.gz");

  string vertexCacheFilename = objFilename + string(".vertices.partial.gz");
  string faceCacheFilename = objFilename + string(".faces.partial.gz");
  string indexCacheFilename = objFilename + string(".index.partial");

  // see if a cache file exists (should just be the index of the last computed slab)
  FILE* indexFile = fopen(indexCacheFilename.c_str(), "r");
  bool cacheExists = (indexFile != NULL);

  if (!cacheExists)
  {
    std::cout << " NO cache found." << std::endl;
    vertexFile = gzopen(vertexFilename.c_str(), "w");
    faceFile = gzopen(faceFilename.c_str(), "w");
    zSlab = 0;
    return;
  }

  std::cout << " Cache FOUND." << std::endl;

  // see what the last slab we saw was
  fscanf(indexFile, "%i\n", &zSlab);
  std::cout << " Got to slab " << zSlab << " last time. " << std::endl;

  // read in pairs seen
  fscanf(indexFile, "%lli\n", &pairsSeen);
  
  // read in triangles seen
  fscanf(indexFile, "%lli\n", &triangleVerticesSeen);

  // read in uniquePairs0
  uniquePairs0.clear();
  int pairsSize0;
  fscanf(indexFile, "%i\n", &pairsSize0);
  for (int x = 0; x < pairsSize0; x++)
  {
    int i,j,k;
    fscanf(indexFile, "%i %i %i\n", &i, &j, &k);
    VEC3I first(i,j,k);

    fscanf(indexFile, "%i %i %i\n", &i, &j, &k);
    VEC3I second(i,j,k);

    pair<VEC3I, VEC3I> entry(first, second);
    uniquePairs0[entry] = true;
  }

  // read in uniquePairs1
  uniquePairs1.clear();
  int pairsSize1;
  fscanf(indexFile, "%i\n", &pairsSize1);
  for (int x = 0; x < pairsSize1; x++)
  {
    int i,j,k;
    fscanf(indexFile, "%i %i %i\n", &i, &j, &k);
    VEC3I first(i,j,k);

    fscanf(indexFile, "%i %i %i\n", &i, &j, &k);
    VEC3I second(i,j,k);

    pair<VEC3I, VEC3I> entry(first, second);
    uniquePairs1[entry] = true;
  }

  // read in _streamingHashCurrent
  _streamingHashCurrent64->clear();
  int hashSize0;
  fscanf(indexFile, "%i\n", &hashSize0);
  for (int x = 0; x < hashSize0; x++)
  {
    int i,j,k;
    fscanf(indexFile, "%i %i %i\n", &i, &j, &k);
    VEC3I first(i,j,k);

    fscanf(indexFile, "%i %i %i\n", &i, &j, &k);
    VEC3I second(i,j,k);

    int64 index;
    fscanf(indexFile, "%lli\n", &index);
    pair<VEC3I, VEC3I> entry(first, second);
    (*_streamingHashCurrent64)[entry] = index;
  }

  // read in _streamingHashCurrent
  _streamingHashPrevious64->clear();
  int hashSize1;
  fscanf(indexFile, "%i\n", &hashSize1);
  for (int x = 0; x < hashSize1; x++)
  {
    int i,j,k;
    fscanf(indexFile, "%i %i %i\n", &i, &j, &k);
    VEC3I first(i,j,k);

    fscanf(indexFile, "%i %i %i\n", &i, &j, &k);
    VEC3I second(i,j,k);

    int64 index;
    fscanf(indexFile, "%lli\n", &index);
    pair<VEC3I, VEC3I> entry(first, second);
    (*_streamingHashPrevious64)[entry] = index;
  }
  
  // read in the slabs
  slab0.read(indexFile);
  slab1.read(indexFile);

  fclose(indexFile);

  // copy over the last saved caches
  std::cout << " Copying over previous results ..." << std::endl;
  string vertexCopy = string("cp -f ") + vertexCacheFilename + string(" ") + vertexFilename;
  system(vertexCopy.c_str());
  std::cout << vertexCopy.c_str() << std::endl;
  string faceCopy = string("cp -f ") + faceCacheFilename + string(" ") + faceFilename;
  system(faceCopy.c_str());
  std::cout << faceCopy.c_str() << std::endl;

  vertexFile = gzopen(vertexFilename.c_str(), "a");
  faceFile = gzopen(faceFilename.c_str(), "a");

  // annotate the difference
  gzprintf(vertexFile, "\n# Restart from here, starting at slab %i\n", zSlab); 
  gzprintf(faceFile, "\n# Restart from here, starting at slab %i\n", zSlab); 
  std::cout << " done." << std::endl;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writeStreamingCheckpoint64gz(const string& objFilename,
                                                 const int z,
                                                 const int64 pairsSeen,
                                                 const int64 triangleVerticesSeen,
                                                 const FIELD_2D& slab0,
                                                 const FIELD_2D& slab1,
                                                 const map<pair<VEC3I, VEC3I>, bool>& uniquePairs0,
                                                 const map<pair<VEC3I, VEC3I>, bool>& uniquePairs1,
                                                 gzFile& vertexFile,
                                                 gzFile& faceFile)
{
  TIMER functionTimer(__FUNCTION__);
  
  // reconstruct the file stream filenames
  string vertexFilename = objFilename + string(".vertices.gz");
  string faceFilename = objFilename + string(".faces.gz");
  string vertexCacheFilename = objFilename + string(".vertices.partial.gz");
  string faceCacheFilename = objFilename + string(".faces.partial.gz");

  // close out the file streams
  gzclose(vertexFile);
  gzclose(faceFile);

  // copy them to partial caches
  TIMER copyTimer("Checkpoint Copying");
  string vertexCopy = string("cp -f ") + vertexFilename + string(" ") + vertexCacheFilename;
  string faceCopy = string("cp -f ") + faceFilename + string(" ") + faceCacheFilename;
  std::cout << " Executing " << vertexCopy.c_str() << std::endl;
  system(vertexCopy.c_str());
  std::cout << " Wrote checkpoint file " << vertexCacheFilename.c_str() << std::endl;
  std::cout << " Executing " << faceCopy.c_str() << std::endl;
  system(faceCopy.c_str());
  std::cout << " Wrote checkpoint file " << faceCacheFilename.c_str() << std::endl;
  copyTimer.stop();

  // write out the last z to a separate file
  string indexCacheFilename = objFilename + string(".index.partial");
  FILE* indexFile = fopen(indexCacheFilename.c_str(), "w");
  fprintf(indexFile, "%i\n", z);

  // write out pairs seen
  fprintf(indexFile, "%lli\n", pairsSeen);

  // write out vertices seen
  fprintf(indexFile, "%lli\n", triangleVerticesSeen);

  // write out uniquePairs0
  int pairsSize0 = uniquePairs0.size();
  fprintf(indexFile, "%i\n", pairsSize0);
  for (auto iter = uniquePairs0.begin(); iter != uniquePairs0.end(); iter++)
  {
    VEC3I first = iter->first.first;
    VEC3I second = iter->first.second;
    fprintf(indexFile, "%i %i %i\n", first[0], first[1], first[2]);
    fprintf(indexFile, "%i %i %i\n", second[0], second[1], second[2]);
  }

  // write out uniquePairs1
  int pairsSize1 = uniquePairs1.size();
  fprintf(indexFile, "%i\n", pairsSize1);
  for (auto iter = uniquePairs1.begin(); iter != uniquePairs1.end(); iter++)
  {
    VEC3I first = iter->first.first;
    VEC3I second = iter->first.second;
    fprintf(indexFile, "%i %i %i\n", first[0], first[1], first[2]);
    fprintf(indexFile, "%i %i %i\n", second[0], second[1], second[2]);
  }
  
  // write out _streamingHashCurrent
  int hashSize0 = _streamingHashCurrent64->size();
  fprintf(indexFile, "%i\n", hashSize0);
  for (auto iter = _streamingHashCurrent64->begin(); iter != _streamingHashCurrent64->end(); iter++)
  {
    VEC3I first = iter->first.first;
    VEC3I second = iter->first.second;
    fprintf(indexFile, "%i %i %i\n", first[0], first[1], first[2]);
    fprintf(indexFile, "%i %i %i\n", second[0], second[1], second[2]);
    int64 index = iter->second;
    fprintf(indexFile, "%lli\n", index);
  }

  // write out _streamingHashPrevious
  int hashSize1 = _streamingHashPrevious64->size();
  fprintf(indexFile, "%i\n", hashSize1);
  for (auto iter = _streamingHashPrevious64->begin(); iter != _streamingHashPrevious64->end(); iter++)
  {
    VEC3I first = iter->first.first;
    VEC3I second = iter->first.second;
    fprintf(indexFile, "%i %i %i\n", first[0], first[1], first[2]);
    fprintf(indexFile, "%i %i %i\n", second[0], second[1], second[2]);
    int64 index = iter->second;
    fprintf(indexFile, "%lli\n", index);
  }

  // write out the slabs
  slab0.write(indexFile);
  slab1.write(indexFile);

  fclose(indexFile);

  // reopen the file streams
  vertexFile = gzopen(vertexFilename.c_str(), "a");
  faceFile = gzopen(faceFilename.c_str(), "a");
  std::cout << " Done checkpointing. Okay to interrupt computation right now." << std::endl;
}

//////////////////////////////////////////////////////////////////////
// add vertex triplets to be interpolated
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::accumulateVertexTriplets(int i, int j, int k, VEC3I index)
{
  pair<VEC3I, VEC3I> toAdd;
  toAdd = getVertexTriplets(i, index);
  _streamingTriplets.push_back(toAdd);
  toAdd = getVertexTriplets(j, index);
  _streamingTriplets.push_back(toAdd);
  toAdd = getVertexTriplets(k, index);
  _streamingTriplets.push_back(toAdd);
}

//////////////////////////////////////////////////////////////////////
// add vertex triplets to be interpolated
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addVertexTriplets(int i, int j, int k, VEC3I index)
{
  pair<VEC3I, VEC3I> toAdd;
  toAdd = getVertexTriplets(i, index);
  _vertexTriplets[toAdd] = true;

  toAdd = getVertexTriplets(j, index);
  _vertexTriplets[toAdd] = true;

  toAdd = getVertexTriplets(k, index);
  _vertexTriplets[toAdd] = true;
}

//////////////////////////////////////////////////////////////////////
// add vertex pairs to be interpolated
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addVertexPairs(int i, int j, int k, int index)
{
  pair<int, int> toAdd;
  toAdd = getVertexPair(i, index);
  _vertexPairs[toAdd] = true;

  toAdd = getVertexPair(j, index);
  _vertexPairs[toAdd] = true;

  toAdd = getVertexPair(k, index);
  _vertexPairs[toAdd] = true;
}

//////////////////////////////////////////////////////////////////////
// get the indices for the first and second vertices in the
// interpolation
//////////////////////////////////////////////////////////////////////
pair<int, int> TRIANGLE_MESH::getVertexPair(int i, int index)
{
  pair<int, int> toAdd;
  toAdd.first = index;
  toAdd.second = index;

  switch (i) {
    case 1: // (x,y,z) to (x,y+1,z)
        toAdd.first += 0;
        toAdd.second += _xRes;
        break;
    case 2: // (x,y,z) to (x,y,z+1)
        toAdd.first += 0;
        toAdd.second += _slabSize;
        break;
    case 3: // (x,y,z+1) to (x,y+1,z+1)
        toAdd.first += _slabSize;
        toAdd.second += _xRes + _slabSize;
        break;
    case 4: // (x,y+1,z) to (x,y+1,z+1)
        toAdd.first += _xRes;
        toAdd.second += _xRes + _slabSize;
        break;
    case 5: // (x+1,y,z) to (x+1,y+1,z)
        toAdd.first += 1;
        toAdd.second += 1 + _xRes;
        break;
    case 6: // (x+1,y,z) to (x+1,y,z+1)
        toAdd.first += 1;
        toAdd.second += 1 + _slabSize;
        break;
    case 7: // (x+1,y,z+1) to (x+1,y+1,z+1)
        toAdd.first += 1 + _slabSize;
        toAdd.second += 1 + _xRes + _slabSize;
        break;
    case 8: // (x+1,y+1,z) to (x+1,y+1,z+1)
        toAdd.first += 1 + _xRes;
        toAdd.second += 1 + _xRes + _slabSize;
        break;
    case 9: // (x,y+1,z) to (x+1,y+1,z)
        toAdd.first += _xRes;
        toAdd.second += 1 + _xRes;
        break;
    case 10: // (x,y,z) to (x+1,y,z)
        toAdd.first += 0;
        toAdd.second += 1;
        break;
    case 11: // (x,y,z+1) to (x+1,y,z+1)
        toAdd.first += _slabSize;
        toAdd.second += 1 + _slabSize;
        break;
    case 12: // (x,y+1,z+1) to (x+1,y+1,z+1)
        toAdd.first += _xRes + _slabSize;
        toAdd.second += 1 + _xRes + _slabSize;
        break;
  }

  return toAdd;  
}

//////////////////////////////////////////////////////////////////////
// get the indices for the first and second vertices in the
// interpolation
//////////////////////////////////////////////////////////////////////
pair<VEC3I, VEC3I> TRIANGLE_MESH::getVertexTriplets(int i, VEC3I v)
{
  //pair<vector<int>, vector<int> > toAdd;
  pair<VEC3I, VEC3I> toAdd;
  toAdd.first = v;
  toAdd.second = v;

  // in general: P = (coord+1), N = (coord)
  switch (i) {
        case 1: // (x,y,z) to (x,y+1,z)
        toAdd.second[1] += 1;
        break;
        case 2: // (x,y,z) to (x,y,z+1)
        toAdd.second[2] += 1;
        break;
        case 3: // (x,y,z+1) to (x,y+1,z+1)
        toAdd.first[2] += 1;
        toAdd.second[1] += 1;
        toAdd.second[2] += 1;
        break;
        case 4: // (x,y+1,z) to (x,y+1,z+1)
        toAdd.first[1] += 1;
        toAdd.second[1] += 1;
        toAdd.second[2] += 1;
        break;
        case 5: // (x+1,y,z) to (x+1,y+1,z)
        toAdd.first[0] += 1;
        toAdd.second[0] += 1;
        toAdd.second[1] += 1;
        break;
        case 6: // (x+1,y,z) to (x+1,y,z+1)
        toAdd.first[0] += 1;
        toAdd.second[0] += 1;
        toAdd.second[2] += 1;
        break;
        case 7: // (x+1,y,z+1) to (x+1,y+1,z+1)
        toAdd.first[0] += 1;
        toAdd.first[2] += 1;
        toAdd.second[0] += 1;
        toAdd.second[1] += 1;
        toAdd.second[2] += 1;
        break;
        case 8: // (x+1,y+1,z) to (x+1,y+1,z+1)
        toAdd.first[0] += 1;
        toAdd.first[1] += 1;
        toAdd.second[0] += 1;
        toAdd.second[1] += 1;
        toAdd.second[2] += 1;
        break;
        case 9: // (x,y+1,z) to (x+1,y+1,z)
        toAdd.first[1] += 1;
        toAdd.second[0] += 1;
        toAdd.second[1] += 1;
        break;
        case 10: // (x,y,z) to (x+1,y,z)
        toAdd.second[0] += 1;
        break;
        case 11: // (x,y,z+1) to (x+1,y,z+1)
        toAdd.first[2] += 1;
        toAdd.second[0] += 1;
        toAdd.second[2] += 1;
        break;
        case 12: // (x,y+1,z+1) to (x+1,y+1,z+1)
        toAdd.first[1] += 1;
        toAdd.first[2] += 1;
        toAdd.second[0] += 1;
        toAdd.second[1] += 1;
        toAdd.second[2] += 1;
        break;
  }

  return toAdd;  
}

//////////////////////////////////////////////////////////////////////
// get the (x,y,z) of an index
//////////////////////////////////////////////////////////////////////
VEC3I TRIANGLE_MESH::getXYZ(const int index) const
{
  VEC3I final;
  final[2] = index / _slabSize;
  final[1] = (index % _slabSize) / _xRes;
  final[0] = (index % _slabSize) % _xRes;

  assert(final[0] + final[1] * _xRes + final[2] * _slabSize == index);
  return final;
}

//////////////////////////////////////////////////////////////////////
// get the nonlinear function value here
//////////////////////////////////////////////////////////////////////
Real TRIANGLE_MESH::nonlinearValue(const VEC3F& position, const bool debug) const
{
  QUATERNION iterate(position[0], position[1], position[2], _quaternionSlice);
 
  // TODO: fractal should really be passing this in ...
  const int maxIterations = _maxIterations;
  const Real escape = _escapeRadius;

  // define a quaternion to set things to when we run into floating point
  // problems and need to just bail
  //const QUATERNION bail(10.0 * escape, 10.0 * escape, 10.0 * escape, 10.0 * escape);
  const Real bailValue = REAL_MAX;
  QUATERNION bail(bailValue, bailValue, bailValue, bailValue);

  Real magnitude = iterate.magnitude();
  int totalIterations = 0;
  bool bailed = false;
  while (magnitude < escape && totalIterations < maxIterations)
  {
    QUATERNION topEval = _top.evaluateScaledPowerFactored(iterate);
    QUATERNION bottomEval;
    
    if (_bottom.totalRoots() > 0)
    {
      bottomEval = _bottom.evaluateScaledPowerFactored(iterate);
      
      // explicitly perform the division here
      const Real magnitudeSq = bottomEval.magnitudeSq();
      const Real magnitudeSqInv = 1.0 / magnitudeSq;
      const QUATERNION bottomInv = bottomEval.conjugate() / magnitudeSq;

      // if the division is tiny, bail
      if (magnitudeSqInv > 1e-200) 
        iterate = topEval * bottomInv;
      else
      {
        iterate = bail;
        bailed = true;
      }
    }
    else
      iterate = topEval;

    if (!bailed)
      iterate *= _expScaling;
    
    if (!bailed)
      magnitude = iterate.magnitude();
    else
      magnitude = REAL_MAX;

    totalIterations++;

    // see if it fell into the black hole at the origin
    if (magnitude < 10.0 * REAL_MIN)
      totalIterations = maxIterations;
  }
  return log(magnitude) - _isosurface;
}

//////////////////////////////////////////////////////////////////////
// initialize CUDA for the edge interpolations
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::mallocEdgeCuda(const vector<pair<VEC3I, VEC3I> >& pairs)
{
  TIMER functionTimer(__FUNCTION__);

  // allocate the vertex fields
  const int totalPairs = pairs.size();
  const size_t widthInt3 = sizeof(int3)*totalPairs;
  const size_t widthDouble3 = sizeof(double3)*totalPairs;
  checkCudaErrors(cudaMalloc((void**)&_firstVerticesCuda, widthInt3));
  checkCudaErrors(cudaMalloc((void**)&_secondVerticesCuda, widthInt3));
  checkCudaErrors(cudaMalloc((void**)&_finalVerticesCuda, widthDouble3));

  // copy over the vertex data
  int3* firstVertices = new int3[totalPairs];
  int3* secondVertices = new int3[totalPairs];

  for (int x = 0; x < totalPairs; x++)
  {
    firstVertices[x].x = pairs[x].first[0];
    firstVertices[x].y = pairs[x].first[1];
    firstVertices[x].z = pairs[x].first[2];
    
    secondVertices[x].x = pairs[x].second[0];
    secondVertices[x].y = pairs[x].second[1];
    secondVertices[x].z = pairs[x].second[2];
  }

  checkCudaErrors(cudaMemcpy(_firstVerticesCuda, firstVertices, widthInt3, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(_secondVerticesCuda, secondVertices, widthInt3, cudaMemcpyHostToDevice));

  delete[] firstVertices;
  delete[] secondVertices;
}

//////////////////////////////////////////////////////////////////////
// free all the CUDA memory for edge interpolations
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::freeEdgeCuda()
{
  TIMER functionTimer(__FUNCTION__);
  checkCudaErrors(cudaFree(_firstVerticesCuda));
  checkCudaErrors(cudaFree(_secondVerticesCuda));
  checkCudaErrors(cudaFree(_finalVerticesCuda));
}

//////////////////////////////////////////////////////////////////////
// do a midpoint search
//////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::midpointSearchForLoop(const VEC3F& positiveVertex, const Real& positiveValue, 
                                           const VEC3F& negativeVertex, const Real& negativeValue) const
{
  VEC3F pVert = positiveVertex;
  VEC3F nVert = negativeVertex;
  Real pVal = positiveValue;
  Real nVal = negativeValue;

  VEC3F midpointVertex;
  Real midpointValue;

  // 6 iterations resolves roughly millimeter level details if a grid edge is a meter
  for (int x = 0; x <= 6; x++)
  {
      midpointVertex = (pVert + nVert) * 0.5;

      if (pVal * nVal >= 0.0)
      {
        // if this gets tripped, and you're using floats or doubles,
        // you may need to bump up the precision of your representation
        cout << " positive: " << pVal << " negative: " << nVal << endl;
        cout << " p vertex: " << pVert << " n vertex: " << nVert << endl;
        cout << " recursion: " << x << endl;
        exit(0);
        return midpointVertex;
      }

      // this resolves roughly millimeter level details if a grid edge is a meter
      //if (x == 6) return midpointVertex;

      midpointValue = nonlinearValue(midpointVertex);
      if (fabs(midpointValue) < 1e-8) return midpointVertex;

      if (midpointValue < 0)
      {
        nVert = midpointVertex;
        nVal = midpointValue;
      }
      else
      {
        pVert = midpointVertex;
        pVal = midpointValue;
      }
  }
  return midpointVertex;
}

//////////////////////////////////////////////////////////////////////
// do a midpoint search
//////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::midpointSearch(const VEC3F& positiveVertex, const Real& positiveValue, 
                                    const VEC3F& negativeVertex, const Real& negativeValue, 
                                    const int recursion) const
{
  VEC3F midpointVertex = positiveVertex + negativeVertex;
  midpointVertex *= 0.5;
  if (positiveValue * negativeValue >= 0.0)
  {
    // if this gets tripped, and you're using floats or doubles,
    // you may need to bump up the precision of your representation
#pragma omp critical
    {
      cout << " positive: " << positiveValue << " negative: " << negativeValue << endl;
      cout << " p vertex: " << positiveVertex << " n vertex: " << negativeVertex << endl;
      cout << " recursion: " << recursion << endl;
      exit(0);
    }

    return midpointVertex;
  }

  // this resolves roughly millimeter level details if a grid edge is a meter
  if (recursion >= 6)
  {
    return midpointVertex;
  }

  Real midpointValue = nonlinearValue(midpointVertex);
  if (fabs(midpointValue) < 1e-8)
  {
    return midpointVertex;
  }

  if (midpointValue < 0)
    return midpointSearch(positiveVertex, positiveValue,
                          midpointVertex, midpointValue, recursion + 1);

  return midpointSearch(midpointVertex, midpointValue,
                        negativeVertex, negativeValue, recursion + 1);
}

//////////////////////////////////////////////////////////////////////
// add a triangle to the list
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addStagedTriangle(int i, int j, int k, int index, const CUBE& cube, const VEC3F& center)
{
  VEC3F p0 = computeStagedVertex(i, index, cube, center);
  VEC3F p1 = computeStagedVertex(j, index, cube, center);
  VEC3F p2 = computeStagedVertex(k, index, cube, center);
  if (p0[0] == _outside || p1[0] == _outside || p2[0] == _outside) return;
 
  // if the vertex has been computed before, don't duplicate it 
  int v0 = storeVertex(p0, index);
  int v1 = storeVertex(p1, index);
  int v2 = storeVertex(p2, index);

  _triangleVertices.push_back(v0); 
  _triangleVertices.push_back(v1); 
  _triangleVertices.push_back(v2); 
}

//////////////////////////////////////////////////////////////////////
// add a triangle to the list
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addTriangle(int i, int j, int k, int index)
{
  VEC3F p0 = computeVertex(i, index);
  VEC3F p1 = computeVertex(j, index);
  VEC3F p2 = computeVertex(k, index);
  if (p0[0] == _outside || p1[0] == _outside || p2[0] == _outside) return;
 
  // if the vertex has been computed before, don't duplicate it 
  int v0 = storeVertex(p0, index);
  int v1 = storeVertex(p1, index);
  int v2 = storeVertex(p2, index);

  _triangleVertices.push_back(v0); 
  _triangleVertices.push_back(v1); 
  _triangleVertices.push_back(v2); 
}

//////////////////////////////////////////////////////////////////////
// add a triangle to the list whose vertices were all precomputed
// by computeEdgeInterpolations()
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addStreamingTriangle64(int i, int j, int k, VEC3I index)
{
  pair<VEC3I,VEC3I> vp0 = getVertexTriplets(i, index);
  pair<VEC3I,VEC3I> vp1 = getVertexTriplets(j, index);
  pair<VEC3I,VEC3I> vp2 = getVertexTriplets(k, index);

  int64 v0 = -1;
  int64 v1 = -1;
  int64 v2 = -1;

  map<pair<VEC3I, VEC3I>, int64>& current = *_streamingHashCurrent64;
  map<pair<VEC3I, VEC3I>, int64>& previous = *_streamingHashPrevious64;
  if (current.find(vp0) != current.end())
      v0 = current[vp0];
  else if (previous.find(vp0) != previous.end())
      v0 = previous[vp0];
  
  if (current.find(vp1) != current.end())
      v1 = current[vp1];
  else if (previous.find(vp1) != previous.end())
      v1 = previous[vp1];

  if (current.find(vp2) != current.end())
      v2 = current[vp2];
  else if (previous.find(vp2) != previous.end())
      v2 = previous[vp2];

  assert(v0 != -1);
  assert(v1 != -1);
  assert(v2 != -1);

  _triangleVertices64.push_back(v0); 
  _triangleVertices64.push_back(v1); 
  _triangleVertices64.push_back(v2); 
}

//////////////////////////////////////////////////////////////////////
// add a triangle to the list whose vertices were all precomputed
// by computeEdgeInterpolations()
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addVertexTripletTriangle(int i, int j, int k, VEC3I index)
{
  pair<VEC3I,VEC3I> vp0 = getVertexTriplets(i, index);
  pair<VEC3I,VEC3I> vp1 = getVertexTriplets(j, index);
  pair<VEC3I,VEC3I> vp2 = getVertexTriplets(k, index);

  int v0 = _vertexTripletHash[vp0];
  int v1 = _vertexTripletHash[vp1];
  int v2 = _vertexTripletHash[vp2];

  _triangleVertices.push_back(v0); 
  _triangleVertices.push_back(v1); 
  _triangleVertices.push_back(v2); 
}

//////////////////////////////////////////////////////////////////////
// add a triangle to the list whose vertices were all precomputed
// by computeEdgeInterpolations()
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addVertexPairTriangle(int i, int j, int k, int index)
{
  pair<int,int> vp0 = getVertexPair(i, index);
  pair<int,int> vp1 = getVertexPair(j, index);
  pair<int,int> vp2 = getVertexPair(k, index);

  int v0 = _vertexPairHash[vp0];
  int v1 = _vertexPairHash[vp1];
  int v2 = _vertexPairHash[vp2];

  _triangleVertices.push_back(v0); 
  _triangleVertices.push_back(v1); 
  _triangleVertices.push_back(v2); 
}

//////////////////////////////////////////////////////////////////////
// get the edge point
//////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::computeStagedVertex(int i, int index, const CUBE& cube, const VEC3F& center) 
{
  int dummy = index;
  dummy++;
  VEC3F point(_outside);
  Real dist[2];

  // the base is the current cell center
  const VEC3F& base = center;

  // in general: P = (coord+1), N = (coord)
  switch (i) {
        case 1: // (x,y,z) to (x,y+1,z)
                point[0] = base[0];
                dist[0] = cube.NNN;
                dist[1] = cube.NPN;
                point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
                point[2] = base[2];
                break;
        case 2: // (x,y,z) to (x,y,z+1)
                point[0] = base[0];
                point[1] = base[1];
                dist[0] = cube.NNN;
                dist[1] = cube.NNP;
                point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
                break;
        case 3: // (x,y,z+1) to (x,y+1,z+1)
                point[0] = base[0];
                dist[0] = cube.NNP;
                dist[1] = cube.NPP;
                point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
                point[2] = base[2] + _fieldDeltas[2];
                break;
        case 4: // (x,y+1,z) to (x,y+1,z+1)
                point[0] = base[0];
                point[1] = base[1] + _fieldDeltas[1];
                dist[0] = cube.NPN;
                dist[1] = cube.NPP;
                point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
                break;
        case 5: // (x+1,y,z) to (x+1,y+1,z)
                point[0] = base[0] + _fieldDeltas[0];
                dist[0] = cube.PNN;
                dist[1] = cube.PPN;
                point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
                point[2] = base[2];
                break;
        case 6: // (x+1,y,z) to (x+1,y,z+1)
                point[0] = base[0] + _fieldDeltas[0];
                point[1] = base[1];
                dist[0] = cube.PNN;
                dist[1] = cube.PNP;
                point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
                break;
        case 7: // (x+1,y,z+1) to (x+1,y+1,z+1)
                point[0] = base[0] + _fieldDeltas[0];
                dist[0] = cube.PNP;
                dist[1] = cube.PPP;
                point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
                point[2] = base[2] + _fieldDeltas[2];
                break;
        case 8: // (x+1,y+1,z) to (x+1,y+1,z+1)
                point[0] = base[0] + _fieldDeltas[0];
                point[1] = base[1] + _fieldDeltas[1];
                dist[0] = cube.PPN;
                dist[1] = cube.PPP;
                point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
                break;
        case 9: // (x,y+1,z) to (x+1,y+1,z)
                dist[0] = cube.NPN;
                dist[1] = cube.PPN;
                point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
                point[1] = base[1] + _fieldDeltas[1];
                point[2] = base[2];
                break;
        case 10: // (x,y,z) to (x+1,y,z)
                dist[0] = cube.NNN;
                dist[1] = cube.PNN;
                point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
                point[1] = base[1];
                point[2] = base[2];
                break;
        case 11: // (x,y,z+1) to (x+1,y,z+1)
                dist[0] = cube.NNP;
                dist[1] = cube.PNP;
                point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
                point[1] = base[1];
                point[2] = base[2] + _fieldDeltas[2];
                break;
        case 12: // (x,y+1,z+1) to (x+1,y+1,z+1)
                dist[0] = cube.NPP;
                dist[1] = cube.PPP;
                point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
                point[1] = base[1] + _fieldDeltas[1];
                point[2] = base[2] + _fieldDeltas[2];
                break;
        }
  return point;
}

//////////////////////////////////////////////////////////////////////
// get the edge point
//////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::computeVertex(int i, int index) 
{
  int dummy = index;
  dummy++;
  VEC3F point(_outside);
  Real dist[2];

  // the base is the current cell center
  const VEC3F& base = _cellCenter;

        switch (i) {
        case 1:
                point[0] = base[0];
                dist[0] = _NNN;
                dist[1] = _NPN;
                point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
                point[2] = base[2];
                break;
        case 2:
                point[0] = base[0];
                point[1] = base[1];
                dist[0] = _NNN;
                dist[1] = _NNP;
                point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
                break;
        case 3:
                point[0] = base[0];
                dist[0] = _NNP;
                dist[1] = _NPP;
                point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
                point[2] = base[2] + _fieldDeltas[2];
                break;
        case 4:
                point[0] = base[0];
                point[1] = base[1] + _fieldDeltas[1];
                dist[0] = _NPN;
                dist[1] = _NPP;
                point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
                break;
        case 5:
                point[0] = base[0] + _fieldDeltas[0];
                dist[0] = _PNN;
                dist[1] = _PPN;
                point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
                point[2] = base[2];
                break;
        case 6:
                point[0] = base[0] + _fieldDeltas[0];
                point[1] = base[1];
                dist[0] = _PNN;
                dist[1] = _PNP;
                point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
                break;
        case 7:
                point[0] = base[0] + _fieldDeltas[0];
                dist[0] = _PNP;
                dist[1] = _PPP;
                point[1] = base[1] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[1];
                point[2] = base[2] + _fieldDeltas[2];
                break;
        case 8:
                point[0] = base[0] + _fieldDeltas[0];
                point[1] = base[1] + _fieldDeltas[1];
                dist[0] = _PPN;
                dist[1] = _PPP;
                point[2] = base[2] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[2];
                break;
        case 9:
                dist[0] = _NPN;
                dist[1] = _PPN;
                point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
                point[1] = base[1] + _fieldDeltas[1];
                point[2] = base[2];
                break;
        case 10:
                dist[0] = _NNN;
                dist[1] = _PNN;
                point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
                point[1] = base[1];
                point[2] = base[2];
                break;
        case 11:
                dist[0] = _NNP;
                dist[1] = _PNP;
                point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
                point[1] = base[1];
                point[2] = base[2] + _fieldDeltas[2];
                break;
        case 12:
                dist[0] = _NPP;
                dist[1] = _PPP;
                point[0] = base[0] + dist[0] / (dist[0] - dist[1]) * _fieldDeltas[0];
                point[1] = base[1] + _fieldDeltas[1];
                point[2] = base[2] + _fieldDeltas[2];
                break;
        }
  return point;
}

//////////////////////////////////////////////////////////////////////
// see if the vertex has been computed before, and if not, store it
//////////////////////////////////////////////////////////////////////
int TRIANGLE_MESH::storeVertex(VEC3F& vertex, int index)
{
  // get vertices that have been created near this cell before
  vector<int>& nearest = _vertexHash[index];

  // check them all to see if any are near
  for (unsigned int x = 0; x < nearest.size(); x++)
  {
    int nearestIndex = nearest[x];
    VEC3F diff = _vertices[nearestIndex] - vertex;
    Real magnitude = norm(diff);

    // if they're really close to each other, just return that.
    if (magnitude < 1e-6)
      return nearestIndex;
  }

  // go ahead and add it
  _vertices.push_back(vertex);
 
  // add it to the hash table as well 
  int vectorIndex = _vertices.size() - 1;

  int fatten = 1;
  for (int z = -fatten; z <= fatten; z++)
    for (int y = -fatten; y <= fatten; y++)
      for (int x = -fatten; x <= fatten; x++)
        _vertexHash[index + x + _xRes * y + _slabSize * z].push_back(vectorIndex);

  return vectorIndex;
}

//////////////////////////////////////////////////////////////////////
//get the range and center of the mesh
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::getBounds(VEC3F& maxVert, VEC3F& minVert, Real& maxLength)
{
  maxVert = minVert = _vertices[0];
  // get bounds
  for (unsigned int x = 0; x < _vertices.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      if (_vertices[x][y] < minVert[y]) minVert[y] = _vertices[x][y];
      if (_vertices[x][y] > maxVert[y]) maxVert[y] = _vertices[x][y];
    }
  VEC3F lengths = maxVert - minVert;
  maxLength = lengths.maxElement();
}

//////////////////////////////////////////////////////////////////////
// Normalize mesh to 1x1x1 cube centered at (0,0,0)
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::normalize()
{
  VEC3F maxVert(_vertices[0]);
  VEC3F minVert(_vertices[0]);

  // get bounds
  for (unsigned int x = 0; x < _vertices.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      if (_vertices[x][y] < minVert[y]) minVert[y] = _vertices[x][y];
      if (_vertices[x][y] > maxVert[y]) maxVert[y] = _vertices[x][y];
    }

  VEC3F half(0.5, 0.5, 0.5);
  VEC3F diff = maxVert - minVert;
  double maxDiff = diff.maxElement();
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x] -= minVert;
    _vertices[x] *= 1.0 / maxDiff;
    _vertices[x] -= half;
  }
}

//////////////////////////////////////////////////////////////////////
// fit the geometry inside a 1x1x1 box
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::normalize(Real padding)
{
  VEC3F currentMin(_vertices[0]);
  VEC3F currentMax(_vertices[0]);

  for (unsigned int i = 0; i < _vertices.size (); i++)
  {
    for (int j = 0; j < 3; j++)
    {
      currentMin[j] = min(currentMin[j],_vertices[i][j]);
      currentMax[j] = max(currentMax[j],_vertices[i][j]);
    }
  }
  cout << " Min: " << currentMin << endl;
  cout << " Max: " << currentMax << endl;

  VEC3F recenter = (currentMin + currentMax) * (Real)0.5f;

  // translate everything to the bounding box center
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] -= recenter;

  // find the maximum magnitude
  double maxVal = 0.0f;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    maxVal = (fabs(_vertices[x][0]) > maxVal) ? fabs(_vertices[x][0]) : maxVal;
    maxVal = (fabs(_vertices[x][1]) > maxVal) ? fabs(_vertices[x][1]) : maxVal;
    maxVal = (fabs(_vertices[x][2]) > maxVal) ? fabs(_vertices[x][2]) : maxVal;
  }
  cout << " Max value found: " << maxVal << endl;

  // add a little padding for the last grid cell
  Real scale = 0.5 - padding;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] *= scale / maxVal;

  // translate everything to 0.5, 0.5, 0.5
  VEC3F half(0.5, 0.5, 0.5);
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] += half;
}

///////////////////////////////////////////////////////////////////////
// center using vertex means
///////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::vertexMean() const
{
  VEC3F sum;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    sum += _vertices[x];

  return sum * (Real)(1.0 / _vertices.size());
}

///////////////////////////////////////////////////////////////////////
// return a bounding box
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  mins = _vertices[0];
  maxs = _vertices[0];

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    if (_vertices[x][0] < mins[0]) mins[0] = _vertices[x][0];
    if (_vertices[x][1] < mins[1]) mins[1] = _vertices[x][1];
    if (_vertices[x][2] < mins[2]) mins[2] = _vertices[x][2];

    if (_vertices[x][0] > maxs[0]) maxs[0] = _vertices[x][0];
    if (_vertices[x][1] > maxs[1]) maxs[1] = _vertices[x][1];
    if (_vertices[x][2] > maxs[2]) maxs[2] = _vertices[x][2];
  }
}

////////////////////////////////////////////////////////////////////////////
// write the raw data to a file stream
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::write(FILE* file)
{
  // write out dims
  int totalVertices = _vertices.size();
  int totalTriangles = _triangles.size();
  fwrite((void*)&totalVertices, sizeof(int), 1, file);
  fwrite((void*)&totalTriangles, sizeof(int), 1, file);

  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x].write(file);
  
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    int indices[3];
    indices[0] = _vertexIndices[_triangles[x].vertex(0)];
    indices[1] = _vertexIndices[_triangles[x].vertex(1)];
    indices[2] = _vertexIndices[_triangles[x].vertex(2)];

    fwrite((void*)&indices[0], sizeof(int), 3, file);
  }
}

////////////////////////////////////////////////////////////////////////////
// read in the raw data from a file stream
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::read(FILE* file)
{
  // read in dims
  int totalVertices;
  int totalTriangles;
  fread((void*)&totalVertices, sizeof(int), 1, file);
  fread((void*)&totalTriangles, sizeof(int), 1, file);

  _vertices.resize(totalVertices);
  for (int x = 0; x < totalVertices; x++)
    _vertices[x].read(file);
 
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

  _triangles.resize(totalTriangles); 
  for (int x = 0; x < totalTriangles; x++)
  {
    int indices[3];
    fread((void*)&indices[0], sizeof(int), 3, file);

    _triangles[x].vertex(0) = &_vertices[indices[0]];
    _triangles[x].vertex(1) = &_vertices[indices[1]];
    _triangles[x].vertex(2) = &_vertices[indices[2]];
  }
}

////////////////////////////////////////////////////////////////////////////
// compute some coarse normals for each vertex
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeNormals()
{
  _normals.clear();
  _normals.resize(_vertices.size());

  // instead of an area weighting, we'll store each normal un-normalized,
  // since the cross product magnitude will reflect the area of the triangle
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    TRIANGLE& triangle = _triangles[x];
    vector<VEC3F*>& vertices = triangle.vertices();
    VEC3F normal = cross(*vertices[1] - *vertices[0], 
                         *vertices[2] - *vertices[0]);

    int indices[] = {_vertexIndices[vertices[0]], _vertexIndices[vertices[1]], 
                     _vertexIndices[vertices[2]]};
    _normals[indices[0]] += normal;
    _normals[indices[1]] += normal;
    _normals[indices[2]] += normal;
  }
  
  // go through and normalize after everything
  for (unsigned int x = 0; x < _normals.size(); x++)
    _normals[x].normalize(); 
}

////////////////////////////////////////////////////////////////////////////
// emulate the cell center lookup for FIELD_3D --
// assumes that _xRes, _yRes, _zRes,
//              _center, _lengths, and _dx are populated
////////////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::cellCenter(int x, int y, int z)
{
  VEC3F halfLengths = (Real)0.5 * _lengths;

  // set it to the lower corner
  VEC3F final = _center - halfLengths;

  // displace to the NNN corner
  final[0] += x * _dxs[0];
  final[1] += y * _dxs[1];
  final[2] += z * _dxs[2];

  // displace it to the cell center
  final[0] += _dxs[0] * 0.5;
  final[1] += _dxs[1] * 0.5;
  final[2] += _dxs[2] * 0.5;

  return final;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSingleSlabEdgeInterpolationsCuda(const vector<pair<VEC3I, VEC3I> >& newPairs,
                                                            vector<VEC3F>& newVertices)
{
  TIMER functionTimer(__FUNCTION__);

  // allocate the space needed for the new pairs
  mallocEdgeCuda(newPairs);

  // do the non-linear interpolations on the GPU
  const int totalPairs = newPairs.size();
  computeEdgesOnGPU(_firstVerticesCuda, _secondVerticesCuda, totalPairs, _finalVerticesCuda);

  // read back the final vertices
  const size_t widthDouble3 = sizeof(double3) * totalPairs;
  double3* finalVertices = new double3[totalPairs];
  checkCudaErrors(cudaMemcpy(finalVertices, _finalVerticesCuda, widthDouble3, cudaMemcpyDeviceToHost));

  // copy them into the palatable form
  newVertices.resize(totalPairs);
  for (int x = 0; x < totalPairs; x++)
  {
    newVertices[x][0] = finalVertices[x].x;
    newVertices[x][1] = finalVertices[x].y;
    newVertices[x][2] = finalVertices[x].z;
  }

  // reclaim the GPU memory
  freeEdgeCuda();
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSingleSlabEdgeInterpolations(const vector<pair<VEC3I, VEC3I> >& newPairs,
                                                        vector<VEC3F>& newVertices)
{
    TIMER functionTimer(__FUNCTION__);

    // do the non-linear interpolations
    newVertices.resize(newPairs.size());

#pragma omp parallel
#pragma omp for  schedule(dynamic)
    for (unsigned int x = 0; x < newPairs.size(); x++)
    {
      const pair<VEC3I,VEC3I>& vertexPair = newPairs[x];

      const VEC3I& firstXYZ = vertexPair.first;
      VEC3F firstVertex = cellCenter(firstXYZ[0], firstXYZ[1], firstXYZ[2]);
      Real firstValue = nonlinearValue(firstVertex, false);
        
      const VEC3I& secondXYZ = vertexPair.second;
      VEC3F secondVertex = cellCenter(secondXYZ[0], secondXYZ[1], secondXYZ[2]);
      Real secondValue = nonlinearValue(secondVertex, false);

      VEC3F positiveVertex = firstVertex;
      Real positiveValue = firstValue;
      VEC3F negativeVertex = secondVertex;
      Real negativeValue = secondValue;

      if (positiveValue * negativeValue >= 0.0)
      {
      #pragma omp critical
          {
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout << " field dims: " << endl;
            cout << " res: " << _res << endl;
            cout << " lengths: " << _lengths << endl;
            cout << " center:  " << _center << endl;
            cout << " dxs:     " << _dxs << endl;
            cout << " p vertex: " << positiveVertex << " n vertex: " << negativeVertex << endl;
            cout << " first XYZ: " << firstXYZ << " second XYZ:" << secondXYZ << endl;

            cout << " positive:           " << positiveValue                                << " negative:           " << negativeValue << endl;

            exit(0);
          }
      }

      if (firstValue < 0)
      {
        positiveVertex = secondVertex;
        positiveValue = secondValue;
        negativeVertex = firstVertex;
        negativeValue = firstValue;
      }

      // this turns the midpoint search on and off. If you want to compare to just traditional
      // marching cubes with linear interpolation, set this to 0.
      //VEC3F finalVertex = midpointSearch(positiveVertex, positiveValue, negativeVertex, negativeValue);
      VEC3F finalVertex = midpointSearchForLoop(positiveVertex, positiveValue, negativeVertex, negativeValue);
      newVertices[x] = finalVertex;
    }
}
