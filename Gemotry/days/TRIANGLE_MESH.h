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
//////////////////////////////////////////////////////////////////////
// Triangle Mesh class
//////////////////////////////////////////////////////////////////////

#ifndef TRIANGLE_MESH_H
#define TRIANGLE_MESH_H

#include <map>
#include <unordered_map>
#include <SETTINGS.h>
#include <TRIANGLE.h>
#include <FIELD_2D.h>
#include <FIELD_3D.h>
#include <POLYNOMIAL_4D.h>
#include <cuda_runtime.h>

typedef long long int int64;

using namespace std;
class TRIANGLE_MESH
{
public:
  TRIANGLE_MESH();

  // do a non-linear marching cubes on just two slabs at a shot,
  // write out the vertices first, delete them, and then output the faces
  //
  // This one can currently handle the largest meshes
  TRIANGLE_MESH(const VEC3F& center, 
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
                const double checkpointFrequency);

  // destructor
  virtual ~TRIANGLE_MESH();

  // Normalize mesh to 1x1x1 cube centered at (0,0,0)
  void normalize();
  void normalize(Real padding);

  //get the range of the original mesh
  void getBounds(VEC3F& maxVert, VEC3F& minVert, Real& maxLength);

  const vector<TRIANGLE>& triangles() const { return _triangles; };
  const vector<VEC3F>& vertices() const { return _vertices; };
  vector<VEC3F>& vertices() { return _vertices; };
  vector<VEC3F>& texcoords() { return _texcoords; };
  int vertexIndex(int triangle, int vertex) { return _vertexIndices[_triangles[triangle].vertex(vertex)]; };

  // center using vertex means
  VEC3F vertexMean() const;

  // return a bounding box
  void boundingBox(VEC3F& mins, VEC3F& maxs);

  // compute marching cubes in stages, with non-linear solve in between
  void computeNonlinearMarchingCubes(const FIELD_3D& field, const bool verbose = false);
  
  ////////////////////////////////////////////////////////////////////////////
  // read/write support
  ////////////////////////////////////////////////////////////////////////////
  int sizeTriangles(){return _triangles.size();}
  int sizeVertices(){return _vertices.size();}

  // write the raw data to a file stream
  void write(FILE* file);
  
  // read from a raw file stream
  void read(FILE* file);

private:
  vector<VEC3F> _vertices;
  vector<VEC3F> _normals;
  vector<VEC3F> _texcoords;
  vector<TRIANGLE> _triangles;

  // hash table allowing lookup of vertex index from address
  map<VEC3F*, int> _vertexIndices;

  // distance field dims from marching cubes
  int _xRes;
  int _yRes;
  int _zRes;
  int _slabSize;

  // hash table that allows quick lookup of if a vertex has been created
  // by marching cubes previously
  map<int, vector<int> > _vertexHash;

  ////////////////////////////////////////////////////////////////////////////
  // Marching Cubes support
  ////////////////////////////////////////////////////////////////////////////

  // cached marching cube distances
  Real _NNN, _NNP, _NPN, _NPP, _PNN, _PNP, _PPN, _PPP; 

  // cached marching cube deltas
  VEC3F _fieldDeltas;

  // cell center of the current grid cell being marched
  VEC3F _cellCenter;

  // what is considered "outside" by grid being marched
  Real _outside;

  // store index triplets initially for the triangles, and then set everything
  // to pointers once the vector is done being resized
  vector<int> _triangleVertices;
  
  // in case we're doing the 64 bit version
  vector<int64> _triangleVertices64;

  // struct needed to pass data during staged marching cubes
  struct CUBE {
    union {
      struct { float NNN, NNP, NPN, NPP, PNN, PNP, PPN, PPP; };
      float data[8];
    };
  };

  // vertex pairs to compute an interpolation point between
  map<pair<int, int>, bool> _vertexPairs;
  map<pair<VEC3I, VEC3I>, bool> _vertexTriplets;
  vector<pair<VEC3I, VEC3I> > _streamingTriplets;

  // where in _vertices is the vertex that corresponds to this pair?
  map<pair<int, int>, int> _vertexPairHash;
  map<pair<VEC3I, VEC3I >, int> _vertexTripletHash;

  // pointers to the hashes for the current and previous slabs
  map<pair<VEC3I, VEC3I>, int64> _streamHash0_64;
  map<pair<VEC3I, VEC3I>, int64> _streamHash1_64;
  map<pair<VEC3I, VEC3I>, int64>* _streamingHashCurrent64;
  map<pair<VEC3I, VEC3I>, int64>* _streamingHashPrevious64;

  // the field being marching cubed
  const FIELD_3D* _toMarch;

  // polynomials for non-linear marching cubes
  POLYNOMIAL_4D _top;
  POLYNOMIAL_4D _bottom;
  Real _expScaling;
  int _maxIterations;
  Real _quaternionSlice;
  Real _isosurface;

  Real _escapeRadius;

  // path to cache intermediate results to
  string _cacheFilename;

  // two-slab non-linear marching cubes vars
  const VEC3I _res;
  const VEC3F _lengths;
  const VEC3F _center;
  const VEC3F _dxs;
  FIELD_2D _slab0;
  FIELD_2D _slab1;

  // rotation from the original field
  const QUATERNION _fieldRotation;

  // CUDA fields
  double* _fieldCurrentCuda;
  double* _fieldOldCuda;
  unsigned char* _flagsCuda;
  int* _occupanciesCuda;
  int* _compactedFlagsCuda;
  int* _indexTranslationCuda;
  size_t _fieldCurrentPitch;
  size_t _fieldOldPitch;

  int3* _firstVerticesCuda;
  int3* _secondVerticesCuda;
  double3* _finalVerticesCuda;

  // pinned host memory
  unsigned char* _flagsPinned;

  // how often should we checkpoint?
  double _checkpointFrequency;

  // add a triangle to the list
  void addTriangle(int i, int j, int k, int index);
  void addStagedTriangle(int i, int j, int k, int index, const CUBE& cube, const VEC3F& center);

  // add a triangle to the list whose vertices were all precomputed
  // by computeEdgeInterpolations()
  void addVertexPairTriangle(int i, int j, int k, int index);
  void addVertexTripletTriangle(int i, int j, int k, VEC3I index);
  void addStreamingTriangle(int i, int j, int k, VEC3I index);
  void addStreamingTriangle64(int i, int j, int k, VEC3I index);

  // add a vertex pairs, to be interpolated later
  void addVertexPairs(int i, int j, int k, int index);
  void addVertexTriplets(int i, int j, int k, VEC3I index);
  pair<int, int> getVertexPair(int i, int index);
  pair<VEC3I, VEC3I> getVertexTriplets(int i, VEC3I v);
 
  // don't check for duplicate collisions, just accumulate
  void accumulateVertexTriplets(int i, int j, int k, VEC3I index);

  // get the edge point
  VEC3F computeVertex(int i, int index);
  VEC3F computeStagedVertex(int i, int index, const CUBE& cube, const VEC3F& center);

  // see if the vertex has been computed before, and if not, store it
  int storeVertex(VEC3F& vertex, int index);

  // compute the non-linear interpolations for matching cubes
  void computeSingleSlabEdgeInterpolations(const vector<pair<VEC3I, VEC3I> >& newPairs,
                                           vector<VEC3F>& newVertices);
  void computeSingleSlabEdgeInterpolationsCuda(const vector<pair<VEC3I, VEC3I> >& newPairs,
                                               vector<VEC3F>& newVertices);
  void mallocEdgeCuda(const vector<pair<VEC3I, VEC3I> >& pairs);
  void freeEdgeCuda();

  // uses 64-bit ints for counting, in case there are more than 2 billion
  // vertices output in the mesh. Writing to a GZipped file as well, because
  // checkpointing is getting really expensive.
  void computeNonlinearMarchingCubesCudaStreamingCheckpointed64gz(const string& objFilename);
  void writeStreamingCheckpoint64gz(const string& objFilename,
                                    const int z,
                                    const int64 pairsSeen,
                                    const int64 triangleVerticesSeen,
                                    const FIELD_2D& slab0,
                                    const FIELD_2D& slab1,
                                    const map<pair<VEC3I, VEC3I>, bool>& uniquePairs0,
                                    const map<pair<VEC3I, VEC3I>, bool>& uniquePairs1,
                                    gzFile& vertexFile,
                                    gzFile& faceFile);
  void readStreamingCheckpoint64gz(const string& objFilename,
                                   int& firstSlab,
                                   int64& pairsSeen,
                                   int64& triangleVerticesSeen,
                                   FIELD_2D& slab0,
                                   FIELD_2D& slab1,
                                   map<pair<VEC3I, VEC3I>, bool>& uniquePairs0,
                                   map<pair<VEC3I, VEC3I>, bool>& uniquePairs1,
                                   gzFile& vertexFile,
                                   gzFile& faceFile);

  // do it on the GPU
  void initSliceCuda();
  void computeNonlinearSliceCuda(const int z);
  void readbackCudaSlice(FIELD_2D& field);
  void swapCudaBuffers();

  // compute flags on the GPU
  void computeFlagsCuda();
  void readbackCudaFlags(unsigned char* flags);

  // compute all the slices for a low memory marching cubes
  void computeAllLowMemorySlices(vector<pair<int, int> >& flags);
  void computeAllLowMemorySlicesHuge(vector<pair<int, VEC3I> >& flags);
  void computeAllSlicesCuda(vector<int>& flags, vector<VEC3I>& indices);
  void computeAllSlicesCudaCompacted(vector<int>& flags, vector<VEC3I>& indices);

  // compute some coarse normals for each vertex
  void computeNormals();

  // read/write an edge cache
  bool readEdgeCache(vector<pair<int, int> >& flags);
  void writeEdgeCache(const vector<pair<int, int> >& flags);
  bool readEdgeCacheHuge(vector<pair<int, VEC3I> >& flags);
  void writeEdgeCacheHuge(const vector<pair<int, VEC3I> >& flags);
  
  bool readEdgeCacheHuge(vector<int>& flags, vector<VEC3I>& indices);
  void writeEdgeCacheHuge(const vector<int>& flags, const vector<VEC3I>& indices);
  void writePartialEdgeCache(const int zCurrent, const vector<int>& flags, const vector<VEC3I>& indices);
  bool readPartialEdgeCache(int& zCurrent, vector<int>& flags, vector<VEC3I>& indices);

  void writePartialInterpolationCache(const int totalChunks, const int lastCompleted, const int verticesCompleted);
  bool readPartialInterpolationCache(int& totalChunks, int& lastCompleted);

  // get the (x,y,z) of an index
  VEC3I getXYZ(const int index) const;
  
  // get the nonlinear function value
  Real nonlinearValue(const VEC3F& position, const bool debug = false) const;
  
  // do a midpoint search
  VEC3F midpointSearch(const VEC3F& positiveVertex, const Real& positiveValue, const VEC3F& negativeVertex, const Real& negativeValue, const int recursion = 0) const;
  VEC3F midpointSearchForLoop(const VEC3F& positiveVertex, const Real& positiveValue, const VEC3F& negativeVertex, const Real& negativeValue) const;
  
  // emulate the cell center lookup for FIELD_3D --
  // assumes that _xRes, _yRes, _zRes,
  //              _center, _lengths, and _dx are populated
  VEC3F cellCenter(int x, int y, int z);
};

#endif
