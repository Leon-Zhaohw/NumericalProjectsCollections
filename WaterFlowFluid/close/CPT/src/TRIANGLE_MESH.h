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

#ifndef TRIANGLE_MESH_H
#define TRIANGLE_MESH_H

#include <SETTINGS.h>
#include <TRIANGLE.h>
#include "FIELD_3D.h"
#include "VECTOR3_FIELD_3D.h"
#include "SURFACE.h"

using namespace std;

class TRIANGLE_MESH : public SURFACE
{
public:
  TRIANGLE_MESH();

  // destructor
  ~TRIANGLE_MESH();

  // draw preview of solid textured version to GL,
  // where the color of each triangle is set to the average texture value
  void drawSolidTexturePreview(Real amp = 1.0, Real shift = 0);

  // center using vertex means
  VEC3F vertexMean() const;

  // return a bounding box
  void boundingBox(VEC3F& mins, VEC3F& maxs);

  // create a set of textures based on a solid texture
  void textureUsingSolidTexture(const FIELD_3D& solidTexture, const int textureRes = 5, float factor = 1.f );
  
  // perform marching cubes
  void computeMarchingCubes(const FIELD_3D& field, const bool verbose = false);

  ////////////////////////////////////////////////////////////////////////////
  // read/write support
  ////////////////////////////////////////////////////////////////////////////
  bool readOBJ(const string& filename);
  bool writeOBJ(const string& filename);
  bool writeOBJ(const string path, const int frame);
  void writePBRT(const char* filename);
  static void readPhysBAMFrame(const int frame, const string path, FIELD_3D& levelSet, VECTOR3_FIELD_3D& velocityField, TRIANGLE_MESH& surface, const bool forceMarch = false, const Real scaleDistance = 1);
  static void readHoudini12(const int frame, const string path, FIELD_3D& levelSet, VECTOR3_FIELD_3D& velocityField, TRIANGLE_MESH& surface, const bool forceMarch = false);

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
  // texturing support
  ////////////////////////////////////////////////////////////////////////////
 
  // the texture, if any
  FIELD_2D _texture;

  // GL handle for the texture
  GLuint _glTextureHandle;

  ////////////////////////////////////////////////////////////////////////////
  // solid texturing support
  ////////////////////////////////////////////////////////////////////////////

  // per triangle texture
  vector<FIELD_2D> _solidTextures;
  
  // per triangle mean solid texture values
  vector<float> _solidTextureMeans;
  
  // GL handle for the texture
  vector<GLuint> _glSolidTextureHandles;

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

  // add a triangle to the list
  void addTriangle(int i, int j, int k, int index);
  
  // get the edge point
  VEC3F computeVertex(int i, int index);
  
  // see if the vertex has been computed before, and if not, store it
  int storeVertex(VEC3F& vertex, int index);
};

#endif
