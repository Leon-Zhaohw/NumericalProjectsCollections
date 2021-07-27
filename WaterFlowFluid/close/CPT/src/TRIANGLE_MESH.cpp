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
//
// TRIANGLE_MESH.h: interface for the TRIANGLE_MESH class.
//
//////////////////////////////////////////////////////////////////////

#include "TRIANGLE_MESH.h"
#include <algorithm>
#include <fstream>
#include <omp.h>
#include <TIMER.h>
#include <cstring>

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::TRIANGLE_MESH() :
  _glTextureHandle(0)
{
}

//////////////////////////////////////////////////////////////////////
// destructor
//////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::~TRIANGLE_MESH()
{
  // hand back any texture handles
  if (_glSolidTextureHandles.size() != 0)
  {
    // to speed things up, try to delete all the handles in one shot
    GLuint* handles = new GLuint[_glSolidTextureHandles.size()];

    for (unsigned int x = 0; x < _glSolidTextureHandles.size(); x++)
      handles[x] = _glSolidTextureHandles[x];
    
    glDeleteTextures(_glSolidTextureHandles.size(), handles);
    delete[] handles;

    _glSolidTextureHandles.clear();
  }
  if (_glTextureHandle)
    glDeleteTextures(1, &_glTextureHandle);
}

//////////////////////////////////////////////////////////////////////
// draw preview of solid textured version to GL,
// where the color of each triangle is set to the average texture value
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::drawSolidTexturePreview(Real amp, Real shift)
{
  assert(_triangles.size() == _solidTextureMeans.size());

  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    // recompute normals every time
    VEC3F normal = cross(*_triangles[x].vertex(1) - *_triangles[x].vertex(0), 
                         *_triangles[x].vertex(2) - *_triangles[x].vertex(0));
    normal.normalize();

    glBegin(GL_TRIANGLES);
      float color = _solidTextureMeans[x];

      color *= amp;
      color += shift;
      glColor4f(color, color, color, 1.0);
      glNormal3f(normal[0], normal[1], normal[2]);

      //int which = 0;
      int index = _vertexIndices[_triangles[x].vertex(0)];
      glVertex3f((*_triangles[x].vertex(0))[0], (*_triangles[x].vertex(0))[1], (*_triangles[x].vertex(0))[2]);

      index = _vertexIndices[_triangles[x].vertex(1)];
      glVertex3f((*_triangles[x].vertex(1))[0], (*_triangles[x].vertex(1))[1], (*_triangles[x].vertex(1))[2]);

      index = _vertexIndices[_triangles[x].vertex(2)];
      glVertex3f((*_triangles[x].vertex(2))[0], (*_triangles[x].vertex(2))[1], (*_triangles[x].vertex(2))[2]);
    glEnd();
  }
}

//////////////////////////////////////////////////////////////////////
// read an obj file using stdlib
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readOBJ(const string& filename)
{
  TIMER functionTimer(__FUNCTION__);

	// clear anything that existed before
	_vertices.clear();
	_normals.clear();
  _triangles.clear();
  _texcoords.clear();

	// open up file
  FILE* file = fopen(filename.c_str(), "r");
	
	if (file == NULL)
	{
		cerr << "Can't read input file " << filename << endl;
		return false;
	}

  // make a type-dependent format string
  string format;
  if (sizeof(Real) == sizeof(double))
    format = string("%lf %lf %lf");
  else
    format = string("%f %f %f");

	// read through line by line
	int lineNumber = 0;
	bool faceSeen = false;
  while (true)
	{
		if (feof(file)) 
      break;

		char type[1024];
		lineNumber++;
    fscanf(file, "%s", type);

		if (feof(file) || ferror(file)) break;

    // see if it's a comment
    if (type[0] == '#')
    {
      fgets(type, 1024, file);
      continue;
    }

		// reading vertices
		if (strcmp(type, "v") == 0)
		{
      if (faceSeen)
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Error! Vertex seen after face read begun!" << endl;
      }

			VEC3F v;
      fscanf(file, format.c_str(), &v[0], &v[1], &v[2]);
			_vertices.push_back(v);
		}

		// vertex normals
		if (strcmp(type, "vn") == 0)
		{
			VEC3F vn;
      fscanf(file, format.c_str(), &vn[0], &vn[1], &vn[2]);
			_normals.push_back(vn);
		}

		// vertex texcoords
		if (strcmp(type, "vt") == 0)
		{
			VEC3F vt;
      fscanf(file, format.c_str(), &vt[0], &vt[1], &vt[2]);
			_texcoords.push_back(vt);
    }

		// reading triangles
	  TRIANGLE f;
		if (type[0] == 'f')
		{
      char indices[256];
      faceSeen = true;
      
      for (int x = 0; x < 3; x++)
      {
        int totalChars = fscanf(file, "%s", indices); 
        if (feof(file) || ferror(file)) break;

        int vertexIndex = -1;
        int texcoordIndex = -1;
        int normalIndex = -1;
        char* texcoordExists = strchr(indices, '/');
        char* normalExists   = strrchr(indices, '/');

        int vertexEnd = texcoordExists - indices;
        int texcoordEnd = normalExists - indices;

        // extract the vertex index
        if (texcoordExists == NULL)
        {
          // if there is nothing but an index
          vertexIndex = atoi(indices);
        }
        else
        {
          // get the right substring
          char vertexString[256];
          strncpy(vertexString, indices, vertexEnd);

          // convert it to an int
          vertexIndex = atoi(vertexString);
        }

        // extract the texture index
        if (texcoordExists)
        {
          if (normalExists == NULL)
          {
            // extract to the end of the string
            char texcoordString[256];
            strncpy(texcoordString, &indices[vertexEnd + 1], totalChars - vertexEnd);
            
            // convert it to an int
            texcoordIndex = atoi(texcoordString);
          }
          else
          {
            // extract to the beginning of the normal index
            char texcoordString[256];
            strncpy(texcoordString, &indices[vertexEnd + 1], texcoordEnd - vertexEnd);
            
            // convert it to an int
            texcoordIndex = atoi(texcoordString);
          }
        }

        // extract the normal index
        if (normalExists)
          {
            // extract to the end of the string
            char normalString[256];
            strncpy(normalString, &indices[texcoordEnd + 1], totalChars - texcoordEnd);
              
            // convert it to an int
            normalIndex = atoi(normalString);
          }

        // subtract one and store
        f.vertex(x) = &_vertices[vertexIndex - 1];
      }

      if (feof(file) || ferror(file)) break;

		  // store the triangle
		  _triangles.push_back(f);
		}
	}

  fclose(file);
  cout << filename.c_str() << " successfully loaded" << endl;
  cout << " Vertices: " << _vertices.size() << endl;
  cout << " Faces: " << _triangles.size() << endl;

  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

	return true;
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

///////////////////////////////////////////////////////////////////////
// create a set of textures based on a solid texture
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::textureUsingSolidTexture(const FIELD_3D& solidTexture, const int textureRes , float factor)
{
  TIMER glTimer("Allocating GL texture handles");
  cout << " Solid texturing ... ";flush(cout);
  unsigned int size = _triangles.size();

  if (_solidTextures.size() != size)
    _solidTextures.resize(size);
  if (_solidTextureMeans.size() != size)
    _solidTextureMeans.resize(size);
  if (_texcoords.size() != 3 * size)
    _texcoords.resize(3 * size);
 
  // return any loaded textures to OGL
  if (_glSolidTextureHandles.size() != 0)
  {
    for (unsigned int x = 0; x < _glSolidTextureHandles.size(); x++)
      glDeleteTextures(1, &(_glSolidTextureHandles[x]));

    _glSolidTextureHandles.clear();
  }

  TIMER meansTimer("Computing GL texture means");
  int sizeCached = _triangles.size();
#pragma omp parallel
#pragma omp for  schedule(dynamic)
  for (int t = 0; t < sizeCached; t++)
  {
    // get the vertices
    vector<VEC3F> verts;
    verts.push_back(*(_triangles[t].vertices()[0]));
    verts.push_back(*(_triangles[t].vertices()[1]));
    verts.push_back(*(_triangles[t].vertices()[2]));
   
    VEC3F vertMean = verts[0] + verts[1] + verts[2];
    vertMean *= 1.0 / 3.0;

    _solidTextureMeans[t] = solidTexture(vertMean) * factor;
  }
  cout << " done. " << endl;

}

///////////////////////////////////////////////////////////////////////
// write out PBRT geometry
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::writePBRT(const char* filename)
{
  FILE* file = fopen(filename, "w");

  // read in the vertices
  fprintf(file, " Shape \"trianglemesh\"\n");
  fprintf(file, "\"point P\" [\n");
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    VEC3F vertex = _vertices[x];
    fprintf(file, "%f %f %f\n", vertex[0], vertex[1], vertex[2]);
  }
  fprintf(file, "]\n");

  fprintf(file, "\"integer indices\" [\n");
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    int indices[3];
    for (int y = 0; y < 3; y++)
      indices[y] = _vertexIndices[_triangles[x].vertex(y)];
    fprintf(file, "%i %i %i\n", indices[0], indices[1], indices[2]);
  }
  fprintf(file, "]\n");

  fclose(file);
}

///////////////////////////////////////////////////////////////////////
// load up a Houdini 12 frame
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::readHoudini12(const int frame, const string path, FIELD_3D& levelSet, VECTOR3_FIELD_3D& velocityField, TRIANGLE_MESH& surface, const bool forceMarch)
{
  // add a slash if needed
  string pathSlash = path;
  if (path[path.length() - 1] != '/')
    pathSlash += string("/");

  // add the frame number
  char buffer[256];

  string filenameLs, filenameVel, filenameTri;

  //filenameLs  = pathSlash + string("surf_out_") + string(buffer) + string(".geo");
  sprintf(buffer, "%i", frame);
  filenameLs  = pathSlash + string("surf") + string(buffer) + string(".geo");
  filenameVel = pathSlash + string("vel_out") + string(buffer) + string(".geo");
  sprintf(buffer, "%02i", frame);
  filenameTri = pathSlash + string("triangles.") + string(buffer) + string(".obj");

  // load the level set
  levelSet.readHoudini12Surf(filenameLs.c_str());

  // load the velocities
  velocityField.readHoudini12(filenameVel);

  // set the lengths the same as the level set
  velocityField.setLengths(levelSet.lengths());

  // 2X padding
  levelSet = levelSet.withAddedPadding(4);
  velocityField = velocityField.withAddedPadding(4);

  // marching cubes the field

  // try to read in a triangle mesh
  bool success = surface.readOBJ(filenameTri.c_str());
  if (!success || forceMarch)
  {
    surface.computeMarchingCubes(levelSet, true);
    surface.writeOBJ(filenameTri);
  }
}

///////////////////////////////////////////////////////////////////////
// load up a PhysBAM frame
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::readPhysBAMFrame(const int frame, const string path, FIELD_3D& levelSet, VECTOR3_FIELD_3D& velocityField, TRIANGLE_MESH& surface, const bool forceMarch, const Real scaleDistance)
{
  TIMER functionTimer(__FUNCTION__);

  // add a slash if needed
  string pathSlash = path;
  if (path[path.length() - 1] != '/')
    pathSlash += string("/");

  // add the frame number
  char buffer[256];

  string filenameLs, filenameVel, filenameTri;
  sprintf(buffer, "%i", frame);
  pathSlash += string(buffer) + string("/");

  //filenameLs  = pathSlash + string("levelset.gz");
  filenameLs  = pathSlash;
  filenameLs += string("levelset.gz");
  filenameVel = pathSlash;
  filenameVel += string("mac_velocities.gz");
  filenameTri = pathSlash;
  filenameTri += string("triangles.obj");

  // load the level set
  levelSet.readPhysBAMGz(filenameLs.c_str());

  levelSet *= scaleDistance;

  // load the velocities
  velocityField.readPhysBAMGz(levelSet, filenameVel.c_str());

  // marching cubes the field

  // try to read in a triangle mesh
  bool success = surface.readOBJ(filenameTri.c_str());
  if (!success || forceMarch)
  {
    surface.computeMarchingCubes(levelSet, true);
    surface.writeOBJ(filenameTri);
  }
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::writeOBJ(const string path, const int frame)
{
  TIMER functionTimer(__FUNCTION__);
  // build the filename
  string filename = path;

  // check for the slash
  if (filename[filename.length() - 1] != '/')
    filename = filename + string("/");

  // add the frame number
  filename = filename + string("surface.");
  char buffer[256];
  sprintf(buffer, "%04i", frame);
  filename = filename + string(buffer) + string(".obj");

  return writeOBJ(filename);
}

//////////////////////////////////////////////////////////////////////
// perform marching cubes
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeMarchingCubes(const FIELD_3D& field, const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);
  if (verbose)
    cout << " Marching cubes ...";flush(cout);

  // clear any previous front
  _vertices.clear();
  _triangles.clear();
  _vertexHash.clear();

  _fieldDeltas[0] = field.dx();
  _fieldDeltas[1] = field.dy();
  _fieldDeltas[2] = field.dz();

  // set "outside" to something a lot bigger than the known grid bounds
  _outside = field.center()[0] + field.lengths()[0] * 10000;
 
  _xRes = field.xRes();
  _yRes = field.yRes();
  _zRes = field.zRes();
  _slabSize = _xRes * _yRes;
  
  for (int z = 0; z < field.zRes() - 1; z++)
  {
    for (int y = 0; y < field.yRes() - 1; y++)
      for (int x = 0; x < field.xRes() - 1; x++) 
      {
        int index = x + y * field.xRes() + z * field.slabSize();

        _cellCenter = field.cellCenter(x,y,z);

        _NNN = field(x,y,z);
        _NNP = field(x,y,z + 1);
        _NPN = field(x,y + 1,z);
        _NPP = field(x,y + 1,z + 1);
        _PNN = field(x + 1,y,z);
        _PNP = field(x + 1,y,z+1);
        _PPN = field(x + 1,y + 1,z);
        _PPP = field(x + 1,y + 1,z + 1);
		
        // construct the flag
        int flag =    ((_NNN > 0) + 2 *   (_NNP > 0) + 4  * (_NPN > 0) +
                   8 * (_NPP > 0) + 16 *  (_PNN > 0) + 32 * (_PNP > 0) +
                   64 *(_PPN > 0) + 128 * (_PPP > 0));
		  
        switch (flag) 
        {
        case 0:  case 255: break;
        case 1: addTriangle(2,1,10,index); break;
        case 2: addTriangle(2,11,3,index); break;
        case 3: addTriangle(10,11,1,index); addTriangle(11,3,1,index); break;
        case 4: addTriangle(1,4,9,index); break;
        case 5: addTriangle(10,2,4,index); addTriangle(4,9,10,index); break;
        case 6: addTriangle(1,4,9,index); addTriangle(2,11,3,index); break;
        case 7: addTriangle(10,4,9,index); addTriangle(10,3,4,index); addTriangle(10,11,3,index); break;
        case 8: addTriangle(3,12,4,index); break;
        case 9: addTriangle(3,12,4,index); addTriangle(2,1,10,index); break;
        case 10: addTriangle(11,12,4,index); addTriangle(4,2,11,index); break;
        case 11: addTriangle(11,1,10,index); addTriangle(11,4,1,index); addTriangle(11,12,4,index); break;
        case 12: addTriangle(9,1,3,index); addTriangle(3,12,9,index); break;
        case 13: addTriangle(9,3,12,index); addTriangle(9,2,3,index); addTriangle(9,10,2,index); break;
        case 14: addTriangle(12,9,1,index); addTriangle(12,1,2,index); addTriangle(12,2,11,index); break;
        case 15: addTriangle(11,12,9,index); addTriangle(9,10,11,index); break;
        case 16: addTriangle(10,5,6,index); break; 
        case 17: addTriangle(5,6,2,index); addTriangle(2,1,5,index); break;
        case 18: addTriangle(10,5,6,index); addTriangle(2,11,3,index); break;
        case 19: addTriangle(1,5,6,index); addTriangle(1,6,11,index); addTriangle(1,11,3,index); break;
        case 20: addTriangle(10,5,6,index); addTriangle(1,4,9,index); break;
        case 21: addTriangle(2,4,9,index); addTriangle(2,9,5,index); addTriangle(2,5,6,index); break;
        case 22: addTriangle(2,11,3,index); addTriangle(1,4,9,index); addTriangle(10,5,6,index); break;
        case 23: addTriangle(5,6,11,index); addTriangle(9,5,11,index); addTriangle(11,3,9,index); addTriangle(3,4,9,index); break;
        case 24: addTriangle(10,5,6,index); addTriangle(3,12,4,index); break;
        case 25: addTriangle(5,6,2,index); addTriangle(2,1,5,index); addTriangle(3,12,4,index); break;
        case 26: addTriangle(11,12,4,index); addTriangle(4,2,11,index); addTriangle(10,5,6,index); break;
        case 27: addTriangle(5,6,11,index); addTriangle(11,1,5,index); addTriangle(1,11,12,index); addTriangle(12,4,1,index); break;
        case 28: addTriangle(9,1,3,index); addTriangle(3,12,9,index); addTriangle(10,5,6,index); break;
        case 29: addTriangle(3,12,9,index); addTriangle(9,2,3,index); addTriangle(2,9,5,index); addTriangle(5,6,2,index); break;
        case 30: addTriangle(12,9,1,index); addTriangle(12,1,2,index); addTriangle(12,2,11,index); addTriangle(10,5,6,index); break;
        case 31: addTriangle(12,9,5,index); addTriangle(12,5,6,index); addTriangle(12,6,11,index); break;
        case 32: addTriangle(6,7,11,index); break;
        case 33: addTriangle(6,7,11,index); addTriangle(2,1,10,index); break;
        case 34: addTriangle(6,7,3,index); addTriangle(3,2,6,index); break;
        case 35: addTriangle(3,1,10,index); addTriangle(3,10,6,index); addTriangle(3,6,7,index); break;
        case 36: addTriangle(6,7,11,index); addTriangle(1,4,9,index); break;
        case 37: addTriangle(10,2,4,index); addTriangle(4,9,10,index); addTriangle(6,7,11,index); break;
        case 38: addTriangle(6,7,3,index); addTriangle(3,2,6,index); addTriangle(1,4,9,index); break;
        case 39: addTriangle(4,9,3,index); addTriangle(3,9,10,index); addTriangle(10,7,3,index); addTriangle(10,6,7,index); break;
        case 40: addTriangle(6,7,11,index); addTriangle(3,12,4,index); break;
        case 41: addTriangle(2,1,10,index); addTriangle(3,12,4,index); addTriangle(6,7,11,index); break;
        case 42: addTriangle(2,6,7,index); addTriangle(2,7,12,index); addTriangle(2,12,4,index); break;
        case 43: addTriangle(1,12,4,index); addTriangle(1,7,12,index); addTriangle(1,10,7,index); addTriangle(10,6,7,index); break;
        case 44: addTriangle(9,1,3,index); addTriangle(3,12,9,index); addTriangle(6,7,11,index); break;
        case 45: addTriangle(9,3,12,index); addTriangle(9,2,3,index); addTriangle(9,10,2,index); addTriangle(6,7,11,index); break;
        case 46: addTriangle(9,6,12,index); addTriangle(12,6,7,index); addTriangle(9,2,6,index); addTriangle(9,1,2,index); break;
        case 47: addTriangle(9,10,6,index); addTriangle(9,6,7,index); addTriangle(9,7,12,index); break;
        case 48: addTriangle(11,10,5,index); addTriangle(5,7,11,index); break; 
        case 49: addTriangle(5,7,11,index); addTriangle(5,11,2,index); addTriangle(5,2,1,index); break; 
        case 50: addTriangle(7,3,2,index); addTriangle(7,2,10,index); addTriangle(7,10,5,index); break; 
        case 51: addTriangle(5,7,3,index); addTriangle(3,1,5,index); break; 
        case 52: addTriangle(11,10,5,index); addTriangle(5,7,11,index); addTriangle(1,4,9,index); break; 
        case 53: addTriangle(11,2,5,index); addTriangle(11,5,7,index); addTriangle(2,4,9,index); addTriangle(9,5,2,index); break; 
        case 54: addTriangle(7,3,2,index); addTriangle(7,2,10,index); addTriangle(7,10,5,index); addTriangle(1,4,9,index); break; 
        case 55: addTriangle(7,3,4,index); addTriangle(7,4,9,index); addTriangle(7,9,5,index); break; 
        case 56: addTriangle(11,10,5,index); addTriangle(5,7,11,index); addTriangle(3,12,4,index); break; 
        case 57: addTriangle(5,7,11,index); addTriangle(5,11,2,index); addTriangle(5,2,1,index); addTriangle(3,12,4,index); break; 
        case 58: addTriangle(2,10,5,index); addTriangle(2,5,7,index); addTriangle(7,4,2,index); addTriangle(7,12,4,index); break; 
        case 59: addTriangle(5,7,12,index); addTriangle(5,12,4,index); addTriangle(5,4,1,index); break; 
        case 60: addTriangle(9,1,3,index); addTriangle(3,12,9,index); addTriangle(11,10,5,index); addTriangle(5,7,11,index); break; 
        case 61: addTriangle(2,9,5,index); addTriangle(9,2,3,index); addTriangle(3,12,9,index); addTriangle(5,7,11,index); addTriangle(11,2,5,index); break; 
        case 62: addTriangle(2,7,12,index); addTriangle(7,2,10,index); addTriangle(10,5,7,index); addTriangle(12,9,1,index); addTriangle(1,2,12,index); break; 
        case 63: addTriangle(12,9,5,index); addTriangle(5,7,12,index); break; 
        case 64: addTriangle(9,8,5,index); break; 
        case 65: addTriangle(9,8,5,index); addTriangle(2,1,10,index); break; 
        case 66: addTriangle(9,8,5,index); addTriangle(2,11,3,index); break; 
        case 67: addTriangle(10,11,1,index); addTriangle(11,3,1,index); addTriangle(9,8,5,index); break; 
        case 68: addTriangle(1,4,8,index); addTriangle(8,5,1,index); break; 
        case 69: addTriangle(4,8,5,index); addTriangle(4,5,10,index); addTriangle(4,10,2,index); break; 
        case 70: addTriangle(1,4,8,index); addTriangle(8,5,1,index); addTriangle(2,11,3,index); break; 
        case 71: addTriangle(3,4,8,index); addTriangle(8,11,3,index); addTriangle(5,10,11,index); addTriangle(11,8,5,index); break; 
        case 72: addTriangle(9,8,5,index); addTriangle(3,12,4,index); break; 
        case 73: addTriangle(2,1,10,index); addTriangle(3,12,4,index); addTriangle(9,8,5,index); break; 
        case 74: addTriangle(11,12,4,index); addTriangle(4,2,11,index); addTriangle(9,8,5,index); break; 
        case 75: addTriangle(11,1,10,index); addTriangle(11,4,1,index); addTriangle(11,12,4,index); addTriangle(9,8,5,index); break; 
        case 76: addTriangle(1,3,12,index); addTriangle(1,12,8,index); addTriangle(1,8,5,index); break; 
        case 77: addTriangle(5,10,2,index); addTriangle(5,2,3,index); addTriangle(5,3,8,index); addTriangle(8,3,12,index); break; 
        case 78: addTriangle(12,1,2,index); addTriangle(2,11,12,index); addTriangle(1,12,8,index); addTriangle(8,5,1,index); break; 
        case 79: addTriangle(11,12,8,index); addTriangle(11,8,5,index); addTriangle(11,5,10,index); break; 
        case 80: addTriangle(6,10,9,index); addTriangle(9,8,6,index); break; 
        case 81: addTriangle(6,2,1,index); addTriangle(6,1,9,index); addTriangle(6,9,8,index); break; 
        case 82: addTriangle(6,10,9,index); addTriangle(9,8,6,index); addTriangle(2,11,3,index); break; 
        case 83: addTriangle(8,6,11,index); addTriangle(11,3,8,index); addTriangle(3,1,9,index); addTriangle(9,8,3,index); break; 
        case 84: addTriangle(8,6,10,index); addTriangle(8,10,1,index); addTriangle(8,1,4,index); break; 
        case 85: addTriangle(2,4,8,index); addTriangle(8,6,2,index); break; 
        case 86: addTriangle(8,6,10,index); addTriangle(8,10,1,index); addTriangle(8,1,4,index); addTriangle(2,11,3,index); break; 
        case 87: addTriangle(8,6,11,index); addTriangle(8,11,3,index); addTriangle(8,3,4,index); break; 
        case 88: addTriangle(6,10,9,index); addTriangle(9,8,6,index); addTriangle(3,12,4,index); break; 
        case 89: addTriangle(6,2,1,index); addTriangle(6,1,9,index); addTriangle(6,9,8,index); addTriangle(3,12,4,index); break; 
        case 90: addTriangle(11,12,4,index); addTriangle(4,2,11,index); addTriangle(6,10,9,index); addTriangle(9,8,6,index); break; 
        case 91: addTriangle(1,6,11,index); addTriangle(11,12,4,index); addTriangle(4,1,11,index); addTriangle(6,1,9,index); addTriangle(9,8,6,index); break; 
        case 92: addTriangle(3,6,1,index); addTriangle(6,3,8,index); addTriangle(3,12,8,index); addTriangle(6,10,1,index); break; 
        case 93: addTriangle(6,2,3,index); addTriangle(6,3,12,index); addTriangle(6,12,8,index); break; 
        case 94: addTriangle(2,11,12,index); addTriangle(2,12,1,index); addTriangle(1,12,8,index); addTriangle(10,1,8,index); addTriangle(10,8,6,index); break; 
        case 95: addTriangle(8,6,11,index); addTriangle(11,12,8,index); break; 
        case 96: addTriangle(9,8,5,index); addTriangle(6,7,11,index); break; 
        case 97: addTriangle(2,1,10,index); addTriangle(6,7,11,index); addTriangle(9,8,5,index); break; 
        case 98: addTriangle(6,7,3,index); addTriangle(3,2,6,index); addTriangle(9,8,5,index); break; 
        case 99: addTriangle(3,1,10,index); addTriangle(3,10,6,index); addTriangle(3,6,7,index); addTriangle(9,8,5,index); break; 
        case 100: addTriangle(1,4,8,index); addTriangle(8,5,1,index); addTriangle(6,7,11,index); break; 
        case 101: addTriangle(4,8,5,index); addTriangle(4,5,10,index); addTriangle(4,10,2,index); addTriangle(6,7,11,index); break; 
        case 102: addTriangle(6,7,3,index); addTriangle(3,2,6,index); addTriangle(1,4,8,index); addTriangle(8,5,1,index); break; 
        case 103: addTriangle(6,7,3,index); addTriangle(10,6,3,index); addTriangle(5,4,8,index); addTriangle(10,4,5,index); addTriangle(10,3,4,index); break; 
        case 104: addTriangle(9,8,5,index); addTriangle(6,7,11,index); addTriangle(3,12,4,index); break; 
        case 105: addTriangle(2,1,10,index); addTriangle(3,12,4,index); addTriangle(6,7,11,index); addTriangle(9,8,5,index); break; 
        case 106: addTriangle(2,6,7,index); addTriangle(2,7,12,index); addTriangle(2,12,4,index); addTriangle(9,8,5,index); break; 
        case 107: addTriangle(9,8,5,index); addTriangle(1,12,4,index); addTriangle(1,7,12,index); addTriangle(1,10,7,index); addTriangle(10,6,7,index); break; 
        case 108: addTriangle(1,3,12,index); addTriangle(1,12,8,index); addTriangle(1,8,5,index); addTriangle(11,6,7,index); break; 
        case 109: addTriangle(2,3,8,index); addTriangle(8,5,2,index); addTriangle(3,12,8,index); addTriangle(5,10,2,index); addTriangle(6,7,11,index); break; 
        case 110: addTriangle(2,6,7,index); addTriangle(2,7,12,index); addTriangle(1,12,8,index); addTriangle(1,8,5,index); addTriangle(1,2,12,index); break; 
        case 111: addTriangle(10,6,7,index); addTriangle(10,7,12,index); addTriangle(12,8,5,index); addTriangle(12,5,10,index); break; 
        case 112: addTriangle(10,9,8,index); addTriangle(10,8,7,index); addTriangle(10,7,11,index); break; 
        case 113: addTriangle(8,7,11,index); addTriangle(9,8,2,index); addTriangle(8,11,2,index); addTriangle(1,9,2,index); break; 
        case 114: addTriangle(3,9,8,index); addTriangle(3,8,7,index); addTriangle(2,9,3,index); addTriangle(2,10,9,index); break; 
        case 115: addTriangle(3,1,9,index); addTriangle(3,9,8,index); addTriangle(3,8,7,index); break; 
        case 116: addTriangle(8,7,11,index); addTriangle(8,11,4,index); addTriangle(10,4,11,index); addTriangle(4,10,1,index); break; 
        case 117: addTriangle(4,8,7,index); addTriangle(4,7,11,index); addTriangle(4,11,2,index); break; 
        case 118: addTriangle(3,2,7,index); addTriangle(2,10,7,index); addTriangle(1,4,8,index); addTriangle(10,8,7,index); addTriangle(10,1,8,index); break; 
        case 119: addTriangle(8,7,3,index); addTriangle(3,4,8,index); break; case 120: addTriangle(10,9,8,index); addTriangle(10,8,7,index); addTriangle(10,7,11,index); addTriangle(3,12,4,index); break;
        case 121: addTriangle(8,7,11,index); addTriangle(9,8,2,index); addTriangle(8,11,2,index); addTriangle(1,9,2,index); addTriangle(3,12,4,index); break; 
        case 122: addTriangle(2,12,4,index); addTriangle(2,7,12,index); addTriangle(10,9,8,index); addTriangle(10,8,7,index); addTriangle(2,10,7,index); break; 
        case 123: addTriangle(1,7,12,index); addTriangle(1,12,4,index); addTriangle(7,1,9,index); addTriangle(7,9,8,index); break; 
        case 124: addTriangle(10,1,8,index); addTriangle(10,8,7,index); addTriangle(7,11,10,index); addTriangle(1,3,12,index); addTriangle(12,8,1,index); break; 
        case 125: addTriangle(2,8,7,index); addTriangle(7,11,2,index); addTriangle(3,12,8,index); addTriangle(8,2,3,index); break; 
        case 126: addTriangle(12,8,7,index); addTriangle(10,1,2,index); break; 
        case 127: addTriangle(12,8,7,index); break; 
        case 128: addTriangle(7,8,12,index); break; 
        case 129: addTriangle(7,8,12,index); addTriangle(2,1,10,index); break; 
        case 130: addTriangle(7,8,12,index); addTriangle(2,11,3,index); break; 
        case 131: addTriangle(10,11,1,index); addTriangle(11,3,1,index); addTriangle(7,8,12,index); break; 
        case 132: addTriangle(1,4,9,index); addTriangle(7,8,12,index); break; 
        case 133: addTriangle(10,2,4,index); addTriangle(4,9,10,index); addTriangle(7,8,12,index); break; 
        case 134: addTriangle(2,11,3,index); addTriangle(1,4,9,index); addTriangle(7,8,12,index); break; 
        case 135: addTriangle(10,4,9,index); addTriangle(10,3,4,index); addTriangle(10,11,3,index); addTriangle(7,8,12,index); break; 
        case 136: addTriangle(3,7,8,index); addTriangle(8,4,3,index); break; 
        case 137: addTriangle(3,7,8,index); addTriangle(8,4,3,index); addTriangle(2,1,10,index); break; 
        case 138: addTriangle(4,2,11,index); addTriangle(4,11,7,index); addTriangle(4,7,8,index); break; 
        case 139: addTriangle(11,7,8,index); addTriangle(8,4,10,index); addTriangle(10,11,8,index); addTriangle(4,1,10,index); break; 
        case 140: addTriangle(3,7,8,index); addTriangle(3,8,9,index); addTriangle(3,9,1,index); break; 
        case 141: addTriangle(3,7,8,index); addTriangle(8,9,3,index); addTriangle(2,3,9,index); addTriangle(9,10,2,index); break; 
        case 142: addTriangle(11,7,8,index); addTriangle(2,11,8,index); addTriangle(2,8,9,index); addTriangle(2,9,1,index); break; 
        case 143: addTriangle(10,11,7,index); addTriangle(10,7,8,index); addTriangle(10,8,9,index); break; 
        case 144: addTriangle(7,8,12,index); addTriangle(10,5,6,index); break; 
        case 145: addTriangle(5,6,2,index); addTriangle(2,1,5,index); addTriangle(7,8,12,index); break; 
        case 146: addTriangle(2,11,3,index); addTriangle(10,5,6,index); addTriangle(7,8,12,index); break; 
        case 147: addTriangle(1,5,6,index); addTriangle(1,6,11,index); addTriangle(1,11,3,index); addTriangle(7,8,12,index); break; 
        case 148: addTriangle(1,4,9,index); addTriangle(10,5,6,index); addTriangle(7,8,12,index); break; 
        case 149: addTriangle(2,4,9,index); addTriangle(2,9,5,index); addTriangle(2,5,6,index); addTriangle(7,8,12,index); break; 
        case 150: addTriangle(7,8,12,index); addTriangle(2,11,3,index); addTriangle(1,4,9,index); addTriangle(10,5,6,index); break; 
        case 151: addTriangle(5,6,11,index); addTriangle(9,5,11,index); addTriangle(11,3,9,index); addTriangle(3,4,9,index); addTriangle(7,8,12,index); break; 
        case 152: addTriangle(3,7,8,index); addTriangle(8,4,3,index); addTriangle(10,5,6,index); break; 
        case 153: addTriangle(5,6,2,index); addTriangle(2,1,5,index); addTriangle(3,7,8,index); addTriangle(8,4,3,index); break;
        case 154: addTriangle(4,2,11,index); addTriangle(4,11,7,index); addTriangle(4,7,8,index); addTriangle(10,5,6,index); break; 
        case 155: addTriangle(11,4,1,index); addTriangle(7,8,4,index); addTriangle(4,11,7,index); addTriangle(1,5,6,index); addTriangle(6,11,1,index); break; 
        case 156: addTriangle(3,7,8,index); addTriangle(3,8,9,index); addTriangle(3,9,1,index); addTriangle(10,5,6,index); break; 
        case 157: addTriangle(9,2,3,index); addTriangle(3,7,8,index); addTriangle(8,9,3,index); addTriangle(5,6,2,index); addTriangle(2,9,5,index); break; 
        case 158: addTriangle(11,7,8,index); addTriangle(2,11,8,index); addTriangle(2,8,9,index); addTriangle(2,9,1,index); break; 
        case 159: addTriangle(8,9,11,index); addTriangle(11,7,8,index); addTriangle(9,6,11,index); addTriangle(9,5,6,index); break; 
        case 160: addTriangle(11,6,8,index); addTriangle(8,12,11,index); break; 
        case 161: addTriangle(11,6,8,index); addTriangle(8,12,11,index); addTriangle(2,1,10,index); break; 
        case 162: addTriangle(6,8,12,index); addTriangle(6,12,3,index); addTriangle(6,3,2,index); break; 
        case 163: addTriangle(8,12,3,index); addTriangle(3,1,6,index); addTriangle(6,8,3,index); addTriangle(1,10,6,index); break; 
        case 164: addTriangle(11,6,8,index); addTriangle(8,12,11,index); addTriangle(1,4,9,index); break; 
        case 165: addTriangle(10,2,4,index); addTriangle(4,9,10,index); addTriangle(11,6,8,index); addTriangle(8,12,11,index); break; 
        case 166: addTriangle(6,8,12,index); addTriangle(6,12,3,index); addTriangle(6,3,2,index); addTriangle(1,4,9,index); break; 
        case 167: addTriangle(6,3,10,index); addTriangle(6,8,12,index); addTriangle(12,3,6,index); addTriangle(4,9,10,index); addTriangle(10,3,4,index); break; 
        case 168: addTriangle(8,4,3,index); addTriangle(8,3,11,index); addTriangle(8,11,6,index); break; 
        case 169: addTriangle(8,4,3,index); addTriangle(8,3,11,index); addTriangle(8,11,6,index); addTriangle(2,1,10,index); break; 
        case 170: addTriangle(2,6,8,index); addTriangle(8,4,2,index); break; 
        case 171: addTriangle(8,4,1,index); addTriangle(8,1,10,index); addTriangle(8,10,6,index); break; 
        case 172: addTriangle(3,11,6,index); addTriangle(6,8,3,index); addTriangle(1,3,8,index); addTriangle(8,9,1,index); break; 
        case 173: addTriangle(9,3,8,index); addTriangle(8,3,11,index); addTriangle(11,6,8,index); addTriangle(2,3,9,index); addTriangle(9,10,2,index); break; 
        case 174: addTriangle(6,8,9,index); addTriangle(6,9,1,index); addTriangle(6,1,2,index); break; 
        case 175: addTriangle(9,10,6,index); addTriangle(6,8,9,index); break; 
        case 176: addTriangle(11,10,5,index); addTriangle(11,5,8,index); addTriangle(11,8,12,index); break; 
        case 177: addTriangle(2,1,5,index); addTriangle(5,11,2,index); addTriangle(11,5,8,index); addTriangle(8,12,11,index); break; 
        case 178: addTriangle(2,10,5,index); addTriangle(2,5,8,index); addTriangle(8,3,2,index); addTriangle(8,12,3,index); break; 
        case 179: addTriangle(1,5,8,index); addTriangle(1,8,12,index); addTriangle(1,12,3,index); break; 
        case 180: addTriangle(11,10,5,index); addTriangle(11,5,8,index); addTriangle(11,8,12,index); addTriangle(1,4,9,index); break; 
        case 181: addTriangle(5,11,2,index); addTriangle(11,5,8,index); addTriangle(8,12,11,index); addTriangle(2,4,9,index); addTriangle(9,5,2,index); break; 
        case 182: addTriangle(2,10,5,index); addTriangle(2,5,8,index); addTriangle(8,3,2,index); addTriangle(8,12,3,index); addTriangle(1,4,9,index); break; 
        case 183: addTriangle(3,4,9,index); addTriangle(9,5,3,index); addTriangle(5,8,12,index); addTriangle(12,3,5,index); break; 
        case 184: addTriangle(11,10,5,index); addTriangle(5,8,11,index); addTriangle(8,4,3,index); addTriangle(3,11,8,index); break; 
        case 185: addTriangle(8,11,5,index); addTriangle(5,11,2,index); addTriangle(2,1,5,index); addTriangle(8,4,3,index); addTriangle(3,11,8,index); break; 
        case 186: addTriangle(4,2,10,index); addTriangle(4,10,5,index); addTriangle(4,5,8,index); break; 
        case 187: addTriangle(8,4,1,index); addTriangle(1,5,8,index); break; 
        case 188: addTriangle(8,3,11,index); addTriangle(8,9,3,index); addTriangle(9,1,3,index); addTriangle(5,8,11,index); addTriangle(11,10,5,index); break; 
        case 189: addTriangle(5,8,9,index); addTriangle(3,11,2,index); break; 
        case 190: addTriangle(2,10,5,index); addTriangle(5,8,2,index); addTriangle(8,9,1,index); addTriangle(1,2,8,index); break; 
        case 191: addTriangle(5,8,9,index); break; 
        case 192: addTriangle(5,9,12,index); addTriangle(12,7,5,index); break; 
        case 193: addTriangle(5,9,12,index); addTriangle(12,7,5,index); addTriangle(2,1,10,index); break; 
        case 194: addTriangle(5,9,12,index); addTriangle(12,7,5,index); addTriangle(2,11,3,index); break; 
        case 195: addTriangle(10,11,1,index); addTriangle(11,3,1,index); addTriangle(5,9,12,index); addTriangle(12,7,5,index); break; 
        case 196: addTriangle(5,1,4,index); addTriangle(5,4,12,index); addTriangle(5,12,7,index); break; 
        case 197: addTriangle(4,12,7,index); addTriangle(7,5,2,index); addTriangle(5,10,2,index); addTriangle(2,4,7,index); break; 
        case 198: addTriangle(5,1,4,index); addTriangle(5,4,12,index); addTriangle(5,12,7,index); addTriangle(2,11,3,index); break; 
        case 199: addTriangle(10,4,5,index); addTriangle(10,11,3,index); addTriangle(3,4,10,index); addTriangle(5,4,12,index); addTriangle(12,7,5,index); break; 
        case 200: addTriangle(7,5,9,index); addTriangle(7,9,4,index); addTriangle(7,4,3,index); break; 
        case 201: addTriangle(7,5,9,index); addTriangle(7,9,4,index); addTriangle(7,4,3,index); addTriangle(2,1,10,index); break; 
        case 202: addTriangle(5,9,4,index); addTriangle(4,7,5,index); addTriangle(4,2,7,index); addTriangle(7,2,11,index); break; 
        case 203: addTriangle(4,11,7,index); addTriangle(11,4,1,index); addTriangle(1,10,11,index); addTriangle(7,5,9,index); addTriangle(9,4,7,index); break; 
        case 204: addTriangle(1,3,7,index); addTriangle(7,5,1,index); break; 
        case 205: addTriangle(7,5,10,index); addTriangle(7,10,2,index); addTriangle(7,2,3,index); break; 
        case 206: addTriangle(5,1,2,index); addTriangle(5,2,11,index); addTriangle(5,11,7,index); break; 
        case 207: addTriangle(11,7,5,index); addTriangle(5,10,11,index); break; 
        case 208: addTriangle(9,12,7,index); addTriangle(9,7,6,index); addTriangle(9,6,10,index); break; 
        case 209: addTriangle(2,1,9,index); addTriangle(9,12,2,index); addTriangle(2,12,6,index); addTriangle(12,7,6,index); break; 
        case 210: addTriangle(9,12,7,index); addTriangle(9,7,6,index); addTriangle(9,6,10,index); addTriangle(2,11,3,index); break; 
        case 211: addTriangle(1,9,6,index); addTriangle(9,12,7,index); addTriangle(7,6,9,index); addTriangle(1,6,11,index); addTriangle(11,3,1,index); break; 
        case 212: addTriangle(4,12,7,index); addTriangle(1,4,7,index); addTriangle(7,6,1,index); addTriangle(6,10,1,index); break; 
        case 213: addTriangle(2,4,12,index); addTriangle(2,12,7,index); addTriangle(2,7,6,index); break; 
        case 214: addTriangle(4,12,7,index); addTriangle(1,4,7,index); addTriangle(7,6,1,index); addTriangle(6,10,1,index); addTriangle(2,11,3,index); break; 
        case 215: addTriangle(3,4,6,index); addTriangle(6,11,3,index); addTriangle(4,12,7,index); addTriangle(7,6,4,index); break; 
        case 216: addTriangle(7,6,10,index); addTriangle(3,7,10,index); addTriangle(10,9,3,index); addTriangle(9,4,3,index); break;
        case 217: addTriangle(9,7,6,index); addTriangle(1,9,6,index); addTriangle(6,2,1,index); addTriangle(4,3,7,index); addTriangle(7,9,4,index); break; 
        case 218: addTriangle(9,4,7,index); addTriangle(6,10,9,index); addTriangle(9,7,6,index); addTriangle(4,2,11,index); addTriangle(11,7,4,index); break; 
        case 219: addTriangle(11,7,6,index); addTriangle(9,4,1,index); break; 
        case 220: addTriangle(3,7,6,index); addTriangle(3,6,10,index); addTriangle(3,10,1,index); break; 
        case 221: addTriangle(6,2,3,index); addTriangle(3,7,6,index); break; 
        case 222: addTriangle(11,7,1,index); addTriangle(1,2,11,index); addTriangle(7,6,10,index); addTriangle(10,1,7,index); break; 
        case 223: addTriangle(11,7,6,index); break;
        case 224: addTriangle(12,11,6,index); addTriangle(12,6,5,index); addTriangle(12,5,9,index); break; 
        case 225: addTriangle(12,11,6,index); addTriangle(12,6,5,index); addTriangle(12,5,9,index); addTriangle(2,1,10,index); break; 
        case 226: addTriangle(6,12,3,index); addTriangle(3,2,6,index); addTriangle(12,6,5,index); addTriangle(5,9,12,index); break; 
        case 227: addTriangle(6,12,3,index); addTriangle(5,9,12,index); addTriangle(12,6,5,index); addTriangle(10,6,3,index); addTriangle(3,1,10,index); break; 
        case 228: addTriangle(1,4,12,index); addTriangle(1,12,11,index); addTriangle(1,11,5,index); addTriangle(11,6,5,index); break; 
        case 229: addTriangle(4,12,5,index); addTriangle(12,11,6,index); addTriangle(6,5,12,index); addTriangle(4,5,10,index); addTriangle(10,2,4,index); break; 
        case 230: addTriangle(12,6,5,index); addTriangle(6,12,3,index); addTriangle(3,2,6,index); addTriangle(5,1,4,index); addTriangle(4,12,5,index); break; 
        case 231: addTriangle(6,5,10,index); addTriangle(4,12,3,index); break; 
        case 232: addTriangle(11,6,5,index); addTriangle(11,5,9,index); addTriangle(9,3,11,index); addTriangle(9,4,3,index); break; 
        case 233: addTriangle(11,6,5,index); addTriangle(11,5,9,index); addTriangle(9,3,11,index); addTriangle(9,4,3,index); addTriangle(2,1,10,index); break; 
        case 234: addTriangle(2,6,5,index); addTriangle(2,5,9,index); addTriangle(2,9,4,index); break; 
        case 235: addTriangle(10,6,4,index); addTriangle(4,1,10,index); addTriangle(9,4,6,index); addTriangle(6,5,9,index); break; 
        case 236: addTriangle(1,3,11,index); addTriangle(1,11,6,index); addTriangle(1,6,5,index); break; 
        case 237: addTriangle(6,5,3,index); addTriangle(3,11,6,index); addTriangle(2,3,5,index); addTriangle(5,10,2,index); break; 
        case 238: addTriangle(2,6,5,index); addTriangle(5,1,2,index); break; 
        case 239: addTriangle(6,5,10,index); break; 
        case 240: addTriangle(11,10,9,index); addTriangle(9,12,11,index); break; 
        case 241: addTriangle(12,11,2,index); addTriangle(12,2,1,index); addTriangle(12,1,9,index); break; 
        case 242: addTriangle(9,12,3,index); addTriangle(9,3,2,index); addTriangle(9,2,10,index); break; 
        case 243: addTriangle(3,1,9,index); addTriangle(9,12,3,index); break; 
        case 244: addTriangle(11,10,1,index); addTriangle(11,1,4,index); addTriangle(11,4,12,index); break; 
        case 245: addTriangle(4,12,11,index); addTriangle(11,2,4,index); break; 
        case 246: addTriangle(2,10,12,index); addTriangle(12,3,2,index); addTriangle(4,12,10,index); addTriangle(10,1,4,index); break; 
        case 247: addTriangle(4,12,3,index); break; 
        case 248: addTriangle(10,9,4,index); addTriangle(10,4,3,index); addTriangle(10,3,11,index); break; 
        case 249: addTriangle(11,9,4,index); addTriangle(4,3,11,index); addTriangle(2,1,9,index); addTriangle(9,11,2,index); break; 
        case 250: addTriangle(4,2,10,index); addTriangle(10,9,4,index); break; 
        case 251: addTriangle(9,4,1,index); break; 
        case 252: addTriangle(1,11,10,index); addTriangle(1,3,11,index); break; 
        case 253: addTriangle(3,11,2,index); break; 
        case 254: addTriangle(10,1,2,index); break;
        }
      }
  }

  // create the final triangles based on the vertex indices -- this
  // couldn't be done in the inner loop because the vector keeps
  // resizing and changing the vertex addresses
  assert(_triangleVertices.size() % 3 == 0);
	
  if (verbose) cout << "computed triangles: " << _triangleVertices.size() << endl;
	
  for (unsigned int x = 0; x < _triangleVertices.size() / 3; x++)
  {
    VEC3F* v0 = &_vertices[_triangleVertices[3 * x]];
    VEC3F* v1 = &_vertices[_triangleVertices[3 * x + 1]];
    VEC3F* v2 = &_vertices[_triangleVertices[3 * x + 2]];

    // ignore any degenerate triangles
    Real dist0 = norm((*v0) - (*v1));
    Real dist1 = norm((*v0) - (*v2));
    Real dist2 = norm((*v1) - (*v2));
    Real eps = 1e-7;
    if (dist0 < eps || dist1 < eps || dist2 < eps)
      continue;
    _triangles.push_back(TRIANGLE(v0, v1, v2));
  }

  // all done -- throw away the indices
  _triangleVertices.clear();

  // rebuild the vertex hash
  _vertexIndices.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexIndices[&(_vertices[x])] = x;

  if (verbose)
    cout << "done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// add a triangle to the list
//////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addTriangle(int i, int j, int k, int index)
{
  //TIMER functionTimer(__FUNCTION__);
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
// get the edge point
//////////////////////////////////////////////////////////////////////
VEC3F TRIANGLE_MESH::computeVertex(int i, int index)
{
  //TIMER functionTimer(__FUNCTION__);
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
  //TIMER functionTimer(__FUNCTION__);

#if 1
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
#else
  // this does a slow linear search, but if it becomes a problem,
  // the obvious thing to do is to put a KD-tree search in instead
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    VEC3F diff = _vertices[x] - vertex;
    Real magnitude = norm(diff);

    // if they're really close to each other, just return that.
    if (magnitude < 1e-7)
      return x;
  }

  _vertices.push_back(vertex);

  return _vertices.size() - 1;
#endif
}

//////////////////////////////////////////////////////////////////////
// read an obj file using stdlib
//////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::writeOBJ(const string& filename)
{
	// try to open out stream
	ofstream out(filename.c_str ());
	if (out.fail())
	{
		cout << "Failed to open " << filename<< " to save OBJ" << endl;
		return false;
	}
	// spew vertices
	for (unsigned int i = 0; i < _vertices.size(); i++)
		out << "v " << _vertices[i][0] << " " << _vertices[i][1] << " " << _vertices[i][2] << endl;

	// normals
	for (unsigned int i = 0; i < _normals.size(); i++)
		out << "vn " << _normals[i] << endl;

	// faces
	for (unsigned int i = 0; i < _triangles.size(); i++)
	{
		out << "f ";
    for (int j = 0; j < 3; j++)
  	  out << _vertexIndices[_triangles[i].vertex(j)] + 1 << " ";
		out << endl;
	}

	// perfunctory error checking
	if (out.fail())
	{
		cerr << "There was an error writing " << filename << endl;
		return false;
	}
	out.close ();
	return true;
}
