/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
//////////////////////////////////////////////////////////////////////
// obj.cxx/h
// 
// Created 2005.08 by christopher.cameron@gmail.com
// Contains code to load and save an OBJ 3d model
//
// Ted Kim started modifications as of 3/14/08. Particularly, support for PLY
// files was added, as well as fast distance queries.
//
// This does the ugly thing for PLY files -- it reads in the PLY into totally
// separate memory and then just copies everything to the OBJ arrays. I assume
// that we will not be using a PLY to gigantic that trying to have two copies
// momentarily in memory will be the limiting factor.
//////////////////////////////////////////////////////////////////////

#include "OBJ.h"
#include <fstream>
#include <sstream>
#include <map>
#include <stdlib.h>
#include <MERSENNETWISTER.h>
#if !USING_OSX
#include <GL/glut.h>
#endif
using namespace std;

//////////////////////////////////////////////////////////////////////
// copy constructor
//////////////////////////////////////////////////////////////////////
OBJ::OBJ(const OBJ& toCopy)
{
  vertices = toCopy.vertices;
  normals = toCopy.normals;
  texcoords = toCopy.texcoords;
  faces = toCopy.faces;
}

//////////////////////////////////////////////////////////////////////
// read an obj file
//////////////////////////////////////////////////////////////////////
bool OBJ::LoadPBRT(string& filename)
{
	vertices.resize (0);
	normals.resize (0);
	texcoords.resize (0);
  faces.resize(0);

  FILE* file = fopen(filename.c_str(), "r");
  cout << " Reading file " << filename.c_str() << endl;

  if (file == NULL)
  {
    cout << " File " << filename.c_str() << " not found! " << endl;
    exit(0);
  }

  // bite off the header
  int success = 0;
  success = fscanf(file, "AttributeBegin\n");
  success = fscanf(file, " Translate 0 0 0\n");
  success = fscanf(file, " Shape \"trianglemesh\"\n");
  success = fscanf(file, "\"point P\" [\n");

  // read in the vertices
  success = 3;
  while (success == 3)
  {
    float in[3];
    VEC3F vertex;
    success = fscanf(file, "%f %f %f\n", &(in[0]), &(in[1]), &(in[2]));
    vertex[0] = in[0];
    vertex[1] = in[1];
    vertex[2] = in[2];
    vertices.push_back(vertex);
  }
  cout << " Read in " << vertices.size() << " vertices " << endl;

  // bite off the face header
  success = fscanf(file, "]\n");
  success = fscanf(file, "\"integer indices\" [\n");
  success = 3;
  while (success == 3)
  {
    int in[3];
    success = fscanf(file, "%i %i %i\n", &(in[0]), &(in[1]), &(in[2]));
    Face face;
    face.vertices.push_back(in[0]);
    face.vertices.push_back(in[1]);
    face.vertices.push_back(in[2]);
    faces.push_back(face);
  }
  cout << " Read in " << faces.size() << " faces " << endl;
  fclose(file);

  return false;
}

//////////////////////////////////////////////////////////////////////
// read an obj file
//////////////////////////////////////////////////////////////////////
bool OBJ::Load(const string& fileName, Real scale)
{
	// clear anything that existed before
	vertices.resize (0);
	normals.resize (0);
	texcoords.resize (0);
  _bccRes = -1;

  // make a copy of the filename
  //_filename = string(fileName);
	
	// open up file
	ifstream in (fileName.c_str ());
	if (in.fail()) 
	{
		cerr << "Can't read input file " << fileName << endl;
		return false;
	}

	// read through line by line
	int lineNumber = 0;
	string type("");
	while (true)
	{
		if (in.eof () ) break;

		string line("");
		lineNumber++;
		getline (in, line);

		if (in.eof () ) break;
		if (in.fail () )
		{
			cerr << "Error reading" << endl;
      exit(1);
		}

		// get the type of command (v, vt, vn, or f supported)
		istringstream is (line);
		is >> type;
		if (is.eof () || is.fail () ) continue;

		// reading vertices
		if (type == "v")
		{
			VEC3 v;
			is >> v;

      // Scale the mesh to the appropriate size
      v *= scale;

			vertices.push_back (v);
		}
		if (type == "g")
		{
      _vertexGroupStarts.push_back(vertices.size());
      _faceGroupStarts.push_back(faces.size());
      _texcoordGroupStarts.push_back(texcoords.size());
      string name("");
      is >> name;

      _groupNames.push_back(name);
    }

    if (type == "usemtl")
    {
      _materialStarts.push_back(faces.size());
      cout << "material line: " << line << endl;
    }

		// vertex normals
		if (type == "vn")
		{
			VEC3 vn;
			is >> vn;
			normals.push_back (vn);
		}

		// texcoords
		if (type == "vt")
		{
			VEC2 vt;
			is >> vt;
			texcoords.push_back (vt);
		}

		// reading faces
		else if (type == "f")
		{
			Face f;
			while (is.eof () == false)
			{
				// read vertex index, texture index, and normal index as v/t/n
				string indexString("");
				is >> indexString;

				// sometimes there will be whitespace after lines, so we'll need
				// to read all the way past the last entry to get an eof, which
				// will result in a fail here
				if (is.fail () ) 
				{
					continue;
				}

				// read the indices of vertex/texcoord/normal
				int indicesIndex = 0;
				int indices[] = {0, 0, 0};
				int sign[] = {1, 1, 1};
				for (unsigned int k = 0; k < indexString.size (); k++)
				{
					char c = indexString[k];
					if (c == '\0') 
					{
						break;
					}
					else if (c == '/') 
					{
						if (indicesIndex == 2)
						{
							cerr << "Malformed face at line " << lineNumber << endl;
							goto Fail;
						}
						indicesIndex++;
					}
					else if (c >= '0' && c <= '9')
					{
						indices[indicesIndex] *= 10;
						indices[indicesIndex] += c - '0';
					}
					else if (c == '-' && indices[indicesIndex] == 0)
					{
						sign[indicesIndex] *= -1;
					}
					else 
					{
						cerr << "Malformed face at line " << lineNumber << endl;
						goto Fail;
					}
				}

				// handle offset-relative indices
				if (sign[0] == -1) indices[0] = vertices.size ()  + 1 - indices[0];
				if (sign[1] == -1) indices[1] = texcoords.size () + 1 - indices[1];
				if (sign[2] == -1) indices[2] = normals.size ()   + 1 - indices[2];

				// verify that indices given are sane
				if (indices[0] < 1 || indices[0] > (int)vertices.size () )
				{
					cerr << "Invalid vertex index at line " << lineNumber << endl;
					goto Fail;
				}
				if (indices[1] < 0 || indices[1] > (int)texcoords.size () )
				{
					cerr << "Invalid texcoord index at line " << lineNumber << endl;
					goto Fail;
				}
				if (indices[2] < 0 || indices[2] > (int)normals.size () )
				{
					cerr << "Invalid normal index at line " << lineNumber << endl;
					goto Fail;
				}

				// subtract one and store
				f.vertices.push_back (indices[0] - 1);
				f.texcoords.push_back   (indices[1] - 1);
				f.normals.push_back  (indices[2] - 1);
			}

		  // store the quad
		  faces.push_back (f);
		}
		if (in.fail () ) 
		{
			cerr << "Failure reading model at line " << lineNumber << endl;
			goto Fail;
		}
	}
  _vertexGroupStarts.push_back(vertices.size());
  _faceGroupStarts.push_back(faces.size());
  _texcoordGroupStarts.push_back(texcoords.size());

	in.close ();
  cout << fileName.c_str() << " successfully loaded" << endl;
  cout << vertices.size() << " total vertices " << endl;
  cout << faces.size() << " total faces " << endl;

  // cache the triangle and edge lengths
  _triangleAreas.clear();
  _maxEdgeLengths.clear();
  for (unsigned int x = 0; x < faces.size(); x++)
  {
    VEC3F* triangleVertices[3];
    for (int y = 0; y < 3; y++)
      triangleVertices[y] = &vertices[faces[x].vertices[y]];

    TRIANGLE triangle(triangleVertices[0], triangleVertices[1], triangleVertices[2]);
    Real area = triangle.area();
    Real maxLength = triangle.maxEdgeLength();

    _triangleAreas.push_back(area);
    _maxEdgeLengths.push_back(maxLength);
  }
  _filteringThreshold = 3.0;

	return true;

Fail:
	in.close ();
	vertices.resize (0);
	normals.resize (0);
	return false;
}

//////////////////////////////////////////////////////////////////////
// read an PLY file
//////////////////////////////////////////////////////////////////////
bool OBJ::LoadPly(const string& fileName)
{
  // Turk's declarations
  PlyFile *ply;
  int nelems;
  char **elist;
  int file_type;
  float version;
  int nprops;
  int num_elems;
  PlyProperty **plist;
  PlyVertex **vlist;
  PlyFace **flist;
  char *elem_name;
  int num_verts;
  int num_faces;

  PlyProperty vert_props[] = { // list of property information for a vertex
    {"x", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,x), 0, 0, 0, 0},
    {"y", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,y), 0, 0, 0, 0},
    {"z", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,z), 0, 0, 0, 0},
  };

  PlyProperty face_props[] = { // list of property information for a vertex
    {"vertex_indices", PLY_INT, PLY_INT, offsetof(PlyFace,verts),
     1, PLY_UCHAR, PLY_UCHAR, offsetof(PlyFace,nverts)},
  };

  // make a copy of the filename
  _filename = string(fileName);

  // open a PLY file for reading
  const char* charName = _filename.c_str();
  ply = 
    (PlyFile*)ply_open_for_reading(
      charName, 
      &nelems, 
      &elist, 
      &file_type, 
      &version);
 
  printf ("version %f\n", version);
  printf ("type %d\n", file_type);

  for (int i = 0; i < nelems; i++) {

    elem_name = elist[i];
    plist = ply_get_element_description (ply, elem_name, &num_elems, &nprops);

    printf ("element %s %d\n", elem_name, num_elems);

    if (equal_strings ("vertex", elem_name)) {
      num_verts = num_elems;
      
      vlist = (PlyVertex **) malloc (sizeof (PlyVertex *) * num_elems);

      ply_get_property (ply, elem_name, &vert_props[0]);
      ply_get_property (ply, elem_name, &vert_props[1]);
      ply_get_property (ply, elem_name, &vert_props[2]);

      // get the file pointer
      FILE* file = ply->fp;

      for (int j = 0; j < num_elems; j++) {
        float data;
        fread(&data, 4,1, file);

        float endian;
        char *lE = (char*) &data;
        char *bE = (char*) &endian;
        lE[0] = bE[3];
        lE[1] = bE[2];
        lE[2] = bE[1];
        lE[3] = bE[0];
      }
    }

    if (equal_strings ("face", elem_name)) {
      num_faces = num_elems;
      
      flist = (PlyFace **) malloc (sizeof (PlyFace *) * num_elems);

      ply_get_property (ply, elem_name, &face_props[0]);
    }
  }
  ply_close(ply);

  return true;
}

//////////////////////////////////////////////////////////////////////
// write a PLY file
/*
ply
format ascii 1.0           { ascii/binary, format version number }
comment made by Greg Turk  { comments keyword specified, like all lines }
comment this file is a cube
element vertex 8           { define "vertex" element, 8 of them in file }
property float x           { vertex contains float "x" coordinate }
property float y           { y coordinate is also a vertex property }
property float z           { z coordinate, too }
element face 6             { there are 6 "face" elements in the file }
property list uchar int vertex_index { "vertex_indices" is a list of ints }
end_header                 { delimits the end of the header }
0 0 0                      { start of vertex list }
0 0 1
0 1 1
0 1 0
1 0 0
1 0 1
1 1 1
1 1 0
4 0 1 2 3                  { start of face list }
4 7 6 5 4
4 0 4 5 1
4 1 5 6 2
4 2 6 7 3
4 3 7 4 0
*/
//////////////////////////////////////////////////////////////////////
bool OBJ::SavePly(const std::string& filename)
{
  FILE* file = fopen(filename.c_str(), "wb");

  fprintf(file, "ply\n");
  fprintf(file, "format ascii 1.0\n");
  fprintf(file, "element vertex %i\n", (int)vertices.size());
  fprintf(file, "property float x\n");
  fprintf(file, "property float y\n");
  fprintf(file, "property float z\n");
  fprintf(file, "element face %i\n", (int)faces.size());
  fprintf(file, "property list uchar int vertex_index\n");
  fprintf(file, "end_header\n");

  for (unsigned int x = 0; x < vertices.size(); x++)
    fprintf(file, "%f %f %f\n", vertices[x][0], vertices[x][1], vertices[x][2]);

  for (unsigned int x = 0; x < faces.size(); x++)
    fprintf(file, "3 %i %i %i\n", faces[x].vertices[0], faces[x].vertices[1], faces[x].vertices[2]);

  fclose(file);

  cout << filename.c_str() << " saved. " << endl;

  return true;
}

//////////////////////////////////////////////////////////////////////
// Save an OBJ
//////////////////////////////////////////////////////////////////////
bool OBJ::SaveFiltered(vector<int>& filteredFaces, const string& fileName)
{
	// try to open out stream
	ofstream  out (fileName.c_str () );
	if (out.fail () )
	{
		cerr << "Failed to open " << fileName << " to save OBJ" << endl;
		return false;
	}
  cout << " Writing out OBJ " << fileName.c_str() << endl;

	// spew vertices
	for (unsigned int i = 0; i < vertices.size (); i++)
	{
		out << "v " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << endl;
	}
	// normals
	for (unsigned int i = 0; i < normals.size (); i++)
	{
		out << "vn " << normals[i] << endl;
	}
	// faces
	for (unsigned int i = 0; i < filteredFaces.size (); i++)
	{
		out << "f ";
		for (unsigned int j = 0; j < faces[filteredFaces[i]].vertices.size (); j++)
		{
			out << faces[filteredFaces[i]].vertices[j] + 1;
			out << " ";

		}
		out << endl;
	}
	// perfunctory error checking
	if (out.fail () )
	{
		cerr << "There was an error writing " << fileName << endl;
		return false;
	}
	out.close ();
	return true;
}

//////////////////////////////////////////////////////////////////////
// Save an OBJ
//////////////////////////////////////////////////////////////////////
bool OBJ::Save (const string& fileName)
{
	// try to open out stream
	ofstream  out (fileName.c_str () );
	if (out.fail () )
	{
		cerr << "Failed to open " << fileName << " to save OBJ" << endl;
		return false;
	}
	// spew vertices
	for (unsigned int i = 0; i < vertices.size (); i++)
	{
		out << "v " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << endl;
	}
	// normals
	for (unsigned int i = 0; i < normals.size (); i++)
	{
		out << "vn " << normals[i] << endl;
	}
	// faces
	for (unsigned int i = 0; i < faces.size (); i++)
	{
		out << "f ";
		for (unsigned int j = 0; j < faces[i].vertices.size (); j++)
		{
			out << faces[i].vertices[j] + 1;
			out << " ";

		}
		out << endl;
	}
	// perfunctory error checking
	if (out.fail () )
	{
		cerr << "There was an error writing " << fileName << endl;
		return false;
	}
	out.close ();
	return true;
}

//////////////////////////////////////////////////////////////////////
// get bounding box
//////////////////////////////////////////////////////////////////////
void OBJ::BoundingBox (VEC3& minP, VEC3& maxP)
{
	if (vertices.size () == 0)
	{
		minP = -1.0;
		maxP =  1.0;
		return;
	}
	minP = maxP = vertices[0];
	for (unsigned int i = 0; i < vertices.size (); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			minP[j] = min (minP[j], vertices[i][j]);
			maxP[j] = max (maxP[j], vertices[i][j]);
		}
	}
}

//////////////////////////////////////////////////////////////////////
// compute smooth vertex normals
//////////////////////////////////////////////////////////////////////
void OBJ::ComputeVertexNormals()
{
	normals.resize (vertices.size () );
	for (unsigned int i = 0; i < normals.size (); i++)
	{
		normals[i] = 0.0;
	}
	for (unsigned int i = 0; i < faces.size (); i++)
	{
		faces[i].normals = faces[i].vertices;
		if (0)
		{
			VEC3 normal = cross (
				vertices[faces[i].vertices[1]] - vertices[faces[i].vertices[0]],
				vertices[faces[i].vertices[2]] - vertices[faces[i].vertices[0]]);
			for (unsigned int j = 0; j < faces[i].normals.size (); j++)
			{
				normals[faces[i].normals[j]] += normal;
			}
		}
		else
		{
			unsigned int n = faces[i].normals.size ();
			for (unsigned int j = 0; j < faces[i].normals.size (); j++)
			{
				const VEC3& v0 = vertices[faces[i].vertices[(j+n-1)%n]];
				const VEC3& v1 = vertices[faces[i].vertices[j]];
				const VEC3& v2 = vertices[faces[i].vertices[(j+1)%n]];
				VEC3 e1 = v0 - v1; unitize(e1);
				VEC3 e2 = v2 - v1; unitize(e2);
				VEC3 normal = cross (e1, e2); unitize (normal);
#ifdef SINGLE_PRECISION
				float cosAngle = e1 * e2;
				float angle = acos (cosAngle);
#else
			  double cosAngle = e1 * e2;
				double angle = acos (cosAngle);
#endif
				normals[faces[i].vertices[j]] -= angle * normal;
			}
		}
	}
	for (unsigned int i = 0; i < normals.size (); i++)
	{
		unitize (normals[i]);
	}
}


//////////////////////////////////////////////////////////////////////
// smooth out vertex normals
//////////////////////////////////////////////////////////////////////
void OBJ::SmoothVertexNormals()
{
  // compute the one ring of each vertex
  vector<map<int, bool> > oneRings;
  oneRings.resize(vertices.size());

  for (unsigned int x = 0; x < faces.size(); x++)
  {
    Face face = faces[x];

    // store each pair
    for (unsigned int y = 0; y < face.vertices.size(); y++)
      for (unsigned int z = 0; z < face.vertices.size(); z++)
      {
        if (y == z) continue;
        oneRings[face.vertices[y]][face.vertices[z]] = true;
        oneRings[face.vertices[z]][face.vertices[y]] = true;
      }
  }

  vector<VEC3> normalsCopy;
  normalsCopy = normals;

  // average each normal with its one ring
  for (unsigned int x = 0; x < normals.size(); x++)
  {
    // get the one ring
    map<int, bool> oneRing = oneRings[x];

    // use the current normal
    VEC3 normal = normalsCopy[x];

    // add all the normals in the one ring too
    map<int,bool>::iterator iter;
    for (iter = oneRing.begin(); iter != oneRing.end(); iter++)
      normal += normalsCopy[iter->first];

    unitize(normal);
    normals[x] = normal;
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void OBJ::FlipVertexNormals()
{
  for (unsigned int x = 0; x < normals.size(); x++)
    normals[x] *= -1.0;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VEC3 OBJ::GetVertex( const Face& f, unsigned int i ) const
{
	unsigned int vertIdx = f.vertices[i];
	return vertices[ vertIdx ];
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VEC3 OBJ::GetCentroid( const Face& f ) const
{
	VEC3 centroid = 0;
	for( unsigned int j = 0; j < f.vertices.size(); j++ )
	{
		centroid += GetVertex( f, j );
	}
	centroid *= 1.0/(double)f.vertices.size();

	return centroid;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
OBJ::Face& OBJ::ClosestFace( const VEC3& p )
{
	double bestDistSq = 0.0;
	int bestFaceIdx = -1;

	for( unsigned int i = 0; i < faces.size(); i++ )
	{
		Face& f = faces[i];

		double distSq = norm2( GetCentroid(f) - p );

		if( distSq < bestDistSq || bestFaceIdx == -1 )
		{
			bestDistSq = distSq;
			bestFaceIdx = i;
		}
	}

	return faces[ bestFaceIdx ];
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
float OBJ::DistanceToClosestFace( const VEC3& p )
{
	double bestDistSq = 0.0;
	int firstSeen = -1;

	for( unsigned int i = 0; i < faces.size(); i++ )
	{
		Face& f = faces[i];

		//double distSq = norm2( GetCentroid(f) - p );
    float distSq = pointFaceDistanceSq(f, p);

		if( distSq < bestDistSq || firstSeen == -1 )
    {
			bestDistSq = distSq;
      firstSeen = 0;
    }
	}

	return bestDistSq;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int OBJ::ClosestVertex( const VEC3& p )
{
	double bestDistSq = 0.0;
	int bestVertIdx = -1;

	for( unsigned int i = 0; i < vertices.size(); i++ )
	{
		VEC3& v = vertices[i];

		double distSq = norm2( v - p );

		if( distSq < bestDistSq || bestVertIdx == -1 )
		{
			bestDistSq = distSq;
			bestVertIdx = i;
		}
	}

	return bestVertIdx;
}

//////////////////////////////////////////////////////////////////////
// fit the geometry inside a 1x1x1 box
//////////////////////////////////////////////////////////////////////
void OBJ::normalize(VEC3F& vertex)
{
  if (_scaling < 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need to normalize the mesh before normalizing other vertices!" << endl;
    return;
  }

  vertex -= _centerOfMass;
  vertex *= _scaling;

  VEC3F half(0.5, 0.5, 0.5);
  vertex += half;
}

//////////////////////////////////////////////////////////////////////
// fit the geometry inside a 1x1x1 box
//////////////////////////////////////////////////////////////////////
void OBJ::normalize(int res)
{
  // first get the center of mass
  VEC3F centerOfMass;
  for (unsigned int x = 0; x < vertices.size(); x++)
    centerOfMass += vertices[x];
  centerOfMass *= 1.0 / vertices.size();

  // translate everything to the center of mass
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] -= centerOfMass;

  // find the maximum magnitude
  double maxVal = 0.0f;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    maxVal = (fabs(vertices[x][0]) > maxVal) ? fabs(vertices[x][0]) : maxVal;
    maxVal = (fabs(vertices[x][1]) > maxVal) ? fabs(vertices[x][1]) : maxVal;
    maxVal = (fabs(vertices[x][2]) > maxVal) ? fabs(vertices[x][2]) : maxVal;
  }

  // scale everything
  double scale = 0.5 - 4.0 / res;
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] *= scale / maxVal;

  // translate everything to 0.5, 0.5, 0.5
  VEC3F half(0.5, 0.5, 0.5);
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] += half;

  // store in case other vertices need to be normalized by the same factors
  _centerOfMass = centerOfMass;
  _scaling = scale / maxVal;

  cout << " center of mass: " << _centerOfMass << endl;
  cout << " scaling: " << _scaling << endl;
  cout << " max val: " << maxVal << endl;
}

//////////////////////////////////////////////////////////////////////
// Do a brute force search in order to get nearest distance
//////////////////////////////////////////////////////////////////////
float OBJ::bruteForceDistance(float* point)
{
  // exhaustive search
  float minDistance = _bccRes * 1000;
  VEC3F vertex(point);

  for (unsigned int x = 0; x < faces.size(); x++)
  {
    Real dist = pointFaceDistanceSq(faces[x], vertex);
    if (dist < minDistance) minDistance = dist;
  }

  return minDistance;
}

//////////////////////////////////////////////////////////////////////
// implements OBJECT virtual function to get distance to nearest
// triangle
//////////////////////////////////////////////////////////////////////
float OBJ::distance(float* point)
{
  // broader bucketed search
  float minDistance = _bccRes * 1000;

  // do explicit tests against triangles in the hash
  int xCoord = (int)(point[0] * _xRes);
  int yCoord = (int)(point[1] * _yRes);
  int zCoord = (int)(point[2] * _zRes);

  int xStart = xCoord - 2;
  xStart = (xCoord == 1) ? xCoord - 1 : xStart;
  xStart = (xCoord == 0) ? xCoord : xStart;
  int xEnd   =  xCoord + 2;
  xEnd = (xCoord == _xRes - 2) ? xCoord + 1 : xEnd;
  xEnd = (xCoord == _xRes - 1) ? xCoord : xEnd;
  int yStart = yCoord - 2;
  yStart = (yCoord == 1) ? yCoord - 1 : yStart;
  yStart = (yCoord == 0) ? yCoord : yStart;
  int yEnd   =  yCoord + 2;
  yEnd = (yCoord == _yRes - 2) ? yCoord + 1 : yEnd;
  yEnd = (yCoord == _yRes - 1) ? yCoord : yEnd;
  int zStart = zCoord - 2;
  zStart = (zCoord == 1) ? zCoord - 1 : zStart;
  zStart = (zCoord == 0) ? zCoord : zStart;
  int zEnd   =  zCoord + 2;
  zEnd = (zCoord == _zRes - 2) ? zCoord + 1 : zEnd;
  zEnd = (zCoord == _zRes - 1) ? zCoord : zEnd;

  for (int z = zStart; z < zEnd; z++)
    for (int y = yStart; y < yEnd; y++)
      for (int x = xStart; x < xEnd; x++)
      {
        int index = x +
                    y * _xRes +
                    z * _xRes  * _yRes;
        vector<int>& bucket = _accelGrid[index];

        // do an explicit test with the triangles in the bucket
        VEC3F vertex(point);
        for (unsigned int i = 0; i < bucket.size(); i++)
        {
          int faceIndex = bucket[i];
          Real dist = sqrt(pointFaceDistanceSq(faces[faceIndex], vertex));
          if (dist < minDistance) minDistance = dist;
        }
      }
  return minDistance;
}

#ifndef EPSILON
#define EPSILON 1e-8
#endif
#ifndef CROSS
#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#endif
#ifndef DOT
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#endif
#ifndef SUB
#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2]; 
#endif

//////////////////////////////////////////////////////////////////////
// ray-triangle test from
// Tomas Möller and Ben Trumbore. Fast, minimum storage ray-triangle intersection. 
// Journal of graphics tools, 2(1):21-28, 1997
//////////////////////////////////////////////////////////////////////
bool OBJ::intersect_triangle(float orig[3], float dir[3], int faceIndex)
{
   float edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
   float det,inv_det;
   float t,u,v;
   float vert0[3], vert1[3], vert2[3];
   vert0[0] = vertices[faces[faceIndex].vertices[0]][0];
   vert0[1] = vertices[faces[faceIndex].vertices[0]][1];
   vert0[2] = vertices[faces[faceIndex].vertices[0]][2];

   vert1[0] = vertices[faces[faceIndex].vertices[1]][0];
   vert1[1] = vertices[faces[faceIndex].vertices[1]][1];
   vert1[2] = vertices[faces[faceIndex].vertices[1]][2];

   vert2[0] = vertices[faces[faceIndex].vertices[2]][0];
   vert2[1] = vertices[faces[faceIndex].vertices[2]][1];
   vert2[2] = vertices[faces[faceIndex].vertices[2]][2];
   
   // find vectors for two edges sharing vert0
   SUB(edge1, vert1, vert0);
   SUB(edge2, vert2, vert0);

   // begin calculating determinant - also used to calculate U parameter
   CROSS(pvec, dir, edge2);

   // if determinant is near zero, ray lies in plane of triangle
   det = DOT(edge1, pvec);

   if (det > -EPSILON && det < EPSILON)
     return false;
   inv_det = 1.0 / det;

   // calculate distance from vert0 to ray origin
   SUB(tvec, orig, vert0);

   // calculate U parameter and test bounds
   u = DOT(tvec, pvec) * inv_det;
   if (u < 0.0 || u > 1.0)
     return false;

   // prepare to test V parameter
   CROSS(qvec, tvec, edge1);

   // calculate V parameter and test bounds
   v = DOT(dir, qvec) * inv_det;
   if (v < 0.0 || u + v > 1.0)
     return false;

   t = DOT(edge2, qvec) * inv_det;

   if (t < 0.0f) return false;
   return true;
}

//////////////////////////////////////////////////////////////////////
// intersection test
//////////////////////////////////////////////////////////////////////
bool OBJ::insideExhaustive(float* point)
{
  unsigned int j;
  float vec[3];

  // try all cardinal directions, take majority
  vec[0] = 1.0f; vec[1] = 1.0f; vec[2] = 0.0f;
  int xPos = 0;
  for (j = 0; j < faces.size(); j++)
    if (intersect_triangle(point, vec, j))
      xPos++;

  // y positive
  vec[0] = 0.0f; vec[1] = 1.0f; vec[2] = 0.0f;
  int yPos = 0;
  for (j = 0; j < faces.size(); j++)
    if (intersect_triangle(point, vec, j))
      yPos++;

  // z positive
  vec[0] = 0.0f; vec[1] = 0.0f; vec[2] = 1.0f;
  int zPos = 0;
  for (j = 0; j < faces.size(); j++)
    if (intersect_triangle(point, vec, j))
      zPos++;

  // x negative
  vec[0] = -1.0f; vec[1] = 0.0f; vec[2] = 0.0f;
  int xNeg = 0;
  for (j = 0; j < faces.size(); j++)
    if (intersect_triangle(point, vec, j))
      xNeg++;
 
  // y negative
  vec[0] = 0.0f; vec[1] = -1.0f; vec[2] = 0.0f;
  int yNeg = 0;
  for (j = 0; j < faces.size(); j++)
    if (intersect_triangle(point, vec, j))
      yNeg++;
  
  // z negative
  vec[0] = 0.0f; vec[1] = 0.0f; vec[2] = -1.0f;
  int zNeg = 0;
  for (j = 0; j < faces.size(); j++)
    if (intersect_triangle(point, vec, j))
      zNeg++;
  
  float sum = (xPos % 2 + yPos % 2 + zPos % 2 + xNeg % 2 + yNeg % 2 + zNeg % 2);
 
  // take a majority vote 
  return (sum > 3);
}

//////////////////////////////////////////////////////////////////////
// translate the entire model
//////////////////////////////////////////////////////////////////////
void OBJ::translate(VEC3F translation)
{
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] += translation;
}

//////////////////////////////////////////////////////////////////////
// do a lookup into the occupancy grid
//////////////////////////////////////////////////////////////////////
bool OBJ::isOccupied(Real* point)
{
  float rotated[] = {point[0], point[1], point[2]};
  rotate(rotated);

  // integer version of the point
  int intPoint[] = {(int)(rotated[0] * _xRes), 
                    (int)(rotated[1] * _yRes), 
                    (int)(rotated[2] * _zRes)};

  // retrieve the faces that hash to that cell
  int index = intPoint[0] + intPoint[1] * _xRes + intPoint[2] * _xRes * _yRes;

  if (_occupancyGrid.find(index) == _occupancyGrid.end())
    return false;
  return true;
}

//////////////////////////////////////////////////////////////////////
// do a lookup into the accel grid
//////////////////////////////////////////////////////////////////////
bool OBJ::gridHit(Real* point)
{
  float rotated[] = {point[0], point[1], point[2]};
  rotate(rotated);

  // integer version of the point
  int intPoint[] = {(int)(rotated[0] * _xRes), 
                    (int)(rotated[1] * _yRes), 
                    (int)(rotated[2] * _zRes)};

  // retrieve the faces that hash to that cell
  int index = intPoint[0] + intPoint[1] * _xRes + intPoint[2] * _xRes * _yRes;
  vector<int>& faces = _accelGrid[index];

  if (faces.size() == 0)
    return false;
  return true;
}

//////////////////////////////////////////////////////////////////////
// intersection test
//////////////////////////////////////////////////////////////////////
bool OBJ::inside(float* point)
{
  float rotated[] = {point[0], point[1], point[2]};
  rotate(rotated);

  float vec[3];

  // integer version of the point
  int intPoint[] = {(int)(rotated[0] * _xRes), 
                    (int)(rotated[1] * _yRes), 
                    (int)(rotated[2] * _zRes)};

  // make a map of hits in the each direction so
  // that we don't double-count
  map<int, bool> xNegHits;
  map<int, bool> yNegHits;
  map<int, bool> zNegHits;
  map<int, bool> xPosHits;
  map<int, bool> yPosHits;
  map<int, bool> zPosHits;
  
  // counters for the total hits in each direction
  map<int, bool>::iterator counter;
  int xNegTotal = 0;
  int yNegTotal = 0;
  int zNegTotal = 0;
  int xPosTotal = 0;
  int yPosTotal = 0;
  int zPosTotal = 0;

  ////////////////////////////////////////////////////////////////////
  // Do the x direction
  ////////////////////////////////////////////////////////////////////

  // form the x direction vector
  vec[0] = -1.0f; vec[1] = 0.0f; vec[2] = 0.0f;

  // for each cell in the negative x direction
  for (int x = 0; x < intPoint[0]; x++)
  {
    // retrieve the faces that hash to that cell
    int index = x + intPoint[1] * _xRes + intPoint[2] * _xRes * _yRes;
    vector<int>& faces = _accelGrid[index];
    
    // do intersection test on the faces
    for (unsigned int i = 0; i < faces.size(); i++)
    {
      if (intersect_triangle(rotated, vec, faces[i]))
        xNegHits[faces[i]] = true;
    }
  }

  // count how many x hits there were
  for(counter = xNegHits.begin(); counter != xNegHits.end(); counter++)
    xNegTotal++;

  // if lots of holes start showing up in the mesh, try activating these
  // tests as well
  // form the x direction vector
  vec[0] = 1.0f; vec[1] = 0.0f; vec[2] = 0.0f;

  // for each cell in the negative x direction
  for (int x = intPoint[0]; x <= _xRes; x++)
  {
    // retrieve the faces that hash to that cell
    int index = x + intPoint[1] * _xRes + intPoint[2] * _xRes * _yRes;
    vector<int>& faces = _accelGrid[index];
    
    // do intersection test on the faces
    for (unsigned int i = 0; i < faces.size(); i++)
    {
      if (intersect_triangle(rotated, vec, faces[i]))
        xPosHits[faces[i]] = true;
    }
  }

  // count how many x hits there were
  for(counter = xPosHits.begin(); counter != xPosHits.end(); counter++)
    xPosTotal++;

  ////////////////////////////////////////////////////////////////////
  // Do the y direction
  ////////////////////////////////////////////////////////////////////

  // form the y direction vector
  vec[0] = 0.0f; vec[1] = -1.0f; vec[2] = 0.0f;

  // for each cell in the negative y direction
  for (int y = 0; y < intPoint[1]; y++)
  {
    // retrieve the faces that hash to that cell
    int index = intPoint[0] + y * _xRes + intPoint[2] * _xRes * _yRes;
    vector<int>& faces = _accelGrid[index];
    
    // do intersection test on the faces
    for (unsigned int i = 0; i < faces.size(); i++)
      if (intersect_triangle(rotated, vec, faces[i]))
        yNegHits[faces[i]] = true;
  }

  // count how many y hits there were
  for(counter = yNegHits.begin(); counter != yNegHits.end(); counter++)
    yNegTotal++;

  // form the y direction vector
  vec[0] = 0.0f; vec[1] = 1.0f; vec[2] = 0.0f;

  // for each cell in the negative y direction
  for (int y = intPoint[0]; y < _yRes; y++)
  {
    // retrieve the faces that hash to that cell
    int index = intPoint[0] + y * _xRes + intPoint[2] * _xRes * _yRes;
    vector<int>& faces = _accelGrid[index];
    
    // do intersection test on the faces
    for (unsigned int i = 0; i < faces.size(); i++)
      if (intersect_triangle(rotated, vec, faces[i]))
        yPosHits[faces[i]] = true;
  }

  // count how many y hits there were
  for(counter = yPosHits.begin(); counter != yPosHits.end(); counter++)
    yPosTotal++;
  
  ////////////////////////////////////////////////////////////////////
  // Do the z direction
  ////////////////////////////////////////////////////////////////////

  // form the z direction vector
  vec[0] = 0.0f; vec[1] = 0.0f; vec[2] = -1.0f;

  // for each cell in the negative z direction
  for (int z = 0; z < intPoint[2]; z++)
  {
    // retrieve the faces that hash to that cell
    int index = intPoint[0] + intPoint[1] * _xRes + z * _xRes * _yRes;
    vector<int>& faces = _accelGrid[index];
    
    // do intersection test on the faces
    for (unsigned int i = 0; i < faces.size(); i++)
      if (intersect_triangle(rotated, vec, faces[i]))
        zNegHits[faces[i]] = true;
  }

  // count how many z hits there were
  for(counter = zNegHits.begin(); counter != zNegHits.end(); counter++)
    zNegTotal++;

  // form the z direction vector
  vec[0] = 0.0f; vec[1] = 0.0f; vec[2] = 1.0f;

  // for each cell in the negative z direction
  for (int z = intPoint[2]; z < _zRes; z++)
  {
    // retrieve the faces that hash to that cell
    int index = intPoint[0] + intPoint[1] * _xRes + z * _xRes * _yRes;
    vector<int>& faces = _accelGrid[index];
    
    // do intersection test on the faces
    for (unsigned int i = 0; i < faces.size(); i++)
      if (intersect_triangle(rotated, vec, faces[i]))
        zPosHits[faces[i]] = true;
  }

  // count how many z hits there were
  for(counter = zPosHits.begin(); counter != zPosHits.end(); counter++)
    zPosTotal++;

  int insideVotes = (xNegTotal % 2 + yNegTotal % 2 + zNegTotal % 2 +
                     xPosTotal % 2 + yPosTotal % 2 + zPosTotal % 2);

  return insideVotes >= 3;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void OBJ::drawOccupancyGrid()
{
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
      {
        int index = x + y * _xRes + z * _xRes * _yRes;

        // make sure it's occupied and on the outside
        if (_occupancyGrid.find(index) == _occupancyGrid.end()) continue;
        if (_occupancyGrid.find(index + 1) != _occupancyGrid.end() &&
            _occupancyGrid.find(index - 1) != _occupancyGrid.end() &&
            _occupancyGrid.find(index + _xRes) != _occupancyGrid.end() &&
            _occupancyGrid.find(index - _xRes) != _occupancyGrid.end() &&
            _occupancyGrid.find(index + _xRes * _yRes) != _occupancyGrid.end() &&
            _occupancyGrid.find(index - _xRes * _yRes) != _occupancyGrid.end())
          continue;

        glPushMatrix();
          glTranslatef(x * 1.0 / _xRes, y * 1.0 / _yRes, z * 1.0 / _zRes);
          glutSolidCube(1.0 / _xRes);
        glPopMatrix();
      }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void OBJ::drawAccelGrid()
{
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++)
      {
        int index = x + y * _xRes + z * _xRes * _yRes;

        if (_accelGrid[index].size() != 0)
        {
          glPushMatrix();
            glTranslatef(x * 1.0 / _xRes, y * 1.0 / _yRes, z * 1.0 / _zRes);
            glutSolidCube(1.0 / _xRes);
          glPopMatrix();
        }
      }
}

//////////////////////////////////////////////////////////////////////
// hash all the triangles to the occupancy grid
//////////////////////////////////////////////////////////////////////
void OBJ::createOccupancyGrid(int fatten)
{
  if (_bccRes == -1)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need to set BCC grid res!!!" << endl;
    return;
  }

  // fatten factor of 3 definitely gives usable results with 128^3 cart
  cout << "Cube mesh fatten factor: " << fatten << endl;
  _xRes = _yRes = _zRes = 4 * _bccRes;
  cout << "occupancyGrid resolution: " << _xRes << endl;
  cout << "Building occupancyGrid ... ";
  flush(cout);
  
  // for each Face
  for (unsigned int i = 0; i < faces.size(); i++)
  {
    if (faces[i].vertices.size() != 3) continue;

    // retrieve the vertices
    VEC3 v0 = vertices[faces[i].vertices[0]];
    VEC3 v1 = vertices[faces[i].vertices[1]];
    VEC3 v2 = vertices[faces[i].vertices[2]];
    
    // calc the bounding box
    float xMin = v0[0];
    float xMax = v0[0];
    float yMin = v0[1];
    float yMax = v0[1];
    float zMin = v0[2];
    float zMax = v0[2];
    xMin = (xMin < v1[0]) ? xMin : v1[0];
    yMin = (yMin < v1[1]) ? yMin : v1[1];
    zMin = (zMin < v1[2]) ? zMin : v1[2];
    xMax = (xMax > v1[0]) ? xMax : v1[0];
    yMax = (yMax > v1[1]) ? yMax : v1[1];
    zMax = (zMax > v1[2]) ? zMax : v1[2];
    
    xMin = (xMin < v2[0]) ? xMin : v2[0];
    yMin = (yMin < v2[1]) ? yMin : v2[1];
    zMin = (zMin < v2[2]) ? zMin : v2[2];
    xMax = (xMax > v2[0]) ? xMax : v2[0];
    yMax = (yMax > v2[1]) ? yMax : v2[1];
    zMax = (zMax > v2[2]) ? zMax : v2[2];
   
    // for each cell in the bounding box
    for (int z = (int)(floor(zMin * _zRes)) - fatten; z <= (int)(ceil(zMax * _zRes)) + fatten; z++)
      for (int y = (int)(floor(yMin * _yRes)) - fatten; y <= (int)(ceil(yMax * _yRes)) + fatten; y++)
        for (int x = (int)(floor(xMin * _xRes)) - fatten; x <= (int)(ceil(xMax * _xRes)) + fatten; x++)
        {
          if (x < 0 || y < 0 || z < 0) continue;
          if (x >= _xRes || y >= _yRes || z >= _zRes) continue;
          int index = x + y * _xRes + z * _xRes * _yRes;
          
          // hash the face into the vector for the cell
          _occupancyGrid[index] = true;
        }

    if (i % (int)(faces.size() / 10) == 0)
    {
      cout << 100 * ((Real)i / faces.size()) << "% ";
      flush(cout);
    }
  }
  cout << "done." << endl;
}

//////////////////////////////////////////////////////////////////////
// hash all the triangles to the accel grid
//////////////////////////////////////////////////////////////////////
void OBJ::createAccelGrid()
{
  if (_bccRes == -1)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need to set BCC grid res!!!" << endl;
    return;
  }

  _xRes = _yRes = _zRes = _bccRes / 2;
  
  // for each Face
  for (unsigned int i = 0; i < faces.size(); i++)
  {
    if (faces[i].vertices.size() != 3) continue;

    // retrieve the vertices
    VEC3 v0 = vertices[faces[i].vertices[0]];
    VEC3 v1 = vertices[faces[i].vertices[1]];
    VEC3 v2 = vertices[faces[i].vertices[2]];
    
    // calc the bounding box
    float xMin = v0[0];
    float xMax = v0[0];
    float yMin = v0[1];
    float yMax = v0[1];
    float zMin = v0[2];
    float zMax = v0[2];
    xMin = (xMin < v1[0]) ? xMin : v1[0];
    yMin = (yMin < v1[1]) ? yMin : v1[1];
    zMin = (zMin < v1[2]) ? zMin : v1[2];
    xMax = (xMax > v1[0]) ? xMax : v1[0];
    yMax = (yMax > v1[1]) ? yMax : v1[1];
    zMax = (zMax > v1[2]) ? zMax : v1[2];
    
    xMin = (xMin < v2[0]) ? xMin : v2[0];
    yMin = (yMin < v2[1]) ? yMin : v2[1];
    zMin = (zMin < v2[2]) ? zMin : v2[2];
    xMax = (xMax > v2[0]) ? xMax : v2[0];
    yMax = (yMax > v2[1]) ? yMax : v2[1];
    zMax = (zMax > v2[2]) ? zMax : v2[2];
   
    // for each cell in the bounding box
    for (int z = (int)(floor(zMin * _zRes)); z <= (int)(ceil(zMax * _zRes)); z++)
      for (int y = (int)(floor(yMin * _yRes)); y <= (int)(ceil(yMax * _yRes)); y++)
        for (int x = (int)(floor(xMin * _xRes)); x <= (int)(ceil(xMax * _xRes)); x++)
        {
          if (x < 0 || y < 0 || z < 0) continue;
          if (x >= _xRes || y >= _yRes || z >= _zRes) continue;
          int index = x + y * _xRes + z * _xRes * _yRes;
          
          // hash the face into the vector for the cell
          _accelGrid[index].push_back(i);
        }

    if (i % (int)(faces.size() / 10) == 0)
    {
      cout << 100 * ((Real)i / faces.size()) << "% ";
      flush(cout);
    }
  }
  cout << "done." << endl;
}

//////////////////////////////////////////////////////////////////////
// hash all the triangles to the distance grid
//////////////////////////////////////////////////////////////////////
void OBJ::createDistanceGrid(int res)
{
  cout << "Building distance grid ...";
  _xResDist = _yResDist = _zResDist = 2 * res;

  // for each Face
  for (unsigned int i = 0; i < faces.size(); i++)
  {
    // retrieve the vertices
    VEC3 v0 = vertices[faces[i].vertices[0]];
    VEC3 v1 = vertices[faces[i].vertices[1]];
    VEC3 v2 = vertices[faces[i].vertices[2]];
    
    // calc the bounding box
    float xMin = v0[0];
    float xMax = v0[0];
    float yMin = v0[1];
    float yMax = v0[1];
    float zMin = v0[2];
    float zMax = v0[2];
    xMin = (xMin < v1[0]) ? xMin : v1[0];
    yMin = (yMin < v1[1]) ? yMin : v1[1];
    zMin = (zMin < v1[2]) ? zMin : v1[2];
    xMax = (xMax > v1[0]) ? xMax : v1[0];
    yMax = (yMax > v1[1]) ? yMax : v1[1];
    zMax = (zMax > v1[2]) ? zMax : v1[2];
    
    xMin = (xMin < v2[0]) ? xMin : v2[0];
    yMin = (yMin < v2[1]) ? yMin : v2[1];
    zMin = (zMin < v2[2]) ? zMin : v2[2];
    xMax = (xMax > v2[0]) ? xMax : v2[0];
    yMax = (yMax > v2[1]) ? yMax : v2[1];
    zMax = (zMax > v2[2]) ? zMax : v2[2];

    // fatten the footprint
    //float fatten = 1.0f;
    float fatten = 3.0f;
    //float fatten = 6.0f;
    if (xMin != 0)         xMin -= fatten / _xResDist;
    if (xMax != _xResDist) xMax += fatten / _xResDist;
    if (yMin != 0)         yMin -= fatten / _yResDist;
    if (yMax != _yResDist) yMax += fatten / _yResDist;
    if (zMin != 0)         zMin -= fatten / _zResDist;
    if (zMax != _xResDist) zMax += fatten / _zResDist;
   
    // for each cell in the bounding box
    for (int z = (int)(floor(zMin * _zResDist)); z <= (int)(ceil(zMax * _zResDist)); z++)
      for (int y = (int)(floor(yMin * _yResDist)); y <= (int)(ceil(yMax * _yResDist)); y++)
        for (int x = (int)(floor(xMin * _xResDist)); x <= (int)(ceil(xMax * _xResDist)); x++)
        {
          int index = x + y * _xResDist + z * _xResDist * _yResDist;

          float xFloat = (float)x / _xResDist;
          float yFloat = (float)y / _yResDist;
          float zFloat = (float)z / _zResDist;
          VEC3F point(xFloat, yFloat, zFloat);
         
          Real distanceSq = pointFaceDistanceSq(faces[i], point);

          map<int,Real>::iterator finder;
          finder = _distanceGrid.find(index);

          if (finder == _distanceGrid.end())
            _distanceGrid[index] = distanceSq;
          else
            _distanceGrid[index] = (distanceSq < _distanceGrid[index]) ? distanceSq : _distanceGrid[index];
        }
  }
  cout << " done." << endl;
}

//////////////////////////////////////////////////////////////////////////////
// use the _distanceGrid to get a distance value
//////////////////////////////////////////////////////////////////////////////
Real OBJ::fastDistToClosestFace(int vertexID)
{
  map<int,Real>::iterator finder;
  finder = _distanceGrid.find(vertexID);
  if (finder == _distanceGrid.end())
  {
    return 1e8;
  }
  return _distanceGrid[vertexID];
}

//////////////////////////////////////////////////////////////////////////////
// This is a modified version of the Wild Magic DistVector3Triangle3 class.
// This is what is used to generate the distance grid.
//
// The license info is as below:
//
// Wild Magic Source Code
// David Eberly
// http://www.geometrictools.com
// Copyright (c) 1998-2008
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or (at
// your option) any later version.  The license is available for reading at
// either of the locations:
//     http://www.gnu.org/copyleft/lgpl.html
//     http://www.geometrictools.com/License/WildMagicLicense.pdf
//
// Version: 4.0.1 (2007/05/06)
//----------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
Real OBJ::pointFaceDistanceSq(Face& f, const VEC3F& point)
{
  VEC3F v0 = vertices[f.vertices[0]];
  VEC3F v1 = vertices[f.vertices[1]];
  VEC3F v2 = vertices[f.vertices[2]];
  VEC3F kDiff = v0 - point;
  VEC3F kEdge0 = v1 - v0;
  VEC3F kEdge1 = v2 - v0;
  Real fA00 = norm2(kEdge0);
  Real fA01 = kEdge0 * kEdge1;
  Real fA11 = norm2(kEdge1);
  Real fB0 = kDiff*kEdge0;
  Real fB1 = kDiff*kEdge1;
  Real fC = norm2(kDiff);
  Real fDet = fabs(fA00*fA11-fA01*fA01);
  Real fS = fA01*fB1-fA11*fB0;
  Real fT = fA01*fB0-fA00*fB1;
  Real fSqrDistance;

  if (fS + fT <= fDet) {
      if (fS < (Real)0.0) {
          if (fT < (Real)0.0)  // region 4
          {
              if (fB0 < (Real)0.0) {
                  fT = (Real)0.0;
                  if (-fB0 >= fA00) {
                      fS = (Real)1.0;
                      fSqrDistance = fA00+((Real)2.0)*fB0+fC;
                  }
                  else {
                      fS = -fB0/fA00;
                      fSqrDistance = fB0*fS+fC;
                  }
              }
              else {
                  fS = (Real)0.0;
                  if (fB1 >= (Real)0.0) {
                      fT = (Real)0.0;
                      fSqrDistance = fC;
                  }
                  else if (-fB1 >= fA11) {
                      fT = (Real)1.0;
                      fSqrDistance = fA11+((Real)2.0)*fB1+fC;
                  }
                  else {
                      fT = -fB1/fA11;
                      fSqrDistance = fB1*fT+fC;
                  }
              }
          }
          else  // region 3
          {
              fS = (Real)0.0;
              if (fB1 >= (Real)0.0) {
                  fT = (Real)0.0;
                  fSqrDistance = fC;
              }
              else if (-fB1 >= fA11) {
                  fT = (Real)1.0;
                  fSqrDistance = fA11+((Real)2.0)*fB1+fC;
              }
              else {
                  fT = -fB1/fA11;
                  fSqrDistance = fB1*fT+fC;
              }
          }
      }
      else if (fT < (Real)0.0)  // region 5
      {
          fT = (Real)0.0;
          if (fB0 >= (Real)0.0) {
              fS = (Real)0.0;
              fSqrDistance = fC;
          }
          else if (-fB0 >= fA00) {
              fS = (Real)1.0;
              fSqrDistance = fA00+((Real)2.0)*fB0+fC;
          }
          else {
              fS = -fB0/fA00;
              fSqrDistance = fB0*fS+fC;
          }
      }
      else  // region 0
      {
          // minimum at interior point
          Real fInvDet = ((Real)1.0)/fDet;
          fS *= fInvDet;
          fT *= fInvDet;
          fSqrDistance = fS*(fA00*fS+fA01*fT+((Real)2.0)*fB0) +
              fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
      }
  }
  else {
      Real fTmp0, fTmp1, fNumer, fDenom;

      if (fS < (Real)0.0)  // region 2
      {
          fTmp0 = fA01 + fB0;
          fTmp1 = fA11 + fB1;
          if (fTmp1 > fTmp0) {
              fNumer = fTmp1 - fTmp0;
              fDenom = fA00-2.0f*fA01+fA11;
              if (fNumer >= fDenom) {
                  fS = (Real)1.0;
                  fT = (Real)0.0;
                  fSqrDistance = fA00+((Real)2.0)*fB0+fC;
              }
              else {
                  fS = fNumer/fDenom;
                  fT = (Real)1.0 - fS;
                  fSqrDistance = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                      fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
              }
          }
          else {
              fS = (Real)0.0;
              if (fTmp1 <= (Real)0.0) {
                  fT = (Real)1.0;
                  fSqrDistance = fA11+((Real)2.0)*fB1+fC;
              }
              else if (fB1 >= (Real)0.0) {
                  fT = (Real)0.0;
                  fSqrDistance = fC;
              }
              else {
                  fT = -fB1/fA11;
                  fSqrDistance = fB1*fT+fC;
              }
          }
      }
      else if (fT < (Real)0.0)  // region 6
      {
          fTmp0 = fA01 + fB1;
          fTmp1 = fA00 + fB0;
          if (fTmp1 > fTmp0) {
              fNumer = fTmp1 - fTmp0;
              fDenom = fA00-((Real)2.0)*fA01+fA11;
              if (fNumer >= fDenom) {
                  fT = (Real)1.0;
                  fS = (Real)0.0;
                  fSqrDistance = fA11+((Real)2.0)*fB1+fC;
              }
              else {
                  fT = fNumer/fDenom;
                  fS = (Real)1.0 - fT;
                  fSqrDistance = fS*(fA00*fS+fA01*fT+((Real)2.0)*fB0) +
                      fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
              }
          }
          else {
              fT = (Real)0.0;
              if (fTmp1 <= (Real)0.0) {
                  fS = (Real)1.0;
                  fSqrDistance = fA00+((Real)2.0)*fB0+fC;
              }
              else if (fB0 >= (Real)0.0) {
                  fS = (Real)0.0;
                  fSqrDistance = fC;
              }
              else {
                  fS = -fB0/fA00;
                  fSqrDistance = fB0*fS+fC;
              }
          }
      }
      else  // region 1
      {
          fNumer = fA11 + fB1 - fA01 - fB0;
          if (fNumer <= (Real)0.0) {
              fS = (Real)0.0;
              fT = (Real)1.0;
              fSqrDistance = fA11+((Real)2.0)*fB1+fC;
          }
          else {
              fDenom = fA00-2.0f*fA01+fA11;
              if (fNumer >= fDenom) {
                  fS = (Real)1.0;
                  fT = (Real)0.0;
                  fSqrDistance = fA00+((Real)2.0)*fB0+fC;
              }
              else {
                  fS = fNumer/fDenom;
                  fT = (Real)1.0 - fS;
                  fSqrDistance = fS*(fA00*fS+fA01*fT+((Real)2.0)*fB0) +
                      fT*(fA01*fS+fA11*fT+((Real)2.0)*fB1)+fC;
              }
          }
      }
  }

  // account for numerical round-off error
  if (fSqrDistance < (Real)0.0) {
      fSqrDistance = (Real)0.0;
  }
  return fSqrDistance;
}

//////////////////////////////////////////////////////////////////////////////
// draw a single vertex
//////////////////////////////////////////////////////////////////////////////
void OBJ::drawVert(int i)
{
#ifdef SINGLE_PRECISION
  glNormal3fv(normals[i]);
  glVertex3fv(vertices[i]);
#else
  glNormal3dv(normals[i]);
  glVertex3dv(vertices[i]);
#endif
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OBJ::transformMesh(MATRIX3 transform)
{
  // do the transform
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] = transform * vertices[x];

  // recompute all the normals
  ComputeVertexNormals();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OBJ::translateMesh(VEC3 translation)
{
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] = vertices[x] + translation;
}

//////////////////////////////////////////////////////////////////////////////
// draw the whole mesh to GL
//////////////////////////////////////////////////////////////////////////////
void OBJ::draw()
{
  if (normals.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need to call ComputeVertexNormals before calling draw! " << endl;
    return;
  }

  glBegin(GL_TRIANGLES);
  for(unsigned int i = 0; i < faces.size(); i++)
  {
    vector<int>& verts = faces[i].vertices;
    TRIANGLE triangle( &vertices[verts[0]], &vertices[verts[1]], &vertices[verts[2]]);

    triangle.draw();
  }
  glEnd();
}

//////////////////////////////////////////////////////////////////////////////
// draw the whole mesh to GL, but filter out the really stretched triangles
//////////////////////////////////////////////////////////////////////////////
void OBJ::drawFilteredSubset()
{
  if (normals.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need to call ComputeVertexNormals before calling draw! " << endl;
    return;
  }

  glBegin(GL_TRIANGLES);
  MERSENNETWISTER twister(123456);
  for(unsigned int i = 0; i < faces.size(); i++)
  {
    vector<int>& verts = faces[i].vertices;

    Real random = twister.rand();

    if (random > 0.2) continue;

    VEC3F* triangleVertices[3];
    for (int y = 0; y < 3; y++)
      triangleVertices[y] = &vertices[faces[i].vertices[y]];
    TRIANGLE triangle(triangleVertices[0], triangleVertices[1], triangleVertices[2]);
    Real currentArea = triangle.area();
    Real currentMaxLength = triangle.maxEdgeLength();

    if (currentArea > _triangleAreas[i] * _filteringThreshold)
      continue;
    if (currentMaxLength > _maxEdgeLengths[i] * _filteringThreshold)
      continue;

    if(verts.size() == 3)
    {
      drawVert(verts[0]);
      drawVert(verts[1]);
      drawVert(verts[2]);
    }
    else if(verts.size() == 4)
    {
      drawVert(verts[0]);
      drawVert(verts[1]);
      drawVert(verts[2]);

      drawVert(verts[2]);
      drawVert(verts[3]);
      drawVert(verts[0]);
    }
  }
  glEnd();
}

//////////////////////////////////////////////////////////////////////////////
// draw the whole mesh to GL, but filter out the really stretched triangles
//////////////////////////////////////////////////////////////////////////////
void OBJ::drawFiltered()
{
  if (normals.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Need to call ComputeVertexNormals before calling draw! " << endl;
    return;
  }

  glBegin(GL_TRIANGLES);
  for(unsigned int i = 0; i < faces.size(); i++)
  {
    vector<int>& verts = faces[i].vertices;

    VEC3F* triangleVertices[3];
    for (int y = 0; y < 3; y++)
      triangleVertices[y] = &vertices[faces[i].vertices[y]];
    TRIANGLE triangle(triangleVertices[0], triangleVertices[1], triangleVertices[2]);
    Real currentArea = triangle.area();
    Real currentMaxLength = triangle.maxEdgeLength();

    if (currentArea > _triangleAreas[i] * _filteringThreshold)
      continue;
    if (currentMaxLength > _maxEdgeLengths[i] * _filteringThreshold)
      continue;

    if(verts.size() == 3)
    {
      drawVert(verts[0]);
      drawVert(verts[1]);
      drawVert(verts[2]);
    }
    else if(verts.size() == 4)
    {
      drawVert(verts[0]);
      drawVert(verts[1]);
      drawVert(verts[2]);

      drawVert(verts[2]);
      drawVert(verts[3]);
      drawVert(verts[0]);
    }
  }
  glEnd();
}

//////////////////////////////////////////////////////////////////////////////
// load in a TET_MESH surface mesh, just to fake out the AOV renderer
//////////////////////////////////////////////////////////////////////////////
void OBJ::LoadTetSurfaceMesh(vector<TRIANGLE*> triangles)
{
	vertices.resize (0);
	normals.resize (0);
	texcoords.resize (0);

  // create a hash of all the vertices in the surface mesh
  map<VEC3F*, int> vertexHash;
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    TRIANGLE& triangle = *triangles[x];

    vertexHash[triangle.vertex(0)] = 1;
    vertexHash[triangle.vertex(1)] = 1;
    vertexHash[triangle.vertex(2)] = 1;
  }

  // push all the vertices onto the vector
  map<VEC3F*, int>::iterator iter;
  for (iter = vertexHash.begin(); iter != vertexHash.end(); iter++)
  {
    VEC3F* vertex = iter->first;
    vertices.push_back(*vertex);
    iter->second = vertices.size() - 1;
  }

  // create the faces
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    Face face;
    TRIANGLE& triangle = *triangles[x];

    face.vertices.push_back(vertexHash[triangle.vertex(0)]);
    face.vertices.push_back(vertexHash[triangle.vertex(1)]);
    face.vertices.push_back(vertexHash[triangle.vertex(2)]);

    faces.push_back(face);
  }

  // populate the vertex normals
  ComputeVertexNormals();
}

//////////////////////////////////////////////////////////////////////////////
// refresh a TET_MESH surface mesh, just to fake out the AOV renderer
//////////////////////////////////////////////////////////////////////////////
void OBJ::RefreshTetSurfaceMesh(vector<TRIANGLE*> triangles)
{
  // create a hash of all the vertices in the surface mesh
  map<VEC3F*, int> vertexHash;
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    TRIANGLE& triangle = *triangles[x];

    vertexHash[triangle.vertex(0)] = 1;
    vertexHash[triangle.vertex(1)] = 1;
    vertexHash[triangle.vertex(2)] = 1;
  }

  // push all the vertices onto the vector
  map<VEC3F*, int>::iterator iter;
  int counter = 0;
  for (iter = vertexHash.begin(); iter != vertexHash.end(); iter++, counter++)
    vertices[counter] = *(iter->first);

  // refresh the vertex normals
  ComputeVertexNormals();
}
