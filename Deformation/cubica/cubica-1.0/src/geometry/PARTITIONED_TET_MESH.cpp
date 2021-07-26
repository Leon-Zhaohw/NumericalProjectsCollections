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
// TET_MESH.h: interface for the TET_MESH class.
//
//////////////////////////////////////////////////////////////////////

#include "PARTITIONED_TET_MESH.h"
#include <SUBMATRIX.h>
#include <fstream>
#include <sstream>
#include <MERSENNETWISTER.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
PARTITIONED_TET_MESH::PARTITIONED_TET_MESH(const char* filename, MATERIAL** materials, int totalMaterials, int partitions, Real springConst, bool simulate, string partitionPath) :
  _partitions(partitions),
  _filename(filename),
  _partitionPath(partitionPath),
  _clonedVertices(NULL),
  _clonedVertexMap(NULL),
  _clonedSurfaceVertices(NULL),
  _clonedTriangles(NULL),
  _clonedTriangleIDs(NULL),
  _interfaceArea(NULL),
  _graphColors(NULL)
{
  // read in the original mesh
  _originalMesh = new TET_MESH(filename, materials, totalMaterials, simulate);

  // allocate per-interface spring consts
  _springConst = springConst;
  _interfaceSpringConst.resizeAndWipe(_partitions, _partitions);
  for (int y = 0; y < _partitions; y++)
    for (int x = 0; x < _partitions; x++)
      _interfaceSpringConst(x,y) = _springConst;

  // read in each separate partition
  _meshes = new TET_MESH*[_partitions];
  for (int x = 0; x < _partitions; x++)
  {
    char buffer[256];
    sprintf(buffer, "%i", x);
    string partitionFilename = partitionPath;
    if (partitionFilename.length() == 0)
      partitionFilename = string(filename);
    partitionFilename += string(".partition.");
    partitionFilename += string(buffer);

    // probe the file to see if it is unconstrained
    FILE* probe = fopen(partitionFilename.c_str(), "rb");
    if (probe == NULL) printf("Filename %s not found!\n", partitionFilename.c_str());
    int unconstrainedSize, constrainedSize;
    fread((void*)&unconstrainedSize, sizeof(int), 1, probe);
    fread((void*)&constrainedSize, sizeof(int), 1, probe);
    fclose(probe);

    if (constrainedSize == 0)
      cout << " Partition " << x << " is unconstrained." << endl;
    
    cout << " Reading partition file: " << partitionFilename.c_str() << endl;
    if (constrainedSize == 0)
    {
      //_meshes[x] = new UNCONSTRAINED_TET_MESH(partitionFilename.c_str(), materials, totalMaterials, simulate);
      _meshes[x] = new TET_MESH(partitionFilename.c_str(), materials, totalMaterials, simulate);
      _unconstrainedPartition.push_back(true);
    }
    else
    {
      _meshes[x] = new TET_MESH(partitionFilename.c_str(), materials, totalMaterials, simulate);
      _unconstrainedPartition.push_back(false);
    }
    cout << endl;
  }
  
  // allocate the clone table
  _clonedVertices = new vector<pair<int, int> >*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _clonedVertices[x] = new vector<pair<int, int> >[_partitions];
  
  // read in the cloned vertices
  _originalIDs = NULL;
  string cloneFilename = partitionPath;
  if (cloneFilename.length() == 0)
    cloneFilename = string(filename);
  cloneFilename += string(".clones");
  cout << " Reading clone table: " << cloneFilename.c_str() << endl;
  FILE* file = fopen(cloneFilename.c_str(), "rb");
  if (file == NULL)
  {
    cout << " No clone table " << cloneFilename.c_str() << " found! " << endl;
    exit(0);
  }
  for (int x = 0; x < _partitions; x++)
  {
    for (int y = 0; y < _partitions; y++)
    {
      // read in how many interactions there are
      int totalClones;
      fread((void*)&totalClones, sizeof(int), 1, file);

      // read in the pairs
      for (int z = 0; z < totalClones; z++)
      {
        int first, second;
        fread((void*)&first, sizeof(int), 1, file);
        fread((void*)&second, sizeof(int), 1, file);
        pair<int, int> clone(first, second);

        _clonedVertices[x][y].push_back(clone);

        VEC3F* leftClone = _meshes[x]->vertices(first);
        VEC3F* rightClone = _meshes[y]->vertices(second);

        _isCloned[leftClone]++;
        _isCloned[rightClone]++;
      }
    }
  }
  fclose(file);

  // sanity check
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
      assert(_clonedVertices[x][y].size() == _clonedVertices[y][x].size());

  // populate the explicit vertex version
  computeClonedVertexMap();

  // read in the correspondences
  cout << " Reading in correspondences ... ";
  _originalIDs = new vector<int>[_partitions];
  for (int x = 0; x < _partitions; x++)
  {
    char buffer[256];
    sprintf(buffer, "%i", x);
    string partitionFilename = partitionPath;
    if (partitionFilename.length() == 0)
      partitionFilename = string(filename);
    partitionFilename += string(".partition.");
    partitionFilename += string(buffer);
    partitionFilename += string(".correspondence");
    FILE* file = fopen(partitionFilename.c_str(), "rb");

    for (int y = 0; y < _meshes[x]->totalNodes(); y++)
    {
      int oldID;
      fread((void*)&oldID, sizeof(int), 1, file);
      _originalIDs[x].push_back(oldID);
    }

    fclose(file);
  }
  cout << " done. " << endl;

  // build inverse correspondences
  cout << " Building inverse correspondences ... "; flush(cout);
  for (int x = 0; x < _partitions; x++)
    for (unsigned int y = 0; y < _originalIDs[x].size(); y++)
    {
      int oldID = _originalIDs[x][y];
      _partitionIDs.insert(pair<int, pair<int, int> >(oldID, pair<int, int>(x,y)));
    }
  cout << " done." << endl;

  // compute the center of mass so we can draw the exploded view
  computeCenterOfMass();

  // allocate the array of graph colors;
  _graphColors = new int[_partitions];
  
  // read in the graph coloring, and if there isn't one, compute one
  string dimacs = partitionPath;
  if (dimacs.length() == 0)
    dimacs = string(filename);
  dimacs += string(".dimacs.res");
  cout << " Reading in graph coloring ... ";
  if (!readDIMACS(dimacs.c_str()))
  {
    cout << " No graph coloring found, computing one ... ";
    computeGraphColoring(partitionPath);
  }
  cout << "done." << endl;

  /*
  // compute the cloned triangles
  computeClonedTriangles();
  computeClonedSurfaceVertices();
  computeBlendedSurfaceMesh();
 
  // compute surface areas of cloned triangles
  computeInterfaceAreas();
  */
 
  // cache the colors - RGB MCYW
  /*
  _colors[0][0] = 1.0f; _colors[0][1] = 0.0f; _colors[0][2] = 0.0f; _colors[0][3] = 1.0f;
  _colors[1][0] = 0.0f; _colors[1][1] = 1.0f; _colors[1][2] = 0.0f; _colors[1][3] = 1.0f;
  _colors[2][0] = 0.0f; _colors[2][1] = 0.0f; _colors[2][2] = 1.0f; _colors[2][3] = 1.0f;
  _colors[3][0] = 1.0f; _colors[3][1] = 0.0f; _colors[3][2] = 1.0f; _colors[3][3] = 1.0f;
  //_colors[4][0] = 0.0f; _colors[4][1] = 1.0f; _colors[4][2] = 1.0f; _colors[4][3] = 1.0f;
  _colors[4][0] = 0.5f; _colors[4][1] = 0.5f; _colors[4][2] = 0.5f; _colors[5][3] = 1.0f;
  _colors[5][0] = 1.0f; _colors[5][1] = 1.0f; _colors[5][2] = 0.0f; _colors[5][3] = 1.0f;
  _colors[6][0] = 1.0f; _colors[6][1] = 1.0f; _colors[6][2] = 1.0f; _colors[6][3] = 1.0f;
  _colors[7][0] = 0.5f; _colors[7][1] = 0.5f; _colors[7][2] = 0.5f; _colors[7][3] = 1.0f;
  */

  _colors[0][0] = 1.0f; _colors[0][1] = 0.7f; _colors[0][2] = 0.7f; _colors[0][3] = 1.0f;
  _colors[1][0] = 0.7f; _colors[1][1] = 1.0f; _colors[1][2] = 0.7f; _colors[1][3] = 1.0f;
  _colors[2][0] = 0.7f; _colors[2][1] = 0.7f; _colors[2][2] = 1.0f; _colors[2][3] = 1.0f;
  _colors[3][0] = 1.0f; _colors[3][1] = 0.7f; _colors[3][2] = 1.0f; _colors[3][3] = 1.0f;
  _colors[4][0] = 0.7f; _colors[4][1] = 1.0f; _colors[4][2] = 1.0f; _colors[4][3] = 1.0f;
  _colors[5][0] = 1.0f; _colors[5][1] = 1.0f; _colors[5][2] = 0.7f; _colors[5][3] = 1.0f;
  _colors[6][0] = 1.0f; _colors[6][1] = 1.0f; _colors[6][2] = 1.0f; _colors[6][3] = 1.0f;
  _colors[7][0] = 0.5f; _colors[7][1] = 0.5f; _colors[7][2] = 0.5f; _colors[7][3] = 1.0f;

  // recompute the mass matrices so that they all add to unity
  if (simulate)
    for (int x = 0; x < _partitions; x++)
    {
      Real trueMass = _originalMesh->mass(0);
      _meshes[x]->TET_MESH::resetMasses(trueMass);

      // recompute center of mass matrix with new masses
      _meshes[x]->computeCenterOfMass();
      _meshes[x]->computeRestCenterOfMass();
    }

  /*
  srand(123456);
  for (int x = 0; x < 8; x++)
  {
    for (int y = 0; y < 3; y++)
    {
      _colors[x][y] = (Real)rand() / RAND_MAX;
    }
    _colors[x][3] = 1.0;
  }
  */

  int totalUnconstrained = 0;
  for (unsigned int x = 0; x < _unconstrainedPartition.size(); x++)
    if (_unconstrainedPartition[x])
      totalUnconstrained++;

  cout << " Unconstrained partitions: " << totalUnconstrained << " of " << _partitions << endl;

  // evenly divide mass
  recomputePartitionedMass();

  computeClonedTriangles();
  computeClonedSurfaceVertices();
  computeInterfaceAreas();
}

PARTITIONED_TET_MESH::PARTITIONED_TET_MESH() :
  _clonedVertices(NULL),
  _clonedSurfaceVertices(NULL),
  _clonedTriangles(NULL),
  _clonedTriangleIDs(NULL),
  _interfaceArea(NULL)
{
}

PARTITIONED_TET_MESH::~PARTITIONED_TET_MESH()
{
  delete _originalMesh;
  delete[] _graphColors;
  
  for (int x = 0; x < _partitions; x++)
    delete _meshes[x];
  delete[] _meshes;

  // clean up clones
  if (_clonedVertices != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedVertices[x];
    delete[] _clonedVertices;
  }

  if (_clonedVertexMap != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedVertexMap[x];
    delete[] _clonedVertexMap;
  }

  // clean up clones
  if (_clonedSurfaceVertices != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedSurfaceVertices[x];
    delete[] _clonedSurfaceVertices;
  }

  if (_clonedTriangles != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedTriangles[x];
    delete[] _clonedTriangles;
  }

  if (_clonedTriangleIDs != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedTriangleIDs[x];
    delete[] _clonedTriangleIDs;
  }

  if (_interfaceArea != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _interfaceArea[x];
    delete[] _interfaceArea;
  }

  if (_blendedTriangles.size() > 0)
    for (unsigned int x = 0; x < _blendedTriangles.size(); x++)
      delete _blendedTriangles[x];

  if (_blendedVertices.size() > 0)
    for (unsigned int x = 0; x < _blendedVertices.size(); x++)
      delete _blendedVertices[x];
}

//////////////////////////////////////////////////////////////////////
// Draw all the partitions, with one highlighted
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawHighlightedPartition(int highlight)
{
#ifdef USING_GLVU
  for (int x = 0; x < _partitions; x++)
  {
    _meshes[x]->updateSurfaceMesh();
    /*
    if (highlight == x)
      glColor4f(1.0, 1.0, 1.0, 1.0);
    else
      glColor4f(0.1, 0.1, 0.1, 0.1);
      */
    if (highlight != x) continue;

    glColor4f(1.0, 1.0, 1.0, 1.0);
    _meshes[x]->drawSurfaceFaces();
    /*
    glDisable(GL_DEPTH_TEST);
    glColor4f(0.0f, 0.0f, 1.0f, 10.0f);
    glPointSize(10.0f);
    _meshes[x]->drawUnconstrainedNodes();
    glEnable(GL_DEPTH_TEST);
    */
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the partitions different colors
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawSurfaceFaces(bool shadowCall)
{
#ifdef USING_GLVU
  // if there's more than 8 partitions just do darker versions
  // of the colors we already have
  for (int x = 0; x < _partitions; x++)
  {
    float color[4];
    //int mod = _graphColors[x] % 8;
    int mod = x % 8;
    //int div = _graphColors[x] / 8;
    //float fraction = 1.0f / (div + 1.0f);
    float fraction = 1.0f;
    color[0] = fraction * _colors[mod][0];
    color[1] = fraction * _colors[mod][1];
    color[2] = fraction * _colors[mod][2];
    color[3] = 1.0f;
    /*
    color[0] = twister.rand();
    color[1] = twister.rand();
    color[2] = twister.rand();
    color[3] = 1.0f;
    */

    // make it a pastel
    float clamp = 0.5;
    //float clamp = 0.25;
    color[0] = (color[0] < clamp) ? clamp : color[0];
    color[1] = (color[1] < clamp) ? clamp : color[1];
    color[2] = (color[2] < clamp) ? clamp : color[2];

    //color[3] = 0.5f;
    if (!shadowCall)
      glColor4f(color[0], color[1], color[2], color[3]);
    _meshes[x]->drawSurfaceFaces();

    /*
    glDisable(GL_DEPTH_TEST);
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
    glPointSize(10.0f);
    _meshes[x]->drawConstrainedNodes();
    glEnable(GL_DEPTH_TEST);
    */

    //glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
    //_meshes[x]->drawUnconstrainedNodes();
  }
  /*
  glEnable(GL_BLEND);
  glEnable(GL_DEPTH);
  */
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawAllTets()
{
  for (int x = 0; x < _partitions; x++)
    _meshes[x]->drawAllTets();
}

//////////////////////////////////////////////////////////////////////
// Draw just the cloned triangles along the interfaces
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawClonedTriangles()
{
#ifdef USING_GLVU  
  // if there's more than 8 partitions just do darker versions
  // of the colors we already have
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      for (unsigned int z = 0; z < _clonedTriangles[x][y].size(); z++)
      {
        float color[4];
        int mod = _graphColors[x] % 8;
        int div = _graphColors[x] / 8;
        float fraction = 1.0f / (div + 1.0f);
        color[0] = fraction * _colors[mod][0];
        color[1] = fraction * _colors[mod][1];
        color[2] = fraction * _colors[mod][2];
        color[3] = 1.0f;
        glColor4f(color[0], color[1], color[2], color[3]);
        _clonedTriangles[x][y][z]->draw();
      }
    }
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw just the cloned triangles along the interfaces
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawClonedTriangles(int partition)
{
#ifdef USING_GLVU  
  // if there's more than 8 partitions just do darker versions
  // of the colors we already have
  for (int y = 0; y < _partitions; y++)
  {
    float color[4];
    int mod = _graphColors[y] % 8;
    int div = _graphColors[y] / 8;
    float fraction = 1.0f / (div + 1.0f);
    color[0] = fraction * _colors[mod][0];
    color[1] = fraction * _colors[mod][1];
    color[2] = fraction * _colors[mod][2];
    color[3] = 1.0f;
    glColor4f(color[0], color[1], color[2], color[3]);

    for (unsigned int z = 0; z < _clonedTriangles[partition][y].size(); z++)
      _clonedTriangles[partition][y][z]->draw();
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw just the cloned points along the interfaces
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawClonedVertices(int partition)
{
#ifdef USING_GLVU  
  glPointSize(4.0);
  glBegin(GL_POINTS);
  for (int x = 0; x < _partitions; x++)
  {
    vector<pair<int,int> >& clonedSurfaceVertices = _clonedVertices[partition][x];
    glColor4f(10.0, 10.0, 10.0, 10.0);
    for (unsigned int y = 0; y < clonedSurfaceVertices.size(); y++)
    {
      int index = clonedSurfaceVertices[y].first;
      VEC3F vertex = *(_meshes[partition]->vertices(index));
      glVertex3f(vertex[0],vertex[1],vertex[2]);
    }
  }
  glEnd();
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw just the cloned vertices along the surface and interface
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawClonedSurfaceVertices(int partition)
{
#ifdef USING_GLVU  
  static bool first = true;
 
  glPointSize(4.0);
  glBegin(GL_POINTS);
  for (int x = 0; x < _partitions; x++)
  {
    vector<pair<int,int> >& clonedSurfaceVertices = _clonedSurfaceVertices[partition][x];
    glColor4f(10.0, 10.0, 10.0, 10.0);
    for (unsigned int y = 0; y < clonedSurfaceVertices.size(); y++)
    {
      int index = clonedSurfaceVertices[y].first;
      VEC3F vertex = *(_meshes[partition]->vertices(index));
      if (first)
        cout << " surface: " << vertex << endl;
      glVertex3f(vertex[0],vertex[1],vertex[2]);
    }
  }
  glEnd();

  first = false;
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the partitions different colors
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawSurfaceFaces(int partition)
{
#ifdef USING_GLVU  
  float color[4];
  int mod = partition % 8;
  int div = partition / 8;
  float fraction = 1.0f / (div + 1.0f);
  color[0] = fraction * _colors[mod][0];
  color[1] = fraction * _colors[mod][1];
  color[2] = fraction * _colors[mod][2];
  color[3] = 1.0f;
  glColor4f(color[0], color[1], color[2], color[3]);
  _meshes[partition]->drawSurfaceFaces();
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
  //glPointSize(10.0f);
  //_meshes[partition]->drawConstrainedNodes();
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw a specific cloned vertex
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawClone(int partition0, int partition1, int clone)
{
#ifdef USING_GLVU  
  vector<pair<int, int> >& clones = _clonedVertices[partition0][partition1];
  if ((unsigned int)clone > clones.size()) return;

  int index0 = clones[clone].first;
  int index1 = clones[clone].second;

  VEC3F vertex0 = *(_meshes[partition0]->vertices(index0));
  VEC3F vertex1 = *(_meshes[partition1]->vertices(index1));

  glBegin(GL_POINTS);
    glVertex3f(vertex0[0], vertex0[1], vertex0[2]);
    glVertex3f(vertex1[0], vertex1[1], vertex1[2]);
  glEnd();
  glBegin(GL_LINES);
    glVertex3f(vertex0[0], vertex0[1], vertex0[2]);
    glVertex3f(vertex1[0], vertex1[1], vertex1[2]);
  glEnd();
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the partitions exploded with respect to the center of mass
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawExplodedCentered(int center)
{
#ifdef USING_GLVU  
  // if there's more than 8 partitions just do darker versions
  // of the colors we already have
  for (int x = 0; x < _partitions; x++)
  {
    if (x != center && !neighbors(center, x)) continue; 

    //float color[4];
    //int mod = _graphColors[x] % 8;
    //int div = _graphColors[x] / 8;
    //float fraction = 1.0f / (div + 1.0f);
    if (x == center)
      glColor4f(1.0, 0.0, 0.0, 1.0);
    else
      glColor4f(1.0, 1.0, 1.0, 1.0);
    glPushMatrix();
      VEC3F explode = _meshes[x]->centerOfMass() - _centerOfMass;
      glTranslatef(explode[0], explode[1], explode[2]);
      _meshes[x]->drawSurfaceFaces();
      glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
      glPointSize(10.0f);
      _meshes[x]->drawConstrainedNodes();
    glPopMatrix();
  }

  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
    {
      if (x != center && y != center) continue;
      if (!neighbors(x,y)) continue;

      for (unsigned int z = 0; z < _clonedVertices[x][y].size(); z++)
      {
        int first = _clonedVertices[x][y][z].first;
        int second = _clonedVertices[x][y][z].second;

        VEC3F firstVertex = *(_meshes[x]->vertices(first));
        VEC3F secondVertex = *(_meshes[y]->vertices(second));

        if (_unconstrainedPartition[x])
          firstVertex += ((UNCONSTRAINED_TET_MESH*)_meshes[x])->rigidTranslation();
        if (_unconstrainedPartition[y])
          secondVertex += ((UNCONSTRAINED_TET_MESH*)_meshes[y])->rigidTranslation();

        VEC3F firstExplode = _meshes[x]->centerOfMass() - _centerOfMass;
        VEC3F secondExplode = _meshes[y]->centerOfMass() - _centerOfMass;

        firstVertex += firstExplode;
        secondVertex += secondExplode;

        glBegin(GL_LINES);
          float color[4];
          int mod = x % 8;
          int div = x / 8;
          float fraction = 1.0f / (div + 1.0f);
          color[0] = fraction * _colors[mod][0];
          color[1] = fraction * _colors[mod][1];
          color[2] = fraction * _colors[mod][2];
          color[3] = 1.0f;
          glColor4f(color[0], color[1], color[2], color[3]);
          glVertex3f(firstVertex[0], firstVertex[1], firstVertex[2]);
          mod = y % 8;
          div = y / 8;
          fraction = 1.0f / (div + 1.0f);
          color[0] = fraction * _colors[mod][0];
          color[1] = fraction * _colors[mod][1];
          color[2] = fraction * _colors[mod][2];
          glColor4f(color[0], color[1], color[2], color[3]);
          glVertex3f(secondVertex[0], secondVertex[1], secondVertex[2]);
        glEnd();
      }
    }
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the partitions exploded with respect to the center of mass
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawExploded(int left, int right)
{
  if (left > right)
  {
    int swap = left;
    left = right;
    right = swap;
  }

#ifdef USING_GLVU  
  // if there's more than 8 partitions just do darker versions
  // of the colors we already have
  for (int x = 0; x < _partitions; x++)
  {
    if (x != left && x != right) continue; 

    float color[4];
    int mod = _graphColors[x] % 8;
    int div = _graphColors[x] / 8;
    float fraction = 1.0f / (div + 1.0f);
    color[0] = fraction * _colors[mod][0];
    color[1] = fraction * _colors[mod][1];
    color[2] = fraction * _colors[mod][2];
    color[3] = 1.0f;
    glColor4f(color[0], color[1], color[2], color[3]);
    glPushMatrix();
      VEC3F explode = _meshes[x]->centerOfMass() - _centerOfMass;
      //explode *= 0.5;
      glTranslatef(explode[0], explode[1], explode[2]);
      _meshes[x]->drawSurfaceFaces();
      glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
      glPointSize(10.0f);
      _meshes[x]->drawConstrainedNodes();
    glPopMatrix();
  }

  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
    {
      if (x != left || y != right) continue;

      for (unsigned int z = 0; z < _clonedVertices[x][y].size(); z++)
      {
        int first = _clonedVertices[x][y][z].first;
        int second = _clonedVertices[x][y][z].second;

        VEC3F firstVertex = *(_meshes[x]->vertices(first));
        VEC3F secondVertex = *(_meshes[y]->vertices(second));

        if (_unconstrainedPartition[x])
          firstVertex += ((UNCONSTRAINED_TET_MESH*)_meshes[x])->rigidTranslation();
        if (_unconstrainedPartition[y])
          secondVertex += ((UNCONSTRAINED_TET_MESH*)_meshes[y])->rigidTranslation();

        VEC3F firstExplode = _meshes[x]->centerOfMass() - _centerOfMass;
        VEC3F secondExplode = _meshes[y]->centerOfMass() - _centerOfMass;

        firstVertex += firstExplode;
        secondVertex += secondExplode;

        glBegin(GL_LINES);
          float color[4];
          int mod = x % 8;
          int div = x / 8;
          float fraction = 1.0f / (div + 1.0f);
          color[0] = fraction * _colors[mod][0];
          color[1] = fraction * _colors[mod][1];
          color[2] = fraction * _colors[mod][2];
          color[3] = 1.0f;
          glColor4f(color[0], color[1], color[2], color[3]);
          glVertex3f(firstVertex[0], firstVertex[1], firstVertex[2]);
          mod = y % 8;
          div = y / 8;
          fraction = 1.0f / (div + 1.0f);
          color[0] = fraction * _colors[mod][0];
          color[1] = fraction * _colors[mod][1];
          color[2] = fraction * _colors[mod][2];
          glColor4f(color[0], color[1], color[2], color[3]);
          glVertex3f(secondVertex[0], secondVertex[1], secondVertex[2]);
        glEnd();
      }
    }
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the cloned vertices and faces only
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawClonesOnly(int left, int right)
{
#ifdef USING_GLVU
  VEC3F firstExplode = _meshes[left]->centerOfMass() - _centerOfMass;
  VEC3F secondExplode = _meshes[right]->centerOfMass() - _centerOfMass;

  float color[4];
  int mod = _graphColors[left] % 8;
  int div = _graphColors[left] / 8;
  float fraction = 1.0f / (div + 1.0f);
  color[0] = fraction * _colors[mod][0];
  color[1] = fraction * _colors[mod][1];
  color[2] = fraction * _colors[mod][2];
  color[3] = 1.0f;
  glColor4f(color[0], color[1], color[2], color[3]);

  glPushMatrix();
  glTranslatef(firstExplode[0], firstExplode[1], firstExplode[2]);
  for (unsigned int z = 0; z < _clonedTriangles[left][right].size(); z++)
    _clonedTriangles[left][right][z]->draw();
  glPopMatrix();

  glPushMatrix();
  glTranslatef(secondExplode[0], secondExplode[1], secondExplode[2]);
  for (unsigned int z = 0; z < _clonedTriangles[right][left].size(); z++)
    _clonedTriangles[right][left][z]->draw();
  glPopMatrix();

  int x = left;
  int y = right;
  for (unsigned int z = 0; z < _clonedVertices[x][y].size(); z++)
  {
    int first = _clonedVertices[x][y][z].first;
    int second = _clonedVertices[x][y][z].second;

    VEC3F firstVertex = *(_meshes[x]->vertices(first));
    VEC3F secondVertex = *(_meshes[y]->vertices(second));

    /*
    if (_unconstrainedPartition[x])
      firstVertex += ((UNCONSTRAINED_TET_MESH*)_meshes[x])->centerOfMass();
    if (_unconstrainedPartition[y])
      secondVertex += ((UNCONSTRAINED_TET_MESH*)_meshes[y])->centerOfMass();
      */

    VEC3F firstExplode = _meshes[x]->centerOfMass() - _centerOfMass;
    VEC3F secondExplode = _meshes[y]->centerOfMass() - _centerOfMass;

    firstVertex += firstExplode;
    secondVertex += secondExplode;

    glBegin(GL_LINES);
      float color[4];
      int mod = x % 8;
      int div = x / 8;
      float fraction = 10;
      color[0] = fraction * _colors[mod][0];
      color[1] = fraction * _colors[mod][1];
      color[2] = fraction * _colors[mod][2];
      color[3] = 10.0f;
      glColor4f(color[0], color[1], color[2], color[3]);
      glVertex3f(firstVertex[0], firstVertex[1], firstVertex[2]);
      mod = y % 8;
      div = y / 8;
      fraction = 1.0f / (div + 1.0f);
      color[0] = fraction * _colors[mod][0];
      color[1] = fraction * _colors[mod][1];
      color[2] = fraction * _colors[mod][2];
      glColor4f(color[0], color[1], color[2], color[3]);
      glVertex3f(secondVertex[0], secondVertex[1], secondVertex[2]);
    glEnd();
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the partitions exploded with respect to the center of mass
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawSprings(int partition)
{
#ifdef USING_GLVU  
  glLineWidth(10.0f);
  MATRIX3 Ri = quaternionRotation(partition).toExplicitMatrix3x3();
  VEC3F ti = rigidTranslation(partition);
  for (int y = 0; y < _partitions; y++)
  {
    MATRIX3 Rj = quaternionRotation(y).toExplicitMatrix3x3();
    VEC3F tj = rigidTranslation(y);
    for (unsigned int z = 0; z < _clonedVertices[partition][y].size(); z++)
    {
      int first = _clonedVertices[partition][y][z].first;
      int second = _clonedVertices[partition][y][z].second;

      VEC3F firstVertex = *(_meshes[partition]->vertices(first));
      VEC3F secondVertex = *(_meshes[y]->vertices(second));

      firstVertex = Ri * firstVertex + ti;
      secondVertex = Rj * secondVertex + tj;

      glBegin(GL_LINES);
        float color[4];
        int mod = partition % 8;
        int div = partition / 8;
        float fraction = 10.0f / (div + 1.0f);
        color[0] = fraction * _colors[mod][0];
        color[1] = fraction * _colors[mod][1];
        color[2] = fraction * _colors[mod][2];
        color[3] = 1.0f;
        glColor4f(color[0], color[1], color[2], color[3]);
        glVertex3f(firstVertex[0], firstVertex[1], firstVertex[2]);
        mod = y % 8;
        div = y / 8;
        fraction = 10.0f / (div + 1.0f);
        color[0] = fraction * _colors[mod][0];
        color[1] = fraction * _colors[mod][1];
        color[2] = fraction * _colors[mod][2];
        glColor4f(color[0], color[1], color[2], color[3]);
        glVertex3f(secondVertex[0], secondVertex[1], secondVertex[2]);
      glEnd();
    }
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the partitions exploded with respect to the center of mass
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawSprings()
{
#ifdef USING_GLVU  
  glLineWidth(4.0f);
  for (int x = 0; x < _partitions; x++)
  {
    MATRIX3 Ri = quaternionRotation(x).toExplicitMatrix3x3();
    VEC3F ti = rigidTranslation(x);
    for (int y = x + 1; y < _partitions; y++)
    {
      MATRIX3 Rj = quaternionRotation(y).toExplicitMatrix3x3();
      VEC3F tj = rigidTranslation(y);
      for (unsigned int z = 0; z < _clonedVertices[x][y].size(); z++)
      {
        int first = _clonedVertices[x][y][z].first;
        int second = _clonedVertices[x][y][z].second;

        VEC3F firstVertex = *(_meshes[x]->vertices(first));
        VEC3F secondVertex = *(_meshes[y]->vertices(second));

        firstVertex = Ri * firstVertex + ti;
        secondVertex = Rj * secondVertex + tj;

        glBegin(GL_LINES);
          float color[4];
          int mod = x % 8;
          int div = x / 8;
          float fraction = 1.0f / (div + 1.0f);
          color[0] = fraction * _colors[mod][0];
          color[1] = fraction * _colors[mod][1];
          color[2] = fraction * _colors[mod][2];
          color[3] = 1.0f;
          //glColor4f(color[0], color[1], color[2], color[3]);
          glVertex3f(firstVertex[0], firstVertex[1], firstVertex[2]);
          mod = y % 8;
          div = y / 8;
          fraction = 1.0f / (div + 1.0f);
          color[0] = fraction * _colors[mod][0];
          color[1] = fraction * _colors[mod][1];
          color[2] = fraction * _colors[mod][2];
          //glColor4f(color[0], color[1], color[2], color[3]);
          glVertex3f(secondVertex[0], secondVertex[1], secondVertex[2]);
        glEnd();
      }
    }
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the partitions exploded with respect to the center of mass
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawExploded()
{
#ifdef USING_GLVU  
  // if there's more than 8 partitions just do darker versions
  // of the colors we already have
  for (int x = 0; x < _partitions; x++)
  {
    float color[4];
    int mod = _graphColors[x] % 8;
    int div = _graphColors[x] / 8;
    float fraction = 1.0f / (div + 1.0f);
    color[0] = fraction * _colors[mod][0];
    color[1] = fraction * _colors[mod][1];
    color[2] = fraction * _colors[mod][2];
    color[3] = 1.0f;
    glColor4f(color[0], color[1], color[2], color[3]);
    //glColor4f(1,1,1,1);
    glPushMatrix();
      VEC3F explode = _meshes[x]->centerOfMass() - _centerOfMass;
      //explode *= 0.5;
      glTranslatef(explode[0], explode[1], explode[2]);
      _meshes[x]->drawSurfaceFaces();
      glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
      glPointSize(10.0f);
      _meshes[x]->drawConstrainedNodes();
    glPopMatrix();
  }

  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
    {
      //cout << " clones " << x << ", " << y << ": " << _clonedVertices[x][y].size() << endl;
      for (unsigned int z = 0; z < _clonedVertices[x][y].size(); z++)
      {
        int first = _clonedVertices[x][y][z].first;
        int second = _clonedVertices[x][y][z].second;

        VEC3F firstVertex = *(_meshes[x]->vertices(first));
        VEC3F secondVertex = *(_meshes[y]->vertices(second));

        /*
        if (_unconstrainedPartition[x])
          firstVertex += ((UNCONSTRAINED_TET_MESH*)_meshes[x])->centerOfMass();
        if (_unconstrainedPartition[y])
          secondVertex += ((UNCONSTRAINED_TET_MESH*)_meshes[y])->centerOfMass();
        */

        VEC3F firstExplode = _meshes[x]->centerOfMass() - _centerOfMass;
        VEC3F secondExplode = _meshes[y]->centerOfMass() - _centerOfMass;

        firstVertex += firstExplode;
        secondVertex += secondExplode;

        VEC3F leftTranslation = rigidTranslation(x); 
        VEC3F rightTranslation = rigidTranslation(y);

        firstVertex += leftTranslation;
        secondVertex += rightTranslation; 

        glBegin(GL_LINES);
          float color[4];
          int mod = x % 8;
          int div = x / 8;
          float fraction = 1.0f / (div + 1.0f);
          color[0] = fraction * _colors[mod][0];
          color[1] = fraction * _colors[mod][1];
          color[2] = fraction * _colors[mod][2];
          color[3] = 1.0f;
          //glColor4f(color[0], color[1], color[2], color[3]);
          glColor4f(1,1,1,1);
          glVertex3f(firstVertex[0], firstVertex[1], firstVertex[2]);
          mod = y % 8;
          div = y / 8;
          fraction = 1.0f / (div + 1.0f);
          color[0] = fraction * _colors[mod][0];
          color[1] = fraction * _colors[mod][1];
          color[2] = fraction * _colors[mod][2];
          //glColor4f(color[0], color[1], color[2], color[3]);
          glColor4f(1,1,1,1);
          glVertex3f(secondVertex[0], secondVertex[1], secondVertex[2]);
        glEnd();
      }
    }
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the connections between the clones
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawCloneConnections()
{
#ifdef USING_GLVU  
  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
      for (unsigned int z = 0; z < _clonedVertices[x][y].size(); z++)
      {
        int first = _clonedVertices[x][y][z].first;
        int second = _clonedVertices[x][y][z].second;

        VEC3F firstVertex = *(_meshes[x]->vertices(first));
        VEC3F secondVertex = *(_meshes[y]->vertices(second));

        glBegin(GL_LINES);
          float color[4];
          int mod = _graphColors[x] % 8;
          int div = _graphColors[x] / 8;
          float fraction = 1.0f / (div + 1.0f);
          color[0] = fraction * _colors[mod][0];
          color[1] = fraction * _colors[mod][1];
          color[2] = fraction * _colors[mod][2];
          color[3] = 1.0f;
          glColor4f(color[0], color[1], color[2], color[3]);
          glVertex3f(firstVertex[0], firstVertex[1], firstVertex[2]);
          mod = y % 8;
          div = y / 8;
          fraction = 1.0f / (div + 1.0f);
          color[0] = fraction * _colors[mod][0];
          color[1] = fraction * _colors[mod][1];
          color[2] = fraction * _colors[mod][2];
          glColor4f(color[0], color[1], color[2], color[3]);
          glVertex3f(secondVertex[0], secondVertex[1], secondVertex[2]);
        glEnd();
      }
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw the blended surface mesh
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawBlendedMesh()
{
#ifdef USING_GLVU  
  // update the blended vertices
  for (unsigned int x = 0; x < _blendedVertices.size(); x++)
  {
    VEC3F* blendee = _blendedVertices[x];
    blendee->clear();

    vector<VEC3F*>& blenders = _verticesToBlend[x];
    for (unsigned int y = 0; y < blenders.size(); y++)
      *(blendee) += *(blenders[y]);

    *(blendee) *= 1.0 / blenders.size();
  }
  
  glColor4f(1.0f, 1.0f, 1.0f, 0.75f);
  for (unsigned int x = 0; x < _blendedSurfaceMesh.size(); x++)
    _blendedSurfaceMesh[x]->draw();

  /*
  // draw the blended vertices
  glColor4f(10.0f, 0.0f, 0.0f, 10.0f);
  glBegin(GL_POINTS);
  for (int x = 0; x < _blendedVertices.size(); x++)
  {
    VEC3F blendee = *(_blendedVertices[x]);
    glVertex3f(blendee[0], blendee[1], blendee[2]);
  }
  glEnd();
  */
#endif
}

//////////////////////////////////////////////////////////////////////
// compute the center of mass for all partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::updateCentersOfMass()
{
  for (int x = 0; x < _partitions; x++)
    _meshes[x]->computeCenterOfMass();
}

//////////////////////////////////////////////////////////////////////
// compute the center of mass over all the partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeCenterOfMass()
{
  Real totalMass = 0;
  _centerOfMass.clear();
  for (int x = 0; x < _partitions; x++)
  {
    _centerOfMass += _meshes[x]->totalMass() * _meshes[x]->centerOfMass();
    totalMass += _meshes[x]->totalMass();
  }

  _centerOfMass /= _partitions * totalMass;
}

//////////////////////////////////////////////////////////////////////
// Find the closest surface node
//////////////////////////////////////////////////////////////////////
int PARTITIONED_TET_MESH::closestPartition(VEC3F point)
{
  Real minDistance = 1e9;
  int minPartition = 0;
  
  for (int x = 0; x < _partitions; x++)
  {
    VEC3F* candidate = _meshes[x]->closestSurfaceNode(point);
    Real distance = norm2(point - *candidate);
    if (x == 0)
    {
      minDistance = distance;
      minPartition = 0;
    }

    if (distance < minDistance)
    {
      minDistance = distance;
      minPartition = x;
    }
  }

  return minPartition;
}

//////////////////////////////////////////////////////////////////////
// compute the lists of triangles along the interfaces
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeClonedTriangles()
{
  // allocate cloned triangle pointers
  if (_clonedTriangles != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedTriangles[x];
    delete[] _clonedTriangles;
  }
  _clonedTriangles = new vector<TRIANGLE*>*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _clonedTriangles[x] = new vector<TRIANGLE*>[_partitions];

  // allocate cloned triangle IDs 
  if (_clonedTriangleIDs != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedTriangleIDs[x];
    delete[] _clonedTriangleIDs;
  }
  _clonedTriangleIDs = new vector<int>*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _clonedTriangleIDs[x] = new vector<int>[_partitions];

  map<VEC3F*, int> triangleCounts;
  for (int x = 0; x < _partitions; x++)
  {
    TET_MESH* leftMesh = _meshes[x];
    vector<TRIANGLE*>& leftFaces = leftMesh->explicitSurfaceFaces();
    //vector<TRIANGLE*>& leftFaces = leftMesh->allFaces();

    for (int y = 0; y < _partitions; y++)
    {
      if (x == y) continue;
      TET_MESH* rightMesh = _meshes[y];
      vector<TRIANGLE*>& rightFaces = rightMesh->explicitSurfaceFaces();
      //vector<TRIANGLE*>& rightFaces = rightMesh->allFaces();
      map<VEC3F*, VEC3F*>& cloneMap = _clonedVertexMap[x][y];

      // build cloned versions of the right faces that use the 
      // indices from this partition
      vector<TRIANGLE> rightFacesThis;
      for (unsigned int z = 0; z < rightFaces.size(); z++)
      {
        TRIANGLE* rightFace = rightFaces[z];

        // if any of the clones don't exist the right face is no
        // longer a candidate
        if (cloneMap.find(rightFace->vertex(0)) == cloneMap.end() ||
            cloneMap.find(rightFace->vertex(1)) == cloneMap.end() ||
            cloneMap.find(rightFace->vertex(2)) == cloneMap.end())
          continue;

        VEC3F* v0 = cloneMap[rightFace->vertex(0)];
        VEC3F* v1 = cloneMap[rightFace->vertex(1)];
        VEC3F* v2 = cloneMap[rightFace->vertex(2)];

        rightFacesThis.push_back(TRIANGLE(v0, v1, v2));
      }

      for (unsigned int i = 0; i < leftFaces.size(); i++)
        for (unsigned int j = 0; j < rightFacesThis.size(); j++)
          if ((*leftFaces[i]) == rightFacesThis[j])
          {
            _clonedTriangles[x][y].push_back(leftFaces[i]);
            _clonedTriangleIDs[x][y].push_back(i);
            break;
          }
    }
  }

  /*
  for (int x = 0; x < _partitions; x++)
    for (int y = x; y < _partitions; y++)
    {
      if (x == y) continue;
      cout << " Cloned triangles between partitions " 
           << x << ", " << y << ": " << _clonedTriangles[x][y].size() << "\t" 
           << y << ", " << x << ": " << _clonedTriangles[y][x].size() << endl;
    }
    */

  /*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " DEACTIVATING SANITY CHECK FOR THE TREE -- SOME TRIANGLES ARE IN FACT CLONED." << endl;
  return;
  */

  // do a sanity check and make sure each triangles along the interface is unique
  for (int x = 0; x < _partitions; x++)
  {
    map<TRIANGLE*, bool> match;
    for (int y = 0; y < _partitions; y++)
      for (unsigned int z = 0; z < _clonedTriangles[x][y].size(); z++)
      {
        TRIANGLE* triangle = _clonedTriangles[x][y][z];
        if (match.find(triangle) != match.end())
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << "Sanity check failed!" << endl;
          exit(1);
        }
        match[triangle] = true;
      }
  }

  // do a sanity check and make sure no triangles with the same positions occur
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
      for (unsigned int z = 0; z < _clonedTriangles[x][y].size(); z++)
      {
        TRIANGLE& leftTriangle = *_clonedTriangles[x][y][z];
        for (unsigned int a = z + 1; a < _clonedTriangles[x][y].size(); a++)
          if (leftTriangle.positionsEqual(*_clonedTriangles[x][y][a]))
          {
            cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
            cout  << " Sanity check failed! " << endl;
            cout << " Triangle: " << leftTriangle;
            cout << " Matches other triangle: " << *_clonedTriangles[x][y][a];
            exit(1);
          }
      }

  // do a sanity check to make that the number of cloned triangles found is
  // symmetric
  vector<TRIANGLE*> noMatch;
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      if (x == y) continue;

      if (_clonedTriangles[x][y].size() == _clonedTriangles[y][x].size()) 
        continue;

      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Sanity check failed!" << endl;
      cout << " x: " << x << endl;
      cout << " y: " << y << endl;
      cout << " x,y cloned triangles: " << _clonedTriangles[x][y].size() << endl;
      cout << " y,x cloned triangles: " << _clonedTriangles[y][x].size() << endl;
      cout << " x,y cloned vertices: " << _clonedVertices[x][y].size() << endl;
      cout << " y,x cloned vertices: " << _clonedVertices[y][x].size() << endl;
      exit(0);

      for (unsigned int i = 0; i < _clonedTriangles[x][y].size(); i++)
      {
        bool matchFound = false;
        for (unsigned int j = 0; j < _clonedTriangles[y][x].size(); j++)
        {
          TRIANGLE leftTriangle = *_clonedTriangles[x][y][i];
          TRIANGLE rightTriangle = *_clonedTriangles[y][x][j];

          if (_unconstrainedPartition[x])
          {
            *(leftTriangle.vertex(0)) += _meshes[x]->centerOfMass();
            *(leftTriangle.vertex(1)) += _meshes[x]->centerOfMass();
            *(leftTriangle.vertex(2)) += _meshes[x]->centerOfMass();
          }
          if (_unconstrainedPartition[y])
          {
            *(rightTriangle.vertex(0)) += _meshes[y]->centerOfMass();
            *(rightTriangle.vertex(1)) += _meshes[y]->centerOfMass();
            *(rightTriangle.vertex(2)) += _meshes[y]->centerOfMass();
          }

          //if (_clonedTriangles[x][y][i]->positionsEqual(*_clonedTriangles[y][x][j]))
          if (leftTriangle.positionsEqual(rightTriangle))
            matchFound = true;
        }
        if (!matchFound)
        {
          cout << " No match found for: " << *_clonedTriangles[x][y][i] << endl;
          noMatch.push_back(_clonedTriangles[x][y][i]);
        }
      }

      cout << " x,y clones and counts: " << endl;
      vector<pair<int,int> >& clones = _clonedVertices[x][y];
      TET_MESH* mesh = _meshes[x];
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        VEC3F* vertex = mesh->vertices(clones[z].first);
        int count = triangleCounts[vertex];
        cout << " Vertex: " << *vertex << " count: " << count << endl;
      }
      cout << endl;

      cout << " y,x clones and counts: " << endl;
      clones = _clonedVertices[y][x];
      mesh = _meshes[y];
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        VEC3F* vertex = mesh->vertices(clones[z].first);
        int count = triangleCounts[vertex];
        cout << " Vertex: " << *vertex << " count: " << count << endl;
      }
    }
}

//////////////////////////////////////////////////////////////////////
// Compute the area between the partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeInterfaceAreas()
{
  if (_interfaceArea != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _interfaceArea[x];
    delete[] _interfaceArea;
  }
  _interfaceArea = new Real*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _interfaceArea[x] = new Real[_partitions];

  _maxInterfaceArea = 0.0;
  // does not give symmetric results?
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      _interfaceArea[x][y] = 0.0;
      for (unsigned int z = 0; z < _clonedTriangles[x][y].size(); z++)
      {
        TRIANGLE* triangle = _clonedTriangles[x][y][z];
        _interfaceArea[x][y] += triangle->area();
      }
      if (_interfaceArea[x][y] > _maxInterfaceArea)
        _maxInterfaceArea = _interfaceArea[x][y];
    }

  /*
  if (_maxInterfaceArea <= 0.0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "Cloned triangle lists are corrupt! Max area found was " << _maxInterfaceArea << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    exit(0);
  }
  */

  // sanity check the symmetry of the interface areas 
  for (int x = 0; x < _partitions; x++)
    for (int y = x+1; y < _partitions; y++)
    {
      Real xyArea = _interfaceArea[x][y];
      Real yxArea = _interfaceArea[y][x];

      if (fabs(xyArea) + fabs(yxArea) < 1e-8) continue;

      assert((fabs(xyArea) - fabs(yxArea)) / fabs(xyArea) < 1e-8);
    }

  /*
  for (int x = 0; x < _partitions; x++)
    for (int y = x+1; y < _partitions; y++)
    {
      _interfaceArea[x][y] = 0.0;

      for (unsigned int z = 0; z < _clonedTriangles[x][y].size(); z++)
      {
        TRIANGLE* triangle = _clonedTriangles[x][y][z];
        _interfaceArea[x][y] += triangle->area();
      }
      _interfaceArea[y][x] = _interfaceArea[x][y];
      if (_interfaceArea[x][y] > _maxInterfaceArea)
        _maxInterfaceArea = _interfaceArea[x][y];
    }
    */
}

//////////////////////////////////////////////////////////////////////
// Compute which cloned vertices are on the surface of the model
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeClonedSurfaceVertices()
{
  if (_originalIDs == NULL)
  {
    cout << " Must read in the correspondences first! " << endl;
    exit(0);
  }

  // clean up if necessary
  if (_clonedSurfaceVertices != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedSurfaceVertices[x];
    delete[] _clonedSurfaceVertices;
  }
  
  // allocate the clone table
  _clonedSurfaceVertices = new vector<pair<int, int> >*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _clonedSurfaceVertices[x] = new vector<pair<int, int> >[_partitions];

  vector<VEC3F>& originalVertices = _originalMesh->vertices();

  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int, int> >& clones = _clonedVertices[x][y];
      if (clones.size() == 0) continue;

      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int vertexID = clones[z].first;
        int originalID = _originalIDs[x][vertexID];

        VEC3F* vertex = &originalVertices[originalID];

        // if it's not on the surface, move on
        if (!_originalMesh->isOnSurface(vertex)) continue;
        
        _clonedSurfaceVertices[x][y].push_back(_clonedVertices[x][y][z]);
      }
    }
  /*
  for (int x = 0; x < _partitions; x++)
  {
    // hash the cloned triangles 
    map<TRIANGLE*, bool> clonedTriangleHash;
    for (int y = 0; y < _partitions; y++)
      for (unsigned int z = 0; z < _clonedTriangles[x][y].size(); z++)
      {
        TRIANGLE* triangle = _clonedTriangles[x][y][z];
        clonedTriangleHash[triangle] = true;
      }

    // hash the vertices of the uncloned triangles
    //vector<TRIANGLE*>& surfaceFaces = _meshes[x]->surfaceFaces();
    vector<TRIANGLE*>& surfaceFaces = _meshes[x]->explicitSurfaceFaces();
    map<VEC3F*, bool> unclonedTriangleVertices;
    for (unsigned int y = 0; y < surfaceFaces.size(); y++)
    {
      // don't hash if it is a cloned triangle
      TRIANGLE* triangle = surfaceFaces[y];
      if (clonedTriangleHash.find(triangle) != clonedTriangleHash.end())
        continue;
      
      for (int z = 0; z < 3; z++)
        unclonedTriangleVertices[triangle->vertex(z)] = true;
    }

    // for every other partition
    for (int y = 0; y < _partitions; y++)
    {
      // scan through all the cloned vertices between
      // the two partitions
      vector<pair<int,int> >& clonedVertices = _clonedVertices[x][y];
      for (unsigned int z = 0; z < clonedVertices.size(); z++)
      {
        // get a pointer to the vertex
        int vertexIndex = clonedVertices[z].first;
        VEC3F* vertex = _meshes[x]->vertices(vertexIndex);
        
        // if the vertex is also a member of an uncloned triangle,
        // it is a surface vertex, so store it
        if (unclonedTriangleVertices.find(vertex) !=
            unclonedTriangleVertices.end())
        {
          pair<int,int> toPush = clonedVertices[z];
          _clonedSurfaceVertices[x][y].push_back(toPush);
        }
      }
    }
  }
  */
}

//////////////////////////////////////////////////////////////////////
// Compute the surface mesh with blended vertices along the partition
// interfaces
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeBlendedSurfaceMesh()
{
  // keep track of the new vertices created
  // if "vertex" is a vertex in a TET_MESH
  // blendedVersion[vertex] returns the new vertex that corresponds
  // to its blended version, if one exists
  map<VEC3F*, VEC3F*> blendedVersion;
 
  // fast lookup of cloned triangles
  map<TRIANGLE*, bool> clonedTriangleHash;
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
      for (unsigned int z = 0; z < _clonedTriangles[x][y].size(); z++)
        clonedTriangleHash[_clonedTriangles[x][y][z]] = true;
  
  // for each partition
  for (int x = 0; x < _partitions; x++)
  {
    // build a table of its cloned *surface* vertices
    map<VEC3F*, vector<VEC3F*> > surfaceClones;
    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int,int> >& clonedVertices = _clonedSurfaceVertices[x][y];
      for (unsigned int z = 0; z < clonedVertices.size(); z++)
      {
        VEC3F* vertex0 = _meshes[x]->vertices(clonedVertices[z].first);
        VEC3F* vertex1 = _meshes[y]->vertices(clonedVertices[z].second);
        surfaceClones[vertex0].push_back(vertex1);
      }
    }

    // build a table of its interior cloned vertices
    map<VEC3F*, bool> interiorClones;
    for (int y = 0; y < _partitions; y++)
    {
      vector<pair<int,int> >& clonedVertices = _clonedVertices[x][y];
      for (unsigned int z = 0; z < clonedVertices.size(); z++)
      {
        VEC3F* vertex0 = _meshes[x]->vertices(clonedVertices[z].first);
        //VEC3F* vertex1 = _meshes[y]->vertices(clonedVertices[z].second);

        if (surfaceClones.find(vertex0) == surfaceClones.end())
          interiorClones[vertex0] = true;
      }
    }
    
    // for each triangle in the surface mesh
    //vector<TRIANGLE*>& surfaceFaces = _meshes[x]->surfaceFaces();
    vector<TRIANGLE*>& surfaceFaces = _meshes[x]->explicitSurfaceFaces();
    for (unsigned int y = 0; y < surfaceFaces.size(); y++)
    {
      TRIANGLE* triangle = surfaceFaces[y];
      
      // see if it's all interior clones
      int totalInteriorClones = 0;
      for (int z = 0; z < 3; z++)
        if (interiorClones.find(triangle->vertex(z)) != interiorClones.end())
          totalInteriorClones++;

      /*
      // if it's all interior clones, move on
      if (totalInteriorClones == 3) continue;
      */
      // if it's all interior clones, see if we are in the special case where
      // all the vertices were cloned but it is still a surface triangle by checking
      // if the whole triangle is cloned
      if (totalInteriorClones == 3)
        if (clonedTriangleHash.find(triangle) != clonedTriangleHash.end())
          continue;

      // see if it has any surface clones
      int totalSurfaceClones = 0;
      for (int z = 0; z < 3; z++)
        if (surfaceClones.find(triangle->vertex(z)) != surfaceClones.end())
          totalSurfaceClones++;

      // if it's all interior plus surface clones, but *not*
      // all surface clones (which is still a surface triangle), move on
      if (totalInteriorClones > 0 &&
          totalSurfaceClones + totalInteriorClones == 3)
        continue;

      // if there are no surface clones, just add the triangle as is
      if (totalSurfaceClones == 0)
      {
        _blendedSurfaceMesh.push_back(triangle);
        continue;
      }

      // create the new triangle, by just doing a copy
      TRIANGLE* newTriangle = new TRIANGLE(triangle);
      _blendedTriangles.push_back(newTriangle);

      // store the new triangle in the blended mesh
      _blendedSurfaceMesh.push_back(newTriangle);

      // detect the clones
      for (int z = 0; z < 3; z++)
      {
        // move on if this one isn't cloned
        if (surfaceClones.find(triangle->vertex(z)) == surfaceClones.end())
          continue;

        // see if this vertex had a blended vertex made for it already
        if (blendedVersion.find(triangle->vertex(z)) != blendedVersion.end())
        {
          // if it did, just set the cloned vertex to that one and
          // move on
          VEC3F* newVertex = blendedVersion[triangle->vertex(z)];
          newTriangle->setVertex(z, newVertex);
          continue;
        }

        // otherwise, this is the first time this cloned vertex has been
        // encountered, so create a new one and store its correspondance
        // with all the other cloned vertices
        VEC3F* oldVertex = triangle->vertex(z);
        VEC3F* newVertex = new VEC3F(*oldVertex);
        _blendedVertices.push_back(newVertex);

        // point the new triangle to this new vertex
        newTriangle->setVertex(z, newVertex);
        
        // hash this new vertex so all the other clones can find it as well
        blendedVersion[oldVertex] = newVertex;
        vector<VEC3F*>& otherClones = surfaceClones[oldVertex];
        for (unsigned int a = 0; a < otherClones.size(); a++)
          blendedVersion[otherClones[a]] = newVertex;

        // create a new vector (STL doesn't seem to like putting it
        // directly in "push_back"
        vector<VEC3F*> dummy;
        _verticesToBlend.push_back(dummy);
        
        // store all the vertices that will blend to form this one
        vector<VEC3F*>& toBlend = _verticesToBlend.back();
        toBlend.push_back(oldVertex);
        for (unsigned int a = 0; a < otherClones.size(); a++)
          toBlend.push_back(otherClones[a]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// write the partition graph out for coloring
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::writeDIMACS(const char* filename)
{
  // first allocate the graph
  bool** edges = new bool*[_partitions];
  for (int x = 0; x < _partitions; x++)
    edges[x] = new bool[_partitions];

  // set everybody initially to false
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
      edges[x][y] = false;

  // find the edges
  int totalEdges = 0;
  for (int x = 0; x < _partitions; x++)
    for (int y = x+1; y < _partitions; y++)
      if (_clonedVertices[x][y].size() > 0)
      {
        edges[x][y] = true;
        totalEdges++;
      }

  // write the graph out
  FILE* file = fopen(filename, "w");

  // write out a comment that identifies the mesh
  string originalFile = _originalMesh->filename();
  fprintf(file, "c DIMACs file for partitioned mesh: %s\n", originalFile.c_str());

  // header for total nodes and edges
  fprintf(file, "p edge %i %i\n", _partitions, totalEdges);

  // find the edges
  for (int x = 0; x < _partitions; x++)
    for (int y = x+1; y < _partitions; y++)
      if (edges[x][y])
        fprintf(file, "e %i %i\n", x+1, y+1);

  // all done, close up file
  fclose(file);
  
  // clean up
  for (int x = 0; x < _partitions; x++)
    delete[] edges[x];
  delete[] edges;
}

//////////////////////////////////////////////////////////////////////
// read in the graph coloring returned by Smallk
//////////////////////////////////////////////////////////////////////
bool PARTITIONED_TET_MESH::readDIMACS(const char* filename)
{
  // read the coloring in
  FILE* file = fopen(filename, "r");
  if (file == NULL) return false;

  float cpu;
  int pid;
  fscanf(file, "CLRS %i FROM DSATUR cpu =  %f pid = %i\n", &_totalGraphColors, &cpu, &pid);
  flush(cout);

  for (int x = 0; x < _partitions; x++)
  {
    int color;
    fscanf(file, "  %i", &color);
    _graphColors[x] = color-1;
  }

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// compute a graph coloring for the partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeGraphColoring(string prefix)
{
  //string tetMeshName = _originalMesh->filename();
  string tetMeshName = prefix;
  if (tetMeshName.length() == 0)
    tetMeshName = _originalMesh->filename();
  
  // output graph for coloring
  string dimacsName(tetMeshName);
  dimacsName += string(".dimacs");
  writeDIMACS(dimacsName.c_str());
  
  // stomp any previous coloring
  string remove("rm ");
  remove += dimacsName;
  remove += string(".res");
  system(remove.c_str());
  
  /*
  int colors = 2;
  int smallKFailed = 1;
  int seed = 12345;
  char buffer[256];

  // starting with 2 color, try to find a coloring up to 9
  // (smallK is only designed for a small number of colors)
  while (colors < 9 && smallKFailed)
  {
    cout << "==================================================" << endl;
    cout << " Trying " << colors << "-coloring" << endl;
    cout << "==================================================" << endl;

    // construct the system call string
    string smallk = string("./smallk/smallk ");
    smallk += string(" ");
    smallk += tetMeshName.c_str();
    smallk += string(".dimacs");
    sprintf(buffer, "%i", seed);
    smallk += string(" ");
    smallk += string(buffer);
    sprintf(buffer, "%i", colors);
    smallk += string(" ");
    smallk += string(buffer);

    // call SmallK
    smallKFailed = system(smallk.c_str());

    // bump up the number of colors to try
    colors++;
  }
  // bad news if we still didn't find a coloring -- need to look into
  // a bigger scale coloring package
  if (smallKFailed)
  {
    cout << "==================================================" << endl;
    cout << " Smallk failed to find any coloring!!!!" << endl;
    cout << "==================================================" << endl;
    exit(1);
  }
  cout << "==================================================" << endl;
  cout << " Smallk found a " << colors - 1 << "-coloring" << endl;
  cout << "==================================================" << endl;
  */
  int seed = 12345;
  char buffer[256];

  cout << "==================================================" << endl;
  cout << " Coloring using DSATUR" << endl;
  cout << "==================================================" << endl;

  // construct the system call string
  string dsatur = string("./bin/dsatur");
  dsatur += string(" ");
  dsatur += tetMeshName.c_str();
  dsatur += string(".dimacs");
  sprintf(buffer, "%i", seed);
  dsatur += string(" ");
  dsatur += string(buffer);

  // call DSATUR
  system(dsatur.c_str());

  // read in the coloring
  string dimacsResult(dimacsName);
  dimacsResult += string(".res");
  readDIMACS(dimacsResult.c_str());
}

//////////////////////////////////////////////////////////////////////
// Update the full mesh of all the partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::updateFullMeshes()
{
  for (int x = 0; x < _partitions; x++)
    _meshes[x]->updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
// normalize the embedding mesh to the size of the tet mesh
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::normalizeEmbedding(OBJ* obj)
{
  // get the bounding box of the embedding
  Real objBoundingBox[6];
  obj->BoundingBox(objBoundingBox);
  VEC3F minObj;
  minObj[0] = objBoundingBox[0];
  minObj[1] = objBoundingBox[2];
  minObj[2] = objBoundingBox[4];
  VEC3F maxObj;
  maxObj[0] = objBoundingBox[1];
  maxObj[1] = objBoundingBox[3];
  maxObj[2] = objBoundingBox[5];
  VEC3F diffObj = maxObj - minObj;
  Real objScale = diffObj.maxElement();

  // get the bounding box of the tet meshes
  VEC3F minTets;
  VEC3F maxTets;
  for (int x = 0; x < _partitions; x++)
  {
    Real tetsBoundingBox[6];
    _meshes[x]->boundingBox(tetsBoundingBox);

    if (x == 0)
    {
      minTets[0] = tetsBoundingBox[0];
      minTets[1] = tetsBoundingBox[2];
      minTets[2] = tetsBoundingBox[4];
      maxTets[0] = tetsBoundingBox[1];
      maxTets[1] = tetsBoundingBox[3];
      maxTets[2] = tetsBoundingBox[5];
    }
    else
    {
      minTets[0] = tetsBoundingBox[0] < minTets[0] ? tetsBoundingBox[0] : minTets[0];
      minTets[1] = tetsBoundingBox[2] < minTets[1] ? tetsBoundingBox[2] : minTets[1];
      minTets[2] = tetsBoundingBox[4] < minTets[2] ? tetsBoundingBox[4] : minTets[2];

      maxTets[0] = tetsBoundingBox[1] > maxTets[0] ? tetsBoundingBox[1] : maxTets[0];
      maxTets[1] = tetsBoundingBox[3] > maxTets[1] ? tetsBoundingBox[3] : maxTets[1];
      maxTets[2] = tetsBoundingBox[5] > maxTets[2] ? tetsBoundingBox[5] : maxTets[2];
    }
  }

  VEC3F diffTets = maxTets - minTets;
  Real tetsScale = diffTets.maxElement();

  // translate the embedding to the origin
  vector<VEC3>& vertices = obj->vertices;
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] -= minObj;

  // scale to unity, then to tet mesh size
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] *= tetsScale / objScale;

  // translate embedding to tet mesh min
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] += minTets;
}

//////////////////////////////////////////////////////////////////////
// compute the embedding of this surface into the tet meshes
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeEmbedding(OBJ* obj)
{
  vector<VEC3F>& vertices = obj->vertices;

  // compute a list of all the surface tetrahedra
  // the first index is a pair of <partition #, tet #>
  map<pair<int, int>, bool> surfaceTetIDsHash;
  for (int i = 0; i < _partitions; i++)
  {
    vector<VEC3F*>& surfaceVertices = _meshes[i]->surfaceVertices();
    cout << " Surface vertices for partition " << i << ": " << surfaceVertices.size() << endl;
    for (unsigned int x = 0; x < surfaceVertices.size(); x++)
    {
      vector<int>& tetMembership = _meshes[i]->tetMembership(surfaceVertices[x]);
      //cout << " Tet membership size: " << tetMembership.size() << endl;
      for (unsigned int y = 0; y < tetMembership.size(); y++)
      {
        pair<int, int> index(i, tetMembership[y]);
        surfaceTetIDsHash[index] = true;
      }
    }
    cout << " Surface hash new size: " << surfaceTetIDsHash.size() << endl;
  }

  vector<pair<int,int> > surfaceTetIDs;
  map<pair<int, int>,bool>::iterator iter;
  for (iter = surfaceTetIDsHash.begin(); iter != surfaceTetIDsHash.end(); iter++)
    surfaceTetIDs.push_back(iter->first);

  cout << " Total surface tets: " << surfaceTetIDs.size() << endl;

  // compute the center of mass of all the tetrahedra
  vector<VEC3F> centers;
  for (unsigned int x = 0; x < surfaceTetIDs.size(); x++)
  {
    int partition = surfaceTetIDs[x].first;
    int tetID = surfaceTetIDs[x].second;
    vector<TET>& tets = _meshes[partition]->tets();
    centers.push_back(tets[tetID].center());
  }
  cout << " Total tet centers: " << surfaceTetIDs.size() << endl;

  cout << " Computing nearest tets ... "; flush(cout);
  // for each point, compute the nearest tet
  _tetEmbeddings.clear();
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    Real minDistance = 0.0;
    pair<int, int> nearestTet;
    for (unsigned int y = 0; y < centers.size(); y++)
    {
      VEC3F diff = vertices[x] - centers[y];
      Real distance = norm(diff);

      if (y == 0 || distance < minDistance)
      {
        minDistance = distance;
        nearestTet = surfaceTetIDs[y];
      }
    }
    _tetEmbeddings.push_back(nearestTet);

    if ((int)(x % (int)(vertices.size() / 10)) == 0)
    {
      cout << 100 * ((Real)x / vertices.size()) << "% ";
      flush(cout);
    }
  }
  cout << endl;

  // for each nearest tet, compute the barycentric coords
  cout << " Computing barycentric coords ... "; flush(cout);
  _barycentricEmbeddings.clear();
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    // get the nearest tet
    int partition = _tetEmbeddings[x].first;
    vector<TET>& tets = _meshes[partition]->tets();
    TET& nearestTet = tets[_tetEmbeddings[x].second];
    VEC3F& v0 = *(nearestTet.vertices[0]);
    VEC3F& v1 = *(nearestTet.vertices[1]);
    VEC3F& v2 = *(nearestTet.vertices[2]);
    VEC3F& v3 = *(nearestTet.vertices[3]);

    // compute the barycenter matrix
    MATRIX3 baryMatrix;
    baryMatrix(0,0) = v0[0] - v3[0];
    baryMatrix(0,1) = v1[0] - v3[0];
    baryMatrix(0,2) = v2[0] - v3[0];
    baryMatrix(1,0) = v0[1] - v3[1];
    baryMatrix(1,1) = v1[1] - v3[1];
    baryMatrix(1,2) = v2[1] - v3[1];
    baryMatrix(2,0) = v0[2] - v3[2];
    baryMatrix(2,1) = v1[2] - v3[2];
    baryMatrix(2,2) = v2[2] - v3[2];

    // invert it
    MATRIX3 inverse = baryMatrix.inverse();

    // solve for the coordinates
    VEC3F diff = vertices[x] - v3;
    VEC3F final = inverse * diff;
    _barycentricEmbeddings.push_back(final);

    if (x % (int)(vertices.size() / 10) == 0)
    {
      cout << 100 * ((Real)x / vertices.size()) << "% ";
      flush(cout);
    }
  }
  cout << endl;

  // save out the embedding
  writeEmbedding();
}

//////////////////////////////////////////////////////////////////////
// write the partitioned embedding to a file
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::writeEmbedding(string filename)
{
  if (filename.length() == 0)
  {
    filename = _filename;
    filename += string(".partitioned.embedding");
  }
  FILE* file = fopen(filename.c_str(), "wb");

  int totalVertices = _tetEmbeddings.size();
  fwrite((void*)&totalVertices, sizeof(int), 1, file);

  // write out the tets
  for (int x = 0; x < totalVertices; x++)
  {
    int partition = _tetEmbeddings[x].first;
    int tetID = _tetEmbeddings[x].second;
    fwrite((void*)&partition, sizeof(int), 1, file);
    fwrite((void*)&tetID, sizeof(int), 1, file);
  }

  // write out the coordinates
  for (int x = 0; x < totalVertices; x++)
  {
    double coords[3];
    coords[0] = _barycentricEmbeddings[x][0];
    coords[1] = _barycentricEmbeddings[x][1];
    coords[2] = _barycentricEmbeddings[x][2];
    fwrite((void*)coords, sizeof(double), 3, file);
  }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// embed this surface into the mesh
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::embedMesh(OBJ* obj)
{
  // store the embedded mesh
  _embeddedMesh = obj;

  // see if a precomputed one exists
  if (!readEmbedding())
  {
    cout << __FILE__ << " " << __LINE__ << " : " << endl;
    cout << "No embedding has been computed for this mesh!!!!" << endl;
    exit(0); 
  }
}

//////////////////////////////////////////////////////////////////////
// try to read an embedding file
//////////////////////////////////////////////////////////////////////
bool PARTITIONED_TET_MESH::readEmbedding(string filename)
{
  if (filename.length() == 0)
  {
    filename = _filename;
    filename += string(".partitioned.embedding");
  }
  FILE* file = fopen(filename.c_str(), "rb");

  if (file == NULL)
  {
    cout << " Tried to open embedding file " << filename.c_str() << endl;
    return false;
  }

  cout << " Found an embedding file: " << filename.c_str() << endl;

  int totalVertices;
  fread((void*)&totalVertices, sizeof(int), 1, file);

  // read in the tets
  _tetEmbeddings.clear();
  for (int x = 0; x < totalVertices; x++)
  {
    int partition;
    int tetID;
    fread((void*)&partition, sizeof(int), 1, file);
    fread((void*)&tetID, sizeof(int), 1, file);
    pair<int, int> newPair(partition, tetID);
    _tetEmbeddings.push_back(newPair);
  }

  // read in the coordinates
  _barycentricEmbeddings.clear();
  for (int x = 0; x < totalVertices; x++)
  {
    double coords[3];
    fread((void*)coords, sizeof(double), 3, file);
    VEC3F baryCoords;
    baryCoords[0] = coords[0];
    baryCoords[1] = coords[1];
    baryCoords[2] = coords[2];
    _barycentricEmbeddings.push_back(baryCoords);
  }

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// update the position of the embedding
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::updateEmbedding()
{
  if (_embeddedMesh == NULL)
  {
    cout << " Embedded mesh has not been set!" << endl;
    return;
  }

  // if _unconstrainedPartition is not populated
  if (_unconstrainedPartition.size() == 0)
    for (int x = 0; x < _partitions; x++)
      if (_meshes[x]->unconstrained())
        _unconstrainedPartition.push_back(true);
      else
        _unconstrainedPartition.push_back(false);

  // recompute positions based on barycentric coords of embedding
  vector<VEC3>& vertices = _embeddedMesh->vertices;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    int partition = _tetEmbeddings[x].first;
    vector<TET>& tets = _meshes[partition]->tets();
    TET& tet = tets[_tetEmbeddings[x].second];

    VEC3F& v0 = *(tet.vertices[0]);
    VEC3F& v1 = *(tet.vertices[1]);
    VEC3F& v2 = *(tet.vertices[2]);
    VEC3F& v3 = *(tet.vertices[3]);

    VEC3F coords = _barycentricEmbeddings[x];
    Real subtract = 1.0 - coords[0] - coords[1] - coords[2];
    VEC3F update = subtract * v3 + v0 * coords[0] + v1 * coords[1] + v2 * coords[2];
    vertices[x] = update;
  }
}

//////////////////////////////////////////////////////////////////////
// Output surface triangles to RenderMan
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawSurfaceToRenderMan()
{
  for (int x = 0; x < _partitions; x++)
  {
#if USING_RENDERMAN
    RiColor(_colors[x % 8]);
#endif
    _meshes[x]->drawSurfaceToRenderMan();
  }
}

//////////////////////////////////////////////////////////////////////
// get a vector of <partition, nodeID in partition> pairs that 
// correspond to the vertexID in the original mesh
//////////////////////////////////////////////////////////////////////
vector<pair<int, int> > PARTITIONED_TET_MESH::partitionIDs(int vertexID)
{
  vector<pair<int, int> > toReturn;

  // argh, so hideous, I remember why I don't like multimaps now
  pair<multimap<int, pair<int, int> >::iterator,
       multimap<int, pair<int, int> >::iterator> hideous;

  hideous = _partitionIDs.equal_range(vertexID);

  multimap<int, pair<int, int> >::iterator iter;
  for (iter = hideous.first; iter != hideous.second; iter++)
  {
    toReturn.push_back(iter->second);
  }
  return toReturn;
}

//////////////////////////////////////////////////////////////////////
// Output embedded triangles to OGL
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawFilteredEmbeddingSubset()
{
  updateEmbedding();
  glColor4f(1,1,1,1);
  //_embeddedMesh->ComputeVertexNormals();
  _embeddedMesh->drawFilteredSubset();
}

//////////////////////////////////////////////////////////////////////
// Output embedded triangles to OGL
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawFilteredEmbedding()
{
  updateEmbedding();
  glColor4f(1,1,1,1);
  //_embeddedMesh->ComputeVertexNormals();
  _embeddedMesh->drawFiltered();
}

//////////////////////////////////////////////////////////////////////
// Output embedded triangles to OGL
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawEmbedding()
{
  updateEmbedding();
  //glColor4f(1,1,1,1);
  _embeddedMesh->ComputeVertexNormals();
  _embeddedMesh->draw();
}

//////////////////////////////////////////////////////////////////////
// Output embedded triangles to RenderMan
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawEmbeddingToRenderMan()
{
#ifdef USING_RENDERMAN  
  _embeddedMesh->ComputeVertexNormals();
  for (int x = 0; x < 10; x++)
    _embeddedMesh->SmoothVertexNormals();

  vector<VEC3F>& vertices = _embeddedMesh->vertices;
  vector<VEC3F>& normals = _embeddedMesh->normals;
  vector<OBJ::Face>& faces = _embeddedMesh->faces;

  int faceStart = 0;
  int vertexStart = 0;
  int totalFaces = faces.size();
  int totalVertices = vertices.size();

  // init point locations
  RtPoint* P = new RtPoint[totalVertices];
  for (int y = 0; y < totalVertices; y++)
  {
    P[y][0] = vertices[vertexStart + y][0];
    P[y][1] = vertices[vertexStart + y][1];
    P[y][2] = vertices[vertexStart + y][2];
  }

  // init normals
  RtPoint* N = new RtPoint[totalVertices];
  for (unsigned int y = 0; y < normals.size(); y++)
  {
    N[y][0] = normals[y][0];
    N[y][1] = normals[y][1];
    N[y][2] = normals[y][2];
  }

  // init all to triangles
  RtInt* nvertices = new RtInt[totalFaces];
  for (int y = 0; y < totalFaces; y++)
    nvertices[y] = 3;

  // init faces
  RtInt* faceIndices= new RtInt[3 * totalFaces];
  for (int y = 0; y < totalFaces; y++)
  {
    OBJ::Face face = faces[faceStart + y];
    faceIndices[y * 3]     = face.vertices[0] - vertexStart;
    faceIndices[y * 3 + 1] = face.vertices[1] - vertexStart;
    faceIndices[y * 3 + 2] = face.vertices[2] - vertexStart;
  }

  RiPointsPolygons(totalFaces, nvertices, faceIndices, 
                   RI_P, P, 
                   RI_N, N, 
                   RI_NULL);
  delete[] faceIndices;
  delete[] nvertices;
  delete[] P;
  delete[] N;
#endif
}

//////////////////////////////////////////////////////////////////////
// compute the full sparse interface stiffness matrix
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeFullInterfaceStiffness(BLOCK_SPARSE_MATRIX& springMatrix)
{
  // populate the matrices
  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
    {
      // get the clones between the two interfaces
      vector<pair<int, int> >& clones = _clonedVertices[x][y];
      if (clones.size() == 0) continue;

      TET_MESH& meshLeft = *_meshes[x];
      TET_MESH& meshRight = *_meshes[y];

      SPARSE_MATRIX diagonalLeft(meshLeft.TET_MESH::rank(), meshLeft.TET_MESH::rank());
      SPARSE_MATRIX diagonalRight(meshRight.TET_MESH::rank(), meshRight.TET_MESH::rank());
      SPARSE_MATRIX offDiagonal(meshLeft.TET_MESH::rank(), meshRight.TET_MESH::rank());

      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int leftIndex = 3 * clones[z].first;
        int rightIndex = 3 * clones[z].second;
        diagonalLeft(leftIndex, leftIndex) += 1.0;
        diagonalLeft(leftIndex + 1, leftIndex + 1) += 1.0;
        diagonalLeft(leftIndex + 2, leftIndex + 2) += 1.0;

        diagonalRight(rightIndex, rightIndex) += 1.0;
        diagonalRight(rightIndex + 1, rightIndex + 1) += 1.0;
        diagonalRight(rightIndex + 2, rightIndex + 2) += 1.0;

        offDiagonal(leftIndex, rightIndex) -= 1.0;
        offDiagonal(leftIndex + 1, rightIndex + 1) -= 1.0;
        offDiagonal(leftIndex + 2, rightIndex + 2) -= 1.0;
      }
    
      diagonalLeft *= _interfaceArea[x][y] / _maxInterfaceArea;
      diagonalRight *= _interfaceArea[y][x] / _maxInterfaceArea;
      offDiagonal *= _interfaceArea[y][x] / _maxInterfaceArea;

      SPARSE_MATRIX offDiagonalTranspose = offDiagonal.transpose();
      springMatrix.add(diagonalLeft, x,x);
      springMatrix.add(diagonalRight, y,y);
      springMatrix.add(offDiagonal, x,y);
      springMatrix.add(offDiagonalTranspose, y,x);
    }
}

//////////////////////////////////////////////////////////////////////
// compute the full diagonal sparse interface stiffness matrix
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeDiagonalInterfaceStiffness(BLOCK_SPARSE_MATRIX& diagonal)
{
  // populate the matrices
  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
    {
      // get the clones between the two interfaces
      vector<pair<int, int> >& clones = _clonedVertices[x][y];
      if (clones.size() == 0) continue;

      TET_MESH& meshLeft = *_meshes[x];

      SPARSE_MATRIX diagonalBlock(meshLeft.TET_MESH::rank(), meshLeft.TET_MESH::rank());

      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int leftIndex = clones[z].first;
        diagonalBlock(leftIndex, leftIndex) += 1.0;
      }
    
      diagonalBlock *= _interfaceArea[x][y] / _maxInterfaceArea;
      diagonal.add(diagonalBlock, x,x);
    }
}

//////////////////////////////////////////////////////////////////////
// compute the full off diagonal sparse interface stiffness matrix
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeOffDiagonalInterfaceStiffness(BLOCK_SPARSE_MATRIX& offDiagonal)
{
  // populate the matrices
  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
    {
      // get the clones between the two interfaces
      vector<pair<int, int> >& clones = _clonedVertices[x][y];
      if (clones.size() == 0) continue;

      TET_MESH& meshLeft = *_meshes[x];
      TET_MESH& meshRight = *_meshes[y];

      SPARSE_MATRIX offDiagonalBlock(meshLeft.TET_MESH::rank(), meshRight.TET_MESH::rank());

      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int leftIndex = clones[z].first;
        int rightIndex = clones[z].second;
        offDiagonalBlock(leftIndex, rightIndex) -= 1.0;
      }
    
      offDiagonalBlock *= _interfaceArea[y][x] / _maxInterfaceArea;

      SPARSE_MATRIX offDiagonalTranspose = offDiagonalBlock.transpose();
      offDiagonal.add(offDiagonalBlock, x,y);
      offDiagonal.add(offDiagonalTranspose, y,x);
    }
}

//////////////////////////////////////////////////////////////////////////////
// extract the portion of this vector that is in the submesh
//////////////////////////////////////////////////////////////////////////////
VECTOR PARTITIONED_TET_MESH::getSubvector(const int partition, const VECTOR& fullVector)
{
  TET_MESH* subMesh= _meshes[partition];
  //TET_MESH* fullMesh = _originalMesh;
  
  // set a matrix to the full vector so we can use SUBMATRIX
  //MATRIX snapshot(fullMesh->x().size(), 1);
  MATRIX snapshot(fullVector.size(), 1);
  snapshot.setColumn(fullVector, 0);

  int partitionNodes = subMesh->unconstrainedNodes();
  MATRIX subSnapshotMatrix(partitionNodes * 3, 1);

  // for each vertex in the partition
  for (int y = 0; y < partitionNodes; y++)
  {
    // get the vertex index in the original mesh
    int ID = originalID(partition, y);

    // get the submatrix corresponding to that vertex
    SUBMATRIX vertexBasis(snapshot, ID * 3, 3);

    // copy the vertex basis into the subbasis
    vertexBasis.copiesInto(subSnapshotMatrix, 3 * y);
  }

  return subSnapshotMatrix.getColumn(0);
}

//////////////////////////////////////////////////////////////////////////////
// compute the area of the interface represented by each cloned vertex
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeClonedVertexAreas()
{
  // stomp any old results
  _clonedVertexAreas.clear();

  // distribute the areas to the vertices
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      vector<TRIANGLE*>& triangles = _clonedTriangles[x][y];

      if (_clonedTriangles[x][y].size() != _clonedTriangles[y][x].size())
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " left: " << x << endl;
        cout << " right: " << y << endl;
        cout << " left  clones: " << _clonedTriangles[x][y].size() << endl;
        cout << " right clones: " << _clonedTriangles[y][x].size() << endl;
      }
      assert(_clonedTriangles[x][y].size() == _clonedTriangles[y][x].size());

      for (unsigned int z = 0; z < triangles.size(); z++)
      {
        Real areaThird = triangles[z]->area() / 3.0;
        _clonedVertexAreas[triangles[z]->vertex(0)] += areaThird;
        _clonedVertexAreas[triangles[z]->vertex(1)] += areaThird;
        _clonedVertexAreas[triangles[z]->vertex(2)] += areaThird;
      }
    }

  // find the max area
  _maxClonedVertexArea = 0.0;
  map<VEC3F*, Real>::iterator iter;
  for (iter = _clonedVertexAreas.begin(); iter != _clonedVertexAreas.end(); iter++)
  {
    if (iter->second > _maxClonedVertexArea)
      _maxClonedVertexArea = iter->second;
  }
}

//////////////////////////////////////////////////////////////////////////////
// return the total number of unconstrianed domains
//////////////////////////////////////////////////////////////////////////////
int PARTITIONED_TET_MESH::totalUnconstrained()
{
  int total = 0;
  for (unsigned int x = 0; x < _unconstrainedPartition.size(); x++)
    if (_unconstrainedPartition[x])
      total++;
  return total;
}

//////////////////////////////////////////////////////////////////////////////
// get the IDs of all the unconstrained partitions
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::getUnconstrainedIDs(vector<int>& unconstrained)
{
  for (unsigned int x = 0; x < _unconstrainedPartition.size(); x++)
    if (_unconstrainedPartition[x])
      unconstrained.push_back(x);
}

//////////////////////////////////////////////////////////////////////
// Return the rigid rotation for a given partition
//////////////////////////////////////////////////////////////////////
QUATERNION PARTITIONED_TET_MESH::quaternionRotation(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_TET_MESH*)_meshes[partition])->rotationQuaternion();

  return QUATERNION(MATRIX3::I());
}

//////////////////////////////////////////////////////////////////////
// Return the rigid rotation for a given partition from the previous
// timestep
//////////////////////////////////////////////////////////////////////
QUATERNION PARTITIONED_TET_MESH::quaternionRotationOld(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_TET_MESH*)_meshes[partition])->rotationQuaternionOld();

  return QUATERNION(MATRIX3::I());
}

//////////////////////////////////////////////////////////////////////////////
// populate _clonedVertexMap
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::computeClonedVertexMap()
{
  // allocate the clone map
  _clonedVertexMap = new map<VEC3F*, VEC3F*>*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _clonedVertexMap[x] = new map<VEC3F*, VEC3F*>[_partitions];
  
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      // check if there are any clones
      vector<pair<int, int> >& clones = _clonedVertices[x][y];
      if (clones.size() == 0) continue;

      // set up for the hash
      map<VEC3F*, VEC3F*>& cloneMap = _clonedVertexMap[x][y];
      TET_MESH* xMesh = _meshes[x];
      TET_MESH* yMesh = _meshes[y];

      // retrieve each vertex and do the hashes
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int xIndex = clones[z].first;
        int yIndex = clones[z].second;

        VEC3F* xVertex = xMesh->vertices(xIndex);
        VEC3F* yVertex = yMesh->vertices(yIndex);

        cloneMap[yVertex] = xVertex;
      }
    }

  // sanity check
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
      assert(_clonedVertexMap[x][y].size() == _clonedVertexMap[y][x].size());
}

//////////////////////////////////////////////////////////////////////////////
// how many unique clone pairs are there?
//////////////////////////////////////////////////////////////////////////////
int PARTITIONED_TET_MESH::totalClones()
{
  int total = 0;
  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
        total += _clonedVertices[x][y].size();

  return total;
}

//////////////////////////////////////////////////////////////////////
// Scatter new q into all partition qs
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::scatterSuperX(BLOCK_VECTOR& superX)
{ 
  for (int x = 0; x < _partitions; x++)
  {
    VECTOR& q = (_meshes[x])->x();
    q.equals(*superX(x));
  }
}

//////////////////////////////////////////////////////////////////////
// Gather all the partitions qs into _superQ
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::gatherSuperX(BLOCK_VECTOR& superX)
{
  for (int x = 0; x < _partitions; x++)
  {
    VECTOR& q = _meshes[x]->x();
    (*superX(x)) = q;
  }
}

//////////////////////////////////////////////////////////////////////
// Compute the stiffness matrix of the interface in block form
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::updateBlockInterfaceStiffness()
{
  _diagonalInterfaceStiffness.resizeAndWipe(_partitions, _partitions);
  _offDiagonalInterfaceStiffness.resizeAndWipe(_partitions, _partitions);

  // do each combination -- full vector version
  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
    {
      vector<pair<int, int> >& clones = clonedVertices(x,y);

      // see if there is even an interface here
      if (clones.size() == 0) continue;

      int leftDofs = _meshes[x]->dofs(); 
      SPARSE_MATRIX springConstLeft(leftDofs, leftDofs);
      int rightDofs = _meshes[y]->dofs(); 
      SPARSE_MATRIX springConstRight(rightDofs, rightDofs);

      Real k = springConst(x,y);
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int leftIndex = clones[z].first;
        int rightIndex = clones[z].second;

        if (_meshes[x]->isConstrained(leftIndex)) continue;

        springConstLeft(3 * leftIndex, 3 * leftIndex) = k;
        springConstLeft(3 * leftIndex + 1, 3 * leftIndex + 1) = k;
        springConstLeft(3 * leftIndex + 2, 3 * leftIndex + 2) = k;
        
        springConstRight(3 * rightIndex, 3 * rightIndex) = k;
        springConstRight(3 * rightIndex + 1, 3 * rightIndex + 1) = k;
        springConstRight(3 * rightIndex + 2, 3 * rightIndex + 2) = k;
      }

      // store self-products
      _diagonalInterfaceStiffness.add(springConstLeft, x, x);
      _diagonalInterfaceStiffness.add(springConstRight, y, y);

      // generate off-diagonal matrix
      SPARSE_MATRIX offDiagonalFull(leftDofs, rightDofs);
      MATRIX3 Ri = quaternionRotation(x).toExplicitMatrix3x3();
      MATRIX3 Rj = quaternionRotation(y).toExplicitMatrix3x3();
      MATRIX3 RiTRj = Ri.transpose() * Rj;
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int leftIndex = clones[z].first;
        int rightIndex = clones[z].second;
        MATRIX3 rotated = -k * RiTRj;
        if (_meshes[x]->isConstrained(leftIndex)) continue;

        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            offDiagonalFull(3 * leftIndex + i, 3 * rightIndex + j) = rotated(i,j);
      }

      // store off-diagonal products
      _offDiagonalInterfaceStiffness.equals(offDiagonalFull, x, y);
      SPARSE_MATRIX transpose = offDiagonalFull.transpose();
      _offDiagonalInterfaceStiffness.equals(transpose, y, x);
    }
}

//////////////////////////////////////////////////////////////////////
// Return the lagrange multiplier from the previous solve
//////////////////////////////////////////////////////////////////////
Real PARTITIONED_TET_MESH::rotationLambda(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_TET_MESH*)_meshes[partition])->rotationLambda();

  return 0.0;
}

//////////////////////////////////////////////////////////////////////
// Return the rigid rotation for a given partition
//////////////////////////////////////////////////////////////////////
VEC3F PARTITIONED_TET_MESH::rigidTranslation(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_TET_MESH*)_meshes[partition])->rigidTranslation();

  VEC3F zero;
  return zero;
}

//////////////////////////////////////////////////////////////////////
// Draw rotated axes for each rigid frame
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawRigidFrames()
{
#if USING_GLVU  
  for (int x = 0; x < _partitions; x++)
    if (_unconstrainedPartition[x])
      ((UNCONSTRAINED_TET_MESH*)_meshes[x])->drawRigidFrame();
#endif
}

//////////////////////////////////////////////////////////////////////
// Return the rigid rotation for a given partition
//////////////////////////////////////////////////////////////////////
MATRIX3 PARTITIONED_TET_MESH::rigidRotation(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_TET_MESH*)_meshes[partition])->rotationQuaternion().toExplicitMatrix3x3();

  return MATRIX3::I();
}

//////////////////////////////////////////////////////////////////////
// Draw constrained nodes on all partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::drawConstrainedNodes()
{
  for (int x = 0; x < _partitions; x++)
    _meshes[x]->drawConstrainedNodes();
}

//////////////////////////////////////////////////////////////////////
// assemble a vector of all the cloned rest vertices along the 
// (x,y) partition boundary
//////////////////////////////////////////////////////////////////////
VECTOR PARTITIONED_TET_MESH::interfaceRestVertices(int x, int y)
{
  vector<pair<int,int> > clones = clonedVertices(x,y);
  VECTOR final(clones.size() * 3);

  vector<VEC3F>& restVertices = _meshes[x]->restPose();

  // if there are no clones, this probably shouldn't be getting called
  assert(clones.size() > 0);

  for (unsigned int i = 0; i < clones.size(); i++)
  {
    int index = clones[i].first;
    VEC3F& restVertex = restVertices[index];
    final[3 * i] = restVertex[0];
    final[3 * i + 1] = restVertex[1];
    final[3 * i + 2] = restVertex[2];
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// are partitions x and y neighbors?
//////////////////////////////////////////////////////////////////////
bool PARTITIONED_TET_MESH::neighbors(int x, int y) { 
  //return _interfaceArea[x][y] > 0.0; 
  return _clonedVertices[x][y].size() != 0;
  //return _clonedVertices[x][y].size() > 3;
}

//////////////////////////////////////////////////////////////////////
// Compute the mass matrix along an interface
//////////////////////////////////////////////////////////////////////
Real PARTITIONED_TET_MESH::interfaceMass(int x, int y)
{
  Real totalMass = 0;
  vector<pair<int,int> > clones = clonedVertices(x,y);

  for (unsigned int i = 0; i < clones.size(); i++)
  {
    int leftIndex = clones[i].first;
    totalMass += _meshes[x]->mass(leftIndex);
  }

  return totalMass;
}

//////////////////////////////////////////////////////////////////////
// get the original tet index for tet 'tetID' in partition 'partition'
//////////////////////////////////////////////////////////////////////
int PARTITIONED_TET_MESH::originalTetID(int partition, int tetID)
{
  vector<int> vertexIDs;
  TET& submeshTet = _meshes[partition]->tets()[tetID];

  // get all the original IDs
  for (int x = 0; x < 4; x++)
  {
    int submeshID = _meshes[partition]->vertexID(submeshTet.vertices[x]);
    vertexIDs.push_back(originalID(partition, submeshID));
  }

  // find the tet in the full mesh
  vector<TET>& tets = _originalMesh->tets();
  int found = -1;
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    int matches = 0;
    for (int y = 0; y < 4; y++)
      for (int z = 0; z < 4; z++)
      {
        int currentID = _originalMesh->vertexID(tets[x].vertices[z]);
        if (currentID == vertexIDs[y])
          matches++;
      }
    if (matches == 4)
      found = x;
  }

  assert(found != -1);
  return found;
}

//////////////////////////////////////////////////////////////////////
// read/write cloned triangles
//////////////////////////////////////////////////////////////////////
bool PARTITIONED_TET_MESH::readClonedTriangles()
{
  string clonesFilename = _partitionPath;
  clonesFilename += string(".cloned.triangles");

  FILE* file = fopen(clonesFilename.c_str(), "rb");

  if (file == NULL) return false;

  // allocate the pointers
  if (_clonedTriangles != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedTriangles[x];
    delete[] _clonedTriangles;
  }
  _clonedTriangles = new vector<TRIANGLE*>*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _clonedTriangles[x] = new vector<TRIANGLE*>[_partitions];

  // allocate cloned triangle IDs 
  if (_clonedTriangleIDs != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedTriangleIDs[x];
    delete[] _clonedTriangleIDs;
  }
  _clonedTriangleIDs = new vector<int>*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _clonedTriangleIDs[x] = new vector<int>[_partitions];

  // read out the IDs
  for (int x = 0; x < _partitions; x++)
  {
    vector<TRIANGLE*>& faces = _meshes[x]->explicitSurfaceFaces();

    for (int y = 0; y < _partitions; y++)
    {
      if (x == y) continue;

      // read in how many there are first
      int size;
      fread((void*)&size, sizeof(int), 1, file);

      // read in the actual indices second
      for (int z = 0; z < size; z++)
      {
        int ID;
        fread((void*)&ID, sizeof(int), 1, file);
        _clonedTriangleIDs[x][y].push_back(ID);
        _clonedTriangles[x][y].push_back(faces[ID]);
      }
    }
  }
  fclose(file);

  return true;
}

//////////////////////////////////////////////////////////////////////
// read/write cloned triangles
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::writeClonedTriangles()
{
  string clonesFilename = _partitionPath;
  clonesFilename += string(".cloned.triangles");

  FILE* file = fopen(clonesFilename.c_str(), "wb");

  if (file == NULL)
  {
    cout << " Failed to write to cloned triangle cache "
         << clonesFilename.c_str() << "!" << endl;
    return;
  }

  /*
  // allocate the pointers
  if (_clonedTriangles != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedTriangles[x];
    delete[] _clonedTriangles;
  }
  _clonedTriangles = new vector<TRIANGLE*>*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _clonedTriangles[x] = new vector<TRIANGLE*>[_partitions];

  // allocate cloned triangle IDs 
  if (_clonedTriangleIDs != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _clonedTriangleIDs[x];
    delete[] _clonedTriangleIDs;
  }
  _clonedTriangleIDs = new vector<int>*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _clonedTriangleIDs[x] = new vector<int>[_partitions];
    */

  // write out the IDs
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      if (x == y) continue;

      // write out how many there are first
      int size = _clonedTriangleIDs[x][y].size();
      fwrite((void*)&size, sizeof(int), 1, file);
      vector<int>& IDs = _clonedTriangleIDs[x][y];

      // write out the actual indices second
      for (unsigned int z = 0; z < IDs.size(); z++)
        fwrite((void*)&IDs[z], sizeof(int), 1, file);
    }

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Dump out the embedded mesh to PBRT
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::exportEmbeddingToPBRT(string renderPath, int frame)
{
  if(_embeddedMesh==NULL)
  {
    cout << " Need to run PartitionedEmbedder!" << endl;
    exit(0);
  }
   
  char buffer[256];
  sprintf(buffer, "%i", frame);
  std::string number = std::string(buffer);
  if (frame < 10) number = std::string("0") + number;
  if (frame < 100) number = std::string("0") + number;
  if (frame < 1000) number = std::string("0") + number;

  FILE* file = NULL;

  string filename = renderPath + string("frame.")  + number + string(".pbrt");
  file = fopen(filename.c_str(), "w");

  if (file == NULL)
  {
    cout << " Couldn't open file " << filename.c_str() << "!!!" << endl;
    exit(0);
  }

  cout << " Writing PBRT file " << filename.c_str() << endl;
  fprintf(file, "AttributeBegin \n Translate 0 0 0\n Shape ");
  fprintf(file, "\"trianglemesh\"\n");
  fprintf(file, "\"point P\" [\n");
  vector<VEC3F>& vertices = _embeddedMesh->vertices;
  for (unsigned int i = 0; i < vertices.size(); i++) {
    VEC3F vec = vertices[i];
    fprintf(file, "%f %f %f\n", vec[0], vec[1], vec[2]);
  }
  fprintf(file, "]\n");

  // smooth out the normals
  _embeddedMesh->SmoothVertexNormals();
  _embeddedMesh->SmoothVertexNormals();

  fprintf(file, "\"normal N\" [\n");
  vector<VEC3F>& normals = _embeddedMesh->normals;
  for (unsigned int i = 0; i < vertices.size(); i++) {
    VEC3F N = normals[i];
    fprintf(file, "%f %f %f\n", N[0], N[1], N[2]);
  }
  fprintf(file, "]\n");

  fprintf(file, "\"integer indices\" [\n");

  vector<OBJ::Face>& faces = _embeddedMesh->faces;
  for (unsigned int i = 0; i < faces.size(); i++) {
    OBJ::Face face = faces[i];
    fprintf(file, "%i %i %i\n", face.vertices[0], face.vertices[1], face.vertices[2]);
  }
  fprintf(file,"]\nAttributeEnd\n");
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Dump out the embedded mesh to an OBJ
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::writeFilteredOBJ(string renderPath, int frame)
{
  if(_embeddedMesh==NULL)
  {
    cout << " Need to run PartitionedEmbedder!" << endl;
    exit(0);
  }
   
  char buffer[256];
  sprintf(buffer, "%i", frame);
  std::string number = std::string(buffer);
  if (frame < 10) number = std::string("0") + number;
  if (frame < 100) number = std::string("0") + number;
  if (frame < 1000) number = std::string("0") + number;

  string filename = renderPath + string("frame.")  + number + string(".obj");

  // get the list of filtered faces to output to OBJ
  Real filteringThreshold = _embeddedMesh->filteringThreshold();
  vector<OBJ::Face>& faces = _embeddedMesh->faces;
  vector<float>& triangleAreas = _embeddedMesh->triangleAreas();
  vector<float>& maxEdgeLengths = _embeddedMesh->maxEdgeLengths();
  vector<VEC3F>& vertices = _embeddedMesh->vertices;
  vector<int> filteredFaces;
  for (unsigned int i = 0; i < faces.size(); i++) 
  {
    OBJ::Face face = faces[i];

    VEC3F* triangleVertices[3];
    for (int y = 0; y < 3; y++)
      triangleVertices[y] = &vertices[faces[i].vertices[y]];
    TRIANGLE triangle(triangleVertices[0], triangleVertices[1], triangleVertices[2]);
    Real currentArea = triangle.area();
    Real currentMaxLength = triangle.maxEdgeLength();

    if (currentArea > triangleAreas[i] * filteringThreshold)
      continue;
    if (currentMaxLength > maxEdgeLengths[i] * filteringThreshold)
      continue;

    filteredFaces.push_back(i);
  }

  _embeddedMesh->SaveFiltered(filteredFaces, filename); 
}

//////////////////////////////////////////////////////////////////////
// Dump out the embedded mesh to PBRT
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::exportFilteredEmbeddingToPBRT(string renderPath, int frame)
{
  if(_embeddedMesh==NULL)
  {
    cout << " Need to run PartitionedEmbedder!" << endl;
    exit(0);
  }
   
  char buffer[256];
  sprintf(buffer, "%i", frame);
  std::string number = std::string(buffer);
  if (frame < 10) number = std::string("0") + number;
  if (frame < 100) number = std::string("0") + number;
  if (frame < 1000) number = std::string("0") + number;

  FILE* file = NULL;

  string filename = renderPath + string("frame.")  + number + string(".pbrt");
  file = fopen(filename.c_str(), "w");

  if (file == NULL)
  {
    cout << " Couldn't open file " << filename.c_str() << "!!!" << endl;
    exit(0);
  }

  cout << " Writing PBRT file " << filename.c_str() << endl;
  fprintf(file, "AttributeBegin \n Translate 0 0 0\n Shape ");
  fprintf(file, "\"trianglemesh\"\n");
  fprintf(file, "\"point P\" [\n");
  vector<VEC3F>& vertices = _embeddedMesh->vertices;
  for (unsigned int i = 0; i < vertices.size(); i++) {
    VEC3F vec = vertices[i];
    fprintf(file, "%f %f %f\n", vec[0], vec[1], vec[2]);
  }
  fprintf(file, "]\n");
  fprintf(file, "\"integer indices\" [\n");

  Real filteringThreshold = _embeddedMesh->filteringThreshold();
  vector<OBJ::Face>& faces = _embeddedMesh->faces;
  vector<float>& triangleAreas = _embeddedMesh->triangleAreas();
  vector<float>& maxEdgeLengths = _embeddedMesh->maxEdgeLengths();
  for (unsigned int i = 0; i < faces.size(); i++) 
  {
    OBJ::Face face = faces[i];

    VEC3F* triangleVertices[3];
    for (int y = 0; y < 3; y++)
      triangleVertices[y] = &vertices[faces[i].vertices[y]];
    TRIANGLE triangle(triangleVertices[0], triangleVertices[1], triangleVertices[2]);
    Real currentArea = triangle.area();
    Real currentMaxLength = triangle.maxEdgeLength();

    if (currentArea > triangleAreas[i] * filteringThreshold)
      continue;
    if (currentMaxLength > maxEdgeLengths[i] * filteringThreshold)
      continue;

    fprintf(file, "%i %i %i\n", face.vertices[0], face.vertices[1], face.vertices[2]);
  }
  fprintf(file,"]\nAttributeEnd\n");
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// recompute the masses, dividing mass evenly between cloned nodes
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::recomputePartitionedMass()
{
  Real trueMass = _originalMesh->mass(0);
  Real summedMass = 0;
  for (int x = 0; x < _partitions; x++)
  {
    _meshes[x]->resetMasses(trueMass);

    SPARSE_MATRIX& massMatrix = _meshes[x]->TET_MESH::massMatrix();
    int totalNodes = _meshes[x]->unconstrainedNodes();

    map<int, int> clonedTimes;
    for (int y = 0; y < _partitions; y++)
    {
      if (x == y) continue;
      vector<pair<int, int> >& clones = _clonedVertices[x][y];
      if (clones.size() == 0) continue;

      for (unsigned int z = 0; z < clones.size(); z++)
        clonedTimes[clones[z].first]++;
    }

    map<int, int>::iterator iter;
    for (iter = clonedTimes.begin(); iter != clonedTimes.end(); iter++)
    {
      int index = 3 * iter->first;
      Real invValence = 1.0 / (iter->second + 1);
      if (index >= 3 * totalNodes) continue;

      massMatrix(index, index) *= invValence;
      massMatrix(index + 1, index + 1) *= invValence;
      massMatrix(index + 2, index + 2) *= invValence;
    }

    _meshes[x]->resetTotalMass();

    summedMass += _meshes[x]->totalMass();

    // recompute centers of mass
    _meshes[x]->computeCenterOfMass();
    _meshes[x]->computeRestCenterOfMass();
      
    cout << " Total mass of partition " << x << ": " << _meshes[x]->totalMass() << endl;
  }
  cout << " Total mass of partitioned mesh: " << summedMass << endl;
}

//////////////////////////////////////////////////////////////////////
// stream read all tet mesh data
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::readState(FILE* file)
{
  for (int x = 0; x < _partitions; x++)
    _meshes[x]->readMeshFile(file, 1.0, false);
}

//////////////////////////////////////////////////////////////////////
// stream write all tet mesh data
//////////////////////////////////////////////////////////////////////
void PARTITIONED_TET_MESH::writeState(FILE* file)
{
  for (int x = 0; x < _partitions; x++)
    _meshes[x]->writeMeshFile(file);
}
