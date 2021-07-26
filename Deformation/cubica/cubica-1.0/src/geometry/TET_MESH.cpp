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

#include <map>
#include <list>
#include "TET_MESH.h"
#include <INVERTIBLE.h>
#ifdef USING_OPENMP
#include <omp.h>
#endif

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
TET_MESH::TET_MESH(const char* filename, MATERIAL** materials, int totalMaterials,
                   bool simulate, Real scale, bool volumeScale) :
  _filename(filename), 
  _totalMaterials(totalMaterials), 
  _meshScaling(scale),
  _embeddedMesh(NULL) 
{
  readMeshFile(filename, scale);

  cout << " Initializing ... " << endl;
  init(simulate, true);

  _materials = new MATERIAL*[_totalMaterials];
  for (int x = 0; x < _totalMaterials; x++)
    _materials[x] = materials[x];

  // allocate a thread-safe copy of material for each thread
  _materialCopies = new MATERIAL**[_totalCores];
  for (int x = 0; x < _totalCores; x++)
  {
    _materialCopies[x] = new MATERIAL*[_totalMaterials];
    for (int y = 0; y < _totalMaterials; y++)
      _materialCopies[x][y] = _materials[y]->copy();
  }

  // if there is more than one material type, look for a material file
  if (totalMaterials > 1)
    readMaterials();

  // Compute tet scales if we are scaling forces, etc. by volume.
  if ( volumeScale )
  {
    _tetScales.resize( _tets.size() );

    for (unsigned int i = 0; i < _tets.size(); i++ )
    {
      _tetScales[ i ] = _tets[ i ].volume();
    }
  }
  else
  {
    _tetScales.resize( _tets.size(), 1.0 );
  }

  cout << " Using " << _totalCores << " cores. " << endl;

  _totalMass = _masses.sum() / 3.0;
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
TET_MESH::~TET_MESH()
{
  // clean up the thread safe copies
  for (int x = 0; x < _totalCores; x++)
  {
    for (int y = 0; y < _totalMaterials; y++)
      delete _materialCopies[x][y];
    delete[] _materialCopies[x];
  }
  delete[] _materialCopies;

  /*
  // clean up the materials
  for (int x = 0; x < _totalMaterials; x++)
    delete _materials[x];
  */
  delete[] _materials;

  delete[] _stiffnessCopies;
  delete[] _RCopies;

  for (unsigned int x = 0; x < _allFaces.size(); x++)
    delete _allFaces[x];

  for (unsigned int x = 0; x < _explicitSurfaceFaces.size(); x++)
    delete _explicitSurfaceFaces[x];
}

//////////////////////////////////////////////////////////////////////
// One-ring Constructor
// this creates a copy of the one-ring of the passed in vertex and TET_MESH,
// constraining all other vertices aside from the one in the middle
//////////////////////////////////////////////////////////////////////
TET_MESH::TET_MESH(VEC3F* point, TET_MESH* mesh,
                   map<VEC3F*, int>& originalToCopy,
                   map<string, double>& timingBreakdown) :
  _filename(""),
  _totalMaterials(mesh->_totalMaterials)
{
  // retrieve the one ring
  vector<TET*> ringTets;
  mesh->oneRing(point, ringTets);

  // add the rest state of the one unconstrained node first
  int pointID = mesh->vertexID(point);
  VEC3F* restPoint = mesh->restVertices(pointID);
  _vertices.push_back(VEC3F(*restPoint));
  originalToCopy[point] = 0;
  
  // clone all the vertices first
  for (unsigned int x = 0; x < ringTets.size(); x++)
    for (int y = 0; y < 4; y++)
    {
      VEC3F* originalVertex = ringTets[x]->vertices[y];
      int originalID = mesh->vertexID(originalVertex);
      VEC3F* originalRestVertex = mesh->restVertices(originalID);

      // if we haven't seen this one yet, make a copy
      if (originalToCopy.find(originalVertex) == originalToCopy.end())
      {
        _vertices.push_back(VEC3F(*originalRestVertex));
        originalToCopy[originalVertex] = _vertices.size() - 1;
      }
    }
  
  // clone all the tets after the vertex addresses have been nailed
  // down in the vector
  for (unsigned int x = 0; x < ringTets.size(); x++)
  {
    // clone the constrained vertices
    VEC3F* newVerts[4];
    for (int y = 0; y < 4; y++)
    {
      VEC3F* originalVertex = ringTets[x]->vertices[y];
      // record the vertex for the tet clone
      newVerts[y] = &(_vertices[originalToCopy[originalVertex]]);
    }

    // clone the tet
    _tets.push_back(TET(newVerts));
  }
  // copy the material type
  for (unsigned int x = 0; x < ringTets.size(); x++)
    _tets[x].materialIndex() = ringTets[x]->materialIndex();

  // copy the skinning matrices (if any)
  for (unsigned int x = 0; x < ringTets.size(); x++)
    for (int y = 0; y < 4; y++)
      _tets[x].skinningMatrix(y) = ringTets[x]->skinningMatrix(y);
  
  // init sizes
  _unconstrainedSize = 1;
  _constrainedSize = _vertices.size() - 1;

  // initialize the mesh
  init(true);

  // allocate copies of the materials
  _materials = new MATERIAL*[_totalMaterials];
  for (int x = 0; x < _totalMaterials; x++)
    _materials[x] = mesh->_materials[x]->copy();

  // allocate a thread-safe copy of material for each thread
  _materialCopies = new MATERIAL**[_totalCores];
  for (int x = 0; x < _totalCores; x++)
  {
    _materialCopies[x] = new MATERIAL*[_totalMaterials];
    for (int y = 0; y < _totalMaterials; y++)
      _materialCopies[x][y] = _materials[y]->copy();
  }

  // fix up the mass matrix
  SPARSE_MATRIX& massMatrix = mesh->massMatrix();
  _masses.clear();
  _masses(0,0) = massMatrix(3 * pointID, 3 * pointID);
  _masses(1,1) = massMatrix(3 * pointID + 1, 3 * pointID + 1);
  _masses(2,2) = massMatrix(3 * pointID + 2, 3 * pointID + 2);
}

//////////////////////////////////////////////////////////////////////
// Initialize mesh after the vertices have been specified
//////////////////////////////////////////////////////////////////////
void TET_MESH::init(bool simulate, bool verbose)
{
  // If we are simulating, allocate the giant matrices.
  // By default SUBSPACE_TET_MESH does not allocate these
  // since allocates its own smaller versions
  if (simulate)
  {
    // initialize Newmark variables
    int rank = _unconstrainedSize * 3;
    int totalTets = _tets.size();
    _x.resizeAndWipe(rank);
    _F.resizeAndWipe(totalTets * 9);
    _R.resizeAndWipe(rank);
  }
  
  if (verbose) 
  {
    cout << "Recording constraint data ...";
    flush(cout);
  }

  // record which nodes are constrained
  for (unsigned int x = 0; x < _vertices.size(); x++)
    if ((int)x < _unconstrainedSize)
      _constrained[&_vertices[x]] = false;
    else
      _constrained[&_vertices[x]] = true;

  if (verbose) 
  {
    cout << "done." << endl;
    cout << "Recording vertex ID data ...";
    flush(cout);
  }

  // record the vertex positions of each vertex
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertexID[&_vertices[x]] = x;

  if (verbose) 
  {
    cout << "done." << endl;
    cout << "Recording rest pose data ...";
    flush(cout);
  }

  // cache a copy of all the vertices for the rest pose
  _restPose.resize(_vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
    //_restPose.push_back(_vertices[x]);
    _restPose[x] = _vertices[x];

  if (verbose) 
  {
    cout << "done." << endl;
    cout << "Recording tet membership ...";
    flush(cout);
  }

  // record what tets each vertex belongs to
  for (unsigned int x = 0; x < _tets.size(); x++)
    for (int y = 0; y < 4; y++)
      _tetMembership[_tets[x].vertices[y]].push_back(x);
  
  if (verbose) 
  {
    cout << "done." << endl;
    cout << "Initializing internal forces ...";
    flush(cout);
  }

  // initialize internal forces to zero
  VEC3F zero;
  _internalForces.resize(_unconstrainedSize);
  for (int x = 0; x < _unconstrainedSize; x++)
    //_internalForces.push_back(zero);
    _internalForces[x] = zero;

  if (verbose) 
  {
    cout << "done." << endl;
    cout << "Computing surface faces ... "; flush(cout);
  }
  computeSurfaceFaces();
  computeExplicitSurfaceFaces();
  cout << "done." << endl;
  //cout << "Computing *all* tet faces ..."; flush(cout);
  //computeAllFaces();
  if (verbose) 
  {
    //cout << "done." << endl;
    cout << "Computing surface vertices ... ";
    flush(cout);
  }
  computeSurfaceVertices();
  if (verbose) 
  {
    cout << "done." << endl;
    flush(cout);
  }

  // create surface vertex hash
  for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
    _vertexOnSurface[_surfaceVertices[x]] = true;
  
  if (verbose) 
  {
    cout << "Computing masses ... ";
    flush(cout);
  }
  computeMasses();
  if (verbose) 
  {
    cout << " done. " << endl;
    flush(cout);
  }

  // compute surface area of each surface node
  computeSurfaceAreas();

  // cache center of mass
  computeCenterOfMass();
  computeRestCenterOfMass();

  // allocate OpenMP copies
#ifdef USING_OPENMP
  //omp_set_num_threads(1);
  _totalCores = omp_get_max_threads();
#else
  _totalCores = 1; 
#endif

  // allocate a thread safe internal force vector for each thread
  _RCopies = new VECTOR[_totalCores];
 
  // allocate a thread safe sparse stiffness matrix for each thread
  _stiffnessCopies = new SPARSE_MATRIX[_totalCores];

  // compute the deformed to rest mapping
  _deformedToRest.clear();
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _deformedToRest[&(_vertices[x])] = &(_restPose[x]);

  // tell all the tets about the deformed to rest mapping
  for (unsigned int x = 0; x < _tets.size(); x++)
    _tets[x].deformedToRest = &_deformedToRest;

  // compute the address to index tet table
  for (unsigned int x = 0; x < _tets.size(); x++)
    _tetID[&(_tets[x])] = x;

  // compute the initial inertia tensor
  _inertiaTensor.resizeAndWipe(3,3);
  refreshInertiaTensor();

  _inertiaTensorOriginal = _inertiaTensor;

  _inertiaTensorDt.resizeAndWipe(3,3);
}

//////////////////////////////////////////////////////////////////////
// Reinitialize the mesh after rescaling
//////////////////////////////////////////////////////////////////////
void TET_MESH::reinitialize()
{
  // reinitialize the traction vectors for all the tets
  for (unsigned int x = 0; x < _tets.size(); x++)
    _tets[x].init();

  // cache a copy of all the vertices for the rest pose
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _restPose[x] = _vertices[x];

  // compute surface area of each surface node
  computeSurfaceAreas();

  // cache center of mass
  computeCenterOfMass();
  computeRestCenterOfMass();
}

//////////////////////////////////////////////////////////////////////
// read in the material types
//////////////////////////////////////////////////////////////////////
void TET_MESH::readMaterials()
{
  // construct the materials filename
  string filename = string(_filename);
  filename += string(".materials");
  FILE* file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " WARNING: No material distribution file found." << endl; 
    return;
  }
  
  int materialIndex;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    fread((void*)&materialIndex, sizeof(int), 1, file);
    _tets[x].materialIndex() = materialIndex;
  }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// write out the material types
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeMaterials()
{
  // construct the materials filename
  string filename = string(_filename);
  filename += string(".materials");
  FILE* file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " WARNING: Could not write out materials file!" << endl; 
    return;
  }
  
  int materialIndex;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    materialIndex = _tets[x].materialIndex();
    fwrite((void*)&materialIndex, sizeof(int), 1, file);
  }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Stream in a tet mesh file from an already open file
//////////////////////////////////////////////////////////////////////
void TET_MESH::readMeshFile(FILE* file, Real scale, bool verbose)
{
  // read in vertex array sizes
  fread((void*)&_unconstrainedSize, sizeof(int), 1, file);
  fread((void*)&_constrainedSize, sizeof(int), 1, file);

  if (verbose)
  {
    printf(" %i unconstrained nodes\n", _unconstrainedSize);
    printf(" %i constrained nodes\n", _constrainedSize);
  }

  _vertices.resize(_unconstrainedSize + _constrainedSize);
  if (verbose)
    cout << " Vertices size: " << (_unconstrainedSize + _constrainedSize) * sizeof(VEC3F) / pow(2.0, 20.0) << " MB " << endl;

  // read in vertex positions
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F node;
#ifdef SINGLE_PRECISION
    double nodeDouble[3];
    fread((void*)nodeDouble, sizeof(double), 3, file);
    node[0] = (float)nodeDouble[0];
    node[1] = (float)nodeDouble[1];
    node[2] = (float)nodeDouble[2];
#else
    fread((void*)&(node), sizeof(Real), 3, file);
#endif
    node *= scale;

    _vertices[x] = node;
  }
  for (int x = 0; x < _constrainedSize; x++)
  {
    VEC3F node;
#ifdef SINGLE_PRECISION
    double nodeDouble[3];
    fread((void*)nodeDouble, sizeof(double), 3, file);
    node[0] = (float)nodeDouble[0];
    node[1] = (float)nodeDouble[1];
    node[2] = (float)nodeDouble[2];
#else
    fread((void*)&(node), sizeof(Real), 3, file);
#endif
    node *= scale;

    _vertices[x + _unconstrainedSize] = node;
  }

  // read in tet vertex lists
  int totalTets;
  fread((void*)&totalTets, sizeof(int), 1, file);

  if (verbose)
  {
    printf(" %i total tets\n", totalTets);
    cout << " Tets size: " << totalTets * sizeof(TET) / pow(2.0, 20.0) << " MB " << endl;
  }
  _tets.resize(totalTets);
  for (int x = 0; x < totalTets; x++)
  {
    // for each node in the tet
    VEC3F* vertices[4];
    int indices[4];
    fread((void*)&indices[0], sizeof(int), 4, file);
    vertices[0] = &(_vertices[indices[0]]);
    vertices[1] = &(_vertices[indices[1]]);
    vertices[2] = &(_vertices[indices[2]]);
    vertices[3] = &(_vertices[indices[3]]);
    TET tet(vertices);
    _tets[x] = tet;

    if (totalTets > 1 && verbose)
      if (x % (int)(totalTets / 10) == 0)
      {
        cout << 100 * ((Real)x / totalTets) << "% ";
        flush(cout);
      }
  }
}

//////////////////////////////////////////////////////////////////////
// Read in a tet mesh file
//////////////////////////////////////////////////////////////////////
void TET_MESH::readMeshFile(const char* filename, Real scale)
{
  FILE* file = fopen(filename, "rb");
  
  if (file == NULL)
    printf("Filename %s not found!\n", filename);

  readMeshFile(file, scale, true);

  /*
  // read in vertex array sizes
  fread((void*)&_unconstrainedSize, sizeof(int), 1, file);
  fread((void*)&_constrainedSize, sizeof(int), 1, file);

  printf(" %i unconstrained nodes\n", _unconstrainedSize);
  printf(" %i constrained nodes\n", _constrainedSize);

  _vertices.resize(_unconstrainedSize + _constrainedSize);
  cout << " Vertices size: " << (_unconstrainedSize + _constrainedSize) * sizeof(VEC3F) / pow(2.0, 20.0) << " MB " << endl;

  // read in vertex positions
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F node;
#ifdef SINGLE_PRECISION
    double nodeDouble[3];
    fread((void*)nodeDouble, sizeof(double), 3, file);
    node[0] = (float)nodeDouble[0];
    node[1] = (float)nodeDouble[1];
    node[2] = (float)nodeDouble[2];
#else
    fread((void*)&(node), sizeof(Real), 3, file);
#endif
    node *= scale;

    //_vertices.push_back(node);
    _vertices[x] = node;
  }
  for (int x = 0; x < _constrainedSize; x++)
  {
    VEC3F node;
#ifdef SINGLE_PRECISION
    double nodeDouble[3];
    fread((void*)nodeDouble, sizeof(double), 3, file);
    node[0] = (float)nodeDouble[0];
    node[1] = (float)nodeDouble[1];
    node[2] = (float)nodeDouble[2];
#else
    fread((void*)&(node), sizeof(Real), 3, file);
#endif
    node *= scale;

    //_vertices.push_back(node);
    _vertices[x + _unconstrainedSize] = node;
  }

  // read in tet vertex lists
  int totalTets;
  fread((void*)&totalTets, sizeof(int), 1, file);

  printf(" %i total tets\n", totalTets);
  cout << " Tets size: " << totalTets * sizeof(TET) / pow(2.0, 20.0) << " MB " << endl;
  _tets.resize(totalTets);
  for (int x = 0; x < totalTets; x++)
  {
    // for each node in the tet
    VEC3F* vertices[4];
    int indices[4];
    fread((void*)&indices[0], sizeof(int), 4, file);
    vertices[0] = &(_vertices[indices[0]]);
    vertices[1] = &(_vertices[indices[1]]);
    vertices[2] = &(_vertices[indices[2]]);
    vertices[3] = &(_vertices[indices[3]]);
    TET tet(vertices);
    _tets[x] = tet;

    if (totalTets > 1)
      if (x % (int)(totalTets / 10) == 0)
      {
        cout << 100 * ((Real)x / totalTets) << "% ";
        flush(cout);
      }
  }
  */

  fclose(file);
  cout << " Done reading mesh file. " << endl;
  flush(cout);

  Real minX = _vertices[0][0];
  Real maxX = _vertices[0][0];
  Real minY = _vertices[0][1];
  Real maxY = _vertices[0][1];
  Real minZ = _vertices[0][2];
  Real maxZ = _vertices[0][2];
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    if (_vertices[x][0] < minX) minX = _vertices[x][0];
    if (_vertices[x][0] > maxX) maxX = _vertices[x][0];
    if (_vertices[x][1] < minY) minY = _vertices[x][1];
    if (_vertices[x][1] > maxY) maxY = _vertices[x][1];
    if (_vertices[x][2] < minZ) minZ = _vertices[x][2];
    if (_vertices[x][2] > maxZ) maxZ = _vertices[x][2];
  }

  cout << "Bounding box: " << endl
       << "(" << minX << ", " << minY << ", " << minZ << ")" << endl
       << "(" << maxX << ", " << maxY << ", " << maxZ << ")" << endl;
}

//////////////////////////////////////////////////////////////////////
// Stream out a mesh file to an already open file
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeMeshFile(FILE* file)
{
  // write out vertex array sizes
  fwrite((void*)&_unconstrainedSize, sizeof(int), 1, file);
  fwrite((void*)&_constrainedSize, sizeof(int), 1, file);

  // write out vertex positions
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    double nodeDouble[3];
    nodeDouble[0] = _vertices[x][0];
    nodeDouble[1] = _vertices[x][1];
    nodeDouble[2] = _vertices[x][2];
    fwrite((void*)nodeDouble, sizeof(double), 3, file);
  }

  // write out tet vertex lists
  int totalTets = _tets.size();
  fwrite((void*)&totalTets, sizeof(int), 1, file);

  for (int x = 0; x < totalTets; x++)
  {
    // for each node in the tet
    int indices[4];
    indices[0] = _vertexID[_tets[x].vertices[0]];
    indices[1] = _vertexID[_tets[x].vertices[1]];
    indices[2] = _vertexID[_tets[x].vertices[2]];
    indices[3] = _vertexID[_tets[x].vertices[3]];
    fwrite((void*)&indices[0], sizeof(int), 4, file);
  }
}

//////////////////////////////////////////////////////////////////////
// Write out a mesh file
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeMeshFile(const char* filename)
{
  FILE* file = fopen(filename, "wb");
  
  if (file == NULL)
    printf("Filename %s could not be opened!\n", filename);

  cout << " Writing file " << filename << " ... ";

  writeMeshFile(file);
  /*
  // write out vertex array sizes
  fwrite((void*)&_unconstrainedSize, sizeof(int), 1, file);
  fwrite((void*)&_constrainedSize, sizeof(int), 1, file);

  // write out vertex positions
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    double nodeDouble[3];
    nodeDouble[0] = _vertices[x][0];
    nodeDouble[1] = _vertices[x][1];
    nodeDouble[2] = _vertices[x][2];
    fwrite((void*)nodeDouble, sizeof(double), 3, file);
  }

  // write out tet vertex lists
  int totalTets = _tets.size();
  fwrite((void*)&totalTets, sizeof(int), 1, file);

  for (int x = 0; x < totalTets; x++)
  {
    // for each node in the tet
    int indices[4];
    indices[0] = _vertexID[_tets[x].vertices[0]];
    indices[1] = _vertexID[_tets[x].vertices[1]];
    indices[2] = _vertexID[_tets[x].vertices[2]];
    indices[3] = _vertexID[_tets[x].vertices[3]];
    fwrite((void*)&indices[0], sizeof(int), 4, file);
  }
  */

  fclose(file);
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// Write out a single tet mesh file (for debugging)
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeSingleTetMesh(const char* filename)
{
  FILE* file = fopen(filename, "wb");
  
  if (file == NULL)
    printf("Filename %s could not be opened!\n", filename);

  cout << " Writing single tet file " << filename << " ... ";

  // write out vertex array sizes
  int unconstrainedSize = 4;
  int constrainedSize = 0;
  fwrite((void*)&unconstrainedSize, sizeof(int), 1, file);
  fwrite((void*)&constrainedSize, sizeof(int), 1, file);

  VEC3F vertices[4];
  /*
  vertices[0][0] = 0;
  vertices[0][1] = 1;
  vertices[0][2] = -0.5;

  vertices[1][0] = 1;
  vertices[1][1] = 0;
  vertices[1][2] = -0.5;

  vertices[2][0] = 0;
  vertices[2][1] = 0;
  vertices[2][2] = 0;

  vertices[3][0] = 0;
  vertices[3][1] = 0;
  vertices[3][2] = -1;
  */
  vertices[0][0] = 0.5;
  vertices[0][1] = sqrt(3.0) / 6;
  vertices[0][2] = sqrt(2.0 / 3.0);

  vertices[1][0] = 0.5;
  vertices[1][1] = sqrt(3.0) / 2.0;
  vertices[1][2] = 0;

  vertices[2][0] = 1.0;
  vertices[2][1] = 0;
  vertices[2][2] = 0;

  vertices[3][0] = 0;
  vertices[3][1] = 0;
  vertices[3][2] = 0;
  
  // write out vertex positions
  for (unsigned int x = 0; x < 4; x++)
  {
    double nodeDouble[3];
    nodeDouble[0] = vertices[x][0];
    nodeDouble[1] = vertices[x][1];
    nodeDouble[2] = vertices[x][2];
    fwrite((void*)nodeDouble, sizeof(double), 3, file);
  }

  // write out tet vertex lists
  int totalTets = 1;
  fwrite((void*)&totalTets, sizeof(int), 1, file);

  for (int x = 0; x < totalTets; x++)
  {
    // for each node in the tet
    int indices[4];
    indices[0] = 0;
    indices[1] = 1;
    indices[2] = 2;
    indices[3] = 3;
    fwrite((void*)&indices[0], sizeof(int), 4, file);
  }

  fclose(file);
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets, including internal faces
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawAllTets()
{
  // colors for the different partitions
  float colors[8][4];

  // cache the colors
  colors[0][0] = 1.0f; colors[0][1] = 0.0f; colors[0][2] = 0.0f; colors[0][3] = 1.0f;
  colors[1][0] = 0.0f; colors[1][1] = 1.0f; colors[1][2] = 0.0f; colors[1][3] = 1.0f;
  colors[2][0] = 0.0f; colors[2][1] = 0.0f; colors[2][2] = 1.0f; colors[2][3] = 1.0f;
  colors[3][0] = 1.0f; colors[3][1] = 0.0f; colors[3][2] = 1.0f; colors[3][3] = 1.0f;
  colors[4][0] = 0.0f; colors[4][1] = 1.0f; colors[4][2] = 1.0f; colors[4][3] = 1.0f;
  colors[5][0] = 1.0f; colors[5][1] = 1.0f; colors[5][2] = 0.0f; colors[5][3] = 1.0f;
  colors[6][0] = 1.0f; colors[6][1] = 1.0f; colors[6][2] = 1.0f; colors[6][3] = 1.0f;
  colors[7][0] = 0.5f; colors[7][1] = 0.5f; colors[7][2] = 0.5f; colors[7][3] = 1.0f;

  srand(123456);
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    int randColor = 8 * ((Real)rand() / RAND_MAX);
    int i = randColor % 8;
    glColor4f(colors[i][0], colors[i][1], colors[i][2], 1.0);
    _tets[x].drawTriangles();
  }
}

//////////////////////////////////////////////////////////////////////
// Draw all the tets, including internal faces
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawTetMaterials(Real z)
{
  // colors for the different partitions
  float colors[8][4];

  // cache the colors
  colors[0][0] = 1.0f; colors[0][1] = 0.0f; colors[0][2] = 0.0f; colors[0][3] = 1.0f;
  colors[1][0] = 0.0f; colors[1][1] = 1.0f; colors[1][2] = 0.0f; colors[1][3] = 1.0f;
  colors[2][0] = 0.0f; colors[2][1] = 0.0f; colors[2][2] = 1.0f; colors[2][3] = 1.0f;
  colors[3][0] = 1.0f; colors[3][1] = 0.0f; colors[3][2] = 1.0f; colors[3][3] = 1.0f;
  colors[4][0] = 0.0f; colors[4][1] = 1.0f; colors[4][2] = 1.0f; colors[4][3] = 1.0f;
  colors[5][0] = 1.0f; colors[5][1] = 1.0f; colors[5][2] = 0.0f; colors[5][3] = 1.0f;
  colors[6][0] = 1.0f; colors[6][1] = 1.0f; colors[6][2] = 1.0f; colors[6][3] = 1.0f;
  colors[7][0] = 0.5f; colors[7][1] = 0.5f; colors[7][2] = 0.5f; colors[7][3] = 1.0f;

  /*
  for (int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE triangle = _tets[_surfaceFaces[x].first].face(_surfaceFaces[x].second);

    float color[4];
    int mod = _tets[_surfaceFaces[x].first].materialIndex() % 8;
    int div = _tets[_surfaceFaces[x].first].materialIndex() / 8;
    color[0] = colors[mod][0];
    color[1] = colors[mod][1];
    color[2] = colors[mod][2];
    color[3] = 1.0f;
    glColor4f(color[0], color[1], color[2], color[3]);
    triangle.draw();
  }
  */
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    /*
    if ((*(_tets[x].vertices[0]))[2] < z ||
        (*(_tets[x].vertices[1]))[2] < z ||
        (*(_tets[x].vertices[2]))[2] < z ||
        (*(_tets[x].vertices[3]))[2] < z)
      continue;
    */

    // DEBUG
    if (_tets[x].materialIndex() == 0) continue;

    float color[4];
    int mod = _tets[x].materialIndex() % 8;
    //int div = _tets[x].materialIndex() / 8;
    color[0] = colors[mod][0];
    color[1] = colors[mod][1];
    color[2] = colors[mod][2];
    color[3] = 1.0f;
    glColor4f(color[0], color[1], color[2], color[3]);
    _tets[x].drawFaces();
  }
}

//////////////////////////////////////////////////////////////////////
// Draw only the surface faces
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawSurfaceFaces()
{
  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE triangle = _tets[_surfaceFaces[x].first].face(_surfaceFaces[x].second);
    triangle.draw();
    //_surfaceFaces[x]->draw();
  }
}

//////////////////////////////////////////////////////////////////////
// Draw only the surface faces of the rest pose
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawRestPose()
{
  VECTOR backup = _x;
  _x *= 0;
  TET_MESH::updateFullMesh();
  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE triangle = _tets[_surfaceFaces[x].first].face(_surfaceFaces[x].second);
    triangle.draw();
  }
  _x = backup;
  TET_MESH::updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
// Draw the tets up to a z slice
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawZSlice(Real z)
{
  // colors for the different partitions
  float colors[8][4];

  // cache the colors
  colors[0][0] = 1.0f; colors[0][1] = 0.0f; colors[0][2] = 0.0f; colors[0][3] = 1.0f;
  colors[1][0] = 0.0f; colors[1][1] = 1.0f; colors[1][2] = 0.0f; colors[1][3] = 1.0f;
  colors[2][0] = 0.0f; colors[2][1] = 0.0f; colors[2][2] = 1.0f; colors[2][3] = 1.0f;
  colors[3][0] = 1.0f; colors[3][1] = 0.0f; colors[3][2] = 1.0f; colors[3][3] = 1.0f;
  colors[4][0] = 0.0f; colors[4][1] = 1.0f; colors[4][2] = 1.0f; colors[4][3] = 1.0f;
  colors[5][0] = 1.0f; colors[5][1] = 1.0f; colors[5][2] = 0.0f; colors[5][3] = 1.0f;
  colors[6][0] = 1.0f; colors[6][1] = 1.0f; colors[6][2] = 1.0f; colors[6][3] = 1.0f;
  colors[7][0] = 0.5f; colors[7][1] = 0.5f; colors[7][2] = 0.5f; colors[7][3] = 1.0f;

  srand(123456);
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    bool visible = true;
    int randColor = 7 * ((Real)rand() / RAND_MAX);
    int i = randColor % 8;
    glColor4f(colors[i][0], colors[i][1], colors[i][2], 1.0);
    //glColor4f(1.0, 1.0, 1.0, 1.0);
    for (int y = 0; y < 4; y++)
    {
      if ((*_tets[x].vertices[0])[2] > z)
        visible = false;
      if ((*_tets[x].vertices[1])[2] > z)
        visible = false;
      if ((*_tets[x].vertices[2])[2] > z)
        visible = false;
      if ((*_tets[x].vertices[3])[2] > z)
        visible = false;
    }
    if (visible)
      _tets[x].drawTriangles();
  }
}

//////////////////////////////////////////////////////////////////////
// Draw the unconstrained nodes as GL points
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawUnconstrainedNodes()
{
  glBegin(GL_POINTS);
  for (int x = 0; x < _unconstrainedSize; x++)
#ifdef SINGLE_PRECISION    
    glVertex3fv(_vertices[x]);
#else
    glVertex3dv(_vertices[x]);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw the constrained nodes as GL points
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawConstrainedNodes()
{
  glBegin(GL_POINTS);
  for (unsigned int x = _unconstrainedSize; x < _vertices.size(); x++)
#ifdef SINGLE_PRECISION    
    glVertex3fv(_vertices[x]);
#else
    glVertex3dv(_vertices[x]);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw the collision nodes as GL points
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawCollisionNodes()
{
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < _collisionNodes.size(); x++)
    if ((*(_collisionNodes[x]))[2] > 0.4)
#ifdef SINGLE_PRECISION    
      glVertex3fv(*(_collisionNodes[x]));
#else
      glVertex3dv(*(_collisionNodes[x]));
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw the surface nodes as GL points
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawSurfaceVertices()
{
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
#ifdef SINGLE_PRECISION    
    glVertex3fv(*_surfaceVertices[x]);
#else
    glVertex3dv(*_surfaceVertices[x]);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw some arbitrary list of points
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawPoints(vector<VEC3F*>& points)
{
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < points.size(); x++)
#ifdef SINGLE_PRECISION    
    glVertex3fv(*points[x]);
#else
    glVertex3dv(*points[x]);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw just one partition
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawPartition(int partition)
{
  for (unsigned int x = 0; x < _tets.size(); x++)
    if (_tets[x].partition == partition)
      _tets[x].drawTriangles();
}

//////////////////////////////////////////////////////////////////////
// Draw all the partitions, but color coded
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawColoredPartitions()
{
  // load up some reasonable colors
  float colors[8][4];
  colors[0][0] = 1.0f; colors[0][1] = 0.0f; colors[0][2] = 0.0f; colors[0][3] = 1.0f;
  colors[1][0] = 0.0f; colors[1][1] = 1.0f; colors[1][2] = 0.0f; colors[1][3] = 1.0f;
  colors[2][0] = 0.0f; colors[2][1] = 0.0f; colors[2][2] = 1.0f; colors[2][3] = 1.0f;
  colors[3][0] = 1.0f; colors[3][1] = 0.0f; colors[3][2] = 1.0f; colors[3][3] = 1.0f;
  colors[4][0] = 0.0f; colors[4][1] = 1.0f; colors[4][2] = 1.0f; colors[4][3] = 1.0f;
  colors[5][0] = 1.0f; colors[5][1] = 1.0f; colors[5][2] = 0.0f; colors[5][3] = 1.0f;
  colors[6][0] = 1.0f; colors[6][1] = 1.0f; colors[6][2] = 1.0f; colors[6][3] = 1.0f;
  colors[7][0] = 0.5f; colors[7][1] = 0.5f; colors[7][2] = 0.5f; colors[7][3] = 1.0f;

  // if there's more than 8 partitions just do darker versions
  // of the colors we already have
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    float color[4];
    int mod = _tets[x].partition % 8;
    int div = _tets[x].partition / 8;
    float fraction = 1.0f / (div + 1.0f);
    color[0] = fraction * colors[mod][0];
    color[1] = fraction * colors[mod][1];
    color[2] = fraction * colors[mod][2];
    color[3] = 1.0f;
    glColor4f(color[0], color[1], color[2], color[3]);
    _tets[x].drawTriangles();
  }
}

//////////////////////////////////////////////////////////////////////
// Compute the surface areas of the surface vertices
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceAreas()
{
  _surfaceArea.clear();

  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    // build the triangle area
    int tetID = _surfaceFaces[x].first;
    int triangleIndex = _surfaceFaces[x].second;
    TRIANGLE triangle = _tets[tetID].face(triangleIndex);
    Real area = triangle.area() / 3.0;

    for (int y = 0; y < 3; y++)
    {
      VEC3F* vertex = triangle.vertex(y);
      if (_surfaceArea.find(vertex) == _surfaceArea.end())
        _surfaceArea[vertex] = area;
      else
        _surfaceArea[vertex] += area;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Construct a vector of the triangle faces on the mesh surface
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceFaces()
{
  if (_filename.size() == 0) return;

  // try to read in a precomputed cache
  if (readSurfaceFaceCache()) 
  {
    cout << " surface cache found ... ";
    return;
  }
  
  // hash the faces according to the sum of their addresses
  //map<long long, vector<TRIANGLE*> > faceHash;
  map<long long, vector<pair<int,int> > > faceHash;

  // hash the tet and triangle index as well so that everything can be dumped
  // to a file after
  map<long long, vector<int> > tetHash;
  map<long long, vector<int> > triangleHash;
  
  // track the addresses generated, make them long longs to support
  // 64 bit architectures
  vector<long long> sums;

  // insert all the tet faces into the face hash
  for (unsigned int x = 0; x < _tets.size(); x++)
    for (int y = 0; y < 4; y++)
    {
      TRIANGLE face = _tets[x].face(y);

      // sum the addresses for use as the hash
      long long sum = 0;
      for (int z = 0; z < 3; z++)
        sum = sum + (long long)(face.vertex(z));

      // hash it
      //faceHash[sum].push_back(&(_tets[x].faces[y]));
      faceHash[sum].push_back(pair<int,int>(x,y));
      tetHash[sum].push_back(x);
      triangleHash[sum].push_back(y);
    }

  // go through all the triangles
  // if more than one hashed in, check for duplicates
  //map<long long, vector<TRIANGLE*> >::iterator iter;
  map<long long, vector<pair<int,int> > >::iterator iter;
  vector<int> faceBelongsToTet;
  vector<int> whichFaceInTet;
  int hashSize = faceHash.size();
  int i = 0;
  for (iter = faceHash.begin(); iter != faceHash.end(); iter++, i++)
  {
    //vector<TRIANGLE*> faces = iter->second;
    vector<pair<int,int> > faces = iter->second;
    if (faces.size() == 1)
    {
      _surfaceFaces.push_back(faces[0]);
      faceBelongsToTet.push_back(tetHash[iter->first][0]);
      whichFaceInTet.push_back(triangleHash[iter->first][0]);
      continue;
    }

    // see if this face matches any other ones
    for (unsigned int x = 0; x < faces.size(); x++)
    {
      bool match = false;
      for (unsigned int y = 0; y < faces.size(); y++)
      {
        if (y == x) continue;

        // exploit overloaded TRIANGLE operator
        TRIANGLE left = _tets[faces[x].first].face(faces[x].second);
        TRIANGLE right = _tets[faces[y].first].face(faces[y].second);
        if (left == right) {
          match = true;
          continue;
        }
      }

      // if there are no matches, it is a surface triangle
      if (!match) 
      {
        _surfaceFaces.push_back(faces[x]);
        faceBelongsToTet.push_back(tetHash[iter->first][x]);
        whichFaceInTet.push_back(triangleHash[iter->first][x]);
      }
    }
    if (i % (int)(hashSize / 10) == 0)
    {
        cout << 100 * ((Real)i / hashSize) << "% ";
        flush(cout);
    }
  }

  // write the surface faces to a file
  string filename = string(_filename);
  filename += string(".surfacefaces");
  FILE* file = fopen(filename.c_str(), "wb");
  int totalFaces = faceBelongsToTet.size();
  fwrite((void*)&totalFaces, sizeof(int), 1, file);
  for (int x = 0; x < totalFaces; x++)
  {
    int whichTet = faceBelongsToTet[x];
    int whichFace = whichFaceInTet[x];
    fwrite((void*)&whichTet, sizeof(int), 1, file);
    fwrite((void*)&whichFace, sizeof(int), 1, file);
  }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// compute all tet faces
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeAllFaces()
{
  // compute all surface faces
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    for (int y = 0; y < 4; y++)
    {
      TRIANGLE triangle = _tets[x].face(y);
      TRIANGLE* newTri = new TRIANGLE(&triangle);
      _allFaces.push_back(newTri); 
    }
  }
}

//////////////////////////////////////////////////////////////////////
// compute explicit surface faces
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeExplicitSurfaceFaces()
{
  // compute explicit surface faces
  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    int tetID = _surfaceFaces[x].first;
    int faceID = _surfaceFaces[x].second;
    TRIANGLE triangle = _tets[tetID].face(faceID);

    TRIANGLE* newTri = new TRIANGLE(&triangle);
    _explicitSurfaceFaces.push_back(newTri); 
  }
}

//////////////////////////////////////////////////////////////////////
// read in the surface faces
//////////////////////////////////////////////////////////////////////
bool TET_MESH::readSurfaceFaceCache()
{
  if (_filename.size() == 0) return false;
  
  // write the surface faces to a file
  string filename = string(_filename);
  filename += string(".surfacefaces");
  FILE* file = fopen(filename.c_str(), "rb");
  if (file == NULL) return false;
  
  int totalFaces;
  fread((void*)&totalFaces, sizeof(int), 1, file);
  for (int x = 0; x < totalFaces; x++)
  {
    int whichTet;
    int whichFace;
    fread((void*)&whichTet, sizeof(int), 1, file);
    fread((void*)&whichFace, sizeof(int), 1, file);

    //TRIANGLE* triangle = &(_tets[whichTet].faces[whichFace]);
    //_surfaceFaces.push_back(triangle);
    _surfaceFaces.push_back(pair<int,int>(whichTet, whichFace));
  }
  fclose(file);

  /*
  // compute explicit surface faces
  for (int x = 0; x < _surfaceFaces.size(); x++)
  {
    int tetID = _surfaceFaces[x].first;
    int faceID = _surfaceFaces[x].second;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    TRIANGLE triangle = _tets[tetID].face(faceID);
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;

    TRIANGLE* newTri = new TRIANGLE(&triangle);
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    _explicitSurfaceFaces.push_back(newTri); 
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  }
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  */
  return true;
}

//////////////////////////////////////////////////////////////////////
// Construct a vector of the vertices on the mesh surface
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceVertices()
{
  if (_filename.size() == 0) return;
  
  _surfaceVertices.clear();

  // try to read in a cached version
  string filename = string(_filename);
  filename += string(".surfacevertices");
  FILE* file = fopen(filename.c_str(), "rb");
  if (file != NULL)
  {
    int totalVertices;
    fread((void*)&totalVertices, sizeof(int), 1, file);
    for (int x = 0; x < totalVertices; x++)
    {
      int whichVertex;
      fread((void*)&whichVertex, sizeof(int), 1, file);
      _surfaceVertices.push_back(&(_vertices[whichVertex]));
    }
    fclose(file);
    //cout << " vertex cache found ... ";
    return;
  }

  // keep track of what vertices we have already added to the vector
  map<VEC3F*, bool> vertexHash;
  map<VEC3F*, bool>::iterator iter;

  // loop through all the surface faces
  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      // check that the vertex was not already added to the list
      TRIANGLE triangle = _tets[_surfaceFaces[x].first].face(_surfaceFaces[x].second);
      //VEC3F* vertex = _surfaceFaces[x]->vertex(y);
      VEC3F* vertex = triangle.vertex(y);
      if (vertexHash.find(vertex) == vertexHash.end())
      {
        // add it to the list and hash
        _surfaceVertices.push_back(vertex);
        vertexHash[vertex] = true;
      }
    }

  // write the surface faces to a file
  filename = string(_filename);
  filename += string(".surfacevertices");
  file = fopen(filename.c_str(), "wb");
  int totalVertices = _surfaceVertices.size();
  fwrite((void*)&totalVertices, sizeof(int), 1, file);
  for (int x = 0; x < totalVertices; x++)
  {
    int whichVertex = _vertexID[_surfaceVertices[x]];
    fwrite((void*)&whichVertex, sizeof(int), 1, file);
  }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// set the masses of all the vertices
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeMasses()
{
  // just set them all to sum to unity
  int size = 3 * _unconstrainedSize;
  //Real mass = 1.0 / _unconstrainedSize;
  
  // DEBUG: scaling the mass
  //Real mass = 10000.0 / _unconstrainedSize;
  Real mass = 1.0 / _unconstrainedSize;

  // populate a sparse matrix
  _masses.resize(size, size);
  for (int x = 0; x < size; x++)
    _masses(x,x) = mass;
}

//////////////////////////////////////////////////////////////////////
// rescale the masses of all the vertices
//////////////////////////////////////////////////////////////////////
void TET_MESH::scaleMasses(Real scale)
{
  // just set them all to sum to unity
  int size = 3 * _unconstrainedSize;

  // populate a sparse matrix
  for (int x = 0; x < size; x++)
    _masses(x,x) *= scale;
}

//////////////////////////////////////////////////////////////////////
// Peek the eigenvalues of each tet and return the smallest found
//////////////////////////////////////////////////////////////////////
Real TET_MESH::inspectStiffnessEigenvalues()
{
  Real minFound = 100000;
  int totalNegative = 0;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    int materialIndex = _tets[x].materialIndex();

    if (_tets[x].invalid())
    {
      cout << "Invalid tet found !!! " << endl;
      exit(0);
    }
    
    MATRIX K = _materials[materialIndex]->stiffness(_tets[x]);
    K *= -1.0;

    VECTOR eigenvalues(K.rows());
    MATRIX eigenvectors(K.rows(), K.cols());

    MATRIX copyK(K);
    copyK.eigensystem(eigenvalues, eigenvectors);

    if (x == 0) minFound = eigenvalues[0];

    if (x == 0)
    {
      cout << " Reference stiffness: " << endl;
      cout << K << endl;

      cout << " eigenvalues: " << eigenvalues << endl;
    }

    for (int y = 0; y < eigenvalues.size(); y++)
      if (eigenvalues[y] < 0.0)
      {
        totalNegative++;
        break;
      }

    for (int y = 0; y < eigenvalues.size(); y++)
    {
      if (eigenvalues[y] < minFound)
      {
        cout << " new min found: " << eigenvalues[y] << endl;
        cout << "K: " << K << endl;
      }

      minFound = (eigenvalues[y] < minFound) ? eigenvalues[y] : minFound;
    }
    
    if (x % (int)(_tets.size() / 10) == 0)
    {
      cout << 100 * ((Real)x / _tets.size()) << "% ";
      flush(cout);
    }
  }

  cout << endl;
  cout << " === " << (Real)totalNegative / _tets.size() * 100.0 << "% negative ===" << endl;
  return minFound;
}

#ifndef IGNORE_PETSC
//////////////////////////////////////////////////////////////////////
// compute the stiffness matrix about the current pose
//////////////////////////////////////////////////////////////////////
void TET_MESH::generateSparseStiffnessMatrix(SPARSE_PETSC_MATRIX& stiffness)
{
  if (stiffness.rows() == 0)
  {
    cout << __FILE__ << " " << __LINE__ << " : " << endl;
    cout << " SPARSE_PETSC_MATRIX sparsity has not been set! Initialize with a SPARSE_MATRIX example!" << endl;
  }
  _stiffness.clear();
  if (_stiffness.rows() != 3 * _unconstrainedSize)
    _stiffness.resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    int materialIndex = _tets[x].materialIndex();
    
    MATRIX K = _materials[materialIndex]->stiffness(_tets[x]);
    int indices[] = {_vertexID[_tets[x].vertices[0]],
                     _vertexID[_tets[x].vertices[1]],
                     _vertexID[_tets[x].vertices[2]],
                     _vertexID[_tets[x].vertices[3]]};
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
          {
            if (indices[i] >= _unconstrainedSize ||
                indices[j] >= _unconstrainedSize) continue;
            int row = indices[i] * 3 + k;
            int col = indices[j] * 3 + l;

            stiffness.add(row, col, -K(i*3 + k, j*3 + l));

            // TESTING
            _stiffness(row, col) -= K(i*3 + k, j*3 + l);
          }
  }

  /*
#ifdef USING_OPENMP
#pragma omp parallel
#endif
  { 
#ifdef USING_OPENMP
    int id  = omp_get_thread_num();
#else
    const int id  = 0;
#endif
#pragma omp for  schedule(static)
    for (int x = 0; x < _tets.size(); x++)
    {
      int materialIndex = _tets[x].materialIndex();
      
      MATRIX K = _materialCopies[id][materialIndex]->stiffness(_tets[x]);
      int indices[] = {_vertexID[_tets[x].vertices[0]],
                       _vertexID[_tets[x].vertices[1]],
                       _vertexID[_tets[x].vertices[2]],
                       _vertexID[_tets[x].vertices[3]]};
      #pragma omp critical
      {
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
            {
              if (indices[i] >= _unconstrainedSize ||
                  indices[j] >= _unconstrainedSize) continue;
              int row = indices[i] * 3 + k;
              int col = indices[j] * 3 + l;

              stiffness.add(row, col, -K(i*3 + k, j*3 + l));
            }
      }
    }
  }
  */
}
#endif

//////////////////////////////////////////////////////////////////////
// compute the stiffness matrix about the current pose
//////////////////////////////////////////////////////////////////////
void TET_MESH::generateSparseStiffnessMatrix(SPARSE_MATRIX& stiffness)
{
  stiffness.clear();
  if (stiffness.rows() != 3 * _unconstrainedSize)
    stiffness.resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    int materialIndex = _tets[x].materialIndex();
    
    MATRIX K = _materials[materialIndex]->stiffness(_tets[x]);
    int indices[] = {_vertexID[_tets[x].vertices[0]],
                     _vertexID[_tets[x].vertices[1]],
                     _vertexID[_tets[x].vertices[2]],
                     _vertexID[_tets[x].vertices[3]]};
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
          {
            if (indices[i] >= _unconstrainedSize ||
                indices[j] >= _unconstrainedSize) continue;
            int row = indices[i] * 3 + k;
            int col = indices[j] * 3 + l;

            stiffness(row, col) -= K(i*3 + k, j*3 + l);
          }

    /*
    if (x % (int)(_tets.size() / 10) == 0)
    {
      cout << 100 * ((Real)x / _tets.size()) << "% ";
      flush(cout);
    }
    */
  }
  //cout << " Size of stiffness matrix: " << stiffness.size() << endl;
}

//////////////////////////////////////////////////////////////////////
// compute the stiffness matrix about the current pose
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& TET_MESH::generateSparseStiffnessMatrix(bool verbose)
{
  _stiffness.clear();
  if (_stiffness.rows() != 3 * _unconstrainedSize)
    _stiffness.resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    int materialIndex = _tets[x].materialIndex();
    
    MATRIX K = _materials[materialIndex]->stiffness(_tets[x]);
    int indices[] = {_vertexID[_tets[x].vertices[0]],
                     _vertexID[_tets[x].vertices[1]],
                     _vertexID[_tets[x].vertices[2]],
                     _vertexID[_tets[x].vertices[3]]};
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
          {
            if (indices[i] >= _unconstrainedSize ||
                indices[j] >= _unconstrainedSize) continue;
            int row = indices[i] * 3 + k;
            int col = indices[j] * 3 + l;

            /*
            if (isnan(K(i*3 + k, j*3 + l)))
            {
              cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
              cout << " NAN FOUND!!! " << endl;

              cout << "Stiffness matrix: " << K << endl;
              cout << " Tet: " << _tets[x] << endl;

              MATRIX3 F = _tets[x].F();
              cout << " F: " << F << endl;
              ((INVERTIBLE*)_materials[0])->debugF(F);
              exit(1);
            }
            */

            _stiffness(row, col) -= K(i*3 + k, j*3 + l);
          }

    if (verbose)
      if (x % (int)(_tets.size() / 10) == 0)
      {
        cout << 100 * ((Real)x / _tets.size()) << "% ";
        flush(cout);
      }
  }
  /*
#ifndef _WIN32    
  cout << " Size of stiffness matrix: " << _stiffness.size() << endl;
#endif
*/

  return _stiffness;

  /*
  // do the resize in serial in case a core is busy
  for (int x = 0; x < _totalCores; x++)
  {
    if (_stiffnessCopies[x].rows() != 3 * _unconstrainedSize)
      _stiffnessCopies[x].resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);
    _stiffnessCopies[x].clear();
  }
  cout << __FILE__ << " " << __LINE__ << " : " << endl; flush(cout);

#ifdef USING_OPENMP
#pragma omp parallel
#endif
  { 
#ifdef USING_OPENMP
    int id  = omp_get_thread_num();
#else
    const int id  = 0;
#endif
#pragma omp for  schedule(static)
    for (unsigned int x = 0; x < _tets.size(); x++)
    {
      int materialIndex = _tets[x].materialIndex();
      
      //MATRIX K = _materialCopies[id]->stiffness(_tets[x]);
      MATRIX K = _materialCopies[id][materialIndex]->stiffness(_tets[x]);
      int indices[] = {_vertexID[_tets[x].vertices[0]],
                       _vertexID[_tets[x].vertices[1]],
                       _vertexID[_tets[x].vertices[2]],
                       _vertexID[_tets[x].vertices[3]]};
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
            {
              if (indices[i] >= _unconstrainedSize ||
                  indices[j] >= _unconstrainedSize) continue;
              int row = indices[i] * 3 + k;
              int col = indices[j] * 3 + l;

              _stiffnessCopies[id](row, col) -= K(i*3 + k, j*3 + l);
            }
    }
  }
  cout << __FILE__ << " " << __LINE__ << " : " << endl; flush(cout);

  // merge results
  for (int x = 0; x < _totalCores; x++)
    _stiffness += _stiffnessCopies[x];
  cout << __FILE__ << " " << __LINE__ << " : " << endl; flush(cout);

  return _stiffness;
  */
}

//////////////////////////////////////////////////////////////////////
// compute the stiffness matrix about the rest pose using invertible
// clamping
//////////////////////////////////////////////////////////////////////
#ifndef IGNORE_PETSC
#if _WIN32
void TET_MESH::generateSparseStiffnessMatrix(vector<MATRIX3>& Us, 
                                             vector<MATRIX3>& Vs, 
                                             vector<MATRIX>& stiffnesses,
                                             SPARSE_PCG_MATRIX& stiffness)
#else
void TET_MESH::generateSparseStiffnessMatrix(vector<MATRIX3>& Us, 
                                             vector<MATRIX3>& Vs, 
                                             vector<MATRIX>& stiffnesses,
                                             SPARSE_PETSC_MATRIX& stiffness)
#endif
{
  stiffness.clear();
  if (stiffness.rows() == 0)
  {
    cout << __FILE__ << " " << __LINE__ << " : " << endl;
    cout << " SPARSE_PETSC_MATRIX sparsity has not been set! Initialize with a SPARSE_MATRIX example!" << endl;
  }
  _stiffness.clear();
  if (_stiffness.rows() != 3 * _unconstrainedSize)
    _stiffness.resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);

  // do the resize in serial in case a core is busy
  for (int x = 0; x < _totalCores; x++)
  {
    if (_stiffnessCopies[x].rows() != 3 * _unconstrainedSize)
      _stiffnessCopies[x].resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);
    _stiffnessCopies[x].clear();
  }

#ifdef USING_OPENMP
#pragma omp parallel
#endif
  { 
#ifdef USING_OPENMP
    int id  = omp_get_thread_num();
#pragma omp for  schedule(static)
#else
    const int id  = 0;
#endif
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    TET& tet = _tets[x];
    int materialIndex = _tets[x].materialIndex();

    // get cached stiffness matrix
    MATRIX& diagonalStiffness = stiffnesses[x];
    MATRIX3& U = Us[x];
    MATRIX3& V = Vs[x];

    // compute Bm matrix
    const VEC3F* b = tet.b();
    
    // get PFPu from the tet - only depends on rest state,
    // so don't need to update tet state at all
    MATRIX pFpu(9,12);
    _materialCopies[id][materialIndex]->computePFPu(tet, pFpu);
    MATRIX stiffnessSubblock(12,12);

    // compute each column of the matrix
    for (int y = 0; y < 12; y++)
    {
      // extract a column from pFpu (ie pretend 'x' is just a 
      // Kronecker delta)
      VECTOR deltaF(9);
      for (int z = 0; z < 9; z++)
        deltaF(z) = pFpu(z, y);

      // rotate deltaF
      MATRIX3 rotated = U.transpose() * TET::repackF(deltaF) * V;
      deltaF = TET::flattenF(rotated);

      VECTOR contraction = diagonalStiffness  * deltaF;
      MATRIX3 deltaP = TET::repackF(contraction);

      // rotate deltaP back
      deltaP = U * deltaP * V.transpose();
      VEC3F forceVecs[4];
      forceVecs[0] = deltaP * b[0];
      forceVecs[1] = deltaP * b[1];
      forceVecs[2] = deltaP * b[2];
      forceVecs[3] = deltaP * b[3];
     
      // copy result into stiffness column
      for (int z = 0; z < 4; z++)
        for (int a = 0; a < 3; a++)
          stiffnessSubblock(z * 3 + a, y) = forceVecs[z][a];
    }
    stiffnessSubblock *= -1.0;
   
    int indices[] = {_vertexID[_tets[x].vertices[0]],
                     _vertexID[_tets[x].vertices[1]],
                     _vertexID[_tets[x].vertices[2]],
                     _vertexID[_tets[x].vertices[3]]};

    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
          {
            if (indices[i] >= _unconstrainedSize ||
                indices[j] >= _unconstrainedSize) continue;
            int row = indices[i] * 3 + k;
            int col = indices[j] * 3 + l;
            //stiffness.add(row, col, -stiffnessSubblock(i*3 + k, j*3 + l));
            _stiffnessCopies[id](row, col) -= stiffnessSubblock(i*3 + k, j*3 + l);
          }
  }
  } // OMP

  // merge results
  for (int x = 0; x < _totalCores; x++)
    _stiffness += _stiffnessCopies[x];

  // copy everything into the PetSc matrix as well
  const map<pair<int,int>, Real>& matrix = _stiffness.matrix();
  map<pair<int,int>, Real>::const_iterator iter;
  for (iter = matrix.begin(); iter != matrix.end(); iter++)
  {
    const int row = iter->first.first;
    const int col = iter->first.second;
    const Real entry = iter->second;
    stiffness.add(row, col, entry);
  }
}
#endif

// SERIAL VERSION
/*
{
  stiffness.clear();
  if (stiffness.rows() == 0)
  {
    cout << __FILE__ << " " << __LINE__ << " : " << endl;
    cout << " SPARSE_PETSC_MATRIX sparsity has not been set! Initialize with a SPARSE_MATRIX example!" << endl;
  }
  _stiffness.clear();
  if (_stiffness.rows() != 3 * _unconstrainedSize)
    _stiffness.resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);

  for (int x = 0; x < _tets.size(); x++)
  {
    TET& tet = _tets[x];
    int materialIndex = _tets[x].materialIndex();

    // get cached stiffness matrix
    MATRIX& diagonalStiffness = stiffnesses[x];
    MATRIX3& U = Us[x];
    MATRIX3& V = Vs[x];

    // compute Bm matrix
    const VEC3F* b = tet.b();
    
    // get PFPu from the tet - only depends on rest state,
    // so don't need to update tet state at all
    MATRIX pFpu(9,12);
    //_material->computePFPu(tet, pFpu);
    _materials[materialIndex]->computePFPu(tet, pFpu);
    MATRIX stiffnessSubblock(12,12);

    // compute each column of the matrix
    for (int y = 0; y < 12; y++)
    {
      // extract a column from pFpu (ie pretend 'x' is just a 
      // Kronecker delta)
      VECTOR deltaF(9);
      for (int z = 0; z < 9; z++)
        deltaF(z) = pFpu(z, y);

      // rotate deltaF
      MATRIX3 rotated = U.transpose() * TET::repackF(deltaF) * V;
      deltaF = TET::flattenF(rotated);

      VECTOR contraction = diagonalStiffness  * deltaF;
      MATRIX3 deltaP = TET::repackF(contraction);

      // rotate deltaP back
      deltaP = U * deltaP * V.transpose();
      VEC3F forceVecs[4];
      forceVecs[0] = deltaP * b[0];
      forceVecs[1] = deltaP * b[1];
      forceVecs[2] = deltaP * b[2];
      forceVecs[3] = deltaP * b[3];
     
      // copy result into stiffness column
      for (int z = 0; z < 4; z++)
        for (int a = 0; a < 3; a++)
          stiffnessSubblock(z * 3 + a, y) = forceVecs[z][a];
    }
    stiffnessSubblock *= -1.0;
   
    int indices[] = {_vertexID[_tets[x].vertices[0]],
                     _vertexID[_tets[x].vertices[1]],
                     _vertexID[_tets[x].vertices[2]],
                     _vertexID[_tets[x].vertices[3]]};

    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
          {
            if (indices[i] >= _unconstrainedSize ||
                indices[j] >= _unconstrainedSize) continue;
            int row = indices[i] * 3 + k;
            int col = indices[j] * 3 + l;
            _stiffness(row, col) -= stiffnessSubblock(i*3 + k, j*3 + l);
            stiffness.add(row, col, -stiffnessSubblock(i*3 + k, j*3 + l));
          }
  }
}
*/


//////////////////////////////////////////////////////////////////////
// compute the stiffness matrix about the rest pose using invertible
// clamping
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& TET_MESH::generateOneRingStiffness(vector<MATRIX3>& Us, 
                                                  vector<MATRIX3>& Vs, 
                                                  vector<MATRIX>& stiffnesses)
{
  int size = _vertices.size();

  _stiffness.clear();
  if (_stiffness.rows() != 3 * size)
    _stiffness.resize(3 * size, 3 * size);

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    TET& tet = _tets[x];
    int materialIndex = _tets[x].materialIndex();

    // get cached stiffness matrix
    MATRIX& diagonalStiffness = stiffnesses[x];
    MATRIX3& U = Us[x];
    MATRIX3& V = Vs[x];

    // compute Bm matrix
    const VEC3F* b = tet.b();
    
    // get PFPu from the tet - only depends on rest state,
    // so don't need to update tet state at all
    MATRIX pFpu(9,12);
    _materials[materialIndex]->computePFPu(tet, pFpu);
    MATRIX stiffnessSubblock(12,12);

    // compute each column of the matrix
    for (int y = 0; y < 12; y++)
    {
      // extract a column from pFpu (ie pretend 'x' is just a 
      // Kronecker delta)
      VECTOR deltaF(9);
      for (int z = 0; z < 9; z++)
        deltaF(z) = pFpu(z, y);

      // rotate deltaF
      MATRIX3 rotated = U.transpose() * TET::repackF(deltaF) * V;
      deltaF = TET::flattenF(rotated);

      VECTOR contraction = diagonalStiffness  * deltaF;
      MATRIX3 deltaP = TET::repackF(contraction);

      // rotate deltaP back
      deltaP = U * deltaP * V.transpose();
      VEC3F forceVecs[4];
      forceVecs[0] = deltaP * b[0];
      forceVecs[1] = deltaP * b[1];
      forceVecs[2] = deltaP * b[2];
      forceVecs[3] = deltaP * b[3];
     
      // copy result into stiffness column
      for (int z = 0; z < 4; z++)
        for (int a = 0; a < 3; a++)
          stiffnessSubblock(z * 3 + a, y) = forceVecs[z][a];
    }
    stiffnessSubblock *= -1.0;
   
    int indices[] = {_vertexID[_tets[x].vertices[0]],
                     _vertexID[_tets[x].vertices[1]],
                     _vertexID[_tets[x].vertices[2]],
                     _vertexID[_tets[x].vertices[3]]};

    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
          {
            int row = indices[i] * 3 + k;
            int col = indices[j] * 3 + l;
            _stiffness(row, col) -= stiffnessSubblock(i*3 + k, j*3 + l);
          }
  }
  return _stiffness;
}


//////////////////////////////////////////////////////////////////////
// compute the stiffness matrix about the rest pose using invertible
// clamping
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& TET_MESH::generateSparseStiffnessMatrix(vector<MATRIX3>& Us, 
                                                       vector<MATRIX3>& Vs, 
                                                       vector<MATRIX>& stiffnesses)
#if USING_OPENMP
{
  _stiffness.clear();
  if (_stiffness.rows() != 3 * _unconstrainedSize)
    _stiffness.resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);

  // do the resize in serial in case a core is busy
  for (int x = 0; x < _totalCores; x++)
  {
    if (_stiffnessCopies[x].rows() != 3 * _unconstrainedSize)
      _stiffnessCopies[x].resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);
    _stiffnessCopies[x].clear();
  }

#pragma omp parallel
  { 
    int id  = omp_get_thread_num();
    //const int id  = 0;
#pragma omp for  schedule(static)
  for (int x = 0; x < _tets.size(); x++)
  {
    TET& tet = _tets[x];
    int materialIndex = _tets[x].materialIndex();

    // get cached stiffness matrix
    MATRIX& diagonalStiffness = stiffnesses[x];
    MATRIX3& U = Us[x];
    MATRIX3& V = Vs[x];

    // compute Bm matrix
    const VEC3F* b = tet.b();
    
    // get PFPu from the tet - only depends on rest state,
    // so don't need to update tet state at all
    MATRIX pFpu(9,12);
    _materialCopies[id][materialIndex]->computePFPu(tet, pFpu);
    MATRIX stiffnessSubblock(12,12);

    // compute each column of the matrix
    for (int y = 0; y < 12; y++)
    {
      // extract a column from pFpu (ie pretend 'x' is just a 
      // Kronecker delta)
      VECTOR deltaF(9);
      for (int z = 0; z < 9; z++)
        deltaF(z) = pFpu(z, y);

      // rotate deltaF
      MATRIX3 rotated = U.transpose() * TET::repackF(deltaF) * V;
      deltaF = TET::flattenF(rotated);

      VECTOR contraction = diagonalStiffness  * deltaF;
      MATRIX3 deltaP = TET::repackF(contraction);

      // rotate deltaP back
      deltaP = U * deltaP * V.transpose();
      VEC3F forceVecs[4];
      forceVecs[0] = deltaP * b[0];
      forceVecs[1] = deltaP * b[1];
      forceVecs[2] = deltaP * b[2];
      forceVecs[3] = deltaP * b[3];
     
      // copy result into stiffness column
      for (int z = 0; z < 4; z++)
        for (int a = 0; a < 3; a++)
          stiffnessSubblock(z * 3 + a, y) = forceVecs[z][a];
    }
    stiffnessSubblock *= -1.0;
   
    int indices[] = {_vertexID[_tets[x].vertices[0]],
                     _vertexID[_tets[x].vertices[1]],
                     _vertexID[_tets[x].vertices[2]],
                     _vertexID[_tets[x].vertices[3]]};

    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
          {
            if (indices[i] >= _unconstrainedSize ||
                indices[j] >= _unconstrainedSize) continue;
            int row = indices[i] * 3 + k;
            int col = indices[j] * 3 + l;
            //_stiffness(row, col) -= K(i*3 + k, j*3 + l);
            //_stiffness(row, col) -= stiffnessSubblock(i*3 + k, j*3 + l);
            _stiffnessCopies[id](row, col) -= stiffnessSubblock(i*3 + k, j*3 + l);
          }
  }
  } // OMP

  // merge results
  for (int x = 0; x < _totalCores; x++)
    _stiffness += _stiffnessCopies[x];

  return _stiffness;
}
#else
{
  _stiffness.clear();
  if (_stiffness.rows() != 3 * _unconstrainedSize)
    _stiffness.resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    TET& tet = _tets[x];
    int materialIndex = _tets[x].materialIndex();

    // get cached stiffness matrix
    MATRIX& diagonalStiffness = stiffnesses[x];
    MATRIX3& U = Us[x];
    MATRIX3& V = Vs[x];

    // compute Bm matrix
    const VEC3F* b = tet.b();
    
    // get PFPu from the tet - only depends on rest state,
    // so don't need to update tet state at all
    MATRIX pFpu(9,12);
    //_material->computePFPu(tet, pFpu);
    _materials[materialIndex]->computePFPu(tet, pFpu);
    MATRIX stiffnessSubblock(12,12);

    // compute each column of the matrix
    for (int y = 0; y < 12; y++)
    {
      // extract a column from pFpu (ie pretend 'x' is just a 
      // Kronecker delta)
      VECTOR deltaF(9);
      for (int z = 0; z < 9; z++)
        deltaF(z) = pFpu(z, y);

      // rotate deltaF
      MATRIX3 rotated = U.transpose() * TET::repackF(deltaF) * V;
      deltaF = TET::flattenF(rotated);

      VECTOR contraction = diagonalStiffness  * deltaF;
      MATRIX3 deltaP = TET::repackF(contraction);

      // rotate deltaP back
      deltaP = U * deltaP * V.transpose();
      VEC3F forceVecs[4];
      forceVecs[0] = deltaP * b[0];
      forceVecs[1] = deltaP * b[1];
      forceVecs[2] = deltaP * b[2];
      forceVecs[3] = deltaP * b[3];
     
      // copy result into stiffness column
      for (int z = 0; z < 4; z++)
        for (int a = 0; a < 3; a++)
          stiffnessSubblock(z * 3 + a, y) = forceVecs[z][a];
    }
    stiffnessSubblock *= -1.0;
   
    int indices[] = {_vertexID[_tets[x].vertices[0]],
                     _vertexID[_tets[x].vertices[1]],
                     _vertexID[_tets[x].vertices[2]],
                     _vertexID[_tets[x].vertices[3]]};

    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
          {
            if (indices[i] >= _unconstrainedSize ||
                indices[j] >= _unconstrainedSize) continue;
            int row = indices[i] * 3 + k;
            int col = indices[j] * 3 + l;
            //_stiffness(row, col) -= K(i*3 + k, j*3 + l);
            _stiffness(row, col) -= stiffnessSubblock(i*3 + k, j*3 + l);
          }
  }

  return _stiffness;
}
#endif

//////////////////////////////////////////////////////////////////////
// is any tet in the mesh inverted?
//////////////////////////////////////////////////////////////////////
bool TET_MESH::inverted(bool verbose)
{
  for (unsigned int x = 0; x < _tets.size(); x++)
    if (_tets[x].inverted())
    {
      /*
      if (verbose)
      {
        cout << " Inverted tet: " << endl;
        MATRIX3 F = _tets[x].F();
        cout << "   F: " << F << endl;
        MATRIX3 firstPK = _materials[0]->firstPiolaKirchhoff(F);
        cout << "   firstPK: " << firstPK << endl;
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        ((INVERTIBLE*)_materials[0])->debugF(F);
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      }
      */
      return true;
    }
  return false;
}

//////////////////////////////////////////////////////////////////////
// clear all internal forces
//////////////////////////////////////////////////////////////////////
void TET_MESH::clearInternalForces()
{
  for (int x = 0; x < _unconstrainedSize; x++)
    _internalForces[x].clear();
}

//////////////////////////////////////////////////////////////////////
// compute internal forces
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeInternalForces()
{
  clearInternalForces();

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    // cache some tet values
    MATRIX3 F = _tets[x].F();
    const VEC3F* b = _tets[x].b();
    int materialIndex = _tets[x].materialIndex();

    // get the second PK for the tet
    MATRIX3 P = _materials[materialIndex]->firstPiolaKirchhoff(F);

    // calculate the final forces
    VEC3F forces[4];
    forces[0] = P * b[0];
    forces[1] = P * b[1];
    forces[2] = P * b[2];
    forces[3] = P * b[3];

    // distribute forces to the nodes
    for (int y = 0; y < 4; y++)
    {
      int index = _vertexID[_tets[x].vertices[y]];

      // make sure it's an unconstrained node
      if (index < _unconstrainedSize)
        _internalForces[index] += forces[y];
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Reset all vertices to the rest pose
//////////////////////////////////////////////////////////////////////
void TET_MESH::resetToRestPose()
{
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] = _restPose[x];
}

//////////////////////////////////////////////////////////////////////
// Given a position, return the closest surface node
//////////////////////////////////////////////////////////////////////
VEC3F* TET_MESH::closestSurfaceNode(VEC3F point)
{
  // make sure a surface list was built 
  if (_surfaceVertices.size() == 0) return NULL;

  // tentatively set it to the first in the list
  VEC3F* closest;
  closest = _surfaceVertices[0];
  Real minDistance = norm2(point - *_surfaceVertices[0]);

  // loop through the rest of the vertices
  for (unsigned int x = 1; x < _surfaceVertices.size(); x++)
  {
    // check if this one is closer
    Real distance = norm2(point - *_surfaceVertices[x]);
    if (distance < minDistance)
    {
      minDistance = distance;
      closest = _surfaceVertices[x];
    }
  }
  return closest;
}

//////////////////////////////////////////////////////////////////////
// Find all surface nodes within a given radius
//////////////////////////////////////////////////////////////////////
void TET_MESH::closestSurfaceNodes( VEC3F point, Real radius,
                                    vector<VEC3F *> &surfaceNodes )
{
  // make sure a surface list was built 
  if (_surfaceVertices.size() == 0) return;

  surfaceNodes.clear();

  Real minDistance = radius * radius;

  // loop through the rest of the vertices
  for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
  {
    // check if this one is closer
    Real distance = norm2(point - *_surfaceVertices[x]);
    if (distance < minDistance)
    {
      surfaceNodes.push_back( _surfaceVertices[x] );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Write the tet mesh to a METIS friendly format for partitioning
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeMETIS(const char* filename)
{
  FILE* file = fopen(filename, "w");

  // total number of tets and the flag telling METIS that they are
  // tets instead of triangles or hexahedra
  fprintf(file, "%i 2\n", (int)(_tets.size()));

  // METIS assumes 1-indexed
  for (unsigned int x = 0; x < _tets.size(); x++)
    fprintf(file, "%i %i %i %i\n",_vertexID[_tets[x].vertices[0]] + 1,
                                  _vertexID[_tets[x].vertices[1]] + 1,
                                  _vertexID[_tets[x].vertices[2]] + 1,
                                  _vertexID[_tets[x].vertices[3]] + 1);
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read in the results of a METIS partitioning and tag the 'partition' 
// field of the tets
//////////////////////////////////////////////////////////////////////
void TET_MESH::readMETIS(const char* filename)
{
  FILE* file = fopen(filename, "r");
  if (file == NULL)
  {
    cout << " No partition file " << filename << " found!" << endl;
    return;
  }

  // keep count of the tets per partition
  map<int,int> totals;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    int partition = 0;
    fscanf(file, "%i\n", &partition);
    _tets[x].partition = partition; 

    totals[partition]++;
  }
  fclose(file);

  // print out some stats
  cout << "==================================================" << endl;
  cout << " Partitioning stats: " << endl;
  cout << "==================================================" << endl;
  map<int,int>::iterator iter;
  for (iter = totals.begin(); iter != totals.end(); iter++)
    cout << " Partition " << iter->first << ": " 
                          << iter->second << " tets " << endl;
  cout << "==================================================" << endl;
}

//////////////////////////////////////////////////////////////////////
// Write the tet mesh to a Scotch friendly format for partitioning
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeScotch(const char* filename)
{
  // compute the total number of arcs (i.e. 2x the number of edges)
  cout << " Building arcs ..."; flush(cout);
  int totalArcs = 0;
  vector<vector<int> > arcs;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    vector<int> oneRing;
    tetOneRing(x, oneRing);
    totalArcs += oneRing.size();

    arcs.push_back(oneRing);
  }
  cout << " done." << endl;

  FILE* file = fopen(filename, "w");
  fprintf(file, "0\n");
  fprintf(file, "%i\t%i\n", (int)(_tets.size()), totalArcs);
  fprintf(file, "0\t000\n");

  cout << " Writing arcs ..."; flush(cout);
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    vector<int>& oneRing = arcs[x];
    fprintf(file, "%i", (int)(oneRing.size()));
    for (unsigned int y = 0; y < oneRing.size(); y++)
      fprintf(file, "\t%i", oneRing[y]);
    fprintf(file, "\n");
  }
  cout << " done." << endl;

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read in the results of a METIS partitioning and tag the 'partition' 
// field of the tets
//////////////////////////////////////////////////////////////////////
void TET_MESH::readScotch(const char* filename)
{
  FILE* file = fopen(filename, "r");
  if (file == NULL)
  {
    cout << " No partition file " << filename << " found!" << endl;
    return;
  }

  int totalLines;
  fscanf(file, "%i\n", &totalLines);
  if (totalLines != (int)_tets.size())
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Number of output tags does not match the number of tets in the mesh! " << endl;
    exit(0);
  }

  // keep count of the tets per partition
  map<int,int> totals;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    int index;
    int partition = 0;
    fscanf(file, "%i\t%i\n", &index, &partition);
    _tets[index].partition = partition; 

    totals[partition]++;
  }
  fclose(file);

  // print out some stats
  cout << "==================================================" << endl;
  cout << " Partitioning stats: " << endl;
  cout << "==================================================" << endl;
  map<int,int>::iterator iter;
  for (iter = totals.begin(); iter != totals.end(); iter++)
    cout << " Partition " << iter->first << ": " 
                          << iter->second << " tets " << endl;
  cout << "==================================================" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generate all the deformation gradients
//////////////////////////////////////////////////////////////////////
void TET_MESH::generateF()
{
  assert(_F.size() == (int)(9 * _tets.size()));
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    MATRIX3 F = _tets[x].F();
    _F(9 * x)     = F(0,0);
    _F(9 * x + 1) = F(1,0);
    _F(9 * x + 2) = F(2,0);
    _F(9 * x + 3) = F(0,1);
    _F(9 * x + 4) = F(1,1);
    _F(9 * x + 5) = F(2,1);
    _F(9 * x + 6) = F(0,2);
    _F(9 * x + 7) = F(1,2);
    _F(9 * x + 8) = F(2,2);
  }
}

//////////////////////////////////////////////////////////////////////
// Evaluate internal forces for Newmark
//////////////////////////////////////////////////////////////////////
VECTOR& TET_MESH::generateInternalForces() 
{
  _R.resizeAndWipe(_unconstrainedSize * 3);

  // wipe previous internal forces
  _R.clear();

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    // compute the forces
    VEC3F forces[4];
    const VEC3F* b = _tets[x].b();
    int materialIndex = _tets[x].materialIndex();

    MATRIX3 F = _tets[x].F();
    MATRIX3 firstPK = _materials[materialIndex]->firstPiolaKirchhoff(F);
    
    forces[0] = firstPK * b[0];
    forces[1] = firstPK * b[1];
    forces[2] = firstPK * b[2];
    forces[3] = firstPK * b[3];

    // add forces to corresponding vertices
    for (int y = 0; y < 4; y++)
    {
      int index = _vertexID[_tets[x].vertices[y]];
      if (index < _unconstrainedSize)
      {
        _R(index * 3)     += forces[y][0];
        _R(index * 3 + 1) += forces[y][1];
        _R(index * 3 + 2) += forces[y][2];
      }
    }
  }

  return _R;

  /*
  _R.resizeAndWipe(_unconstrainedSize * 3);

  // wipe previous internal forces
  _R.clear();

  for (int x = 0; x < _totalCores; x++)
  {
    _RCopies[x].resizeAndWipe(_unconstrainedSize * 3);
    _RCopies[x].clear();
  }

#ifdef USING_OPENMP
#pragma omp parallel
#endif
  { 
#ifdef USING_OPENMP
    const int id  = omp_get_thread_num();
#pragma omp for  schedule(static)
#else
    const int id  = 0;
#endif
    for (int x = 0; x < _tets.size(); x++)
    {
      // compute the forces
      VEC3F forces[4];
      const VEC3F* b = _tets[x].b();
      int materialIndex = _tets[x].materialIndex();

      MATRIX3 F = _tets[x].F();
      MATRIX3 firstPK = _materialCopies[id][materialIndex]->firstPiolaKirchhoff(F);
      
      forces[0] = firstPK * b[0];
      forces[1] = firstPK * b[1];
      forces[2] = firstPK * b[2];
      forces[3] = firstPK * b[3];

      // add forces to corresponding vertices
      for (int y = 0; y < 4; y++)
      {
        int index = _vertexID[_tets[x].vertices[y]];
        if (index < _unconstrainedSize)
        {
          _RCopies[id](index * 3)     += forces[y][0];
          _RCopies[id](index * 3 + 1) += forces[y][1];
          _RCopies[id](index * 3 + 2) += forces[y][2];
        }
      }
    }
  }

  // merge all the copies
  for (int x = 0; x < _totalCores; x++)
    _R += _RCopies[x];

  return _R;
  */
}

//////////////////////////////////////////////////////////////////////
// Evaluate internal forces for Newmark, but with precached
// diagonalizations
//////////////////////////////////////////////////////////////////////
VECTOR& TET_MESH::generateInternalForces(vector<MATRIX3>& Us,
                                         vector<MATRIX3>& Fhats,
                                         vector<MATRIX3>& Vs) 
{
  // wipe previous internal forces
  _R.clear();

  // populate the forces
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    // compute the forces
    VEC3F forces[4];
    const VEC3F* b = _tets[x].b();
    int materialIndex = _tets[x].materialIndex();
    
    // cast explicitly to invertible
    INVERTIBLE* material = (INVERTIBLE*)_materials[materialIndex];
  
    MATRIX3 firstPK = material->firstPiolaKirchhoff(Us[x], Fhats[x], Vs[x]);
    forces[0] = firstPK * b[0];
    forces[1] = firstPK * b[1];
    forces[2] = firstPK * b[2];
    forces[3] = firstPK * b[3];

    // add forces to corresponding vertices
    for (int y = 0; y < 4; y++)
    {
      int index = _vertexID[_tets[x].vertices[y]];
      if (index < _unconstrainedSize)
      {
        _R(index * 3)     += forces[y][0];
        _R(index * 3 + 1) += forces[y][1];
        _R(index * 3 + 2) += forces[y][2];
      }
    }
  }
  return _R;
}

//////////////////////////////////////////////////////////////////////
// Evaluate stiffness for Newmark
//////////////////////////////////////////////////////////////////////
MATRIX& TET_MESH::generateStiffnessMatrix()
{
  // Check that this isn't the first call and _K has already been
  // allocated. This check has been moved here because _K is a
  // giant memory hog, and explicit integration doesn't even need it.
  if (_K.rows() != this->TET_MESH::rank())
  {
    int rank = this->TET_MESH::rank();
    _K.resizeAndWipe(rank, rank);
  }
  
  // stomp old stiffness matrix
  _K.clear();

  for (unsigned int x = 0; x < _tets.size(); x++)
	{
    int materialIndex = _tets[x].materialIndex();
    
    // get the stiffness
    MATRIX stiffness = _materials[materialIndex]->stiffness(_tets[x]);

    // get the IDs for the vertices
		int nodeIDs[4];
    for (int i = 0; i < 4; i++)
    {
      nodeIDs[i] = _vertexID[_tets[x].vertices[i]];
      if (nodeIDs[i] >= _unconstrainedSize)
        nodeIDs[i] = -1;
    }

    // assemble the final stiffness matrix
		for (int i = 0; i < 4; i++)
		{
			if (nodeIDs[i] == -1) continue;
			for (int j = 0; j < 4; j++)
			{
        // copy 3x3 interactions (ie jacobian of one vertex wrt another)
				if (nodeIDs[j] == -1) continue;
				for (int id = 0; id < 3; id++)
					for (int jd = 0; jd < 3; jd++)
            _K(nodeIDs[i] * 3 + id, nodeIDs[j] * 3 + jd) -= stiffness(i * 3 + id, j * 3 + jd);
			}
		}
	}
  
  return _K;
}

//////////////////////////////////////////////////////////////////////
// update the vertices on the mesh surface using the full update
// vector _x
//////////////////////////////////////////////////////////////////////
void TET_MESH::updateSurfaceMesh()
{
  // total surface vertices
  int totalVertices = _surfaceVertices.size();
  for (int x = 0; x < totalVertices; x++)
  {
    VEC3F& restPose = restVertex(_surfaceVertices[x]);
    VEC3F& deformed = (*_surfaceVertices[x]);

    int index = 3 * _vertexID[_surfaceVertices[x]];

    if (index < _unconstrainedSize)
    {
      deformed[0] = restPose[0] + _x(index);
      deformed[1] = restPose[1] + _x(index + 1);
      deformed[2] = restPose[2] + _x(index + 2);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// update all the vertices on the mesh 
//////////////////////////////////////////////////////////////////////
void TET_MESH::updateFullMesh()
{
  assert(_unconstrainedSize * 3 == _x.size());

  // total surface vertices
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F& restPose = _restPose[x];
    VEC3F& deformed = _vertices[x];
    int index = 3 * x;

    deformed[0] = restPose[0] + _x(index);
    deformed[1] = restPose[1] + _x(index + 1);
    deformed[2] = restPose[2] + _x(index + 2);
  }
}

//////////////////////////////////////////////////////////////////////
// update all the vertices on the mesh with a temporary position
//////////////////////////////////////////////////////////////////////
void TET_MESH::updateFullMesh(VECTOR& position)
{
  // total surface vertices
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F& restPose = _restPose[x];
    VEC3F& deformed = _vertices[x];
    int index = 3 * x;

    deformed[0] = restPose[0] + position(index);
    deformed[1] = restPose[1] + position(index + 1);
    deformed[2] = restPose[2] + position(index + 2);
  }
}

//////////////////////////////////////////////////////////////////////
// reset the mass matrix to a specific mass
//////////////////////////////////////////////////////////////////////
void TET_MESH::resetMasses(Real mass)
{
  for (int x = 0; x < _masses.rows(); x++)
    _masses(x,x) = mass;

  _totalMass = _masses.sum() / 3;
}

//////////////////////////////////////////////////////////////////////
// Translate all vertices so that the center of mass is at the
// origin
//////////////////////////////////////////////////////////////////////
void TET_MESH::centerOfMassToOrigin()
{
  // compute the rest center of mass
  VEC3F restCenterOfMass;
  Real totalMass = 0;
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    restCenterOfMass += mass(x) * _restPose[x];
    totalMass += mass(x);
  }
  
  for (int x = 0; x < _unconstrainedSize; x++)
    _restPose[x] -= restCenterOfMass;

  // propagate the centered rest mesh
  TET_MESH::updateFullMesh();

  // compute the center of mass
  VEC3F centerOfMass;
  for (int x = 0; x < _unconstrainedSize; x++)
    centerOfMass += mass(x) * _vertices[x];
  centerOfMass *= 1.0 / totalMass;

  // translate the vertices and the rest pose to the center of mass
  for (int x = 0; x < _unconstrainedSize; x++)
    _vertices[x] -= centerOfMass;

  // propagate change to the deformation vector as well
  TET_MESH::recoverX();
 
  assert(TET_MESH::centerOfMassIsZero());
}

//////////////////////////////////////////////////////////////////////
// dump the stiffness matrix to a Matlab file
//////////////////////////////////////////////////////////////////////
void TET_MESH::stiffnessToMatlab(const char* filename)
{
  SPARSE_MATRIX& K = generateSparseStiffnessMatrix();
  K.writeToMatlab(filename);
}

//////////////////////////////////////////////////////////////////////
// calculate and return the center of mass
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeCenterOfMass()
{
  _centerOfMass.clear();
  Real totalMass = 0.0;
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    _centerOfMass += mass(x) * _vertices[x];
    totalMass += mass(x);
  }
  _centerOfMass *= 1.0 / totalMass;
}

//////////////////////////////////////////////////////////////////////
// calculate and return the center of mass
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeRestCenterOfMass()
{
  _restCenterOfMass.clear();
  Real totalMass = 0.0;
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    _restCenterOfMass += _restPose[x] * mass(x);
    totalMass += mass(x);
  }
  _restCenterOfMass *= 1.0 / totalMass;
}

//////////////////////////////////////////////////////////////////////
// compute and return the one-ring of a vertex
//////////////////////////////////////////////////////////////////////
void TET_MESH::oneRing(VEC3F* vertex, vector<VEC3F*>& result)
{
  // get all the tets the vertex is a member of
  vector<int>& membership = _tetMembership[vertex];

  // hash all the vertices in the tets -- repeated
  // vertices will still only be counter once.
  map<VEC3F*, bool> hash;
  for (unsigned int x = 0; x < membership.size(); x++)
  {
    TET& tet = _tets[membership[x]];

    for (int y = 0; y < 4; y++)
      // make sure it is not the current vertex
      if (tet.vertices[y] != vertex)
        hash[tet.vertices[y]] = true;
  }

  // scan the hash out into a flat array
  map<VEC3F*, bool>::iterator iter;
  result.clear();
  for (iter = hash.begin(); iter != hash.end(); iter++)
    result.push_back(iter->first);
}

//////////////////////////////////////////////////////////////////////
// compute and return the one-ring of a tet
//////////////////////////////////////////////////////////////////////
void TET_MESH::oneRing(int tetID, vector<VEC3F*>& result)
{
  // get the tets that the vertices are a member of
  vector<int> membership[4];
  TET& tet = _tets[tetID];
  membership[0] = _tetMembership[tet.vertices[0]];
  membership[1] = _tetMembership[tet.vertices[1]];
  membership[2] = _tetMembership[tet.vertices[2]];
  membership[3] = _tetMembership[tet.vertices[3]];
  
  // get the outlier for the 0,1,3 face
  VEC3F* outlier;
  outlier = outlierVertex(membership, tetID, 0, 1, 3);
  if (outlier != NULL) result.push_back(outlier);

  // get the tet for the 0,2,3 face
  outlier = outlierVertex(membership, tetID, 0, 2, 3);
  if (outlier != NULL) result.push_back(outlier);

  // get the tet for the 0,1,2 face
  outlier = outlierVertex(membership, tetID, 0, 1, 2);
  if (outlier != NULL) result.push_back(outlier);

  // get the tet for the 1,2,3 face
  outlier = outlierVertex(membership, tetID, 1, 2, 3);
  if (outlier != NULL) result.push_back(outlier);
}

//////////////////////////////////////////////////////////////////////
// compute and return the one-ring of a tet
//////////////////////////////////////////////////////////////////////
void TET_MESH::vertexOneRing(VEC3F* vertex, vector<int>& result)
{
  // get the tets that the vertices are a member of
  vector<int> membership;
  membership = _tetMembership[vertex];
  
  for (unsigned int x = 0; x < membership.size(); x++)
    result.push_back(membership[x]);
}

//////////////////////////////////////////////////////////////////////
// compute and return the one-ring of a tet
//////////////////////////////////////////////////////////////////////
void TET_MESH::oneRing(VEC3F* vertex, vector<TET*>& result)
{
  vector<int>& tetMembers = _tetMembership[vertex];
  for (unsigned int x = 0; x < tetMembers.size(); x++)
    result.push_back(&_tets[tetMembers[x]]);
}

//////////////////////////////////////////////////////////////////////
// helper function for oneRing
//
// Given the tetMembership lists of all the vertices in in a tet,
// find the other tet that also contains v0, v1, and v2, and return the
// vertex it contains that is not v0, v1, or v2.
//////////////////////////////////////////////////////////////////////
VEC3F* TET_MESH::outlierVertex(vector<int>* membership, int currentTet, int v0, int v1, int v2)
{
  // get the tet for the 0,1,3 face
  int foundID = -1;
  for (unsigned int x = 0; x < membership[v0].size(); x++)
  {
    int tetID = membership[v0][x];

    // they will have the current tet in common, so make sure to
    // skip it
    if (tetID == currentTet) continue;
  
    // if the tet is a member of the two other vertices as well.
    // found should equal 2 by the end.
    int found = 0;
    for (unsigned int y = 0; y < membership[v1].size(); y++)
    {
      if (membership[v1][y] == tetID)
      {
        found++;
        break;
      }
      // list is sorted, so break if we've gone past the index number
      if (membership[v1][y] > tetID)
        break;
    }
    
    // if the tet is a member of the two other vertices as well.
    // found should equal 2 by the end.
    for (unsigned int y = 0; y < membership[v2].size(); y++)
    {
      if (membership[v2][y] == tetID)
      {
        found++;
        break;
      }
      // list is sorted, so break if we've gone past the index number
      if (membership[v2][y] > tetID)
        break;
    }
    if (found == 2)
    {
      foundID = tetID;
      break;
    }
  }
  // go through all the vertices in the found tet, and store
  // the one that isn't on the current tet
  if (foundID != -1)
    for (int x = 0; x < 4; x++)
    {
      // check if the vertex is in the current tet
      if (_tets[foundID].vertices[x] == _tets[currentTet].vertices[0]) continue;
      if (_tets[foundID].vertices[x] == _tets[currentTet].vertices[1]) continue;
      if (_tets[foundID].vertices[x] == _tets[currentTet].vertices[2]) continue;
      if (_tets[foundID].vertices[x] == _tets[currentTet].vertices[3]) continue;

      // store the current vertex
      return _tets[foundID].vertices[x];
    }

  return NULL;
}

//////////////////////////////////////////////////////////////////////
// Write out a partitioned version of this mesh -- assumes METIS
// has been run on the mesh and the partition field of all the tets
// has been assigned
//////////////////////////////////////////////////////////////////////
void TET_MESH::writePartitions(string prefix)
{
  // find out how many partitions there are
  int maxPartition = 0;
  for (unsigned int x = 0; x < _tets.size(); x++)
    maxPartition = (_tets[x].partition > maxPartition) ? _tets[x].partition : maxPartition;

  int totalPartitions = maxPartition + 1;

  // return if there is only one
  if (totalPartitions == 0)
  {
    cout << " There is only one partition! Run METIS!" << endl;
    return;
  }

  for (unsigned int x = 0; x < _tets.size(); x++)
    if (_tets[x].partition < 0)
    {
      cout << " Some tets have not been assigned a partition! " << endl;
      exit(0);
    }
  
  // keep around new vertex IDs for clone detection later
  map<VEC3F*, int>* newVertexIDs;
  newVertexIDs = new map<VEC3F*, int>[totalPartitions];

  // for each partition
  for (int partition = 0; partition < totalPartitions; partition++)
  {
    // create vector of tets in this partition
    vector<TET*> newTets;
    for (unsigned int x = 0; x < _tets.size(); x++)
      if (_tets[x].partition == partition)
        newTets.push_back(&_tets[x]);

    // create a list of vertices used in this partition
    list<VEC3F*> newVerticesList;
    map<VEC3F*, bool> vertexUsed;
    int totalUnconstrained = 0;
    int totalConstrained = 0;
    for (unsigned int x = 0; x < newTets.size(); x++)
      for (int y = 0; y < 4; y++)
      {
        // see if this vertex has been added to the list already
        VEC3F* vertex = newTets[x]->vertices[y];
        if (vertexUsed.find(vertex) == vertexUsed.end())
        {
          // make sure the constrained vertices always go to the
          // back of the list by pushing them to the back and
          // the unconstrained to the front
          if (_constrained[vertex])
          {
            newVerticesList.push_back(vertex);
            totalConstrained++;
          }
          else
          {
            newVerticesList.push_front(vertex);
            totalUnconstrained++;
          }
          // hash the vertex so it isn't double-counter later
          vertexUsed[vertex] = true;
        }
      }

    // copy the list to a vector, and create a hash between the
    // vertex addresses and their position in this new vector
    vector<VEC3F*> newVertices;
    list<VEC3F*>::iterator iter;
    int index = 0;
    for (iter = newVerticesList.begin(); iter != newVerticesList.end(); iter++)
    {
      newVertices.push_back(*iter);
      newVertexIDs[partition][*iter] = index;
      index++;
    }

    // open file for dump
    char buffer[256];
    sprintf(buffer, "%i", partition);

    //string partitionFilename = _filename;
    string partitionFilename = prefix;
    if (prefix.length() == 0)
      partitionFilename = _filename;
    partitionFilename += string(".partition.");
    partitionFilename += string(buffer);
    FILE* file = fopen(partitionFilename.c_str(), "wb");
    if (file == NULL)
      printf("Partition file %s could not be opened!\n", partitionFilename.c_str());
    cout << "Writing partition file: " << partitionFilename.c_str() << endl;

    // output the vertices
    fwrite((void*)&totalUnconstrained, sizeof(int), 1, file);
    fwrite((void*)&totalConstrained, sizeof(int), 1, file);

    // write out vertex positions
    for (unsigned int x = 0; x < newVertices.size(); x++)
    {
      VEC3F node = *(newVertices[x]);
      fwrite((void*)&(node), sizeof(Real), 3, file);
    }
    // write out the tet indices
    int totalTets = newTets.size();
    fwrite((void*)&totalTets, sizeof(int), 1, file);
    for (unsigned int x = 0; x < newTets.size(); x++)
      for (int y = 0; y < 4; y++)
      {
        int vertexID = newVertexIDs[partition][newTets[x]->vertices[y]];
        fwrite((void*)&vertexID, sizeof(int), 1, file);
      }

    fclose(file);
  }

  // write out the cloned vertices
  writeClones(totalPartitions, newVertexIDs, prefix);

  // write the correspondence between the old and new vertices
  writeVertexCorrespondence(totalPartitions, newVertexIDs, prefix);

  delete[] newVertexIDs;
}

//////////////////////////////////////////////////////////////////////
// helper function for writePartitions -- writes out the file
// that describes which vertices are cloned across partitions
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeClones(int totalPartitions, map<VEC3F*, int>* newVertexIDs, string prefix)
{
  // initialize clones -- an n x n grid of pairs
  vector<pair<int, int> >** clones;
  clones = new vector<pair<int, int> >*[totalPartitions];
  for (int x = 0; x < totalPartitions; x++)
    clones[x] = new vector<pair<int, int> >[totalPartitions];
  
  // record which vertices span more than one partition
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    // vertex we are looking at
    VEC3F* vertex = &_vertices[x];
    
    // get its tet membership list
    vector<int>& membership = _tetMembership[vertex];

    // count how many partitions the vertex belongs to
    map<int, bool> count;
    for (unsigned int y = 0; y < membership.size(); y++)
    {
      TET& tet = _tets[membership[y]];
      count[tet.partition] = true;
    }

    // if it just belongs to one partition, move on
    if (count.size() <= 1) continue;

    // copy partition list to vector for convenient random access
    vector<int> partitions;
    map<int, bool>::iterator iter;
    for (iter = count.begin(); iter != count.end(); iter++)
      partitions.push_back(iter->first);

    // store all the pairwise interactions
    for (unsigned int y = 0; y < partitions.size(); y++)
      for (unsigned int z = y + 1; z < partitions.size(); z++)
      {
        // call the two partitions "first" and "second" just so
        // we have something to call them
        int first = partitions[y];
        int second = partitions[z];

        // what index # does this vertex have in the first partition?
        int firstID = newVertexIDs[first][vertex];

        // what index # does this vertex have in the second partition?
        int secondID = newVertexIDs[second][vertex];
       
        // store the interaction symmetrically
        clones[first][second].push_back(pair<int,int>(firstID, secondID));
        clones[second][first].push_back(pair<int,int>(secondID, firstID));
      }
  }

  // sanity check -- make sure that the number of clones is symmetric
  // this should be trivially true
  for (int x = 0; x < totalPartitions; x++)
    for (int y = 0; y < totalPartitions; y++)
      assert(clones[x][y].size() == clones[y][x].size());

  // write out the clones file
  string cloneFilename = prefix;
  if (cloneFilename.length() == 0)
    cloneFilename = _filename;
  cloneFilename += string(".clones");
  FILE* file = fopen(cloneFilename.c_str(), "wb");
  for (int x = 0; x < totalPartitions; x++)
  {
    for (int y = 0; y < totalPartitions; y++)
    {
      // write out how many interactions there are
      int totalClones = clones[x][y].size();
      fwrite((void*)&totalClones, sizeof(int), 1, file);

      // write out the pairs
      for (int z = 0; z < totalClones; z++)
      {
        pair<int, int> clone = clones[x][y][z];
        int first = clone.first;
        int second = clone.second;
        fwrite((void*)&first, sizeof(int), 1, file);
        fwrite((void*)&second, sizeof(int), 1, file);
      }
    }
  }
  fclose(file);

  // clean up clones
  for (int x = 0; x < totalPartitions; x++)
    delete[] clones[x];
  delete[] clones;
}

//////////////////////////////////////////////////////////////////////
// helper function for writePartitions -- writes out the file
// that translates the partition vertex index to an index in the
// original unpartitioned mesh
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeVertexCorrespondence(int totalPartitions, map<VEC3F*, int>* newVertexIDs, string prefix)
{
  for (int x = 0; x < totalPartitions; x++)
  {
    // lookup what the original vertexID is 
    map<int,int> inverseID;
    map<VEC3F*, int>::iterator iter;
    for (iter = newVertexIDs[x].begin(); iter != newVertexIDs[x].end(); iter++)
    {
      VEC3F* vertex = iter->first;
      int newID = iter->second;
      int oldID = _vertexID[vertex];
      inverseID[newID] = oldID;
    }

    // open the correspondence file
    char buffer[256];
    sprintf(buffer, "%i", x);
    string partitionFilename = prefix;
    if (partitionFilename.length() == 0)
      partitionFilename = _filename;
    partitionFilename += string(".partition.");
    partitionFilename += string(buffer);
    partitionFilename += string(".correspondence");
    FILE* file = fopen(partitionFilename.c_str(), "wb");

    // dump out the correspondence
    map<int,int>::iterator dumpIter;
    int sanity = 0;
    for (dumpIter = inverseID.begin(); dumpIter != inverseID.end(); dumpIter++)
    {
      int newID = dumpIter->first;
      int oldID = dumpIter->second;
      fwrite((void*)&oldID, sizeof(int), 1, file);
      if (newID != sanity)
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << "The correspondence table is missing a vertex!" << endl;
        cout << " Sanity check is on vertex: " << sanity << endl;
        cout << " Map is on vertex: " << newID << endl;
      }
      sanity++;
    }
    
    fclose(file);
  }
}

//////////////////////////////////////////////////////////////////////
// Write out the deformed mesh
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeDeformedMesh(int timestep, string path)
{
  char buffer[256];
  sprintf(buffer, "%i", timestep);
  string filename = path;
  filename += string("deformed.");
  filename += string(buffer);
  filename += string(".tetmesh");
  cout << " Writing deformation file " << filename.c_str() << endl;
  
  FILE* file = fopen(filename.c_str(), "wb");
  if (file == NULL)
    printf("Filename %s could not be opened!\n", filename.c_str());

  // output vertex array sizes
  fwrite((void*)&_unconstrainedSize, sizeof(int), 1, file);
  fwrite((void*)&_constrainedSize, sizeof(int), 1, file);

  // write out vertex positions
  for (unsigned int x = 0; x < _vertices.size(); x++)
    fwrite((void*)&(_vertices[x]), sizeof(Real), 3, file);

  // output tet vertex lists
  int totalTets = _tets.size();
  fwrite((void*)&totalTets, sizeof(int), 1, file);

  for (int x = 0; x < totalTets; x++)
  {
    TET& tet = _tets[x];
    for (int y = 0; y < 4; y++)
    {
      int index = _vertexID[tet.vertices[y]];
      fwrite((void*)&index, sizeof(int), 1, file);
    }
  }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Read in the deformed mesh
//////////////////////////////////////////////////////////////////////
bool TET_MESH::readDeformedMesh(int timestep, string path, bool verbose)
{
  char buffer[256];
  sprintf(buffer, "%i", timestep);
  string filename = path;
  filename += string("deformed.");
  filename += string(buffer);
  filename += string(".tetmesh");
  if (verbose)
    cout << " Reading deformation file " << filename.c_str() << endl;
  
  FILE* file = fopen(filename.c_str(), "rb");
  if (file == NULL)
    printf("Filename %s could not be opened!\n", filename.c_str());

  size_t check;

  // read in dummy dimensions
  int dummy;
  check = fread((void*)&dummy, sizeof(int), 1, file);
  if (check != 1)
  {
    cout << " Error encountered reading first size of " << filename.c_str() << endl;
    cout << " File pointer: " << file << endl;
    fclose(file);
    return false;
  }
  check = fread((void*)&dummy, sizeof(int), 1, file);
  if (check != 1)
  {
    cout << " Error encountered reading second size of " << filename.c_str() << endl;
    cout << " File pointer: " << file << endl;
    fclose(file);
    return false;
  }

  // read in vertex positions
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    check = fread((void*)&(_vertices[x][0]), sizeof(Real), 1, file);
    if (check != 1)
    {
      cout << " Error encountered reading vertex x of " << filename.c_str() << endl;
      fclose(file);
      return false;
    }
    check = fread((void*)&(_vertices[x][1]), sizeof(Real), 1, file);
    if (check != 1)
    {
      cout << " Error encountered reading vertex y of " << filename.c_str() << endl;
      fclose(file);
      return false;
    }
    check = fread((void*)&(_vertices[x][2]), sizeof(Real), 1, file);
    if (check != 1) 
    {
      cout << " Error encountered reading vertex z of " << filename.c_str() << endl;
      fclose(file);
      return false;
    }
  }
  if (verbose)
    cout << " File read successfully." << endl;

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// displacements of the constrained nodes
//////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::constraintDisplacements()
{
  int totalConstraints = 3 * (_vertices.size() - _unconstrainedSize);
  VECTOR displacements(totalConstraints);

  int i = 0;
  for (unsigned int x = _unconstrainedSize; x < _vertices.size(); x++, i++)
  {
    displacements[3 * i] = _vertices[x][0] - _restPose[x][0];
    displacements[3 * i + 1] = _vertices[x][1] - _restPose[x][1];
    displacements[3 * i + 2] = _vertices[x][2] - _restPose[x][2];
  }

  return displacements;
}

//////////////////////////////////////////////////////////////////////
// animate the mesh somehow
//////////////////////////////////////////////////////////////////////
void TET_MESH::animate(int timestep)
{
  // increment the max y nodes
  vector<VEC3F*> topNodes;
  for (unsigned int x = _unconstrainedSize; x < _vertices.size(); x++)
  {
    // do nothing if it is a bottom node
    if (_vertices[x][1] < 0.5) continue;
    VEC3F& vertex = _vertices[x];

    // the rotation
    //if (timestep < 72)
    if (timestep < 36)
    {
      vertex -= _restCenterOfMass;
      MATRIX3 rotate;
      if (timestep < 18 || timestep >= 54)
        rotate = MATRIX3::rotation(VEC3F(0,1,0), 0.1);
      if (timestep >= 18 && timestep < 54)
        rotate = MATRIX3::rotation(VEC3F(0,1,0), -0.1);
      vertex = rotate * vertex;
      vertex += _restCenterOfMass;
      continue;
    }

    // the extension
    if (timestep < 56)
    {
      if (timestep < 46)
      {
        vertex[1] += 0.1;
        continue;
      }
      vertex[1] -= 0.1;
    }
    
    /*
    if (timestep < 36 * 2)
    {
      vertex -= _restCenterOfMass;
      MATRIX3 rotate;
      if (timestep < 18 * 2 || timestep >= 54 * 2)
        rotate = MATRIX3::rotation(VEC3F(0,1,0), 0.1/2);
      if (timestep >= 18 * 2 && timestep < 54 * 2)
        rotate = MATRIX3::rotation(VEC3F(0,1,0), -0.1/2);
      vertex = rotate * vertex;
      vertex += _restCenterOfMass;
      continue;
    }
    */

  }
}

//////////////////////////////////////////////////////////////////////
// get the bounding box
//////////////////////////////////////////////////////////////////////
void TET_MESH::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  boundingBox(mins[0], maxs[0],
              mins[1], maxs[1],
              mins[2], maxs[2]);
}

//////////////////////////////////////////////////////////////////////
// get the bounding box
//////////////////////////////////////////////////////////////////////
void TET_MESH::boundingBox(Real& xMin, Real& xMax, 
                           Real& yMin, Real& yMax, 
                           Real& zMin, Real& zMax)
{
  if (_vertices.size() == 0) return;
  
  xMin = _vertices[0][0];
  xMax = _vertices[0][0];
  yMin = _vertices[0][1];
  yMax = _vertices[0][1];
  zMin = _vertices[0][2];
  zMax = _vertices[0][2];

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    if (_vertices[x][0] < xMin) xMin = _vertices[x][0];
    if (_vertices[x][0] > xMax) xMax = _vertices[x][0];
    if (_vertices[x][1] < yMin) yMin = _vertices[x][1];
    if (_vertices[x][1] > yMax) yMax = _vertices[x][1];
    if (_vertices[x][2] < zMin) zMin = _vertices[x][2];
    if (_vertices[x][2] > zMax) zMax = _vertices[x][2];
  }
}

//////////////////////////////////////////////////////////////////////
// Compute collision nodes
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeCollisionNodes()
{
  // create a hash of all the surface nodes
  map<VEC3F*, bool> surfaceHash;
  for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
    surfaceHash[_surfaceVertices[x]] = true;
 
  // track which nodes have been added
  map<VEC3F*, bool> addedNodes;

  // wipe any previous results
  _collisionNodes.clear();
  
  // for each surface vertex, find all the tets it is a member of
  // and see if it has a sub-dermal node
  for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
  {
    vector<int>& tetMembership = _tetMembership[_surfaceVertices[x]];

    for (unsigned int y = 0; y < tetMembership.size(); y++)
    {
      TET& tet = _tets[tetMembership[y]];

      for (int z = 0; z < 4; z++)
      {
        if (surfaceHash.find(tet.vertices[z]) == surfaceHash.end() &&
            addedNodes.find(tet.vertices[z]) == addedNodes.end())
        {
          _collisionNodes.push_back(tet.vertices[z]);
          addedNodes[tet.vertices[z]] = true;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Output surface triangles to RenderMan
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawSurfaceToRenderMan()
{
  //int x;

#ifdef USING_RENDERMAN
  int nfaces = _surfaceFaces.size();
  computeSurfaceNormals();

  RtPoint* P = new RtPoint[_surfaceVertices.size()];
  for (unsigned int y = 0; y < _surfaceVertices.size(); y++)
  {
    P[y][0] = (*_surfaceVertices[y])[0];
    P[y][1] = (*_surfaceVertices[y])[1];
    P[y][2] = (*_surfaceVertices[y])[2];
  }

  /*
  RtPoint* N = new RtPoint[_surfaceNormals.size()];
  for (int y = 0; y < _surfaceNormals.size(); y++)
  {
    N[y][0] = _surfaceNormals[y][0];
    N[y][1] = _surfaceNormals[y][1];
    N[y][2] = _surfaceNormals[y][2];
  }
  */

  RtInt* faceIndices= new RtInt[3 * nfaces];
  for (int y = 0; y < nfaces; y++)
  {
    // construct the triangle
    TRIANGLE triangle = _tets[_surfaceFaces[y].first].face(_surfaceFaces[y].second);
  
    faceIndices[3 * y] = _surfaceVertexID[triangle.vertex(0)]; 
    faceIndices[3 * y + 1] = _surfaceVertexID[triangle.vertex(1)]; 
    faceIndices[3 * y + 2] = _surfaceVertexID[triangle.vertex(2)]; 
  }

  // init all to triangles
  RtInt* nvertices = new RtInt[nfaces];
  for (int x = 0; x < nfaces; x++)
    nvertices[x] = 3;
  
  //RiPointsPolygons(nfaces, nvertices, faceIndices, RI_P, P, RI_N, (RtPointer)N, RI_NULL);
  RiPointsPolygons(nfaces, nvertices, faceIndices, RI_P, P, RI_NULL);

  delete[] nvertices;
  delete[] P;
  //delete[] N;
  delete[] faceIndices;
#endif

  /*
  int x;
  int nfaces = _surfaceFaces.size();

#ifdef USING_RENDERMAN
  // init point locations and normals
  RtPoint* P = new RtPoint[nfaces * 3];
  RtPoint* N = new RtPoint[nfaces * 3];
  for (x = 0; x < nfaces; x++)
  {
    // construct the triangle
    TRIANGLE triangle = _tets[_surfaceFaces[x].first].face(_surfaceFaces[x].second);
   
    // get the normal
    VEC3F& normal = triangle.normal();
    
    for (int y = 0; y < 3; y++)
    { 
      VEC3F& vertex = *(triangle.vertex(y));
      P[3 * x + y][0] = vertex[0];
      P[3 * x + y][1] = vertex[1];
      P[3 * x + y][2] = vertex[2];

      N[3 * x + y][0] = normal[0];
      N[3 * x + y][1] = normal[1];
      N[3 * x + y][2] = normal[2];
    }
  }

  // init all to triangles
  RtInt* nvertices = new RtInt[nfaces];
  for (x = 0; x < nfaces; x++)
    nvertices[x] = 3;
  
  // init faces
  RtInt* rifaces = new RtInt[3 * nfaces];
  for (x = 0; x < nfaces; x++)
  {
    rifaces[x * 3]     = 3 * x;
    rifaces[x * 3 + 1] = 3 * x + 1;
    rifaces[x * 3 + 2] = 3 * x + 2;
  }
  
  RiPointsPolygons(nfaces, nvertices, rifaces, RI_P, P, RI_N, (RtPointer)N, RI_NULL);

  delete[] nvertices;
  delete[] P;
  delete[] N;
  delete[] rifaces;
#endif
*/
}

//////////////////////////////////////////////////////////////////////
// Output surface forces to RenderMan
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawSurfaceForcesToRenderMan()
{
#ifdef USING_RENDERMAN
  // draw outlines on a per-tet basis
  RtColor black;
  black[0] = 1.0;
  black[1] = 0.0;
  black[2] = 0.0;
  RiColor(black);
  for (unsigned int i = 0; i < _vertices.size(); i++)
  {
    //if (_vertexOnSurface.find(&_vertices[i]) == _vertexOnSurface.end())
    //  continue;

    int vertexID = _vertexID[&_vertices[i]];

    // get the force
    VEC3F force;
    if (vertexID < _unconstrainedSize)
    {
      force[0] = _R[3 * vertexID];
      force[1] = _R[3 * vertexID + 1];
      force[2] = _R[3 * vertexID + 2];

      force *= 20.0;
    }

    // init all to triangles
    int nfaces = 1;
    RtInt* nvertices = new RtInt[nfaces];
    for (int x = 0; x < nfaces; x++)
      nvertices[x] = 2;

    // init point locations and normals
    RtPoint* P = new RtPoint[2];
    P[0][0] = _vertices[i][0];
    P[0][1] = _vertices[i][1];
    P[0][2] = _vertices[i][2];
    P[1][0] = _vertices[i][0] + force[0];
    P[1][1] = _vertices[i][1] + force[1];
    P[1][2] = _vertices[i][2] + force[2];

    RtFloat width = 0.002;
    RtToken type = "linear";
    RtToken wrap = "periodic";
    RiDeclare("width", "constant float");
    RiCurves(type, nfaces, nvertices, wrap, RI_P, (RtPointer)P, RI_CONSTANTWIDTH, &width ,RI_NULL);
     
    delete[] P; 
    delete[] nvertices;
  }

#endif
}

//////////////////////////////////////////////////////////////////////
// Output outlines of triangles to renderman
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawOutlinesToRenderMan()
{
#ifdef USING_RENDERMAN
  int x;

  // init all to triangles
  int nfaces = _tets.size() * 4;
  RtInt* nvertices = new RtInt[nfaces];
  for (x = 0; x < nfaces; x++)
    nvertices[x] = 3;
  
  // init point locations and normals
  RtPoint* P = new RtPoint[nfaces * 3];
  int i = 0;
  for (unsigned int x = 0; x < _tets.size(); x++)
    for (int y = 0; y < 4; y++, i++)
    {
      TRIANGLE triangle = _tets[x].face(y);
      VEC3F v0 = *(triangle.vertex(0));
      VEC3F v1 = *(triangle.vertex(1));
      VEC3F v2 = *(triangle.vertex(2));

      P[i * 3][0] = v0[0];
      P[i * 3][1] = v0[1];
      P[i * 3][2] = v0[2];
      P[i * 3 + 1][0] = v1[0];
      P[i * 3 + 1][1] = v1[1];
      P[i * 3 + 1][2] = v1[2];
      P[i * 3 + 2][0] = v2[0];
      P[i * 3 + 2][1] = v2[1];
      P[i * 3 + 2][2] = v2[2];
    }

  RtFloat width = 0.002;
  RtToken type = "linear";
  RtToken wrap = "periodic";
  RiDeclare("width", "constant float");
  RiCurves(type, nfaces, nvertices, wrap, RI_P, (RtPointer)P, RI_CONSTANTWIDTH, &width ,RI_NULL);

  delete[] P;
  delete[] nvertices;
#endif
}

//////////////////////////////////////////////////////////////////////
// rainbow ramp
//////////////////////////////////////////////////////////////////////
void TET_MESH::rainbowRamp(Real input, float color[3])
{
  VEC3F red;
  VEC3F green;
  VEC3F blue;
  VEC3F yellow;
  VEC3F cyan;
  VEC3F magenta;
  VEC3F black;

  red[0] = 1.0;
  red[1] = 0.0;
  red[2] = 0.0;

  yellow[0] = 1.0;
  yellow[1] = 1.0;
  yellow[2] = 0.0;

  magenta[0] = 1.0;
  magenta[1] = 0.0;
  magenta[2] = 1.0;

  black[0] = 0.0;
  black[1] = 0.0;
  black[2] = 0.0;

  cyan[0] = 0.0;
  cyan[1] = 1.0;
  cyan[2] = 1.0;

  green[0] = 0.0;
  green[1] = 1.0;
  green[2] = 0.0;

  blue[0] = 0.0;
  blue[1] = 0.0;
  blue[2] = 1.0;

  VEC3F left;
  VEC3F right;
  /*
  if (input < 0.5)
  {
    input *= 2.0;
    left = red;
    right = green;
  }
  else if (input >= 0.5)
  {
    input -= 0.5;
    input *= 2.0;
    left = green;
    right = blue;
  }
  */
  if (input > 0.8)
  {
    input -= 0.8;
    input *= 5.0;
    right = red;
    left = yellow;
  }
  else if (input > 0.6)
  {
    input -= 0.6;
    input *= 5.0;
    right = yellow;
    left = green;
  }
  else if (input > 0.4)
  {
    input -= 0.4;
    input *= 5.0;
    right = green;
    left = cyan;
  }
  else if (input > 0.2)
  {
    input -= 0.2;
    input *= 5.0;
    right = cyan;
    left = blue;
  }
  else
  {
    input *= 5.0;
    right = blue;
    left = black;
  }

  color[0] = input * right[0] + (1.0 - input) * left[0];
  color[1] = input * right[1] + (1.0 - input) * left[1];
  color[2] = input * right[2] + (1.0 - input) * left[2];
}

//////////////////////////////////////////////////////////////////////
// Output colored tets to OGL
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawFirstPK()
{
  // find the maximum density
  Real maxMagnitude = 0.0;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    // get the force density
    MATERIAL* material = _materials[_tets[i].materialIndex()];
    MATRIX3 F = _tets[i].F();
    MATRIX3 firstPK = material->firstPiolaKirchhoff(F);
    Real magnitude = 0.0;
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++)
        magnitude += firstPK(x,y) * firstPK(x,y);
    magnitude = sqrt(magnitude);

    if (magnitude > maxMagnitude)
      maxMagnitude = magnitude;
  }

  // draw on a per-tet basis
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    // get the force density
    MATERIAL* material = _materials[_tets[i].materialIndex()];
    MATRIX3 F = _tets[i].F();
    MATRIX3 firstPK = material->firstPiolaKirchhoff(F);
    Real magnitude = 0.0;
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++)
        magnitude += firstPK(x,y) * firstPK(x,y);
    magnitude = sqrt(magnitude);
    magnitude *= 1.0 / maxMagnitude;
    magnitude += 0.00001;

    float densityColor[3];
    rainbowRamp(magnitude, densityColor);
    
    glColor4f(densityColor[0], densityColor[1], densityColor[2], 1.0);
    _tets[i].drawTriangles();
  }
}

//////////////////////////////////////////////////////////////////////
// Output colored test to RenderMan
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawColoredToRenderMan()
{
#ifdef USING_RENDERMAN
  Real sliceZ = 1.0;
  
  // find the maximum density
  Real maxMagnitude = 0.0;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    // draw a slice, not all the tets
    bool found = false;
    for (int x = 0; x < 4; x++)
      if ((*_tets[i].vertices[x])[2] > sliceZ)
        found = true; 
    if (found) continue;

    // get the force density
    MATERIAL* material = _materials[_tets[i].materialIndex()];
    MATRIX3 F = _tets[i].F();
    MATRIX3 firstPK = material->firstPiolaKirchhoff(F);
    Real magnitude = 0.0;
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++)
        magnitude += firstPK(x,y) * firstPK(x,y);
    magnitude = sqrt(magnitude);

    if (magnitude > maxMagnitude)
      maxMagnitude = magnitude;
  }

  // draw on a per-tet basis
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    // draw a slice, not all the tets
    bool found = false;
    for (int x = 0; x < 4; x++)
      if ((*_tets[i].vertices[x])[2] > sliceZ)
        found = true; 
    if (found) continue;

    // get the force density
    MATERIAL* material = _materials[_tets[i].materialIndex()];
    MATRIX3 F = _tets[i].F();
    MATRIX3 firstPK = material->firstPiolaKirchhoff(F);
    Real magnitude = 0.0;
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++)
        magnitude += firstPK(x,y) * firstPK(x,y);
    magnitude = sqrt(magnitude);
    magnitude *= 1.0 / maxMagnitude;
    magnitude += 0.00001;

    // get a black run
    //magnitude = 0.0;

    RtColor densityColor;
    rainbowRamp(magnitude, densityColor);
    RiColor(densityColor);

    // init point locations and normals
    RtPoint* P = new RtPoint[4];

    for (int x = 0; x < 4; x++)
    {
      VEC3F vertex = *(_tets[i].vertices[x]);

      P[x][0] = vertex[0];
      P[x][1] = vertex[1];
      P[x][2] = vertex[2];
    }

    // init all to triangles
    int nfaces = 4;
    RtInt* nvertices = new RtInt[nfaces];
    for (int x = 0; x < nfaces; x++)
      nvertices[x] = 3;
    
    // init faces
    RtInt* rifaces = new RtInt[3 * nfaces];
    rifaces[0] = 0;
    rifaces[1] = 1;
    rifaces[2] = 3;
    rifaces[3] = 1;
    rifaces[4] = 2;
    rifaces[5] = 3;
    rifaces[6] = 0;
    rifaces[7] = 3;
    rifaces[8] = 2;
    rifaces[9] = 1;
    rifaces[10] = 0;
    rifaces[11] = 2;
    
    RiPointsPolygons(nfaces, nvertices, rifaces, RI_P, P,  RI_NULL);

    delete[] nvertices;
    delete[] P;
    delete[] rifaces;
  }

  /*
  // draw outlines on a per-tet basis
  RtColor black;
  black[0] = 0.0;
  black[1] = 0.0;
  black[2] = 0.0;
  RiColor(black);
  for (int i = 0; i < _tets.size(); i++)
  {
    // draw a slice, not all the tets
    bool found = false;
    for (int x = 0; x < 4; x++)
      if ((*_tets[i].vertices[x])[2] > sliceZ)
        found = true; 
    if (found) continue;

    // init all to triangles
    int nfaces = 4;
    RtInt* nvertices = new RtInt[nfaces];
    for (int x = 0; x < nfaces; x++)
      nvertices[x] = 3;

    // init point locations and normals
    RtPoint* P = new RtPoint[nfaces * 3];
    Real offset = 0.001;
    for (int y = 0; y < 4; y++)
    {
      TRIANGLE triangle = _tets[i].face(y);
      VEC3F v0 = *(triangle.vertex(0));
      VEC3F v1 = *(triangle.vertex(1));
      VEC3F v2 = *(triangle.vertex(2));

      P[0][0] = v0[0];
      P[0][1] = v0[1] + offset;
      P[0][2] = v0[2];
      P[1][0] = v1[0];
      P[1][1] = v1[1] + offset;
      P[1][2] = v1[2];
      P[2][0] = v2[0];
      P[2][1] = v2[1] + offset;
      P[2][2] = v2[2];
    }

    RtFloat width = 0.002;
    RtToken type = "linear";
    RtToken wrap = "periodic";
    RiDeclare("width", "constant float");
    RiCurves(type, nfaces, nvertices, wrap, RI_P, (RtPointer)P, RI_CONSTANTWIDTH, &width ,RI_NULL);
     
    delete[] P; 
    delete[] nvertices;
  }
  */

#endif
}

//////////////////////////////////////////////////////////////////////
// Output surface triangles to RenderMan
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawExhaustiveToRenderMan()
{
#ifdef USING_RENDERMAN
  // init point locations and normals
  RtPoint* P = new RtPoint[_vertices.size()];

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    P[x][0] = _vertices[x][0];
    P[x][1] = _vertices[x][1];
    P[x][2] = _vertices[x][2];
  }

  // init all to triangles
  int nfaces = _tets.size() * 4;
  RtInt* nvertices = new RtInt[nfaces];
  for (int x = 0; x < nfaces; x++)
    nvertices[x] = 3;
  
  // init faces
  RtInt* rifaces = new RtInt[3 * nfaces];
  int i = 0;
  for (unsigned int x = 0; x < _tets.size(); x++)
    for (int y = 0; y < 4; y++, i++)
    {
      TRIANGLE triangle = _tets[x].face(y);
      rifaces[i * 3] = _vertexID[triangle.vertex(0)];
      rifaces[i * 3 + 1] = _vertexID[triangle.vertex(1)];
      rifaces[i * 3 + 2] = _vertexID[triangle.vertex(2)];
    }
  
  RiPointsPolygons(nfaces, nvertices, rifaces, RI_P, P,  RI_NULL);

  delete[] nvertices;
  delete[] P;
  delete[] rifaces;
#endif
}

//////////////////////////////////////////////////////////////////////
// Output surface triangles to RenderMan
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawZSliceToRenderMan(float zSlice)
{
//#ifndef _WIN32
#ifdef USING_RENDERMAN
  // init point locations and normals
  RtPoint* P = new RtPoint[_vertices.size()];

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    P[x][0] = _vertices[x][0];
    P[x][1] = _vertices[x][1];
    P[x][2] = _vertices[x][2];
  }

  // count the tets
  int totalTets = 0;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    bool add = false;
    for (int y = 0; y < 4; y++)
    {
      VEC3F& vertex = *(_tets[x].vertices[y]);
      if (vertex[0] > zSlice)
        add = true;
    }
    if (add) totalTets++;
  }

  // init all to triangles
  int nfaces = totalTets * 4;
  RtInt* nvertices = new RtInt[nfaces];
  for (int x = 0; x < nfaces; x++)
    nvertices[x] = 3;
  
  // init faces
  RtInt* rifaces = new RtInt[3 * nfaces];
  int i = 0;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    bool add = false;
    for (int y = 0; y < 4; y++)
    {
      VEC3F& vertex = *(_tets[x].vertices[y]);
      if (vertex[0] > zSlice)
        add = true;
    }
    if (!add) continue;

    for (int y = 0; y < 4; y++, i++)
    {
      TRIANGLE triangle = _tets[x].face(y);
      rifaces[i * 3] = _vertexID[triangle.vertex(0)];
      rifaces[i * 3 + 1] = _vertexID[triangle.vertex(1)];
      rifaces[i * 3 + 2] = _vertexID[triangle.vertex(2)];
    }
  }
  
  RiPointsPolygons(nfaces, nvertices, rifaces, RI_P, P,  RI_NULL);

  delete[] nvertices;
  delete[] P;
  delete[] rifaces;
#endif
}

//////////////////////////////////////////////////////////////////////
// constrain nodes instide this surface
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeNewConstraints(vector<VEC3F*>& newConstraints, const char* filename)
{
  // make a new hash of the constrained nodes
  map<VEC3F*, bool> constrainedHash;
  for (unsigned int x = 0; x < newConstraints.size(); x++)
    constrainedHash[newConstraints[x]] = true;

  // record the translation from the old to new vertexIDs
  map<VEC3F*, int> newID;

  // total unconstrained and constrained nodes
  int unconstrainedSize = 0;
  int constrainedSize = 0;

  // go through all the vertices, record the unconstrained ones first
  vector<VEC3F*> newVertices;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    if (constrainedHash.find(&_vertices[x]) == constrainedHash.end())
    {
      newVertices.push_back(&_vertices[x]);
      newID[&_vertices[x]] = newVertices.size() - 1;
      unconstrainedSize++;
    }

  // now record the constrained ones
  for (unsigned int x = 0; x < _vertices.size(); x++)
    if (constrainedHash.find(&_vertices[x]) != constrainedHash.end())
    {
      newVertices.push_back(&_vertices[x]);
      newID[&_vertices[x]] = newVertices.size() - 1;
      constrainedSize++;
    }

  // start dumping to a new file
  FILE* file = fopen(filename, "wb");
  
  if (file == NULL)
    printf("Filename %s could not be opened!\n", filename);

  cout << " Writing file " << filename << " ... ";

  // write out vertex array sizes
  fwrite((void*)&unconstrainedSize, sizeof(int), 1, file);
  fwrite((void*)&constrainedSize, sizeof(int), 1, file);

  // write out vertex positions
  for (unsigned int x = 0; x < newVertices.size(); x++)
  {
    double nodeDouble[3];
    nodeDouble[0] = (*newVertices[x])[0];
    nodeDouble[1] = (*newVertices[x])[1];
    nodeDouble[2] = (*newVertices[x])[2];
    fwrite((void*)nodeDouble, sizeof(double), 3, file);
  }

  // read in tet vertex lists
  int totalTets = _tets.size();
  fwrite((void*)&totalTets, sizeof(int), 1, file);

  for (int x = 0; x < totalTets; x++)
  {
    // for each node in the tet
    int indices[4];
    indices[0] = newID[_tets[x].vertices[0]];
    indices[1] = newID[_tets[x].vertices[1]];
    indices[2] = newID[_tets[x].vertices[2]];
    indices[3] = newID[_tets[x].vertices[3]];
    fwrite((void*)&indices[0], sizeof(int), 4, file);
  }

  fclose(file);
  cout << " done. " << endl;

  string faceFile = string(filename) + string(".surfacefaces");
  string vertexFile = string(filename) + string(".surfacevertices");
  faceFile = string("rm ") + faceFile;
  vertexFile = string("rm ") + vertexFile;
  system(faceFile.c_str());
  system(vertexFile.c_str());
}

//////////////////////////////////////////////////////////////////////
// embed this surface into the mesh
//////////////////////////////////////////////////////////////////////
void TET_MESH::embedMesh(OBJ* obj)
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
// normalize the embedding mesh to the size of the tet mesh
//////////////////////////////////////////////////////////////////////
void TET_MESH::normalizeEmbedding(OBJ* obj)
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

  // get the bounding box of the tet mesh
  Real tetsBoundingBox[6];
  this->boundingBox(tetsBoundingBox);
  VEC3F minTets;
  minTets[0] = tetsBoundingBox[0];
  minTets[1] = tetsBoundingBox[2];
  minTets[2] = tetsBoundingBox[4];
  VEC3F maxTets;
  maxTets[0] = tetsBoundingBox[1];
  maxTets[1] = tetsBoundingBox[3];
  maxTets[2] = tetsBoundingBox[5];
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
// compute the embedding of this surface into the tet mesh
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeEmbedding(OBJ* obj)
{
  vector<VEC3F>& vertices = obj->vertices;

  // compute a list of all the surface tetrahedra
  map<int, bool> surfaceTetIDsHash;
  for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
  {
    vector<int>& tetMembership = _tetMembership[_surfaceVertices[x]];
    for (unsigned int y = 0; y < tetMembership.size(); y++)
      surfaceTetIDsHash[tetMembership[y]] = true;
  }
  vector<int> surfaceTetIDs;
  map<int,bool>::iterator iter;
  for (iter = surfaceTetIDsHash.begin(); iter != surfaceTetIDsHash.end(); iter++)
    surfaceTetIDs.push_back(iter->first);

  // compute the center of mass of all the tetrahedra
  vector<VEC3F> centers;
  for (unsigned int x = 0; x < surfaceTetIDs.size(); x++)
    centers.push_back(_tets[surfaceTetIDs[x]].center());

  cout << " Computing nearest tets ... "; flush(cout);
  // for each point, compute the nearest tet
  _tetEmbeddings.clear();
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    Real minDistance = 0.0;
    int nearestTet = 0;
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
    TET& nearestTet = _tets[_tetEmbeddings[x]];
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
// compute area weighted normals of the surface vertices
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceNormals()
{
  if (_surfaceNormals.size() != _surfaceVertices.size())
    _surfaceNormals.resize(_surfaceVertices.size());

  // stomp old normals
  for (unsigned int x = 0; x < _surfaceNormals.size(); x++)
    _surfaceNormals[x].clear();

  // hash surface vertices (just the first time)
  if (_surfaceVertexID.size() == 0)
    for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
      _surfaceVertexID[_surfaceVertices[x]] = x;

  // compute normal for all faces the vertex is a member of --
  // opposing internal faces will cancel each other out
  for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
  {
    vector<int>& tetIDs = _tetMembership[_surfaceVertices[x]];

    for (unsigned int y = 0; y < tetIDs.size(); y++)
    {
      TET& tet = _tets[tetIDs[y]];

      // for each face of the tet
      for (int z = 0; z < 4; z++)
      {
        TRIANGLE face = tet.face(z);

        // see if the current vertex is a member of the face
        bool isMember = false;
        for (int a = 0; a < 3; a++)
          if (face.vertex(a) == _surfaceVertices[x])
            isMember = true;

        // if not, move on
        if (!isMember) continue;

        // add the normal to the current sum
        _surfaceNormals[x] += face.area() * face.normal();
      }
    }
  }
  /*
  // calculate area weighted normals for surface vertices
  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    pair<int, int> face = _surfaceFaces[x];
    TET& tet = _tets[face.first];
    TRIANGLE triangle = tet.face(face.second);

    VEC3F& normal = triangle.normal();
    Real area = 1.0;

    for (int y = 0; y < 3; y++)
    {
      int vertexID = _surfaceVertexID[triangle.vertex(y)];
      _surfaceNormals[vertexID] += area * normal;
    }
  }
  */

  // normalize the normals
  for (unsigned int x = 0; x < _surfaceNormals.size(); x++)
    if (norm(_surfaceNormals[x]) > 0.0)
      _surfaceNormals[x].normalize();
  /*
    else
    {
      cout << __FILE__ << " " << __LINE__ << " : " << endl;
      cout << "Surface vertex has no normal! " << endl;
    }
    */
}

//////////////////////////////////////////////////////////////////////
// Interpolate low res mesh normals to embedded mesh
//////////////////////////////////////////////////////////////////////
void TET_MESH::updateEmbeddingNormals()
{
  // make sure mesh normals are up to date
  computeSurfaceNormals();

  vector<VEC3>& normals = _embeddedMesh->normals;
  vector<VEC3>& vertices = _embeddedMesh->vertices;

  if (normals.size() != vertices.size())
    normals.resize(vertices.size());

  // stomp old normals
  for (unsigned int x = 0; x < normals.size(); x++)
    normals[x].clear();

  // for each embedded vertex
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    // get its tet and barycentric coords
    TET& tet = _tets[_tetEmbeddings[x]];
    VEC3F coords = _barycentricEmbeddings[x];
    Real weights[] = {coords[0], coords[1], coords[2], 0.0};
    weights[3] = 1.0 - coords[0] - coords[1] - coords[2];
  
    // for each vertex in the tet
    int totalNormals = 0;
    Real totalWeights = 0.0;
    for (int y = 0; y < 4; y++)
    {
      // see if it has a surface normal associated with it
      VEC3F* vertex = tet.vertices[y];
      if (_surfaceVertexID.find(vertex) != _surfaceVertexID.end())
      {
        int vertexID = _surfaceVertexID[vertex];
        //normals[x] += fabs(weights[y]) * _surfaceNormals[vertexID];
        normals[x] += weights[y] * _surfaceNormals[vertexID];

        totalNormals++;
        totalWeights += fabs(weights[y]);
      }
    } 
    if (totalNormals != 0)
      normals[x].normalize();
    else
    {
      cout << __FILE__ << " " << __LINE__ << " : " << endl;
      cout << " No normals found! " << endl;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// write the embedding to a file
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeEmbedding()
{
  string filename = _filename;
  filename += string(".embedding");
  FILE* file = fopen(filename.c_str(), "wb");

  int totalVertices = _tetEmbeddings.size();
  fwrite((void*)&totalVertices, sizeof(int), 1, file);

  // write out the tets
  for (int x = 0; x < totalVertices; x++)
    fwrite((void*)&_tetEmbeddings[x], sizeof(int), 1, file);

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
// try to read an embedding file
//////////////////////////////////////////////////////////////////////
bool TET_MESH::readEmbedding()
{
  string filename = _filename;
  filename += string(".embedding");
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
    int which;
    fread((void*)&which, sizeof(int), 1, file);
    _tetEmbeddings.push_back(which);
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
void TET_MESH::updateEmbedding()
{
  if (_embeddedMesh == NULL)
  {
    cout << " Embedded mesh has not been set!" << endl;
    return;
  }

  // recompute positions based on barycentric coords of embedding
  vector<VEC3>& vertices = _embeddedMesh->vertices;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    TET& tet = _tets[_tetEmbeddings[x]];
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
// update the position of the head embedding
//////////////////////////////////////////////////////////////////////
void TET_MESH::updateHeadEmbedding()
{
  if (_embeddedMesh == NULL)
  {
    cout << " Embedded mesh has not been set!" << endl;
    return;
  }

  // recompute positions based on barycentric coords of embedding
  vector<VEC3>& vertices = _embeddedMesh->vertices;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    /*
    // see if the vertex is in the right group
    vector<int>& groups = _embeddedMesh->groups();
    bool faceVertex = false;
    if (x >= groups[1] && x < groups[2]) faceVertex = true;
    */

    TET& tet = _tets[_tetEmbeddings[x]];
    /*
    VEC3F& v0 = (faceVertex) ? *(tet.vertices[0]) : _restPose[_vertexID[tet.vertices[0]]];
    VEC3F& v1 = (faceVertex) ? *(tet.vertices[1]) : _restPose[_vertexID[tet.vertices[1]]];
    VEC3F& v2 = (faceVertex) ? *(tet.vertices[2]) : _restPose[_vertexID[tet.vertices[2]]];
    VEC3F& v3 = (faceVertex) ? *(tet.vertices[3]) : _restPose[_vertexID[tet.vertices[3]]];
    */
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
// Output embedded triangles to RenderMan
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawEmbeddingToRenderMan()
{
#ifdef USING_RENDERMAN  
  //updateEmbeddingNormals();
  _embeddedMesh->ComputeVertexNormals();
  for (int x = 0; x < 10; x++)
    _embeddedMesh->SmoothVertexNormals();

  vector<VEC3F>& vertices = _embeddedMesh->vertices;
  vector<VEC3F>& normals = _embeddedMesh->normals;
  //vector<VEC2>& texcoords = _embeddedMesh->texcoords;
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
// Output embedded triangles to RenderMan
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawHeadEmbeddingToRenderMan(bool bakingSurface, bool subsurface, string filename)
{
#ifdef USING_RENDERMAN
  _embeddedMesh->ComputeVertexNormals();
  //cout << __FILE__ << " " << __LINE__ << " : " << endl;
  //cout << " NORMAL SMOOTHING DISABLED." << endl;
  for (int x = 0; x < 5; x++)
    _embeddedMesh->SmoothVertexNormals();

  vector<VEC3F>& vertices = _embeddedMesh->vertices;
  vector<VEC3F>& normals = _embeddedMesh->normals;
  vector<VEC2>& texcoords = _embeddedMesh->texcoords;
  vector<OBJ::Face>& faces = _embeddedMesh->faces;
  vector<int>& vertexGroups = _embeddedMesh->vertexGroups();
  vector<int>& faceGroups = _embeddedMesh->faceGroups();
  //vector<int>& texcoordGroups = _embeddedMesh->texcoordGroups();
  //vector<string>& groupNames = _embeddedMesh->groupNames();

  /*
  cout << " total vertices: " << vertices.size() << endl;
  cout << " total texcoords: " << texcoords.size() << endl;
  cout << " total faces: " << faces.size() << endl << endl;
  cout << " total vertex groups: " << vertexGroups.size() << endl;
  cout << " total face groups: " << faceGroups.size() << endl;
  cout << " total texcoord groups: " << texcoordGroups.size() << endl;
  */

  RtColor colors[5];
  colors[0][0] = 1.0; colors[0][1] = 1.0; colors[0][2] = 1.0;
  colors[1][0] = 1.0; colors[1][1] = 0.0; colors[1][2] = 0.0;
  colors[2][0] = 0.0; colors[2][1] = 1.0; colors[2][2] = 0.0;
  colors[3][0] = 0.0; colors[3][1] = 0.0; colors[3][2] = 1.0;
  colors[4][0] = 1.0; colors[4][1] = 0.0; colors[4][2] = 1.0;

  RtString specmapname = "/hydra_nfs/tedkim/svnrepos/head_model/specMap.tex";

  RtString textureFiles[] = {
    "",
    "/hydra_nfs/tedkim/svnrepos/head_model/eye.tex",
    "",
    "/hydra_nfs/tedkim/svnrepos/head_model/eye.tex",
    "/hydra_nfs/tedkim/svnrepos/head_model/colorMap.tex"
  };

  RtColor surfaceColors[] = {
    {0.552000, 0.320000, 0.320000}, // mouth
    //{0,0,0},                        // eye
    {1,1,1},                        // eye
    {0.400000, 0.384000, 0.312000}, // teeth
    //{0,0,0},                        // eye
    //{0,0,0}                         // head
    {1,1,1},                        // eye
    {1,1,1}                         // head
  };

  RtColor specularColors[] = {
    {0.360000, 0.360000, 0.360000}, // mouth
    {0.250000, 0.250000, 0.250000}, // eye
    {0.250000, 0.230000, 0.200000}, // teeth
    {0.250000, 0.250000, 0.250000}, // eye
    //{0,0,0}                         // head
    {1,1,1}                         // head
  };

  RtFloat roughness[] = {
    96.078431,
    96.078431,
    15.686275,
    96.078431,
    96.078431
  };

  RtFloat eccentricity[] = {
    0.074, // mouth
    0.074, // eye
    0.074, // teeth
    0.074, // eye
    0.231 // head
  };

  RtFloat rolloff[] = {
    0.7, // mouth
    0.7, // eye
    0.7, // teeth
    0.7, // eye
    0.719 // head
  };

  RtFloat reflectivity[] = {
    0.5, // mouth
    0.116, // eye
    0, // teeth
    0.116, // eye
    0// head
  };

  for (unsigned int x = 0; x < faceGroups.size() - 1; x++)
  {
    if (!bakingSurface && !subsurface)
    {
      //RtString textureFile = "/hydra_nfs/tedkim/svnrepos/head_model/colorMap.tex"; 
#define TEXTURED 1
#if TEXTURED
      RiColor(surfaceColors[x]);

      // setup the surface shader
      RiDeclare("texname", "uniform string");
      RiDeclare("specularcolor", "uniform color");
      RiDeclare("roughness", "uniform float");
      //RiSurface("faceshader", 
      RiSurface("./shaders/faceshader", 
                (RtToken)"texname", (RtPointer)&textureFiles[x], 
                (RtToken)"specularcolor", (RtPointer)&specularColors[x],
                (RtToken)"roughness", (RtPointer)&roughness[x],
                RI_NULL);
#else
      cout << __FILE__ << " " << __LINE__ << " : " << endl;
      cout << " Face shader disabled!" << endl;
#endif
    }
    if (bakingSurface)
    {
      RiColor(surfaceColors[x]);
      //RtColor white = {10,10,10};
      //RiColor(white);
    }

    if (subsurface)
    {
      RiColor(surfaceColors[x]);
      char buffer[256];
      sprintf(buffer, "%s.sss.ptc", filename.c_str());
      RtString sssFile = buffer;
      RiDeclare("unitlength", "uniform float");
      RtFloat unitlength = 0.01;

      /*
      RiSurface("render_ssdiffusion", 
                (RtToken)"filename", (RtPointer)&sssFile, 
                (RtToken)"unitlength", (RtPointer)&unitlength,
                RI_NULL);
                */
      if (x == 0 || x == 4)
      { 
        RiDeclare("texname", "uniform string");
        RiDeclare("specmapname", "uniform string");
        RiDeclare("ptcname", "uniform string");
        RiDeclare("specularcolor", "uniform color");
        RiDeclare("roughness", "uniform float");
        RiDeclare("eccentricity", "uniform float");
        RiDeclare("specularRollOff", "uniform float");
        RiDeclare("reflectivity", "uniform float");
        RiSurface("./shaders/faceshader", 
                  (RtToken)"texname", (RtPointer)&textureFiles[x], 
                  (RtToken)"specmapname", (RtPointer)&specmapname, 
                  (RtToken)"ptcname", (RtPointer)&sssFile,
                  (RtToken)"specularcolor", (RtPointer)&specularColors[x],
                  //(RtToken)"roughness", (RtPointer)&roughness[x],
                  (RtToken)"unitlength", (RtPointer)&unitlength,
                  (RtToken)"eccentricity", (RtPointer)&eccentricity[x],
                  (RtToken)"specularRollOff", (RtPointer)&rolloff[x],
                  (RtToken)"reflectivity", (RtPointer)&reflectivity[x],
                  RI_NULL);
      }
      else
      {
        RiDeclare("texname", "uniform string");
        RiDeclare("ptcname", "uniform string");
        RiDeclare("reflectivity", "uniform float");
        RiDeclare("diffuseCoeff", "uniform float");
        RiDeclare("cosPower", "uniform float");

        RtFloat diffuseCoeff = 0.8;
        RtFloat cosPower = 15;
        RtFloat reflectivity = 1;
        diffuseCoeff = (x == 2) ? 0.901 : 1.0;
        cosPower = (x == 2) ? 20.0 : 15.0;
        reflectivity = (x == 2) ? 0.5 : 0.116; 
        RiSurface("./shaders/phongshader",
                  (RtToken)"texname", (RtPointer)&textureFiles[x], 
                  (RtToken)"ptcname", (RtPointer)&sssFile,
                  (RtToken)"unitlength", (RtPointer)&unitlength,
                  (RtToken)"reflectivity", (RtPointer)&reflectivity,
                  (RtToken)"diffuseCoeff", (RtPointer)&diffuseCoeff,
                  (RtToken)"cosPower", (RtPointer)&cosPower,
                  RI_NULL);
      }
    }

    int faceStart = faceGroups[x];
    int vertexStart = vertexGroups[x];
    //int texcoordStart = texcoordGroups[x];
    int totalFaces = faceGroups[x+1] - faceGroups[x];
    int totalVertices = vertexGroups[x+1] - vertexGroups[x];
    //int totalTexcoords = texcoordGroups[x+1] - texcoordGroups[x];

    // init point locations
    RtPoint* P = new RtPoint[totalVertices];
    for (int y = 0; y < totalVertices; y++)
    {
      P[y][0] = vertices[vertexStart + y][0];
      P[y][1] = vertices[vertexStart + y][1];
      P[y][2] = vertices[vertexStart + y][2];
    }

    RtPoint* N = new RtPoint[totalVertices];
    //for (int y = 0; y < normals.size(); y++)
    for (int y = 0; y < totalVertices; y++)
    {
      N[y][0] = normals[vertexStart + y][0];
      N[y][1] = normals[vertexStart + y][1];
      N[y][2] = normals[vertexStart + y][2];
    }

    // init vertex surface colors
    RtColor* vertexColors = new RtColor[totalVertices];
    for (int y = 0; y < totalVertices; y++)
    {
      vertexColors[y][0] = surfaceColors[x][0];
      vertexColors[y][1] = surfaceColors[x][1];
      vertexColors[y][2] = surfaceColors[x][2];
    }

    // init all to triangles
    RtInt* nvertices = new RtInt[totalFaces];
    for (int y = 0; y < totalFaces; y++)
      nvertices[y] = 3;

    // init faces
    RtInt* faceIndices= new RtInt[3 * totalFaces];
    RtFloat* s = new RtFloat[3 * totalFaces];
    RtFloat* t = new RtFloat[3 * totalFaces];
    for (int y = 0; y < totalFaces; y++)
    {
      OBJ::Face face = faces[faceStart + y];
      faceIndices[y * 3]     = face.vertices[0] - vertexStart;
      faceIndices[y * 3 + 1] = face.vertices[1] - vertexStart;
      faceIndices[y * 3 + 2] = face.vertices[2] - vertexStart;

      s[y * 3]     = texcoords[face.texcoords[0]][0];
      s[y * 3 + 1] = texcoords[face.texcoords[1]][0];
      s[y * 3 + 2] = texcoords[face.texcoords[2]][0];

      t[y * 3]     = 1.0 - texcoords[face.texcoords[0]][1];
      t[y * 3 + 1] = 1.0 - texcoords[face.texcoords[1]][1];
      t[y * 3 + 2] = 1.0 - texcoords[face.texcoords[2]][1];
    }
    RtToken svarying = "facevarying float s";
    RtToken tvarying = "facevarying float t";

    RiPointsPolygons(totalFaces, nvertices, faceIndices, 
                     RI_P, P, 
                     RI_N, N, 
                     svarying, s, 
                     tvarying, t, 
                     //RI_CS, vertexColors,
                     RI_NULL);
    delete[] faceIndices;
    delete[] nvertices;
    delete[] P;
    delete[] s;
    delete[] t;
    delete[] vertexColors;
    delete[] N;
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// compute the Hessian product with respect to two vectors --
// H * first * second
//////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::hessianProduct(VECTOR& first, VECTOR& second)
{
  VECTOR result(this->TET_MESH::rank());
  result.clear();

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    int materialIndex = _tets[x].materialIndex();
    
    // get the Hessian
    TENSOR3 hessian = _materials[materialIndex]->hessian(_tets[x]);

    //if (hessian.sum2() > 1e-4)
    //  cout << " Hessian: " << hessian << endl;

    // gather the vertex indices for this tet
    int indices[4];
    for (int y = 0; y < 4; y++)
      indices[y] = _vertexID[_tets[x].vertices[y]];

    // gather the vector info for this tet
    VECTOR subFirst(12);
    VECTOR subSecond(12);
    for (int y = 0; y < 4; y++)
      // if it's not a constrained node
      if (indices[y] < _unconstrainedSize)
      {
        int index3 = indices[y] * 3;
        subFirst[3 * y]     = first[index3];
        subFirst[3 * y + 1] = first[index3 + 1];
        subFirst[3 * y + 2] = first[index3 + 2];
        
        subSecond[3 * y] = second[index3];
        subSecond[3 * y + 1] = second[index3 + 1];
        subSecond[3 * y + 2] = second[index3 + 2];
      }

    // do the product
    MATRIX firstProduct = hessian.modeThreeProduct(subFirst);
    VECTOR secondProduct = firstProduct * subSecond;

    // scatter the vector info into the final result
    for (int y = 0; y < 4; y++)
      // if it's not a constrained node
      if (indices[y] < _unconstrainedSize)
      {
        int index3 = indices[y] * 3;
        result[index3] += secondProduct[3 * y];
        result[index3 + 1] += secondProduct[3 * y + 1];
        result[index3 + 2] += secondProduct[3 * y + 2];
      }
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// get the normal to this vertex
//////////////////////////////////////////////////////////////////////
VEC3F TET_MESH::surfaceNormal(VEC3F* vertex)
{
  if (_surfaceNormals.size() == 0)
    computeSurfaceNormals();

  VEC3F normal;
  if (_vertexOnSurface.find(vertex) == _vertexOnSurface.end())
    return normal;

  int id = _surfaceVertexID[vertex];
  return _surfaceNormals[id];
}

//////////////////////////////////////////////////////////////////////
// Get the mean stretch ratio
//////////////////////////////////////////////////////////////////////
Real TET_MESH::meanStretchRatio()
{
  Real accum = 0.0;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    MATRIX3 F3 = _tets[x].F();
    MATRIX F(F3);
    VECTOR eigenvalues(3);
    MATRIX eigenvectors(3,3);
#ifdef USING_MKL
    F.eigensystem(eigenvalues, eigenvectors);
#else
    F.eigensystem3x3(eigenvalues, eigenvectors);
#endif

    //Real incompressibility = eigenvalues[0] * eigenvalues[1] * eigenvalues[2];
    //cout << "eigenvalues: " << eigenvalues << endl;
    //cout << "Incompressbility: " << incompressibility << endl;

    //accum += eigenvalues.maxValue();
    accum += eigenvalues.minValue();
  }
  return accum / _tets.size();
}

//////////////////////////////////////////////////////////////////////
// min tet volume
//////////////////////////////////////////////////////////////////////
Real TET_MESH::minTetVolume()
{
  Real minFound = _tets[0].volume();

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    if (_tets[x].volume() < minFound)
      minFound = _tets[x].volume();
  }
  return minFound;
}

//////////////////////////////////////////////////////////////////////
// max tet volume
//////////////////////////////////////////////////////////////////////
Real TET_MESH::maxTetVolume()
{
  Real maxFound = _tets[0].volume();

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    if (_tets[x].volume() > maxFound)
      maxFound = _tets[x].volume();
  }
  return maxFound;
}

//////////////////////////////////////////////////////////////////////
// Total fully connected nodes
//////////////////////////////////////////////////////////////////////
int TET_MESH::testFullyConnected()
{
  map<VEC3F*, bool> verticesCrawled;
  map<int, bool> tetsCrawled;

  // seed it with the last vertex, which should be constrained if
  // there are any
  vector<VEC3F*> connections;
  oneRing(&_vertices[_vertices.size() - 1], connections);
  vector<VEC3F*> toCrawl;

  // prime the stack
  for (unsigned int x = 0; x < connections.size(); x++)
    toCrawl.push_back(connections[x]);

  while (toCrawl.size() > 0)
  {
    VEC3F* connected = toCrawl.back();
    toCrawl.pop_back();
    verticesCrawled[connected] = true;

    // tag the tets this vertex is a member of as crawled as well
    vector<int>& tetMembers = _tetMembership[connected];
    for (unsigned int x = 0; x < tetMembers.size(); x++)
      tetsCrawled[tetMembers[x]] = true;

    oneRing(connected, connections);
    for (unsigned int x = 0; x < connections.size(); x++)
      if (verticesCrawled.find(connections[x]) == verticesCrawled.end())
        toCrawl.push_back(connections[x]);
  }

  cout << verticesCrawled.size() << " vertices crawled of " << _vertices.size() << endl;
  cout << tetsCrawled.size() << " tets crawled of " << _tets.size() << endl;

  // check the dihedrals
  Real minDihedral = _tets[0].minDihedral();
  Real maxDihedral = _tets[0].maxDihedral();
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    if (_tets[x].minDihedral() < minDihedral)
      minDihedral = _tets[x].minDihedral();
    if (_tets[x].maxDihedral() > maxDihedral)
      maxDihedral = _tets[x].maxDihedral();
  }
  cout << " Min dihedral: " << minDihedral / (2.0 * M_PI) * 360.0 << " Max dihedral: " << maxDihedral / (2.0 * M_PI) * 360.0 << endl;

  return verticesCrawled.size();
}

//////////////////////////////////////////////////////////////////////
// get the face-based one ring of a tet
//////////////////////////////////////////////////////////////////////
void TET_MESH::tetOneRing(int tetIndex, vector<int>& oneRing)
{
  TET& tet = _tets[tetIndex];

  // get the one rings of its vertices
  map<int, bool> allTets;
  for (int x = 0; x < 4; x++)
  {
    VEC3F* vertex = tet.vertices[x];
    vector<int>& membership = _tetMembership[vertex];

    for (unsigned int y = 0; y < membership.size(); y++)
      if (membership[y] != tetIndex)
        allTets[membership[y]] = true;
  }

  // see which tets share a face
  map<int, bool>::iterator i;
  for (i = allTets.begin(); i != allTets.end(); i++)
  {
    // get the index of the tet to compare to
    int otherIndex = i->first;
    TET& otherTet = _tets[otherIndex];

    // if they do in fact share a face, store it
    int sharedFace = tet.sharedFace(otherTet);
    if (sharedFace >= 0)
    {
      TRIANGLE face = tet.face(sharedFace);
      // removed shared face membership recording here
      // (see getTetOneRing in ISO_STUFFER)
      oneRing.push_back(otherIndex);
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Recompute x based on current deformation
//////////////////////////////////////////////////////////////////////
void TET_MESH::recoverX()
{
  _x.resizeAndWipe(_unconstrainedSize * 3);
  
  // total surface vertices
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F& restPose = _restPose[x];
    VEC3F& deformed = _vertices[x];
    int index = 3 * x;

    _x(index)     = deformed[0] - restPose[0];
    _x(index + 1) = deformed[1] - restPose[1];
    _x(index + 2) = deformed[2] - restPose[2];
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
bool TET_MESH::centerOfMassIsZero()
{
  // unambiguously call the TET_MESH version
  TET_MESH::updateFullMesh();
 
  // is the center of mass for _vertices being used?
  VEC3F vertexSum;
  int size = _vertices.size();
  Real totalMass = 0;
  for (int x = 0; x < size; x++)
  {
    vertexSum += mass(x) * _vertices[x];
    totalMass += mass(x);
  }
  vertexSum *= 1.0 / totalMass;

  // if this is tripped, then the center of mass aren't centered at zero
  if (vertexSum * vertexSum > 1e-4)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Center of mass test failed! " << endl;
    cout << " vertexSum 2 norm: " << vertexSum * vertexSum << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////
// Return all the rest vertices stacked into a vector
//////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::restVector()
{
  VECTOR final(3 * _unconstrainedSize);

  for (int x = 0; x < _unconstrainedSize; x++)
  {
    final[3 * x] = _restPose[x][0];
    final[3 * x + 1] = _restPose[x][1];
    final[3 * x + 2] = _restPose[x][2];
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a block identity matrix with the same rank as this mesh
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_MESH::sparseIdentity()
{
  SPARSE_MATRIX final(3 * _unconstrainedSize, 3);

  for (int x = 0; x < _unconstrainedSize; x++)
  {
    final(3 * x, 0) = 1.0;
    final(3 * x + 1, 1) = 1.0;
    final(3 * x + 2, 2) = 1.0;
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a block identity matrix with the same rank as this mesh
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX TET_MESH::blockIdentity()
{
  BLOCK_MATRIX final(_unconstrainedSize, 1);

  MATRIX eye(3,3);
  eye(0,0) = 1.0;
  eye(1,1) = 1.0;
  eye(2,2) = 1.0;

  for (int x = 0; x < _unconstrainedSize; x++)
    final.add(eye, x, 0);

  return final;
}

//////////////////////////////////////////////////////////////////////
// refresh the inertia tensor
//////////////////////////////////////////////////////////////////////
const MATRIX& TET_MESH::refreshInertiaTensor()
{
  _inertiaTensor.clear();
  vector<VEC3F>& vertices = _vertices;
  //for (unsigned int x = 0; x < vertices.size(); x++)
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    Real mass = this->mass(x);
    Real xSq = vertices[x][0] * vertices[x][0];
    Real ySq = vertices[x][1] * vertices[x][1];
    Real zSq = vertices[x][2] * vertices[x][2];

    MATRIX entry(3,3);

    entry(0,0) = (ySq + zSq);
    entry(0,1) = -vertices[x][1] * vertices[x][0];
    entry(0,2) = -vertices[x][2] * vertices[x][0];
    entry(1,0) = entry(0,1);
    entry(1,1) = (xSq + zSq);
    entry(1,2) = -vertices[x][2] * vertices[x][1];
    entry(2,0) = entry(0,2);
    entry(2,1) = entry(1,2);
    entry(2,2) = (xSq + ySq);
    entry *= mass;
    _inertiaTensor += entry;
    //cout << " mass: " << mass << endl;
    //cout << " vertex: " << vertices[x] << endl;
    //cout << " current tensor: " << _inertiaTensor << endl;
  }

  // center of mass when computing inertia tensor
  VEC3F center;
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    Real mass = this->mass(x);
    center += mass * _vertices[x];
  }
  center *= 1.0 / _totalMass;

  return _inertiaTensor;
}

//////////////////////////////////////////////////////////////////////
// refresh the inertia tensor time derivative
//////////////////////////////////////////////////////////////////////
const MATRIX& TET_MESH::refreshInertiaTensorDt(VECTOR& velocity)
{
  _inertiaTensorDt.clear();
  vector<VEC3F>& vertices = _vertices;
  for (int x = 0; x < velocity.size() / 3; x++)
  {
    Real mass = this->mass(x);
    VEC3F uDot;
    uDot[0] = velocity[3 * x];
    uDot[1] = velocity[3 * x + 1];
    uDot[2] = velocity[3 * x + 2];

    VEC3F& u = vertices[x];

    MATRIX entry(3,3);

    entry(0,0) = (u[1] * uDot[1] + u[2] * uDot[2]);
    entry(0,0) = 2.0 * entry(0,0);

    entry(0,1) = -(uDot[1] * u[0] + u[1] * uDot[0]);
    entry(0,2) = -(uDot[2] * u[0] + u[2] * uDot[0]);

    entry(1,1) = (u[0] * uDot[0] + u[2] * uDot[2]);
    entry(1,1) = 2.0 * entry(1,1);

    entry(1,2) = -(uDot[2] * u[1] + u[2] * uDot[1]);
    
    entry(2,2) = (u[0] * uDot[0] + u[1] * uDot[1]);
    entry(2,2) = 2.0 * entry(2,2);

    entry(1,0) = entry(0,1);
    entry(2,0) = entry(0,2);
    entry(2,1) = entry(1,2);

    entry *= mass;
    _inertiaTensorDt += entry;
  }

  return _inertiaTensorDt;
}

//////////////////////////////////////////////////////////////////////
// recompute current center of mass
//////////////////////////////////////////////////////////////////////
VEC3F TET_MESH::SitBar()
{
  updateFullMesh();
  
  VEC3F centerOfMass;
  Real totalMass = 0.0;
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    centerOfMass += mass(x) * _vertices[x];
    totalMass += mass(x);
  }
  centerOfMass *= 1.0 / totalMass;

  return centerOfMass;
}

//////////////////////////////////////////////////////////////////////
// recompute current center of mass, minus rest pose
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_MESH::SiBar()
{
  SPARSE_MATRIX final(3, this->rank());

  for (int x = 0; x < _unconstrainedSize; x++)
  {
    final(0, 3 * x) = mass(x);
    final(1, 3 * x + 1) = mass(x);
    final(2, 3 * x + 2) = mass(x);
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// draw vector field on top of the existing mesh
//////////////////////////////////////////////////////////////////////
void TET_MESH::drawVectorField(VECTOR field)
{
  assert(field.size() == _unconstrainedSize * 3);

  // factor the scale the field down by
  Real scale = 0.1;
 
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F& vertex = _vertices[x];

    VEC3F endPoint;
    endPoint[0] = field[3 * x];
    endPoint[1] = field[3 * x + 1];
    endPoint[2] = field[3 * x + 2];
    endPoint *= scale;

    glLineWidth(1);
    glPushMatrix();
      glTranslatef(vertex[0], vertex[1], vertex[2]);
      glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(endPoint[0], endPoint[1], endPoint[2]);
      glEnd();
    glPopMatrix();
  }
}

//////////////////////////////////////////////////////////////////////
// compute the shape matching transformation
// make sure to call updateFullMesh() before this!
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeShapeMatching(VEC3F& translation, MATRIX3& rotation) const
{
  // recompute the center of mass, just in case
  //computeCenterOfMass();
  //computeRestCenterOfMass();

  // compute centers of mass from scratch, but don't store it in order
  // to enforce const
  /*
  VEC3F centerOfMass;
  VEC3F restCenterOfMass;
  Real totalMass = 0.0;
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    Real mass = _masses.constEntry(3 * x, 3 * x);
    centerOfMass += mass * _vertices[x];
    restCenterOfMass += mass * _restPose[x];
    totalMass += mass;
  }
  centerOfMass *= 1.0 / totalMass;
  restCenterOfMass *= 1.0 / totalMass;
  */

  // find geometric center, not center of mass, as uneven mass distribution
  // can just introduce rigid translation modes
  VEC3F geometricCenter;
  VEC3F geometricRestCenter;
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    geometricCenter += _vertices[x];
    geometricRestCenter += _restPose[x];
  }
  geometricCenter *= 1.0 / _unconstrainedSize;
  geometricRestCenter *= 1.0 / _unconstrainedSize;

  // compute the matching matrix
  MATRIX3 Apq;
  for (unsigned int x = 0; x < _unconstrainedSize; x++)
  {
    // center the vertex
    VEC3F vertex = _vertices[x];
    vertex -= geometricCenter;

    // center the rest pose
    VEC3F rest = _restPose[x];
    rest -= geometricRestCenter;

    Real mass = _masses.constEntry(3 * x, 3 * x);
    Real vertexMass = ((int)x < _unconstrainedSize) ?  mass : 0;
    Apq += vertexMass * MATRIX3::outer_product(vertex, rest);
  }

  // diagonalize
  MATRIX fullApq(Apq);
  VECTOR fullEigenvalues(3);
  MATRIX fullEigenvectors(3,3);
  MATRIX ApqTApq = fullApq.transpose() * fullApq;

  ApqTApq.eigensystem(fullEigenvalues, fullEigenvectors);

  // get the inverse square root
  for (int x = 0; x < 3; x++)
    fullEigenvalues[x] = 1.0 / sqrt(fullEigenvalues[x]);

  MATRIX3 diag;
  diag(0,0) = fullEigenvalues[0];
  diag(1,1) = fullEigenvalues[1];
  diag(2,2) = fullEigenvalues[2];
 
  MATRIX3 eigenvectors;
  MATRIX3 eigenvectorsInv;
  fullEigenvectors.copiesInto(eigenvectors);

  MATRIX3 Sinv = eigenvectors * diag * eigenvectors.inverse();

  rotation = Apq * Sinv;

  // for the translation being returned, compute it in terms of the original
  // rest center of mass
  translation = geometricCenter;
}

//////////////////////////////////////////////////////////////////////
// compute the deformation minus any rigid shape matching component
//////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::computeDefoMinusShapeMatching() const
{
  VEC3F translation;
  MATRIX3 rotation;

  // compute the shape matching
  computeShapeMatching(translation, rotation);

  // be careful to use the geometric center, not the center of mass.
  // All we care about it the kinematics, not the mass distribution
  VEC3F geometricRestCenter;
  for (int x = 0; x < _unconstrainedSize; x++)
    geometricRestCenter += _restPose[x];
  geometricRestCenter *= 1.0 / _unconstrainedSize;

  VECTOR defoOnly(_unconstrainedSize * 3);

  // undo the shape matching transform
  MATRIX3 RT = rotation.transpose();
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F vertex = _vertices[x];
    vertex -= translation;
    vertex = RT * vertex;

    // center the rest pose
    VEC3F rest = _restPose[x];
    rest -= geometricRestCenter;

    VEC3F defo = vertex - rest;

    defoOnly[3 * x] = defo[0];
    defoOnly[3 * x + 1] = defo[1];
    defoOnly[3 * x + 2] = defo[2];
  }

  return defoOnly;
}

//////////////////////////////////////////////////////////////////////
// compute the angular vector induced by a vector field
//////////////////////////////////////////////////////////////////////
VEC3F TET_MESH::computeAngularVector(const VECTOR& field, const VEC3F& center, bool verbose)
{
  //VEC3F angularCenter = computeAngularCenter();

  // build the cross product matrix
  vector<MATRIX> matrices;
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F vertex = _vertices[x];
    //vertex -= angularCenter;
    vertex -= center;

    /*
    if (x == 0 && verbose)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Original vertex position: " << _vertices[x] << endl;
      cout << " Vertex position used: " << vertex << endl;
      cout << " translation used: " << center << endl;
    }
    */

    matrices.push_back(MATRIX::cross(vertex));
  }

  // build the stacked matrix
  MATRIX stacked = MATRIX::columnOfMatrices(matrices);

  // solve the least squares problem
  VECTOR fieldCopy = field;
  fieldCopy *= -1;

  stacked.solveLeastSquares(fieldCopy);

  VEC3F foundAngular;
  foundAngular[0] = fieldCopy[0];
  foundAngular[1] = fieldCopy[1];
  foundAngular[2] = fieldCopy[2];

  return foundAngular;
}

//////////////////////////////////////////////////////////////////////
// compute the angular vector induced by a vector field
//////////////////////////////////////////////////////////////////////
VECTOR TET_MESH::computeAngularField(const VEC3F& angular, const VEC3F& center)
{
  //VEC3F angularCenter = computeAngularCenter();

  VECTOR field(_unconstrainedSize * 3);

  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F vertex = _vertices[x];
    //vertex -= angularCenter;
    vertex -= center;

    VEC3F angularCross = cross(angular, vertex);
    field[3 * x] = angularCross[0];
    field[3 * x + 1] = angularCross[1];
    field[3 * x + 2] = angularCross[2];
  }

  return field;
}

//////////////////////////////////////////////////////////////////////
// subtract the rigid components from a field
//////////////////////////////////////////////////////////////////////
void TET_MESH::subtractRigidComponent(const VECTOR& position, const VEC3F& rigidCenter, 
                                      VECTOR& field, VEC3F& linearVector, VEC3F& angularVector)
{
  VECTOR backupX = _x;

  _x = position;
  TET_MESH::updateFullMesh();

  // compute and subtract angular component
  //  order is important here! If the translation is subtracted off
  //  first, the field will obtain a new translation component after
  //  the angular component is factored off
  angularVector = computeAngularVector(field, rigidCenter, true);
  field -= computeAngularField(angularVector, rigidCenter);

  // compute and subtract translation component
  linearVector *= 0;
  for (int x = 0; x < field.size() / 3; x++)
  {
    linearVector[0] += field[3 * x];
    linearVector[1] += field[3 * x + 1];
    linearVector[2] += field[3 * x + 2];
  }
  linearVector *= 1.0 / (field.size() / 3.0);
  
  for (int x = 0; x < field.size() / 3; x++)
  {
    field[3 * x] -= linearVector[0];
    field[3 * x + 1] -= linearVector[1];
    field[3 * x + 2] -= linearVector[2];
  }

  _x = backupX;
  TET_MESH::updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
// compute geometric center (i.e. non mass-weighted):
//////////////////////////////////////////////////////////////////////
VEC3F TET_MESH::computeGeometricCenter() const
{
  VEC3F geometricCenter;
  for (int x = 0; x < _unconstrainedSize; x++)
    geometricCenter += _vertices[x];
  geometricCenter *= 1.0 / _unconstrainedSize;

  return geometricCenter;
}

//////////////////////////////////////////////////////////////////////
// compute geometric center (i.e. non mass-weighted):
//////////////////////////////////////////////////////////////////////
VEC3F TET_MESH::computeGeometricRestCenter() const
{
  VEC3F geometricRestCenter;
  for (int x = 0; x < _unconstrainedSize; x++)
    geometricRestCenter += _restPose[x];
  geometricRestCenter *= 1.0 / _unconstrainedSize;

  return geometricRestCenter;
}

//////////////////////////////////////////////////////////////////////
// compute the center that is being used to compute any angular fields
// this is only here so that we can switch easily between the
// geometric and mass centers
//////////////////////////////////////////////////////////////////////
VEC3F TET_MESH::computeAngularCenter()
{
  //computeCenterOfMass();
  //return _centerOfMass;

  //return computeGeometricCenter();
  
  return computeGeometricRestCenter();

  /*
  VEC3F translation;
  MATRIX3 rotation;
  computeShapeMatching(translation, rotation);

  translation -= rotation * computeGeometricRestCenter();

  return translation;
  */

  //computeRestCenterOfMass();
  //return _restCenterOfMass;
}
