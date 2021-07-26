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
// SUBSPACE_TET_MESH.h: interface for the SUBSPACE_TET_MESH class.
//
//////////////////////////////////////////////////////////////////////

#include "PARTITIONED_SUBSPACE_TET_MESH.h"

#if !USING_OSX
#include <GL/glut.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
PARTITIONED_SUBSPACE_TET_MESH::PARTITIONED_SUBSPACE_TET_MESH(const char* filename, MATERIAL** materials, int totalMaterials, int partitions, Real springConst, bool simulate, bool loadOriginal, string partitionPath) :
  _interfaceU(NULL),
  _restDiffs(NULL),
  _restSandwich(NULL),
  _basisSandwich(NULL),
  _translationSandwich(NULL),
  _constraintRestSandwich(NULL),
  _constraintBasisSandwich(NULL),
  _meanRestInterface(NULL),
  _meanBasisInterface(NULL)
{
  _partitions = partitions;
  _filename = string(filename);
  _partitionPath = partitionPath;
  
  // allocate per-interface spring consts
  _interfaceSpringConst.resizeAndWipe(_partitions, _partitions);
  _springConst = springConst;
  for (int y = 0; y < _partitions; y++)
    for (int x = 0; x < _partitions; x++)
      _interfaceSpringConst(x,y) = _springConst;

  // read in the original mesh
  if (loadOriginal)
    _originalMesh = new SUBSPACE_TET_MESH(filename, materials, totalMaterials);
  else
    _originalMesh = new SUBSPACE_TET_MESH(filename, materials, totalMaterials, false, "dontReadAnything", "dontReadAnything");

  // read in each separate partition
  _meshes = new TET_MESH*[_partitions];
  _ranks = new int[_partitions];
  _superRank = 0;
  _superRankPlusRigids = 0;
  for (int x = 0; x < _partitions; x++)
  {
    cout << " ==============================" << endl;
    cout << " Loading partition " << x << endl;
    cout << " ==============================" << endl;
    char buffer[256];
    sprintf(buffer, "%i", x);
    string partitionFilename = partitionPath;
    if (partitionFilename.length() == 0)
      partitionFilename = string(filename);
    partitionFilename += string(".partition.");
    partitionFilename += string(buffer);

    // if we are simulating, attach frames. Otherwise, we are training,
    // and extra translation modes should *not* be added to the basis
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
      bool simulateFullspace = true;
      if (simulate)
        simulateFullspace = false;
      _meshes[x] = new UNCONSTRAINED_SUBSPACE_TET_MESH(partitionFilename.c_str(), materials, totalMaterials, simulateFullspace);
      _unconstrainedPartition.push_back(true);
      if (simulate)
        ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x])->cacheSubspaceCenterOfMass();
    }
    else
    {
      _meshes[x] = new SUBSPACE_TET_MESH(partitionFilename.c_str(), materials, totalMaterials);
      _unconstrainedPartition.push_back(false);
    }
    cout << endl;

    int rank = ((SUBSPACE_TET_MESH*)_meshes[x])->rank();
    cout << " Submesh rank: " << rank << endl;
    _ranks[x] = rank;
    _superRank += rank;
    _superRankPlusRigids += rank;
    if (constrainedSize == 0)
      _superRankPlusRigids += 6;
  }

  cout << " ==============================" << endl;
  cout << " Done loading all partitions." << endl;
  cout << " Super rank: " << _superRank << endl;
  cout << " Super rank with rigids: " << _superRankPlusRigids << endl;
  cout << " ==============================" << endl;

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

        // only store the clone if it is unconstrained, otherwise
        // it will just screw up the indexing of the interface nodes
        if (first < _meshes[x]->unconstrainedNodes())
          _clonedVertices[x][y].push_back(clone);

        VEC3F* leftClone = _meshes[x]->vertices(first);
        VEC3F* rightClone = _meshes[y]->vertices(second);
        _isCloned[leftClone]++;
        _isCloned[rightClone]++;
      }
    }
  }
  fclose(file);

  // recompute the mass matrix based on the clones
  recomputePartitionedMass();

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
  cout << " Inverse correspondence size: " << _partitionIDs.size() << endl;
  cout << " Tet mesh size: " << _originalMesh->totalNodes() << endl;

  // compute the center of mass so we can draw the exploded view
  computeCenterOfMass();

  // compute the cloned triangles
  cout << " Computing cloned triangles ..."; flush(cout);
  if (!readClonedTriangles())
  {
    cout << " no cache found! ... "; flush(cout);
    computeClonedTriangles();
    writeClonedTriangles();
  }
  else
    cout << " cache found! ... "; flush(cout);
  computeClonedSurfaceVertices();
  computeInterfaceAreas();
  computeBlendedSurfaceMesh();
  cout << "done." << endl;

  // compute the Us for just the nodes between partitions
  cout << " Computing interface bases ... ";
  flush(cout);
  computeInterfaceUs();
  cout << "done." << endl;

  // if we are simulating, compute a bunch of stuff
  if (simulate)
  {
    // compute the Us for just the nodes between partitions
    cout << " Computing interface rest pose differences ... ";
    flush(cout);
    computeRestDiffs();
    cout << "done." << endl;
   
    // compute the sandwiches between all the interfaces
    cout << " Computing interface sandwiches ... ";
    flush(cout);
    computeSandwiches();
    cout << "done." << endl;

    // compute projected interface rest poses
    cout << " Computing projected rest interfaces ... ";
    flush(cout);
    computeProjectedRestInterfaces();
    cout << "done." << endl;
  }
  
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
    cout << " No graph coloring found, computing one ... " << endl;
    computeGraphColoring(partitionPath);

    // try reading in again
    readDIMACS(dimacs.c_str());
  }
  cout << "done." << endl;

  // compute the total number of interfaces
  _totalInterfaces = 0;
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < x; y++)
      if (_interfaceU[x][y].rows() != 0)
        _totalInterfaces++;

  // recompute inertial variables due to mass adjustments
  cout << " Renormalizing masses ... "; flush(cout);
  for (int x = 0; x < _partitions; x++)
  {
    string cacheName = _meshes[x]->filename();
    cacheName = cacheName + string(".normalized.inertia");
    SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_meshes[x];

    // if there's no cubature yet, don't cache an incorrect matrix
    if (mesh->totalKeyTets() == 0) continue;

    if (!mesh->readInertiaCache(cacheName))
    {
      cout << " No cache found for partition " << x << " ... "; flush(cout);
      ((SUBSPACE_TET_MESH*)_meshes[x])->cacheMassMatrixVars();
      mesh->writeInertiaCache(cacheName);
    }
    else
      cout << " Cache found for partition " << x << "! ... "; flush(cout);
  }
  cout << " done." << endl;

  // cache the colors
  _colors[0][0] = 1.0f; _colors[0][1] = 0.7f; _colors[0][2] = 0.7f; _colors[0][3] = 1.0f;
  _colors[1][0] = 0.7f; _colors[1][1] = 1.0f; _colors[1][2] = 0.7f; _colors[1][3] = 1.0f;
  _colors[2][0] = 0.7f; _colors[2][1] = 0.7f; _colors[2][2] = 0.7f; _colors[2][3] = 1.0f;
  _colors[3][0] = 1.0f; _colors[3][1] = 0.7f; _colors[3][2] = 1.0f; _colors[3][3] = 1.0f;
  _colors[4][0] = 0.7f; _colors[4][1] = 1.0f; _colors[4][2] = 1.0f; _colors[4][3] = 1.0f;
  _colors[5][0] = 1.0f; _colors[5][1] = 1.0f; _colors[5][2] = 0.7f; _colors[5][3] = 1.0f;
  _colors[6][0] = 1.0f; _colors[6][1] = 1.0f; _colors[6][2] = 1.0f; _colors[6][3] = 1.0f;
  _colors[7][0] = 0.5f; _colors[7][1] = 0.5f; _colors[7][2] = 0.5f; _colors[7][3] = 1.0f;
}

PARTITIONED_SUBSPACE_TET_MESH::PARTITIONED_SUBSPACE_TET_MESH() : 
  _interfaceU(NULL),
  _restDiffs(NULL),
  _restSandwich(NULL),
  _basisSandwich(NULL),
  _translationSandwich(NULL),
  _constraintRestSandwich(NULL),
  _constraintBasisSandwich(NULL),
  _meanRestInterface(NULL),
  _meanBasisInterface(NULL)
{
}

PARTITIONED_SUBSPACE_TET_MESH::~PARTITIONED_SUBSPACE_TET_MESH()
{
  // clean up interface Us
  if (_interfaceU != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _interfaceU[x];
    delete[] _interfaceU;
  }

  // clean up interface diffs
  if (_restDiffs != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _restDiffs[x];
    delete[] _restDiffs;
  }
 
  // clean up any old sandwiches
  if (_restSandwich != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _restSandwich[x];
    delete[] _restSandwich;
  }
  if (_basisSandwich != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _basisSandwich[x];
    delete[] _basisSandwich;
  }
  if (_translationSandwich != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _translationSandwich[x];
    delete[] _translationSandwich;
  }
 
  if (_constraintRestSandwich != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _constraintRestSandwich[x];
    delete[] _constraintRestSandwich;
  }
  if (_constraintBasisSandwich != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _constraintBasisSandwich[x];
    delete[] _constraintBasisSandwich;
  }

  if (_meanRestInterface != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _meanRestInterface[x];
    delete[] _meanRestInterface;
  }
  if (_meanBasisInterface != NULL)
  {
    for (int x = 0; x < _partitions; x++)
      delete[] _meanBasisInterface[x];
    delete[] _meanBasisInterface;
  }
  delete[] _originalIDs;

  delete[] _ranks;
}

//////////////////////////////////////////////////////////////////////
// Distribute the _UBasis of _originalMesh to the partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::distributeSubbases()
{
  // get the U basis
  MATRIX& UBasis = ((SUBSPACE_TET_MESH*)_originalMesh)->U();
  VECTOR& eigenvalues = ((SUBSPACE_TET_MESH*)_originalMesh)->eigenvalues();

  // for each partition, extract the rows corresponding to each
  // vertex and concatenate them into a basis
  for (int x = 0; x < _partitions; x++)
  {
    // get the basis of the partition
    MATRIX& subbasis = ((SUBSPACE_TET_MESH*)_meshes[x])->U();
    VECTOR& subvalues = ((SUBSPACE_TET_MESH*)_meshes[x])->eigenvalues();

    // resize it appropriately
    int partitionNodes = _meshes[x]->unconstrainedNodes();
    subbasis.resizeAndWipe(partitionNodes * 3, UBasis.cols());

    // for each vertex in the partition
    for (int y = 0; y < partitionNodes; y++)
    {
      // get the vertex index in the original mesh
      int originalID = _originalIDs[x][y];

      // get the submatrix corresponding to that vertex
      SUBMATRIX vertexBasis(UBasis, originalID * 3, 3);

      // copy the vertex basis into the subbasis
      vertexBasis.copiesInto(subbasis, 3 * y);
    }

    // copy all the eigenvalues
    subvalues = eigenvalues;
  }
}

//////////////////////////////////////////////////////////////////////
// Read in existing cubature files
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::readCubature()
{
  for (int x = 0; x < _partitions; x++)
    ((SUBSPACE_TET_MESH*)_meshes[x])->readCubature();
}

//////////////////////////////////////////////////////////////////////
// Compute the basis between partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::computeInterfaceUs()
{
  // allocate the matrices
  _interfaceU = new MATRIX*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _interfaceU[x] = new MATRIX[_partitions];

  // populate the matrices
  for (int x = 0; x < _partitions; x++)
    // construct the matrices _interfaceU[x][y] and
    // _interfaceU[y][x] simulteously in order to guarantee
    // that the vertex ordering is consistent
    for (int y = x + 1; y < _partitions; y++)
    {
      // get the clones between the two interfaces
      vector<pair<int, int> >& clones = _clonedVertices[x][y];

      int totalClones = clones.size();
      //int rank = ((SUBSPACE_TET_MESH*)(_meshes[0]))->rank();
      int rank;
      
      // set the sizes of the matrices
      rank = ((SUBSPACE_TET_MESH*)(_meshes[x]))->rank();
      _interfaceU[x][y].resizeAndWipe(3 * totalClones, rank);
      rank = ((SUBSPACE_TET_MESH*)(_meshes[y]))->rank();
      _interfaceU[y][x].resizeAndWipe(3 * totalClones, rank);

      // if there's nothing to do, move on
      if (totalClones == 0) continue;

      // get the basis for partition x
      MATRIX& xU = ((SUBSPACE_TET_MESH*)(_meshes[x]))->U();

      // get the basis for partition y
      MATRIX& yU = ((SUBSPACE_TET_MESH*)(_meshes[y]))->U();
      
      // copy in the subbases for each clone
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int xIndex = clones[z].first;
        int yIndex = clones[z].second;

        // if the node is constrained, move on. The U basis does not
        // have rows to accomodate these.
        if (xIndex >=  _meshes[x]->unconstrainedNodes()) continue;

        // get the subbasis for the x vertex
        SUBMATRIX xVertexBasis(xU, xIndex * 3, 3);
       
        // copy into the interface basis
        xVertexBasis.copiesInto(_interfaceU[x][y], 3 * z);
        
        // get the subbasis for the y vertex
        SUBMATRIX yVertexBasis(yU, yIndex * 3, 3);
       
        // copy into the interface basis
        yVertexBasis.copiesInto(_interfaceU[y][x], 3 * z);
      }
    }
}

//////////////////////////////////////////////////////////////////////
// properly ordered rest pose vertices along interface x,y
//////////////////////////////////////////////////////////////////////
VECTOR PARTITIONED_SUBSPACE_TET_MESH::interfaceRests(int x, int y)
{
  // get the clones between the two interfaces
  vector<pair<int, int> >& clones = _clonedVertices[x][y];

  // set the sizes of the final vector
  int totalClones = clones.size();
  VECTOR final(3 * totalClones);

  // if there's nothing to do, move on
  if (totalClones == 0) return final;

  // copy in the subbases for each clone
  for (unsigned int z = 0; z < clones.size(); z++)
  {
    int index = clones[z].first;
    VEC3F& rest = *_meshes[x]->restVertices(index);

    final[3 * z] = rest[0];
    final[3 * z + 1] = rest[1];
    final[3 * z + 2] = rest[2];
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// Compute the difference between rest poses between partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::computeRestDiffs()
{
  // allocate the matrices
  _restDiffs = new VECTOR*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _restDiffs[x] = new VECTOR[_partitions];

  // populate the vectors
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      // get the clones between the two interfaces
      vector<pair<int, int> >& clones = _clonedVertices[x][y];

      // see if there is even an interface here
      if (clones.size() == 0 || x == y) continue;

      VECTOR diffs(3 * clones.size());
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int xIndex = clones[z].first;
        int yIndex = clones[z].second;

        // get the rest pose vertex for the first domain
        VEC3F* xVertex = _meshes[x]->restVertices(xIndex);
        VEC3F* yVertex = _meshes[y]->restVertices(yIndex);

        VEC3F diff = (*xVertex) - (*yVertex);
        //diff *= _springConst;
        diff *= _interfaceSpringConst(x,y);
        diff *= _interfaceArea[x][y] / _maxInterfaceArea;
        diffs[3 * z] = diff[0];
        diffs[3 * z + 1] = diff[1];
        diffs[3 * z + 2] = diff[2];
      }
      _restDiffs[x][y] = _interfaceU[x][y] ^ diffs;
    }
}

//////////////////////////////////////////////////////////////////////
// Verify the basis
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::verifyInterfaceBasis(int index)
{
  // strip off the first interface clone
  vector<pair<int, int> >& clones = _clonedVertices[0][1];
  cout << " cloned nodes: " << clones.size() << endl;

  cout << " clones: " << endl;
  for (unsigned int x = 0; x < clones.size(); x++)
    cout << clones[x].first << " " << clones[x].second << endl;

  int xIndex = clones[index].first;
  int yIndex = clones[index].second;

  VEC3F* xVertex = _meshes[0]->vertices(xIndex);
  VEC3F* yVertex = _meshes[1]->vertices(yIndex);
  VEC3F* xRest = _meshes[0]->restVertices(xIndex);
  VEC3F* yRest = _meshes[1]->restVertices(yIndex);

  cout << " Vertex " << index << " in mesh 0: " << *xVertex << endl;
  cout << " Vertex " << index << " in mesh 1: " << *yVertex << endl;
  cout << " Displacement " << index << " in mesh 0: " << *xVertex - *xRest << endl;
  cout << " Displacement " << index << " in mesh 1: " << *yVertex - *yRest << endl;
  VEC3F diff = (*xVertex) - (*yVertex);
  cout << " distance: " << diff << " " << norm2(diff) << endl;

  // get the subbasis for this vertex
  MATRIX& xInterface = _interfaceU[0][1];
  MATRIX& yInterface = _interfaceU[1][0];
  SUBMATRIX xVertexBasis(xInterface, 3 * index, 3);

  // this is going to get really confusing...
  SUBMATRIX yVertexBasis(yInterface, 3 * index, 3);

  // get the q's for each mesh
  VECTOR& xQ = ((SUBSPACE_TET_MESH*)_meshes[0])->q();
  VECTOR& yQ = ((SUBSPACE_TET_MESH*)_meshes[1])->q();

  // unproject
  VECTOR xFinal = xVertexBasis * xQ;
  VECTOR yFinal = yVertexBasis * yQ;
  cout << " Displacement " << index << " according to mesh 0 interface basis: " << xFinal << endl;
  cout << " Displacement " << index << " according to mesh 1 interface basis: " << yFinal << endl;
  xFinal(0) += (*xRest)[0];
  xFinal(1) += (*xRest)[1];
  xFinal(2) += (*xRest)[2];
  yFinal(0) += (*yRest)[0];
  yFinal(1) += (*yRest)[1];
  yFinal(2) += (*yRest)[2];

  cout << " Vertex " << index << " according to mesh 0 interface basis: " << xFinal << endl;
  cout << " Vertex " << index << " according to mesh 1 interface basis: " << yFinal << endl;
}

//////////////////////////////////////////////////////////////////////
// Gather all the partitions qs into _superQ
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::gatherSuperQ()
{ 
  if (_superQ.size() != _superRank)
    _superQ.resizeAndWipe(_superRank);

  int i = 0;
  for (int x = 0; x < _partitions; x++)
  {
    VECTOR& q = ((SUBSPACE_TET_MESH*)_meshes[x])->q();
    int rank = q.size();
    for (int y = 0; y < rank; y++, i++)
      _superQ(i) = q(y);
  }
}

//////////////////////////////////////////////////////////////////////
// Gather all the partitions qs into _superQ
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::gatherSuperQ(BLOCK_VECTOR& superQ)
{
  for (int x = 0; x < _partitions; x++)
  {
    VECTOR& q = ((SUBSPACE_TET_MESH*)_meshes[x])->q();
    (*superQ(x)) = q;
  }
}
 
//////////////////////////////////////////////////////////////////////
// Gather all the partitions qOlds into _superQ
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::gatherSuperQOld(BLOCK_VECTOR& superQOld)
{
  for (int x = 0; x < _partitions; x++)
  {
    VECTOR& qOld = ((SUBSPACE_TET_MESH*)_meshes[x])->qOld();
    (*superQOld(x)) = qOld;
  }
}

//////////////////////////////////////////////////////////////////////
// Scatter _superQ into all partition qs
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::scatterSuperQ()
{ 
  int i = 0;
  for (int x = 0; x < _partitions; x++)
  {
    VECTOR& q = ((SUBSPACE_TET_MESH*)_meshes[x])->q();
    int rank = q.size();
    for (int y = 0; y < rank; y++, i++)
      q(y) = _superQ(i);
  }
}

//////////////////////////////////////////////////////////////////////
// Scatter new qOld into all partition qs
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::scatterSuperQOld(BLOCK_VECTOR& superQOld)
{ 
  for (int x = 0; x < _partitions; x++)
  {
    VECTOR& qOld = ((SUBSPACE_TET_MESH*)_meshes[x])->qOld();
    qOld.equals(*superQOld(x));
  }
}

//////////////////////////////////////////////////////////////////////
// Scatter new q into all partition qs
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::scatterSuperQ(BLOCK_VECTOR& superQ)
{ 
  for (int x = 0; x < _partitions; x++)
  {
    VECTOR& q = ((SUBSPACE_TET_MESH*)_meshes[x])->q();
    q.equals(*superQ(x));
  }
}

//////////////////////////////////////////////////////////////////////
// Compute the stiffness matrix of the interface
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::computeInterfaceStiffness()
{
  if (_interfaceStiffness.rows() != _superRank ||
      _interfaceStiffness.cols() != _superRank)
    _interfaceStiffness.resize(_superRank, _superRank);

  // get the prefix sum of the partition ranks
  int* prefixRanks = new int[_partitions];
  prefixRanks[0] = 0;
  if (_partitions > 1)
    prefixRanks[1] = _ranks[0];
  for (int x = 2; x < _partitions; x++)
    prefixRanks[x] = _ranks[x-1] + prefixRanks[x-1];

  // do each combination
  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
    {
      // see if there is even an interface here
      if (_interfaceU[x][y].rows() == 0) continue;

      // create the spring constant matrices
      int leftRows = _interfaceU[x][y].rows();
      MATRIX springConstLeft(leftRows, leftRows);
      springConstLeft.setToIdentity();
      springConstLeft *= _interfaceSpringConst(x,y);
      
      int rightRows = _interfaceU[y][x].rows();
      MATRIX springConstRight(rightRows, rightRows);
      springConstRight.setToIdentity();
      springConstRight *= _interfaceSpringConst(y,x);

      // weight by the area
      springConstLeft *= _interfaceArea[x][y] / _maxInterfaceArea;
      springConstRight *= _interfaceArea[y][x] / _maxInterfaceArea;

      // store self-products
      MATRIX spring = springConstLeft * _interfaceU[x][y];
      MATRIX self = _interfaceU[x][y] ^ spring;

      _interfaceStiffness.add(self, prefixRanks[x], prefixRanks[x]);

      spring = springConstRight * _interfaceU[y][x];
      self = _interfaceU[y][x] ^ spring;
      _interfaceStiffness.add(self, prefixRanks[y], prefixRanks[y]);

      // store off-diagonal products
      spring = springConstRight * _interfaceU[y][x];
      MATRIX offDiagonal = _interfaceU[x][y] ^ spring;
      _interfaceStiffness.subtract(offDiagonal, prefixRanks[x], prefixRanks[y]);
      MATRIX transpose = offDiagonal.transpose();
      _interfaceStiffness.subtract(transpose, prefixRanks[y], prefixRanks[x]);
    }

  delete[] prefixRanks;
}

//////////////////////////////////////////////////////////////////////
// Compute the stiffness matrix of the interface in block form
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::computeBlockInterfaceStiffness()
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

        springConstLeft(3 * leftIndex, 3 * leftIndex) = k;
        springConstLeft(3 * leftIndex + 1, 3 * leftIndex + 1) = k;
        springConstLeft(3 * leftIndex + 2, 3 * leftIndex + 2) = k;
        
        springConstRight(3 * rightIndex, 3 * rightIndex) = k;
        springConstRight(3 * rightIndex + 1, 3 * rightIndex + 1) = k;
        springConstRight(3 * rightIndex + 2, 3 * rightIndex + 2) = k;
      }

      // store self-products
      MATRIX& Ux = ((SUBSPACE_TET_MESH*)_meshes[x])->U();
      MATRIX spring = springConstLeft * Ux;
      MATRIX self = Ux ^ spring;
      _diagonalInterfaceStiffness.add(self, x, x);

      MATRIX& Uy = ((SUBSPACE_TET_MESH*)_meshes[y])->U();
      spring = springConstRight * Uy;
      self = Uy ^ spring;
      _diagonalInterfaceStiffness.add(self, y, y);

      // generate off-diagonal matrix
      SPARSE_MATRIX offDiagonalFull(leftDofs, rightDofs);
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int leftIndex = clones[z].first;
        int rightIndex = clones[z].second;

        offDiagonalFull(3 * leftIndex, 3 * rightIndex)         = -k;
        offDiagonalFull(3 * leftIndex + 1, 3 * rightIndex + 1) = -k;
        offDiagonalFull(3 * leftIndex + 2, 3 * rightIndex + 2) = -k;
      }

      // store off-diagonal products
      spring = offDiagonalFull * Uy;
      MATRIX offDiagonal = Ux ^ spring;
      _offDiagonalInterfaceStiffness.equals(offDiagonal, x, y);
      MATRIX transpose = offDiagonal.transpose();
      _offDiagonalInterfaceStiffness.equals(transpose, y, x);
    }
}

//////////////////////////////////////////////////////////////////////
// Compute the stiffness matrix of the interface in block form
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::updateBlockInterfaceStiffness()
{
  // do each combination
  for (int x = 0; x < _partitions; x++)
    for (int y = x + 1; y < _partitions; y++)
    {
      // see if there is even an interface here
      if (_interfaceU[x][y].rows() == 0) continue;

      MATRIX3 Ri = quaternionRotation(x).toExplicitMatrix3x3();
      MATRIX3 Rj = quaternionRotation(y).toExplicitMatrix3x3();

      MATRIX3 RiTRj = Ri.transpose() * Rj;
      MATRIX sandwich = _basisSandwich[x][y].transform(RiTRj);
      sandwich *= -springConst(x,y);

      _offDiagonalInterfaceStiffness.equals(sandwich, x, y);
      
      MATRIX transpose = sandwich.transpose();
      _offDiagonalInterfaceStiffness.equals(transpose, y, x);
    }
}

//////////////////////////////////////////////////////////////////////
// combine the vertex positions in all the partitions and 
// write out to Matlab format
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::verticesToMatlab(const char* filename, const char* varName, int column)
{
  // total nodes in the original mesh;
  int totalNodes = _originalMesh->totalNodes();
  
  // make sure all the meshes are updated
  for (int x = 0; x < _partitions; x++)
    _meshes[x]->updateFullMesh();

  // allocate a vector large enough to store all the vertex positions
  vector<VEC3F> vertices;
  vertices.resize(totalNodes);

  // track how many values were stored to each vertex so they can
  // be averaged later
  vector<int> totalSeen;
  totalSeen.resize(totalNodes);
  for (int x = 0; x < totalNodes; x++)
    totalSeen[x] = 0;

  // for each partition
  for (int x = 0; x < _partitions; x++)
  {
    // get its vertex list
    vector<VEC3F>& partitionVertices = _meshes[x]->vertices();
    for (unsigned int y = 0; y < partitionVertices.size(); y++)
    {
      // see what vertex this was in the original mesh
      int originalID = _originalIDs[x][y];

      // store it
      vertices[originalID] += partitionVertices[y];
        
      // increment the number of vertices seen in this slot
      totalSeen[originalID]++;
    }
  }

  // average the vertices that saw more than one instance
  for (int x = 0; x < totalNodes; x++)
    if (totalSeen[x] > 1)
      vertices[x] *= 1.0 / totalSeen[x];
    else if (totalSeen[x] == 0)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " No instances of a vertex seen!!!!!" << endl;
    }

  // open a file or writing
  FILE* file = NULL;
  if (column == 0)
    file = fopen(filename, "w");
  else
    file = fopen(filename, "a");
  fprintf(file, "%s(:,%i) = [\n", varName, column + 1);
  for (int x = 0; x < totalNodes; x++)
    fprintf(file, "%f\n%f\n%f\n", vertices[x][0], vertices[x][1], vertices[x][2]);
  fprintf(file, "];\n");
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// combine the vertex positions in all the partitions and 
// write out to binary
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::verticesToBinary(const char* filename, int column)
{
  // total nodes in the original mesh;
  int totalNodes = _originalMesh->totalNodes();
  
  // make sure all the meshes are updated
  for (int x = 0; x < _partitions; x++)
    _meshes[x]->updateFullMesh();

  // allocate a vector large enough to store all the vertex positions
  vector<VEC3F> vertices;
  vertices.resize(totalNodes);

  // track how many values were stored to each vertex so they can
  // be averaged later
  vector<int> totalSeen;
  totalSeen.resize(totalNodes);
  for (int x = 0; x < totalNodes; x++)
    totalSeen[x] = 0;

  // for each partition
  for (int x = 0; x < _partitions; x++)
  {
    // get its vertex list
    vector<VEC3F>& partitionVertices = _meshes[x]->vertices();
    for (unsigned int y = 0; y < partitionVertices.size(); y++)
    {
      // see what vertex this was in the original mesh
      int originalID = _originalIDs[x][y];

      // store it
      vertices[originalID] += partitionVertices[y];
        
      // increment the number of vertices seen in this slot
      totalSeen[originalID]++;
    }
  }

  // average the vertices that saw more than one instance
  for (int x = 0; x < totalNodes; x++)
    if (totalSeen[x] > 1)
      vertices[x] *= 1.0 / totalSeen[x];
    else if (totalSeen[x] == 0)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " No instances of a vertex seen!!!!!" << endl;
    }

  // open a file or writing
  FILE* file = NULL;
  if (column == 0)
  {
    file = fopen(filename, "wb");
    int totalScalars = totalNodes * 3;
    fwrite((void*)&totalScalars, sizeof(int), 1, file);
  }
  else
    file = fopen(filename, "ab");
  for (int x = 0; x < totalNodes; x++)
    for (int y = 0; y < 3; y++)
    {
      Real scalar = vertices[x][y];
      fwrite((void*)&scalar, sizeof(Real), 1, file);
    }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::resetRigidTranslations()
{
  for (unsigned int x = 0; x < _unconstrainedPartition.size(); x++)
    if (_unconstrainedPartition[x])
      ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x])->resetRigidTranslation();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::updateRotations()
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " updateRotations() disabled! " << endl;
}

//////////////////////////////////////////////////////////////////////
// Draw centers of mass for each rigid frame
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::drawCentersOfMass()
{
#if USING_GLVU  
  for (int x = 0; x < _partitions; x++)
    if (_unconstrainedPartition[x])
      ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x])->drawCenterOfMass();
#endif
}

//////////////////////////////////////////////////////////////////////
// Draw rotated axes for each rigid frame
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::drawRigidFrames()
{
#if USING_GLVU  
  for (int x = 0; x < _partitions; x++)
    if (_unconstrainedPartition[x])
      ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x])->drawRigidFrame();
#endif
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::drawRigidDebug()
{
#if USING_GLVU  
  for (int x = 0; x < _partitions; x++)
  {
    float color[4];
    color[0] = 1;
    color[1] = 1;
    color[2] = 1;
    color[3] = 1.0f;
    glColor4f(color[0], color[1], color[2], color[3]);

    if (_unconstrainedPartition[x])
      ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x])->drawRigidDebug();
  }
#endif
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::drawRigidOnly()
{
#if USING_GLVU  
  for (int x = 0; x < _partitions; x++)
  {
    float color[4];
    color[0] = 1.0; color[1] = 1.0; color[2] = 1.0;
    color[3] = 1.0f;
    glColor4f(color[0], color[1], color[2], color[3]);

    if (_unconstrainedPartition[x])
      ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x])->drawRigidOnly();
  }
#endif
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::computeSandwiches()
{
  _restSandwich  = new SANDWICH_TRANSFORM*[_partitions];
  _basisSandwich = new SANDWICH_TRANSFORM*[_partitions];
  _translationSandwich = new SANDWICH_TRANSFORM*[_partitions];
  _constraintRestSandwich = new SANDWICH_TRANSFORM*[_partitions];
  _constraintBasisSandwich = new SANDWICH_TRANSFORM*[_partitions];
  _meanRestInterface = new VEC3F*[_partitions];
  _meanBasisInterface = new MATRIX*[_partitions];

  for (int x = 0; x < _partitions; x++)
  {
    _restSandwich[x] = new SANDWICH_TRANSFORM[_partitions];
    _basisSandwich[x] = new SANDWICH_TRANSFORM[_partitions];
    _translationSandwich[x] = new SANDWICH_TRANSFORM[_partitions];
    _constraintRestSandwich[x] = new SANDWICH_TRANSFORM[_partitions];
    _constraintBasisSandwich[x] = new SANDWICH_TRANSFORM[_partitions];
    _meanRestInterface[x] = new VEC3F[_partitions];
    _meanBasisInterface[x] = new MATRIX[_partitions];
  }

  // do each combination
  for (int x = 0; x < _partitions; x++)
  {
    for (int y = 0; y < _partitions; y++)
    {
      // see if there is even an interface here
      if (_interfaceU[x][y].rows() == 0 || x == y) continue;

      // do the basis sandwich first (easiest)
      _basisSandwich[x][y].init3x3(_interfaceU[x][y], _interfaceU[y][x]);

      // build a vector of the interface rest positions
      vector<pair<int, int> >& clones = _clonedVertices[x][y];
      VECTOR rests(3 * clones.size());

      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int index = clones[z].second;

        // get the rest pose vertex for the first domain
        VEC3F rest = *(_meshes[y]->restVertices(index));
        rests[3 * z] = rest[0];
        rests[3 * z + 1] = rest[1];
        rests[3 * z + 2] = rest[2];
      }

      // build the rest pose sandwich
      Real springConst = _interfaceSpringConst(x,y);
      springConst *= _interfaceArea[x][y] / _maxInterfaceArea;
      MATRIX springInterface = springConst * _interfaceU[x][y];
      _restSandwich[x][y].init3x3(springInterface, rests);

      // build the blocked identity matrix for the translation
      MATRIX identities(3 * clones.size(), 3);
      for (unsigned int z = 0; z < 3 * clones.size(); z++)
        identities(z, z % 3) = 1.0;

      _translationSandwich[x][y].init3x3(_interfaceU[x][y], identities);

      // for the hard constraints, need the specific vertices for *this*
      // interface (x), not the other one (y), because the mean may already
      // have been subtracted out
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int index = clones[z].first;

        // get the rest pose vertex for the first domain
        VEC3F rest = *(_meshes[x]->restVertices(index));
        rests[3 * z] = rest[0];
        rests[3 * z + 1] = rest[1];
        rests[3 * z + 2] = rest[2];
      }

      // for the hard constraints, normalize for total number of vertices
      identities *= 1.0 / clones.size();

      // compute the hard constraint versions
      _constraintRestSandwich[x][y].init3x3(identities, rests);
      _constraintBasisSandwich[x][y].init3x3(identities, _interfaceU[x][y]);

      // compute the mean rest pose position for partitions that are
      // constrained (just do it for everybody -- we'll need it eventually
      // for the very stiff case anyway)
      VECTOR meanRest = identities ^ rests;
      _meanRestInterface[x][y][0] = meanRest[0];
      _meanRestInterface[x][y][1] = meanRest[1];
      _meanRestInterface[x][y][2] = meanRest[2];
    
      // Do the same for the interface U
      _meanBasisInterface[x][y] = identities ^ _interfaceU[x][y];
    }
  }
}

//////////////////////////////////////////////////////////////////////
// return a block identity matrix the dimensions of the interface at (x,y)
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX PARTITIONED_SUBSPACE_TET_MESH::interfaceIdentity(int x, int y)
{
  if (!interfaceExists(x,y))
    return SPARSE_MATRIX(0,0);

  int totalClones = _clonedVertices[x][y].size();
  SPARSE_MATRIX final(3 * totalClones, 3);

  for (int i = 0; i < totalClones; i++)
  {
    final(3 * i, 0) = 1.0;
    final(3 * i + 1, 1) = 1.0;
    final(3 * i + 2, 2) = 1.0;
  }
  return final;
}

//////////////////////////////////////////////////////////////////////
// Return the rigid rotation for a given partition
//////////////////////////////////////////////////////////////////////
MATRIX3 PARTITIONED_SUBSPACE_TET_MESH::rigidRotation(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[partition])->rotationQuaternion().toExplicitMatrix3x3();

  return MATRIX3::I();
}

//////////////////////////////////////////////////////////////////////
// Return the rigid rotation for a given partition
//////////////////////////////////////////////////////////////////////
QUATERNION PARTITIONED_SUBSPACE_TET_MESH::quaternionRotation(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[partition])->rotationQuaternion();

  return QUATERNION(MATRIX3::I());
}

//////////////////////////////////////////////////////////////////////
// Return the rigid rotation for a given partition from the previous
// timestep
//////////////////////////////////////////////////////////////////////
QUATERNION PARTITIONED_SUBSPACE_TET_MESH::quaternionRotationOld(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[partition])->rotationQuaternionOld();

  return QUATERNION(MATRIX3::I());
}

//////////////////////////////////////////////////////////////////////
// Return the lagrange multiplier from the previous solve
//////////////////////////////////////////////////////////////////////
Real PARTITIONED_SUBSPACE_TET_MESH::rotationLambda(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[partition])->rotationLambda();

  return 0.0;
}

//////////////////////////////////////////////////////////////////////
// Return the rigid rotation for a given partition
//////////////////////////////////////////////////////////////////////
MATRIX3 PARTITIONED_SUBSPACE_TET_MESH::rigidRotationOld(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[partition])->rotationQuaternionOld().toExplicitMatrix3x3();

  return MATRIX3::I();
}

//////////////////////////////////////////////////////////////////////
// Return the rigid rotation for a given partition
//////////////////////////////////////////////////////////////////////
VEC3F PARTITIONED_SUBSPACE_TET_MESH::rigidTranslation(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[partition])->rigidTranslation();

  VEC3F zero;
  return zero;
}

//////////////////////////////////////////////////////////////////////
// Return the rigid rotation for a given partition
//////////////////////////////////////////////////////////////////////
VEC3F PARTITIONED_SUBSPACE_TET_MESH::rigidTranslationOld(int partition)
{
  if (_unconstrainedPartition[partition])
    return ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[partition])->rigidTranslationOld();

  VEC3F zero;
  return zero;
}

//////////////////////////////////////////////////////////////////////
// Compute projected rest interfaces
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::computeProjectedRestInterfaces()
{
  _projectedRestInterfaces.resizeAndWipe(_partitions);
  
  // populate the vectors
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      // get the clones between the two interfaces
      vector<pair<int, int> >& clones = _clonedVertices[x][y];

      // see if there is even an interface here
      if (clones.size() == 0 || x == y) continue;

      VECTOR rests(3 * clones.size());
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int index = clones[z].first;

        // get the rest pose vertex for the first domain
        VEC3F vertex = *(_meshes[x]->restVertices(index));

        rests[3 * z]     = vertex[0];
        rests[3 * z + 1] = vertex[1];
        rests[3 * z + 2] = vertex[2];
      }
      VECTOR projectedRest = _interfaceU[x][y] ^ rests;
      projectedRest *= _interfaceSpringConst(x,y);
      projectedRest *= _interfaceArea[x][y] / _maxInterfaceArea;
      _projectedRestInterfaces.add(projectedRest, x);
    }
}

//////////////////////////////////////////////////////////////////////
// Find the closest surface node
//////////////////////////////////////////////////////////////////////
int PARTITIONED_SUBSPACE_TET_MESH::closestPartition(VEC3F point)
{
  Real minDistance = 1e9;
  int minPartition = 0;
  
  for (int x = 0; x < _partitions; x++)
  {
    VEC3F* candidate = _meshes[x]->closestSurfaceNode(point);
    VEC3F transformed = *candidate;

    if (_unconstrainedPartition[x])
      transformed = rigidRotation(x) * transformed + rigidTranslation(x);
    Real distance = norm2(point - transformed);
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
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::writeSubbasis(int x)
{
  SUBSPACE_TET_MESH* tetMesh = (SUBSPACE_TET_MESH*)_meshes[x];
  string filename = tetMesh->filename();
  filename = filename + string(".eigenvectors.matrix");
  MATRIX& basis = tetMesh->U();
  basis.write(filename.c_str());

  filename = tetMesh->filename() + string(".eigenvalues.vector");
  VECTOR& values = tetMesh->eigenvalues();
  values.write(filename.c_str());
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::writeSubbases()
{
  // save out the partition bases
  for (int x = 0; x < _partitions; x++)
    writeSubbasis(x);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VEC3F PARTITIONED_SUBSPACE_TET_MESH::computeInterfaceMean(int x, int y)
{
  VEC3F zero;

  // if it's on the diagonal, skip it
  if (x == y) return zero;

  // if there are no clones, skip it
  vector<pair<int,int> >& clonedVertices = _clonedVertices[x][y];
  if (clonedVertices.size() == 0) return zero;

  // just to be safe
  _meshes[x]->updateFullMesh();

  // take the mean of all the vertices
  VEC3F mean;
  for (unsigned int z = 0; z < clonedVertices.size(); z++)
  {
    VEC3F& vertex = *(_meshes[x]->vertices(clonedVertices[z].first));
    mean += vertex;
  }
  mean *= 1.0 / clonedVertices.size();
  mean = rigidRotation(x) * mean;

  // add the interface translation
  VEC3F translation = rigidTranslation(x);
  mean += translation;

  return mean;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VEC3F PARTITIONED_SUBSPACE_TET_MESH::computeInterfaceMeanRest(int x, int y)
{
  VEC3F zero;

  // if it's on the diagonal, skip it
  if (x == y) return zero;

  // if there are no clones, skip it
  vector<pair<int,int> >& clonedVertices = _clonedVertices[x][y];
  if (clonedVertices.size() == 0) return zero;

  // just to be safe
  _meshes[x]->updateFullMesh();

  // take the mean of all the vertices
  VEC3F mean;
  for (unsigned int z = 0; z < clonedVertices.size(); z++)
  {
    VEC3F& vertex = *(_meshes[x]->restVertices(clonedVertices[z].first));
    mean += vertex;
  }
  mean *= 1.0 / clonedVertices.size();
  mean = rigidRotation(x) * mean;

  // add the interface translation
  VEC3F translation = rigidTranslation(x);
  mean += translation;

  return mean;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::drawInterfaceMeans()
{
  glDisable(GL_DEPTH);
  glDisable(GL_BLEND);

  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < x; y++)
    {
      if (x == y) continue;
      vector<pair<int,int> >& clonedVertices = _clonedVertices[x][y];
      if (clonedVertices.size() == 0) continue;

      VEC3F meanLeft = computeInterfaceMean(x,y);

      // give it the same color as the domain
      float colorLeft[4];
      int mod = y % 8;
      colorLeft[0] = _colors[mod][0];
      colorLeft[1] = _colors[mod][1];
      colorLeft[2] = _colors[mod][2];
      colorLeft[3] = 1.0;

      glColor4f(colorLeft[0], colorLeft[1], colorLeft[2], colorLeft[3]);
      glPushMatrix();
        glTranslatef(meanLeft[0], meanLeft[1], meanLeft[2]);
        glutSolidSphere(0.025, 10, 10);
      glPopMatrix();

      VEC3F meanRight = computeInterfaceMean(y,x);

      // give it the same color as the domain
      float colorRight[4];
      mod = x % 8;
      colorRight[0] = _colors[mod][0];
      colorRight[1] = _colors[mod][1];
      colorRight[2] = _colors[mod][2];
      colorRight[3] = 1.0;

      glColor4f(colorRight[0], colorRight[1], colorRight[2], colorRight[3]);
      glPushMatrix();
        glTranslatef(meanRight[0], meanRight[1], meanRight[2]);
        glutSolidSphere(0.025, 10, 10);
      glPopMatrix();

      glLineWidth(10.0);
      glBegin(GL_LINES);
        glVertex3f(meanLeft[0], meanLeft[1], meanLeft[2]);
        glVertex3f(meanRight[0], meanRight[1], meanRight[2]);
      glEnd();

      VEC3F diff = meanRight - meanLeft;
    }

  glEnable(GL_BLEND);
  glEnable(GL_DEPTH);
}

//////////////////////////////////////////////////////////////////////
// Update the partitions in the mesh using a subspace update
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::updateSubspaceMeshes()
{
  for (int x = 0; x < _partitions; x++)
  {
    SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_meshes[x];
    mesh->SUBSPACE_TET_MESH::updateFullMesh();
  }
}

//////////////////////////////////////////////////////////////////////
// add the rigid transforms to the surface vertices
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::addRigidsToSurface()
{
  // update all the partitions, applying the rigid transforms
  // as appropriate
  for (int x = 0; x < _partitions; x++)
  {
    SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_meshes[x];
    mesh->updateSurfaceMesh();

    if (_unconstrainedPartition[x])
    {
      UNCONSTRAINED_SUBSPACE_TET_MESH* rigidMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)mesh;
      MATRIX3 rotation = rigidMesh->rotationQuaternion().toExplicitMatrix3x3();
      VEC3F translation = rigidMesh->rigidTranslation();

      vector<VEC3F*>& surfaceVertices = rigidMesh->surfaceVertices();

      for (unsigned int y = 0; y < surfaceVertices.size(); y++)
      {
        VEC3F& vertex = *(surfaceVertices[y]);
        vertex = rotation * vertex;
        vertex += translation;
      }
    }
  }
  
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
}

//////////////////////////////////////////////////////////////////////
// undo the rigid transforms from the surface vertices
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::subtractRigidsFromSurface()
{
  // undo the rigid transforms, just in case
  for (int x = 0; x < _partitions; x++)
  {
    SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_meshes[x];
    mesh->updateSurfaceMesh();
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::drawBlendedMesh()
{
  // update all the partitions, applying the rigid transforms
  // as appropriate
  for (int x = 0; x < _partitions; x++)
  {
    SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_meshes[x];
    mesh->updateSurfaceMesh();

    if (_unconstrainedPartition[x])
    {
      UNCONSTRAINED_SUBSPACE_TET_MESH* rigidMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)mesh;
      MATRIX3 rotation = rigidMesh->rotationQuaternion().toExplicitMatrix3x3();
      VEC3F translation = rigidMesh->rigidTranslation();

      vector<VEC3F*>& surfaceVertices = rigidMesh->surfaceVertices();

      for (unsigned int y = 0; y < surfaceVertices.size(); y++)
      {
        VEC3F& vertex = *(surfaceVertices[y]);
        vertex = rotation * vertex;
        vertex += translation;
      }
    }
  }
  
#if USING_GLVU  
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
#endif

  // undo the rigid transforms, just in case
  for (int x = 0; x < _partitions; x++)
  {
    SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_meshes[x];
    mesh->updateSurfaceMesh();
  }
}

//////////////////////////////////////////////////////////////////////
// update the position of the embedding
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::updateEmbedding()
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

    VEC3F v0 = *(tet.vertices[0]);
    VEC3F v1 = *(tet.vertices[1]);
    VEC3F v2 = *(tet.vertices[2]);
    VEC3F v3 = *(tet.vertices[3]);

    if (_unconstrainedPartition[partition])
    {
      UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[partition];

      VEC3F translation = mesh->rigidTranslation();
      MATRIX3 rotation = mesh->rotationQuaternion().toExplicitMatrix3x3();
      v0 = rotation * v0 + translation;
      v1 = rotation * v1 + translation;
      v2 = rotation * v2 + translation;
      v3 = rotation * v3 + translation;
    }
   
    VEC3F coords = _barycentricEmbeddings[x];
    Real subtract = 1.0 - coords[0] - coords[1] - coords[2];
    VEC3F update = subtract * v3 + v0 * coords[0] + v1 * coords[1] + v2 * coords[2];
    vertices[x] = update;
  }
}

//////////////////////////////////////////////////////////////////////
// set the spring constant for a specific interface to something new, 
// and update all the quantities that need to be updated
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::updateSpringConst(int x, int y, Real newSpringConst)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << "NEED TO ADD COM SANDWICHES HERE" << endl;

  _interfaceSpringConst(x,y) = newSpringConst;
  _interfaceSpringConst(y,x) = newSpringConst;

  // recompute _restSandwich
  for (int i = 0; i < _partitions; i++)
    for (int j = 0; j < _partitions; j++)
    {
      // see if there is even an interface here
      if (_interfaceU[i][j].rows() == 0 || i == j) continue;

      // build the rest pose sandwich
      vector<pair<int, int> >& clones = _clonedVertices[i][j];
      VECTOR rests(3 * clones.size());
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int index = clones[z].second;

        // get the rest pose vertex for the first domain
        VEC3F rest = *(_meshes[j]->restVertices(index));
        rests[3 * z] = rest[0];
        rests[3 * z + 1] = rest[1];
        rests[3 * z + 2] = rest[2];
      }
      //Real springConst = _springConst;
      Real springConst = _interfaceSpringConst(i,j);
      springConst *= _interfaceArea[i][j] / _maxInterfaceArea;
      MATRIX springInterface = springConst * _interfaceU[i][j];
      _restSandwich[i][j].init3x3(springInterface, rests);
    }

  // recompute _projectedRestInterfaces
  computeProjectedRestInterfaces();

  // recompute _interfaceStiffness
  _interfaceStiffness *= 0.0;
  computeInterfaceStiffness();

  // recompute _diagonalInterfaceStiffness 
  // and _offDiagonalInterfaceStiffness
  computeBlockInterfaceStiffness();
}

//////////////////////////////////////////////////////////////////////
// set the spring constant to something new, and update all the
// quantities that need to be updated
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::updateSpringConst(Real newSpringConst)
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << "NEED TO ADD COM SANDWICHES HERE" << endl;

  _springConst = newSpringConst;
  for (int y = 0; y < _partitions; y++)
    for (int x = 0; x < _partitions; x++)
      _interfaceSpringConst(x,y) = newSpringConst;

  // recompute _restSandwich
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
    {
      // see if there is even an interface here
      if (_interfaceU[x][y].rows() == 0 || x == y) continue;

      // build the rest pose sandwich
      vector<pair<int, int> >& clones = _clonedVertices[x][y];
      VECTOR rests(3 * clones.size());
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int index = clones[z].second;

        // get the rest pose vertex for the first domain
        VEC3F rest = *(_meshes[y]->restVertices(index));
        rests[3 * z] = rest[0];
        rests[3 * z + 1] = rest[1];
        rests[3 * z + 2] = rest[2];
      }
      //Real springConst = _springConst;
      Real springConst = _interfaceSpringConst(x,y);
      springConst *= _interfaceArea[x][y] / _maxInterfaceArea;
      MATRIX springInterface = springConst * _interfaceU[x][y];
      _restSandwich[x][y].init3x3(springInterface, rests);
    }

  // recompute _projectedRestInterfaces
  computeProjectedRestInterfaces();

  // recompute _interfaceStiffness
  _interfaceStiffness *= 0.0;
  computeInterfaceStiffness();

  // recompute _diagonalInterfaceStiffness 
  // and _offDiagonalInterfaceStiffness
  computeBlockInterfaceStiffness();
}

//////////////////////////////////////////////////////////////////////
// reset everything to zero
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::reset()
{
  // reset all the all constituent meshes
  for (int x = 0; x < _partitions; x++)
  {
    if (_unconstrainedPartition[x])
    {
      UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x];
      mesh->UNCONSTRAINED_SUBSPACE_TET_MESH::reset();
    }
    else
    {
      SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_meshes[x];
      mesh->SUBSPACE_TET_MESH::reset();
    }
  }

  // recompute the center of mass
  computeCenterOfMass();

  // update the vertex positions
  updateFullMeshes();
}

//////////////////////////////////////////////////////////////////////
// write out an interface spring const matrix
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::writeInterfaceSpringConsts(string directory)
{
  string springFile = directory + string("/interface.spring.matrix");

  cout << " Writing out " << _interfaceSpringConst << endl;
  _interfaceSpringConst.write(springFile.c_str());
}

//////////////////////////////////////////////////////////////////////
// Update centers of mass
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::updateCentersOfMass()
{
  for (int x = 0; x < _partitions; x++)
  {
    _meshes[x]->computeCenterOfMass();

    if (x == 2)
    {
      UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x];
      VEC3F COM = mesh->centerOfMass();
      MATRIX3 rigidRotation = mesh->rotationQuaternion().toExplicitMatrix3x3();
      VEC3F rigidTranslation = mesh->rigidTranslation();
    }
  }
}

//////////////////////////////////////////////////////////////////////
// write out an interface spring const matrix
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::readInterfaceSpringConsts(string directory)
{
  string springFile = directory + string("/interface.spring.matrix");

  // probe for the file
  FILE* file = fopen(springFile.c_str(), "rb");
  if (file == NULL)
  {
    fclose(file);
    cout << " No interface spring constant matrix found! " << endl;
    return;
  }

  // close it so the MATRIX can reopen it
  fclose(file);

  // read it in finally
  _interfaceSpringConst.read(springFile.c_str());

  // sanity check
  assert(_interfaceSpringConst.rows() == _partitions &&
         _interfaceSpringConst.cols() == _partitions);

  // recompute _restSandwich
  for (int i = 0; i < _partitions; i++)
    for (int j = 0; j < _partitions; j++)
    {
      // see if there is even an interface here
      if (_interfaceU[i][j].rows() == 0 || i == j) continue;

      // build the rest pose sandwich
      vector<pair<int, int> >& clones = _clonedVertices[i][j];
      VECTOR rests(3 * clones.size());
      for (unsigned int z = 0; z < clones.size(); z++)
      {
        int index = clones[z].second;

        // get the rest pose vertex for the first domain
        VEC3F rest = *(_meshes[j]->restVertices(index));
        rests[3 * z] = rest[0];
        rests[3 * z + 1] = rest[1];
        rests[3 * z + 2] = rest[2];
      }
      //Real springConst = _springConst;
      Real springConst = _interfaceSpringConst(i,j);
      springConst *= _interfaceArea[i][j] / _maxInterfaceArea;
      MATRIX springInterface = springConst * _interfaceU[i][j];
      _restSandwich[i][j].init3x3(springInterface, rests);
    }

  // recompute _projectedRestInterfaces
  computeProjectedRestInterfaces();

  // recompute _interfaceStiffness
  _interfaceStiffness *= 0.0;
  computeInterfaceStiffness();

  // recompute _diagonalInterfaceStiffness 
  // and _offDiagonalInterfaceStiffness
  computeBlockInterfaceStiffness();
}

//////////////////////////////////////////////////////////////////////
// return a matrix of all the cross product matrices produced by the 
// rest vertices along interface (x,y)
//////////////////////////////////////////////////////////////////////
MATRIX PARTITIONED_SUBSPACE_TET_MESH::interfaceRestTildeU(int x, int y)
{
  vector<pair<int,int> > clones = clonedVertices(x, y);
  vector<VEC3F>& restVertices = _meshes[x]->restPose();

  MATRIX final(3 * clones.size(), 3);
  for (unsigned int i = 0; i < clones.size(); i++)
  {
    int index = clones[i].first;
    VEC3F restVertex = restVertices[index];

    MATRIX restCross = MATRIX::cross(restVertex);
    final.setSubmatrix(restCross, 3 * i);
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a matrix of all the cross product matrices produced by the 
// rest vertices in partition x
//////////////////////////////////////////////////////////////////////
MATRIX PARTITIONED_SUBSPACE_TET_MESH::restTildeU(int x)
{
  vector<VEC3F>& restVertices = _meshes[x]->restPose();

  int totalVertices = ((SUBSPACE_TET_MESH*) _meshes[x])->U().rows() / 3;
  MATRIX final(3 * totalVertices, 3);
  for (int i = 0; i < totalVertices; i++)
  {
    VEC3F restVertex = restVertices[i];

    MATRIX restCross = MATRIX::cross(restVertex);
    final.setSubmatrix(restCross, 3 * i);
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a tensor of all the cross product tensors produced by the 
// displacements in partition x
//////////////////////////////////////////////////////////////////////
TENSOR3 PARTITIONED_SUBSPACE_TET_MESH::tildeU(int x)
{
  MATRIX& U = ((SUBSPACE_TET_MESH*)_meshes[x])->U();
  TENSOR3 final(U.rows(), 3, U.cols());

  for (int i = 0; i < U.rows() / 3; i++)
  {
    SUBMATRIX subU(U, 3 * i, 3);

    TENSOR3 subUTilde = TENSOR3::cross(subU);

    final.setSubtensor(subUTilde, 3 * i);
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a tensor of all the cross product tensors produced by the 
// displacements along interface (x,y)
//////////////////////////////////////////////////////////////////////
TENSOR3 PARTITIONED_SUBSPACE_TET_MESH::interfaceTildeU(int x, int y)
{
  MATRIX& Uij = _interfaceU[x][y];
  TENSOR3 final(Uij.rows(), 3, Uij.cols());

  for (int i = 0; i < Uij.rows() / 3; i++)
  {
    SUBMATRIX subU(Uij, 3 * i, 3);

    TENSOR3 subUTilde = TENSOR3::cross(subU);

    final.setSubtensor(subUTilde, 3 *i);
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// return a tensor of all the cross product tensors produced by the 
// vertices along interface (x,y)
//////////////////////////////////////////////////////////////////////
MATRIX PARTITIONED_SUBSPACE_TET_MESH::interfaceTildeVertices(int x, int y)
{
  vector<pair<int,int> > clones = clonedVertices(x, y);
  vector<VEC3F>& vertices = _meshes[x]->vertices();

  MATRIX final(3 * clones.size(), 3);
  for (unsigned int i = 0; i < clones.size(); i++)
  {
    int index = clones[i].first;
    VEC3F vertex = vertices[index];

    MATRIX vertexCross = MATRIX::cross(vertex);
    final.setSubmatrix(vertexCross, 3 * i);
  }

  return final;
}

//////////////////////////////////////////////////////////////////////
// reduced mass matrix of only the nodes along the (x,y) interface
//////////////////////////////////////////////////////////////////////
MATRIX PARTITIONED_SUBSPACE_TET_MESH::interfaceMassMatrix(int x, int y)
{
  // make a copy
  MATRIX Uij = interfaceU(x,y);
  MATRIX M(Uij.rows(), Uij.rows());

  vector<pair<int,int> > clones = clonedVertices(x, y);

  assert(Uij.rows() == clones.size() * 3);
  for (unsigned int i = 0; i < clones.size(); i++)
  {
    int index = clones[i].first;
    Real mass = _meshes[x]->mass(index);
    M(3 * i, 3 * i) = mass; 
    M(3 * i + 1, 3 * i + 1) = mass; 
    M(3 * i + 2, 3 * i + 2) = mass; 
  }

  return Uij * (M * Uij.transpose());
}

//////////////////////////////////////////////////////////////////////////////
// extract the portion of this vector that is in the submesh
//////////////////////////////////////////////////////////////////////////////
VECTOR PARTITIONED_SUBSPACE_TET_MESH::getSubvector(const int partition, const VECTOR& fullVector)
{
  SUBSPACE_TET_MESH* subMesh= (SUBSPACE_TET_MESH*)_meshes[partition];
  TET_MESH* fullMesh = _originalMesh;

  if (fullVector.size() != fullMesh->x().size())
    fullMesh->x().resizeAndWipe(fullVector.size());

  // set a matrix to the full vector so we can use SUBMATRIX
  MATRIX snapshot(fullMesh->x().size(), 1);
  snapshot.setColumn(fullVector, 0);

  int partitionNodes = subMesh->unconstrainedNodes();
  MATRIX subSnapshotMatrix(partitionNodes * 3, 1);

  // for each vertex in the partition
  for (int y = 0; y < partitionNodes; y++)
  {
    // get the vertex index in the original mesh
    int originalID = this->originalID(partition, y);

    // get the submatrix corresponding to that vertex
    SUBMATRIX vertexBasis(snapshot, originalID * 3, 3);

    // copy the vertex basis into the subbasis
    vertexBasis.copiesInto(subSnapshotMatrix, 3 * y);
  }

  return subSnapshotMatrix.getColumn(0);
}

//////////////////////////////////////////////////////////////////////
// recompute the masses, dividing mass evenly between cloned nodes
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::recomputePartitionedMass()
{
  Real trueMass = _originalMesh->mass(0);
  for (int x = 0; x < _partitions; x++)
  {
    // if no cubature scheme has been computed, don't bother
    if (((SUBSPACE_TET_MESH*)_meshes[x])->totalKeyTets() == 0)
      continue;

    ((SUBSPACE_TET_MESH*)_meshes[x])->resetMasses(trueMass);

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

      massMatrix(index, index) *= invValence;
      massMatrix(index + 1, index + 1) *= invValence;
      massMatrix(index + 2, index + 2) *= invValence;
    }

    if (_meshes[x]->unconstrained())
    {
      UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x];
      if (mesh->U().rows() > 0)
      {
        mesh->cacheSubspaceCenterOfMass();
        mesh->recenterMesh();
      }
    }

    if (((SUBSPACE_TET_MESH*)_meshes[x])->U().rows() > 0)
      _meshes[x]->refreshInertiaTensor();
    _meshes[x]->resetTotalMass();
   
    SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_meshes[x]; 
    mesh->recomputeReducedMass();

    mesh->TET_MESH::refreshInertiaTensor();
  }
}

//////////////////////////////////////////////////////////////////////
// write out state for regression testing
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SUBSPACE_TET_MESH::writeState(const char* filename)
{
  BLOCK_VECTOR state(_partitions);

  for (int x = 0; x < _partitions; x++)
    state.set(*q(x), x);

  VECTOR fullState = state.full();

  fullState.write(filename);
}
