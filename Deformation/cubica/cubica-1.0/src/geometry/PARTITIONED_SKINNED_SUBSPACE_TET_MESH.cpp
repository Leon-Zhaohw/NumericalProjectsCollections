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

#include "PARTITIONED_SKINNED_SUBSPACE_TET_MESH.h"

#if !USING_OSX
#include <GL/glut.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
PARTITIONED_SKINNED_SUBSPACE_TET_MESH::PARTITIONED_SKINNED_SUBSPACE_TET_MESH(SKELETON* skeleton, const char* filename, MATERIAL** materials, int totalMaterials, int partitions, Real springConst, bool simulate, bool loadOriginal, string partitionPath) :
  _skeleton(skeleton)
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

    // force everything to unconstrained -- this is the main deviation from
    // the PARTITIONED_SUBSPACE_TET_MESH constructor   
    cout << " Reading partition file: " << partitionFilename.c_str() << endl;
    bool simulateFullspace = true;
    if (simulate)
      simulateFullspace = false;

    _meshes[x] = new UNCONSTRAINED_SUBSPACE_TET_MESH(partitionFilename.c_str(), materials, totalMaterials, simulateFullspace, NULL, NULL, false);
    _unconstrainedPartition.push_back(true);

    if (simulate)
      ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x])->cacheSubspaceCenterOfMass();
    cout << endl;

    // recompute the mass matrix to make the mesh matrix unity, not just each submesh
    if (simulate)
    {
      Real trueMass = _originalMesh->mass(0);
      ((SUBSPACE_TET_MESH*)_meshes[x])->resetMasses(trueMass);

      // recompute center of mass matrix with new masses
      if (constrainedSize == 0)
        ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x])->cacheSubspaceCenterOfMass();

      // recompute the inertia tensor accordingly
      _meshes[x]->refreshInertiaTensor();
    }

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

  // allocate the clone table
  _clonedVertices = new vector<pair<int, int> >*[_partitions];
  for (int x = 0; x < _partitions; x++)
    _clonedVertices[x] = new vector<pair<int, int> >[_partitions];
  
  // read in the cloned vertices
  string cloneFilename = partitionPath;
  if (cloneFilename.length() == 0)
    cloneFilename = string(filename);
  cloneFilename += string(".clones");
  cout << " Reading clone table: " << cloneFilename.c_str() << endl;
  FILE* file = fopen(cloneFilename.c_str(), "rb");
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
      }
    }
  }
  fclose(file);

  // sanity check
  cout << " Sanity checking clone table ... "; flush(cout);
  for (int x = 0; x < _partitions; x++)
    for (int y = 0; y < _partitions; y++)
      assert(_clonedVertices[x][y].size() == _clonedVertices[y][x].size());
  cout << " done. " << endl;

  // compute the surface vertices
  computeClonedSurfaceVertices();

  // populate the explicit vertex version
  computeClonedVertexMap();

  // compute the center of mass so we can draw the exploded view
  computeCenterOfMass();

  // compute the Us for just the nodes between partitions
  cout << " Computing interface bases ... ";
  flush(cout);
  computeInterfaceUs();
  cout << "done." << endl;

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
  computeInterfaceAreas();
  computeBlendedSurfaceMesh();
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
  _colors[0][0] = 1.0f; _colors[0][1] = 0.5f; _colors[0][2] = 0.5f; _colors[0][3] = 1.0f;
  _colors[1][0] = 0.5f; _colors[1][1] = 1.0f; _colors[1][2] = 0.5f; _colors[1][3] = 1.0f;
  _colors[2][0] = 0.5f; _colors[2][1] = 0.5f; _colors[2][2] = 1.0f; _colors[2][3] = 1.0f;
  _colors[3][0] = 1.0f; _colors[3][1] = 0.5f; _colors[3][2] = 1.0f; _colors[3][3] = 1.0f;
  _colors[4][0] = 0.5f; _colors[4][1] = 1.0f; _colors[4][2] = 1.0f; _colors[4][3] = 1.0f;
  _colors[5][0] = 1.0f; _colors[5][1] = 1.0f; _colors[5][2] = 0.5f; _colors[5][3] = 1.0f;
  _colors[6][0] = 1.0f; _colors[6][1] = 1.0f; _colors[6][2] = 1.0f; _colors[6][3] = 1.0f;
  _colors[7][0] = 0.5f; _colors[7][1] = 0.5f; _colors[7][2] = 0.5f; _colors[7][3] = 1.0f;
}

PARTITIONED_SKINNED_SUBSPACE_TET_MESH::~PARTITIONED_SKINNED_SUBSPACE_TET_MESH()
{
}

//////////////////////////////////////////////////////////////////////
// load a mocap frame into the skeleton, and propagate the bone transforms
// to the partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::loadInterpolatedMocapFrame(float frame)
{
  _skeleton->loadInterpolatedMocapFrame(frame);

  scatterBoneTransforms();
}

//////////////////////////////////////////////////////////////////////
// load a mocap frame into the skeleton, and propagate the bone transforms
// to the partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::loadMocapFrame(int frame)
{
  _skeleton->loadMocapFrame(frame);

  scatterBoneTransforms();
}

//////////////////////////////////////////////////////////////////////
// propagate bone transforms to partitions
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::scatterBoneTransforms()
{
  vector<BONE*> bones = _skeleton->bones();

  for (unsigned int x = 0; x < bones.size() - 1; x++)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x];

    int whichBone = x + 1;
    QUATERNION quaternion = bones[whichBone]->quaternion();
    MATRIX3 rotation = quaternion.toExplicitMatrix3x3();
    VEC3F translation = bones[whichBone]->translation();
    VEC3F& centerOfMass = mesh->originalCenterOfMass();

    VEC3F domainTranslation = rotation * centerOfMass + translation;

    mesh->rigidTranslation() = domainTranslation;
    mesh->rotationQuaternion() = quaternion;
  }
}

//////////////////////////////////////////////////////////////////////
// scatter the constrained node positions to the submeshes -- in the case of a skinned
// animation, these positions will actually change along with the skeleton
//////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::scatterConstrainedNodes()
{
  vector<VEC3F>& originalVertices = _originalMesh->vertices();

  for (int x = 0; x < _partitions; x++)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x];

    // get the rigid components to undo
    VEC3F translation = mesh->rigidTranslation();
    VEC3F centerOfMass = mesh->originalCenterOfMass();
    MATRIX3 rotation = mesh->rotationQuaternion().toExplicitMatrix3x3();
    MATRIX3 RT = rotation.transpose();

    int start = mesh->unconstrainedNodes();
    vector<VEC3F>& vertices = mesh->vertices();

    for (unsigned int y = start; y < vertices.size(); y++)
    {
      int originalID = this->originalID(x,y);
      vertices[y] = RT * (originalVertices[originalID] - translation);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// load skinned simulation data snapshot, including the mocap data
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::loadInterpolatedSkeletonFrame(string dataPath, float frame)
{
  // scatters bone transforms to the partitions as well
  loadInterpolatedMocapFrame(frame);

  // load up just the skinning
  _skeleton->updatePinocchioSkinning(true);

  // does it scatter the cubature points only?
  scatterConstrainedNodes();
}

//////////////////////////////////////////////////////////////////////////////
// scatter the current skeleton bone transforms to the meshes
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::scatterSkeletonTransforms()
{
  // scatter the bone transform data to the meshes -- don't call scatterBoneTransforms(),
  // as that uses Pinocchio's offset indices
  vector<BONE*> bones = _skeleton->bones();
  for (unsigned int x = 0; x < bones.size(); x++)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x];

    QUATERNION quaternion = bones[x]->quaternion();
    MATRIX3 rotation = quaternion.toExplicitMatrix3x3();
    VEC3F translation = bones[x]->translation();
    VEC3F& centerOfMass = mesh->originalCenterOfMass();

    mesh->rigidTranslation() = translation;
    mesh->rotationQuaternion() = quaternion;
  }
}

//////////////////////////////////////////////////////////////////////////////
// load skinned simulation data snapshot, including the mocap data
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::loadOdeSkeletonFrame(string dataPath, int frame)
{
  // if this is the first frame, cache the rigids
  static bool called = false;
  if (frame == 0 && !called) 
  {
    setOdeRestRigids(dataPath);
    called = true;
  }

  // update the frame
  char buffer[256];
  sprintf(buffer, "%04i", frame);
  string skeletonFile = dataPath + string("ode.motion.") + string(buffer) + string(".skeleton");
  _skeleton->loadOdeFrame(skeletonFile.c_str());

  scatterSkeletonTransforms();

  // Don't need to update the constrained nodes because of the way ODE handles
  // the rigid components
  //scatterConstrainedNodes();
}

//////////////////////////////////////////////////////////////////////////////
// load skinned simulation data snapshot, including the mocap data
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::loadSkeletonFrame(string dataPath, int frame)
{
  // scatters bone transforms to the partitions as well
  loadMocapFrame(frame);

  // load up just the skinning
  _skeleton->updatePinocchioSkinning(true);

  // does it scatter the cubature points only?
  scatterConstrainedNodes();
}

//////////////////////////////////////////////////////////////////////////////
// load skinned simulation data snapshot, including the mocap data
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::loadOdeFrame(string dataPath, int frame, VECTOR& position)
{
  // if this is the first frame, cache the rigids
  if (frame == 0)
  {
    vector<BONE*> bones = _skeleton->bones();
    for (unsigned int x = 0; x < bones.size(); x++)
    {
      bones[x]->translationOriginal() = bones[x]->translation();
      bones[x]->quaternionOriginal() = bones[x]->quaternion();
    }
  }

  TET_MESH* originalMesh = _originalMesh;
  originalMesh->x() = position;
  originalMesh->TET_MESH::updateFullMesh();
  originalMesh->readDeformedMesh(frame, dataPath, false);

  // update defo nodes
  vector<BONE*> bones = _skeleton->bones();
  for (int x = 0; x < _partitions; x++)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* subMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x];

    // retrieve the rigid component from the skeleton
    QUATERNION quaternion = bones[x]->quaternion();
    VEC3F translation = bones[x]->translation();

    subMesh->rigidTranslation() = translation;
    subMesh->rotationQuaternion() = quaternion;

    VECTOR restVector = subMesh->restVector();
    VECTOR subPosition = getSubvector(x, position);

    subPosition += restVector;
    subtractTranslation(x, subPosition);
    subtractRotation(x, subPosition);
    subPosition -= restVector;
    subMesh->x() = subPosition;

    subMesh->TET_MESH::updateFullMesh();
  }

  // apply the transform to the constrained nodes as well
  for (int x = 0; x < _partitions; x++)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* subMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x];
    QUATERNION quaternion = bones[x]->quaternion();
    QUATERNION quaternionOriginal = bones[x]->quaternionOriginal();
    VEC3F translation = bones[x]->translation();
    VEC3F translationOriginal = bones[x]->translationOriginal();

    MATRIX3 rotationOriginal = quaternionOriginal.toExplicitMatrix3x3();

    vector<VEC3F>& restPose = subMesh->restPose();
    vector<VEC3F>& vertices = subMesh->vertices();
    int unconstrainedNodes = subMesh->unconstrainedNodes();
    MATRIX3 R = quaternion.toExplicitMatrix3x3();
    MATRIX3 RT = R.transpose();
    VEC3F centerOfMass = subMesh->originalCenterOfMass();

    for (unsigned int y = unconstrainedNodes; y < restPose.size(); y++)
    {
      int original = originalID(x, y);
      vertices[y] = _originalMesh->vertices()[original];

      // subtracting off to local frame would then be
      vertices[y] = (vertices[y] - translation);
      vertices[y] = RT * vertices[y];
    }
    subMesh->TET_MESH::updateFullMesh();
  }
}

//////////////////////////////////////////////////////////////////////////////
// load skinned simulation data snapshot, including the mocap data
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::loadSimulationFrame(string dataPath, int frame, VECTOR& position)
{
  loadMocapFrame(frame);

  TET_MESH* originalMesh = _originalMesh;
  originalMesh->x() = position;
  originalMesh->TET_MESH::updateFullMesh();
  originalMesh->readDeformedMesh(frame, dataPath, false);

  for (int x = 0; x < _partitions; x++)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* subMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x];
    VECTOR restVector = subMesh->restVector();
    VECTOR subPosition = getSubvector(x, position);
    subPosition += restVector;
    subtractTranslation(x, subPosition);
    subtractRotation(x, subPosition);
    subPosition -= restVector;

    subMesh->x() = subPosition;
    subMesh->TET_MESH::updateFullMesh();
  }

  scatterConstrainedNodes();
}

//////////////////////////////////////////////////////////////////////////////
// undo the translation in a snapshot for a given partition
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::subtractTranslation(int partition, VECTOR& snapshot)
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[partition];
  VEC3F translation = mesh->rigidTranslation();

  VEC3F centerOfMass = mesh->originalCenterOfMass();
  translation -= centerOfMass;

  // subtract out the translation
  for (int x = 0; x < snapshot.size() / 3; x++)
  {
    snapshot[3 * x] -= translation[0];
    snapshot[3 * x + 1] -= translation[1];
    snapshot[3 * x + 2] -= translation[2];
  }
}

//////////////////////////////////////////////////////////////////////////////
// undo the rotation in a snapshot for a given partition
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::subtractRotation(int partition, VECTOR& snapshot)
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[partition];
  MATRIX3 rotation = mesh->rotationQuaternion().toExplicitMatrix3x3();
  MATRIX3 transpose = rotation.transpose();

  for (int x = 0; x < snapshot.size() / 3; x++)
  {
    // copy into VEC3s
    VEC3F position3;
    position3[0] = snapshot[3 * x];
    position3[1] = snapshot[3 * x + 1];
    position3[2] = snapshot[3 * x + 2];

    // perform the undo
    position3 = transpose * position3;

    // copy back into the bigger vector
    snapshot[3 * x] = position3[0];
    snapshot[3 * x + 1] = position3[1];
    snapshot[3 * x + 2] = position3[2];
  }
}

//////////////////////////////////////////////////////////////////////////////
// set the rest poses of the meshes according to the first ODE frame
//////////////////////////////////////////////////////////////////////////////
void PARTITIONED_SKINNED_SUBSPACE_TET_MESH::setOdeRestRigids(string dataPath)
{
  string skeletonFile = dataPath + string("ode.motion.0000.skeleton");
  _skeleton->loadOdeFrame(skeletonFile.c_str());

  vector<BONE*> bones = _skeleton->bones();
  for (unsigned int x = 0; x < bones.size(); x++)
  {
    UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_meshes[x];
    vector<VEC3F>& restPose = mesh->restPose();
    vector<VEC3F>& vertices = mesh->vertices();

    QUATERNION quaternion = bones[x]->quaternion();
    MATRIX3 rotation = quaternion.toExplicitMatrix3x3();
    VEC3F translation = bones[x]->translation();
    VEC3F& centerOfMass = mesh->originalCenterOfMass();
    MATRIX3 RT = rotation.transpose();

    for (unsigned int y = 0; y < restPose.size(); y++)
    {
      restPose[y] = RT * ((restPose[y] + centerOfMass) - translation);
      vertices[y] = RT * ((vertices[y] + centerOfMass) - translation);
    }

    mesh->rigidTranslation() = translation;
    mesh->rotationQuaternion() = quaternion;
    mesh->originalCenterOfMass() = translation;

    bones[x]->quaternionOriginal() = quaternion;
    bones[x]->translationOriginal() = translation;
  }
  updateFullMeshes();
}
