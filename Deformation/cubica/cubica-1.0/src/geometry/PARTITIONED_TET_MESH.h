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
// PARTITIONED_TET_MESH.h: interface for the PARTITIONED_TET_MESH class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PARTITIONED_TET_MESH_H
#define PARTITIONED_TET_MESH_H

#include "TET_MESH.h"
#include "UNCONSTRAINED_TET_MESH.h"
#include <BLOCK_SPARSE_MATRIX.h>

class PARTITIONED_TET_MESH {
public:
  PARTITIONED_TET_MESH(const char* filename, MATERIAL** materials, int totalMaterials, int partitions, Real springConst = 2.0, bool simulate = false, string partitionPath = string(""));
  PARTITIONED_TET_MESH();
  virtual ~PARTITIONED_TET_MESH();

  // draw partitions by color
  void drawSurfaceFaces(bool shadowCall = false);
  void drawHighlightedPartition(int highlight);
  void drawSurfaceFaces(int partition);
  void drawExploded();
  void drawExploded(int left, int right);
  void drawExplodedCentered(int center);
  void drawSprings();
  void drawSprings(int partition);
  void drawClonedTriangles();
  void drawClonedTriangles(int partition);
  void drawClonedVertices(int partition);
  void drawClonedSurfaceVertices(int partition);
  virtual void drawBlendedMesh();
  void drawCloneConnections();
  void drawSurfaceToRenderMan();
  void drawEmbeddingToRenderMan();
  void drawConstrainedNodes();
  void drawClonesOnly(int left, int right);
  virtual void drawRigidFrames();
  void drawEmbedding();
  void drawFilteredEmbedding();
  void drawFilteredEmbeddingSubset();

  void drawAllTets();

  // draw cloned vertices
  void drawClone(int partition0, int partition1, int clone);
  
  // find closest node to an arbitrary point
  virtual int closestPartition(VEC3F point);

  // update all meshes
  void updateFullMeshes();
  
  // accessors
  TET_MESH* originalMesh() { return _originalMesh; };
  TET_MESH* mesh(int x)    { return _meshes[x]; };
  TET_MESH** meshes()      { return _meshes; };
  int partitions()         { return _partitions; };
  int originalID(int partition, int nodeID) { return _originalIDs[partition][nodeID]; };
  string filename()        { return _filename; };
  string partitionPath()   { return _partitionPath; };
  int* graphColors()       { return _graphColors; };
  int& totalGraphColors()  { return _totalGraphColors; };
  bool unconstrained(int partition) { return _unconstrainedPartition[partition]; };
  bool constrained(int partition) { return !_unconstrainedPartition[partition]; };
  vector<pair<int, int> >& tetEmbeddings() { return _tetEmbeddings; };
  vector<VEC3F>& barycentricEmbeddings() { return _barycentricEmbeddings; };
  vector<TRIANGLE*>& blendedSurfaceMesh() { return _blendedSurfaceMesh; };
  int totalUnconstrained();
  void getUnconstrainedIDs(vector<int>& unconstrained);
  VEC3F centerOfMass(int partition) { return _meshes[partition]->centerOfMass(); };
  BLOCK_SPARSE_MATRIX& diagonalInterfaceStiffness() { return _diagonalInterfaceStiffness; };
  BLOCK_SPARSE_MATRIX& offDiagonalInterfaceStiffness() { return _offDiagonalInterfaceStiffness; };
  virtual QUATERNION quaternionRotation(int partition);
  virtual QUATERNION quaternionRotationOld(int partition);
  virtual Real rotationLambda(int partition);
  virtual VEC3F rigidTranslation(int partition);
  int dofs(int partition) { return _meshes[partition]->dofs(); };
  Real interfaceArea(int x, int y) { return _interfaceArea[x][y]; };
  int totalClones(int x, int y) { return _clonedVertices[x][y].size(); };
  Real interfaceMass(int x, int y);
  int rank(int partition) { return _meshes[partition]->rank(); };
  Real totalMass(int partition) { return _meshes[partition]->totalMass(); };
 
  // get a vector of <partition, nodeID in partition> pairs that correspond to the
  // vertexID in the original mesh
  vector<pair<int,int> > partitionIDs(int vertexID);

  // are partitions x and y neighbors?
  //bool neighbors(int x, int y) { return _clonedVertices[x][y].size() != 0; };
  //bool neighbors(int x, int y) { return _interfaceArea[x][y] > 0.0; };
  bool neighbors(int x, int y);
  
  // do partitions x and y share faces?
  int sharedFaces(int x, int y) { return _clonedTriangles[x][y].size(); };

  // return the cloned vertex list
  vector<pair<int,int> >& clonedVertices(int x, int y) { return _clonedVertices[x][y]; };

  // is this vertex cloned?
  bool isCloned(VEC3F* vertex) { return _isCloned.find(vertex) == _isCloned.end() ? false : true; };

  // how many unique clone pairs are there?
  int totalClones();

  // embedding support functions
  void normalizeEmbedding(OBJ* obj);
  void computeEmbedding(OBJ* obj);
  void writeEmbedding(string filename = string(""));
  bool readEmbedding(string filename = string(""));
  virtual void updateEmbedding();

  // embed this surface into the mesh
  void embedMesh(OBJ* obj);

  // compute the full sparse interface stiffness matrix
  void computeFullInterfaceStiffness(BLOCK_SPARSE_MATRIX& springMatrix);
  void computeDiagonalInterfaceStiffness(BLOCK_SPARSE_MATRIX& diagonal);
  void computeOffDiagonalInterfaceStiffness(BLOCK_SPARSE_MATRIX& offDiagonal);

  // get the subvector corresponding to this submesh
  VECTOR getSubvector(const int partition, const VECTOR& fullVector);
  VECTOR getSubposition(const int partition) { return getSubvector(partition, _originalMesh->x()); };

  // compute the center of mass over all partitions
  virtual void updateCentersOfMass();

  // gather/scatter all the partition xs into _superX
  void scatterSuperX(BLOCK_VECTOR& superX);
  void gatherSuperX(BLOCK_VECTOR& superX);

  // compute interface stiffness with rotations
  void updateBlockInterfaceStiffness();

  Real springConst(int x, int y) { return _interfaceSpringConst(x,y) * _interfaceArea[x][y] / _maxInterfaceArea; };

  virtual MATRIX3 rigidRotation(int partition);

  // assemble a vector of all the cloned rest vertices along the (x,y) partition boundary
  VECTOR interfaceRestVertices(int x, int y);

  // get the original tet index for tet 'tetID' in partition 'partition'
  int originalTetID(int partition, int tetID);
 
  // write embedded file out to pbrt 
  void exportEmbeddingToPBRT(string renderPath, int frame);
  void exportFilteredEmbeddingToPBRT(string renderPath, int frame);

  void writeFilteredOBJ(string renderPath, int frame);

  // stream read/write
  void readState(FILE* file);
  void writeState(FILE* file);

protected:
  // total number of partitions
  int _partitions;
  TET_MESH** _meshes;
  TET_MESH* _originalMesh;
  string _filename;
  string _partitionPath;

  // center of mass
  VEC3F _centerOfMass;
  
  // the list of cloned vertices between the different
  // partitions. The declaration is ugly but the
  // interpretation is straightforward
  //
  // _clonedVertices[0][1] is the list of vertices that are cloned
  // between partition 0 and partition 1
  // _clonedVertices[0][1].first is the index of the vertex in partition 0
  // _clonedVertices[0][1].second is the index of the same vertex in partition 1
  vector<pair<int, int> >** _clonedVertices;

  // _clonedVertexMap[0][1] is a map of vertices that are cloned
  // between partition 0 and partition 1
  //
  // Performing a lookup for some VEC3F* vertex in partition 1
  // _clonedVertexMap[0][1][vertex] will return the address of the clone
  // in partition 0 (note the indices are reversed!)
  map<VEC3F*, VEC3F*>** _clonedVertexMap;

  // cloned vertex areas - fraction of the interface area represeneted
  // by each cloned vertex. Used to fine-tune the spring constant for
  // the interface stiffness matrix
  map<VEC3F*, Real> _clonedVertexAreas;

  // maximum cloned vertex area
  Real _maxClonedVertexArea;

  // subset of _clonedVertices, where these are the ones that are
  // on the surface of the original mesh
  vector<pair<int, int> >** _clonedSurfaceVertices;
  
  // lists what the index of the vertices in the new partition
  // were in the original mesh -- this is necessary to carve
  // out the portion of the U matrix that corresponds to the partition
  //
  // For example, _originalID[x][y] is the index in _vertices of the
  // y-th vertex in partition x
  vector<int>* _originalIDs;

  // the inverse of _originalIDs -- given some vertex in the original mesh,
  // what is the vertex in the partitioned mesh?
  multimap<int, pair<int, int> > _partitionIDs;
 
  // which triangles are along partition interfaces
  //
  // _interfaceTriangles[0][1] is the list of triangles in partition
  // 0 that are also cloned in partition 1
  // between partitions 0 and 1
  vector<TRIANGLE*>** _clonedTriangles;

  // indices of the triangles in _clonedTriangles
  vector<int>** _clonedTriangleIDs;

  // Interface surface areas
  Real** _interfaceArea;

  // maximum of all the interface areas
  Real _maxInterfaceArea;

  // the surface triangles, with cloned surface vertex positions
  // blended between partitions
  vector<TRIANGLE*> _blendedSurfaceMesh;

  // vertices of the surface mesh that are blended --
  // the VEC3F pointers in the triangles should point here.
  vector<VEC3F*> _blendedVertices;

  // lists of partitions to blend between for each entry of _blendedVertices
  // -- should be the same size as _blendedVertices, and each entry should
  // have a vector of at least size 2.
  vector<vector<VEC3F*> > _verticesToBlend;
 
  // storage location of triangles with blended vertices. *Only* the triangles
  // with blended vertices are stored here, the other triangles are just the
  // ones from the underlying TET_MESH* objects.
  vector<TRIANGLE*> _blendedTriangles;
 
  // is a partition an unconstrained mesh?
  vector<bool> _unconstrainedPartition;
  
  // colors for the different partitions
  float _colors[8][4];

  // graph coloring returned by DSATUR
  int* _graphColors;

  // total colors in the graph coloring
  int _totalGraphColors;

  // track which vertices are cloned
  map<VEC3F*, int> _isCloned;

  // embedded mesh variables --
  // the pairs are <partition #, tet #> 
  vector<pair<int, int> > _tetEmbeddings;
  vector<VEC3F> _barycentricEmbeddings;
  OBJ* _embeddedMesh;

  // split into diagonal and off-diagonal parts
  BLOCK_SPARSE_MATRIX _diagonalInterfaceStiffness;
  BLOCK_SPARSE_MATRIX _offDiagonalInterfaceStiffness;

  // spring const in interface stiffness
  Real _springConst;

  // make the spring const a per-interface constant
  MATRIX _interfaceSpringConst;
  
  // compute the center of mass over all partitions
  void computeCenterOfMass();

  // compute triangles along the interface
  void computeClonedTriangles();

  // compute the area of all the interfaces
  void computeInterfaceAreas();
  
  // compute the area of the interface represented by each cloned vertex
  void computeClonedVertexAreas();

  // compute the cloned surface vertex list
  void computeClonedSurfaceVertices();

  // compute the blended surface mesh
  void computeBlendedSurfaceMesh();

  // write the partition graph out for coloring
  void writeDIMACS(const char* filename);

  // read a graph coloring in
  bool readDIMACS(const char* filename);

  // compute a graph coloring for the partitions
  void computeGraphColoring(string prefix);

  // populate _clonedVertexMap
  void computeClonedVertexMap();

  // read/write cloned triangles
  bool readClonedTriangles();
  void writeClonedTriangles();
  
  // recompute the masses, dividing mass evenly between cloned nodes
  virtual void recomputePartitionedMass();
};

#endif
