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
// The sparse eigensystem solve uses Slepc, which in turn uses PetSc.
// It may not be lightweight, but it is robust, well-documented, and
// relatively easy to build.
//
//////////////////////////////////////////////////////////////////////

#ifndef TET_MESH_H
#define TET_MESH_H

#include <SETTINGS.h>
#include <TET.h>
#include <TRIANGLE.h>
#include <VEC3.h>
#include <SPARSE_MATRIX.h>
#include <BLOCK_MATRIX.h>
#include <MATERIAL.h>
#include <MATRIX.h>
#include <VECTOR.h>
#include <cstdio>
#include <map>
#include <vector>
#include <algorithm>
#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#ifdef USING_RENDERMAN
#include <ri.h>
#endif
#include <OBJ.h>
#ifndef IGNORE_PETSC
#include <SPARSE_PETSC_MATRIX.h>
#endif
#include <SPARSE_PCG_MATRIX.h>
#include <TIMER.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Tetrahedron mesh class
//
// Note these use a LEFT handed coordinate system (like in screenspace)
//////////////////////////////////////////////////////////////////////
class TET_MESH {

public:
  // "filename" is the tet mesh file
  // "material" is the constitutive model
  // "simulate" decides if we want to allocate the gigantic matrices for
  //            direct simulation
  // "volumeScale" decides whether or not we want to scale tet stiffnesses,
  //               forces, etc. by volume
  TET_MESH(const char* filename, MATERIAL** materials, int totalMaterials,
           bool simulate = false, Real scale = 1.0,
           bool volumeScale = false);

  // one-ring tet mesh constructor --
  // 
  // this creates a copy of the one-ring of the passed in vertex and TET_MESH,
  // constraining all other vertices aside from the one in the middle
  //
  // originalToCopy returns the correspondance between the vertices passed
  // in and the ones in the new mesh
  TET_MESH(VEC3F* point, TET_MESH* mesh, map<VEC3F*, int>& originalToCopy,
           map<string, double>& timingBreakdown);
  virtual ~TET_MESH(); 

  // bunch of drawing options
  virtual void drawAllTets();
  virtual void drawZSlice(Real z);
  virtual void drawConstrainedNodes();
  virtual void drawUnconstrainedNodes();
  virtual void drawCollisionNodes();
  virtual void drawSurfaceFaces();
  virtual void drawSurfaceVertices();
  virtual void drawPartition(int partition);
  virtual void drawColoredPartitions();
  virtual void drawTetMaterials(Real z = 0.5);
  virtual void drawRestPose();
  void drawFirstPK();
  void drawPoints(vector<VEC3F*>& points);
  void drawVectorField(VECTOR field);

  void drawSurfaceToRenderMan();
  void drawEmbeddingToRenderMan();
  void drawHeadEmbeddingToRenderMan(bool bakingSurface = false,
                                    bool subsurface = false,
                                    string filename = string(""));
  void drawExhaustiveToRenderMan();
  void drawColoredToRenderMan();
  void drawOutlinesToRenderMan();
  void drawSurfaceForcesToRenderMan();
  void drawZSliceToRenderMan(float zSlice = 0.5);
  
  // mesh degrees of freedom
  virtual int rank() { return _unconstrainedSize * 3; };
  
  // is any tet in the mesh inverted?
  bool inverted(bool verbose = false);

  // reset all vertex positions to rest pose
  void resetToRestPose();
 
  // reinitialize the mesh (after a scaling for example)
  void reinitialize();

  // accessors
  MATERIAL** materials()       { return _materials; };
  MATERIAL*** materialCopies() { return _materialCopies; };
  int& totalMaterials()          { return _totalMaterials; };
  int totalMaterials() const     { return _totalMaterials; };
  vector<TET>& tets()            { return _tets; };
  int totalTets() const          { return (int)(_tets.size()); };
  int unconstrainedNodes() const { return _unconstrainedSize; };
  int constrainedNodes() const   { return _constrainedSize; };
  int dofs() const               { return 3 * _unconstrainedSize; };
  int totalNodes() const         { return (int)(_vertices.size()); };
  string filename()        { return _filename; };
  VECTOR& x()              { return _x; };
  SPARSE_MATRIX& massMatrix()      { return _masses; };
  SPARSE_MATRIX& stiffnessMatrix() { return _stiffness; };
  int vertexID(VEC3F* vertex)      { return (_vertexID.find(vertex) != _vertexID.end()) ? _vertexID[vertex] : -1; };

  int tetID(TET* tet)              { return _tetID[tet]; };
  vector<VEC3F>& vertices()        { return _vertices; };
  vector<VEC3F>& restPose()        { return _restPose; };
  VEC3F* vertices(int index)       { return &_vertices[index]; };
  VEC3F* restVertices(int index)   { return &_restPose[index]; };
  Real mass(int vertex)            { return _masses(3 * vertex, 3 * vertex); };
  VEC3F& centerOfMass()             { return _centerOfMass; };
  VEC3F& restCenterOfMass()         { return _restCenterOfMass; };
  vector<int>& tetMembership(VEC3F* vertex) { return _tetMembership[vertex]; };
  vector<pair<int,int> >& surfaceFaces() { return _surfaceFaces; };
  vector<TRIANGLE*>& explicitSurfaceFaces() { return _explicitSurfaceFaces; };
  vector<TRIANGLE*>& allFaces() { return _allFaces; };
  vector<VEC3F*>& surfaceVertices() { return _surfaceVertices; };
  vector<VEC3F*>& collisionNodes()  { return _collisionNodes; };
  virtual VECTOR& F()               { return _F; };
  virtual VECTOR& internalForce()   { return _R; };
  vector<VEC3F>& forceVectors()     { return _internalForces; };
  vector<int>& tetEmbeddings()      { return _tetEmbeddings; };
  const MATRIX& inertiaTensor()     { return _inertiaTensor; };
  const MATRIX& inertiaTensorOriginal() { return _inertiaTensorOriginal; };
  MATRIX& inertiaTensorOld()        { return _inertiaTensorOld; };
  const MATRIX& inertiaTensorDt()   { return _inertiaTensorDt; };
  MATRIX& inertiaTensorDtOld()      { return _inertiaTensorDtOld; };
  Real surfaceArea(int vertexID)    { return surfaceArea(&(_vertices[vertexID])); };
  Real surfaceArea(VEC3F* vertex)   { return (_surfaceArea.find(vertex) == _surfaceArea.end()) ? 0.0 : _surfaceArea[vertex]; };
  vector<VEC3F>& barycentricEmbeddings() { return _barycentricEmbeddings; };
  Real& meshScaling()               { return _meshScaling; };
  VEC3F surfaceNormal(VEC3F* vertex);

  // queries
  bool constrained()                { return _constrainedSize != 0; };
  bool unconstrained()              { return _constrainedSize == 0; };
  bool isOnSurface(VEC3F* vertex)   { return _vertexOnSurface[vertex]; };
  bool isConstrained(int vertexID)  { return vertexID >= _unconstrainedSize; };
  bool isConstrained(VEC3F* vertex) { return _vertexID[vertex] >= _unconstrainedSize; };

  // lookup rest vertex version of this vertex
  VEC3F* restVersion(VEC3F* vertex) { return &(_restPose[_vertexID[vertex]]); };

  // find closest node to an arbitrary point
  VEC3F* closestSurfaceNode(VEC3F point);

  // Find all surface nodes within a given radius
  void closestSurfaceNodes( VEC3F point, Real radius,
                            std::vector<VEC3F *> &surfaceNodes );

  // Newmark Integrator support functions
  virtual void generateF();
  virtual VECTOR& generateInternalForces();   //< assumes generateF() was already called!
  virtual MATRIX& generateStiffnessMatrix();  //< assumes generateF() was already called!
  SPARSE_MATRIX& generateSparseStiffnessMatrix(bool verbose = false);
  void generateSparseStiffnessMatrix(SPARSE_MATRIX& stiffness);
#ifndef IGNORE_PETSC
  void generateSparseStiffnessMatrix(SPARSE_PETSC_MATRIX& stiffness);
#endif
  SPARSE_MATRIX& generateSparseStiffnessMatrix(vector<MATRIX3>& Us, 
                                               vector<MATRIX3>& Vs, 
                                               vector<MATRIX>& stiffnesses);
  SPARSE_MATRIX& generateOneRingStiffness(vector<MATRIX3>& Us, 
                                               vector<MATRIX3>& Vs, 
                                               vector<MATRIX>& stiffnesses);
#ifndef IGNORE_PETSC
#if _WIN32  
  void generateSparseStiffnessMatrix(vector<MATRIX3>& Us, 
                                     vector<MATRIX3>& Vs, 
                                     vector<MATRIX>& stiffnesses, 
                                     SPARSE_PCG_MATRIX& stiffness);
#else
  void generateSparseStiffnessMatrix(vector<MATRIX3>& Us, 
                                     vector<MATRIX3>& Vs, 
                                     vector<MATRIX>& stiffnesses, 
                                     SPARSE_PETSC_MATRIX& stiffness);
#endif
#endif
  VECTOR& generateInternalForces(vector<MATRIX3>& Us,
                                 vector<MATRIX3>& Fhats,
                                 vector<MATRIX3>& Vs);

  // Peek the eigenvalues of each tet and return the smallest found
  Real inspectStiffnessEigenvalues();

  // compute the Hessian product with respect to two vectors --
  // H * first * second
  VECTOR hessianProduct(VECTOR& first, VECTOR& second);
  
  virtual void updateSurfaceMesh();
  virtual void updateFullMesh();
  virtual void updateFullMesh(VECTOR& position);

  // write the mesh out to a METIS-friendly partitioning format
  void writeMETIS(const char* filename);

  // read in the results of a METIS partitioning and tag the 
  // 'partition' field of the tets
  void readMETIS(const char* filename);

  // write the mesh out to a Scotch-friendly partitioning format
  void writeScotch(const char* filename);

  // read in the results of a Scotch partitioning and tag the 
  // 'partition' field of the tets
  void readScotch(const char* filename);

  // write out a partitioned version of this mesh
  void writePartitions(string prefix = string(""));
  
  // write the stiffness matrix to matlab
  void stiffnessToMatlab(const char* filename = "stiffness.m");
  
  // populate _centerOfMass
  void computeCenterOfMass();
  void computeRestCenterOfMass();

  // return the one-ring of a vertex
  void oneRing(VEC3F* vertex, vector<VEC3F*>& result);
  void oneRing(int tetID, vector<VEC3F*>& result);
  void oneRing(VEC3F* vertex, vector<TET*>& result);

  // return the one-ring of a tet
  void tetOneRing(int tetIndex, vector<int>& oneRing);
  void vertexOneRing(VEC3F* vertex, vector<int>& oneRing);

  // write out the deformed mesh
  void writeDeformedMesh(int timestep, string path);
  bool readDeformedMesh(int timestep, string path, bool verbose = true);

  // animate the nodes of the mesh somehow
  void animate(int timestep);

  // displacements of the constrained nodes
  VECTOR constraintDisplacements();
  
  // read in tet mesh materials if there is more than one type
  void readMaterials();
  
  // write out tet mesh materials if there is more than one type
  void writeMaterials();

  // get the bounding box
  void boundingBox(Real& xMin, Real& xMax, Real& yMin, Real& yMax, Real& zMin, Real& zMax);
  void boundingBox(Real* bbox) { boundingBox(bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]); };
  void boundingBox(VEC3F& mins, VEC3F& maxs);

  // compute collision nodes
  void computeCollisionNodes();

  // constrain nodes inside this surface
  void writeNewConstraints(vector<VEC3F*>& newConstraints, const char* filename);

  // embed this surface into the mesh
  void embedMesh(OBJ* obj);

  // compute the embedding for this mesh
  void normalizeEmbedding(OBJ* obj);
  void computeEmbedding(OBJ* obj);

  // update the current embedded mesh
  void updateEmbedding();

  // update embedded mesh normals
  void updateEmbeddingNormals();

  // update the current embedded mesh
  void updateHeadEmbedding();

  // scale the mass of the mesh
  void scaleMasses(Real scale);

  // get the mean stretch ratio of the mesh
  Real meanStretchRatio();

  // write out the embedding
  void writeEmbedding();

  // get some tet stats
  Real minTetVolume();
  Real maxTetVolume();

  // find the number of fully connected vertices
  int testFullyConnected();

  // recompute x by subtracting the restVertices for vertices
  void recoverX();
  
  // write out a tet mesh containing just a single tet (for debugging purposes)
  static void writeSingleTetMesh(const char* filename);

  // check that the center of mass is in fact centered at zero
  virtual bool centerOfMassIsZero();

  // reset mesh so that the center of mass is at the origin
  void centerOfMassToOrigin();

  // return a vector of all the rest positions
  VECTOR restVector();

  // return a block identity matrix with the same rank as this mesh
  SPARSE_MATRIX sparseIdentity();
  BLOCK_MATRIX blockIdentity();

  // return the total mass in the mesh (div 3 because the mass is repeated
  // three times per vertex)
  //Real totalMass() { return _masses.sum() / 3; };
  Real totalMass() { return _totalMass; };

  int totalCores() const { return _totalCores; }

  // refresh the inertia tensor
  virtual const MATRIX& refreshInertiaTensor();

  // refresh the inertia tensor time derivative
  virtual const MATRIX& refreshInertiaTensorDt(VECTOR& velocity);

  // reset the mass matrix to a different weight
  virtual void resetMasses(Real mass);

  // reset the total mass
  void resetTotalMass() { _totalMass = _masses.sum() / 3; };

  // recache the original inertia tensor
  void recacheInertiaTensorOriginal() { _inertiaTensorOriginal = _inertiaTensor; };

  // recompute current center of mass
  VEC3F SitBar();
  SPARSE_MATRIX SiBar();

  // stream in/out a mesh file to an already open file
  void writeMeshFile(FILE* file);
  void readMeshFile(FILE* file, Real scale = 1.0, bool verbose = true);

  // compute the shape matching transformation
  // make sure to call updateFullMesh() before this!
  void computeShapeMatching(VEC3F& translation, MATRIX3& rotation) const;

  // compute the deformation minus any rigid shape matching component
  VECTOR computeDefoMinusShapeMatching() const;

  // subtract the rigid components from a field
  void subtractRigidComponent(const VECTOR& position, const VEC3F& rigidCenter, VECTOR& field, VEC3F& linearVector, VEC3F& angularVector);

  // compute geometric center (i.e. non mass-weighted):
  VEC3F computeGeometricCenter() const;
  VEC3F computeGeometricRestCenter() const;

  // compute the center that is being used to compute any angular fields
  VEC3F computeAngularCenter();

protected:
  vector<VEC3F> _vertices;
  vector<VEC3F> _restPose;
  vector<TET> _tets;
  vector<VEC3F> _internalForces;
  string _filename;

  // vertex masses -- _masses.entry(i,i) = mass of _vertices[i]
  SPARSE_MATRIX _masses;

  // stiffness matrix
  SPARSE_MATRIX _stiffness;
  
  // list of surface faces -- the int pair is <tetID, faceIndex>
  vector<pair<int,int> > _surfaceFaces;
  vector<VEC3F*> _surfaceVertices;
  map<VEC3F*, bool> _vertexOnSurface;
  vector<VEC3F> _surfaceNormals;
  map<VEC3F*, int> _surfaceVertexID;
  
  // list of surface faces
  vector<TRIANGLE*> _explicitSurfaceFaces;

  // list of all faces
  vector<TRIANGLE*> _allFaces;

  // _vertices[_unconstrainedSize] to 
  // _vertices[_unconstrainedSize + _constrainedSize - 1], 
  // are constrained
  int _constrainedSize;

  // _vertices[0] to _vertices[_unconstrainedSize - 1]
  // are unconstrained
  int _unconstrainedSize;

  // tracks if a vertex is constrained
  map<VEC3F*, bool> _constrained;

  // maps vertex address to its index in _vertices
  map<VEC3F*, int> _vertexID;

  // maps tet address to its index in _tets
  map<TET*, int> _tetID;

  // track which tets a vertex is a member of
  map<VEC3F*, vector<int> > _tetMembership;

  // track the surface area associated with each vertex
  map<VEC3F*, Real> _surfaceArea;

  // constitutive model
  //MATERIAL* _material;
  MATERIAL** _materials;
  int _totalMaterials;

  // current center of mass
  VEC3F _centerOfMass;
  VEC3F _restCenterOfMass;

  // sub-dermal nodes to check for collision
  vector<VEC3F*> _collisionNodes;

  // inertia tensor of mesh
  MATRIX _inertiaTensor;
  MATRIX _inertiaTensorOld;
  MATRIX _inertiaTensorDt;
  MATRIX _inertiaTensorDtOld;
  MATRIX _inertiaTensorOriginal;

  Real _totalMass;

  Real _meshScaling;

  // initialize the mesh after the vertices have been specified
  void init(bool simulate, bool verbose = false);

  // read in the mesh file
  void readMeshFile(const char* filename, Real scale = 1.0);

  // write out a mesh file
  void writeMeshFile(const char* filename);

  // read in the surface faces
  bool readSurfaceFaceCache();

  // compute which faces are on the surface of the tet mesh
  void computeSurfaceFaces();
  void computeExplicitSurfaceFaces();
  void computeAllFaces();
  void computeSurfaceVertices();

  // compute surface area of the vertices on the surface
  void computeSurfaceAreas();

  // compute the masses - currently just makes all the vertices sum
  // to unity
  void computeMasses();

  // given a deformed vertex, return the corrsponding rest pose vertex
  // this assumes the rest pose vertices are in the same array positions
  // as the tet vertices
  VEC3F& restVertex(VEC3F* vertex) { return _restPose[_vertexID[vertex]]; };

  // internal force functions
  void clearInternalForces();
  void computeInternalForces();

  // helper function for oneRing
  VEC3F* outlierVertex(vector<int>* membership, int currentTet, int v0, int v1, int v2);

  // helper function for writePartitions -- writes out the file
  // that describes which vertices are cloned across partitions
  void writeClones(int totalPartitions, map<VEC3F*, int>* newVertexIDs, string prefix);
 
  // helper function for writePartitions -- writes out the file
  // that translates the partition vertex index to an index in the
  // original unpartitioned mesh
  void writeVertexCorrespondence(int totalPartitions, map<VEC3F*, int>* newVertexIDs, string prefix);
  
  // compute area weighted normals of the surface vertices
  void computeSurfaceNormals();

  ///////////////////////////////////////////////////////////////////
  // Newmark variables
  ///////////////////////////////////////////////////////////////////
  
  // current deformation state
  VECTOR _x;

  // force vector for Newmark
  VECTOR _R;

  // stiffness matrix for Newmark
  MATRIX _K;

  // stacked deformation gradients, so that they can all be computed
  // in one big matrix-vector multiply
  VECTOR _F;

  // DEBUG ONLY!
  VECTOR _debug;

  ///////////////////////////////////////////////////////////////////
  // OpenMP support variables -- per thread versions of existings
  // variables
  ///////////////////////////////////////////////////////////////////
  int _totalCores;
  MATERIAL*** _materialCopies;
  SPARSE_MATRIX* _stiffnessCopies;
  VECTOR* _RCopies;

  ///////////////////////////////////////////////////////////////////
  // Embedded mesh variables
  ///////////////////////////////////////////////////////////////////
  OBJ* _embeddedMesh;
  vector<int> _tetEmbeddings;
  vector<VEC3F> _barycentricEmbeddings;

  void computeEmbedding();
  bool readEmbedding();

  // Lookup rest version of a defomed vertex
  map<VEC3F*, VEC3F*> _deformedToRest;

  ///////////////////////////////////////////////////////////////////
  // Scales for each key tet.  These store the tet volumes if
  // we are scaling by volume, or 1 otherwise.
  ///////////////////////////////////////////////////////////////////
  vector<Real> _tetScales;

  ///////////////////////////////////////////////////////////////////
  // RenderMan support function
  ///////////////////////////////////////////////////////////////////
  void rainbowRamp(Real input, float color[3]);
  
  // compute the angular vector induced by a vector field
  VEC3F computeAngularVector(const VECTOR& field, const VEC3F& center, bool verbose = false);
  
  // compute the angular vector induced by a vector field
  VECTOR computeAngularField(const VEC3F& angular, const VEC3F& center);
};

#endif
