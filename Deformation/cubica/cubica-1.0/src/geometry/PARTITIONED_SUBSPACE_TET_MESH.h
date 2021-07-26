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
// PARTITIONED_SUBSPACE_TET_MESH.h: interface for the PARTITIONED_SUBSPACE_TET_MESH class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PARTITIONED_SUBSPACE_TET_MESH_H
#define PARTITIONED_SUBSPACE_TET_MESH_H

#include "PARTITIONED_TET_MESH.h"
#include "SUBSPACE_TET_MESH.h"
#include "UNCONSTRAINED_SUBSPACE_TET_MESH.h"
#include <BLOCK_MATRIX.h>
#include <SANDWICH_TRANSFORM.h>
#include <TENSOR3.h>

class PARTITIONED_SUBSPACE_TET_MESH : public PARTITIONED_TET_MESH {

public:
  PARTITIONED_SUBSPACE_TET_MESH(const char* filename, MATERIAL** materials, int totalMaterials, int partitions, Real springConst = 2.0, bool simulate = false, bool loadOriginal = true, string partitionPath = string(""));
  PARTITIONED_SUBSPACE_TET_MESH();
  virtual ~PARTITIONED_SUBSPACE_TET_MESH();

  SPARSE_MATRIX& interfaceStiffness() { return _interfaceStiffness; };
  BLOCK_MATRIX& diagonalInterfaceStiffness() { return _diagonalInterfaceStiffness; };
  BLOCK_MATRIX& offDiagonalInterfaceStiffness() { return _offDiagonalInterfaceStiffness; };
  VECTOR& superQ() { return _superQ; };
  int rank(int partition) const { return _ranks[partition]; };
  int totalInterfaces() const { return _totalInterfaces; };
  int superRank() const { return _superQ.size(); };
  
  // distribute the _UBasis of _originalMesh to the partitions
  void distributeSubbases();

  // read in the cubature for all the partitions
  void readCubature();
  
  // gather/scatter all the partition qs into _superQ
  void gatherSuperQ();
  void gatherSuperQ(BLOCK_VECTOR& superQ);
  void gatherSuperQOld(BLOCK_VECTOR& superQ);
  void scatterSuperQ();
  void scatterSuperQ(BLOCK_VECTOR& superQ);
  void scatterSuperQOld(BLOCK_VECTOR& superQ);

  // get the q for a specific partition
  VECTOR* q(int partition) { return &((SUBSPACE_TET_MESH*)_meshes[partition])->q(); };
  VECTOR* qOld(int partition) { return &((SUBSPACE_TET_MESH*)_meshes[partition])->qOld(); };
  
  // get the U for a specific partition
  MATRIX& U(int partition) { return ((SUBSPACE_TET_MESH*)_meshes[partition])->U(); };
  
  // get the M for a specific partition
  SPARSE_MATRIX& massMatrix(int partition) { return ((SUBSPACE_TET_MESH*)_meshes[partition])->massMatrix(); };
  
  // was the basis of the original mesh loaded?
  bool basisLoaded() { return ((SUBSPACE_TET_MESH*)_originalMesh)->basisLoaded(); };

  int totalKeyTets(int partition) { return((SUBSPACE_TET_MESH*)_meshes[partition])->totalKeyTets(); };
  VECTOR& keyWeights(int partition) { return((SUBSPACE_TET_MESH*)_meshes[partition])->keyWeights(); };
  
  // interface spring constant
  //Real& interfaceSpringConst() { return _springConst; };
 
  // draw partitions with only rigid transform -- no deformation
  void drawRigidOnly();
  
  // draw partitions with debug rigid transform -- no deformation
  void drawRigidDebug();

  // draw coordinate axes for the rigid frames
  virtual void drawRigidFrames();

  // draw the centers of mass for the unconstrained meshes
  void drawCentersOfMass();

  // draw the surface mesh with blended vertices
  virtual void drawBlendedMesh();
 
  // draw a point for the mean interface position
  void drawInterfaceMeans();
  VEC3F computeInterfaceMean(int x, int y);
  VEC3F computeInterfaceMeanRest(int x, int y);

  // combine the vertex positions in all the partitions and 
  // write out to Matlab format
  virtual void verticesToMatlab(const char* filename, const char* varName, int column = 0);
  virtual void verticesToBinary(const char* filename, int column = 0);

  // DEBUG only
  void verifyInterfaceBasis(int index);

  void resetRigidTranslations();

  // return whether this partition is constrained or not
  bool unconstrained(int partition) { return _unconstrainedPartition[partition]; };

  VECTOR& restDiff(int x, int y) { return _restDiffs[x][y]; };
  MATRIX& interfaceU(int x, int y) { return _interfaceU[x][y]; };
  bool interfaceExists(int x, int y) const { return (_interfaceU[x][y].rows() != 0); };
  int interfaceSize(int x, int y)    const { return _interfaceU[x][y].rows(); };
  MATRIX interfaceMassMatrix(int x, int y);

  // return a matrix of all the cross product matrices produced by the 
  // rest vertices along interface (x,y)
  MATRIX interfaceRestTildeU(int x, int y);

  // return a matrix of all the cross product matrices produced by the 
  // rest vertices in partition x
  MATRIX restTildeU(int x);

  // return a tensor of all the cross product tensors produced by the 
  // displacements along interface (x,y)
  TENSOR3 interfaceTildeU(int x, int y);
  
  // return a tensor of all the cross product tensors produced by the 
  // displacements in partition x
  TENSOR3 tildeU(int x);
  
  // return a tensor of all the cross product tensors produced by the 
  // vertices along interface (x,y)
  MATRIX interfaceTildeVertices(int x, int y);

  // properly ordered rest pose vertices along interface x,y
  VECTOR interfaceRests(int x, int y);
 
  // properly ordered rest pose vertices for an entire partition
  VECTOR restVector(int x) { return _meshes[x]->restVector(); };

  // update rotations of all unconstrained meshes
  void updateRotations();  
 
  // compute interface stiffness with rotations
  void updateBlockInterfaceStiffness();
 
  // return the spring constant being used across the x,y interface
  //Real interfaceStiffness(int x, int y) { return _springConst * _interfaceArea[x][y]; };
  
  // return the spring constant being used across the x,y interface
  Real interfaceSpringConst(int x, int y) { return _interfaceSpringConst(x,y); };
  const MATRIX& interfaceSpringConsts() const { return _interfaceSpringConst; };

  // get the rigid components of a partition
  MATRIX3 rigidRotation(int partition);
  virtual QUATERNION quaternionRotation(int partition);
  virtual QUATERNION quaternionRotationOld(int partition);
  virtual Real rotationLambda(int partition);
  MATRIX3 rigidRotationOld(int partition);
  virtual VEC3F rigidTranslation(int partition);
  VEC3F rigidTranslationOld(int partition);

  SANDWICH_TRANSFORM& translationSandwich(int x, int y) { return _translationSandwich[x][y]; };
  SANDWICH_TRANSFORM& restSandwich(int x, int y) { return _restSandwich[x][y]; };

  SANDWICH_TRANSFORM& constraintRestSandwich(int x, int y) { return _constraintRestSandwich[x][y]; };
  SANDWICH_TRANSFORM& constraintBasisSandwich(int x, int y) { return _constraintBasisSandwich[x][y]; };

  BLOCK_VECTOR& projectedRestInterfaces() { return _projectedRestInterfaces; };
  VEC3F meanRestInterface(int x, int y) { return _meanRestInterface[x][y]; };
  MATRIX& meanBasisInterface(int x, int y) { return _meanBasisInterface[x][y]; };
  void updateSpringConst(Real newSpringConst);
  void updateSpringConst(int x, int y, Real newSpringConst);
  
  // take into account transforms when picking
  virtual int closestPartition(VEC3F point);
 
  // write out the subbases
  void writeSubbases();

  // write out a single subbasis
  void writeSubbasis(int partition);
 
  // update meshes with the subspace
  void updateSubspaceMeshes();

  // update the mesh embedding
  void updateEmbedding();

  // reset everything to zero
  void reset();

  // read/write out the interface spring const matrix
  void writeInterfaceSpringConsts(string directory);
  void readInterfaceSpringConsts(string directory);

  // compute the center of mass over all partitions
  virtual void updateCentersOfMass();

  // get the reduced mass for a partition
  MATRIX& reducedMass(int partition) { return ((SUBSPACE_TET_MESH*)_meshes[partition])->reducedMass(); };

  // extract the portion of this vector that is in the submesh
  VECTOR getSubvector(const int partition, const VECTOR& fullVector);

  // hack to determine if this is a PARTITIONED_SKINNED_SUBSPACE_TET_MESH or not
  virtual bool hasSkeleton() { return false; };

  // add the rigid transforms to the surface vertices
  void addRigidsToSurface();

  // undo the rigid transforms from the surface vertices
  void subtractRigidsFromSurface();

  // write out state for regression testing
  void writeState(const char* filename);

protected:

  // U matrices for each partition corresponding to the
  // shared nodes between partitions
  //
  // ie _interfaceU[x][y] is the U for the surface nodes between
  // partition x and y
  MATRIX** _interfaceU;
  
  // differences between the rest poses across interfaces
  VECTOR** _restDiffs;
  
  // stiffness of the entire interface -- assumes that _superQ will be
  // the multiplier.
  SPARSE_MATRIX _interfaceStiffness;

  // block matrix version of the stiffness of the interface,
  // split into diagonal and off-diagonal parts
  BLOCK_MATRIX _diagonalInterfaceStiffness;
  BLOCK_MATRIX _offDiagonalInterfaceStiffness;

  // the size of the _superQ
  int _superRank;

  // the size of the _superQ, plus rigid components
  int _superRankPlusRigids;

  // rank of each of the partitions
  int* _ranks;
  
  // all the qs from the partitions stacked into one
  VECTOR _superQ;

  // sandwich matrices
  SANDWICH_TRANSFORM** _restSandwich;
  SANDWICH_TRANSFORM** _basisSandwich;
  SANDWICH_TRANSFORM** _translationSandwich;

  // mean interfae sandwich matrices
  SANDWICH_TRANSFORM** _constraintRestSandwich;
  SANDWICH_TRANSFORM** _constraintBasisSandwich;

  // mean interface position for partitions that are constrained
  VEC3F** _meanRestInterface;
  
  // mean of basis for partitions that are constrained
  MATRIX** _meanBasisInterface;

  // projected rest poses
  BLOCK_VECTOR _projectedRestInterfaces;

  // the total number of interfaces
  int _totalInterfaces;

  // compute the U for the interfaces
  void computeInterfaceUs();

  // compute the differences between the rest poses
  // across interfaces
  void computeRestDiffs();

  // compute the stiffness along the interface
  void computeInterfaceStiffness();

  // compute the stiffness along the interface
  void computeBlockInterfaceStiffness();

  // compute the sandwiches for frames
  void computeSandwiches();

  // compute the projected rest interfaces
  void computeProjectedRestInterfaces();

  // return a block identity matrix the dimensions of the interface at (x,y)
  SPARSE_MATRIX interfaceIdentity(int x, int y);

  // recompute the masses, dividing mass evenly between cloned nodes
  virtual void recomputePartitionedMass();
};

#endif
