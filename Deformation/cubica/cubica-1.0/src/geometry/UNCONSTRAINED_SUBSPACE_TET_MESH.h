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
// UNCONSTRAINED_SUBSPACE_TET_MESH.h: interface for the UNCONSTRAINED_SUBSPACE_TET_MESH class.
//
//////////////////////////////////////////////////////////////////////

#ifndef UNCONSTRAINED_SUBSPACE_TET_MESH_H
#define UNCONSTRAINED_SUBSPACE_TET_MESH_H

#include <SUBSPACE_TET_MESH.h>
#include <SANDWICH_TRANSFORM.h>
#include <QUATERNION.h>

//////////////////////////////////////////////////////////////////////
// Unconstrained subspace tet mesh (tracks rigid modes)
//////////////////////////////////////////////////////////////////////
class UNCONSTRAINED_SUBSPACE_TET_MESH : public SUBSPACE_TET_MESH {

public:
  UNCONSTRAINED_SUBSPACE_TET_MESH(const char* filename, MATERIAL** materials,
                    int totalMaterials, bool simulateFullspace = false,
                    const char* eigenvectors = NULL, const char* eigenvalues = NULL,
                    bool projectTranslation = true);
  ~UNCONSTRAINED_SUBSPACE_TET_MESH();

  // accessors
  VEC3F& rigidTranslation() { return _translation; };
  VEC3F& rigidTranslationOld() { return _translationOld; };
  VEC3F& rigidTranslationNewtonOld() { return _translationNewtonOld; };
  VEC3F& originalCenterOfMass() { return _originalCenterOfMass; };
  //MATRIX3& rigidRotation() { return _rotation; };
  //MATRIX3& rigidRotationOld() { return _rotationOld; };
  QUATERNION& rotationQuaternion() { return _rotationQuaternion; };
  QUATERNION& rotationQuaternionOld() { return _rotationQuaternionOld; };
  Real& rotationLambda() { return _rotationLambda; };
  //MATRIX3& rigidRotationNewtonOld() { return _rotationNewtonOld; };
  
  void stepRigid(VEC3F& velocity, VEC3F& omega, Real dt);

  // drawing routines, taking into account the rigid components
  void drawRigidOnly();
  void drawRigidDebug();
  virtual void drawSurfaceFaces();
  virtual void drawAllTets();
  /*
  void drawFacesWithoutRigidTransform();
  */
  void drawRigidFrame();
  void drawCenterOfMass();
  virtual void drawConstrainedNodes();

  // find closest node to an arbitrary point, taking into account
  // the rigis components
  virtual VEC3F* closestSurfaceNode(VEC3F point);

  void resetRigidTranslation();

  // compute the rigid rotation
  //void updateRotation(int debug);
  bool updateRotation(int debug, VECTOR& velocity, VECTOR& acceleration);

  // precache subspace center of mass matrix
  void cacheSubspaceCenterOfMass();

  // compute the center of mass in the subspace
  VEC3F computeSubspaceCenterOfMass();

  SANDWICH_TRANSFORM& UTRU()  { return _UTRU; };
  SANDWICH_TRANSFORM& UTRp()  { return _UTRp; };
  VECTOR& projectedRestPose() { return _projectedRestPose; };
  SANDWICH_TRANSFORM& translationSandwich() { return _translationSandwich; };
  int rotationsKept() { return _totalRotationsKept; };
  int totalGeorgiiIterations() { return _totalGeorgiiIterations; };

  // compute the full rank shape matching rotation
  MATRIX3 computeFullRankShapeMatchingRotation(bool centerCheck = true);
  
  // compute the full rank Georgii and Westermann rotation
  QUATERNION computeFullRankGeorgiiRotation(unsigned int digits = 6, QUATERNION guess = QUATERNION(0,0,0,1), bool useQ = true);

  // compute the reduced Georgii and Westermann rotation
  QUATERNION computeGeorgiiRotation(int debug = 0);

  // compute and update the translation component
  void updateRigidTranslation();

  // compute rotation using shape matching
  MATRIX3 computeShapeMatchingRotation();

  // check that the center of mass is in fact centered at zero
  bool centerOfMassIsZero();

  // reset state to zero
  virtual void reset();

  // recenter the mesh after modifying the mass matrix
  void recenterMesh();

protected:
  ///////////////////////////////////////////////////////////////////
  // Hybrid Rigid/Deformable variables and functions
  ///////////////////////////////////////////////////////////////////
  VEC3F _translationNewtonOld;    // previous rigid translation
                                  // from previous Newton iteration,
                                  // not full time step
  VEC3F _translationOld;    // rigid translation
  VEC3F _translation;       // rigid translation
  //MATRIX3 _rotationOld;        // rigid rotation;
  MATRIX3 _rotationNewtonOld;  // rigid rotation
                               // from previous Newton iteration,
                               // not full time step
  //MATRIX3 _rotation;        // rigid rotation;
  QUATERNION _rotationQuaternion; // unnormalized rigid rotation quaternion
  QUATERNION _rotationQuaternionOld; // unnormalized rigid rotation quaternion
  Real _rotationLambda;
  
  MATRIX3 _rotationDebug;        // rigid rotation;

  VECTOR _qNewtonOld;

  VEC3F _originalCenterOfMass;  // keep the original center of mass around,
                                // before it's subtracted out in the constructor

  ///////////////////////////////////////////////////////////////////
  // Support vars for rigid rotation computation
  ///////////////////////////////////////////////////////////////////
  MATRIX _crossMatrix;   // \II_{sum} \CC_{all} \UU

  // shape matching needs to compute a bulk outer product:
  // A (p + Uq), where A is the outer product operator.
  // Ap and AU can be precomputed, as they are here.
  MATRIX3 _Ap;
  MATRIX _AU;
 
  MATRIX _subspaceCenterOfMass;
  VECTOR _centerOfMassRest;

  // rotation update support variables
  SANDWICH_TRANSFORM _restSandwich;
  SANDWICH_TRANSFORM _basisSandwich;
  SANDWICH_TRANSFORM _translationSandwich;
  VECTOR _projectedRestPose;

  // error estimation support variables
  SANDWICH_TRANSFORM _UTRU;
  SANDWICH_TRANSFORM _UTRp;
  SANDWICH_TRANSFORM _ITRp;
  SANDWICH_TRANSFORM _ITRU;
  SANDWICH_TRANSFORM _pTRp;
  Real _pDot;

  SANDWICH_TRANSFORM _pTripleTRU;
  SANDWICH_TRANSFORM _pTripleTRp;

  VECTOR _UTp;
  VEC3F _ITp;

  // DEBUG
  int _totalRotationsSeen;
  int _totalRotationsKept;
  Real _meanRotationError;

  // initialize different rotation calculation types
  void initializeCrossProductRotation();
  void initializeShapeMatchingRotation();

  // compute rotation using a monolithic cross product
  void computeCrossProductRotation();

  // compute the error norm for the new rotation
  Real computeRotationErrorNorm(VECTOR& q,    MATRIX3& rotation,    VEC3& translation,
                                VECTOR& qOld, MATRIX3& rotationOld, VEC3& translationOld, int debug);

  ///////////////////////////////////////////////////////////////////
  // Support functions for computeGeorgiiRotation
  ///////////////////////////////////////////////////////////////////
  QUATERNION fullRotationDistance(QUATERNION& q, VEC3F& rest, VEC3F& displacement);
  QUATERNION fullDistanceHessian(QUATERNION& q, QUATERNION& p, int i, int j);
  QUATERNION fullDistanceGradient(QUATERNION& q, QUATERNION& p, int which);
  MATRIX fullEnergyHessian(QUATERNION& q, Real lambda);
  VECTOR fullEnergyGradient(QUATERNION& q, Real lambda);
  VECTOR projectedEnergyGradient(QUATERNION& q, Real lambda);
  Real fullRotationEnergy(QUATERNION& q, Real lambda = 0);
  int _georgiiIterations; // track how many Newton iterations the rotation took
  int _totalGeorgiiIterations; // track how many Newton iterations the rotation took
  int _maxGeorgiiIterations;
  int _totalSteps;

  MATRIX hessianMatrix(QUATERNION& q, int i, int j);
  MATRIX gradientMatrix(QUATERNION& q, int which);
  MATRIX energyHessian(QUATERNION& q, Real lambda);
  VECTOR energyGradient(QUATERNION& q, Real lambda);
  Real rotationEnergy(QUATERNION& q, Real lambda = 0);

  // sandwiches for Georgii rotations
  SANDWICH_TRANSFORM _IBU_U;
  SANDWICH_TRANSFORM _IBU_IV;
  SANDWICH_TRANSFORM _IBU_IBU;
  SANDWICH_TRANSFORM _IV_U;
  SANDWICH_TRANSFORM _IV_IV;
  SANDWICH_TRANSFORM _IV_IBU;

  // for debugging only -- unpack, project, expand, and repack a vector
  VECTOR projectQuaternionVector(VECTOR& d);
  
  // for debugging only -- make a padded U for quaternions
  void buildQU(BLOCK_MATRIX& qU);
  void buildBlockU(BLOCK_MATRIX& blockU);

  // read/write basis with translation projected out
  bool readTranslationlessBasis();
  void writeTranslationlessBasis();

  // read/write variable mass matrix vars
  bool readVariableMassMatrix(string filename = string(""));
  void writeVariableMassMatrix(string filename = string(""));

};

#endif
