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
// UNCONSTRAINED_SUBSPACE_INTEGRATOR.h: interface for the UNCONSTRAINED_SUBSPACE_INTEGRATOR_H class.
//
//////////////////////////////////////////////////////////////////////

#ifndef UNCONSTRAINED_SUBSPACE_INTEGRATOR_H
#define UNCONSTRAINED_SUBSPACE_INTEGRATOR_H

#include <SETTINGS.h>
#include <SUBSPACE_INTEGRATOR.h>
#include <SANDWICH_TRANSFORM.h>
#include <QUATERNION.h>

//////////////////////////////////////////////////////////////////////
// Timestepping class for the tet mesh, follows the notation
// of [Barbic and James], Appendix C. The variables names
// adhere as closely as possible to those therein.
//
// This does not explicitly subclass FULLSPACE_INTEGRATOR, it just
// mirrors the functionality. The two classes are sufficiently
// different that subclassing would just introduce compatibility 
// headaches. You wouldn't swap one with the other at runtime anyway,
// so there's no reason to keep their pointers interchangeable.
//////////////////////////////////////////////////////////////////////
class UNCONSTRAINED_SUBSPACE_INTEGRATOR : public SUBSPACE_INTEGRATOR {

public:
  UNCONSTRAINED_SUBSPACE_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, Real dt = 1.0 / 60.0, Real alpha = 0.01, Real beta = 0.01, Real gravity = 9.8);
  virtual ~UNCONSTRAINED_SUBSPACE_INTEGRATOR() {};

  void initializeImplicitStep();
  void generateImplicitMatrices(map<string, double>& timingBreakdown);
  void generateImplicitMatricesInvertible(map<string, double>& timingBreakdown);
  void finalizeImplicitStep();
  void finalizeImplicitStepPlusTranslations();
  void generateQuasistaticMatrices(map<string,double>& timingBreakdown);
  // UNIMPLEMENTED. Only included here to throw a flag.
  virtual void initializeQuasistaticStep();
  virtual void finalizeQuasistaticStep();
  virtual void finalizeImplicitAccelerationStep();
  virtual void finalizeImplicitVelocityStep();

  VEC3F& translation() { return _translation; };
  VEC3F& translationOld() { return _translationOld; };
  VEC3F& translationVelocity() { return _translationVelocity; };
  VEC3F& translationVelocityOld() { return _translationVelocityOld; };
  VEC3F& translationAcceleration() { return _translationAcceleration; };
  VEC3F& translationAccelerationOld() { return _translationAccelerationOld; };
  VEC3F& translationForce() { return _translationForce; };

  QUATERNION& rotation() { return _rotation; };
  QUATERNION& rotationOld() { return _rotationOld; };
  VEC3F& angularVelocity() { return _angularVelocity; };
  VEC3F& angularVelocityOld() { return _angularVelocityOld; };
  VEC3F& angularAcceleration() { return _angularAcceleration; };
  VEC3F& angularAccelerationOld() { return _angularAccelerationOld; };
  VEC3F& angularDelta() { return _angularDelta; };
  BLOCK_VECTOR& quadraticForces() { return _quadraticForces; };
  TENSOR3& partialDelta() { return _partialDelta; };
  TENSOR3& partialDeltaT() { return _partialDeltaT; };
  VEC3F& torque() { return _torque; };

  //MATRIX3& rigidRotation() { return _rotation; };
  //MATRIX3& rigidRotationOld() { return _rotationOld; };
  VECTOR& rotatedExternalForces() { return _reducedRotatedExternals; };

  virtual void reset();

  // compute and return full coordinate velocity
  virtual VECTOR fullCoordinateVelocity();
  
  // compute and return full coordinate velocity
  virtual VECTOR computeNewestReducedVelocity();

  // Add gravity only to the translation component
  //virtual void addGravity() { _translationForce += _gravityMagnitude * _gravityDown; };
  virtual void addGravity();
  virtual void addBodyForce(VEC3F bodyForce);
  virtual void clearForces() { 
    _forceVectors.clear(); 
    _forceNodes.clear(); 
    _forceNodeIDs.clear(); 
    _externalForces.clear(); 
    _translationForce *= 0.0; 
    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //cout << " EXTERNAL CLEAR DISABLED " << endl;
    _quadraticForces *= 0.0; 
  };

  // update rigid body components
  void updateRigids();

  // update the current state using acceleration -- does not update olds however
  virtual void updateStateUsingAcceleration();
  virtual void updateStateUsingVelocity();

  // update the quadratic forces
  void computeQuadraticForces();
  void computeQuadraticForcesReduced();

  // computes only the defo quadractic velocity -- assumes the rigid components
  // are specified kinematically
  void computeSkinnedQuadraticForcesReduced();

  virtual void drawClickedNode();
  virtual void drawForceVector();
  virtual bool click(VEC3F& point, Real maxRadius = -1.0);

  // unconstrained timestepping, acceleration-level only support for now
  virtual int stepImplicitAcceleration();
  virtual int stepImplicitWithExplicitCollisions(vector<pair<SURFACE*, int> >& collisions);
  virtual int stepImplicitWithImplicitCollisions(vector<pair<SURFACE*, int> >& collisions);

  // add implicit collision forces directly to the residual
  //BLOCK_VECTOR blockImplicitCollisionResiduals(vector<pair<SURFACE*, int> >& collisions);
  BLOCK_VECTOR blockImplicitCollisionResiduals();
  BLOCK_VECTOR blockImplicitCollisionResidualsDebug();

  // add implicit collision forces directly to the Jacobian
  //BLOCK_MATRIX blockImplicitCollisionJacobians(vector<pair<SURFACE*, int> >& collisions);
  //BLOCK_MATRIX blockImplicitCollisionJacobians(map<string, double>& timingBreakdown);
  //BLOCK_MATRIX blockImplicitCollisionJacobiansDebug();
  void blockImplicitCollisionJacobians(BLOCK_MATRIX& systemMatrix, map<string, double>& timingBreakdown);

  // build the list of vertices in collision -- cull away those that are a member
  // of a colliding triangle but actually outside the surface
  void buildCollisionList(vector<pair<SURFACE*, int> >& collisions);

  vector<VEC3F>& collisionForces() { return _collisionForces; };

  // compute the kinetic energy of the entire mesh
  virtual Real kineticEnergy();

  // compute the rigid derivatives based on the current rigid positions
  void computeRigidDerivatives();

protected:
  VEC3F _translationForce;

  VEC3F _translation;
  VEC3F _translationOld;

  VEC3F _translationVelocity;
  VEC3F _translationVelocityOld;

  VEC3F _translationAcceleration;
  VEC3F _translationAccelerationOld;

  QUATERNION _rotation;
  QUATERNION _rotationOld;

  VEC3F _angularVelocity;
  VEC3F _angularVelocityOld;
  VEC3F _angularDelta;

  VEC3F _angularAcceleration;
  VEC3F _angularAccelerationOld;

  VEC3F _torque;

  VECTOR _reducedRotatedExternals;

  BLOCK_VECTOR _quadraticForces;

  VECTOR _externalRotationForces;

  // sandwiches needed for quadratic force computation
  SANDWICH_TRANSFORM _MUi_R_Ui;
  SANDWICH_TRANSFORM _MUi_R_uiBarO;

  TENSOR3 _partialDelta;
  TENSOR3 _partialDeltaT;
  vector<pair<VEC3F*, SURFACE*> > _collisionVertices;
  vector<SURFACE*> _collisionSurfaces;

  // cache the collision forces for Jacobian computation
  vector<VEC3F> _collisionForces;

  void updatePartialDelta();

  void computeMassMatrix(BLOCK_MATRIX& blockMass);
  void computeResidual(BLOCK_VECTOR& blockResidual);

  // add explicit collision forces
  //virtual void addExplicitCollisions(vector<pair<SURFACE*, int> >& collisions);
  virtual void addExplicitCollisions();

public:  
  // verify the correctness of the collision jacobian
  //void verifyCollisionJacobian(vector<pair<SURFACE*, int> >& collisions);
  void verifyCollisionJacobian(vector<pair<SURFACE*, int> >& collisions);

public:
  // compute versions of all the old quantities transformed into the new frame
  void oldsTransformed(VECTOR& positionOldRotated, 
                       VECTOR& velocityOldRotated, 
                       VECTOR& accelerationOldRotated);
  void newtonOldsTransformed(VECTOR& positionOldRotated, 
                       VECTOR& velocityOldRotated, 
                       VECTOR& accelerationOldRotated);

  void computeExternalForces();

  // compute in full coordinates the position in the local frame
  // taking into account the translation
  VECTOR projectedPositionPlusTranslation(VECTOR& position, VEC3F& translation);
};

#endif
