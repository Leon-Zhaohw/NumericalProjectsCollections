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
// UNCONSTRAINED_FULLSPACE_INTEGRATOR.h: interface for the UNCONSTRAINED_FULLSPACE_INTEGRATOR class.
//
//////////////////////////////////////////////////////////////////////

#ifndef UNCONSTRAINED_FULLSPACE_INTEGRATOR_H
#define UNCONSTRAINED_FULLSPACE_INTEGRATOR_H

#include "FULLSPACE_INTEGRATOR.h"
#include "UNCONSTRAINED_TET_MESH.h"

//////////////////////////////////////////////////////////////////////
// This child class of FULLSPACE_INTEGRATOR handles unconstrained meshes by
// tracking the rigid and deformable components separately. The
// original FULLSPACE_INTEGRATOR class handles unconstrained meshes just fine,
// but this class is a lead-up to handling reduced unconstrained meshes.
//
// The methodology is from:
//
// "Physically Based Models with Rigid and Deformable Components"
// Terzopoulos and Witkin, IEEE CG&A, 1988
//
// Their notation is followed wherever possible. Unless otherwise
// noted, the equation numbers refer to the numbers from that paper
//
//////////////////////////////////////////////////////////////////////
class UNCONSTRAINED_FULLSPACE_INTEGRATOR : public FULLSPACE_INTEGRATOR {

public:
  UNCONSTRAINED_FULLSPACE_INTEGRATOR(UNCONSTRAINED_TET_MESH* tetMesh, Real dt = 1.0 / 60.0, Real alpha = 0.01, Real beta = 0.01);
  ~UNCONSTRAINED_FULLSPACE_INTEGRATOR() {};

  // various solver options are available
  void stepFullImplicit();
  void stepFullExplicit();
  void stepSparseImplicit();
  void stepSparseExplicit();

  // same as stepSparseImplicitAcceleration, but with acceleration as the primary variable
  virtual void stepImplicitAcceleration();

  /*
  VEC3F& rigidTranslation() { return _translation; };
  VEC3F& rigidTranslationOld() { return _translationOld; };
  MATRIX3& rigidRotation() { return _rotation; };
  MATRIX3& rigidRotationOld() { return _rotationOld; };
  VEC3F&   omega()         { return _omega; };
  */

  VEC3F& translation() { return _translation; };
  VEC3F& translationOld() { return _translationOld; };
  VEC3F& translationVelocity() { return _translationVelocity; };
  VEC3F& translationVelocityOld() { return _translationVelocityOld; };
  VEC3F& translationAcceleration() { return _translationAcceleration; };
  VEC3F& translationAccelerationOld() { return _translationAccelerationOld; };
  VEC3F& translationForce() { return _translationForce; };

  QUATERNION& rotation() { return _rotation; };
  QUATERNION& rotationOld() { return _rotationOld; };
  VEC3F& angularDelta() { return _angularDelta; };
  VEC3F& angularDeltaOld() { return _angularDeltaOld; };
  VEC3F& angularVelocity() { return _angularVelocity; };
  VEC3F& angularVelocityOld() { return _angularVelocityOld; };
  VEC3F& angularAcceleration() { return _angularAcceleration; };
  VEC3F& angularAccelerationOld() { return _angularAccelerationOld; };
  TENSOR3& partialDelta() { return _partialDelta; };

  // the collision hack, taking into account rigid components
  bool collideWithFloor(VEC3F down);

  virtual void printRigidComponents();
  virtual void initializeQuasistaticStep();
  virtual void initializeImplicitStep();
  virtual void finalizeQuasistaticStep();
  virtual void finalizeImplicitAccelerationStep();
  virtual void addGravity();
  virtual void generateQuasistaticMatrices(map<string,double>& timingBreakdown);
  virtual void clearForces() { _translationForce.clear(); FULLSPACE_INTEGRATOR::clearForces(); };
  void generateImplicitAccelerationMatrices(SPARSE_MATRIX& A, map<string, double>& timingBreakdown);
  
  // update the current state using acceleration -- does not update olds however
  virtual void updateStateUsingAcceleration();

  // update rigid body components
  void updateRigids();
  
protected:
  VEC3F _translationForce;
  VEC3F _angularForce;

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
  VEC3F _angularDeltaOld;

  VEC3F _angularAcceleration;
  VEC3F _angularAccelerationOld;

  Real _meanMass;   // average mass of a node
  
  VECTOR _xDot;      // aggregate velocity - v + omega ^ q + edot
  
  VEC3F _gV;     // rigid translation force
  VEC3F _gOmega; // rigid rotation force

  VEC3F _meanF;
  VEC3F _meanTorque;
  
  MATRIX3 _inertiaTensor;     // (t)
  MATRIX3 _inertiaTensorOld;  // (t - \delta t)

  // compute external forces 
  void computeXDot();
  void computeGV();
  void computeGOmega();
  void computeInertiaTensor();

  void computeMultibodyExternalForces();

  TENSOR3 _partialDelta;
  TENSOR3 _partialDeltaT;

  // update the exponential derivative
  void updatePartialDelta();
};

#endif
