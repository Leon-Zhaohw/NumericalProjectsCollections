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
// UNCONSTRAINED_FULLSPACE_INTEGRATOR.h: interface for the 
// UNCONSTRAINED_FULLSPACE_INTEGRATOR class.
//
//////////////////////////////////////////////////////////////////////

#include "UNCONSTRAINED_FULLSPACE_INTEGRATOR.h"
#include "TIMER.h"

//////////////////////////////////////////////////////////////////////
// Constructor for newmark integrator
//////////////////////////////////////////////////////////////////////
UNCONSTRAINED_FULLSPACE_INTEGRATOR::UNCONSTRAINED_FULLSPACE_INTEGRATOR(UNCONSTRAINED_TET_MESH* tetMesh, Real dt, Real alpha, Real beta) :
  FULLSPACE_INTEGRATOR(tetMesh, dt, alpha, beta)
{
  _rotation = MATRIX3::I();
  _rotationOld = MATRIX3::I();

  _translationOld = tetMesh->rigidTranslationOld();
  _translation = tetMesh->rigidTranslation();

  /*
  // cache the rest state inertial tensor
  computeInertiaTensor();
  _inertiaTensorOld = _inertiaTensor;

  // compute the mean mass once and for all
  _meanMass = 0;
  for (int x = 0; x < _tetMesh->unconstrainedNodes(); x++)
    _meanMass += _tetMesh->mass(x);
  _meanMass /= _tetMesh->unconstrainedNodes();
  */
}

//////////////////////////////////////////////////////////////////////
// Explicit Integration solved with a full LU factorization
//
// This follows pages 324-325 of: 
// Computational Contact Mechanics by Peter Wriggers, 2006
//
// The equation to solve is:
// (M + (dt/2) * C) u_{n+1} = 
// (dt)^2 * [P_n - R(u_n)] + (dt/2) * C * u_{n-1} + M* (2 * u_n - u_{n-1})
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::stepFullExplicit()
{
  cout << " UNCONSTRAINED_FULLSPACE_INTEGRATOR::stepFullExplicit() integrator is untested! " << endl;

  // this will probably work if you just cut and paste the code
  // from stepSparseExplicit() here.
  
  /*
  // call FULLSPACE_INTEGRATOR, but computeExternalForces will
  // be called from this class
  //
  // rigid translation and rotation are computed inside 
  // computeExternalForces as well
  FULLSPACE_INTEGRATOR::stepFullExplicit();
 
  // update rigid components
  ((UNCONSTRAINED_TET_MESH*)_tetMesh)->stepRigid(_v, _omega, _dt);
  */
}

//////////////////////////////////////////////////////////////////////
// Explicit Integration
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::stepSparseExplicit()
{
  /*
  // call FULLSPACE_INTEGRATOR, but computeExternalForces will
  // be called from this class
  //
  // rigid translation and rotation are computed inside 
  // computeExternalForces as well
  FULLSPACE_INTEGRATOR::stepSparseExplicit();

  // compute rigid translation
  _vOld = _v;
  _v = _vOld + _dt * _gV / _meanMass;

  // compute rigid rotation
  // this appears to use the inertiaTensor from one step before,
  // but it seems to give the correct answer wrt the full integrator?
  _omega = _inertiaTensor.inverse() * (_inertiaTensorOld * _omegaOld + 
                                       _dt * _gOmega);
  _inertiaTensorOld = _inertiaTensor;
  computeInertiaTensor();
  _omegaOld = _omega;

  // update rigid components
  ((UNCONSTRAINED_TET_MESH*)_tetMesh)->stepRigid(_v, _omega, _dt);
  */
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::stepFullImplicit()
{
  cout << " UNCONSTRAINED_FULLSPACE_INTEGRATOR::stepFullImplicit() integrator is untested! " << endl;

  // Not as straightforward as explicit --
  // the fact the Terzopoulos does only explicit integration
  // will probably wreak havoc here.
  
  /*
  // call FULLSPACE_INTEGRATOR, but computeExternalForces will
  // be called from this class
  //
  // rigid translation and rotation are computed inside 
  // computeExternalForces as well
  FULLSPACE_INTEGRATOR::stepFullImplicit();

  // update rigid components
  ((UNCONSTRAINED_TET_MESH*)_tetMesh)->stepRigid(_v, _omega, _dt);
  */
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::stepSparseImplicit()
{
  cout << " UNCONSTRAINED_FULLSPACE_INTEGRATOR::stepSparseImplicit() integrator is untested! " << endl;

  // Not as straightforward as explicit --
  // the fact the Terzopoulos does only explicit integration
  // will probably wreak havoc here.
  //
  /*
  // call FULLSPACE_INTEGRATOR, but computeExternalForces will
  // be called from this class
  //
  // rigid translation and rotation are computed inside 
  // computeExternalForces as well
  FULLSPACE_INTEGRATOR::stepSparseImplicit();

  // get the new inertia tensor -- we need to do this after the deformable
  // component has been updated so we can compute I based on q_{t+\delta t}
  computeInertiaTensor();

  // update rigid components
  ((UNCONSTRAINED_TET_MESH*)_tetMesh)->stepRigid(_v, _omega, _dt);
  */
}

//////////////////////////////////////////////////////////////////////
// Compute rigid translation components
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::computeGV()
{
  int size = _tetMesh->unconstrainedNodes();

  // get the rigid external force - Eqn. 33a
  _gV.clear();
  for (int x = 0; x < size; x++)
  {
    int index = 3 * x;
    VEC3F force;
    force[0] = _externalForces(index); 
    force[1] = _externalForces(index + 1); 
    force[2] = _externalForces(index + 2); 
    _gV += force;
  }
  _gV *= 1.0 / size;
  _meanF = _gV;

  // compute the \frac{d}{dt} \mu \dot{e} term (same as M\ddot{x})
  SPARSE_MATRIX& M = _tetMesh->massMatrix();
  VECTOR product = M * _acceleration;
  VEC3F sum;
  for (int x = 0; x < size; x++)
  {
    int index = 3 * x;
    sum[0] += product(index);
    sum[1] += product(index + 1);
    sum[2] += product(index + 2);
  }

  if (_CDampingSparse.rows() != _rank)
  {
    _CDampingSparse.resize(_rank, _rank);
    cacheStaticDampingSparse();
  }
  product = _CDampingSparse * _xDot;
  VEC3F dampingSum;
  for (int x = 0; x < size; x++)
  {
    int index = 3 * x;
    VEC3F damping(product(index), product(index + 1), product(index + 2));
    sum += damping;
    dampingSum += damping;
  }
  sum *= 1.0 / size;
  
  // subtract M\ddot{x} + C\dot{x} from f
  _gV -= sum;
}

//////////////////////////////////////////////////////////////////////
// compute the inertia tensor
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::computeInertiaTensor()
{
  _inertiaTensor.clear();
  for (int x = 0; x < _tetMesh->unconstrainedNodes(); x++)
  {
    VEC3F q = *(_tetMesh->vertices(x));
    Real squared = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
    MATRIX3 temp;
    temp(0,0) = squared - q[0] * q[0];
    temp(1,1) = squared - q[1] * q[1];
    temp(2,2) = squared - q[2] * q[2];
    temp(0,1) = temp(1,0) = -q[0] * q[1];
    temp(0,2) = temp(2,0) = -q[0] * q[2];
    temp(1,2) = temp(2,1) = -q[1] * q[2];

    _inertiaTensor += _tetMesh->mass(x) * temp;
  }
}

//////////////////////////////////////////////////////////////////////
// compute rigid rotation components
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::computeGOmega()
{
  int size = _tetMesh->unconstrainedNodes();
 
  // compute q cross f
  VEC3F qCrossF;
  for (int x = 0; x < size; x++)
  {
    int index = 3 * x;
    VEC3F q = *(_tetMesh->vertices(x));
    VEC3F f(_externalForces(index), 
            _externalForces(index + 1), 
            _externalForces(index + 2));
    VEC3F cross = q ^ f;
    qCrossF += cross;
  }
  _meanTorque = qCrossF;

  // compute q cross \mu \ddot{e}
  SPARSE_MATRIX& M = _tetMesh->massMatrix();
  VECTOR product = M * _acceleration;
  VEC3F sum;
  for (int x = 0; x < size; x++)
  {
    int index = 3 * x;
    VEC3F q = *(_tetMesh->vertices(x));
    VEC3F e(product(index),
            product(index + 1),
            product(index + 2));
    VEC3F cross = q ^ e;
    sum += q ^ e;
  }

  // compute q cross velocity first
  if (_CDampingSparse.rows() != _rank)
  {
    _CDampingSparse.resize(_rank, _rank);
    cacheStaticDampingSparse();
  }
  product.clear();
  for (int x = 0; x < size; x++)
  {
    VEC3F q = *(_tetMesh->vertices(x));
    int index = 3 * x;
    VEC3F xDot(_xDot(index), _xDot(index + 1), _xDot(index + 2));
    VEC3F cross = q ^ xDot;
    product(index)     = cross[0];
    product(index + 1) = cross[1];
    product(index + 2) = cross[2];
  }
  VECTOR final = _CDampingSparse * _xDot;
  for (int x = 0; x < size; x++)
  {
    int index = 3 * x;
    VEC3F torque(final(index), final(index+1), final(index+2));
    VEC3F q = *(_tetMesh->vertices(x));
    VEC3F damping = q ^ torque;
    
    sum += damping;
  }

  _gOmega = qCrossF - sum;
}

//////////////////////////////////////////////////////////////////////
// Compute the sum of all the velocity components
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::computeXDot()
{
  /*
  // compute xdot
  int size = _tetMesh->unconstrainedNodes();
  if (_xDot.size() != 3 * size)
    _xDot.resizeAndWipe(3 * size);

  for (int x = 0; x < size; x++)
  {
    VEC3F q = *(_tetMesh->vertices(x));
    VEC3F cross = _omega ^ q;
   
    int index = 3 * x;
    _xDot(index)     = _v[0] + cross[0] + _velocity(index);
    _xDot(index + 1) = _v[1] + cross[1] + _velocity(index + 1);
    _xDot(index + 2) = _v[2] + cross[2] + _velocity(index + 2);
  }
  */
}

// no longer valid, and the overloading causes all sorts of problems
/*
//////////////////////////////////////////////////////////////////////
// Compute all external forces
//////////////////////////////////////////////////////////////////////  
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::computeExternalForces()
{
  // call superclass
  FULLSPACE_INTEGRATOR::computeExternalForces();

  // compute the sum of all the velocity components
  computeXDot();
  
  // compute _gV - rigid translation
  computeGV();

  // compute _gOmega - rigid rotation
  computeGOmega();

  // subtract the rigid translation from the overall force
  VEC3F vDot = _meanF / _meanMass;
  for (int x = 0; x < _externalForces.size() / 3; x++)
  {
    int index = 3 * x;
    Real mass = _tetMesh->mass(x);

    _externalForces(index) -= mass * vDot[0];
    _externalForces(index + 1) -= mass * vDot[1];
    _externalForces(index + 2) -= mass * vDot[2];
  }

  // subtract the rigid rotation from the overall force
  VEC3F omegaDot = _inertiaTensor.inverse() * _meanTorque;
  for (int x = 0; x < _tetMesh->unconstrainedNodes(); x++)
  {
    VEC3F subtract;
    Real mass = _tetMesh->mass(x);

    // centripetal force
    // \mu \omega \cross \omega \q
    VEC3F q = (*_tetMesh->vertices(x));
    subtract += mass * _omega ^ _omega ^ q;
      
    // 2 \mu \omega \cross \dot{e}
    int index = 3 * x;
    VEC3F velocity;
    velocity[0] = _velocity(index);
    velocity[1] = _velocity(index + 1);
    velocity[2] = _velocity(index + 2);
    subtract += 2.0 * mass * _omega ^ velocity;
 
    // angular acceleration
    // \mu \dot{omega} \cross q
    subtract += mass * (omegaDot ^ q);

    // subtract from the overall force
    _externalForces(index)     -= subtract[0];
    _externalForces(index + 1) -= subtract[1];
    _externalForces(index + 2) -= subtract[2];
  }
}
*/

//////////////////////////////////////////////////////////////////////
// Collision detection and response for a floor
// (Hacky -- for debugging purposes only)
//////////////////////////////////////////////////////////////////////
bool UNCONSTRAINED_FULLSPACE_INTEGRATOR::collideWithFloor(VEC3F down)
{
  float floorPosition = -0.25;
  bool collided = false;
  
  VEC3F up = down;
  VEC3F translate = ((UNCONSTRAINED_TET_MESH*)_tetMesh)->rigidTranslation();
  MATRIX3 rotation = ((UNCONSTRAINED_TET_MESH*)_tetMesh)->rotationQuaternion().toExplicitMatrix3x3();
  MATRIX3 rotationInverse = rotation.inverse();
  up *= -1.0;
  up = rotation * up;
  for (int x = 0; x < _tetMesh->unconstrainedNodes(); x++)
  {
    VEC3F vertex = (*_tetMesh->vertices(x));

    vertex = rotation * vertex;
    vertex += translate;
    if (vertex[1] < floorPosition)
    {
      float diff = floorPosition - vertex[1];
      VEC3F response = up;
      response *= diff * 100.0;
      _forceVectors.push_back(response);
      _forceNodes.push_back(_tetMesh->vertices(x));
      collided = true;
    }
  }
  return collided;
}

//////////////////////////////////////////////////////////////////////
// Compute rigid translation, rotation and velocities
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::printRigidComponents()
{
  //cout << " rigid velocity: " << _v << endl;
  //cout << " rigid rotation: " << _omega << endl;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::initializeQuasistaticStep()
{
  TIMER preamble;

  // accumulate external forces
  computeExternalForces();
  
  // copy the current values into the old values
  _positionOld = _position;
  _translationOld = _translation;
  _rotationOld = _rotation;

  _timingBreakdown["Preamble"] += preamble.timing();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::finalizeQuasistaticStep()
{
  _tetMesh->updateFullMesh();
  _externalForces.clear();
  _totalSteps++;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::addGravity()
{
  // need to retain the option of multiplying by the rotation derivative instead
  // of the pure rotation, so the unrotated version should be stored here.
  /*
  VECTOR rotatedGravity(_gravity.size());
 
  MATRIX3 Ri = _rotation;
  MATRIX3 RiT = _rotation.transpose();
  for (int x = 0; x < _gravity.size() / 3; x++)
  {
    VEC3F nodeGravity;
    nodeGravity[0] = _gravity[3 * x];
    nodeGravity[1] = _gravity[3 * x + 1];
    nodeGravity[2] = _gravity[3 * x + 2];

    nodeGravity = RiT * nodeGravity;
    //nodeGravity = Ri * nodeGravity;

    rotatedGravity[3 * x] = nodeGravity[0];
    rotatedGravity[3 * x + 1] = nodeGravity[1];
    rotatedGravity[3 * x + 2] = nodeGravity[2];
  }

  _externalForces += rotatedGravity;
  */

  _translationForce += _gravityMagnitude * _gravityDown;
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an quasistatic step -- should be called 
// right after initializeQuasistaticStep
//
// this has been forked out into its own function in case we want
// to do more than one Newton step
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::generateQuasistaticMatrices(map<string, double>& timingBreakdown)
{
  /*
  VECTOR& x = _tetMesh->x();
  //TIMER updateMesh;

  // update q vector with current tet mesh state
  // an assignment is overly expensive -- the "=" forces an allocation
  x.equals(_position);
  //timingBreakdown["Update Mesh"] += updateMesh.timing();
  
  // compute the new deformation gradient F (needed by both R and K)
  _tetMesh->generateF();
  
  // precompute diagonalizations --
  // this call is not timed out because timers are inside the function call
  cacheDiagonalizations();
  
  //TIMER internalTimer;
  VECTOR& R = _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);
  //timingBreakdown["Internal Forces"] += internalTimer.timing();

  _tetMesh->generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);

  VECTOR rotatedExternals(_externalForces);
  for (int x = 0; x < _externalForces.size() / 3; x++)
  {
    VEC3F force;
    force[0] = _externalForces[3 * x];
    force[1] = _externalForces[3 * x + 1];
    force[2] = _externalForces[3 * x + 2];

    VEC3F rotatedForce = _rotation.transpose() * force;
    rotatedExternals[3 * x] = rotatedForce[0];
    rotatedExternals[3 * x + 1] = rotatedForce[1];
    rotatedExternals[3 * x + 2] = rotatedForce[2];
  }

  //TIMER RHSTimer;
  //_residual.equals(_externalForces);
  _residual.equals(rotatedExternals);
  _residual += R;

  _totalMeshForces = R;
  _totalMeshForces *= -1.0;

  // account for the fact that the implicit solver does not negate
  // the gradient
  //
  // No, this is just to push the R + f_ext term to the LHS
  _residual *= -1.0;
  //timingBreakdown["RHS Assembly"] += RHSTimer.timing();
  */
}

//////////////////////////////////////////////////////////////////////
// Build matrices for an invertible implicit step
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::initializeImplicitStep()
{
  TIMER preamble;

  // accumulate external forces
  computeExternalForces();

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;
  _translationOld = _translation;
  _translationAccelerationOld = _translationAcceleration;
  _translationVelocityOld = _translationVelocity;

  _rotationOld = _rotation;
  _angularVelocityOld = _angularVelocity;
  _angularAccelerationOld = _angularAcceleration;

  _tetMesh->inertiaTensorOld() = _tetMesh->inertiaTensor();
  _tetMesh->inertiaTensorDtOld() = _tetMesh->inertiaTensorDt();

  UNCONSTRAINED_TET_MESH* mesh = (UNCONSTRAINED_TET_MESH*)_tetMesh;
  mesh->rotationQuaternionOld() = _rotation;
  mesh->rigidTranslationOld() = _translation;

  _timingBreakdown["Preamble"] += preamble.timing();
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an implicit step -- should be called 
// right after initializeImplicitStep
//
// this has been forked out into its own function in case we want
// to do more than one Newton step
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::generateImplicitAccelerationMatrices(SPARSE_MATRIX& A, map<string, double>& timingBreakdown)
{
  Real* alpha = _accelerationAlpha;
  VECTOR& position = _tetMesh->x();
  SPARSE_MATRIX& M = _tetMesh->massMatrix();
  _position = _positionOld + alpha[2] * _velocityOld +
                             alpha[3] * _accelerationOld + 
                             alpha[4] * _acceleration;
  _velocity = _velocityOld + alpha[0] * _accelerationOld +
                             alpha[1] * _acceleration;

  // update tet mesh with the new q vector
  // an assignment is overly expensive -- the "=" forces an allocation
  TIMER updateMesh;
  position.equals(_position);
  timingBreakdown["Update Mesh"] += updateMesh.timing();

  // compute the new deformation gradient F (needed by both R and K)
  TIMER timerF;
  _tetMesh->generateF();
  timingBreakdown["Generate F"] += timerF.timing();

  // precompute diagonalizations --
  // this call is not timed out because timers are inside the function call
  TIMER diagonalizeTimer;
  cacheDiagonalizations();
  timingBreakdown["Diagonalize F"] += diagonalizeTimer.timing();

  // get the reduced internal forces
  TIMER internalTimer;
  VECTOR& R = _tetMesh->TET_MESH::generateInternalForces(_Us, _Fhats, _Vs);
  timingBreakdown["Internal Forces"] += internalTimer.timing();

  // get the reduced stiffness matrix -- do this even if we are 
  // doing the conjugate gradient solve, because the damping matrix
  // needs it
  SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);
  
  // compute the damping matrix, but only if a new K is available
  TIMER RHSTimer;
  _CDampingSparse = _rayleighAlpha * M;
  _CDampingSparse += _rayleighBeta * K;

  _temp = M * _acceleration;
  _residual = _CDampingSparse * _velocity;
  _residual += _temp;
  _residual -= R;
  _residual -= _externalForces;
  timingBreakdown["RHS Assembly"] += RHSTimer.timing();

  // use the acceleration level Newmark consts
  TIMER ATimer;
  if (A.rows() == 0)
  {
    A.clear();
    A.resize(M.rows(), M.cols());
  }
  A.equals(M);
  A.axpy(alpha[4], K);
  A.axpy(alpha[1], _CDampingSparse);
  timingBreakdown["A Assembly"] += ATimer.timing();
}

//////////////////////////////////////////////////////////////////////
// Update the mesh after the partitioned implicit step
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::finalizeImplicitAccelerationStep()
{
  FULLSPACE_INTEGRATOR::finalizeImplicitAccelerationStep();
  Real* alpha = _accelerationAlpha;

  // velocity update
  _angularVelocity = _angularVelocityOld + 
                       alpha[0] * _angularAccelerationOld + 
                       alpha[1] * _angularAcceleration;

  _angularDelta = alpha[2] * _angularVelocityOld + 
                  alpha[3] * _angularAccelerationOld + 
                  alpha[4] * _angularAcceleration;

  _translationVelocity = _translationVelocityOld + 
                   alpha[0] * _translationAccelerationOld + 
                   alpha[1] * _translationAcceleration;

  // position update
  QUATERNION update = QUATERNION::fromAxisAngle(_angularDelta);

  _rotation = update * _rotationOld;

  _translation = _translationOld + 
                alpha[2] * _translationVelocityOld + 
                alpha[3] * _translationAccelerationOld + 
                alpha[4] * _translationAcceleration;

  // commit to mesh
  UNCONSTRAINED_TET_MESH* mesh = (UNCONSTRAINED_TET_MESH*)_tetMesh;
  mesh->rotationQuaternion() = _rotation;
  mesh->rigidTranslation() = _translation;
}

//////////////////////////////////////////////////////////////////////
// Update the integrator and mesh state using and acceleration
// level update. Does not update the "olds" however -- this should
// have been done at the beginning of a time step.
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::updateStateUsingAcceleration()
{
  FULLSPACE_INTEGRATOR::updateStateUsingAcceleration();
  updateRigids();

  UNCONSTRAINED_TET_MESH* mesh = (UNCONSTRAINED_TET_MESH*)_tetMesh;
  
  mesh->TET_MESH::refreshInertiaTensor();
  mesh->TET_MESH::refreshInertiaTensorDt(_velocity);

  updatePartialDelta();
}

//////////////////////////////////////////////////////////////////////
// update rigid body components
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::updateRigids()
{
  Real* alpha = _accelerationAlpha;

  // rotation update
  _angularVelocity.equals(_angularVelocityOld);
  _angularVelocity.axpy(alpha[0], _angularAccelerationOld);
  _angularVelocity.axpy(alpha[1], _angularAcceleration);

  _angularDelta.clearingAxpy(alpha[2], _angularVelocityOld);
  _angularDelta.axpy(alpha[3], _angularAccelerationOld);
  _angularDelta.axpy(alpha[4], _angularAcceleration); 

  QUATERNION update = QUATERNION::fromAxisAngle(_angularDelta);
  _rotation = update * _rotationOld;

  UNCONSTRAINED_TET_MESH* mesh = (UNCONSTRAINED_TET_MESH*)_tetMesh;
  mesh->rotationQuaternion().equals(_rotation);

  // translation update
  _translationVelocity.equals(_translationVelocityOld);
  _translationVelocity.axpy(alpha[0], _translationAccelerationOld);
  _translationVelocity.axpy(alpha[1], _translationAcceleration);

  _translation.equals(_translationOld);
  _translation.axpy(alpha[2], _translationVelocityOld);
  _translation.axpy(alpha[3], _translationAccelerationOld);
  _translation.axpy(alpha[4], _translationAcceleration);

  mesh->rigidTranslation().equals(_translation);
}

//////////////////////////////////////////////////////////////////////
// update the exponential derivative
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::updatePartialDelta()
{
  vector<MATRIX> partials;
  Real halfDt2 = _accelerationAlpha[4];
  for (int x = 0; x < 3; x++)
  {
    // depending on the component, the derivative of the
    // cross product (tilde) operator
    MATRIX partial(3,3);
    
    if (x == 0)
    {
      partial(1,2) = -halfDt2;
      partial(2,1) = halfDt2;
    }
    if (x == 1)
    {
      partial(0,2) = halfDt2;
      partial(2,0) = -halfDt2;
    }
    if (x == 2)
    {
      partial(0,1) = -halfDt2;
      partial(1,0) = halfDt2;
    }
    partials.push_back(partial);
  }
  MATRIX deltaATilde = MATRIX::cross(_angularDelta);
  _partialDelta = TENSOR3(deltaATilde, partials);
  _partialDeltaT = _partialDelta.transpose();
}

//////////////////////////////////////////////////////////////////////
// Compute external forces, but the multibody way
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::computeMultibodyExternalForces()
{
  computeExternalForces();

  int unconstrained = _externalForces.size() / 3;

  // get the current rotation
  MATRIX3 R = _rotation.toExplicitMatrix3x3();
  MATRIX3 RT = R.transpose();

  cout << " rotation: " << R << endl;

  // get the current vertices
  vector<VEC3F>& vertices = _tetMesh->vertices();

  // project out the rigid forces
  _translationForce *= 0;
  _angularForce *= 0;
  for (int x = 0; x < unconstrained; x++)
  {
    VEC3F force;
    force[0] = _externalForces[3 * x];
    force[1] = _externalForces[3 * x + 1];
    force[2] = _externalForces[3 * x + 2];

    _translationForce += force;

    VEC3F vertex = vertices[x];
    VEC3F rotatedForce = RT * force;

    _angularForce += cross(vertex, rotatedForce);
  }

  // subtract the rigid forces from the external force
  VEC3F meanTranslation = _translationForce;
  meanTranslation *= 1.0 / unconstrained;
  VEC3F meanAngular = _angularForce;
  meanAngular *= 1.0 / unconstrained;

  for (int x = 0; x < unconstrained; x++)
  {
    VEC3F angular = R * cross(meanAngular, vertices[x]);
    VEC3F subtract = meanTranslation + angular;

    _externalForces[3 * x] -= subtract[0];
    _externalForces[3 * x + 1] -= subtract[1];
    _externalForces[3 * x + 2] -= subtract[2];
  }

  // sanity check -- compute again, see if they are zero
  VEC3F translationCheck;
  VEC3F angularCheck;
  for (int x = 0; x < unconstrained; x++)
  {
    VEC3F force;
    force[0] = _externalForces[3 * x];
    force[1] = _externalForces[3 * x + 1];
    force[2] = _externalForces[3 * x + 2];

    translationCheck += force;

    VEC3F vertex = vertices[x];
    VEC3F rotatedForce = RT * force;
    angularCheck += cross(vertex, rotatedForce);
  }

  cout << " translation check: " << translationCheck << endl;
  cout << " angular check: " << angularCheck << endl;
  cout << " angular force: " << _angularForce << endl;
  cout << " diff: " << _angularForce - angularCheck << endl;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration, using acceleration as the primary variable
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_FULLSPACE_INTEGRATOR::stepImplicitAcceleration()
{
  TIMER total;

  // allocate matrices if necessary -- should only fire
  // on the first step
  if (_CDampingSparse.rows() != _rank)
  {
    _CDampingSparse.resize(_rank, _rank);
    cacheStaticDampingSparse();

    _ASparse.setSparsity(_CDampingSparse);
  }

  TIMER preamble;
  // get the reduced mass matrix
  SPARSE_MATRIX& M = _tetMesh->TET_MESH::massMatrix();

  // get the state vector
  VECTOR& position = _tetMesh->x();

  // accumulate external forces
  //computeExternalForces();
  computeMultibodyExternalForces();
  
  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld =  _velocity;
  _accelerationOld = _acceleration;
  _timingBreakdown["Preamble"] += preamble.timing();

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  Real initialR = maxR;
  int step = 0;
  Real* alpha = _accelerationAlpha;

  Real PCGEps = pow(10.0, -_PCGDigits);
  _ASparse.eps() = PCGEps;
 
  cout << "==================================================" << endl;
  cout << " Solving multibody acceleration-level implicit step " << _totalSteps << endl;
  cout << "==================================================" << endl;
  
  while (step < _maxNewtonSteps && maxR > eps)
  {
    // update position and velocity as well
    _position = _positionOld + alpha[2] * _velocityOld +
                               alpha[3] * _accelerationOld + 
                               alpha[4] * _acceleration;
    _velocity = _velocityOld + alpha[0] * _accelerationOld +
                               alpha[1] * _acceleration;

    // update tet mesh with the new q vector
    // an assignment is overly expensive -- the "=" forces an allocation
    TIMER updateMesh;
    position.equals(_position);
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    _tetMesh->TET_MESH::updateFullMesh();

    // precompute diagonalizations --
    // this call is not timed out because timers are inside the function call
    cacheDiagonalizations();

    // get the reduced internal forces
    TIMER internalTimer;
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces(_Us, _Fhats, _Vs);
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the reduced stiffness matrix -- do this even if we are 
    // doing the conjugate gradient solve, because the damping matrix
    // needs it
    SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);
    
    // compute the damping matrix, but only if a new K is available
    TIMER RHSTimer;
    _temp = M * _acceleration;
    _residual = _CDampingSparse * _velocity;
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();
    maxR = _residual.norm2();
    if (step == 0) 
      initialR = maxR; 
    step++;

    // use the acceleration level Newmark consts
    TIMER ATimer;
    if (_ASparse.rows() == 0)
      _ASparse.setSparsity(K);
    _ASparse.equals(M);
    _ASparse.axpy(alpha[4], K);
    _ASparse.axpy(alpha[1], _CDampingSparse);
    _ASparse.useJacobi();
    _ASparse.usePCG();
    _timingBreakdown["A Assembly"] += ATimer.timing();

    TIMER solverTimer;
    cout << " residual " << step << ": " << maxR << endl;

    _solution.resizeAndWipe(_residual.size());
    _ASparse.solveCG(_solution, _residual);

    // update acceleration according to solution
    _acceleration -= _solution;

    _timingBreakdown["Solver"] += solverTimer.timing();
  }

	// update velocity
  TIMER finalUpdatePosition;
  _tetMesh->x() = _position;
  _timingBreakdown["Final Update Position"] += finalUpdatePosition.timing();
  
  // update node positions
  TIMER finalUpdateForces;
  // pushing to main
  //_tetMesh->updateSurfaceMesh();
  _tetMesh->updateFullMesh();
  _externalForces.clear();

  _timingBreakdown["Final Update Forces"] += finalUpdateForces.timing();
  _totalSteps++;
  _totalTime += total.timing();
}
