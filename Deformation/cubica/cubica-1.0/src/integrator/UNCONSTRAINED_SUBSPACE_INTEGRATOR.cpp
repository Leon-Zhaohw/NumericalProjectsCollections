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
// UNCONSTRAINED_SUBSPACE_INTEGRATOR.h: interface for the UNCONSTRAINED_SUBSPACE_INTEGRATOR class.
//
//////////////////////////////////////////////////////////////////////

#include "UNCONSTRAINED_SUBSPACE_INTEGRATOR.h"
#include <UNCONSTRAINED_SUBSPACE_TET_MESH.h>
#ifdef USING_OPENMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor for newmark integrator
//////////////////////////////////////////////////////////////////////
UNCONSTRAINED_SUBSPACE_INTEGRATOR::UNCONSTRAINED_SUBSPACE_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, Real dt, Real alpha, Real beta, Real gravity) :
  SUBSPACE_INTEGRATOR(tetMesh, dt, alpha, beta, gravity),
  _quadraticForces(3)
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)tetMesh;

  _translation = mesh->rigidTranslation();
  _translationOld = mesh->rigidTranslationOld();
  _externalRotationForces.resizeAndWipe(3);

  // compute the matrices needed for quadratic force sandwiches
  MATRIX& U = mesh->U();
  SPARSE_MATRIX& M = mesh->massMatrix();
  VECTOR uiBarO = mesh->restVector();
  MATRIX MU = M * U;
  
  // compute and store the quadratic force sandwiches
  string filename = mesh->filename();
  filename = filename + string(".quadratic.sandwiches");
  FILE* file = fopen(filename.c_str(), "rb");

  cout << " Computing quadratic force sandwiches ... ";flush(cout);
  if (file == NULL)
  {
    cout << " no cache found ... ";flush(cout);
    //fclose(file);
    file = fopen(filename.c_str(), "wb");

    SANDWICH_TRANSFORM MUi_R_Ui(MU, U);
    _MUi_R_Ui = MUi_R_Ui;
    SANDWICH_TRANSFORM MUi_R_uiBarO(MU, uiBarO);
    _MUi_R_uiBarO = MUi_R_uiBarO;

    _MUi_R_Ui.write(file);
    _MUi_R_uiBarO.write(file);
    fclose(file);
  }
  else
  {
    cout << " cache found! ... ";flush(cout);
    _MUi_R_Ui.read(file);
    _MUi_R_uiBarO.read(file);
    fclose(file);
  }
  cout << "done." << endl;

  updatePartialDelta();
}

//////////////////////////////////////////////////////////////////////
// Form the A and b so that they can be passed out for solution
// by the partitioned integrator -- implicit dynamics case
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::initializeImplicitStep()
{
  TIMER preamble;

  // accumulate external forces
  computeExternalForces();
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " EXTERNAL FORCES SPECIFIED FROM FILE " << endl;

  //cout << " external force magnitude: " << _externalForces.norm2() << endl;

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

  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
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
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::generateImplicitMatrices(map<string, double>& timingBreakdown)
{
  TIMER oldsTimer;
  // Compute transformed versions of old vectors
  VECTOR positionOldRotated(_positionOld.size());
  VECTOR velocityOldRotated(_velocityOld.size());
  VECTOR accelerationOldRotated(_accelerationOld.size());
  oldsTransformed(positionOldRotated, velocityOldRotated, accelerationOldRotated);
  timingBreakdown["Transform Olds"] += oldsTimer.timing();

  // get the reduced mass matrix
  TIMER updateMesh;
  MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // update tet mesh with the new q vector
  // an assignment is overly expensive -- the "=" forces an allocation
  q.equals(_position);
  timingBreakdown["Update Mesh"] += updateMesh.timing();

  // compute the new deformation gradient F (needed by both R and K)
  TIMER timerF;
  _tetMesh->generateF();
  timingBreakdown["Generate F"] += timerF.timing();

  // get the reduced internal forces
  TIMER internalTimer;
  _tetMesh->generateInternalForces();
  VECTOR& R = _tetMesh->reducedInternalForce();
  timingBreakdown["Internal Forces"] += internalTimer.timing();

  // get the reduced stiffness matrix -- do this even if we are 
  // doing the conjugate gradient solve, because the damping matrix
  // needs it
  TIMER stiffnessTimer;
  _tetMesh->generateStiffnessMatrix();
  MATRIX& K = _tetMesh->stiffness();
  timingBreakdown["Generate Stiffness Matrix"] += stiffnessTimer.timing();
  
  TIMER RHSTimer;
  _CDamping.clearingAxpy(_rayleighAlpha, M);
  _CDamping.axpy(_rayleighBeta, K);

  // get the positions taking into account translations
  VECTOR projectedPosition = projectedPositionPlusTranslation(_position, _translation);
  VECTOR projectedPositionOld = projectedPositionPlusTranslation(positionOldRotated, _translationOld);

  // compute the LHS of the residual:
  // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
  _temp.clearingAxpy(-_alpha[2], accelerationOldRotated);
  _temp.axpy(-_alpha[1], velocityOldRotated);
  //_temp.axpy(-_alpha[0], positionOldRotated);
  //_temp.axpy(_alpha[0], _position);
  _temp.axpy(-_alpha[0], projectedPositionOld);
  _temp.axpy(_alpha[0], projectedPosition);
  _temp = M * _temp;

  // compute the RHS of the residual:
  // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
  _totalMeshForces.clearingAxpy(_alpha[5], _accelerationOld);
  _totalMeshForces.axpy(_alpha[4], _velocityOld);
  //_totalMeshForces.axpy(-_alpha[3], _positionOld);
  //_totalMeshForces.axpy(_alpha[3], _position);
  _totalMeshForces.axpy(-_alpha[3], projectedPositionOld);
  _totalMeshForces.axpy(_alpha[3], projectedPosition);
  _totalMeshForces = _CDamping * _totalMeshForces;

  // assemble the mesh forces LHS + RHS + R
  _totalMeshForces += _temp;
  _totalMeshForces -= R;

  VECTOR fullExternal = _tetMesh->U() * _externalForces;

  for (int x = 0; x < fullExternal.size() / 3; x++)
  {
    VEC3F force;
    force[0] = fullExternal[3 * x];
    force[1] = fullExternal[3 * x + 1];
    force[2] = fullExternal[3 * x + 2];

    MATRIX3 RT = _rotation.toExplicitMatrix3x3().transpose();
    VEC3F rotatedForce = RT * force;
    fullExternal[3 * x] = rotatedForce[0];
    fullExternal[3 * x + 1] = rotatedForce[1];
    fullExternal[3 * x + 2] = rotatedForce[2];
  }
  _reducedRotatedExternals = _tetMesh->U() ^ fullExternal;

  // assemble full residual: LHS + RHS + R - F
  _residual = _totalMeshForces;
  //_residual -= _externalForces;
  _residual -= _reducedRotatedExternals;
  timingBreakdown["RHS Assembly"] += RHSTimer.timing();

  // form A
  TIMER timerA;
  _A = K;
  _A.axpy(_alpha[0], M);
  _A.axpy(_alpha[3], _CDamping);
  timingBreakdown["Form A"] += timerA.timing();
 
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << "RESIDUAL IS NOT NEGATED. " << endl;
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an implicit step -- should be called 
// right after initializeImplicitStep
//
// this has been forked out into its own function in case we want
// to do more than one Newton step
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::generateImplicitMatricesInvertible(map<string, double>& timingBreakdown)
{
  TIMER oldsTimer;
  // Compute transformed versions of old vectors
  VECTOR positionOldRotated(_positionOld.size());
  VECTOR velocityOldRotated(_velocityOld.size());
  VECTOR accelerationOldRotated(_accelerationOld.size());
  oldsTransformed(positionOldRotated, velocityOldRotated, accelerationOldRotated);
  timingBreakdown["Transform Olds"] += oldsTimer.timing();

  // get the reduced mass matrix
  MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // update tet mesh with the new q vector
  // an assignment is overly expensive -- the "=" forces an allocation
  TIMER updateMesh;
  q.equals(_position);
  timingBreakdown["Update Mesh"] += updateMesh.timing();

  // compute the new deformation gradient F (needed by both R and K)
  TIMER timerF;
  _tetMesh->generateF();
  timingBreakdown["Generate F"] += timerF.timing();

  // precompute diagonalizations --
  // this call is not timed out because timers are inside the function call
  TIMER timerCache;
  cacheDiagonalizations();
  timingBreakdown["Diagonalize F"] += timerCache.timing();

  // get the reduced internal forces
  TIMER internalTimer;
  VECTOR& R = _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);
  timingBreakdown["Internal Forces"] += internalTimer.timing();

  // get the reduced stiffness matrix -- do this even if we are 
  // doing the conjugate gradient solve, because the damping matrix
  // needs it
  MATRIX& K = _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, timingBreakdown);
  
  TIMER RHSTimer;
  _CDamping.clearingAxpy(_rayleighAlpha, M);
  _CDamping.axpy(_rayleighBeta, K);

  // get the positions taking into account translations
  VECTOR projectedPosition = projectedPositionPlusTranslation(_position, _translation);
  VECTOR projectedPositionOld = projectedPositionPlusTranslation(positionOldRotated, _translationOld);

  // compute the LHS of the residual:
  // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
  _temp.clearingAxpy(-_alpha[2], accelerationOldRotated);
  _temp.axpy(-_alpha[1], velocityOldRotated);
  _temp.axpy(-_alpha[0], projectedPositionOld);
  _temp.axpy(_alpha[0], projectedPosition);
  _temp = M * _temp;

  // compute the RHS of the residual:
  // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
  _totalMeshForces.clearingAxpy(_alpha[5], _accelerationOld);
  _totalMeshForces.axpy(_alpha[4], _velocityOld);
  _totalMeshForces.axpy(-_alpha[3], projectedPositionOld);
  _totalMeshForces.axpy(_alpha[3], projectedPosition);
  _totalMeshForces = _CDamping * _totalMeshForces;

  // assemble the mesh forces LHS + RHS + R
  _totalMeshForces += _temp;
  _totalMeshForces -= R;

  // assemble full residual: LHS + RHS + R - F
  _residual = _totalMeshForces;
  _residual -= _externalForces;
  timingBreakdown["RHS Assembly"] += RHSTimer.timing();

  // form A
  TIMER timerA;
  _A = K;
  _A.axpy(_alpha[0], M);
  _A.axpy(_alpha[3], _CDamping);
  timingBreakdown["Form A"] += timerA.timing();

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << "RESIDUAL IS NOT NEGATED. " << endl;
}

//////////////////////////////////////////////////////////////////////
// Update the mesh after the partitioned implicit step
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::finalizeImplicitStep()
{
  // Compute transformed versions of old vectors
  VECTOR positionOldRotated(_positionOld.size());
  VECTOR velocityOldRotated(_velocityOld.size());
  VECTOR accelerationOldRotated(_accelerationOld.size());
  oldsTransformed(positionOldRotated, velocityOldRotated, accelerationOldRotated);

	// update velocity
  _velocity.clearingAxpy(_alpha[5], accelerationOldRotated);
  _velocity.axpy(_alpha[4], velocityOldRotated);
  _velocity.axpy(-_alpha[3], positionOldRotated);
  _velocity.axpy(_alpha[3], _position);
  
	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], accelerationOldRotated);
  _acceleration.axpy(-_alpha[1], velocityOldRotated);
  _acceleration.axpy(-_alpha[0], positionOldRotated);
  _acceleration.axpy(_alpha[0], _position);

  /*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " DERIVATIVES DEACTIVATED" << endl;
  _velocity *= 0.0;
  _acceleration *= 0.0;
  */

  /*
  // update node positions
  //_tetMesh->updateSurfaceMesh();
  _externalForces.clear();
  _totalSteps++;
  */
}

//////////////////////////////////////////////////////////////////////
// Update the mesh after the partitioned implicit step
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::finalizeImplicitStepPlusTranslations()
{
  // Compute transformed versions of old vectors
  VECTOR positionOldRotated(_positionOld.size());
  VECTOR velocityOldRotated(_velocityOld.size());
  VECTOR accelerationOldRotated(_accelerationOld.size());
  oldsTransformed(positionOldRotated, velocityOldRotated, accelerationOldRotated);

	// update velocity
  _velocity.clearingAxpy(_alpha[5], accelerationOldRotated);
  _velocity.axpy(_alpha[4], velocityOldRotated);

  // get the positions taking into account translations
  VECTOR projectedPosition = projectedPositionPlusTranslation(_position, _translation);
  VECTOR projectedPositionOld = projectedPositionPlusTranslation(positionOldRotated, _translationOld);

  //_velocity.axpy(-_alpha[3], positionOldRotated);
  //_velocity.axpy(_alpha[3], _position);
  _velocity.axpy(-_alpha[3], projectedPositionOld);
  _velocity.axpy(_alpha[3], projectedPosition);
  
	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], accelerationOldRotated);
  _acceleration.axpy(-_alpha[1], velocityOldRotated);
  //_acceleration.axpy(-_alpha[0], positionOldRotated);
  //_acceleration.axpy(_alpha[0], _position);
  _acceleration.axpy(-_alpha[0], projectedPositionOld);
  _acceleration.axpy(_alpha[0], projectedPosition);

  // update node positions
  //_tetMesh->updateSurfaceMesh();
  _externalForces.clear();
  _totalSteps++;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::initializeQuasistaticStep()
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
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::finalizeQuasistaticStep()
{
  _externalForces.clear();
  _totalSteps++;
}

//////////////////////////////////////////////////////////////////////
// compute versions of all the old quantities transformed into 
// the new frame
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::oldsTransformed(VECTOR& positionOldRotated, 
                                                        VECTOR& velocityOldRotated, 
                                                        VECTOR& accelerationOldRotated)
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;

  VECTOR& projectedRestPose = mesh->projectedRestPose();

  MATRIX3 RT = _rotation.toExplicitMatrix3x3().transpose();
  MATRIX3 RnewRold = RT * _rotationOld.toExplicitMatrix3x3();
  MATRIX UTRU = mesh->UTRU().transform(RnewRold);
  VECTOR UTRp = mesh->UTRp().vectorTransform(RnewRold);

  VEC3F diff3 = _translationOld - _translation;
  VECTOR diff(3); diff[0] = diff3[0]; diff[1] = diff3[1]; diff[2] = diff3[2];
  MATRIX3 RnewT = RT;
  MATRIX UTRI = mesh->translationSandwich().transform(RnewT);

  positionOldRotated = UTRp;
  positionOldRotated += UTRU * _positionOld;
  positionOldRotated -= projectedRestPose;
  positionOldRotated += UTRI * diff;

  velocityOldRotated = UTRU * _velocityOld;
  accelerationOldRotated = UTRU * _accelerationOld;
}

//////////////////////////////////////////////////////////////////////
// Compute external forces with local rotations added
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::computeExternalForces()
{
  // process mouse input
  if (_clickedNode != NULL)
  {
    // forces are defined in world space
    VEC3F worldClick = *_clickedNode;
    worldClick = _rotationOld.toExplicitMatrix3x3() * worldClick + _translation;
    VEC3F dragForce = _draggedPosition - worldClick;
    dragForce.normalize();
    _forceVectors.push_back(dragForce);
    _forceNodes.push_back(_clickedNode);
  }

  // make sure some forces are in the list
  if (_forceVectors.size() == 0) return;

  // clear the force vector on the first sample
  MATRIX subbasis = _tetMesh->vertexSubBasis(_forceNodes[0]);

  // if there is no basis yet, do nothing
  if (subbasis.cols() == 0 || subbasis.rows() == 0)
    return;

  // rotate all of the force vectors with respect to the local frame
  MATRIX3 rotationT = _rotationOld.toExplicitMatrix3x3().transpose();
  VEC3F totalForce;
  for (unsigned int x = 0; x < _forceVectors.size(); x++)
  {
    //_forceVectors[x] = rotationT * _forceVectors[x];
    totalForce += _forceVectors[x];
  }

  // add the total force to the overall translation force
  _translationForce += totalForce;

  subbasis = subbasis.transpose();
  VEC3F rotated = rotationT * _forceVectors[0];
  _externalForces += subbasis.gemv(_forceMultiplier, rotated);

  // process the rest of the forces
  for (unsigned int x = 1; x < _forceVectors.size(); x++)
  {
    subbasis = _tetMesh->vertexSubBasis(_forceNodes[x]);
    subbasis = subbasis.transpose();
    //_externalForces += subbasis.gemv(_forceMultiplier, _forceVectors[x]);

    VEC3F rotated = rotationT * _forceVectors[x];
    _externalForces += subbasis.gemv(_forceMultiplier, rotated);
  }

  // compute the net torque
  _torque = 0;
  for (unsigned int x = 0; x < _forceVectors.size(); x++)
  {
    VEC3F displace = _tetMesh->vertexSubBasis(_forceNodes[x]) * _tetMesh->q();
    VEC3F vertex = _tetMesh->restPose()[x] + displace;

    // rotate force into local frame
    VEC3F localForce = rotationT * _forceVectors[x];
    localForce *= _forceMultiplier;
    _torque += cross(vertex, localForce);
  }

  _forceVectors.clear();
  _forceNodes.clear();
}

//////////////////////////////////////////////////////////////////////
// reset everything to zero
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::reset()
{
  SUBSPACE_INTEGRATOR::reset();

  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
  _translation = mesh->originalCenterOfMass();
  _translationOld = mesh->originalCenterOfMass();
  _rotation = MATRIX3::I(); 
  _rotationOld = MATRIX3::I(); 
}

//////////////////////////////////////////////////////////////////////
// compute velocity using the current _tetMesh q
//////////////////////////////////////////////////////////////////////
VECTOR UNCONSTRAINED_SUBSPACE_INTEGRATOR::computeNewestReducedVelocity()
{
  VECTOR reducedVelocity(_rank);

  VECTOR positionOldRotated(_positionOld.size());
  VECTOR velocityOldRotated(_velocityOld.size());
  VECTOR accelerationOldRotated(_accelerationOld.size());
  oldsTransformed(positionOldRotated, velocityOldRotated, accelerationOldRotated);

  // get the most recent displacement
  VECTOR& qNewest = _tetMesh->q();

  // get the positions taking into account translations
  VECTOR projectedPosition = projectedPositionPlusTranslation(qNewest, _translation);
  VECTOR projectedPositionOld = projectedPositionPlusTranslation(positionOldRotated, _translationOld);

  reducedVelocity.clearingAxpy(_alpha[5], accelerationOldRotated);
  reducedVelocity.axpy(_alpha[4], velocityOldRotated);
  reducedVelocity.axpy(-_alpha[3], projectedPositionOld);
  reducedVelocity.axpy(_alpha[3], projectedPosition);

  return reducedVelocity;
}

//////////////////////////////////////////////////////////////////////
// compute and return full coordinate velocity
//////////////////////////////////////////////////////////////////////
VECTOR UNCONSTRAINED_SUBSPACE_INTEGRATOR::fullCoordinateVelocity()
{
  VECTOR reducedVelocity(_rank);
  MATRIX& U = _tetMesh->U();

  VECTOR positionOldRotated(_positionOld.size());
  VECTOR velocityOldRotated(_velocityOld.size());
  VECTOR accelerationOldRotated(_accelerationOld.size());
  oldsTransformed(positionOldRotated, velocityOldRotated, accelerationOldRotated);

  // get the most recent displacement
  VECTOR& qNewest = _tetMesh->q();

  // get the positions taking into account translations
  VECTOR projectedPosition = projectedPositionPlusTranslation(qNewest, _translation);
  VECTOR projectedPositionOld = projectedPositionPlusTranslation(positionOldRotated, _translationOld);

  reducedVelocity.clearingAxpy(_alpha[5], accelerationOldRotated);
  reducedVelocity.axpy(_alpha[4], velocityOldRotated);
  reducedVelocity.axpy(-_alpha[3], projectedPositionOld);
  reducedVelocity.axpy(_alpha[3], projectedPosition);

  return U * reducedVelocity;
}

//////////////////////////////////////////////////////////////////////
// compute in full coordinates the position in the local frame
// taking into account the translation
//////////////////////////////////////////////////////////////////////
VECTOR UNCONSTRAINED_SUBSPACE_INTEGRATOR::projectedPositionPlusTranslation(VECTOR& position, VEC3F& translation)
{
  MATRIX& U = _tetMesh->U();
  VECTOR fullPosition = U * position;
  VEC3F rotatedTranslation = _rotation.toExplicitMatrix3x3().transpose() * translation;
  for (int x = 0; x < fullPosition.size() / 3; x++)
  {
    fullPosition[3 * x] += rotatedTranslation[0];
    fullPosition[3 * x + 1] += rotatedTranslation[1];
    fullPosition[3 * x + 2] += rotatedTranslation[2];
  }
  return U ^ fullPosition;
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an quasistatic step -- should be called 
// right after initializeQuasistaticStep
//
// this has been forked out into its own function in case we want
// to do more than one Newton step
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::generateQuasistaticMatrices(map<string, double>& timingBreakdown)
{
  VECTOR& q = _tetMesh->q();
  //TIMER updateMesh;

  // update q vector with current tet mesh state
  // an assignment is overly expensive -- the "=" forces an allocation
  q.equals(_position);
  //timingBreakdown["Update Mesh"] += updateMesh.timing();
  
  // compute the new deformation gradient F (needed by both R and K)
  //_tetMesh->generateF();
  _tetMesh->generateF();
  
  // precompute diagonalizations --
  // this call is not timed out because timers are inside the function call
  cacheDiagonalizations();
  
  //TIMER internalTimer;
  // get the reduced internal forces
  VECTOR& R = _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);
  //timingBreakdown["Internal Forces"] += internalTimer.timing();

  MATRIX& K = _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
  _A = K;

  VECTOR fullExternal = _tetMesh->U() * _externalForces;
  for (int x = 0; x < fullExternal.size() / 3; x++)
  {
    VEC3F force;
    force[0] = fullExternal[3 * x];
    force[1] = fullExternal[3 * x + 1];
    force[2] = fullExternal[3 * x + 2];

    VEC3F rotatedForce = _rotation.toExplicitMatrix3x3().transpose() * force;
    fullExternal[3 * x] = rotatedForce[0];
    fullExternal[3 * x + 1] = rotatedForce[1];
    fullExternal[3 * x + 2] = rotatedForce[2];
  }
  VECTOR reducedRotatedExternals = _tetMesh->U() ^ fullExternal;

  //TIMER RHSTimer;
  //_residual.equals(_externalForces);
  _residual.equals(reducedRotatedExternals);
  _residual += R;

  _totalMeshForces = R;
  _totalMeshForces *= -1.0;

  // account for the fact that the implicit solver does not negate
  // the gradient
  //
  // No, this is just to push the R + f_ext term to the LHS
  _residual *= -1.0;
  //timingBreakdown["RHS Assembly"] += RHSTimer.timing();
}

//////////////////////////////////////////////////////////////////////
// Update the mesh after the partitioned implicit step
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::finalizeImplicitVelocityStep()
{
  SUBSPACE_INTEGRATOR::finalizeImplicitVelocityStep();
  Real* alpha = _velocityAlpha;

  // position = position_old + alpha[3] * v_old + alpha[4] * v_new + alpha[5] * accel_old
  _angularDelta = alpha[3] * _angularVelocityOld + 
                  alpha[4] * _angularVelocity + 
                  alpha[5] * _angularAccelerationOld; 
  // acceleration = alpha[0] * v_new + alpha[1] * v_old + alpha[2] * accel_old
  _angularAcceleration = 
                       alpha[0] * _angularVelocity + 
                       alpha[1] * _angularVelocityOld + 
                       alpha[2] * _angularAccelerationOld; 
 
  // compute new translation quantities
  _translation = _translationOld + 
                 alpha[3] * _translationVelocityOld + 
                 alpha[4] * _translationVelocity + 
                 alpha[5] * _translationAccelerationOld;
  // acceleration = alpha[0] * v_new + alpha[1] * v_old + alpha[2] * accel_old
  _translationAcceleration = 
                 alpha[0] * _translationVelocity + 
                 alpha[1] * _translationVelocityOld + 
                 alpha[2] * _translationAccelerationOld;

  // position update
  QUATERNION update = QUATERNION::fromAxisAngle(_angularDelta);
  _rotation = update * _rotationOld;

  // commit to mesh
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
  mesh->rotationQuaternion() = _rotation;
  mesh->rigidTranslation() = _translation;
}

//////////////////////////////////////////////////////////////////////
// Update the mesh after the partitioned implicit step
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::finalizeImplicitAccelerationStep()
{
  SUBSPACE_INTEGRATOR::finalizeImplicitAccelerationStep();
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

  /*
  VEC3F translationUpdate = 
                alpha[2] * _translationVelocityOld + 
                alpha[3] * _translationAccelerationOld + 
                alpha[4] * _translationAcceleration;

  // pack velocities into twist matrix
  MATRIX twist(4,4);
  twist(3,3) = 1.0;
  MATRIX3 cross = MATRIX3::cross(_angularDelta);
  for (int x = 0; x < 3; x++)
  {
    for (int y = 0; y < 3; y++)
      twist(x,y) = cross(x,y);
    twist(x,3) = translationUpdate[x];
  }

  // take the exponential
  MATRIX exp = twist.exp();

  // multiply through to get the new values
  MATRIX olds(4,4);
  MATRIX3 rotationOld = _rotationOld.toExplicitMatrix3x3();
  for (int x = 0; x < 3; x++)
  {
    for (int y = 0; y < 3; y++)
      olds(x,y) = rotationOld(x,y);
    olds(x,3) = _translationOld[x];
  }
  MATRIX expUpdate = exp * olds;

  // unpack the new values
  MATRIX3 rotationNew;
  VEC3F translationNew;
  for (int x = 0; x < 3; x++)
  {
    for (int y = 0; y < 3; y++)
      rotationNew(x,y) = expUpdate(x,y);
    translationNew[x] = expUpdate(x,3);
  }
  cout << " exponential update: " << expUpdate << endl;

  QUATERNION quatNew(rotationNew);
  */

  // commit to updates
  //_rotation = quatNew;
  //_translation = translationNew;

  /*
  // DEBUG: see if the odeUpdate is right
  _angularVelocity = VEC3F(0, 0, -1.879895561357683);
  _rotationOld = QUATERNION(0,0,0,1);
  */

  // position update
  QUATERNION update = QUATERNION::fromAxisAngle(_angularDelta);

  // ONLY APPLIES TO IMPLICIT EULER!!!
  //QUATERNION update = QUATERNION::odeUpdate(_dt, _angularVelocity);

  _rotation = update * _rotationOld;

  //_rotation.normalize();

  _translation = _translationOld + 
                alpha[2] * _translationVelocityOld + 
                alpha[3] * _translationAccelerationOld + 
                alpha[4] * _translationAcceleration;

  //cout << " exp translation: " << translationNew << endl;
  //cout << " new translation: " << _translation << endl;
  //cout.precision(16);
  //cout << " exp rotation: " << quatNew << endl;
  //cout << " new rotation: " << update * _rotationOld << endl;

  /*
  VEC3F translationUpdate = alpha[2] * _translationVelocityOld + alpha[3] * _translationAccelerationOld + alpha[4] * _translationAcceleration;
  MATRIX3 dexp = MATRIX3::dexp(_angularDelta);
  MATRIX3 rotationOld = _rotationOld.toExplicitMatrix3x3();
  cout << " translation before: " << translationUpdate << endl;
  translationUpdate = dexp * translationUpdate;
  MATRIX3 hat = MATRIX3::cross(_angularDelta);
  _translation = _translationOld + translationUpdate;

  cout << " angular delta: " << _angularDelta << endl;
  cout << " translation after: " << translationUpdate << endl;
  */

  /*
  //_translation = _translationOld + _dt * MATRIX3::dexp(_dt * _angularVelocity) * _translationVelocity;
  MATRIX3 dexp = MATRIX3::dexp(_angularDelta);
  cout << " dexp: " << dexp << endl;
  dexp = _rotationOld.toExplicitMatrix3x3() * dexp;
  _translation = _translationOld + _dt * dexp * _translationVelocityOld;
  */

  // commit to mesh
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
  mesh->rotationQuaternion() = _rotation;
  mesh->rigidTranslation() = _translation;
}

//////////////////////////////////////////////////////////////////////
// update rigid body components
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::updateRigids()
{
  Real* alpha = _accelerationAlpha;

  // rotation update
  _angularVelocity.equals(_angularVelocityOld);
  _angularVelocity.axpy(alpha[0], _angularAccelerationOld);
  _angularVelocity.axpy(alpha[1], _angularAcceleration);
  //_angularVelocity = _angularVelocityOld + 
  //                     alpha[0] * _angularAccelerationOld + 
  //                     alpha[1] * _angularAcceleration;

  _angularDelta.clearingAxpy(alpha[2], _angularVelocityOld);
  _angularDelta.axpy(alpha[3], _angularAccelerationOld);
  _angularDelta.axpy(alpha[4], _angularAcceleration); 
  //_angularDelta = alpha[2] * _angularVelocityOld + 
  //                alpha[3] * _angularAccelerationOld + 
  //                alpha[4]* _angularAcceleration; 

  QUATERNION update = QUATERNION::fromAxisAngle(_angularDelta);
  _rotation = update * _rotationOld;

  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
  mesh->rotationQuaternion().equals(_rotation);

  // translation update
  //_translationVelocity = _translationVelocityOld + 
  //                 alpha[0] * _translationAccelerationOld + 
  //                 alpha[1] * _translationAcceleration;
  _translationVelocity.equals(_translationVelocityOld);
  _translationVelocity.axpy(alpha[0], _translationAccelerationOld);
  _translationVelocity.axpy(alpha[1], _translationAcceleration);

  //_translation = _translationOld + 
  //              alpha[2] * _translationVelocityOld + 
  //              alpha[3] * _translationAccelerationOld + 
  //              alpha[4] * _translationAcceleration;

  _translation.equals(_translationOld);
  _translation.axpy(alpha[2], _translationVelocityOld);
  _translation.axpy(alpha[3], _translationAccelerationOld);
  _translation.axpy(alpha[4], _translationAcceleration);

  mesh->rigidTranslation().equals(_translation);
}

//////////////////////////////////////////////////////////////////////
// Update the integrator and mesh state using and acceleration
// level update. Does not update the "olds" however -- this should
// have been done at the beginning of a time step.
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::updateStateUsingAcceleration()
{
  SUBSPACE_INTEGRATOR::updateStateUsingAcceleration();
  updateRigids();

  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
  
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " DEACTIVATING VARIABLE MASS " << endl;
  mesh->SUBSPACE_TET_MESH::refreshInertiaTensor();
  mesh->SUBSPACE_TET_MESH::refreshInertiaTensorDt(_velocity);

  updatePartialDelta();
}

//////////////////////////////////////////////////////////////////////
// Update the integrator and mesh state using a velocity
// level update. Does not update the "olds" however -- this should
// have been done at the beginning of a time step.
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::updateStateUsingVelocity()
{
  SUBSPACE_INTEGRATOR::updateStateUsingVelocity();

  Real* alpha = _velocityAlpha;

  // position = position_old + alpha[3] * v_old + alpha[4] * v_new + alpha[5] * accel_old
  _angularDelta = alpha[3] * _angularVelocityOld + 
                  alpha[4] * _angularVelocity + 
                  alpha[5] * _angularAccelerationOld; 
  // acceleration = alpha[0] * v_new + alpha[1] * v_old + alpha[2] * accel_old
  _angularAcceleration = 
                       alpha[0] * _angularVelocity + 
                       alpha[1] * _angularVelocityOld + 
                       alpha[2] * _angularAccelerationOld; 
 
  // compute new translation quantities
  _translation = _translationOld + 
                 alpha[3] * _translationVelocityOld + 
                 alpha[4] * _translationVelocity + 
                 alpha[5] * _translationAccelerationOld;
  // acceleration = alpha[0] * v_new + alpha[1] * v_old + alpha[2] * accel_old
  _translationAcceleration = 
                 alpha[0] * _translationVelocity + 
                 alpha[1] * _translationVelocityOld + 
                 alpha[2] * _translationAccelerationOld; 

  // position update
  QUATERNION update = QUATERNION::fromAxisAngle(_angularDelta);
  _rotation = update * _rotationOld;

  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
  mesh->rotationQuaternion().equals(_rotation);
  mesh->rigidTranslation().equals(_translation);
  
  mesh->SUBSPACE_TET_MESH::refreshInertiaTensor();
  mesh->SUBSPACE_TET_MESH::refreshInertiaTensorDt(_velocity);

  //updatePartialDelta();
}

//////////////////////////////////////////////////////////////////////
// compute the quadratic velocity vectors
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::computeQuadraticForces()
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;

  /////////////////////////////////////////////////////////////////////////
  // Add angular quadratics
  /////////////////////////////////////////////////////////////////////////

  // First term: -\tilde{\bar{\omega}}^i (\bar \II^i_{\theta \theta} \bar\omega^i)
  MATRIX3 omegaTilde = MATRIX3::cross(_angularVelocity);
  VECTOR quadratic = omegaTilde * mesh->inertiaTensor() * _angularVelocity.toVector();
  VECTOR QTheta = -1.0 * quadratic;

  // second term: {\dot{\bar{\II}}^i_{\theta \theta}} \bar\omega^i 
  quadratic = mesh->inertiaTensorDt() * _angularVelocity.toVector();
  QTheta -= quadratic;

  // third term: \tilde{\bar \omega}({\bar\II}^i_{\theta f} \dot \qq^i_f)
  MATRIX& U = mesh->U();
  vector<VEC3F>& vertices = mesh->vertices();
  MATRIX Ithetaf(3, U.cols());
  for (int x = 0; x < U.rows() / 3; x++)
  {
    Real mass = mesh->mass(x);
    MATRIX uBarTilde = MATRIX::cross(vertices[x]);
    SUBMATRIX submatrix(U, 3 * x, 3);

    Ithetaf += mass * uBarTilde * submatrix;
  }
  quadratic = omegaTilde * Ithetaf * _velocity;
  QTheta -= quadratic;

  _quadraticForces.clear();
  _quadraticForces.add(QTheta, 1);

  /////////////////////////////////////////////////////////////////////////
  // Add translation quadratics
  /////////////////////////////////////////////////////////////////////////
  VEC3F SitBar;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    Real mass = mesh->mass(x);
    SitBar += mass * vertices[x];
  }
  MATRIX SiBar(3, U.cols()); 
  for (int x = 0; x < U.rows() / 3; x++)
  {
    Real mass = mesh->mass(x);
    SUBMATRIX submatrix(U, 3 * x, 3);
    SiBar += mass * submatrix;
  }
  MATRIX3 omegaTilde2 = omegaTilde * omegaTilde;

  VEC3F translationFinal = omegaTilde2 * SitBar;
  translationFinal += 2.0 * omegaTilde * SiBar * _velocity;
  MATRIX3 A = _rotation.toExplicitMatrix3x3();
  translationFinal = -A * translationFinal;
  _quadraticForces.add(translationFinal.toVector(), 0);

  /////////////////////////////////////////////////////////////////////////
  // Add deformation quadratics
  /////////////////////////////////////////////////////////////////////////
  VECTOR fullVelocity = U * _velocity;
  VECTOR deformationFinal(U.cols());
  for (int x = 0; x < U.rows() / 3; x++)
  {
    Real mass = mesh->mass(x);
    VEC3F uBar = vertices[x];
    VEC3F uDotBar;
    uDotBar[0] = fullVelocity[3 * x];
    uDotBar[1] = fullVelocity[3 * x + 1];
    uDotBar[2] = fullVelocity[3 * x + 2];

    VEC3F sum = omegaTilde2 * uBar + 2.0 * omegaTilde2  * uDotBar;
    SUBMATRIX submatrix(U, 3 * x, 3);
    deformationFinal -= mass * submatrix ^ sum.toVector();
  }
  _quadraticForces.add(deformationFinal, 2);
}

//////////////////////////////////////////////////////////////////////
// compute the quadratic velocity vectors in reduced coordinates
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::computeSkinnedQuadraticForcesReduced()
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;

  /*
  /////////////////////////////////////////////////////////////////////////
  // Add angular quadratics
  /////////////////////////////////////////////////////////////////////////

  // First term: -\tilde{\bar{\omega}}^i (\bar \II^i_{\theta \theta} \bar\omega^i)
  //
  // this is already reduced, since inertiaTensor is computed in reduced
  VECTOR quadratic = omegaTilde * mesh->inertiaTensor() * _angularVelocity.toVector();
  VECTOR QTheta = -1.0 * quadratic;

  // second term: {\dot{\bar{\II}}^i_{\theta \theta}} \bar\omega^i 
  //
  // inertiaTensorDt is already computed in reduced
  quadratic = mesh->inertiaTensorDt() * _angularVelocity.toVector();
  QTheta -= quadratic;

  // third term: \tilde{\bar \omega}({\bar\II}^i_{\theta f} \dot \qq^i_f)
  //
  // rotation defo tensor already computed in reduced
  quadratic = omegaTilde * mesh->rotationDefoTensor() * _velocity;
  QTheta -= quadratic;

  _quadraticForces.clear();
  _quadraticForces.add(QTheta, 1);

  /////////////////////////////////////////////////////////////////////////
  // Add translation quadratics
  /////////////////////////////////////////////////////////////////////////
  // in explicit, SitBar doesn't get updated until updateStateUsingAcceleration(),
  // which happens after the explicit force call
  VEC3F SitBar = mesh->SitBar();
  MATRIX SiBar = mesh->SiBar();
  VEC3F translationFinal = omegaTilde2 * SitBar;
  translationFinal += 2.0 * omegaTilde * SiBar * _velocity;
  MATRIX3 A = _rotation.toExplicitMatrix3x3();
  translationFinal = -A * translationFinal;
  _quadraticForces.add(translationFinal.toVector(), 0);
  */
  MATRIX3 omegaTilde = MATRIX3::cross(_angularVelocity);
  MATRIX3 omegaTilde2 = omegaTilde * omegaTilde;

  /////////////////////////////////////////////////////////////////////////
  // Add deformation quadratics
  /////////////////////////////////////////////////////////////////////////
  //VECTOR deformationFinal = _MUi_R_Ui.transform(omegaTilde2) * (_position + 2.0 * _velocity);
  //VECTOR deformationFinal = _MUi_R_Ui.transform(omegaTilde) * (_position + 2.0 * _velocity);
  //deformationFinal += _MUi_R_uiBarO.vectorTransform(omegaTilde2);
  
  VECTOR deformationFinal = _MUi_R_Ui.transform(omegaTilde2) * (_position);
  deformationFinal += _MUi_R_Ui.transform(omegaTilde) * (2.0 * _velocity);
  deformationFinal += _MUi_R_uiBarO.vectorTransform(omegaTilde2);
  deformationFinal *= -1.0;

  _quadraticForces.add(deformationFinal, 2);
}
  
//////////////////////////////////////////////////////////////////////
// compute the quadratic velocity vectors in reduced coordinates
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::computeQuadraticForcesReduced()
{
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;

  /////////////////////////////////////////////////////////////////////////
  // Add angular quadratics
  /////////////////////////////////////////////////////////////////////////

  // First term: -\tilde{\bar{\omega}}^i (\bar \II^i_{\theta \theta} \bar\omega^i)
  //
  // this is already reduced, since inertiaTensor is computed in reduced
  MATRIX3 omegaTilde = MATRIX3::cross(_angularVelocity);
  VECTOR quadratic = omegaTilde * mesh->inertiaTensor() * _angularVelocity.toVector();
  VECTOR QTheta = -1.0 * quadratic;

  // second term: {\dot{\bar{\II}}^i_{\theta \theta}} \bar\omega^i 
  //
  // inertiaTensorDt is already computed in reduced
  quadratic = mesh->inertiaTensorDt() * _angularVelocity.toVector();
  QTheta -= quadratic;

  // third term: \tilde{\bar \omega}({\bar\II}^i_{\theta f} \dot \qq^i_f)
  //
  // rotation defo tensor already computed in reduced
  quadratic = omegaTilde * mesh->rotationDefoTensor() * _velocity;
  QTheta -= quadratic;

  _quadraticForces.clear();
  _quadraticForces.add(QTheta, 1);

  /////////////////////////////////////////////////////////////////////////
  // Add translation quadratics
  /////////////////////////////////////////////////////////////////////////
  // in explicit, SitBar doesn't get updated until updateStateUsingAcceleration(),
  // which happens after the explicit force call
  VEC3F SitBar = mesh->SitBar();
  MATRIX SiBar = mesh->SiBar();
  MATRIX3 omegaTilde2 = omegaTilde * omegaTilde;
  VEC3F translationFinal = omegaTilde2 * SitBar;
  translationFinal += 2.0 * omegaTilde * SiBar * _velocity;
  MATRIX3 A = _rotation.toExplicitMatrix3x3();
  translationFinal = -A * translationFinal;
  _quadraticForces.add(translationFinal.toVector(), 0);

  /////////////////////////////////////////////////////////////////////////
  // Add deformation quadratics
  /////////////////////////////////////////////////////////////////////////
  //VECTOR deformationFinal = _MUi_R_Ui.transform(omegaTilde2) * (_position + 2.0 * _velocity);
  VECTOR deformationFinal = _MUi_R_Ui.transform(omegaTilde) * (_position + 2.0 * _velocity);
  deformationFinal += _MUi_R_uiBarO.vectorTransform(omegaTilde2);
  deformationFinal *= -1.0;

  _quadraticForces.add(deformationFinal, 2);
}

//////////////////////////////////////////////////////////////////////
// update the exponential derivative
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::updatePartialDelta()
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
// Draw the clicked node
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::drawClickedNode()
{
  if (_clickedNode == NULL) return;

  VEC3F clickedPosition = _clickedPosition;
  //MATRIX3 rotation = _rotation.toExplicitMatrix3x3();
  //clickedPosition = rotation * clickedPosition + _translation;

  glBegin(GL_POINTS);
#ifdef SINGLE_PRECISION    
    glVertex3fv(clickedPosition);
#else
    glVertex3dv(clickedPosition);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw a line from the dragged to clicked node
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::drawForceVector()
{
  if (_clickedNode == NULL) return;
  glPointSize(10.0f);
  glColor4f(0.0f, 10.0f, 0.0f, 1.0f);

  VEC3F localPosition = *_clickedNode;
  VEC3F clickedPosition = localPosition;
  MATRIX3 rotation = _rotation.toExplicitMatrix3x3();
  clickedPosition = rotation * clickedPosition + _translation;

  VEC3F draggedPosition = _draggedPosition;
  //draggedPosition = rotation * draggedPosition + _translation;

  glBegin(GL_POINTS);
#ifdef SINGLE_PRECISION    
    glVertex3fv(clickedPosition);
    glVertex3fv(draggedPosition);
#else
    glVertex3dv(clickedPosition);
    glVertex3dv(draggedPosition);
#endif
  glEnd();

  glLineWidth(4.0f);
  glBegin(GL_LINES);
#ifdef SINGLE_PRECISION    
    glVertex3fv(clickedPosition);
    glVertex3fv(draggedPosition);
#else
    glVertex3dv(clickedPosition);
    glVertex3dv(draggedPosition);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// process a mouse click event
//////////////////////////////////////////////////////////////////////
bool UNCONSTRAINED_SUBSPACE_INTEGRATOR::click(VEC3F& point, Real maxRadius)
{
  // find and store the nearest surface point
  _clickedNode = ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh)->closestSurfaceNode(point);

  Real clickDistance = norm2(point - *_clickedNode);
  if (maxRadius > 0 && clickDistance > maxRadius * maxRadius)
  {
    // The point is too far away, negate the interaction
    _clickedNode = NULL;
    return false;
  }

  _clickedID = _tetMesh->vertexID( _clickedNode );
  //_clickedPosition = *_clickedNode;
  //_draggedPosition = *_clickedNode;
  VEC3F clicked = *_clickedNode;
  _clickedPosition = _rotation.toExplicitMatrix3x3() * clicked + _translation;
  _draggedPosition = point;

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " clicked ID: " << _clickedID << endl;

  return true;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration, using acceleration as the primary variable
//////////////////////////////////////////////////////////////////////
int UNCONSTRAINED_SUBSPACE_INTEGRATOR::stepImplicitWithImplicitCollisions(vector<pair<SURFACE*, int> >& collisions)
{
  TIMER total;
  TIMER preamble;
  // get the reduced mass matrix
  //MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // if there isn't a basis yet, do nothing
  if (q.size() == 0 || 
      (!_useKryslStiffness && _tetMesh->totalKeyTets() == 0))
    return 0;

  initializeImplicitStep();

  // compute the quadratic forces
  computeQuadraticForcesReduced();
  _quadraticForces *= 0;

  // build a list of the vertices in collision
  buildCollisionList(collisions);

  // allocate a mass matrix that includes rigid components
  BLOCK_MATRIX blockMass(3,3);
  blockMass.resizeAndWipeBlock(0,0,3,3);
  blockMass.resizeAndWipeBlock(1,1,3,3);
  int rank = _position.size();
  blockMass.resizeAndWipeBlock(2,2,rank, rank);

  BLOCK_VECTOR blockResidual(3);
  blockResidual.resizeAndWipeBlock(0,3);
  blockResidual.resizeAndWipeBlock(1,3);
  blockResidual.resizeAndWipeBlock(2,rank);

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  Real initialR = maxR;
  int step = 0;
  _timingBreakdown["Preamble"] += preamble.timing();
 
  cout << "=========================" << endl;
  cout << " Timestep " << _totalSteps << endl;
  cout << "=========================" << endl;
  cout << " total collisions: " << collisions.size() << endl;
  int maxNewtonIterations = 10;

  while (step < maxNewtonIterations && maxR > initialR * eps)
  {
    // update position and velocity as well
    TIMER updateTimer;
    updateStateUsingAcceleration();
    _timingBreakdown["Update state"] += updateTimer.timing();

    generateImplicitAccelerationMatrices(_timingBreakdown);

    TIMER residualTimer;
    blockResidual.clear();
    computeResidual(blockResidual);
    _timingBreakdown["Residual computation"] += residualTimer.timing();

    ///////////////////////////////////////////////////////////
    // full multibody solve
    ///////////////////////////////////////////////////////////
    // build the system matrix
    TIMER massTimer;
    computeMassMatrix(blockMass);
    _timingBreakdown["Mass computation"] += massTimer.timing();

    // add MA terms to residual
    BLOCK_VECTOR blockAcceleration(3);
    blockAcceleration.add(_translationAcceleration.toVector(), 0);
    blockAcceleration.add(_angularAcceleration.toVector(), 1);
    blockAcceleration.add(_acceleration, 2);
    BLOCK_VECTOR MA = blockMass * blockAcceleration;
    blockResidual += MA;

    // add collision forces to the residual
    //blockResidual += blockImplicitCollisionResiduals(collisions);
    TIMER collisionResidualTimer;
    blockResidual += blockImplicitCollisionResiduals();
    _timingBreakdown["Collision residuals"] += collisionResidualTimer.timing();
    
    // add collision forces to the Jacobian
    //TIMER collisionJacobianTimer;
    //BLOCK_MATRIX collisionJacobian = blockImplicitCollisionJacobians(_timingBreakdown);
    BLOCK_MATRIX collisionJacobian;
    blockImplicitCollisionJacobians(collisionJacobian, _timingBreakdown);
    blockMass += collisionJacobian;
    //_timingBreakdown["Collision Jacobian"] += collisionJacobianTimer.timing();

    // DEBUG: check the collision jacobian
    //if (collisions.size() > 0 && _totalSteps > 160) verifyCollisionJacobian(collisions);
    //if (collisions.size() > 0 && _totalSteps >= 50) verifyCollisionJacobian(collisions);

    // add internal force matrix afterwards -- the residual is already baked in to b
    *blockMass(2,2) += _A;

    maxR = blockResidual.norm2();
    Real transNorm = blockResidual[0].norm2();
    Real angularNorm = blockResidual[1].norm2();
    Real defoNorm = blockResidual[2].norm2();
    cout << " residual " << step << ": " << maxR << " [" << transNorm << " " << angularNorm << " " << defoNorm << "]" << endl;

    TIMER solverTimer;
    MATRIX A = blockMass.full();
    VECTOR residual = blockResidual.full();
    //bool success = A.factorCholesky();
    //A.solveCholesky(residual);
    A.factorLU();
    A.solveLU(residual);
    _timingBreakdown["Solve time"] += solverTimer.timing();

    for (int x = 0; x < 3; x++)
    {
      _translationAcceleration[x] -= residual[x];
      _angularAcceleration[x] -= residual[3 + x];
    }
    for (int x = 0; x < _acceleration.size(); x++)
      _acceleration[x] -= residual[6 + x];
    step++;
  }

  // update node positions
  TIMER finalUpdateForces;
  // pushing to main
  //_tetMesh->updateSurfaceMesh();
  finalizeImplicitAccelerationStep();

  clearForces();

  _timingBreakdown["Final Update Forces"] += finalUpdateForces.timing();
  _totalSteps++;
  _totalTime += total.timing();

  return step;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration, using acceleration as the primary variable
//////////////////////////////////////////////////////////////////////
int UNCONSTRAINED_SUBSPACE_INTEGRATOR::stepImplicitWithExplicitCollisions(vector<pair<SURFACE*, int> >& collisions)
{
  TIMER total;
  TIMER preamble;
  // get the reduced mass matrix
  //MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // if there isn't a basis yet, do nothing
  if (q.size() == 0 || 
      (!_useKryslStiffness && _tetMesh->totalKeyTets() == 0))
    return 0;

  initializeImplicitStep();

  // compute the quadratic forces
  computeQuadraticForcesReduced();
  _quadraticForces *= 0;

  // build a list of the vertices in collision
  buildCollisionList(collisions);

  // add collisions to the external force
  addExplicitCollisions();

  // allocate a mass matrix that includes rigid components
  BLOCK_MATRIX blockMass(3,3);
  blockMass.resizeAndWipeBlock(0,0,3,3);
  blockMass.resizeAndWipeBlock(1,1,3,3);
  int rank = _position.size();
  blockMass.resizeAndWipeBlock(2,2,rank, rank);

  BLOCK_VECTOR blockResidual(3);
  blockResidual.resizeAndWipeBlock(0,3);
  blockResidual.resizeAndWipeBlock(1,3);
  blockResidual.resizeAndWipeBlock(2,rank);

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  Real initialR = maxR;
  int step = 0;
  //Real* alpha = _accelerationAlpha;
  _timingBreakdown["Preamble"] += preamble.timing();
 
  cout << "=========================" << endl;
  cout << " Timestep " << _totalSteps << endl;
  cout << "=========================" << endl;
  
  while (step < 10 && maxR > initialR * eps)
  {
    // update position and velocity as well
    TIMER updateTimer;
    updateStateUsingAcceleration();
    _timingBreakdown["Update state"] += updateTimer.timing();

    generateImplicitAccelerationMatrices(_timingBreakdown);

    TIMER residualTimer;
    blockResidual.clear();
    computeResidual(blockResidual);
    _timingBreakdown["Residual computation"] += residualTimer.timing();

    ///////////////////////////////////////////////////////////
    // full multibody solve
    ///////////////////////////////////////////////////////////
    // build the system matrix
    TIMER massTimer;
    computeMassMatrix(blockMass);
    _timingBreakdown["Mass computation"] += massTimer.timing();

    // add MA terms to residual
    BLOCK_VECTOR blockAcceleration(3);
    blockAcceleration.add(_translationAcceleration.toVector(), 0);
    blockAcceleration.add(_angularAcceleration.toVector(), 1);
    blockAcceleration.add(_acceleration, 2);
    BLOCK_VECTOR MA = blockMass * blockAcceleration;
    blockResidual += MA;

    // add internal force matrix afterwards -- the residual is already baked in to b
    *blockMass(2,2) += _A;

    maxR = blockResidual.norm2();
    Real transNorm = blockResidual[0].norm2();
    Real angularNorm = blockResidual[1].norm2();
    Real defoNorm = blockResidual[2].norm2();
    cout << " residual " << step << ": " << maxR << " [" << transNorm << " " << angularNorm << " " << defoNorm << "]" << endl;

    TIMER solverTimer;
    MATRIX A = blockMass.full();
    VECTOR residual = blockResidual.full();
    A.factorCholesky();
    A.solveCholesky(residual);
    _timingBreakdown["Solve time"] += solverTimer.timing();

    for (int x = 0; x < 3; x++)
    {
      _translationAcceleration[x] -= residual[x];
      _angularAcceleration[x] -= residual[3 + x];
    }
    for (int x = 0; x < _acceleration.size(); x++)
      _acceleration[x] -= residual[6 + x];
    step++;
  }

  // update node positions
  TIMER finalUpdateForces;
  // pushing to main
  //_tetMesh->updateSurfaceMesh();
  finalizeImplicitAccelerationStep();

  clearForces();

  _timingBreakdown["Final Update Forces"] += finalUpdateForces.timing();
  _totalSteps++;
  _totalTime += total.timing();

  return step;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration, using acceleration as the primary variable
//////////////////////////////////////////////////////////////////////
int UNCONSTRAINED_SUBSPACE_INTEGRATOR::stepImplicitAcceleration()
{
  TIMER total;
  TIMER preamble;
  // get the reduced mass matrix
  //MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // if there isn't a basis yet, do nothing
  if (q.size() == 0 || 
      (!_useKryslStiffness && _tetMesh->totalKeyTets() == 0))
    return 0;

  initializeImplicitStep();

  // compute the quadratic forces
  computeQuadraticForcesReduced();
  _quadraticForces *= 0;

  // allocate a mass matrix that includes rigid components
  BLOCK_MATRIX blockMass(3,3);
  blockMass.resizeAndWipeBlock(0,0,3,3);
  blockMass.resizeAndWipeBlock(1,1,3,3);
  int rank = _position.size();
  blockMass.resizeAndWipeBlock(2,2,rank, rank);

  BLOCK_VECTOR blockResidual(3);
  blockResidual.resizeAndWipeBlock(0,3);
  blockResidual.resizeAndWipeBlock(1,3);
  blockResidual.resizeAndWipeBlock(2,rank);

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  Real initialR = maxR;
  int step = 0;
  //Real* alpha = _accelerationAlpha;
  _timingBreakdown["Preamble"] += preamble.timing();
 
  cout << "=========================" << endl;
  cout << " Timestep " << _totalSteps << endl;
  cout << "=========================" << endl;
  
  //while (step < _maxNewtonSteps && maxR > initialR * eps)
  while (step < 10 && maxR > initialR * eps)
  {
    // update position and velocity as well
    TIMER updateTimer;
    updateStateUsingAcceleration();
    _timingBreakdown["Update state"] += updateTimer.timing();

    generateImplicitAccelerationMatrices(_timingBreakdown);

    /*
    ///////////////////////////////////////////////////////////
    // defo only solve
    ///////////////////////////////////////////////////////////
    _A.axpy(1.0, M);
    _residual += M * _acceleration;
    maxR = _residual.norm2();

    TIMER solverTimer;
    
    // solve with LU factorization
    //bool success = _A.factorLU();
    bool success = _A.factorCholesky();
    cout << " residual " << step << ": " << maxR << endl;
    _A.solveCholesky(_residual);

    // update acceleration according to solution
    _acceleration -= _residual;

    _timingBreakdown["Solver"] += solverTimer.timing();
    */

    TIMER residualTimer;
    blockResidual.clear();
    computeResidual(blockResidual);
    _timingBreakdown["Residual computation"] += residualTimer.timing();

    ///////////////////////////////////////////////////////////
    // full multibody solve
    ///////////////////////////////////////////////////////////
    // build the system matrix
    TIMER massTimer;
    computeMassMatrix(blockMass);
    _timingBreakdown["Mass computation"] += massTimer.timing();

    // add MA terms to residual
    BLOCK_VECTOR blockAcceleration(3);
    blockAcceleration.add(_translationAcceleration.toVector(), 0);
    blockAcceleration.add(_angularAcceleration.toVector(), 1);
    blockAcceleration.add(_acceleration, 2);
    BLOCK_VECTOR MA = blockMass * blockAcceleration;
    blockResidual += MA;

    // add internal force matrix afterwards -- the residual is already baked in to b
    *blockMass(2,2) += _A;

    maxR = blockResidual.norm2();
    Real transNorm = blockResidual[0].norm2();
    Real angularNorm = blockResidual[1].norm2();
    Real defoNorm = blockResidual[2].norm2();
    cout << " residual " << step << ": " << maxR << " [" << transNorm << " " << angularNorm << " " << defoNorm << "]" << endl;

    TIMER solverTimer;
    MATRIX A = blockMass.full();
    VECTOR residual = blockResidual.full();
    A.factorCholesky();
    A.solveCholesky(residual);
    _timingBreakdown["Solve time"] += solverTimer.timing();

    for (int x = 0; x < 3; x++)
    {
      _translationAcceleration[x] -= residual[x];
      _angularAcceleration[x] -= residual[3 + x];
    }
    for (int x = 0; x < _acceleration.size(); x++)
      _acceleration[x] -= residual[6 + x];
    step++;
  }

  // update node positions
  TIMER finalUpdateForces;
  // pushing to main
  //_tetMesh->updateSurfaceMesh();
  finalizeImplicitAccelerationStep();

  clearForces();

  _timingBreakdown["Final Update Forces"] += finalUpdateForces.timing();
  _totalSteps++;
  _totalTime += total.timing();

  return step;
}

//////////////////////////////////////////////////////////////////////
// compute the non-linear mass matrix
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::computeMassMatrix(BLOCK_MATRIX& massMatrix)
{
  massMatrix.clear();
  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;

  // get the interia tensor in anticipation of computing rotation A
  const MATRIX& Ithetatheta = mesh->inertiaTensor();

  /////////////////////////////////////////////////////////////////////////
  // copy in rotation A
  /////////////////////////////////////////////////////////////////////////
  
  // shouldn't assume that the quadratics are implicit -- make a separate function
  // to add these in if necessary
  MATRIX rotationA = Ithetatheta;
  massMatrix.add(rotationA, 1,1);

  /////////////////////////////////////////////////////////////////////////
  // Assemble translation A
  /////////////////////////////////////////////////////////////////////////

  // get the translation mass matrix
  MATRIX translationA = mesh->totalMass() * MATRIX3::I();

  // finalize the mass matrix
  massMatrix.add(translationA, 0,0);

  /////////////////////////////////////////////////////////////////////////
  // Assemble rigid-rotation coupling terms
  /////////////////////////////////////////////////////////////////////////
  mesh->refreshSitBar();
  MATRIX3 SitBarTilde = MATRIX3::cross(mesh->SitBar());

  // get the current rotation
  MATRIX3 A = _rotation.toExplicitMatrix3x3();

  MATRIX AST = A * SitBarTilde.transpose();
  massMatrix.add(AST, 0, 1);
  MATRIX SAT = AST.transpose();
  massMatrix.add(SAT, 1, 0);

  /////////////////////////////////////////////////////////////////////////
  // Assemble deformable-translation couplings
  /////////////////////////////////////////////////////////////////////////
  MATRIX& SiBar = mesh->SiBar();
  MATRIX ASiBar = A * SiBar;
  massMatrix.add(ASiBar, 0, 2);
  MATRIX SiBarTAT = ASiBar.transpose();
  massMatrix.add(SiBarTAT, 2, 0);

  /////////////////////////////////////////////////////////////////////////
  // Assemble deformable-rotation coupling
  /////////////////////////////////////////////////////////////////////////
  MATRIX Ithetaf = mesh->rotationDefoTensor();

  massMatrix.add(Ithetaf, 1, 2);
  MATRIX IthetafT = Ithetaf.transpose();
  massMatrix.add(IthetafT, 2, 1);

  /////////////////////////////////////////////////////////////////////////
  // add deformable components
  /////////////////////////////////////////////////////////////////////////
  massMatrix.add(mesh->reducedMass(), 2,2);
}

//////////////////////////////////////////////////////////////////////
// compute the block residual, including the rigid components
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::computeResidual(BLOCK_VECTOR& blockResidual)
{
  //VECTOR& qDotDot = _acceleration;
  //VEC3F& translationDotDot = _translationAcceleration;
  //VEC3F& angularAcceleration = _angularAcceleration;

  /////////////////////////////////////////////////////////////////////////
  // Assemble rotation b
  /////////////////////////////////////////////////////////////////////////

  // compute the current center of mass
  _tetMesh->updateFullMesh();
  //vector<VEC3F>& vertices = _tetMesh->vertices();
  VEC3F SitBar = _tetMesh->SitBar();

  MATRIX A = _rotation.toExplicitMatrix3x3();
  VECTOR QTheta(3);
  VECTOR rotationForces(3);
  // TODO: do this in reduced
  /*
  VECTOR fullExternals(_tetMesh->U().rows());
  fullExternals = _tetMesh->U() * _externalForces;
  for (unsigned int y = 0; y < vertices.size(); y++)
  {
    // was updated with most recent deformation above
    VEC3F& vertex = vertices[y];
    VEC3F centered = vertex - SitBar;
    MATRIX uBar = MATRIX::cross(centered);

    MATRIX AuBar = A * uBar;
    MATRIX AuBarT = AuBar.transpose();

    VECTOR F(3);

    // TODO: put rotation forces here
    F(0) = fullExternals[3 * y];
    F(1) = fullExternals[3 * y + 1];
    F(2) = fullExternals[3 * y + 2];

    F += _translationForce.toVector();
    rotationForces += -1.0 * AuBarT * F;
  }
  */

  // get the spring forces
  vector<int> blockSizes;
  blockSizes.push_back(3);
  blockSizes.push_back(3);
  blockSizes.push_back(_tetMesh->rank());
  QTheta = rotationForces;

  BLOCK_VECTOR& quadratics = _quadraticForces;
  QTheta += quadratics[1];

  // compute the rotational residual
  blockResidual.set(-1.0 * QTheta, 1);
  blockResidual.subtract(_externalRotationForces, 1);

  /////////////////////////////////////////////////////////////////////////
  // Assemble translation b
  /////////////////////////////////////////////////////////////////////////
  VEC3F translationForce = _translationForce;
  VEC3F b;
  b -= translationForce;
  b -= quadratics[0];

  // compute the translation residual
  blockResidual.set(b.toVector(), 0);

  /////////////////////////////////////////////////////////////////////////
  // Assemble deformation b
  /////////////////////////////////////////////////////////////////////////
  blockResidual.set(_residual, 2);
  blockResidual.subtract(quadratics[2], 2);
}

//////////////////////////////////////////////////////////////////////
// add implicit collision forces
//////////////////////////////////////////////////////////////////////
BLOCK_VECTOR UNCONSTRAINED_SUBSPACE_INTEGRATOR::blockImplicitCollisionResidualsDebug()
{
  BLOCK_VECTOR blockResidual(3);
  blockResidual.resizeAndWipeBlock(0,3);
  blockResidual.resizeAndWipeBlock(1,3);
  blockResidual.resizeAndWipeBlock(2,_position.size());

  // cache the rotations
  MATRIX3 R = ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh)->rotationQuaternion().toExplicitMatrix3x3();
  MATRIX3 RT = R.transpose();
  VEC3F& translation = _translation;

  // get the force for each point
  VEC3F summedForces;
  vector<VEC3F> forces;
  VECTOR rotationForce(3);
  VEC3F SitBar = _tetMesh->SitBar();
  for (unsigned int x = 0; x < _collisionVertices.size(); x++)
  {
    VEC3F* vertex = _collisionVertices[x].first;
    SURFACE* plane = _collisionVertices[x].second;
    int vertexID = _tetMesh->vertexID(vertex);

    // retrieve the velocity for damping
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);
    VEC3F vertexVelocity = VEC3F(vertexSubBasis * _velocity);
    vertexVelocity = R * vertexVelocity;

    VEC3F& restVertex= *(_tetMesh->restVertices(vertexID));
    VECTOR displace = vertexSubBasis * _position;
    VEC3F updatedPosition = restVertex + VEC3F(displace);
    updatedPosition = R * updatedPosition + translation;

    // retrieve the damped force
    VEC3F force = plane->force(updatedPosition, vertexVelocity);
    summedForces += force;

    // store the forces so we don't have to recompute later
    forces.push_back(force);

    // compute the angular term
    VEC3F localPosition = restVertex + VEC3F(displace);
    VEC3F centered = localPosition - SitBar;
    MATRIX uBarTilde = MATRIX::cross(centered);
    MATRIX AuBar = R * uBarTilde;
    MATRIX AuBarT = AuBar.transpose();

    VECTOR F(3);
    F(0) = force[0];
    F(1) = force[1];
    F(2) = force[2];

    rotationForce += -1.0 * AuBarT * F;
  }
  blockResidual[0] -= summedForces.toVector();
  blockResidual[1] -= rotationForce;

  // get the defo force
  for (unsigned int x = 0; x < _collisionVertices.size(); x++)
  {
    //VEC3F* vertex = iter->first.first;
    VEC3F* vertex = _collisionVertices[x].first;
    //int vertexID = _tetMesh->vertexID(vertex);

    // retrieve the damped force
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);
    VEC3F force = forces[x];

    // rotate into the local frame
    force = RT * force;
    
    VECTOR projectedForce = vertexSubBasis ^ force.toVector();
    blockResidual[2] -= projectedForce;
  }

  return blockResidual;  
}

//////////////////////////////////////////////////////////////////////
// add implicit collision forces
//////////////////////////////////////////////////////////////////////
//BLOCK_VECTOR UNCONSTRAINED_SUBSPACE_INTEGRATOR::blockImplicitCollisionResiduals(vector<pair<SURFACE*, int> >& collisions)
BLOCK_VECTOR UNCONSTRAINED_SUBSPACE_INTEGRATOR::blockImplicitCollisionResiduals()
{
  BLOCK_VECTOR blockResidual(3);
  blockResidual.resizeAndWipeBlock(0,3);
  blockResidual.resizeAndWipeBlock(1,3);
  blockResidual.resizeAndWipeBlock(2,_position.size());

  // clear the previous collision force cache
  _collisionForces.clear();

  // cache the rotations
  MATRIX3 R = ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh)->rotationQuaternion().toExplicitMatrix3x3();
  MATRIX3 RT = R.transpose();
  VEC3F& translation = _translation;
  MATRIX3 crossAngularVelocity = MATRIX3::cross(_angularVelocity);

  // get the force for each point
  VEC3F summedForces;
  VECTOR rotationForce(3);
  VEC3F SitBar = _tetMesh->SitBar();
  for (unsigned int x = 0; x < _collisionVertices.size(); x++)
  {
    VEC3F* vertex = _collisionVertices[x].first;
    SURFACE* surface = _collisionVertices[x].second;
    int vertexID = _tetMesh->vertexID(vertex);

    // NEW: get the surface area of the collision vertex
    Real surfaceArea = _tetMesh->surfaceArea(vertexID);

    // retrieve the velocity for damping
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);

    VEC3F& restVertex= *(_tetMesh->restVertices(vertexID));
    VECTOR displace = vertexSubBasis * _position;
    VEC3F localPosition = restVertex + VEC3F(displace);
    VEC3F updatedPosition = R * localPosition + translation;

    VEC3F vertexVelocity = VEC3F(vertexSubBasis * _velocity);
    vertexVelocity = R * vertexVelocity;
    vertexVelocity += _translationVelocity;
    vertexVelocity += R * crossAngularVelocity * localPosition;

    // retrieve the damped force
    VEC3F force = surface->force(updatedPosition, vertexVelocity);

    // cache the unscaled force
    _collisionForces.push_back(force);

    // NEW: weight the force by the surface area
    force *= surfaceArea;
    summedForces += force;

    // compute the angular term
    VEC3F centered = localPosition - SitBar;
    MATRIX uBarTilde = MATRIX::cross(centered);
    MATRIX AuBar = R * uBarTilde;
    MATRIX AuBarT = AuBar.transpose();

    VECTOR F(3);
    F(0) = force[0];
    F(1) = force[1];
    F(2) = force[2];

    rotationForce += -1.0 * AuBarT * F;

    // compute the defo term
    VEC3F rotatedForce = RT * force;
    VECTOR projectedForce = vertexSubBasis ^ rotatedForce.toVector();
    blockResidual[2] -= projectedForce;
  }
  blockResidual[0] -= summedForces.toVector();
  blockResidual[1] -= rotationForce;

  return blockResidual;  
}

//////////////////////////////////////////////////////////////////////
// add explicit collision forces
//////////////////////////////////////////////////////////////////////
//void UNCONSTRAINED_SUBSPACE_INTEGRATOR::addExplicitCollisions(vector<pair<SURFACE*, int> >& collisions)
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::addExplicitCollisions()
{
  // cache the rotations
  MATRIX3 R = ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh)->rotationQuaternionOld().toExplicitMatrix3x3();
  MATRIX3 RT = R.transpose();

  /*
  // build a list of the unique vertices in collision
  vector<TRIANGLE*> faces = _tetMesh->explicitSurfaceFaces();
  map<pair<VEC3F*, SURFACE*>, bool> vertexHash;
  for (int x = 0; x < collisions.size(); x++)
  {
    int index = collisions[x].second;
    SURFACE* plane = collisions[x].first;

    TRIANGLE* face = faces[index];
    vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(0), plane)] = true;
    vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(1), plane)] = true;
    vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(2), plane)] = true;
  }
  */

  // get the force for each point and get the mean
  VEC3F summedForces;
  vector<VEC3F> forces;
  //map<pair<VEC3F*, SURFACE*>, bool>::iterator iter;
  //for (iter = vertexHash.begin(); iter != vertexHash.end(); iter++)
  for (unsigned int x = 0; x < _collisionVertices.size(); x++)
  {
    //VEC3F* vertex = iter->first.first;
    //SURFACE* plane = iter->first.second;
    VEC3F* vertex = _collisionVertices[x].first;
    SURFACE* plane = _collisionVertices[x].second;
    //int vertexID = _tetMesh->vertexID(vertex);

    // retrieve the velocity for damping
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);
    VEC3F vertexVelocity = VEC3F(vertexSubBasis * _velocity);
    vertexVelocity = R * vertexVelocity;

    // retrieve the damped force
    VEC3F force = plane->force(*vertex, vertexVelocity);
    summedForces += force;

    // store the forces so we don't have to recompute later
    forces.push_back(force);
  }
  _translationForce += summedForces;

  // compute the rotation component
  VEC3F SitBar = _tetMesh->SitBar();
  VECTOR rotationForce(3);
  //for (iter = vertexHash.begin(); iter != vertexHash.end(); iter++, x++)
  for (unsigned int x = 0; x < _collisionVertices.size(); x++)
  {
    //VEC3F& vertex = *(iter->first.first);
    VEC3F& vertex = *(_collisionVertices[x].first);
    VEC3F centered = vertex - SitBar;
    MATRIX uBarTilde = MATRIX::cross(centered);

    MATRIX AuBar = R * uBarTilde;
    MATRIX AuBarT = AuBar.transpose();

    VECTOR F(3);
    F(0) = forces[x][0];
    F(1) = forces[x][1];
    F(2) = forces[x][2];

    rotationForce += -1.0 * AuBarT * F;
  }
  _externalRotationForces = rotationForce;

  // get the defo force, but remember to subtract out the translation and rotation components
  //for (iter = vertexHash.begin(), x = 0; iter != vertexHash.end(); iter++, x++)
  for (unsigned int x = 0; x < _collisionVertices.size(); x++)
  {
    //VEC3F* vertex = iter->first.first;
    VEC3F* vertex = _collisionVertices[x].first;
    //int vertexID = _tetMesh->vertexID(vertex);

    // retrieve the damped force
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);
    VEC3F force = forces[x];

    // rotate into the local frame
    force = RT * force;
    
    VECTOR projectedForce = vertexSubBasis ^ force.toVector();
    _externalForces += projectedForce;
  }
}

/*
//////////////////////////////////////////////////////////////////////
// add implicit collision forces directly to the residual
//////////////////////////////////////////////////////////////////////
BLOCK_MATRIX UNCONSTRAINED_SUBSPACE_INTEGRATOR::blockImplicitCollisionJacobiansDebug()
{
  BLOCK_MATRIX systemMatrix(3,3);
  int rank = _position.size();
  systemMatrix.resizeAndWipeBlock(0,0,3,3);
  systemMatrix.resizeAndWipeBlock(1,1,3,3);
  systemMatrix.resizeAndWipeBlock(2,2,rank,rank);

  systemMatrix.resizeAndWipeBlock(0,1,3,3);
  systemMatrix.resizeAndWipeBlock(1,0,3,3);
  systemMatrix.resizeAndWipeBlock(0,2,3,rank);
  systemMatrix.resizeAndWipeBlock(1,2,3,rank);
  systemMatrix.resizeAndWipeBlock(2,0,rank,3);
  systemMatrix.resizeAndWipeBlock(2,1,rank,3);

  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
  Real* alpha = _accelerationAlpha;
  MATRIX rotationOld  = _rotationOld.toExplicitMatrix3x3(); 
  VEC3F& angularDelta = _angularDelta;
  QUATERNION update = QUATERNION::fromAxisAngle(angularDelta);
  MATRIX3 leftRotation = _rotation.toExplicitMatrix3x3();

  MATRIX3 R = mesh->rotationQuaternion().toExplicitMatrix3x3();
  MATRIX3 RT = R.transpose();
  VEC3F translation = mesh->rigidTranslation();

  // build the partial of Sit
  MATRIX SitPartial(3, rank);
  MATRIX& U = mesh->U();
  for (int x = 0; x < U.rows() / 3; x++)
  {
    SUBMATRIX Uij(U, 3 * x, 3);
    Real mass = mesh->mass(x);
    SitPartial += mass * alpha[4] * Uij;
  }
  SitPartial *= 1.0 / mesh->totalMass();

  TENSOR3& partialDelta = _partialDelta;
  TENSOR3 partialRotation = partialDelta * rotationOld;

  // compute the center of mass
  VEC3F SitBar = mesh->SitBar();

  // partial with respect to the spring terms
  MATRIX springPartials(3,3);
  MATRIX springPartialsFiniteDiff(3,3);

  // partial with respect to the cross product term
  MATRIX crossPartials(3,3);

  MATRIX translationPartialTranslation(3,3);
  MATRIX translationPartialAngular(3,3);
  MATRIX angularPartialTranslation(3,3);
  MATRIX translationPartialDefo(3, rank);
  MATRIX translationPartialDefoDamping(3, rank);
  MATRIX defoPartialAngular(rank, 3);
  MATRIX angularPartialDefo(3, rank);
  MATRIX angularPartialAngular(3,3);
  MATRIX defoPartialDefo(rank,rank);

  // get all the spring forces
  for (unsigned int x = 0; x < _collisionVertices.size(); x++)
  {
    VEC3F* vertex = _collisionVertices[x].first;
    SURFACE* surface = _collisionVertices[x].second;
    int vertexID = _tetMesh->vertexID(vertex);

    // retrieve the velocity for damping
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);
    VEC3F localVelocity = VEC3F(vertexSubBasis * _velocity);
    VEC3F vertexVelocity = R * localVelocity;

    // retrieve position for stiffness
    VEC3F& restVertex= *(_tetMesh->restVertices(vertexID));
    VECTOR displace = vertexSubBasis * _position;
    VEC3F localPosition = restVertex + VEC3F(displace);
    VEC3F updatedPosition = R * localPosition + translation;

    // if the solves have pushed this vertex outside, its jacobian is zero
    if (!surface->inside(updatedPosition)) continue;

    // force partial collision point coordinates
    MATRIX subbasis = ((SUBSPACE_TET_MESH*)_tetMesh)->vertexSubBasis(vertex);
    MATRIX springJacobian = surface->springJacobian(updatedPosition);
    springJacobian *= _accelerationAlpha[4];

    // add to easiest term, translation wrt itself
    translationPartialTranslation += springJacobian;

    // defo-defo
    defoPartialDefo += _accelerationAlpha[4] * surface->defoSpringJacobian(updatedPosition, subbasis, R, translation);
    defoPartialDefo += _accelerationAlpha[1] * surface->defoDampingJacobian(updatedPosition, subbasis, R);

    // angular-angular
    angularPartialAngular += surface->angularSpringJacobian(updatedPosition, vertexVelocity, localPosition, R, partialRotation);

    // trans-defo and transpose
    MATRIX currentTransPartialDefo = surface->translationDefoJacobian(subbasis, R);
    currentTransPartialDefo *= _accelerationAlpha[4];
    translationPartialDefo += currentTransPartialDefo;

    MATRIX currentTransPartialDefoDamping = surface->translationDefoDampingJacobian(subbasis, R);
    currentTransPartialDefoDamping *= _accelerationAlpha[1];
    translationPartialDefoDamping += currentTransPartialDefoDamping;
    
    // trans-angular and transpose
    translationPartialAngular += surface->translationAngularJacobian(localPosition, partialRotation);
    translationPartialAngular += surface->translationAngularDampingJacobian(localVelocity, partialRotation);
    angularPartialTranslation += _accelerationAlpha[4] * surface->angularTranslationJacobian(localPosition, R);

    // defo-angular
    defoPartialAngular += surface->defoAngularJacobian(updatedPosition, vertexVelocity, localPosition, R, partialRotation, subbasis);

    // angular-defo
    MATRIX currentAngularPartialDefo = surface->angularDefoJacobian(updatedPosition, vertexVelocity, localPosition, subbasis, R);
    currentAngularPartialDefo *= _accelerationAlpha[4];
    angularPartialDefo += currentAngularPartialDefo;
    
    MATRIX currentAngularPartialDefoDamping = surface->angularDefoJacobianDamping(updatedPosition, vertexVelocity, localPosition, subbasis, R);
    currentAngularPartialDefoDamping *= _accelerationAlpha[1];
    angularPartialDefo += currentAngularPartialDefoDamping;
  }

  // cache the symmetric part
  MATRIX defoPartialTranslation = translationPartialDefo.transpose();

  // add in the asymmetric part
  translationPartialDefo += translationPartialDefoDamping;

  // add translation diagonal term
  systemMatrix.add(translationPartialTranslation, 0,0);
  systemMatrix.add(angularPartialTranslation, 1,0);
  systemMatrix.add(defoPartialTranslation, 2,0);

  systemMatrix.add(translationPartialAngular, 0,1);
  systemMatrix.add(angularPartialAngular, 1,1);
  systemMatrix.add(defoPartialAngular, 2,1);

  systemMatrix.add(angularPartialDefo, 1,2);
  systemMatrix.add(translationPartialDefo, 0,2);
  systemMatrix.add(defoPartialDefo, 2,2);

  return systemMatrix;
}
*/

//////////////////////////////////////////////////////////////////////
// add implicit collision forces directly to the residual
//////////////////////////////////////////////////////////////////////
//BLOCK_MATRIX UNCONSTRAINED_SUBSPACE_INTEGRATOR::blockImplicitCollisionJacobians(map<string, double>& timingBreakdown)
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::blockImplicitCollisionJacobians(BLOCK_MATRIX& systemMatrix, map<string, double>& timingBreakdown)
{
  TIMER jacobianPreamble;
  systemMatrix.resizeAndWipe(3,3);
  int rank = _position.size();
  systemMatrix.resizeAndWipeBlock(0,0,3,3);
  systemMatrix.resizeAndWipeBlock(1,1,3,3);
  systemMatrix.resizeAndWipeBlock(2,2,rank,rank);

  systemMatrix.resizeAndWipeBlock(0,1,3,3);
  systemMatrix.resizeAndWipeBlock(1,0,3,3);
  systemMatrix.resizeAndWipeBlock(0,2,3,rank);
  systemMatrix.resizeAndWipeBlock(1,2,3,rank);
  systemMatrix.resizeAndWipeBlock(2,0,rank,3);
  systemMatrix.resizeAndWipeBlock(2,1,rank,3);

  // cache a single entry so it can be area-scaled
  BLOCK_MATRIX singleEntry(systemMatrix);

  UNCONSTRAINED_SUBSPACE_TET_MESH* mesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh;
  Real* alpha = _accelerationAlpha;
  MATRIX rotationOld  = _rotationOld.toExplicitMatrix3x3(); 
  VEC3F& angularDelta = _angularDelta;
  QUATERNION update = QUATERNION::fromAxisAngle(angularDelta);
  MATRIX3 leftRotation = _rotation.toExplicitMatrix3x3();

  MATRIX3 R = mesh->rotationQuaternion().toExplicitMatrix3x3();
  MATRIX3 RT = R.transpose();
  VEC3F translation = mesh->rigidTranslation();
  timingBreakdown["Collision Jacobian preamble"] += jacobianPreamble.timing();

  TIMER partialDeltaTimer;
  TENSOR3& partialDelta = _partialDelta;
  TENSOR3 partialRotation = partialDelta * rotationOld;
  timingBreakdown["Collision Partial Delta"] += partialDeltaTimer.timing();

  // set the surface consts
  for (unsigned int x = 0; x < _collisionSurfaces.size(); x++)
  {
    _collisionSurfaces[x]->collisionAlpha() = _accelerationAlpha;
    _collisionSurfaces[x]->collisionTranslation() = translation;
    _collisionSurfaces[x]->collisionRotation() = R;
    _collisionSurfaces[x]->collisionRotationPartial() = partialRotation;
  }

  // get all the spring forces
  MATRIX3 crossAngularVelocity = MATRIX3::cross(_angularVelocity);
  for (unsigned int x = 0; x < _collisionVertices.size(); x++)
  {
    TIMER collisionSetupTimer;
    VEC3F* vertex = _collisionVertices[x].first;
    SURFACE* surface = _collisionVertices[x].second;
    int vertexID = _tetMesh->vertexID(vertex);

    // retrieve the velocity for damping
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);
    
    // retrieve position for stiffness
    VEC3F& restVertex= *(_tetMesh->restVertices(vertexID));
    VECTOR displace = vertexSubBasis * _position;
    VEC3F localPosition = restVertex + VEC3F(displace);
    VEC3F updatedPosition = R * localPosition + translation;

    VEC3F localVelocity = VEC3F(vertexSubBasis * _velocity);
    VEC3F vertexVelocity = R * localVelocity;
    vertexVelocity += _translationVelocity;
    vertexVelocity += R * crossAngularVelocity * localPosition;

    timingBreakdown["Collision Jacobian setup"] += collisionSetupTimer.timing();

    // if the solves have pushed this vertex outside, its jacobian is zero
    if (surface->inside(updatedPosition)) 
    {
      // stomp the previous entry
      TIMER insidePreamble;
      singleEntry.clear();
      MATRIX subbasis = ((SUBSPACE_TET_MESH*)_tetMesh)->vertexSubBasis(vertex);
      timingBreakdown["Inside preamble"] += insidePreamble.timing();

#if 1
      TIMER surfaceTimer;
      surface->springJacobian(updatedPosition, 
                              vertexVelocity, 
                              _collisionForces[x], 
                              localPosition,
                              subbasis,
                              singleEntry);
      timingBreakdown["Collision Surface Call"] += surfaceTimer.timing();
#else    
      surface->springJacobianDebug(updatedPosition, 
                              vertexVelocity,
                              _collisionForces[x], 
                              localPosition,
                              subbasis,
                              singleEntry, timingBreakdown);
#endif

      // NEW: area weight the collision
      TIMER collisionCommitTimer;
      Real surfaceArea = _tetMesh->surfaceArea(vertexID);
      singleEntry *= surfaceArea;
      systemMatrix += singleEntry;
      timingBreakdown["Collision commit"] += collisionCommitTimer.timing();
    }
  }
}
/*
{
  // cache the rotations
  MATRIX3 R = ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh)->rotationQuaternion().toExplicitMatrix3x3();
  MATRIX3 RT = R.transpose();
  VEC3F translation = ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh)->rigidTranslation();

  // build a list of the unique vertices in collision
  vector<TRIANGLE*> faces = _tetMesh->explicitSurfaceFaces();
  map<pair<VEC3F*, SURFACE*>, bool> vertexHash;
  for (int x = 0; x < collisions.size(); x++)
  {
    int index = collisions[x].second;
    SURFACE* surface = collisions[x].first;

    TRIANGLE* face = faces[index];
    vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(0), surface)] = true;
    vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(1), surface)] = true;
    vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(2), surface)] = true;
  }

  // get the Jacobian for each point and add it to the block Jacobian
  MATRIX translationJacobian(3,3);
  map<pair<VEC3F*, SURFACE*>, bool>::iterator iter;
  for (iter = vertexHash.begin(); iter != vertexHash.end(); iter++)
  {
    VEC3F* vertex = iter->first.first;
    SURFACE* surface = iter->first.second;
    int vertexID = _tetMesh->vertexID(vertex);

    // retrieve the velocity for damping
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);
    VEC3F vertexVelocity = VEC3F(vertexSubBasis * _velocity);

    // retrieve position for stiffness
    VEC3F& restVertex= *(_tetMesh->restVertices(vertexID));
    VECTOR displace = vertexSubBasis * _position;
    VEC3F updatedPosition = restVertex + VEC3F(displace);
    updatedPosition = R * updatedPosition + translation;

    MATRIX springJacobian = surface->springJacobian(updatedPosition);
    MATRIX dampingJacobian = surface->dampingJacobian(updatedPosition, vertexVelocity);
    springJacobian *= _accelerationAlpha[4];
    dampingJacobian *= _accelerationAlpha[1];
    MATRIX jacobian = springJacobian - dampingJacobian;

    MATRIX defoJacobian = R * jacobian;
    translationJacobian += jacobian;

    MATRIX reducedJacobian = vertexSubBasis ^ defoJacobian * vertexSubBasis;
    *blockMass(2,2) += reducedJacobian;
  }
  *blockMass(0,0) += translationJacobian;

  // do the rotation Jacobian, which uses the force directly instead of the force Jacobian
  int x = 0;
  MATRIX rotationJacobian(3,3);
  VEC3F SitBar = _tetMesh->SitBar();
  MATRIX rotationOld = ((UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh)->rotationQuaternionOld().toExplicitMatrix3x3();
  for (iter = vertexHash.begin(), x = 0; iter != vertexHash.end(); iter++, x++)
  {
    VEC3F* vertex = iter->first.first;
    SURFACE* plane = iter->first.second;
    int vertexID = _tetMesh->vertexID(vertex);

    // retrieve the velocity for damping
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);
    VEC3F vertexVelocity = VEC3F(vertexSubBasis * _velocity);

    VEC3F force = plane->force(*vertex, vertexVelocity);
    VEC3F centered = *vertex - SitBar;
    MATRIX uBarTilde = MATRIX::cross(centered);

    TENSOR3 partialDeltaRuTilde = _partialDelta * rotationOld * uBarTilde;
    TENSOR3 partialT = partialDeltaRuTilde.transpose();

    MATRIX final = partialT.modeOneProduct(force.toVector());
    rotationJacobian -= final;
  }
  if (vertexHash.size() > 0)
  {
    //rotationJacobian *= 1.0 / vertexHash.size();
    *blockMass(1,1) += rotationJacobian;
  }
}
*/

//////////////////////////////////////////////////////////////////////
// verify the correctness of the collision Jacobian
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::verifyCollisionJacobian(vector<pair<SURFACE*, int> >& collisions)
{
  cout << " total collisions: " << _collisionVertices.size() << endl;

  // wipe all vertices but the first one
  pair<VEC3F*, SURFACE*> cached = _collisionVertices[0];
  _collisionVertices.clear();
  _collisionVertices.push_back(cached);

  //BLOCK_VECTOR originalResidual = blockImplicitCollisionResidualsDebug();
  BLOCK_VECTOR originalResidual = blockImplicitCollisionResiduals();

  VEC3F angularAccelerationOriginal = _angularAcceleration;
  VEC3F translationAccelerationOriginal = _translationAcceleration;
  VECTOR accelerationOriginal = _acceleration;

  //BLOCK_MATRIX computedJacobian = blockImplicitCollisionJacobians(_timingBreakdown);
  BLOCK_MATRIX computedJacobian;
  blockImplicitCollisionJacobians(computedJacobian, _timingBreakdown);
  SUBSPACE_TET_MESH* mesh = (SUBSPACE_TET_MESH*)_tetMesh;

  vector<VECTOR> differences;
  //Real delta = pow(10.0f, -7.0f);
  //Real delta = pow(10.0f, -6.0f);
  Real delta = pow(10.0f, -5.0f);
  //Real delta = pow(10.0f, -4.0f);
  int rank = _position.size();

  // step through translation components
  for (int x = 0; x < 3; x++)
  {
    // perturb the angular component
    _translationAcceleration = translationAccelerationOriginal;
    _translationAcceleration[x] += delta;
    updateStateUsingAcceleration();
    mesh->updateFullMesh();

    // get force at the perturbation
    BLOCK_VECTOR newForceVector(3);
    newForceVector.resizeAndWipeBlock(0, 3);
    newForceVector.resizeAndWipeBlock(1, 3);
    newForceVector.resizeAndWipeBlock(2, rank);
    //newForceVector = blockImplicitCollisionResidualsDebug();
    newForceVector = blockImplicitCollisionResiduals();

    newForceVector -= originalResidual;
    newForceVector *= 1.0 / delta;

    differences.push_back(newForceVector.full());
  }
  _translationAcceleration = translationAccelerationOriginal;
  updateStateUsingAcceleration();
  mesh->updateFullMesh();
  
  // step through angular components
  for (int x = 0; x < 3; x++)
  {
    // perturb the angular component
    _angularAcceleration = angularAccelerationOriginal;
    _angularAcceleration[x] += delta;
    updateStateUsingAcceleration();
    mesh->updateFullMesh();

    // get force at the perturbation
    BLOCK_VECTOR newForceVector(3);
    newForceVector.resizeAndWipeBlock(0, 3);
    newForceVector.resizeAndWipeBlock(1, 3);
    newForceVector.resizeAndWipeBlock(2, rank);
    //newForceVector = blockImplicitCollisionResidualsDebug();
    newForceVector = blockImplicitCollisionResiduals();

    newForceVector -= originalResidual;
    newForceVector *= 1.0 / delta;

    differences.push_back(newForceVector.full());
  }
  _angularAcceleration = angularAccelerationOriginal;
  updateStateUsingAcceleration();
  mesh->updateFullMesh();

  // step through defo components
  for (int x = 0; x < rank; x++)
  {
    // perturb the angular component
    _acceleration = accelerationOriginal;
    _acceleration[x] += delta;
    updateStateUsingAcceleration();
    mesh->updateFullMesh();

    // get force at the perturbation
    BLOCK_VECTOR newForceVector(3);
    newForceVector.resizeAndWipeBlock(0, 3);
    newForceVector.resizeAndWipeBlock(1, 3);
    newForceVector.resizeAndWipeBlock(2, rank);
    //newForceVector = blockImplicitCollisionResidualsDebug();
    newForceVector = blockImplicitCollisionResiduals();

    newForceVector -= originalResidual;
    newForceVector *= 1.0 / delta;

    differences.push_back(newForceVector.full());
  }
  _acceleration = accelerationOriginal;
  updateStateUsingAcceleration();
  mesh->updateFullMesh();

  MATRIX finiteDiffs(differences);

  cout << " computed: " << computedJacobian.full() << endl;
  cout << " finite diff: " << finiteDiffs << endl;

  MATRIX diff = computedJacobian.full() - finiteDiffs;
  cout << " diff: " << diff << endl;

  MATRIX relative = diff;
  for (int x = 0; x < relative.rows(); x++)
    for (int y = 0; y < relative.cols(); y++)
    {
      if (fabs(finiteDiffs(x,y)) > 0.0)
        relative(x,y) *= 1.0 / finiteDiffs(x,y);
      relative(x,y) = fabs(relative(x,y));
    }
  cout << " relative: " << relative << endl;
  cout << " max absolute entry: " << relative.maxAbsEntry() << endl;

  MATRIX ratios = finiteDiffs;
  MATRIX full = computedJacobian.full();
  for (int x = 0; x < relative.rows(); x++)
    for (int y = 0; y < relative.cols(); y++)
    {
      if (fabs(full(x,y)) > 0.0)
        ratios(x,y) *= 1.0 / full(x,y);
    }
  cout << " ratios: " << ratios << endl;

  cout << " absolute sum sq: " << diff.sum2() << endl;
  cout << " relative sum sq: " << relative.sum2() << endl;
  /*
  cout << " SitBar:" << _tetMesh->SitBar() << endl;


  MATRIX signs = full;
  for (int x = 0; x < signs.rows(); x++)
    for (int y = 0; y < signs.cols(); y++)
    {
      if (finiteDiffs(x,y) * full(x,y) < 0.0)
        signs(x,y) = 1;
      else
        signs(x,y) = 0;
    }
  cout << " signs match?: " << signs << endl;
  */
  exit(0);
}

//////////////////////////////////////////////////////////////////////
// build the list of vertices in collision -- cull away those that are a member
// of a colliding triangle but actually outside the surface
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::buildCollisionList(vector<pair<SURFACE*, int> >& collisions)
{
  // get the rigid transforms
  VEC3F& translation = _translation;
  MATRIX3 rotation = _rotation.toExplicitMatrix3x3();

  // build a list of unique collision vertices
  vector<TRIANGLE*> faces = _tetMesh->explicitSurfaceFaces();
  map<pair<VEC3F*, SURFACE*>, bool> vertexHash;
  map<SURFACE*, bool> surfaceHash;
  map<pair<VEC3F*, SURFACE*>, bool>::iterator iter;
  for (unsigned int x = 0; x < collisions.size(); x++)
  {
    int index = collisions[x].second;
    SURFACE* surface = collisions[x].first;
    TRIANGLE* face = faces[index];

    VEC3F v0 = rotation * (*(face->vertex(0))) + translation;
    VEC3F v1 = rotation * (*(face->vertex(1))) + translation;
    VEC3F v2 = rotation * (*(face->vertex(2))) + translation;

    surfaceHash[surface] = true;

    if (surface->inside(v0))
      vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(0), surface)] = true;
    if (surface->inside(v1))
      vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(1), surface)] = true;
    if (surface->inside(v2))
      vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(2), surface)] = true;
  }

  // clear the old collision list
  _collisionVertices.clear();

  // build a new one
  for (iter = vertexHash.begin(); iter != vertexHash.end(); iter++)
    _collisionVertices.push_back(iter->first);

  // clear the old collision surface list
  _collisionSurfaces.clear();

  // build a new one
  map<SURFACE*, bool>::iterator surfaceIter;
  for (surfaceIter = surfaceHash.begin(); surfaceIter != surfaceHash.end(); surfaceIter++)
    _collisionSurfaces.push_back(surfaceIter->first);
}

//////////////////////////////////////////////////////////////////////
// add in gravity both to the translation component and as a 
// body force to the deformable mesh. In the case of a skinned
// unconstrained mesh, there might still be some translation
// component in the basis. Otherwise, it projects to nothing and
// does no harm.
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::addGravity() 
{
  _translationForce += _gravityMagnitude * _gravityDown; 

  MATRIX3 R = _rotation.toExplicitMatrix3x3();
  MATRIX3 RT = R.transpose();
  VEC3F rotatedGravity = RT * _gravityDown;
  rotatedGravity *= _gravityMagnitude;

  //VECTOR projectedGravity = _bodyForceU * rotatedGravity.toVector();
  //_externalForces += projectedGravity;
}

void UNCONSTRAINED_SUBSPACE_INTEGRATOR::addBodyForce(VEC3F bodyForce)
{
  _translationForce += bodyForce;

  MATRIX3 R = _rotation.toExplicitMatrix3x3();
  MATRIX3 RT = R.transpose();
  VEC3F rotatedForce = RT * bodyForce;

  VECTOR projectedForce = _bodyForceU * bodyForce.toVector();
  _externalForces += projectedForce;
}

//////////////////////////////////////////////////////////////////////
// compute the kinetic energy of the entire mesh
//////////////////////////////////////////////////////////////////////
Real UNCONSTRAINED_SUBSPACE_INTEGRATOR::kineticEnergy()
{
  // update all the vertex positions
  _tetMesh->updateFullMesh();

  // rotate the linear velocity into the local frame
  QUATERNION rotation = _rotation;
  rotation.negateIm();
  VEC3F localTranslationVelocity = rotation.rotates(_translationVelocity);
  //VEC3F localTranslationVelocity = _rotation.toExplicitMatrix3x3().transpose() * _translationVelocity;

  // unproject the velocity
  VECTOR velocity = _tetMesh->U() * _velocity;

  Real totalEnergy = 0.0;

  for (unsigned int x = 0; x < velocity.size() / 3; x++)
  {
    VEC3F vertexVelocity;
    vertexVelocity[0] = velocity[3 * x];
    vertexVelocity[1] = velocity[3 * x + 1];
    vertexVelocity[2] = velocity[3 * x + 2];

    VEC3F angularComponent = cross(_angularVelocity, *_tetMesh->vertices(x));
    Real mass = _tetMesh->mass(x);

    //totalEnergy += mass * (vertexVelocity * vertexVelocity);
    //totalEnergy += 0.5 * mass * (angularComponent * angularComponent);
  }

  Real totalMass = _tetMesh->totalMass();
  totalEnergy += 0.5 * (_angularVelocity.toVector() * (_tetMesh->inertiaTensor() * _angularVelocity.toVector()));
  totalEnergy += 0.5 * totalMass * (localTranslationVelocity * localTranslationVelocity);

  return totalEnergy;
}

//////////////////////////////////////////////////////////////////////
// compute the rigid velocites based on the current rigids
//////////////////////////////////////////////////////////////////////
void UNCONSTRAINED_SUBSPACE_INTEGRATOR::computeRigidDerivatives()
{
  _translationVelocityOld = _translationVelocity;
  _translationAccelerationOld = _translationAcceleration;

  // update translation derivatives
  /*
  _translationVelocity = _alpha[3] * _translation;
  _translationVelocity -= _alpha[3] * _translationOld;
  _translationVelocity += _alpha[4] * _translationVelocityOld;
  _translationVelocity += _alpha[5] * _translationAccelerationOld;

  _translationAcceleration = _alpha[0] * _translation;
  _translationAcceleration -= _alpha[0] * _translationOld;
  _translationAcceleration -= _alpha[1] * _translationVelocityOld;
  _translationAcceleration -= _alpha[2] * _translationAccelerationOld;
  */
  //_translationVelocity = (_translation - _translationOld) * 0.5;
  //_translationAcceleration = (_translationVelocity - _translationVelocityOld) * 0.5;
  _translationVelocity = (_translation - _translationOld)  / _dt;
  _translationAcceleration = (_translationVelocity - _translationVelocityOld)  / _dt;

  //cout << " translation acceleration: " << _translationAcceleration << endl;
 
  _angularVelocityOld = _angularVelocity;
  _angularAccelerationOld = _angularAcceleration;

  MATRIX3 change = _rotation.toExplicitMatrix3x3() * _rotationOld.toExplicitMatrix3x3().transpose();

  QUATERNION qChange(change);
  //cout << " angular change: " << qChange << endl;
  VEC3F axis;
  Real angle = 0;
  qChange.axisAngle(axis, angle);

  if (fabs(qChange[0]) + fabs(qChange[1]) + fabs(qChange[2]) < 1e-8)
    angle = 0.0;
  angle *= 1.0 / _dt;

  _angularVelocity = axis * angle;

  _angularVelocity = _rotation.toExplicitMatrix3x3().transpose() * _angularVelocity;
 
  // velocity computation is still flaky, so dial it back a bit
  _angularVelocity *= 0.01;

  //cout << " angular velocity: " << _angularVelocity << endl;

  /*
  VECTOR eigenvalues;
  MATRIX eigenvectors;
  MATRIX logMatrix(preLog);
  logMatrix.eigensystem(eigenvalues, eigenvectors);
  //eigenvectors = eigenvectors.transpose();

  MATRIX invEigenvectors(eigenvectors);
  invEigenvectors.invert();

  cout << " original: " << logMatrix << endl;
  cout << " combined: " << (eigenvectors * MATRIX(eigenvalues) * invEigenvectors) << endl;
  //cout << " combined: " << (invEigenvectors * MATRIX(eigenvalues) * eigenvectors) << endl;
  cout << " eigenvectors: " << eigenvectors << endl;
  cout << " values: " << eigenvalues << endl;
  cout << endl;
  */
}
