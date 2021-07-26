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
// FULLSPACE_INTEGRATOR.h: interface for the FULLSPACE_INTEGRATOR class.
//
//////////////////////////////////////////////////////////////////////

#include "FULLSPACE_INTEGRATOR.h"
#include <TIMER.h>
#include <INVERTIBLE.h>
#ifdef USING_OPENMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor for newmark integrator
//////////////////////////////////////////////////////////////////////
FULLSPACE_INTEGRATOR::FULLSPACE_INTEGRATOR(TET_MESH* tetMesh, Real dt, Real alpha, Real beta, Real gravity) :
  _tetMesh(tetMesh),
  _maxNewtonSteps(3),
  _dt(dt),
  _rayleighAlpha(alpha),
  _rayleighBeta(beta),
  _time(0.0),
  _preconditioner(NULL),
  _totalTime(0),
  _totalSteps(0),
  _PCGDigits(5),
  _clickedNode(NULL),
  _forceMultiplier(0.1),
  _gravityMagnitude(gravity),
  _gravityDown(0,0,-1),
  _solverEps(1e-8),
  _substeppingEnabled(false),
  _totalNewtonStepsSeen(0)
{
  cout << "got rayleigh alpha = " << _rayleighAlpha << endl;
  cout << "got rayleigh beta = " << _rayleighBeta << endl;
  _rank = _tetMesh->TET_MESH::dofs();
  _position.resizeAndWipe(_rank);
  _velocity.resizeAndWipe(_rank);
  _acceleration.resizeAndWipe(_rank);
  _positionOld.resizeAndWipe(_rank);
  _velocityOld.resizeAndWipe(_rank);
  _accelerationOld.resizeAndWipe(_rank);
  _residual.resizeAndWipe(_rank);
  _temp.resizeAndWipe(_rank);
  _externalForces.resizeAndWipe(_rank);
  _externalForcesOld.resizeAndWipe(_rank);
  _positionSample.resizeAndWipe(_rank);
  _gravity.resizeAndWipe(_rank);

  _solverEps = 1e-2;
  
  int vectorSize = sizeof(Real) * _rank * 10;
  cout << " Integrator size: " << vectorSize / pow(2.0, 20.0) << " MB " << endl;

  // implicit Newmark
	//_beta = 0.25;
	//_gamma = 0.50;

  // fully implicit Euler
  _beta = 0.5;
	_gamma = 1.0;

	_alpha[0] = 1.0 / (_beta * _dt * _dt);
	_alpha[1] = 1.0 / (_beta * _dt);
	_alpha[2] = (1.0 - 2.0 * _beta) / (2.0 * _beta);
	_alpha[3] = _gamma / (_beta * _dt);
	_alpha[4] = 1.0 - _gamma / _beta;
	_alpha[5] = (1.0 - _gamma / (2.0 * _beta)) * _dt;

  // alphas for acceleration-level update
  _accelerationAlpha[0] = _dt * (1.0 - _gamma);
  _accelerationAlpha[1] = _dt * _gamma;
  _accelerationAlpha[2] = _dt;
  _accelerationAlpha[3] = _dt * _dt * (0.5 - _beta);
  _accelerationAlpha[4] = _dt * _dt * _beta;

  /*
  // For fully implicit Euler, these have to be set explicitly. 
  // It is not possible to get fully implicit Euler using the Newmark consts.
  _accelerationAlpha[3] = 0.0;
  _accelerationAlpha[4] = _dt * _dt;
  */

  // precompute the gravity vector
  for (int x = 0; x < _tetMesh->unconstrainedNodes(); x++)
  {
    VEC3F downCopy = _gravityDown;
    downCopy.normalize();
    downCopy *= _gravityMagnitude;
    downCopy[0] *= _tetMesh->massMatrix()(3 * x, 3 * x);
    downCopy[1] *= _tetMesh->massMatrix()(3 * x + 1, 3 * x + 1);
    downCopy[2] *= _tetMesh->massMatrix()(3 * x + 2, 3 * x + 2);

    _gravity[3 * x] = downCopy[0];
    _gravity[3 * x + 1] = downCopy[1];
    _gravity[3 * x + 2] = downCopy[2];
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::computeAlphas()
{
	_beta = 0.25;
	_gamma = 0.50;
	_alpha[0] = 1.0 / (_beta * _dt * _dt);
	_alpha[1] = 1.0 / (_beta * _dt);
	_alpha[2] = (1.0 - 2.0 * _beta) / (2.0 * _beta);
	_alpha[3] = _gamma / (_beta * _dt);
	_alpha[4] = 1.0 - _gamma / _beta;
	_alpha[5] = (1.0 - _gamma / (2.0 * _beta)) * _dt;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::stepFullImplicit()
{
  TIMER total;

  // allocate matrices if necessary -- should only fire
  // on the first step
  if (_CDamping.rows() != _rank)
  {
    _CDamping.resizeAndWipe(_rank, _rank);
    cacheStaticDampingFull();
  }
  if (_A.rows() != _rank)
    _A.resizeAndWipe(_rank, _rank);

  // get the reduced mass matrix
  SPARSE_MATRIX& M = _tetMesh->massMatrix();

  // accumulate external forces
  computeExternalForces();
 
  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  int step = 0;
  while (step < 100 && maxR > eps)
  //for (int step = 0; step < 3; step++)
  {
    // update tet mesh with the new q vector
    TIMER updateMesh;
    VECTOR& x = _tetMesh->x();

    // an assignment is overly expensive -- the "=" forces an allocation
    //q = _position;
    x.copyInplace(_position);
    _positionSample.equals(_position);
    _tetMesh->TET_MESH::updateFullMesh();
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // get the reduced internal forces
    TIMER internalTimer;
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the reduced stiffness matrix
    TIMER stiffnessTimer;
    MATRIX& K = _tetMesh->TET_MESH::generateStiffnessMatrix();
    _timingBreakdown["Stiffness Assembly"] += stiffnessTimer.timing();

    // compute the LHS of the residual:
    // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
    TIMER RHSTimer;
    _temp.clearingAxpy(-_alpha[2], _accelerationOld);
    _temp.axpy(-_alpha[1], _velocityOld);
    _temp.axpy(-_alpha[0], _positionOld);
    _temp.axpy(_alpha[0], _position);

    // do an explicit multiply because of the sparse matrix
    for (int i = 0; i < M.rows(); i++)
      _temp(i) *= M(i,i);

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    _residual.clearingAxpy(_alpha[5], _accelerationOld);
    _residual.axpy(_alpha[4], _velocityOld);
    _residual.axpy(-_alpha[3], _positionOld);
    _residual.axpy(_alpha[3], _position);
    _residual = _CDamping * _residual;

    // assemble full residual: LHS + RHS + R - F
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();
   
    // assemble system matrix A
    TIMER ATimer;
    _A = K;
    for (int i = 0; i < M.rows(); i++)
      _A(i,i) += _alpha[0] * M(i,i);
    _A.axpy(_alpha[3], _CDamping);
    _timingBreakdown["A Assembly"] += ATimer.timing();

    // solve with LU factorization
    TIMER solverTimer;
    _A.solve(_residual);
    //_A.solveICCG(_solution, _residual);
    _timingBreakdown["Solver"] += solverTimer.timing();

    // update positions
    _position -= _residual;
    //_position -= _solution;
    step++;
  }

	// update velocity
  TIMER finalUpdate;
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);

	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], _accelerationOld);
  _acceleration.axpy(-_alpha[1], _velocityOld);
  _acceleration.axpy(-_alpha[0], _positionOld);
  _acceleration.axpy(_alpha[0], _position);

  // update node positions
  _tetMesh->TET_MESH::updateFullMesh();
  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment counters
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// trying out a support routine for online cubature generation
//////////////////////////////////////////////////////////////////////
VECTOR& FULLSPACE_INTEGRATOR::computeInternalForces(bool invertible)
{
  /*
  cout << __FILE__ << " " << __LINE__ << " : " << endl;
  cout << "DEBUG: INVERSION DEACTIVATED" << endl;
  cout << "       REMEMBER IT'S DEACTIVATED IN INVERTIBLE AS WELL." << endl;
  */
  if (invertible)
  {
    cacheDiagonalizations();
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces(_Us, _Fhats, _Vs);
    return R;
  }
  VECTOR& R = _tetMesh->TET_MESH::generateInternalForces();
  return R;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::stepSparseImplicit(bool computeExternals)
{
  TIMER total;

  TIMER preambleTimer;
  // allocate matrices if necessary -- should only fire
  // on the first step
  if (_CDampingSparse.rows() != _rank)
  {
    cout << "Caching damping stuff" << endl;
    _CDampingSparse.resize(_rank, _rank);
    cacheStaticDampingSparse();
    cout << "Done" << endl;
    
    cout << "Setting sparsity in _ASparse" << endl;
    _ASparse.setSparsity(_CDampingSparse);
    cout << "Done" << endl;
  }
  
  //if (_ASparse.rows() != _rank)
  //  _ASparse.resize(_rank, _rank);
  
  if (_solution.size() != _rank)
    _solution.resizeAndWipe(_rank);

  // get the mass matrix
  SPARSE_MATRIX& M = _tetMesh->TET_MESH::massMatrix();

  // accumulate external forces
  if (computeExternals)
    computeExternalForces();
 
  // add in frictional forces
  computeFrictionForces();
  _externalForces += _frictionForces;

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;
  _timingBreakdown["Preamble"] += preambleTimer.timing();

  cout << "==================================================" << endl;
  cout << " Solving implicit timestep " << _totalSteps << endl;
  cout << "==================================================" << endl;

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  int step = 0;
  //while (step < 10 && maxR > eps)
  //while (step < 100 && maxR > eps)
#ifndef _WIN32
  while (step < _maxNewtonSteps && maxR > eps)
#else
  while (step < _maxNewtonSteps && maxR > eps)
#endif
  {
    // update tet mesh with the new q vector
    TIMER updateMesh;
    VECTOR& x = _tetMesh->TET_MESH::x();

    // an assignment is overly expensive -- the "=" forces an allocation
    x.equals(_position);
    _positionSample.equals(_position);

    _tetMesh->TET_MESH::updateFullMesh();
    _timingBreakdown["Update Mesh"] += updateMesh.timing();

    // get the internal forces
    TIMER internalTimer;
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();
    
    // compute the LHS of the residual:
    // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
    TIMER RHSTimer;
    _temp.clearingAxpy(-_alpha[2], _accelerationOld);
    _temp.axpy(-_alpha[1], _velocityOld);
    _temp.axpy(-_alpha[0], _positionOld);
    _temp.axpy(_alpha[0], _position);

    // do an explicit multiply because of the sparse matrix
    for (int i = 0; i < M.rows(); i++)
      _temp(i) *= M(i,i);

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    _residual.clearingAxpy(_alpha[5], _accelerationOld);
    _residual.axpy(_alpha[4], _velocityOld);
    _residual.axpy(-_alpha[3], _positionOld);
    _residual.axpy(_alpha[3], _position);
    _residual = _CDampingSparse * _residual;

    // assemble full residual: LHS + RHS + R - F
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;

    // add collision responses
    updateVelocities();
    computeCollisionForces();
    computeCollisionForceJacobians();
    _residual -= _collisionForces;

    // apply constraints
    applyResidualConstraints();

    maxR = _residual.norm2();
    step++;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

#if 1
    cout << " Newton-Raphson residual " << step << ": " << maxR;
#ifndef _WIN32    
    cout << endl; flush(cout);
#endif
#endif
    if (maxR < eps) break;

#ifndef _WIN32    
    // assemble system matrix A
    TIMER ACopy;
    _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_ASparse);
    for (int i = 0; i < M.rows(); i++)
      //_ASparse(i,i) += _alpha[0] * M(i,i);
      _ASparse.add(i,i,_alpha[0] * M(i,i));
    _ASparse.axpy(_alpha[3], _CDampingSparse);
    _timingBreakdown["Stiffness Assembly"] += ACopy.timing();

    TIMER addForceTimer;
    // add force jacocbians
    _ASparse.axpy(1.0, _collisionForceJacobians);

    // apply constraints
    applyMatrixConstraints();
    _timingBreakdown["Adding collision forces to matrix"] += addForceTimer.timing();

    // solve with sparse solver
    TIMER solverTimer;
    Real PCGEps = pow(10.0, -_PCGDigits);
    _ASparse.eps() = PCGEps;
    bool converged = _ASparse.solveCG(_solution, _residual);

    // DEBUG
    if (!converged)
    {
      cout << "Nan encountered!" << endl;
      cout << " Is mesh inverted? " << _tetMesh->inverted(true) << endl;
    }

    //_ASparse.solvePCG(_solution, _residual, _preconditioner);
    _timingBreakdown["Solver"] += solverTimer.timing();
#else
    SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix();

    // assemble system matrix A
    TIMER ACopy;
    if (_ASparse.rows() != _rank)
      _ASparse.resize(_rank, _rank);
    _ASparse.equals(K);
    _timingBreakdown["Stiffness Assembly"] += ACopy.timing();
    TIMER AMultiply;
    for (int i = 0; i < M.rows(); i++)
      _ASparse(i,i) += _alpha[0] * M(i,i);
    _timingBreakdown["A Multiply"] += AMultiply.timing();
    TIMER ATimer;
    _ASparse.axpy(_alpha[3], _CDampingSparse);
    _timingBreakdown["A Assembly"] += ATimer.timing();

    if (_preconditioner == NULL)
      _preconditioner = new INCOMPLETE_CHOLESKY(_ASparse);

    // solve with sparse solver
    TIMER solverTimer;
    Real PCGEps = pow(10.0, -_PCGDigits);
    //_ASparse.eps() = _solverEps / 100000.0;
    //_ASparse.eps() = _solverEps / PCGEps;
    _ASparse.eps() = PCGEps;
    //_ASparse.solveCG(_solution, _residual);
    _ASparse.solvePCG(_solution, _residual, _preconditioner);

    //_ASparse.solvePCG(_solution, _residual, _preconditioner);
    _timingBreakdown["Solver"] += solverTimer.timing();
#endif

    // update positions
    _position -= _solution;
  }
#ifdef _WIN32
  cout << endl << " Step complete. " << endl;
#endif

	// update velocity
  TIMER finalUpdate;
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);

	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], _accelerationOld);
  _acceleration.axpy(-_alpha[1], _velocityOld);
  _acceleration.axpy(-_alpha[0], _positionOld);
  _acceleration.axpy(_alpha[0], _position);

  // update node positions
  VECTOR& x = _tetMesh->TET_MESH::x();
  x.copyInplace(_position);
  _tetMesh->TET_MESH::updateFullMesh();

  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment counters
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration
//////////////////////////////////////////////////////////////////////
bool FULLSPACE_INTEGRATOR::stepSparseImplicitInvertible(bool computeExternals, bool substepping)
{
  TIMER total;

  TIMER preambleTimer;
  // allocate matrices if necessary -- should only fire
  // on the first step
  if (_CDampingSparse.rows() != _rank)
  {
    _CDampingSparse.resize(_rank, _rank);
    cacheStaticDampingSparse();

    _ASparse.setSparsity(_CDampingSparse);
  }

  //if (_ASparse.rows() != _rank)
  //  _ASparse.resize(_rank, _rank);

  if (_solution.size() != _rank)
    _solution.resizeAndWipe(_rank);
 
  // get the mass matrix
  SPARSE_MATRIX& M = _tetMesh->TET_MESH::massMatrix();

  // accumulate external forces
  if (computeExternals)
    computeExternalForces();

  // add in frictional forces
  computeFrictionForces();
  _externalForces += _frictionForces;

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;
  _timingBreakdown["Preamble"] += preambleTimer.timing();

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  int step = 0;
  //while (step < 100 && maxR > eps)
  while (step < _maxNewtonSteps && maxR > eps)
  {
    // update tet mesh with the new q vector
    TIMER updateMesh;
    VECTOR& x = _tetMesh->TET_MESH::x();

    // an assignment is overly expensive -- the "=" forces an allocation
    x.equals(_position);
    _positionSample.equals(_position);

    _tetMesh->TET_MESH::updateFullMesh();
    _timingBreakdown["Update Mesh"] += updateMesh.timing();

    // diagonalize
    TIMER diagonal;
    cacheDiagonalizations();
    _timingBreakdown["Diagonalization and Stiffness Density"] += diagonal.timing();

    // get the internal forces
    TIMER internalTimer;
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces(_Us, _Fhats, _Vs);
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // compute the LHS of the residual:
    // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
    TIMER RHSTimer;
    _temp.clearingAxpy(-_alpha[2], _accelerationOld);
    _temp.axpy(-_alpha[1], _velocityOld);
    _temp.axpy(-_alpha[0], _positionOld);
    _temp.axpy(_alpha[0], _position);

    // do an explicit multiply because of the sparse matrix
    for (int i = 0; i < M.rows(); i++)
      _temp(i) *= M(i,i);

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    _residual.clearingAxpy(_alpha[5], _accelerationOld);
    _residual.axpy(_alpha[4], _velocityOld);
    _residual.axpy(-_alpha[3], _positionOld);
    _residual.axpy(_alpha[3], _position);
    _residual = _CDampingSparse * _residual;

    // assemble full residual: LHS + RHS + R - F
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;

    // add collision responses
    updateVelocities();
    computeCollisionForces();
    computeCollisionForceJacobians();
    _residual -= _collisionForces;

    // apply constraints
    applyResidualConstraints();
    
    maxR = _residual.norm2();
    step++;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

    cout << " Newton-Raphson residual " << step << ": " << maxR;
    cout << endl; flush(cout);
    if (maxR < eps) break;
   
    // assemble system matrix A
    TIMER genSparseTimer;
    //_tetMesh->TET_MESH::generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses, _ASparse);
    SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);
    _timingBreakdown["Generate Sparse Stiffness Matrix"] += genSparseTimer.timing();
    TIMER copyPetscTimer;
    _ASparse.equals(K);
    _timingBreakdown["Petsc Stiffness copy"] += copyPetscTimer.timing();

    TIMER addMassTimer;
    for (int i = 0; i < M.rows(); i++)
      _ASparse.add(i,i, _alpha[0] * M(i,i));
    _timingBreakdown["Adding mass to stiffness matrix"] += addMassTimer.timing();
    TIMER addDampingTimer;
    _ASparse.axpy(_alpha[3], _CDampingSparse);
    _timingBreakdown["Adding damping to stiffness matrix"] += addDampingTimer.timing();

    TIMER addForceTimer;
    // add force jacocbians
    _ASparse.axpy(1.0, _collisionForceJacobians);

    // apply constraints
    applyMatrixConstraints();
    _timingBreakdown["Adding collision forces to matrix"] += addForceTimer.timing();

    // solve with CG
    TIMER solverTimer;
    Real PCGEps = pow(10.0, -_PCGDigits);
    _ASparse.eps() = PCGEps;

    bool converged = _ASparse.solveCG(_solution, _residual);
    _timingBreakdown["Solver"] += solverTimer.timing();

    if (!converged)
    {
      cout << " Is mesh inverted? " << _tetMesh->inverted(true) << endl;
      //if (_substeppingEnabled && !substepping) break;
      return false;
    }

    // update positions
    _position -= _solution;
  }

  if (maxR > eps)
    return false;

  if (maxR > eps && !substepping && _substeppingEnabled)
  {
    // roll back the state
    _position = _positionOld;
    _velocityOld = _velocityOld;
    _acceleration = _accelerationOld;
    _totalTime += total.timing();
    return false;
  }

	// update velocity
  TIMER finalUpdate;
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);

	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], _accelerationOld);
  _acceleration.axpy(-_alpha[1], _velocityOld);
  _acceleration.axpy(-_alpha[0], _positionOld);
  _acceleration.axpy(_alpha[0], _position);

  // update node positions
  VECTOR& x = _tetMesh->TET_MESH::x();
  x.copyInplace(_position);
  _tetMesh->TET_MESH::updateFullMesh();

  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment counters
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();

  return true;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration
//////////////////////////////////////////////////////////////////////
bool FULLSPACE_INTEGRATOR::stepSparseImplicitWithCollisions(vector<pair<SURFACE*, int> > collisionPairs)
{
  TIMER total;

  // build collision lists
  vector<pair<VEC3F*, SURFACE*> > collisionVertices;
  vector<SURFACE*> collisionSurfaces;
  buildCollisionList(collisionPairs, collisionVertices, collisionSurfaces);

  TIMER preambleTimer;
  // allocate matrices if necessary -- should only fire
  // on the first step
  if (_CDampingSparse.rows() != _rank)
  {
    _CDampingSparse.resize(_rank, _rank);
    cacheStaticDampingSparse();

    _ASparse.setSparsity(_CDampingSparse);
  }

  //if (_ASparse.rows() != _rank)
  //  _ASparse.resize(_rank, _rank);

  if (_solution.size() != _rank)
    _solution.resizeAndWipe(_rank);
 
  // get the mass matrix
  SPARSE_MATRIX& M = _tetMesh->TET_MESH::massMatrix();

  // accumulate external forces
  computeExternalForces();

  cout << "==================================================" << endl;
  cout << " Solving implicit with collisions timestep " << _totalSteps << endl;
  cout << "==================================================" << endl;

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;
  _timingBreakdown["Preamble"] += preambleTimer.timing();

  // cache the smallest seen position
  VECTOR smallest(_position);
  Real smallestResidualSeen = 0.0;

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  int step = 0;
  while (step < _maxNewtonSteps && maxR > eps)
  {
    // update tet mesh with the new q vector
    TIMER updateMesh;
    VECTOR& x = _tetMesh->TET_MESH::x();

    // an assignment is overly expensive -- the "=" forces an allocation
    x.equals(_position);
    _positionSample.equals(_position);

    _tetMesh->TET_MESH::updateFullMesh();
    _timingBreakdown["Update Mesh"] += updateMesh.timing();

    // diagonalize
    TIMER diagonal;
    cacheDiagonalizations();
    _timingBreakdown["Diagonalization and Stiffness Density"] += diagonal.timing();

    // get the internal forces
    TIMER internalTimer;
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces(_Us, _Fhats, _Vs);
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // compute the LHS of the residual:
    // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
    TIMER RHSTimer;
    _temp.clearingAxpy(-_alpha[2], _accelerationOld);
    _temp.axpy(-_alpha[1], _velocityOld);
    _temp.axpy(-_alpha[0], _positionOld);
    _temp.axpy(_alpha[0], _position);

    // do an explicit multiply because of the sparse matrix
    for (int i = 0; i < M.rows(); i++)
      _temp(i) *= M(i,i);

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    _residual.clearingAxpy(_alpha[5], _accelerationOld);
    _residual.axpy(_alpha[4], _velocityOld);
    _residual.axpy(-_alpha[3], _positionOld);
    _residual.axpy(_alpha[3], _position);
    _residual = _CDampingSparse * _residual;

    // assemble full residual: LHS + RHS + R - F
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;

    // add collision responses
    updateVelocities();
   
    // NEW: add in collision residuals
    _residual += computeCollisionResidual(collisionVertices);

    maxR = _residual.norm2();
    step++;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

    cout << " Newton-Raphson residual " << step << ": " << maxR;
    cout << endl; flush(cout);
    if (maxR < eps) break;
   
#if !USING_PETSC
    cout << " Need Petsc installed to run " << __FUNCTION__ << "!!!" << endl;
    exit(0);
#endif

    // assemble system matrix A
    TIMER genSparseTimer;
    SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);
    _timingBreakdown["Generate Sparse Stiffness Matrix"] += genSparseTimer.timing();
    TIMER copyPetscTimer;
    _ASparse.equals(K);
    _timingBreakdown["Petsc Stiffness copy"] += copyPetscTimer.timing();

    TIMER addMassTimer;
    for (int i = 0; i < M.rows(); i++)
      _ASparse.add(i,i, _alpha[0] * M(i,i));
    _timingBreakdown["Adding mass to stiffness matrix"] += addMassTimer.timing();
    TIMER addDampingTimer;
    _ASparse.axpy(_alpha[3], _CDampingSparse);
    _timingBreakdown["Adding damping to stiffness matrix"] += addDampingTimer.timing();

    // NEW: add in collision Jacobians
    addCollisionForceJacobians(collisionVertices);

    // solve with CG
    TIMER solverTimer;
    Real PCGEps = pow(10.0, -_PCGDigits);
    _ASparse.eps() = PCGEps;
    //_ASparse.usePCG();
    _ASparse.useMINRES();
    bool converged = _ASparse.solveCG(_solution, _residual);
    _timingBreakdown["Solver"] += solverTimer.timing();

    if (!converged)
    {
      cout << " Is mesh inverted? " << _tetMesh->inverted(true) << endl;
      return false;
    }

    // update positions
    _position -= _solution;

    // update the smallest found
    if (step == 1 || maxR < smallestResidualSeen)
    {
      smallest = _position;
      smallestResidualSeen = maxR; 
    }
  }
  _totalNewtonStepsSeen += step;

  // if the final solution is bad, use the best one seen
  if (maxR > smallestResidualSeen)
  {
    cout << " Rolled solution back to iteration with residual " << smallestResidualSeen << endl;
    _position = smallest;
  }

	// update velocity
  TIMER finalUpdate;
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);

	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], _accelerationOld);
  _acceleration.axpy(-_alpha[1], _velocityOld);
  _acceleration.axpy(-_alpha[0], _positionOld);
  _acceleration.axpy(_alpha[0], _position);

  // update node positions
  VECTOR& x = _tetMesh->TET_MESH::x();
  x.copyInplace(_position);
  _tetMesh->TET_MESH::updateFullMesh();

  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment counters
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();

  return true;
}

//////////////////////////////////////////////////////////////////////
// Quasistatic implicit
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::stepSparseQuasistatic()
{
  TIMER total;

  TIMER preambleTimer;
  if (_solution.size() != _rank)
    _solution.resizeAndWipe(_rank);
 
  // accumulate external forces
  computeExternalForces();
  _timingBreakdown["Preamble"] += preambleTimer.timing();
 
  // copy the current values into the old values
  _positionOld = _position;

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  int step = 0;
  while (step < _maxNewtonSteps && maxR > eps)
  {
    TIMER updateMesh;
    // update tet mesh with the new q vector
    VECTOR& x = _tetMesh->x();

    // an assignment is overly expensive -- the "=" forces an allocation
    //q = _position;
    x.equals(_position);
    _positionSample.equals(_position);

    _tetMesh->TET_MESH::updateFullMesh();
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // get the reduced internal forces
    TIMER internalTimer;
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the sparse stiffness matrix
    TIMER stiffnessTimer;
    SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix();
    _timingBreakdown["Stiffness Assembly"] += stiffnessTimer.timing();

    // if this is the first time, set the sparsity of the PETSC matrix
    if (_ASparse.rows() != _rank)
      _ASparse.setSparsity(K);

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    TIMER RHSTimer;
    _residual = _externalForces;
    _residual += R;

    // add collition responses
    _velocity.clear();
    computeCollisionForces();
    computeCollisionForceJacobians();
    _residual -= _collisionForces;
   
    // apply constraints
    applyResidualConstraints();

    maxR = _residual.norm2();
    step++;
    cout << " Newton-Raphson residual: " << _residual.norm2() << endl;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();
    if (maxR < eps) break;
  
    // solve with explicit stiffness matrix A
    TIMER ACopy;
    _ASparse.equals(K);
    _timingBreakdown["A Copy"] += ACopy.timing();

    // add force jacobians and constraints
    TIMER addForceTimer;
    _ASparse.axpy(1.0, _collisionForceJacobians);
    applyMatrixConstraints();
    _timingBreakdown["Adding collision forces to matrix"] += addForceTimer.timing();
      
    TIMER solverTimer;
    Real PCGEps = pow(10.0, -_PCGDigits);
    _ASparse.eps() = PCGEps;
    //bool converged = _ASparse.solveCG(_solution, _residual);
    _ASparse.solveCG(_solution, _residual);
    _timingBreakdown["Solver"] += solverTimer.timing();

    // update positions
    _position += _solution;
  }

  TIMER finalUpdate;
  // update node positions
  VECTOR& x = _tetMesh->x();
  x.copyInplace(_position);
  _tetMesh->TET_MESH::updateFullMesh();
  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment counters
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// Quasistatic implicit
//////////////////////////////////////////////////////////////////////
bool FULLSPACE_INTEGRATOR::stepSparseQuasistaticInvertible()
{
  TIMER total;

  TIMER preambleTimer;
  if (_solution.size() != _rank)
    _solution.resizeAndWipe(_rank);
 
  // accumulate external forces
  computeExternalForces();
  _timingBreakdown["Preamble"] += preambleTimer.timing();
 
  // copy the current values into the old values
  _positionOld = _position;

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  int step = 0;
  cout << "==================================================" << endl;
  cout << " Solving quasistatic invertible timestep " << _totalSteps << endl;
  cout << "==================================================" << endl;
  while (step < _maxNewtonSteps && maxR > eps)
  {
    TIMER updateMesh;
    // update tet mesh with the new q vector
    VECTOR& x = _tetMesh->x();

    // an assignment is overly expensive -- the "=" forces an allocation
    //q = _position;
    x.equals(_position);
    _positionSample.equals(_position);

    _tetMesh->TET_MESH::updateFullMesh();
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // caching diagonalizations
    INVERTIBLE::inversions() = 0;
    TIMER diagonalTimer;
    cacheDiagonalizations();
    _timingBreakdown["Diagonalization"] += diagonalTimer.timing();  
    Real percent = 100.0 * INVERTIBLE::inversions() / (_tetMesh->tets().size() * 9);
    static Real maxSeen = percent;
    static Real minSeen = percent;
    if (percent > maxSeen) maxSeen = percent;
    if (percent < minSeen) minSeen = percent;

    /*
    cout << " tet inversions: " << INVERTIBLE::inversions() << " of " << _tetMesh->dofs() * 9 << "(" << percent << "%)" << endl;
    cout << " max seen: " << maxSeen << endl;
    cout << " min seen: " << minSeen << endl;
    */
      
    // get the reduced internal forces
    TIMER internalTimer;
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces(_Us, _Fhats, _Vs);
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the sparse stiffness matrix
    TIMER stiffnessTimer;
    SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);
    _timingBreakdown["Stiffness Assembly"] += stiffnessTimer.timing();

    // if this is the first time, set the sparsity of the PETSC matrix
    if (_ASparse.rows() != _rank)
      _ASparse.setSparsity(K);

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    TIMER RHSTimer;
    _residual = _externalForces;
    _residual += R;

    // add collision responses
    _velocity.clear();
    computeCollisionForces();
    computeCollisionForceJacobians();
    _residual -= _collisionForces;
   
    // apply constraints
    applyResidualConstraints();

    maxR = _residual.norm2();
    step++;
    cout << " Newton-Raphson residual " << step << ": " << maxR << endl;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

    /*
    // DEBUG -- breakdown detector
    if (step == 0) oldResidual = maxR;
    Real relativeProgress = fabs((oldResidual - maxR) / oldResidual);
    oldResidual = maxR;
    cout << " Relative progress: " << relativeProgress << endl;
    if (relativeProgress < 1e-2)
      badProgress++;
    else
      badProgress = 0;

    if (badProgress == 10)
    {
      cout << " Solver breakdown!!!!!! " << endl;
      break;
    }
    */

    if (maxR < eps) break;

    // solve with explicit stiffness matrix A
    TIMER ACopy;
    _ASparse.equals(K);
    _timingBreakdown["A Copy"] += ACopy.timing();

    // add force jacobians and constraints
    TIMER addForceTimer;
    _ASparse.axpy(1.0, _collisionForceJacobians);
    applyMatrixConstraints();
    _timingBreakdown["Adding collision forces to matrix"] += addForceTimer.timing();

    TIMER solverTimer;
    Real PCGEps = pow(10.0, -_PCGDigits);
    _ASparse.eps() = PCGEps;
    //bool converged = _ASparse.solveCG(_solution, _residual);
    _ASparse.solveCG(_solution, _residual);
    _timingBreakdown["Solver"] += solverTimer.timing();

    // update positions
    _position += _solution;
  }

  // DEBUG
  static int maxStepsSeen = 0;
  if (step > maxStepsSeen)
    maxStepsSeen = step;
  cout << "Max Newton steps seen: " << maxStepsSeen << endl;
  static int totalStepsSeen = 0;
  totalStepsSeen += step;
  cout << "Mean Newton steps seen: " << (float)totalStepsSeen / (_totalSteps + 1) << endl;
  static Real totalResidual = 0;
  totalResidual += maxR;
  cout << "Mean residual seen: " << totalResidual / (_totalSteps + 1) << endl;

  TIMER finalUpdate;
  // update node positions
  VECTOR& x = _tetMesh->x();
  x.copyInplace(_position);
  _tetMesh->TET_MESH::updateFullMesh();
  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment counters
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();

  return (step == _maxNewtonSteps);
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
void FULLSPACE_INTEGRATOR::stepFullExplicit()
{
  TIMER total;

  TIMER preambleTimer;
  if (_A.rows() != _rank)
    _A.resizeAndWipe(_rank, _rank);

  // allocate matrices if necessary -- should only fire
  // on the first step
  if (_CDamping.rows() != _rank)
  {
    _CDamping.resizeAndWipe(_rank, _rank);
    cacheStaticDampingFull();

    // prefactor M and C
    SPARSE_MATRIX& M = _tetMesh->massMatrix();
    _A.clear();
    for (int i = 0; i < M.rows(); i++)
      _A(i,i) = M(i,i);
    _A.axpy(_dt / 2.0, _CDamping);

    cout << " Factoring system matrix (just once) ... ";
    _A.factorCholesky();
    cout << " done." << endl;
  }

  // explicit does not use the acceleration vectors,
  // so use them as scratch space
  VECTOR& scratch = _acceleration;
  _timingBreakdown["Preamble"] += preambleTimer.timing();
  
  // update tet mesh with the new q vector
  TIMER updateTimer;
  VECTOR& x = _tetMesh->x();
  x.copyInplace(_position);
  _positionSample.equals(_position);
  _tetMesh->TET_MESH::updateFullMesh();
  _timingBreakdown["Update Mesh"] += updateTimer.timing();

  // get the forces
  // R = (P_n - R(u_n))
  TIMER internalTimer;
  computeExternalForces();
  VECTOR& R = _externalForces;
  R += _tetMesh->TET_MESH::generateInternalForces();
  _timingBreakdown["Internal Forces"] += internalTimer.timing();

  // start accum'ing in _temp,
  // switch to 'scratch' eventually to avoid a dynamic allocation
  // scratch = M * (2 * u_n - u_{n-1})
  TIMER RHSTimer;
  SPARSE_MATRIX& M = _tetMesh->massMatrix();
  _temp.clearingAxpy(2.0, _position);
  _temp.axpy(-1.0, _positionOld);
  // do explicit multiply because of sparse matrix
  for (int i = 0; i < scratch.size(); i++)
    scratch(i) = M(i,i) * _temp(i);

  // scratch += dt^2 * (P_n - R(u_n))
  scratch.axpy(_dt * _dt, R);

  // add damping
  // scratch += dt/2 * C * u_{n-1}
  VECTOR damping = _CDamping * _positionOld;
  damping *= _dt / 2.0;
  scratch.axpy(1.0, damping);
  _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

  // A = M + dt/2 * C
  /*
  TIMER ATimer;
  _A.clear();
  for (int i = 0; i < M.rows(); i++)
    _A(i,i) = M(i,i);
  _A.axpy(_dt / 2.0, _CDamping);
  _timingBreakdown["A Assembly"] += ATimer.timing();
  TIMER solve;
  _A.solve(scratch);
  _timingBreakdown["Solver"] += solve.timing();
  */
  _A.solveCholesky(scratch);

  TIMER finalUpdate;
  _positionOld.swap(_position);
  _position.swap(scratch);

  _velocity = _position - _positionOld;
  _velocity *= 1.0 / _dt;

  // update node positions
  _tetMesh->TET_MESH::updateFullMesh();
  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment simulation time
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();
}
/*
{
  TIMER total;

  TIMER preambleTimer;
  // allocate matrices if necessary -- should only fire
  // on the first step
  if (_CDamping.rows() != _rank)
  {
    _CDamping.resizeAndWipe(_rank, _rank);
    cacheStaticDampingFull();
  }
  if (_A.rows() != _rank)
    _A.resizeAndWipe(_rank, _rank);

  // explicit does not use the acceleration vectors,
  // so use them as scratch space
  VECTOR& scratch = _acceleration;
  _timingBreakdown["Preamble"] += preambleTimer.timing();
  
  // update tet mesh with the new q vector
  TIMER updateTimer;
  VECTOR& x = _tetMesh->x();
  x.copyInplace(_position);
  _positionSample.equals(_position);
  _tetMesh->TET_MESH::updateFullMesh();
  _timingBreakdown["Update Mesh"] += updateTimer.timing();

  // get the forces
  // R = (P_n - R(u_n))
  TIMER internalTimer;
  computeExternalForces();
  VECTOR& R = _externalForces;
  R += _tetMesh->TET_MESH::generateInternalForces();
  _timingBreakdown["Internal Forces"] += internalTimer.timing();

  // start accum'ing in _temp,
  // switch to 'scratch' eventually to avoid a dynamic allocation
  // scratch = M * (2 * u_n - u_{n-1})
  TIMER RHSTimer;
  SPARSE_MATRIX& M = _tetMesh->massMatrix();
  _temp.clearingAxpy(2.0, _position);
  _temp.axpy(-1.0, _positionOld);
  // do explicit multiply because of sparse matrix
  for (int i = 0; i < scratch.size(); i++)
    scratch(i) = M(i,i) * _temp(i);

  // scratch += dt^2 * (P_n - R(u_n))
  scratch.axpy(_dt * _dt, R);

  // add damping
  // scratch += dt/2 * C * u_{n-1}
  VECTOR damping = _CDamping * _positionOld;
  damping *= _dt / 2.0;
  scratch.axpy(1.0, damping);
  _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

  // A = M + dt/2 * C
  TIMER ATimer;
  _A.clear();
  for (int i = 0; i < M.rows(); i++)
    _A(i,i) = M(i,i);
  _A.axpy(_dt / 2.0, _CDamping);
  _timingBreakdown["A Assembly"] += ATimer.timing();
  TIMER solve;
  _A.solve(scratch);
  _timingBreakdown["Solver"] += solve.timing();

  TIMER finalUpdate;
  _positionOld.swap(_position);
  _position.swap(scratch);

  _velocity = _position - _positionOld;
  _velocity *= 1.0 / _dt;

  // update node positions
  _tetMesh->TET_MESH::updateFullMesh();
  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment simulation time
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();
}
*/

//////////////////////////////////////////////////////////////////////
// Explicit Integration solved with Pardiso sparse solver
//
// This follows pages 324-325 of: 
// Computational Contact Mechanics by Peter Wriggers, 2006
//
// The equation to solve is:
// (M + (dt/2) * C) u_{n+1} = 
// (dt)^2 * [P_n - R(u_n)] + (dt/2) * C * u_{n-1} + M* (2 * u_n - u_{n-1})
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::stepSparseExplicit()
{
  TIMER total;

  TIMER preambleTimer;
  // allocate matrices if necessary -- should only fire
  // on the first step
  if (_CDampingSparse.rows() != _rank)
  {
    _CDampingSparse.resize(_rank, _rank);
    cacheStaticDampingSparse();

    _ASparse.setSparsity(_CDampingSparse);
  }
  if (_solution.size() != _rank)
    _solution.resizeAndWipe(_rank);

  _timingBreakdown["Preamble"] += preambleTimer.timing();

  // backup the position for sampling
  _positionSample.equals(_tetMesh->x());
  
  // get the forces
  // R = (P_n - R(u_n))
  TIMER internalTimer;
  computeExternalForces();
  VECTOR R = _externalForces;
  R = _externalForces;
  R += _tetMesh->TET_MESH::generateInternalForces();
  VECTOR internal = _tetMesh->TET_MESH::generateInternalForces();
  _timingBreakdown["Internal Forces"] += internalTimer.timing();

  // start accum'ing in _temp,
  // switch to 'scratch' eventually to avoid a dynamic allocation
  // scratch = M * (2 * u_n - u_{n-1})
  TIMER RHSTimer;
  SPARSE_MATRIX M = _tetMesh->massMatrix();
  _temp.clearingAxpy(2.0, _position);
  _temp.axpy(-1.0, _positionOld);
  VECTOR scratch = M * _temp;

  // scratch += dt^2 * (P_n - R(u_n))
  scratch.axpy(_dt * _dt, R);

  // add damping
  // scratch += dt/2 * C * u_{n-1}
  VECTOR damping = _CDampingSparse * _positionOld;
  damping *= _dt / 2.0;

  // activate for FULLSPACE_INTEGRATOR only sim
  scratch.axpy(1.0, damping);
  _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

  // A = M + dt/2 * C
  TIMER ATimer;
  _ASparse.clear();
  _ASparse.axpy(1.0, M);
  _ASparse.axpy(_dt / 2.0, _CDampingSparse);
  Real PCGEps = pow(10.0, -_PCGDigits);
  _ASparse.eps() = PCGEps;
  _timingBreakdown["A Assembly"] += ATimer.timing();
  TIMER solve;
  //_ASparse.solveSPD(_solution, scratch);
  _ASparse.solveCG(_solution, scratch);
  //_ASparse.solveSPI(_solution, scratch);
  //_ASparse.solveGeneral(_solution, scratch);
  _timingBreakdown["Solver"] += solve.timing();

  /*
  // PARDISO version
  // A = M + dt/2 * C
  TIMER ATimer;

  // initialize the pardiso solver if necessary
  if (_APardiso.rows() == 0)
  {
    _APardiso.clear();
    _APardiso.axpy(1.0, M);
    _APardiso.axpy(_dt / 2.0, _CDampingSparse);
  }
  _timingBreakdown["A Assembly"] += ATimer.timing();
  TIMER solve;
  _APardiso.solveSPD(_solution, scratch);
  _timingBreakdown["Solver"] += solve.timing();
  */

  TIMER finalUpdate;

  //_accelerationOld = _acceleration;
  _acceleration = _position;
  _acceleration *= -2.0;
  _acceleration += _solution;
  _acceleration += _positionOld;
  _acceleration *= 1.0 / (_dt * _dt);

  //_velocityOld = _velocity;
  _velocity = _solution - _positionOld;
  _velocity *= 1.0 / (2.0 * _dt);
  
  _positionOld.swap(_position);
  _position.swap(_solution);

  // update node positions
  TIMER updateTimer;
  VECTOR& x = _tetMesh->x();
  x.copyInplace(_position);
  _tetMesh->TET_MESH::updateFullMesh();
  _timingBreakdown["Update Mesh"] += updateTimer.timing();

  // increment simulation time
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// Print out a detailed timing breakdown
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::printTimingBreakdown(string name, map<string,double>& timingBreakdown, double totalTime)
{
  // create an inverse map so that it will sort by time
  map<double, string> inverseMap;
  map<string, double>::iterator forwardIter;
  for (forwardIter = timingBreakdown.begin(); forwardIter != timingBreakdown.end(); forwardIter++)
    inverseMap[forwardIter->second] = forwardIter->first;

  // print the map out backwards since it sorts from least to greatest
  cout << name.c_str() << " TIMING BREAKDOWN: " << endl;
  cout << "===============================================================================" << endl;
  map<double,string>::reverse_iterator backwardIter;
  double totalSeen = 0.0;
  for (backwardIter = inverseMap.rbegin(); backwardIter != inverseMap.rend(); backwardIter++)
  {
    string name = (*backwardIter).second + string("                         ");
    name = name.substr(0,30);

    cout << "[" << (*backwardIter).first / totalTime * 100.0 << "%\t]: "
         << name.c_str() << "\t" << (*backwardIter).first / _totalSteps << "s / frame -- \t" << (*backwardIter).first << "s total" << endl;
    totalSeen += (*backwardIter).first;
  }
  Real misc = ((_totalTime - totalSeen) / _totalTime) * 100.0;
  cout << "[" << misc << "%\t]: " << "Misc. " << endl;
  cout << "===============================================================================" << endl;
  //cout << " Mean iterations: " << _ASparse.meanIterations() << endl;
  //cout << " Mean residual: " << _ASparse.meanResidual() << endl;
  cout << " Current FPS: " << _totalSteps / totalTime << endl;
  cout << " Current seconds / frame: " << totalTime / _totalSteps << endl;
  cout << " Total running time: " << totalTime << endl;
  cout << " Total Newton iterations: " << _totalNewtonStepsSeen << endl;
  cout << "===============================================================================" << endl;

  if (misc < 0.0)
  {
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << " BREAKDOWN ADDS UP TO MORE THAN 100! TIMERS ARE OVERLAPPING! " << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// cache the full static damping matrix
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::cacheStaticDampingFull()
{
  // make sure the matrix is the right size
  _CDamping.resizeAndWipe(_rank, _rank);

  // backup the current tet mesh state
  _solution = _tetMesh->x();
  
  // get the rest pose damping matrix in case we want to do just
  // linear Rayleigh damping
  _tetMesh->x().clear();
  _tetMesh->TET_MESH::generateF();

  // construct the constant damping matrix
  SPARSE_MATRIX massMatrix = _tetMesh->massMatrix();
  _CDamping.clear();
  for (int x = 0; x < massMatrix.rows(); x++)
    _CDamping(x,x) = _rayleighAlpha * massMatrix(x,x);
  //MATRIX& stiff = _tetMesh->TET_MESH::generateStiffnessMatrix();
  _tetMesh->TET_MESH::generateStiffnessMatrix();
  
  _CDamping.axpy(_rayleighBeta, _tetMesh->TET_MESH::generateStiffnessMatrix());

  // restore the tet mesh state
  _tetMesh->x() = _solution;
}

//////////////////////////////////////////////////////////////////////
// cache the sparse static damping matrix
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::cacheStaticDampingSparse()
{
  // make sure the matrix is the right size
  _CDampingSparse.resize(_rank, _rank);
  
  // backup the current tet mesh state
  _solution = _tetMesh->x();
  
  // get the rest pose damping matrix in case we want to do just
  // linear Rayleigh damping
  _tetMesh->x().clear();
  _tetMesh->TET_MESH::updateFullMesh();
  _tetMesh->TET_MESH::generateF();

  _CDampingSparse.clear();
  _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_CDampingSparse);
  _CDampingSparse *= _rayleighBeta;

  // construct the constant damping matrix
  SPARSE_MATRIX massMatrix = _tetMesh->massMatrix();
  for (int x = 0; x < massMatrix.rows(); x++)
    _CDampingSparse(x,x) += _rayleighAlpha * massMatrix(x,x);

  //_CDampingSparse *= _dt * 0.5f;

  /*
  // construct the constant damping matrix
  SPARSE_MATRIX massMatrix = _tetMesh->massMatrix();
  _CDampingSparse.clear();
  for (int x = 0; x < massMatrix.rows(); x++)
    _CDampingSparse(x,x) = _rayleighAlpha * massMatrix(x,x);
  cout << __FILE__ << " " << __LINE__ << " : " << endl; flush(cout);
  
  _CDampingSparse.axpy(_rayleighBeta, _tetMesh->TET_MESH::generateSparseStiffnessMatrix());
  cout << __FILE__ << " " << __LINE__ << " : " << endl; flush(cout);
  */

  // restore the tet mesh state
  _tetMesh->x() = _solution;
  _tetMesh->TET_MESH::updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
// process a mouse click event
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::click(VEC3F& point)
{
  // find and store the nearest surface point
  _clickedNode = _tetMesh->closestSurfaceNode(point);
  _clickedPosition = *_clickedNode;
  _draggedPosition = *_clickedNode;
}

//////////////////////////////////////////////////////////////////////
// Compute all external forces
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::computeExternalForces()
{
  // backup previous forces
  _externalForcesOld = _externalForces;

  // wipe previous force vectors
  _externalForces.clear();
 
  // process mouse input
  if (_clickedNode != NULL)
  {
    VEC3F dragForce = _draggedPosition - _clickedPosition;
    //dragForce.normalize();
    dragForce *= _forceMultiplier;
    _forceVectors.push_back(dragForce);
    _forceNodes.push_back(_clickedNode);
  }

  // make sure some forces are in the list
  if (_forceVectors.size() == 0) return;

  // accumulate the forces into the full force vector
  for (unsigned int x = 0; x < _forceVectors.size(); x++)
  {
    int index = _tetMesh->vertexID(_forceNodes[x]);
    if (index < _tetMesh->unconstrainedNodes())
    {
      _externalForces(3 * index)     += _forceMultiplier * _forceVectors[x][0];
      _externalForces(3 * index + 1) += _forceMultiplier * _forceVectors[x][1];
      _externalForces(3 * index + 2) += _forceMultiplier * _forceVectors[x][2];
    }
  }
  _forceVectors.clear();
  _forceNodes.clear();
}

//////////////////////////////////////////////////////////////////////
// poke the mesh in a deterministic way for debugging
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::poke()
{
  /*
  //VEC3F* vertex = _tetMesh->vertices(0);
  VEC3F* vertex = _tetMesh->vertices(1);
  //VEC3F direction(0.0f,1.0f,0.0f);
  //VEC3F direction(0.0f,100.0f,0.0f);
  //VEC3F direction(100.0f,0.0f,0.0f);
  //VEC3F direction(123.0f,456.0f,789.0f);
  VEC3F direction(1.254f, 0.234f, 1000.0f);
  */
  /*
  {
    VEC3F* vertex = _tetMesh->vertices(0);
    VEC3F direction(0.0, 100.0f, 0.0f);
    //VEC3F direction(0.0, 300.0f, 0.0f);
    //VEC3F direction(100.0, 0, 0.0f);

    _forceVectors.push_back(direction);
    _forceNodes.push_back(vertex);
  }

  {
    VEC3F* vertex = _tetMesh->vertices(9);
    VEC3F direction(0.0, -100.0f, 0.0f);
    //VEC3F direction(-100.0f, 0, 0.0f);

    _forceVectors.push_back(direction);
    _forceNodes.push_back(vertex);
  }
  */

  // try to induce a uniform rotation
  int size = _tetMesh->unconstrainedNodes();
  VEC3F omega;
  omega[0] = 0.0f;
  omega[1] = 1.0f;
  omega[2] = 0.0f;
  _tetMesh->computeCenterOfMass();
  VEC3F centerOfMass = _tetMesh->centerOfMass();
  Real theta = _dt * norm(omega);
  for (int x = 0; x < size; x++)
  {
    VEC3F* vertex = _tetMesh->vertices(x);
    VEC3F q = *vertex - centerOfMass;
    VEC3F cross = omega ^ q;
    Real magnitude = norm(cross);
    cross /= magnitude;
    
    Real qMagnitude = norm(q);
    VEC3F force = cross * qMagnitude * tan(theta);
    force *= _tetMesh->mass(x) / (_dt * _dt);

    _forceVectors.push_back(force);
    _forceNodes.push_back(vertex);
  }
}

//////////////////////////////////////////////////////////////////////
// Draw the clicked node
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::drawClickedNode()
{
  if (_clickedNode == NULL) return;
  glBegin(GL_POINTS);
#ifdef SINGLE_PRECISION    
    glVertex3fv(_clickedPosition);
#else
    glVertex3dv(_clickedPosition);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw a line from the dragged to clicked node
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::drawForceVector()
{
  if (_clickedNode == NULL) return;
  glPointSize(10.0f);
  glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
  glBegin(GL_POINTS);
#ifdef SINGLE_PRECISION    
    glVertex3fv(_clickedPosition);
    glVertex3fv(_draggedPosition);
#else
    glVertex3dv(_clickedPosition);
    glVertex3dv(_draggedPosition);
#endif
  glEnd();

  glLineWidth(4.0f);
  glBegin(GL_LINES);
#ifdef SINGLE_PRECISION    
    glVertex3fv(_clickedPosition);
    glVertex3fv(_draggedPosition);
#else
    glVertex3dv(_clickedPosition);
    glVertex3dv(_draggedPosition);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Add a gravity force to all the nodes
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::addGravity(VEC3F down, float gravity)
{
  for (int x = 0; x < _tetMesh->unconstrainedNodes(); x++)
  {
    VEC3F downCopy = down;
    downCopy.normalize();
    downCopy *= gravity;
    downCopy[0] *= _tetMesh->massMatrix()(3 * x, 3 * x);
    downCopy[1] *= _tetMesh->massMatrix()(3 * x + 1, 3 * x + 1);
    downCopy[2] *= _tetMesh->massMatrix()(3 * x + 2, 3 * x + 2);
    _forceVectors.push_back(downCopy);
    _forceNodes.push_back(_tetMesh->vertices(x));
  }
  /*
  down.normalize();
  down *= gravity;
  for (int x = 0; x < _tetMesh->unconstrainedNodes(); x++)
  {
    _forceVectors.push_back(downCopy);
    _forceNodes.push_back(_tetMesh->vertices(x));
  }
  */
}

//////////////////////////////////////////////////////////////////////
// Collision detection and response for a floor
// (Hacky -- for debugging purposes only)
//////////////////////////////////////////////////////////////////////
bool FULLSPACE_INTEGRATOR::collideWithFloor(VEC3F down, float floorPosition)
{
  bool collided = false;
    
  VEC3F up = down;
  up *= -1.0;
  for (int x = 0; x < _tetMesh->unconstrainedNodes(); x++)
  {
    VEC3F* vertex = _tetMesh->vertices(x);
    if ((*vertex)[1] < floorPosition)
    {
      float diff = floorPosition - (*vertex)[1];
      VEC3F response = up;
      //response *= diff * 100.0;
      response *= diff * 10.0;
      _forceVectors.push_back(response);
      _forceNodes.push_back(vertex);

      collided = true;
    }
  }
  return collided;
}

//////////////////////////////////////////////////////////////////////
// Compute rigid translation, rotation and velocities
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::printRigidComponents()
{
  // compute rigid velocity
  VEC3F rigidVelocity;
  int size = _tetMesh->unconstrainedNodes();
  for (int x = 0; x < size; x++)
  {
    int index = 3 * x;
    rigidVelocity[0] += _velocity(index);
    rigidVelocity[1] += _velocity(index + 1);
    rigidVelocity[2] += _velocity(index + 2);
  }
  rigidVelocity /= size;
  cout << " rigid velocity: " << rigidVelocity << endl;

  // compute rigid rotation
  VEC3F rigidRotation;
  _tetMesh->computeCenterOfMass();
  VEC3F centerOfMass = _tetMesh->centerOfMass();
  for (int x = 0; x < size; x++)
  {
    int index = 3 * x;
    VEC3F q = *(_tetMesh->vertices(x)) - centerOfMass;
    VEC3F v(_velocity(index), _velocity(index+1), _velocity(index+2));

    VEC3F rotate = q ^ v;
    rotate /= norm2(q);
    rigidRotation += rotate;
  }
  rigidRotation /= size;
  cout << " rigid rotation: " << rigidRotation << " " << norm(rigidRotation) << endl;
}

//////////////////////////////////////////////////////////////////////
// solve using the invertible finite element stiffness matrix
// from "Robust Quasisatic Finite Elements and Flesh Simulation"
// Teran, Sifakis, Irving, Fedkiw, SCA 2005
//
// There is a lot of magic that has to occur to do the matrix-vector
// multiply, so the whole CG solve is done explicitly here.
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::solveInvertibleStiffness(VECTOR& x, VECTOR& b, bool quasistatic)
{
  /*
  // get the temp CG arrays
  VECTOR& residual = _ASparse.residual();
  VECTOR& direction = _ASparse.direction();
  VECTOR& q = _ASparse.q();
  Real eps = _ASparse.eps();
  Real iterations = _ASparse.maxIterations();
 
  // resize if necessary
  if (residual.size() != _rank)
  {
    residual.resizeAndWipe(_rank);
    direction.resizeAndWipe(_rank);
    q.resizeAndWipe(_rank);
  }

  // non-invertible case
  x.clear();
  residual = b;

	// d = r
  direction.copyInplace(residual);
  
	// deltaNew = transpose(r) * r
  Real deltaNew = residual ^ residual;
  
	// delta0 = deltaNew
  Real delta0 = deltaNew;

	// While deltaNew > (eps^2) * delta0
	Real maxR = 2.0f * eps;
  int i = 0;
	while ((i < iterations) && (maxR > eps))
	{
    // multiply by K
    invertibleStiffnessMultiply(direction, q);

    // multiply and add the M and C components too
    if (!quasistatic)
      q += _ASparse * direction;

		// alpha = deltaNew / (transpose(d) * q)
    Real alpha = direction ^ q;
		if (fabs(alpha) > 0.0)
      alpha = deltaNew / alpha;

		// x = x + alpha * d
    x.axpy(alpha, direction);

		// r = r - alpha * q
    residual.axpy(-alpha, q);

		// deltaOld = deltaNew
		Real deltaOld = deltaNew;

		// deltaNew = transpose(r) * r
		deltaNew = residual ^ residual;

		// beta = deltaNew / deltaOld
		Real beta = deltaNew / deltaOld;

		// d = r + beta * d
    direction *= beta;
    direction += residual;

    // maxR = max(r);
    maxR = residual.maxValue();
    
		// i = i + 1
		i++;
  }
  cout << " iterations: " << i << " residual: " << maxR << endl;
  */
}

//////////////////////////////////////////////////////////////////////
// Do the matrix-vector multiply for solveInvertibleQuasistatic()
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::stiffnessMultiply(VECTOR& direction, VECTOR& answer)
{
  // wipe any old answer
  answer.clear();
  
  // if the material isn't INVERTIBLE, you did something really bad
  int totalTets = _tetMesh->totalTets();
  vector<TET>& tets = _tetMesh->tets();

  for (int x = 0; x < totalTets; x++)
  {
    TET& tet = tets[x];
    int materialIndex = tet.materialIndex();
    MATRIX3 F = tet.F();
    Real stiffness[81];
    VECTOR flatF = TET::flattenF(F);

    INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materials()[materialIndex]);
    material->stiffnessDensity(flatF.data(), stiffness);
    MATRIX stiffnessMatrix(9,9,stiffness);

    // copy current tet state
    int indices[4];
    for (int y = 0; y < 4; y++)
    {
      indices[y] = _tetMesh->vertexID(tet.vertices[y]);
      if (_tetMesh->isConstrained(indices[y]))
        indices[y] = -1;
    }

    VECTOR deltaX(12);
    for (int y = 0; y < 4; y++)
      if (indices[y] != -1)
      {
        int index = 3 * indices[y];
        deltaX(3 * y)     = direction(index);
        deltaX(3 * y + 1) = direction(index + 1);
        deltaX(3 * y + 2) = direction(index + 2);
      }

    MATRIX pFpu(9,12);
    material->computePFPu(tet, pFpu);

    VECTOR deltaF = pFpu * deltaX;
    VECTOR contraction = stiffnessMatrix * deltaF;
    MATRIX3 deltaP = TET::repackF(contraction);

    // compute the forces
    VEC3F forces[4];
    const VEC3F* b = tet.b();
    forces[0] = deltaP * b[0];
    forces[1] = deltaP * b[1];
    forces[2] = deltaP * b[2];
    forces[3] = deltaP * b[3];

    // add forces to corresponding vertices
    for (int y = 0; y < 4; y++)
    {
      if (indices[y] != -1)
      {
        int index = 3 * indices[y];
        answer(index)     += forces[y][0];
        answer(index + 1) += forces[y][1];
        answer(index + 2) += forces[y][2];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Do the matrix-vector multiply for solveInvertibleQuasistatic()
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::invertibleStiffnessMultiply(VECTOR& direction, VECTOR& answer)
{
  // wipe any old answer
  answer.clear();
  
  // if the material isn't INVERTIBLE, you did something really bad
  int totalTets = _tetMesh->totalTets();
  vector<TET>& tets = _tetMesh->tets();

  for (int x = 0; x < totalTets; x++)
  {
    TET& tet = tets[x];
    int materialIndex = tet.materialIndex();
    INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materials()[materialIndex]);

    // get cached stiffness matrix
    MATRIX& stiffnessMatrix = _stiffnesses[x];
    MATRIX3& U = _Us[x];
    MATRIX3& V = _Vs[x];

    // cache current tet vertex indices
    int indices[4];
    for (int y = 0; y < 4; y++)
    {
      indices[y] = _tetMesh->vertexID(tet.vertices[y]);
      if (_tetMesh->isConstrained(indices[y]))
        indices[y] = -1;
    }

    VECTOR deltaX(12);
    for (int y = 0; y < 4; y++)
      if (indices[y] != -1)
      {
        int index = 3 * indices[y];
        deltaX(3 * y)     = direction(index);
        deltaX(3 * y + 1) = direction(index + 1);
        deltaX(3 * y + 2) = direction(index + 2);
      }

    MATRIX pFpu(9,12);
    material->computePFPu(tet, pFpu);

    VECTOR deltaF = pFpu * deltaX;
    
    // rotate deltaF
    MATRIX3 repacked = TET::repackF(deltaF);
    repacked = U.transpose() * repacked * V;
    deltaF = TET::flattenF(repacked);

    VECTOR contraction = stiffnessMatrix * deltaF;
    MATRIX3 deltaP = TET::repackF(contraction);

    // rotate deltaP back
    deltaP = U * deltaP * V.transpose();

    // compute the forces
    VEC3F forces[4];
    const VEC3F* b = tet.b();
    forces[0] = deltaP * b[0];
    forces[1] = deltaP * b[1];
    forces[2] = deltaP * b[2];
    forces[3] = deltaP * b[3];

    // add forces to corresponding vertices
    for (int y = 0; y < 4; y++)
    {
      if (indices[y] != -1)
      {
        int index = 3 * indices[y];
        answer(index)     += forces[y][0];
        answer(index + 1) += forces[y][1];
        answer(index + 2) += forces[y][2];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Cache the deformation gradient diagonalizations
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::cacheDiagonalizations()
#if USING_OPENMP
{
  int totalTets = _tetMesh->totalTets();
  if (_Us.size() != totalTets)
  {
    _Us.resize(totalTets);
    _Vs.resize(totalTets);
    _Fhats.resize(totalTets);
    _stiffnesses.resize(totalTets);
  }

  // if the material isn't INVERTIBLE, you did something really bad
  vector<TET>& tets = _tetMesh->tets();
#pragma omp parallel
  { 
  int id  = omp_get_thread_num();
#pragma omp for  schedule(static)
  for (int x = 0; x < totalTets; x++)
  {
    int materialIndex = tets[x].materialIndex();
    //INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materials()[materialIndex]);
    INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materialCopies()[id][materialIndex]);
    
    MATRIX3 F = tets[x].F();
    MATRIX3 U;
    MATRIX3 Fhat;
    MATRIX3 V;
    Real stiffnessData[81];
    material->diagonalizeF(F, U, Fhat, V);
    material->stiffnessDensity(U, Fhat, V, stiffnessData);

    _Us[x] = U;
    _Vs[x] = V;
    _Fhats[x] = Fhat;

    MATRIX stiffness(9,9,stiffnessData);
    _stiffnesses[x] = stiffness;

#if 0
    // DEBUG -- sanity check
    MATRIX stiffnessSanity(9,9,stiffnessData);
    MATRIX eigenvectors9x9(9,9);
    VECTOR eigenvalues9x9(9);
    stiffnessSanity.eigensystem(eigenvalues9x9, eigenvectors9x9);

    bool anynegative = false;
    for ( int ii = 0; ii < 9; ii++ )
    {
      if ( eigenvalues9x9( ii ) > 0.0 ) anynegative = true;
    }
    if ( anynegative )
    {
      cout << "Stiffness eigenvalues: " << endl;
      cout << eigenvalues9x9 << endl << endl;
    }
#endif
  }
  } //OMP
}
#else
{
  int totalTets = _tetMesh->totalTets();
  if (_Us.size() != (unsigned int)totalTets)
  {
    _Us.resize(totalTets);
    _Vs.resize(totalTets);
    _Fhats.resize(totalTets);
    _stiffnesses.resize(totalTets);
  }

  // if the material isn't INVERTIBLE, you did something really bad
  vector<TET>& tets = _tetMesh->tets();
  for (int x = 0; x < totalTets; x++)
  {
    int materialIndex = tets[x].materialIndex();
    INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materials()[materialIndex]);
    
    MATRIX3 F = tets[x].F();
    MATRIX3 U;
    MATRIX3 Fhat;
    MATRIX3 V;
    Real stiffnessData[81];
    material->diagonalizeF(F, U, Fhat, V);
    material->stiffnessDensity(U, Fhat, V, stiffnessData);

    _Us[x] = U;
    _Vs[x] = V;
    _Fhats[x] = Fhat;

    MATRIX stiffness(9,9,stiffnessData);
    _stiffnesses[x] = stiffness;
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// DEBUG only for stiffness assembly
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::debugStiffness()
{
  TET& tet = (_tetMesh->tets())[28];
  int materialIndex = tet.materialIndex();
  INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materials()[materialIndex]);
  VECTOR deltaX(12);
  for (int x = 0; x < 12; x++)
    deltaX(x) = x + 1;
  VEC3F forces[4];
  const VEC3F* b = tet.b();
  cout << "Tet under consideration: " << endl;
  for (int x = 0; x < 4; x++)
    cout << "Vertex " << x << ": " << *(tet.vertices[x]) << endl;

  // do stiffness the gold standard way
  MATRIX goldStiffness = material->stiffness(tet);
  VECTOR goldF = goldStiffness* deltaX;
  cout << " gold standard deltaF: " << goldF << endl;

  // do stiffness the Teran way
  MATRIX3 F = tet.F();
  VECTOR flatF = TET::flattenF(F);
  Real stiffness[81];
  material->stiffnessDensity(flatF.data(), stiffness);
  MATRIX stiffnessMatrix(9,9,stiffness);
  stiffnessMatrix *= -1.0;
  MATRIX pFpu(9,12);
  material->computePFPu(tet, pFpu);
  VECTOR deltaF = pFpu * deltaX;
  VECTOR contraction = stiffnessMatrix * deltaF;
  MATRIX3 deltaP = TET::repackF(contraction);

  forces[0] = deltaP * b[0];
  forces[1] = deltaP * b[1];
  forces[2] = deltaP * b[2];
  forces[3] = deltaP * b[3];
  cout << " Teran deltaF: " << endl;
  for (int x = 0; x < 4; x++)
    cout << forces[x] << endl;
  exit(0);
}

//////////////////////////////////////////////////////////////////////
// Stream state to already open file
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::writeState(FILE* file)
{
  _position.write(file);
  _positionOld.write(file);
  _velocity.write(file);
  _velocityOld.write(file);
  _acceleration.write(file);
  _accelerationOld.write(file);
  _residual.write(file);
  _temp.write(file);
  _solution.write(file);
  _externalForces.write(file);
}

//////////////////////////////////////////////////////////////////////
// Dump state to a file
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::writeState(string filename)
{
  cout << " Dumping state to file " << filename.c_str() << endl;
  FILE* file = fopen(filename.c_str(), "wb");
  writeState(file);
  /*
  _position.write(file);
  _positionOld.write(file);
  _velocity.write(file);
  _velocityOld.write(file);
  _acceleration.write(file);
  _accelerationOld.write(file);
  _residual.write(file);
  _temp.write(file);
  _solution.write(file);
  _externalForces.write(file);
  */
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Read state from already open file
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::readState(FILE* file)
{
  _position.read(file);
  _positionOld.read(file);
  _velocity.read(file);
  _velocityOld.read(file);
  _acceleration.read(file);
  _accelerationOld.read(file);
  _residual.read(file);
  _temp.read(file);
  _solution.read(file);
  _externalForces.read(file);
}

//////////////////////////////////////////////////////////////////////
// Read state from file
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::readState(string filename)
{
  //cout << " Reading state from file " << filename.c_str() << endl;
  FILE* file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << " Error opening file: " << filename.c_str() << endl;
    exit(0);
  }
  readState(file);
  /*
  _position.read(file);
  _positionOld.read(file);
  _velocity.read(file);
  _velocityOld.read(file);
  _acceleration.read(file);
  _accelerationOld.read(file);
  _residual.read(file);
  _temp.read(file);
  _solution.read(file);

  // see if it's a newer file
  if (!feof(file))
    _externalForces.read(file);
    */
  fclose(file);

  // update tet state
  VECTOR& x = _tetMesh->TET_MESH::x();

  if (x.size() != _position.size())
    x.resizeAndWipe(_position.size());

  //cout << " x size: " << x.size() << endl;
  //cout << " position size: " << _position.size() << endl;

  x.equals(_position);
  _tetMesh->TET_MESH::updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
// compute the current force residual
//////////////////////////////////////////////////////////////////////
VECTOR& FULLSPACE_INTEGRATOR::computeForceResidual(bool invertible)
{
  SPARSE_MATRIX dummy;

  _tetMesh->x() = this->position();
  generateImplicitMatrices(dummy, invertible);

  SPARSE_MATRIX M = _tetMesh->massMatrix();
  SPARSE_MATRIX C = this->dampingMatrix();
  VECTOR acceleration = this->acceleration();
  VECTOR velocity = this->velocity();
  VECTOR position = this->position();
  VECTOR R = _tetMesh->TET_MESH::internalForce();
  VECTOR F = this->externalForces();

  /*
  cout << " M sum2: " << M.sum2() << endl;
  cout << " C sum2: " << C.sum2() << endl;
  cout << " R norm: " << R.norm2() << endl;
  cout << " F norm: " << F.norm2() << endl;
  cout << " acceleration norm: " << acceleration.norm2() << endl;
  cout << " velocity norm:     " << velocity.norm2() << endl;
  cout << " position norm:     " << position.norm2() << endl;
  */

  _residual = (M * acceleration) + (C * velocity) - R - F;
  return _residual;
}

//////////////////////////////////////////////////////////////////////
// compute the friction forces and store them in _frictionForces
//////////////////////////////////////////////////////////////////////
VECTOR& FULLSPACE_INTEGRATOR::computeFrictionForces()
{
  if (_frictionForces.size() != _rank)
    _frictionForces.resizeAndWipe(_rank);
  else
    _frictionForces.clear();

  // get the force for each vertex
  int frictionVertices = 0;
  for (int x = 0; x < _position.size() / 3; x++)
  {
    VEC3F& restVertex = *(_tetMesh->restVertices(x));
    
    VEC3F vertex = restVertex;
    vertex[0] += _position[3 * x];
    vertex[1] += _position[3 * x + 1];
    vertex[2] += _position[3 * x + 2];

    VEC3F velocity;
    velocity[0] = _velocity[3 * x];
    velocity[1] = _velocity[3 * x + 1];
    velocity[2] = _velocity[3 * x + 2];

    VEC3F force;
    if (_collisionForces.size() == _rank)
    {
      force[0] = _collisionForces[3 * x];
      force[1] = _collisionForces[3 * x + 1];
      force[2] = _collisionForces[3 * x + 2];
    }

    VEC3F frictionForce;
    for (unsigned int y = 0; y < _collisionResponses.size(); y++)
      frictionForce += _collisionResponses[y]->friction(vertex, velocity, force);

    if (norm(frictionForce) > 0.0)
      frictionVertices++;

    Real area = _tetMesh->surfaceArea(x);
    _frictionForces[3 * x] += frictionForce[0] * area / _dt;
    _frictionForces[3 * x + 1] += frictionForce[1] * area / _dt;
    _frictionForces[3 * x + 2] += frictionForce[2] * area / _dt;
  }

  return _frictionForces;
}

//////////////////////////////////////////////////////////////////////
// compute the collision forces and store them in _collisionForces
//////////////////////////////////////////////////////////////////////
VECTOR FULLSPACE_INTEGRATOR::computeCollisionResidual(vector<pair<VEC3F*, SURFACE*> >& collisionVertices)
{
  VECTOR collisionResidual(_rank);

  // get the force for each point
  int unconstrainedSize = _tetMesh->unconstrainedNodes();
  for (unsigned int x = 0; x < collisionVertices.size(); x++)
  {
    VEC3F* vertex = collisionVertices[x].first;
    SURFACE* surface = collisionVertices[x].second;
    int vertexID = _tetMesh->vertexID(vertex);
    if (vertexID >= unconstrainedSize) continue;

    // NEW: get the surface area of the collision vertex
    Real surfaceArea = _tetMesh->surfaceArea(vertexID);

    // retrieve the velocity for damping
    VEC3F vertexVelocity;
    vertexVelocity[0] = _velocity[3 * vertexID];
    vertexVelocity[1] = _velocity[3 * vertexID + 1];
    vertexVelocity[2] = _velocity[3 * vertexID + 2];

    VEC3F& restVertex= *(_tetMesh->restVertices(vertexID));
    VEC3F displace;
    displace[0] = _position[3 * vertexID];
    displace[1] = _position[3 * vertexID + 1];
    displace[2] = _position[3 * vertexID + 2];
    VEC3F updatedPosition = restVertex + displace;

    // retrieve the damped force
    VEC3F force = surface->force(updatedPosition, vertexVelocity);
    force *= surfaceArea;
    force *= -1.0;

    collisionResidual[3 * vertexID] = force[0];
    collisionResidual[3 * vertexID + 1] = force[1];
    collisionResidual[3 * vertexID + 2] = force[2];
  }

  return collisionResidual;  
}

//////////////////////////////////////////////////////////////////////
// compute the collision forces and store them in _collisionForces
//////////////////////////////////////////////////////////////////////
VECTOR& FULLSPACE_INTEGRATOR::computeCollisionForces()
{
  if (_collisionForces.size() != _rank)
    _collisionForces.resizeAndWipe(_rank);
  else
    _collisionForces.clear();

  // get the force for each vertex
  for (int x = 0; x < _position.size() / 3; x++)
  {
    VEC3F& restVertex = *(_tetMesh->restVertices(x));
    
    VEC3F vertex = restVertex;
    vertex[0] += _position[3 * x];
    vertex[1] += _position[3 * x + 1];
    vertex[2] += _position[3 * x + 2];

    VEC3F velocity;
    velocity[0] = _velocity[3 * x];
    velocity[1] = _velocity[3 * x + 1];
    velocity[2] = _velocity[3 * x + 2];

    VEC3F force;
    for (unsigned int y = 0; y < _collisionResponses.size(); y++)
      force += _collisionResponses[y]->force(vertex, velocity);

    Real mass = 1.0 / _tetMesh->mass(x);
    _collisionForces[3 * x] += force[0] * mass;
    _collisionForces[3 * x + 1] += force[1] * mass;
    _collisionForces[3 * x + 2] += force[2] * mass;
  }

  return _collisionForces;
}

//////////////////////////////////////////////////////////////////////
// compute the collision Jacobians and store them in 
// _collisionForceJacobians
//////////////////////////////////////////////////////////////////////
SPARSE_MATRIX& FULLSPACE_INTEGRATOR::computeCollisionForceJacobians()
{
  if (_collisionForceJacobians.rows() != _rank ||
      _collisionForceJacobians.cols() != _rank)
    _collisionForceJacobians.resize(_rank, _rank);
  else
    _collisionForceJacobians.clear();

  // set the newmark alphas
  for (unsigned int x = 0; x < _collisionResponses.size(); x++)
    //_collisionResponses[x]->newmarkAlpha() = fabs(_alpha[3]);
    _collisionResponses[x]->newmarkAlpha() = _alpha[3];

  // make sure the velocities aren't stale
  updateVelocities();

  // get the Jacobian for each vertex
  for (int x = 0; x < _position.size() / 3; x++)
  {
    VEC3F& restVertex = *(_tetMesh->restVertices(x));
   
    VEC3F displacement;
    displacement[0] = _position[3 * x];
    displacement[1] = _position[3 * x + 1];
    displacement[2] = _position[3 * x + 2];
    VEC3F vertex = restVertex + displacement;

    VEC3F velocity;
    velocity[0] = _velocity[3 * x];
    velocity[1] = _velocity[3 * x + 1];
    velocity[2] = _velocity[3 * x + 2];

    MATRIX3 jacobian;
    for (unsigned int y = 0; y < _collisionResponses.size(); y++)
      jacobian += _collisionResponses[y]->forceJacobian(vertex, velocity);

    Real mass = 1.0 / _tetMesh->mass(x);
    jacobian *= mass;

    // place it along the diagonal in the sparse matrix
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        _collisionForceJacobians(3 * x + i, 3 * x + j) += jacobian(i,j);
  }

  return _collisionForceJacobians;
}

//////////////////////////////////////////////////////////////////////
// compute collision jacobian based on collision pairs
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::addCollisionForceJacobians(vector<pair<VEC3F*, SURFACE*> >& collisionVertices)
{
  // build a list of the unique vertices in collision
  vector<TRIANGLE*> faces = _tetMesh->explicitSurfaceFaces();

  if (_systemMatrix.rows() == 0)
    _systemMatrix.resize(_rank, _rank);
  else
    _systemMatrix.clear();

  // get the force for each point and add it to the residual
  int unconstrainedSize = _tetMesh->unconstrainedNodes();
  for (unsigned int x = 0; x < collisionVertices.size(); x++)
  {
    VEC3F* vertex = collisionVertices[x].first;
    SURFACE* surface = collisionVertices[x].second;
    int vertexID = _tetMesh->vertexID(vertex);
    if (vertexID >= unconstrainedSize) continue;

    // need to debug the sphere implementation at some point
    if (surface->type().compare("SPHERE") == 0)
      continue;

    // retrieve the velocity for damping
    VEC3F vertexVelocity;
    vertexVelocity[0] = _velocity[3 * vertexID];
    vertexVelocity[1] = _velocity[3 * vertexID + 1];
    vertexVelocity[2] = _velocity[3 * vertexID + 2];

    // retrieve position for stiffness
    VEC3F& restVertex= *(_tetMesh->restVertices(vertexID));
    VEC3F displace;
    displace[0] = _position[3 * vertexID];
    displace[1] = _position[3 * vertexID + 1];
    displace[2] = _position[3 * vertexID + 2];
    VEC3F updatedPosition = restVertex + displace;

    MATRIX springJacobian = surface->springJacobian(updatedPosition);
    MATRIX dampingJacobian = surface->dampingJacobian(updatedPosition, vertexVelocity);
    //springJacobian *= _alpha[4];
    dampingJacobian *= _alpha[4];
    MATRIX jacobian = springJacobian - dampingJacobian;
    //MATRIX jacobian = springJacobian;

    /*
    // DEBUG: see if the Jacobian looks correct
    //surface->verifySpringJacobian(*vertex);

    MATRIX reducedJacobian = vertexSubBasis ^ jacobian * vertexSubBasis;
    _A += reducedJacobian;
    */

    // place it along the diagonal in the sparse matrix
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
        if (_ASparse.rows() > 0)
          _ASparse.add(3 * vertexID + i, 3 * vertexID + j, jacobian(i,j));
        _systemMatrix.add(jacobian(i,j), 3 * vertexID + i, 3 * vertexID + j);
      }
  }
}

//////////////////////////////////////////////////////////////////////
// update the velocity from the current position --
// this is necessary for the proper collision forces to be computed
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::updateVelocities()
{
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);
}

//////////////////////////////////////////////////////////////////////
// apply constraints -- ie "frozen vertices"
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::applyMatrixConstraints()
{
  for (unsigned int x = 0; x < _frozenVertices.size(); x++)
  {
    // get the index of the vertex
    int index = _tetMesh->vertexID(_frozenVertices[x]);
    
    // zero out the matrix entries for that vertex
    _ASparse.eraseRowColumn(3 * index);
    _ASparse.eraseRowColumn(3 * index + 1);
    _ASparse.eraseRowColumn(3 * index + 2);

    // apply the identity along the diagonal
#if _WIN32    
    _ASparse(3 * index, 3 * index) = 1.0;
    _ASparse(3 * index + 1, 3 * index + 1) = 1.0;
    _ASparse(3 * index + 2, 3 * index + 2) = 1.0;
#else
    _ASparse.set(3 * index, 3 * index, 1.0);
    _ASparse.set(3 * index + 1, 3 * index + 1, 1.0);
    _ASparse.set(3 * index + 2, 3 * index + 2, 1.0);
#endif
  }
}

//////////////////////////////////////////////////////////////////////
// apply constraints -- ie "frozen vertices"
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::applyResidualConstraints()
{
  for (unsigned int x = 0; x < _frozenVertices.size(); x++)
  {
    // get the index of the vertex
    int index = _tetMesh->vertexID(_frozenVertices[x]);
    
    // fix up the RHS so that it works out to zero
    _residual[3 * index] = 0.0;
    _residual[3 * index + 1] = 0.0;
    _residual[3 * index + 2] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////
// Reset everything to zero
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::reset()
{
  _position *= 0;
  _velocity *= 0;
  _acceleration *= 0;
  _positionOld *= 0;
  _velocityOld *= 0;
  _accelerationOld *= 0;
  _residual *= 0;
  _temp *= 0;
  _externalForces *= 0;
  _externalForcesOld *= 0;
  _collisionForces *= 0;
  _frictionForces *= 0;
  _positionSample *= 0;
}

//////////////////////////////////////////////////////////////////////
// needed when calling the interface stiffness auto-tuner
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::generateImplicitMatrices(SPARSE_MATRIX& ASparse, bool invertible)
{
  int rank = _tetMesh->TET_MESH::rank();
  int totalTets = _tetMesh->tets().size();
  _tetMesh->x().resizeAndWipe(rank);
  _tetMesh->TET_MESH::F().resizeAndWipe(totalTets * 9);
  _tetMesh->TET_MESH::internalForce().resizeAndWipe(rank);

  // allocate matrices if necessary -- should only fire
  // on the first step
  if (_CDampingSparse.rows() != _rank)
  {
    _CDampingSparse.resize(_rank, _rank);
    cacheStaticDampingSparse();
  }

  // get the mass matrix
  SPARSE_MATRIX& M = _tetMesh->TET_MESH::massMatrix();

  // update tet mesh with the new q vector
  VECTOR& x = _tetMesh->TET_MESH::x();

  // an assignment is overly expensive -- the "=" forces an allocation
  x.equals(_position);
  _tetMesh->TET_MESH::updateFullMesh();

  // diagonalize
  _tetMesh->TET_MESH::generateF();
  if (invertible)
    cacheDiagonalizations();

  // get the internal forces
  //VECTOR& R = (invertible) ? _tetMesh->TET_MESH::generateInternalForces(_Us, _Fhats, _Vs) :
  //                           _tetMesh->TET_MESH::generateInternalForces();
  VECTOR& R = _tetMesh->TET_MESH::generateInternalForces();

  // compute the LHS of the residual:
  // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
  _temp.clearingAxpy(-_alpha[2], _accelerationOld);
  _temp.axpy(-_alpha[1], _velocityOld);
  _temp.axpy(-_alpha[0], _positionOld);
  _temp.axpy(_alpha[0], _position);

  // do an explicit multiply because of the sparse matrix
  for (int i = 0; i < M.rows(); i++)
    _temp(i) *= M(i,i);
  
  // compute the RHS of the residual:
  // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
  _residual.clearingAxpy(_alpha[5], _accelerationOld);
  _residual.axpy(_alpha[4], _velocityOld);
  _residual.axpy(-_alpha[3], _positionOld);
  _residual.axpy(_alpha[3], _position);
  _residual = _CDampingSparse * _residual;

  // assemble full residual: LHS + RHS + R - F
  _residual += _temp;
  _residual -= R;
  _residual -= _externalForces;

  SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);

  // assemble system matrix A
  if (ASparse.rows() != _rank)
    ASparse.resize(_rank, _rank);
  ASparse.equals(K);
  for (int i = 0; i < M.rows(); i++)
    ASparse(i,i) += _alpha[0] * M(i,i);
  ASparse.axpy(_alpha[3], _CDampingSparse);

  /*
  VECTOR residual = _CDampingSparse * _velocity;
  residual -= R;
  residual -= _externalForces;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " full acceleration residual: " << residual << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " full position residual: " << _residual << endl;
  */
}

//////////////////////////////////////////////////////////////////////
// Form the A and b so that they can be passed out for solution
// by the partitioned integrator -- quasistatic case
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::initializeQuasistaticStep()
{
  TIMER preamble;

  // accumulate external forces
  computeExternalForces();
  
  // copy the current values into the old values
  _positionOld = _position;

  _timingBreakdown["Preamble"] += preamble.timing();
}

//////////////////////////////////////////////////////////////////////
// Update the mesh after the partitioned quasistatic step
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::finalizeQuasistaticStep()
{
  // update node positions
  //_tetMesh->updateSurfaceMesh();
  _tetMesh->updateFullMesh();
  _externalForces.clear();
  _totalSteps++;
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an quasistatic step -- should be called 
// right after initializeQuasistaticStep
//
// this has been forked out into its own function in case we want
// to do more than one Newton step
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::generateQuasistaticMatrices(map<string, double>& timingBreakdown)
{
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

  //TIMER RHSTimer;
  _residual.equals(_externalForces);
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
// Implicit Integration, using acceleration as the primary variable
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::stepImplicitAcceleration()
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
  computeExternalForces();
  
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
  cout << " Solving acceleration-level implicit step " << _totalSteps << endl;
  cout << "==================================================" << endl;
  
  while (step < _maxNewtonSteps && maxR > eps)
  //while (step < _maxNewtonSteps && maxR > initialR * eps)
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
    //_CDampingSparse = _rayleighAlpha * M;
    //_CDampingSparse += _rayleighBeta * K;

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

//////////////////////////////////////////////////////////////////////
// change the direction of the gravity vector
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::changeGravity(VEC3F gravityDown, Real gravityMagnitude) 
{ 
  _gravityDown = gravityDown; 
  _gravityMagnitude = gravityMagnitude;

  VEC3F gravity = gravityMagnitude * gravityDown;
  _gravity *= 0.0;
  for (int x = 0; x < _gravity.size() / 3; x++)
  {
    _gravity[3 * x] = gravity[0];
    _gravity[3 * x + 1] = gravity[1];
    _gravity[3 * x + 2] = gravity[2];
  }
}

//////////////////////////////////////////////////////////////////////
// Implicit integration with a skinned mesh
//////////////////////////////////////////////////////////////////////
bool FULLSPACE_INTEGRATOR::stepSkinnedImplicit()
{
  TIMER total;

  TIMER preambleTimer;
  // allocate matrices if necessary -- should only fire
  // on the first step
  if (_CDampingSparse.rows() != _rank)
  {
    _CDampingSparse.resize(_rank, _rank);
    cacheStaticDampingSparse();

    _ASparse.setSparsity(_CDampingSparse);
  }

  if (_solution.size() != _rank)
    _solution.resizeAndWipe(_rank);
 
  // get the mass matrix
  SPARSE_MATRIX& M = _tetMesh->TET_MESH::massMatrix();

  cout << " Mass: " << M.sum() << endl;

  // accumulate external forces
  computeExternalForces();

  cout << "=================================================" << endl;
  cout << " Solving skinned implicit timestep " << _totalSteps << endl;
  cout << "==================================================" << endl;
  cout << " External forces: " << _externalForces.norm2() << endl;

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;
  _timingBreakdown["Preamble"] += preambleTimer.timing();

  // load in the skinned guess from the tet mesh
  _position = _tetMesh->TET_MESH::x();

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  int step = 0;
  while (step < _maxNewtonSteps && maxR > eps)
  {
    // update tet mesh with the new q vector
    TIMER updateMesh;
    VECTOR& x = _tetMesh->TET_MESH::x();

    // an assignment is overly expensive -- the "=" forces an allocation
    x.equals(_position);
    _positionSample.equals(_position);

    _tetMesh->TET_MESH::updateFullMesh();
    _timingBreakdown["Update Mesh"] += updateMesh.timing();

    // diagonalize
    TIMER diagonal;
    cacheDiagonalizations();
    _timingBreakdown["Diagonalization and Stiffness Density"] += diagonal.timing();

    // get the internal forces
    TIMER internalTimer;
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces(_Us, _Fhats, _Vs);
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // compute the LHS of the residual:
    // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
    TIMER RHSTimer;
    _temp.clearingAxpy(-_alpha[2], _accelerationOld);
    _temp.axpy(-_alpha[1], _velocityOld);
    _temp.axpy(-_alpha[0], _positionOld);
    _temp.axpy(_alpha[0], _position);

    // do an explicit multiply because of the sparse matrix
    for (int i = 0; i < M.rows(); i++)
      _temp(i) *= M(i,i);

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    _residual.clearingAxpy(_alpha[5], _accelerationOld);
    _residual.axpy(_alpha[4], _velocityOld);
    _residual.axpy(-_alpha[3], _positionOld);
    _residual.axpy(_alpha[3], _position);
    _residual = _CDampingSparse * _residual;

    // assemble full residual: LHS + RHS + R - F
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;

    maxR = _residual.norm2();
    step++;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

    cout << " Newton-Raphson residual " << step << ": " << maxR;
    cout << endl; flush(cout);
    if (maxR < eps) break;
   
    // assemble system matrix A
    TIMER genSparseTimer;
    SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);
    _timingBreakdown["Generate Sparse Stiffness Matrix"] += genSparseTimer.timing();
    TIMER copyPetscTimer;
    _ASparse.equals(K);
    _timingBreakdown["Petsc Stiffness copy"] += copyPetscTimer.timing();

    TIMER addMassTimer;
    for (int i = 0; i < M.rows(); i++)
      _ASparse.add(i,i, _alpha[0] * M(i,i));
    _timingBreakdown["Adding mass to stiffness matrix"] += addMassTimer.timing();
    TIMER addDampingTimer;
    _ASparse.axpy(_alpha[3], _CDampingSparse);
    _timingBreakdown["Adding damping to stiffness matrix"] += addDampingTimer.timing();

    // solve with CG
    TIMER solverTimer;
    Real PCGEps = pow(10.0, -_PCGDigits);
    _ASparse.eps() = PCGEps;
    _ASparse.usePCG();
    bool converged = _ASparse.solveCG(_solution, _residual);
    _timingBreakdown["Solver"] += solverTimer.timing();

    if (!converged)
    {
      cout << " Is mesh inverted? " << _tetMesh->inverted(true) << endl;
      //if (_substeppingEnabled && !substepping) break;
      return false;
    }

    // update positions
    _position -= _solution;
  }

  if (maxR > eps)
  {
    _time += _dt;
    _totalSteps++;
    _totalTime += total.timing();
    return false;
  }

	// update velocity
  TIMER finalUpdate;
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);

	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], _accelerationOld);
  _acceleration.axpy(-_alpha[1], _velocityOld);
  _acceleration.axpy(-_alpha[0], _positionOld);
  _acceleration.axpy(_alpha[0], _position);

  // update node positions
  VECTOR& x = _tetMesh->TET_MESH::x();
  x.copyInplace(_position);
  _tetMesh->TET_MESH::updateFullMesh();

  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment counters
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();

  return true;
}

//////////////////////////////////////////////////////////////////////
// Quasistatic implicit
//////////////////////////////////////////////////////////////////////
bool FULLSPACE_INTEGRATOR::stepSkinnedQuasistaticWithCollisions(vector<pair<SURFACE*, int> > collisionPairs)
{
  TIMER total;

  // build collision lists
  vector<pair<VEC3F*, SURFACE*> > collisionVertices;
  vector<SURFACE*> collisionSurfaces;
  buildCollisionList(collisionPairs, collisionVertices, collisionSurfaces);

  TIMER preambleTimer;
  if (_solution.size() != _rank)
    _solution.resizeAndWipe(_rank);
 
  // accumulate external forces
  computeExternalForces();
  _timingBreakdown["Preamble"] += preambleTimer.timing();
 
  // copy the current values into the old values
  _positionOld = _position;

  _position = _tetMesh->x();

  VECTOR bestPosition = _position;

  // do Newton-Raphson
  Real eps = _solverEps;
  Real bestResidual = eps * 100;
  Real maxR = eps * 10;
  int step = 0;
  cout << "==================================================" << endl;
  cout << " Solving skinned quasistatic with collisions timestep " << _totalSteps << endl;
  cout << "==================================================" << endl;
  while (step < _maxNewtonSteps && maxR > eps)
  {
    TIMER updateMesh;
    // update tet mesh with the new q vector
    VECTOR& x = _tetMesh->x();

    // an assignment is overly expensive -- the "=" forces an allocation
    //q = _position;
    x.equals(_position);
    _positionSample.equals(_position);

    _tetMesh->TET_MESH::updateFullMesh();
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // caching diagonalizations
    INVERTIBLE::inversions() = 0;
    TIMER diagonalTimer;
    cacheDiagonalizations();
    _timingBreakdown["Diagonalization"] += diagonalTimer.timing();  
    Real percent = 100.0 * INVERTIBLE::inversions() / (_tetMesh->tets().size() * 9);
    static Real maxSeen = percent;
    static Real minSeen = percent;
    if (percent > maxSeen) maxSeen = percent;
    if (percent < minSeen) minSeen = percent;

    /*
    cout << " tet inversions: " << INVERTIBLE::inversions() << " of " << _tetMesh->dofs() * 9 << "(" << percent << "%)" << endl;
    cout << " max seen: " << maxSeen << endl;
    cout << " min seen: " << minSeen << endl;
    */
      
    // get the reduced internal forces
    TIMER internalTimer;
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces(_Us, _Fhats, _Vs);
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the sparse stiffness matrix
    TIMER stiffnessTimer;
    SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);
    _timingBreakdown["Stiffness Assembly"] += stiffnessTimer.timing();

    // if this is the first time, set the sparsity of the PETSC matrix
    if (_ASparse.rows() != _rank)
      _ASparse.setSparsity(K);

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    TIMER RHSTimer;
    _residual = _externalForces;
    _residual += R;

    // add collision responses (old version)
    _velocity.clear();
    //computeCollisionForces();
    //computeCollisionForceJacobians();
    //_residual -= _collisionForces;
   
    // NEW: add in collision residuals
    _residual += computeCollisionResidual(collisionVertices);

    // apply constraints
    applyResidualConstraints();

    maxR = _residual.norm2();
    cout << " Newton-Raphson residual " << step << ": " << maxR << endl;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

    /*
    // DEBUG -- breakdown detector
    if (step == 0) oldResidual = maxR;
    Real relativeProgress = fabs((oldResidual - maxR) / oldResidual);
    oldResidual = maxR;
    cout << " Relative progress: " << relativeProgress << endl;
    if (relativeProgress < 1e-2)
      badProgress++;
    else
      badProgress = 0;

    if (badProgress == 10)
    {
      cout << " Solver breakdown!!!!!! " << endl;
      break;
    }
    */

    if (maxR < eps) break;

    // solve with explicit stiffness matrix A
    TIMER ACopy;
    _ASparse.equals(K);
    _timingBreakdown["A Copy"] += ACopy.timing();

    // NEW: add in collision Jacobians
    addCollisionForceJacobians(collisionVertices);

    /*
    // add force jacobians and constraints (old version)
    TIMER addForceTimer;
    _ASparse.axpy(1.0, _collisionForceJacobians);
    applyMatrixConstraints();
    _timingBreakdown["Adding collision forces to matrix"] += addForceTimer.timing();
    */

    TIMER solverTimer;
    Real PCGEps = pow(10.0, -_PCGDigits);
    _ASparse.eps() = PCGEps;

    // DEBUG: force to PCG
    _ASparse.usePCG();
    _ASparse.useICC();
    //bool converged = _ASparse.solveCG(_solution, _residual);
    _ASparse.solveCG(_solution, _residual);
    _timingBreakdown["Solver"] += solverTimer.timing();

    // update positions
    _position += _solution;

    if (step == 0 || maxR < bestResidual)
    {
      bestResidual = maxR;
      bestPosition = _position;
    }
    step++;
  }

  if (maxR > bestResidual)
  {
    cout << " Rolling solution back to " << bestResidual << endl;
    _position = bestPosition;
  }

  // DEBUG
  static int maxStepsSeen = 0;
  if (step > maxStepsSeen)
    maxStepsSeen = step;
  cout << "Max Newton steps seen: " << maxStepsSeen << endl;
  static int totalStepsSeen = 0;
  totalStepsSeen += step;
  cout << "Mean Newton steps seen: " << (float)totalStepsSeen / (_totalSteps + 1) << endl;
  static Real totalResidual = 0;
  totalResidual += maxR;
  cout << "Mean residual seen: " << totalResidual / (_totalSteps + 1) << endl;

  TIMER finalUpdate;
  // update node positions
  VECTOR& x = _tetMesh->x();
  x.copyInplace(_position);
  _tetMesh->TET_MESH::updateFullMesh();
  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment counters
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();

  return (step == _maxNewtonSteps);
}

//////////////////////////////////////////////////////////////////////
// Quasistatic implicit
//////////////////////////////////////////////////////////////////////
bool FULLSPACE_INTEGRATOR::stepSkinnedQuasistatic()
{
  TIMER total;

  TIMER preambleTimer;
  if (_solution.size() != _rank)
    _solution.resizeAndWipe(_rank);
 
  // accumulate external forces
  computeExternalForces();
  _timingBreakdown["Preamble"] += preambleTimer.timing();
 
  // copy the current values into the old values
  _positionOld = _position;

  _position = _tetMesh->x();

  VECTOR bestPosition = _position;

  // do Newton-Raphson
  Real eps = _solverEps;
  Real bestResidual = eps * 100;
  Real maxR = eps * 10;
  int step = 0;
  cout << "==================================================" << endl;
  cout << " Solving skinned quasistatic timestep " << _totalSteps << endl;
  cout << "==================================================" << endl;
  while (step < _maxNewtonSteps && maxR > eps)
  {
    TIMER updateMesh;
    // update tet mesh with the new q vector
    VECTOR& x = _tetMesh->x();

    // an assignment is overly expensive -- the "=" forces an allocation
    //q = _position;
    x.equals(_position);
    _positionSample.equals(_position);

    _tetMesh->TET_MESH::updateFullMesh();
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // caching diagonalizations
    INVERTIBLE::inversions() = 0;
    TIMER diagonalTimer;
    cacheDiagonalizations();
    _timingBreakdown["Diagonalization"] += diagonalTimer.timing();  
    Real percent = 100.0 * INVERTIBLE::inversions() / (_tetMesh->tets().size() * 9);
    static Real maxSeen = percent;
    static Real minSeen = percent;
    if (percent > maxSeen) maxSeen = percent;
    if (percent < minSeen) minSeen = percent;

    /*
    cout << " tet inversions: " << INVERTIBLE::inversions() << " of " << _tetMesh->dofs() * 9 << "(" << percent << "%)" << endl;
    cout << " max seen: " << maxSeen << endl;
    cout << " min seen: " << minSeen << endl;
    */
      
    // get the reduced internal forces
    TIMER internalTimer;
    VECTOR& R = _tetMesh->TET_MESH::generateInternalForces(_Us, _Fhats, _Vs);
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the sparse stiffness matrix
    TIMER stiffnessTimer;
    SPARSE_MATRIX& K = _tetMesh->TET_MESH::generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);
    _timingBreakdown["Stiffness Assembly"] += stiffnessTimer.timing();

    // if this is the first time, set the sparsity of the PETSC matrix
    if (_ASparse.rows() != _rank)
      _ASparse.setSparsity(K);

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    TIMER RHSTimer;
    _residual = _externalForces;
    _residual += R;

    // add collision responses
    _velocity.clear();
    computeCollisionForces();
    computeCollisionForceJacobians();
    _residual -= _collisionForces;
   
    // apply constraints
    applyResidualConstraints();

    maxR = _residual.norm2();
    step++;
    cout << " Newton-Raphson residual " << step << ": " << maxR << endl;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

    /*
    // DEBUG -- breakdown detector
    if (step == 0) oldResidual = maxR;
    Real relativeProgress = fabs((oldResidual - maxR) / oldResidual);
    oldResidual = maxR;
    cout << " Relative progress: " << relativeProgress << endl;
    if (relativeProgress < 1e-2)
      badProgress++;
    else
      badProgress = 0;

    if (badProgress == 10)
    {
      cout << " Solver breakdown!!!!!! " << endl;
      break;
    }
    */

    if (maxR < eps) break;

    // solve with explicit stiffness matrix A
    TIMER ACopy;
    _ASparse.equals(K);
    _timingBreakdown["A Copy"] += ACopy.timing();

    // add force jacobians and constraints
    TIMER addForceTimer;
    _ASparse.axpy(1.0, _collisionForceJacobians);
    applyMatrixConstraints();
    _timingBreakdown["Adding collision forces to matrix"] += addForceTimer.timing();

    TIMER solverTimer;
    Real PCGEps = pow(10.0, -_PCGDigits);
    _ASparse.eps() = PCGEps;

    // DEBUG: force to PCG
    _ASparse.usePCG();
    _ASparse.useICC();
    //bool converged = _ASparse.solveCG(_solution, _residual);
    _ASparse.solveCG(_solution, _residual);
    _timingBreakdown["Solver"] += solverTimer.timing();

    // update positions
    _position += _solution;

    if (step == 0 || maxR < bestResidual)
    {
      bestResidual = maxR;
      bestPosition = _position;
    }
  }

  if (maxR > bestResidual)
  {
    cout << " Rolling solution back to " << bestResidual << endl;
    _position = bestPosition;
  }

  // DEBUG
  static int maxStepsSeen = 0;
  if (step > maxStepsSeen)
    maxStepsSeen = step;
  cout << "Max Newton steps seen: " << maxStepsSeen << endl;
  static int totalStepsSeen = 0;
  totalStepsSeen += step;
  cout << "Mean Newton steps seen: " << (float)totalStepsSeen / (_totalSteps + 1) << endl;
  static Real totalResidual = 0;
  totalResidual += maxR;
  cout << "Mean residual seen: " << totalResidual / (_totalSteps + 1) << endl;

  TIMER finalUpdate;
  // update node positions
  VECTOR& x = _tetMesh->x();
  x.copyInplace(_position);
  _tetMesh->TET_MESH::updateFullMesh();
  _timingBreakdown["Final Update"] += finalUpdate.timing();

  // increment counters
  _time += _dt;
  _totalSteps++;
  _totalTime += total.timing();

  return (step == _maxNewtonSteps);
}

//////////////////////////////////////////////////////////////////////
// build the list of vertices in collision -- cull away those that are a member
// of a colliding triangle but actually outside the surface
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::buildCollisionList(vector<pair<SURFACE*, int> >& collisions,
                                              vector<pair<VEC3F*, SURFACE*> >& collisionVertices,
                                              vector<SURFACE*>& collisionSurfaces)
{
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

    VEC3F v0 = (*(face->vertex(0)));
    VEC3F v1 = (*(face->vertex(1)));
    VEC3F v2 = (*(face->vertex(2)));

    surfaceHash[surface] = true;

    if (surface->inside(v0))
      vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(0), surface)] = true;
    if (surface->inside(v1))
      vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(1), surface)] = true;
    if (surface->inside(v2))
      vertexHash[pair<VEC3F*, SURFACE*>(face->vertex(2), surface)] = true;
  }

  // clear the old collision list
  collisionVertices.clear();

  // build a new one
  for (iter = vertexHash.begin(); iter != vertexHash.end(); iter++)
    collisionVertices.push_back(iter->first);

  // clear the old collision surface list
  collisionSurfaces.clear();

  // build a new one
  map<SURFACE*, bool>::iterator surfaceIter;
  for (surfaceIter = surfaceHash.begin(); surfaceIter != surfaceHash.end(); surfaceIter++)
    collisionSurfaces.push_back(surfaceIter->first);
}

//////////////////////////////////////////////////////////////////////
// Build matrices for an invertible implicit step
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::initializeImplicitStep()
{
  TIMER preamble;

  // accumulate external forces
  computeExternalForces();
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " EXTERNAL FORCES LOADED FROM A FILE " << endl;

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;

  _tetMesh->inertiaTensorOld() = _tetMesh->inertiaTensor();
  _tetMesh->inertiaTensorDtOld() = _tetMesh->inertiaTensorDt();

  _timingBreakdown["Preamble"] += preamble.timing();
}

//////////////////////////////////////////////////////////////////////
// Update the mesh after the partitioned implicit step
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::finalizeImplicitAccelerationStep()
{
  _velocity = _velocityOld + 
              _accelerationAlpha[0] * _accelerationOld +
              _accelerationAlpha[1] * _acceleration;

  _position = _positionOld +
              _accelerationAlpha[2] * _velocityOld +
              _accelerationAlpha[3] * _accelerationOld +
              _accelerationAlpha[4] * _acceleration;

  _tetMesh->x() = _position;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::updateState()
{
	// update velocity
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);

	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], _accelerationOld);
  _acceleration.axpy(-_alpha[1], _velocityOld);
  _acceleration.axpy(-_alpha[0], _positionOld);
  _acceleration.axpy(_alpha[0], _position);

  _tetMesh->x() = _position;
}

//////////////////////////////////////////////////////////////////////
// Update the integrator and mesh state using and acceleration
// level update. Does not update the "olds" however -- this should
// have been done at the beginning of a time step.
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::updateStateUsingAcceleration()
{
  Real* alpha = _accelerationAlpha;

  _position.equals(_positionOld);
  _position.axpy(alpha[2], _velocityOld);
  _position.axpy(alpha[3], _accelerationOld);
  _position.axpy(alpha[4], _acceleration);

  _velocity.equals(_velocityOld); 
  _velocity.axpy(alpha[0], _accelerationOld);
  _velocity.axpy(alpha[1], _acceleration);

  VECTOR& x = _tetMesh->x();
  x.equals(_position);
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an implicit step -- should be called 
// right after initializeImplicitStep
//
// Note that this does not load the mass matrix into A, or add
// the MA term to the residual!
//////////////////////////////////////////////////////////////////////
void FULLSPACE_INTEGRATOR::generateImplicitAccelerationMatrices(map<string,double>& timingBreakdown)
{
  Real* alpha = _accelerationAlpha;
  SPARSE_MATRIX& M = _tetMesh->massMatrix();

  _tetMesh->generateF();

  // precompute diagonalizations --
  // this call is not timed out because timers are inside the function call
  cacheDiagonalizations();

  VECTOR& R = _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);
  const SPARSE_MATRIX& K = _tetMesh->generateSparseStiffnessMatrix(_Us, _Vs, _stiffnesses);
  
  // compute the damping matrix, but only if a new K is available
  _CDampingSparse = _rayleighAlpha * M;
  _CDampingSparse += _rayleighBeta * K;

  //_temp = M * _acceleration;
  _residual = _CDampingSparse * _velocity;

  _residual -= R;
  _residual -= _externalForces;
  // MA is left out because multibody solver adds in MA terms later  
  //_residual += M * _acceleration;
  //timingBreakdown["RHS Assembly"] += RHSTimer.timing();

  // use the acceleration level Newmark consts
  //TIMER ATimer;
  _systemMatrix = alpha[4] * K;
  _systemMatrix += alpha[1] * _CDampingSparse;

  // M is left out because multibody solver adds in MA terms later  
  //_systemMatrix += M;

  _ASparse.setSparsity(_systemMatrix);
  _ASparse.equals(_systemMatrix);
  //timingBreakdown["A Assembly"] += ATimer.timing();
}
