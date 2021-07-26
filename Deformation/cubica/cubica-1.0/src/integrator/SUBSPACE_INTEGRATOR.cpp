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
// SUBSPACE_INTEGRATOR.h: interface for the SUBSPACE_INTEGRATOR class.
//
//////////////////////////////////////////////////////////////////////

#include "SUBSPACE_INTEGRATOR.h"
#include "INVERTIBLE.h"
#ifdef USING_OPENMP
#include <omp.h>
#endif

#define DUMP(x)	" " << #x << "=[ " << x << " ] "

Real SUBSPACE_INTEGRATOR::MAX_EXTERNAL_FORCE = -1.0;

//////////////////////////////////////////////////////////////////////
// Constructor for newmark integrator
//////////////////////////////////////////////////////////////////////
SUBSPACE_INTEGRATOR::SUBSPACE_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, Real dt,
                                         Real alpha, Real beta, Real gravity,
                                         bool inexactNewton,
                                         Real deformationTolerance) :
  _tetMesh(tetMesh),
  _maxNewtonSteps(100),
  _dt(dt),
  _rayleighAlpha(alpha),
  _rayleighBeta(beta),
  _time(0.0),
  _gravityMagnitude(gravity),
  _gravityDown(0,0,-1),
  _clickedNode(NULL),
  _clickedID(-1),
  _forceMultiplier(0.1f),
  _inexactNewton( inexactNewton ),
  _deformationTolerance( deformationTolerance ),
  _totalSolveSteps( 0.0 ),
  _totalSkipped( 0.0 ),
  _totalTime(0.0),
  _solverEps(1e-8),
  _totalSteps(0),
  _useKryslStiffness(false),
  _pulledNode(NULL)
{
  _rank = _tetMesh->rank();
  _position.resizeAndWipe(_rank);
  _velocity.resizeAndWipe(_rank);
  _acceleration.resizeAndWipe(_rank);
  _positionOld.resizeAndWipe(_rank);
  _velocityOld.resizeAndWipe(_rank);
  _accelerationOld.resizeAndWipe(_rank);
  _residual.resizeAndWipe(_rank);
  _totalMeshForces.resizeAndWipe(_rank);
  _temp.resizeAndWipe(_rank);
  _externalForces.resizeAndWipe(_rank);

  _gravity.resizeAndWipe(_rank);

  _CDamping.resizeAndWipe(_tetMesh->rank(), _tetMesh->rank());
  _A.resizeAndWipe(_tetMesh->rank(), _tetMesh->rank());

  // implicit Newmark  
	_beta = 0.25;
	_gamma = 0.50;
	
  // explicit Newmark  
  //_beta = 0;
	//_gamma = 0.50;

  // explicit Euler
  //_beta = 0;
	//_gamma = 0;
 
  // fully implicit Euler
  //_beta = 0.5;
	//_gamma = 1.0;

  // Damped Newmark consts from original HHT (Hilber-Hughes-Taylor) paper
	//_beta = 0.3025;
	//_gamma = 0.6;

  // compute the Newmark consts
  setTimestep(_dt);

  if (_tetMesh->rank() == 0)
    return;

  // get the rest pose damping matrix in case we want to do just
  // linear Rayleigh damping
  //_tetMesh->q().clear();
  _tetMesh->reset();
  _tetMesh->generateF();
  _CDamping.clearingAxpy(_rayleighAlpha, _tetMesh->reducedMass());

  // check that a cubature was in fact loaded by checking the size
  // of the stiffness matrix
  //MATRIX& stiffness = _tetMesh->generateStiffnessMatrix();
  //
  // Only doing this once, so go ahead and do it the Krysl way for
  // maximum accuracy
  MATRIX& stiffness = _tetMesh->generateKryslStiffnessMatrix();
  if (stiffness.rows() > 0)
  {
    _CDamping.axpy(_rayleighBeta, stiffness);
    //_CDamping *= _dt * 0.5f;
  }

  // precompute the gravity vector
  int nodes = _tetMesh->unconstrainedNodes();
  VECTOR fullGravity(nodes * 3);

  VEC3F gravityVector = _gravityDown * _gravityMagnitude;
  for (int x = 0; x < nodes; x++)
  {
    fullGravity(3 * x) = gravityVector[0];
    fullGravity(3 * x + 1) = gravityVector[1];
    fullGravity(3 * x + 2) = gravityVector[2];
  }

  VECTOR massGravity = (_tetMesh->massMatrix()) * fullGravity;
  if (_tetMesh->U().rows() != 0)
    _gravity = _tetMesh->U() ^ massGravity;

  _updateKeyTetDeformation.resize( _tetMesh->totalKeyTets() );

  // compute the body force matrix
  MATRIX I3nx3 = MATRIX::columnOfIdentities(_tetMesh->unconstrainedNodes());

  _bodyForceU = _tetMesh->U() ^ (_tetMesh->massMatrix() * I3nx3);
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::stepImplicit(bool finalize)
{
  TIMER total;
  TIMER preamble;
  // get the reduced mass matrix
  MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // if there isn't a basis yet, do nothing
  if (q.size() == 0 ||
      (!_useKryslStiffness && _tetMesh->totalKeyTets() == 0))
    return;

  // accumulate external forces
  computeExternalForces();
  
  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;
  _tetMesh->qOld() = q;
  _timingBreakdown["Preamble"] += preamble.timing();

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  int step = 0;
  while (step < _maxNewtonSteps && maxR > eps)
  //for (int step = 0; step < 1; step++)
  {
    // an assignment is overly expensive -- the "=" forces an allocation
    TIMER updateMesh;
    q.equals(_position);
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // compute the new deformation gradient F (needed by both R and K)
    TIMER internalTimer;
    _tetMesh->generateF();

    // get the reduced internal forces
    if (_useKryslStiffness)
      _tetMesh->generateKryslInternalForces();
    else
      _tetMesh->generateInternalForces();

    VECTOR& R = _tetMesh->reducedInternalForce();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the reduced stiffness matrix
    TIMER stiffnessTimer;
    if (_useKryslStiffness)
      _tetMesh->generateKryslStiffnessMatrix();
    else
      _tetMesh->generateStiffnessMatrix();
    MATRIX& K = _tetMesh->stiffness();
    _timingBreakdown["Stiffness Assembly"] += stiffnessTimer.timing();

    // compute the damping matrix
    TIMER RHSTimer;
    // DO NOT DO DYNAMIC DAMPING
    _CDamping.clearingAxpy(_rayleighAlpha, M);
    _CDamping.axpy(_rayleighBeta, K);

    // compute the LHS of the residual:
    // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
    _temp.clearingAxpy(-_alpha[2], _accelerationOld);
    _temp.axpy(-_alpha[1], _velocityOld);
    _temp.axpy(-_alpha[0], _positionOld);
    _temp.axpy(_alpha[0], _position);
    _temp = M * _temp;

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
    maxR = _residual.norm2();
    step++;
    //cout << " Subspace Newton-Raphson residual " << step << ": " << maxR << endl;
   
    // assemble system matrix A
    TIMER ATimer;
    _A = K;
    _A.axpy(_alpha[0], M);
    _A.axpy(_alpha[3], _CDamping);
    _timingBreakdown["A Assembly"] += ATimer.timing();

    TIMER solverTimer;
    // solve with LU factorization
    //_A.solve(_residual);
    //_position -= _residual;
    bool success = _A.factorLU();
    if (!success)
    {
      cout << __FILE__ << " " << __LINE__ << " : " << endl;
      _A = K;
      _A.axpy(_alpha[0], M);
      _A.axpy(_alpha[3], _CDamping);
      cout << " Bad A: " << _A << endl;
      cout << " key weights: " << _tetMesh->keyWeights() << endl;
      cout << " K: " << K << endl;
      cout << " M: " << M << endl;
      cout << " CDamping: " << _CDamping << endl;
      cout << " F: " << _tetMesh->F() << endl;
      cout << " left H: " << _tetMesh->leftH() << endl;
      cout << " force density: " << _tetMesh->forceDensity() << endl;
      cout << " key inverted? " << _tetMesh->keyInverted() << endl;

      // factor it again so the solveLU bombs too
      _A.factorLU();
    }
    _A.solveLU(_residual);
    _position -= _residual;

    /*
    // solve with CG
    _temp.clear();
    _A.maxIterations() = 10000;
    _A.eps() = _solverEps;
    //_A.solveCG(_temp, _residual);
    _A.solveICCG(_temp, _residual);
    _position -= _temp;
    */
    _timingBreakdown["Solver"] += solverTimer.timing();

    //if (_tetMesh->keyInverted())
    //  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : Key tet inverted! " << endl;
  }
  //cout << " Subspace Newton solver steps: " << step << endl;

	// NEW NEW NEW update velocity
  TIMER finalUpdate;
  q.equals(_position);

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

  // update node positions
  if (finalize)
  {
    // pushing this to main
    //_tetMesh->updateSurfaceMesh();
    _externalForces.clear();
  }
  _timingBreakdown["Final Update"] += finalUpdate.timing();
  _totalSteps++;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::stepInvertibleImplicit(bool finalize, bool debugPrint)
{
  cout << "=========================================" << endl;
  cout << " Dynamics step " << _totalSteps << endl;
  cout << "=========================================" << endl;

  TIMER total;
  TIMER preamble;
  // get the reduced mass matrix
  MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // if there isn't a basis yet, do nothing
  if (q.size() == 0 || 
      (!_useKryslStiffness && _tetMesh->totalKeyTets() == 0))
    return;

  // accumulate external forces
  computeExternalForces();

  // copy the current values into the old values
#if 0
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;
  _tetMesh->qOld() = q;
#endif
  _positionOld.copyInplace( _position, true /* resize if needed */ );
  _velocityOld.copyInplace( _velocity, true /* resize if needed */ );
  _accelerationOld.copyInplace( _acceleration, true /* resize if needed */ );
  _tetMesh->qOld().copyInplace( q, true /* resize if needed */ );
  _timingBreakdown["Preamble"] += preamble.timing();

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  int step = 0;
  if ( debugPrint ) 
    cout << "_maxNewtonSteps = " << _maxNewtonSteps << endl;
  while (step < _maxNewtonSteps && maxR > eps)
  //while (step < 100 && maxR > eps)
  //for (int step = 0; step < 3; step++)
  {
    // update tet mesh with the new q vector
    // an assignment is overly expensive -- the "=" forces an allocation
    TIMER updateMesh;
    q.equals(_position);
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // compute the new deformation gradient F (needed by both R and K)
    if (_useKryslStiffness)
    {
      _tetMesh->updateFullMesh();
      _tetMesh->TET_MESH::generateF();
    }
    else  
      _tetMesh->SUBSPACE_TET_MESH::generateF();

    // precompute diagonalizations --
    // this call is not timed out because timers are inside the function call.
    // Indicate to cacheDiagonalizations whether or not we want it to
    // possibly ignore certain mildly deformed tets (ie. if we are forming
    // an inexact stiffness matrix).
    if (_useKryslStiffness)
      cacheKryslDiagonalizations();
    else
      cacheDiagonalizations( _inexactNewton );

    // get the reduced internal forces
    TIMER internalTimer;
    if (_useKryslStiffness)
      _tetMesh->generateKryslInternalForces(_Us, _Fhats, _Vs);
    else
      _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);

    VECTOR& R = _tetMesh->reducedInternalForce();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the reduced stiffness matrix -- do this even if we are 
    // doing the conjugate gradient solve, because the damping matrix
    // needs it
    //TIMER stiffnessTimer;
    if (_useKryslStiffness)
    {
      _tetMesh->generateKryslStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
    }
    else
    {
      //_tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);

      // If we are doing inexact newton solves, we just need to update the
      // reduced stiffness matrix using the appropriate key tets
      if ( _inexactNewton )
      {
        _tetMesh->updateStiffnessMatrixFast(_Us, _Vs, _stiffnesses,
                                            _updateKeyTetDeformation,
                                            _timingBreakdown);

        if ( debugPrint )
        {
          int numSkipped = 0;
          _totalSolveSteps += 1.0;

          for (unsigned int i = 0; i < _updateKeyTetDeformation.size(); i++ )
          {
            if ( !_updateKeyTetDeformation[i] )
              numSkipped++;
          }
          _totalSkipped += (Real)numSkipped;

          cout << "Ignoring " << numSkipped << " key tets" << endl;
        }
      }
      else
      {
        _tetMesh->generateStiffnessMatrixFast(_Us, _Vs, _stiffnesses,
                                              _timingBreakdown);
      }
    }
    MATRIX& K = _tetMesh->stiffness();
    //_timingBreakdown["Stiffness Assembly"] += stiffnessTimer.timing();
    
    // compute the damping matrix, but only if a new K is available
    TIMER RHSTimer;
    // DO NOT DO DYNAMIC DAMPING
    //_CDamping.clearingAxpy(_rayleighAlpha, M);
    //_CDamping.axpy(_rayleighBeta, K);

    // compute the LHS of the residual:
    // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
    _temp.clearingAxpy(-_alpha[2], _accelerationOld);
    _temp.axpy(-_alpha[1], _velocityOld);
    _temp.axpy(-_alpha[0], _positionOld);
    _temp.axpy(_alpha[0], _position);
    _temp = M * _temp;

    VECTOR debugMa = _temp;

    if ( debugPrint )
    {
      cout << DUMP( step ) << endl;
      cout << DUMP( _temp.norm2() ) << endl;
      cout << DUMP( _accelerationOld.norm2() ) << endl;
      cout << DUMP( _velocityOld.norm2() ) << endl;
      cout << DUMP( _positionOld.norm2() ) << endl;
      cout << DUMP( _position.norm2() ) << endl;
    }

    // compute the RHS of the residual:
    // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
    _residual.clearingAxpy(_alpha[5], _accelerationOld);
    _residual.axpy(_alpha[4], _velocityOld);
    _residual.axpy(-_alpha[3], _positionOld);
    _residual.axpy(_alpha[3], _position);
    _residual = _CDamping * _residual;

#if 0
    if ( debugPrint )
    {
      cout << "After damping..." << endl;
      cout << DUMP( _residual.norm2() ) << endl;
      cout << DUMP( R.norm2() ) << endl;
      cout << DUMP( _externalForces.norm2() ) << endl;
      cout << endl << endl;
    }
#endif

    VECTOR debugCv= _residual;

    // assemble full residual: LHS + RHS + R - F
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;
    
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();
    maxR = _residual.norm2();
    step++;
   
    // assemble system matrix A
    TIMER ATimer;
    _A = K;
    _A.axpy(_alpha[0], M);
    _A.axpy(_alpha[3], _CDamping);
    _timingBreakdown["A Assembly"] += ATimer.timing();

    //cout << " temp: " << _temp.norm2() << endl;
    //cout << " R: " << R.norm2() << endl;
    //cout << " external: " << _externalForces.norm2() << endl;
    //cout << " Residual: " << _residual.norm2() << endl;

    TIMER solverTimer;
    
    // solve with LU factorization
    //bool success = _A.factorLU();
    bool success = _A.factorCholesky();
    if (!success)
    {
      _A = K;
      _A.axpy(_alpha[0], M);
      _A.axpy(_alpha[3], _CDamping);
      cout << " Bad A: " << _A << endl;
      cout << " key weights: " << _tetMesh->keyWeights() << endl;
      cout << " K: " << K << endl;
      cout << " M: " << M << endl;
      cout << " CDamping: " << _CDamping << endl;
      cout << " F: " << _tetMesh->F() << endl;
      cout << " left H: " << _tetMesh->leftH() << endl;
      cout << " force density: " << _tetMesh->forceDensity() << endl;
      cout << " key inverted? " << _tetMesh->keyInverted() << endl;

      // factor it again so the solveLU bombs too
      _A.factorLU();
    }
    //_A.solveLU(_residual);
    _A.solveCholesky(_residual);
    _position -= _residual;

    /*
    if (directSolve)
    {
      // solve with LU
      //_A.solve(_residual);
      //_position -= _residual;

      _A.factorCholesky();
      _A.solveCholesky(_residual);
      _position -= _residual;
    }
    else
    {
      // solve with CG
      _temp.clear();
      solveInvertibleStiffness(_temp, _residual, false);
      _position -= _temp;
    }
    */

    _timingBreakdown["Solver"] += solverTimer.timing();
  }

  if ( debugPrint )
  {
    cout << " Subspace Newton solver steps: " << step << endl;
    cout << " Solver eps: " << _solverEps << endl;
    cout << " Simulation dt: " << _dt << endl;
    cout << " \t\t\t\t\t\t\tMaxR: " << maxR << endl;
    cout << " epx: " << eps << endl;
    cout << " Skipping an average of " << ( _totalSkipped / _totalSolveSteps ) << endl;
#if 0
    cout << DUMP( ( _position - _positionOld ).norm2() ) << endl;
#endif
    if ( maxR > 1.0 )
    {
      cout << "Large residual: " << DUMP( maxR ) << endl;
    }
  }

	// update velocity
  TIMER finalUpdatePosition;
  q.equals(_position);
  _timingBreakdown["Final Update Position"] += finalUpdatePosition.timing();
  
  TIMER finalUpdateVelocity;
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);
  _timingBreakdown["Final Update Velocity"] += finalUpdateVelocity.timing();

	// update acceleration
  TIMER finalUpdateAcceleration;
  _acceleration.clearingAxpy(-_alpha[2], _accelerationOld);
  _acceleration.axpy(-_alpha[1], _velocityOld);
  _acceleration.axpy(-_alpha[0], _positionOld);
  _acceleration.axpy(_alpha[0], _position);
  _timingBreakdown["Final Update Acceleration"] += finalUpdateAcceleration.timing();

  //cout << " position: " << _position << endl;
  //cout << " velocity: " << _velocity << endl;
  //cout << " acceleration: " << _acceleration << endl;

#if 0
  if ( debugPrint )
  {
#endif
    Real diff = ( _acceleration - _accelerationOld ).norm2() / ( _accelerationOld.norm2() );

    if ( diff > 5.0 )
    {
      cout << "Large acceleration increase!" << endl;
      cout << DUMP( _accelerationOld.norm2() ) << endl;
      cout << DUMP( _velocityOld.norm2() ) << endl;
      cout << DUMP( _positionOld.norm2() ) << endl;
      cout << DUMP( _position.norm2() ) << endl;
      cout << DUMP( ( _position - _positionOld ).norm2() ) << endl;
      cout << DUMP( _acceleration.norm2() ) << endl;
      cout << endl;
    }
#if 0
  }
#endif

  // update node positions
  TIMER finalUpdateForces;
  if (finalize)
  {
    // pushing to main
    //_tetMesh->updateSurfaceMesh();
    _externalForces.clear();
  }
  _timingBreakdown["Final Update Forces"] += finalUpdateForces.timing();
  _totalSteps++;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration, using acceleration as the primary variable
//////////////////////////////////////////////////////////////////////
int SUBSPACE_INTEGRATOR::stepImplicitAcceleration()
{
  TIMER total;
  TIMER preamble;
  // get the reduced mass matrix
  MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // if there isn't a basis yet, do nothing
  if (q.size() == 0 || 
      (!_useKryslStiffness && _tetMesh->totalKeyTets() == 0))
    return 0;

  // accumulate external forces
  computeExternalForces();
  
  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld =  _velocity;
  _accelerationOld = _acceleration;
  _tetMesh->qOld() = q;
  _timingBreakdown["Preamble"] += preamble.timing();

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  Real initialR = maxR;
  int step = 0;
  Real* alpha = _accelerationAlpha;
 
  cout << "=========================" << endl;
  cout << " Timestep " << _totalSteps << endl;
  cout << "=========================" << endl;
  
  //while (step < _maxNewtonSteps && maxR > eps)
  while (step < _maxNewtonSteps && maxR > initialR * eps && maxR > eps)
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
    q.equals(_position);
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // compute the new deformation gradient F (needed by both R and K)
    _tetMesh->SUBSPACE_TET_MESH::generateF();

    // precompute diagonalizations --
    // this call is not timed out because timers are inside the function call
    cacheDiagonalizations();

    // get the reduced internal forces
    TIMER internalTimer;
    _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);

    VECTOR& R = _tetMesh->reducedInternalForce();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the reduced stiffness matrix -- do this even if we are 
    // doing the conjugate gradient solve, because the damping matrix
    // needs it
    _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
    MATRIX& K = _tetMesh->stiffness();
    
    // compute the damping matrix, but only if a new K is available
    TIMER RHSTimer;
    _CDamping.clearingAxpy(_rayleighAlpha, M);
    _CDamping.axpy(_rayleighBeta, K);

    _temp = M * _acceleration;
    _residual = _CDamping * _velocity;
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();
    maxR = _residual.norm2();
    if (step == 0) 
      initialR = maxR; 
    step++;

    //cout << " residual: " << _residual << endl;

    // use the acceleration level Newmark consts
    TIMER ATimer;
    _A = M;
    _A.axpy(alpha[4], K);
    _A.axpy(alpha[1], _CDamping);
    _timingBreakdown["A Assembly"] += ATimer.timing();

    TIMER solverTimer;
    
    // solve with LU factorization
    //bool success = _A.factorLU();
    bool success = _A.factorCholesky();
    if (!success)
    {
      _A = K;
      _A.axpy(_alpha[0], M);
      _A.axpy(_alpha[3], _CDamping);
      cout << " Bad A: " << _A << endl;
      cout << " key weights: " << _tetMesh->keyWeights() << endl;
      cout << " K: " << K << endl;
      cout << " M: " << M << endl;
      cout << " CDamping: " << _CDamping << endl;
      cout << " F: " << _tetMesh->F() << endl;
      cout << " left H: " << _tetMesh->leftH() << endl;
      cout << " force density: " << _tetMesh->forceDensity() << endl;
      cout << " key inverted? " << _tetMesh->keyInverted() << endl;

      // factor it again so the solveLU bombs too
      _A.factorLU();
      exit(0);
    }
    cout << " residual " << step << ": " << maxR << endl;
    _A.solveCholesky(_residual);

    // update acceleration according to solution
    _acceleration -= _residual;

    _timingBreakdown["Solver"] += solverTimer.timing();
  }

	// update velocity
  TIMER finalUpdatePosition;
  _tetMesh->q() = _position;
  _timingBreakdown["Final Update Position"] += finalUpdatePosition.timing();
  
  // update node positions
  TIMER finalUpdateForces;
  // pushing to main
  //_tetMesh->updateSurfaceMesh();
  _externalForces.clear();

  _timingBreakdown["Final Update Forces"] += finalUpdateForces.timing();
  _totalSteps++;
  _totalTime += total.timing();

  return step;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration, using acceleration as the primary variable,
// and taking into account accelerations, with collisions
// handled implicitly
//////////////////////////////////////////////////////////////////////
int SUBSPACE_INTEGRATOR::stepImplicitWithImplicitCollisions(vector<pair<SURFACE*, int> >& collisions)
{
  TIMER total;
  TIMER preamble;
  // get the reduced mass matrix
  MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // if there isn't a basis yet, do nothing
  if (q.size() == 0 || 
      (!_useKryslStiffness && _tetMesh->totalKeyTets() == 0))
    return 0;

  // accumulate external forces
  computeExternalForces();
 
  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld =  _velocity;
  _accelerationOld = _acceleration;
  _tetMesh->qOld() = q;
  _timingBreakdown["Preamble"] += preamble.timing();

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  Real initialR = maxR;
  int step = 0;
  Real* alpha = _accelerationAlpha;
 
  cout << "=========================" << endl;
  cout << " Timestep " << _totalSteps << endl;
  cout << "=========================" << endl;
  
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
    q.equals(_position);
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // compute the new deformation gradient F (needed by both R and K)
    TIMER timerF;
    _tetMesh->SUBSPACE_TET_MESH::generateF();
    _timingBreakdown["Generate F"] += timerF.timing();

    // precompute diagonalizations --
    // this call is not timed out because timers are inside the function call
    cacheDiagonalizations();

    // get the reduced internal forces
    TIMER internalTimer;
    _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);

    VECTOR& R = _tetMesh->reducedInternalForce();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the reduced stiffness matrix -- do this even if we are 
    // doing the conjugate gradient solve, because the damping matrix
    // needs it
    _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
    MATRIX& K = _tetMesh->stiffness();
    
    // compute the damping matrix, but only if a new K is available
    TIMER RHSTimer;
    _CDamping.clearingAxpy(_rayleighAlpha, M);
    _CDamping.axpy(_rayleighBeta, K);

    _temp = M * _acceleration;
    _residual = _CDamping * _velocity;
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;

    // add collision forces to the residual
    addImplicitCollisionResiduals(collisions);
    
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();
    maxR = _residual.norm2();
    if (step == 0) 
      initialR = maxR; 
    step++;

    // use the acceleration level Newmark consts
    TIMER ATimer;
    _A = M;
    _A.axpy(alpha[4], K);
    _A.axpy(alpha[1], _CDamping);

    // add collision forces to the Jacobian
    addImplicitCollisionJacobians(collisions);

    _timingBreakdown["A Assembly"] += ATimer.timing();

    TIMER solverTimer;
    
    // solve with LU factorization
    bool success = _A.factorLU();
    //bool success = _A.factorCholesky();
    if (!success)
    {
      _A = K;
      _A.axpy(_alpha[0], M);
      _A.axpy(_alpha[3], _CDamping);
      /*
      cout << " Bad A: " << _A << endl;
      cout << " key weights: " << _tetMesh->keyWeights() << endl;
      cout << " K: " << K << endl;
      cout << " M: " << M << endl;
      cout << " CDamping: " << _CDamping << endl;
      cout << " F: " << _tetMesh->F() << endl;
      cout << " left H: " << _tetMesh->leftH() << endl;
      cout << " force density: " << _tetMesh->forceDensity() << endl;
      cout << " key inverted? " << _tetMesh->keyInverted() << endl;

      // factor it again so the solveLU bombs too
      _A.factorLU();
      */
      exit(0);
    }
    cout << " residual " << step << ": " << maxR << endl;
    //_A.solveCholesky(_residual);
    _A.solveLU(_residual);

    // update acceleration according to solution
    _acceleration -= _residual;

    _timingBreakdown["Solver"] += solverTimer.timing();
  }

	// update velocity
  TIMER finalUpdatePosition;
  _tetMesh->q() = _position;
  _timingBreakdown["Final Update Position"] += finalUpdatePosition.timing();
  
  // update node positions
  TIMER finalUpdateForces;
  // pushing to main
  //_tetMesh->updateSurfaceMesh();
  _externalForces.clear();

  _timingBreakdown["Final Update Forces"] += finalUpdateForces.timing();
  _totalSteps++;
  _totalTime += total.timing();

  return step;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration, using acceleration as the primary variable,
// and taking into account accelerations, with collisions
// handled explicitly
//////////////////////////////////////////////////////////////////////
int SUBSPACE_INTEGRATOR::stepImplicitWithExplicitCollisions(vector<pair<SURFACE*, int> >& collisions)
{
  TIMER total;
  TIMER preamble;
  // get the reduced mass matrix
  MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // if there isn't a basis yet, do nothing
  if (q.size() == 0 || 
      (!_useKryslStiffness && _tetMesh->totalKeyTets() == 0))
    return 0;

  // accumulate external forces
  computeExternalForces();
 
  // add explicit collision forces
  addExplicitCollisions(collisions);

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld =  _velocity;
  _accelerationOld = _acceleration;
  _tetMesh->qOld() = q;
  _timingBreakdown["Preamble"] += preamble.timing();

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  Real initialR = maxR;
  int step = 0;
  Real* alpha = _accelerationAlpha;
 
  cout << "=========================" << endl;
  cout << " Timestep " << _totalSteps << endl;
  cout << "=========================" << endl;
  
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
    q.equals(_position);
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // compute the new deformation gradient F (needed by both R and K)
    _tetMesh->SUBSPACE_TET_MESH::generateF();

    // precompute diagonalizations --
    // this call is not timed out because timers are inside the function call
    cacheDiagonalizations();

    // get the reduced internal forces
    TIMER internalTimer;
    _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);

    VECTOR& R = _tetMesh->reducedInternalForce();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the reduced stiffness matrix -- do this even if we are 
    // doing the conjugate gradient solve, because the damping matrix
    // needs it
    _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
    MATRIX& K = _tetMesh->stiffness();
    
    // compute the damping matrix, but only if a new K is available
    TIMER RHSTimer;
    _CDamping.clearingAxpy(_rayleighAlpha, M);
    _CDamping.axpy(_rayleighBeta, K);

    _temp = M * _acceleration;
    _residual = _CDamping * _velocity;
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();
    maxR = _residual.norm2();
    if (step == 0) 
      initialR = maxR; 
    step++;

    //cout << " residual: " << _residual << endl;

    // use the acceleration level Newmark consts
    TIMER ATimer;
    _A = M;
    _A.axpy(alpha[4], K);
    _A.axpy(alpha[1], _CDamping);
    _timingBreakdown["A Assembly"] += ATimer.timing();

    TIMER solverTimer;
    
    // solve with LU factorization
    //bool success = _A.factorLU();
    bool success = _A.factorCholesky();
    if (!success)
    {
      _A = K;
      _A.axpy(_alpha[0], M);
      _A.axpy(_alpha[3], _CDamping);
      cout << " Bad A: " << _A << endl;
      cout << " key weights: " << _tetMesh->keyWeights() << endl;
      cout << " K: " << K << endl;
      cout << " M: " << M << endl;
      cout << " CDamping: " << _CDamping << endl;
      cout << " F: " << _tetMesh->F() << endl;
      cout << " left H: " << _tetMesh->leftH() << endl;
      cout << " force density: " << _tetMesh->forceDensity() << endl;
      cout << " key inverted? " << _tetMesh->keyInverted() << endl;

      // factor it again so the solveLU bombs too
      _A.factorLU();
      exit(0);
    }
    cout << " residual " << step << ": " << maxR << endl;
    _A.solveCholesky(_residual);

    // update acceleration according to solution
    _acceleration -= _residual;

    _timingBreakdown["Solver"] += solverTimer.timing();
  }

	// update velocity
  TIMER finalUpdatePosition;
  _tetMesh->q() = _position;
  _timingBreakdown["Final Update Position"] += finalUpdatePosition.timing();
  
  // update node positions
  TIMER finalUpdateForces;
  // pushing to main
  //_tetMesh->updateSurfaceMesh();
  _externalForces.clear();

  _timingBreakdown["Final Update Forces"] += finalUpdateForces.timing();
  _totalSteps++;
  _totalTime += total.timing();

  return step;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::stepInvertibleSemiImplicit(bool directSolve)
{
  TIMER total;
  TIMER preamble;
  // get the reduced mass matrix
  MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // accumulate external forces
  computeExternalForces();

  bool verbose = (_externalForces.norm2() > 1e-2) ? true : false;
  
  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;
  _timingBreakdown["Preamble"] += preamble.timing();

  if (!directSolve)
  {
    // set CG solve parameters
    _A.maxIterations() = 10000;
    _A.eps() = 1e-7;
    if (_A.rows() != _rank)
      _A.resizeAndWipe(_rank, _rank);
  }
  
  // do one Newton-Raphson step
  
  // update tet mesh with the new q vector
  // an assignment is overly expensive -- the "=" forces an allocation
  TIMER updateMesh;
  q.equals(_position);
  _timingBreakdown["Update Mesh"] += updateMesh.timing();
  
  // compute the new deformation gradient F (needed by both R and K)
  _tetMesh->generateF();

  // precompute diagonalizations --
  // this call is not timed out because timers are inside the function call
  cacheDiagonalizations();

  // get the reduced internal forces
  TIMER internalTimer;
  VECTOR& R = _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);
  _timingBreakdown["Internal Forces"] += internalTimer.timing();

  // get the reduced stiffness matrix -- do this even if we are 
  // doing the conjugate gradient solve, because the damping matrix
  // needs it
  //TIMER stiffnessTimer;
  MATRIX& K = _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
  
  // compute the damping matrix, but only if a new K is available
  TIMER RHSTimer;
  _CDamping.clearingAxpy(_rayleighAlpha, M);
  _CDamping.axpy(_rayleighBeta, K);

  // compute the LHS of the residual:
  // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
  _temp.clearingAxpy(-_alpha[2], _accelerationOld);
  _temp.axpy(-_alpha[1], _velocityOld);
  _temp.axpy(-_alpha[0], _positionOld);
  _temp.axpy(_alpha[0], _position);
  _temp = M * _temp;

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
  if (directSolve)
    _A = K;
  else
    _A.clear();
  _A.axpy(_alpha[0], M);
  _A.axpy(_alpha[3], _CDamping);
  _timingBreakdown["A Assembly"] += ATimer.timing();

  if (verbose)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "A: " << _A << endl;
    cout << "b: " << _residual << endl;
  }

  TIMER solverTimer;
  if (directSolve)
  {
    // solve with LU
    //_A.solve(_residual);
    //_position -= _residual;

    _A.factorCholesky();
    _A.solveCholesky(_residual);
    _position -= _residual;
  }
  else
  {
    // solve with CG
    _temp.clear();
    solveInvertibleStiffness(_temp, _residual, false);
    _position -= _temp;
  }
  _timingBreakdown["Solver"] += solverTimer.timing();

	// update velocity
  TIMER finalUpdate;
  q.equals(_position);
  
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);

	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], _accelerationOld);
  _acceleration.axpy(-_alpha[1], _velocityOld);
  _acceleration.axpy(-_alpha[0], _positionOld);
  _acceleration.axpy(_alpha[0], _position);
  
  if (verbose)
  {
    cout << "solution: " << _residual << endl;
    cout << "velocity: " << _velocity << endl;
    cout << "acceleration: " << _acceleration << endl;
  }

  // update node positions
  _tetMesh->updateSurfaceMesh();
  _externalForces.clear();
  _timingBreakdown["Final Update"] += finalUpdate.timing();
  _totalSteps++;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// Quasistatic solver
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::stepInvertibleQuasistatic(bool directSolve)
{
  cout << "=====================================" << endl;
  cout << " Subspace quasistatic step " << _totalSteps << endl;
  cout << "=====================================" << endl;

  TIMER total;
  TIMER preamble;

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // if there isn't a basis yet, do nothing
  if (q.size() == 0 || 
      (!_useKryslStiffness && _tetMesh->totalKeyTets() == 0))
    return;

  // accumulate external forces
  computeExternalForces();
  
  // copy the current values into the old values
  _positionOld = _position;

  _timingBreakdown["Preamble"] += preamble.timing();

  /*
  _position[0] = 0.00532371;
  _position[1] = -0.00582856;
  _position[2] = -0.00226783;
  _position[3] = 0.00156537;
  _position[4] = -0.00735499;
  _position[5] = 0.0165258;
  _position[6] = 0.00248066;
  _position[7] = -0.00334456;
  _position[8] = 0.0130519;
  _position[9] = 0.0158753;
  */

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = 10.0 * eps;
  int iterations = 0;
  //for (int step = 0; step < 3; step++)
  int step = 0;
  while (step < _maxNewtonSteps && maxR > eps)
  {
    TIMER updateMesh;
    // update tet mesh with the new q vector
    // an assignment is overly expensive -- the "=" forces an allocation
    q.equals(_position);
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
   
    // compute the new deformation gradient F (needed by both R and K)
    if (_useKryslStiffness)
    {
      _tetMesh->updateFullMesh();
      _tetMesh->TET_MESH::generateF();
    }
    else  
      _tetMesh->SUBSPACE_TET_MESH::generateF();

    // precompute diagonalizations --
    // this call is not timed out because timers are inside the function call
    if (_useKryslStiffness)
      cacheKryslDiagonalizations();
    else
      cacheDiagonalizations();
    
    TIMER internalTimer;
    // get the reduced internal forces
    if (_useKryslStiffness)
      _tetMesh->generateKryslInternalForces(_Us, _Fhats, _Vs);
    else
      _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);
    VECTOR& R = _tetMesh->reducedInternalForce();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // only do stiffness assembly if we're doing a direct solve
    if (directSolve)
    {
      if (_useKryslStiffness)
        _tetMesh->generateKryslStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
      else
        _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
      _A = _tetMesh->stiffness();
    }
    
    TIMER RHSTimer;
    _residual.equals(_externalForces);
    _residual += R;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();
    maxR = _residual.norm2();
    //cout << " Energy: " << (0.5 * (_residual ^ _residual)) << endl;
    //cout << " Norm: " << _residual.norm2() << endl;

    cout << " Subspace Newton-Raphson residual " << step << ": " << maxR << endl;
    if (maxR < eps) continue;
    iterations++;

    TIMER solverTimer;
    if (directSolve)
    {
      // solve with Cholesky
      _A.factorCholesky();
      _A.solveCholesky(_residual);
      _position += _residual;

      //cout << " solution norm: " << _residual.norm2() << endl;
    }
    else
    {
      // solve with conjugate gradient
      _temp.clear();
      solveInvertibleStiffness(_temp, _residual, true);
      _position += _temp;
    }
    _timingBreakdown["Solver"] += solverTimer.timing();
    step++;
  }

  //cout << " Subspace Newton-Raphson iterations: " << iterations << " step: " << _totalSteps << endl;
  TIMER finalUpdate;
  q.equals(_position);

  // update node positions
  _tetMesh->updateSurfaceMesh();
  _externalForces.clear();
  _timingBreakdown["Final Update"] += finalUpdate.timing();
  _totalSteps++;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// Quasistatic solver
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::stepQuasistatic()
{
  TIMER total;
  TIMER preamble;

  // accumulate external forces
  computeExternalForces();
  
  // copy the current values into the old values
  _positionOld = _position;

  VECTOR& q = _tetMesh->q();
  _timingBreakdown["Preamble"] += preamble.timing();

  // do Newton-Raphson
  Real eps = 1e-15;
  Real rNorm = 10.0 * eps;
  int iterations = 0;
  //for (int step = 0; step < 3; step++)
  while (rNorm > eps)
  {
    TIMER updateMesh;
    // update tet mesh with the new q vector
    // an assignment is overly expensive -- the "=" forces an allocation
    q.equals(_position);
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    TIMER internalTimer;

    // compute the new deformation gradient F (needed by both R and K)
    _tetMesh->generateF();
    
    // get the reduced internal forces
    VECTOR& R = _tetMesh->generateInternalForces();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    TIMER stiffnessAssembly;
    MATRIX& K = _tetMesh->generateStiffnessMatrix();
    _A = K;
    _timingBreakdown["Stiffness Assembly"] += stiffnessAssembly.timing();
    
    TIMER RHSTimer;
    _residual.equals(_externalForces);
    _residual += R;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();
    rNorm = _residual.norm2();
    if (rNorm < eps) continue;
    iterations++;
  
    TIMER solverTimer;
    // solve with LU factorization
    _A.solve(_residual);
    _position += _residual;
    _timingBreakdown["Solver"] += solverTimer.timing();
  }
  cout << " Newton-Raphson iterations: " << iterations << endl;
  TIMER finalUpdate;
  q.equals(_position);

  // update node positions
  _tetMesh->updateSurfaceMesh();
  _externalForces.clear();
  _timingBreakdown["Final Update"] += finalUpdate.timing();
  _totalSteps++;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// Explicit Integration
//
// This follows pages 324-325 of: 
// Computational Contact Mechanics by Peter Wriggers, 2006
//
// The equation to solve is:
// (M + \frac{\delta t}{2} C) u_{n+1} = 
//    (\delta t)^2 * [P_n - R * u_n] + 
//    \frac{\delta t}{2} * C * u_{n-1} + 
//    M* (2 * u_n - u_{n-1})
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::stepExplicit()
{
  // explicit does not use the acceleration vectors,
  // so use them as scratch space
  VECTOR& scratch = _acceleration;
  
  // update tet mesh with the new q vector
  VECTOR& q = _tetMesh->q();
  q.equals(_position);

  // get the forces
  _tetMesh->generateF();
  computeExternalForces();
  VECTOR& R = _externalForces;
  R += _tetMesh->generateInternalForces();

  // start accum'ing in _temp,
  // switch to 'scratch' to avoid a dynamic allocation
  MATRIX& M = _tetMesh->reducedMass();
  _temp.clearingAxpy(2.0, _position);
  _temp.axpy(-1.0, _positionOld);
  M.multiplyInplace(_temp, scratch);
  scratch.axpy(_dt * _dt, R);

  /*
  // use non-linear Rayleigh damping
  // This breaks invertible finite elements
  MATRIX& K = _tetMesh->generateStiffnessMatrix();
  _CDamping.clearingAxpy(_rayleighAlpha, M);
  _CDamping.axpy(_rayleighBeta, K);
  _CDamping *= _dt * 0.5f;
  */
  
  // add damping
  VECTOR damping = _CDamping * _positionOld;
  scratch.axpy(1.0, damping);
  _A = _tetMesh->reducedMass();
  _A += _CDamping;
  _A.solve(scratch);

  _positionOld.swap(_position);
  _position.swap(scratch);

  _velocity = _position - _positionOld;
  _velocity *= 1.0 / _dt;

  // update node positions
  _tetMesh->updateSurfaceMesh();
  _externalForces.clear();

  // increment simulation time
  _time += _dt;
}

//////////////////////////////////////////////////////////////////////
// process a mouse click event
//////////////////////////////////////////////////////////////////////
bool SUBSPACE_INTEGRATOR::click(VEC3F& point, Real maxRadius)
{
  // find and store the nearest surface point
  _clickedNode = _tetMesh->closestSurfaceNode(point);

  Real clickDistance = norm2( point - *_clickedNode );
  if ( maxRadius > 0 && clickDistance > maxRadius * maxRadius )
  {
    // The point is too far away, negate the interaction
    _clickedNode = NULL;
    return false;
  }

  _clickedID = _tetMesh->vertexID( _clickedNode );
  _clickedPosition = *_clickedNode;
  _draggedPosition = *_clickedNode;

  cout << " clicked ID: " << _clickedID << endl;

  return true;
}

//////////////////////////////////////////////////////////////////////
// Compute all external forces
// FIXME: Lots of dynamic allocations here!  Very slow probably!
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::computeExternalForces()
{
  // process mouse input
  if (_clickedNode != NULL)
  {
    VEC3F dragForce = _draggedPosition - _clickedPosition;

#if 0
    // Correct the force magnitude, if necessary
    Real forceMagnitude = norm( dragForce * _forceMultiplier );
    if ( MAX_EXTERNAL_FORCE > 0.0 && forceMagnitude > MAX_EXTERNAL_FORCE )
    {
      // Clamp it!
      Real forceRatio = forceMagnitude / MAX_EXTERNAL_FORCE;
      dragForce /= forceRatio;
    }
#endif

    //dragForce.normalize();
    _forceVectors.push_back(dragForce);
    _forceNodes.push_back(_clickedNode);
    _forceNodeIDs.push_back(_clickedID);
  }

  // make sure some forces are in the list
  if (_forceVectors.size() == 0) return;

#if 0
  // clear the force vector on the first sample
  // FIXME: Repeatedly looking up node by its vertex pointer
  // is probably slow.  We should just get the index once and
  // be done with it.
  MATRIX subbasis = _tetMesh->vertexSubBasis(_forceNodes[0]);
#endif

  // TODO - try this out
  if ( _forceNodes.size() == 0 || _forceNodeIDs.size() == 0)
    return;

  MATRIX &subbasis = _subBasisWorkspace.get();

  if ( subbasis.rows() != 3 || subbasis.cols() != _rank )
  {
    subbasis.resizeAndWipe( 3, _rank );
  }

  //_tetMesh->vertexSubBasis( _forceNodeIDs[0], subbasis, true /* transpose */ );
  _tetMesh->vertexSubBasis( _forceNodeIDs[0], subbasis );
  // TODO - end

#if 0
  // if there is no basis yet, do nothing
  if (subbasis.cols() == 0 || subbasis.rows() == 0)
    return;

  subbasis = subbasis.transpose();
  _externalForces += subbasis.gemv(_forceMultiplier, _forceVectors[0]);
#endif

  // TODO - do multiplication in place
  // _externalForces += _forceMultiplier * subbasis * _forceVectors[0]
  subbasis.gemvInplace( _forceMultiplier, _forceVectors[0],
                        //_externalForces, 1.0 );
                        _externalForces, 1.0, true /* transpose */ );
  // TODO - end
  
  // process the rest of the forces
  for (unsigned int x = 1; x < _forceVectors.size(); x++)
  {
#if 0
    subbasis = _tetMesh->vertexSubBasis(_forceNodes[x]);
    subbasis = subbasis.transpose();
    _externalForces += subbasis.gemv(_forceMultiplier, _forceVectors[x]);
#endif

    // TODO - try to do everything in place
    //_tetMesh->vertexSubBasis( _forceNodeIDs[x], subbasis, true /* transpose */ );
    _tetMesh->vertexSubBasis( _forceNodeIDs[x], subbasis );

    // _externalForces += _forceMultiplier * subbasis * _forceVectors[x]
    subbasis.gemvInplace( _forceMultiplier, _forceVectors[x],
                          //_externalForces, 1.0 );
                          _externalForces, 1.0, true /* transpose */ );
    // TODO - end
  }

  // Correct the force magnitude, if necessary
  // Try doing this by forming a weightes norm
  Real forceMagnitude = 0.0;
  VECTOR &eigs = _tetMesh->eigenvalues();
  for ( int i = 0; i < _externalForces.size(); i++ )
  {
    forceMagnitude += _externalForces[i] * _externalForces[i] * eigs[i];
  }

  if ( MAX_EXTERNAL_FORCE > 0.0 && forceMagnitude > MAX_EXTERNAL_FORCE )
  {
    // Clamp it!
    Real forceRatio = MAX_EXTERNAL_FORCE / forceMagnitude;
    _externalForces *= forceRatio;
  }

  _forceVectors.clear();
  _forceNodes.clear();
  _forceNodeIDs.clear();
}

//////////////////////////////////////////////////////////////////////
// poke the mesh in a deterministic way for debugging
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::poke()
{
  /*
  VEC3F* vertex = _tetMesh->vertices(100);
  VEC3F direction(0.0f,100.0f,0.0f);

  _forceVectors.push_back(direction);
  _forceNodes.push_back(vertex);
  */
  VEC3F clicked(0.462357, 0.567332, 0.883899);
  VEC3F dragged(0.540074, -0.359055, 0.654308);
  VEC3F direction = dragged - clicked;
  VEC3F* node = _tetMesh->closestSurfaceNode(clicked);
  _forceVectors.push_back(direction);
  _forceNodes.push_back(node);
  _forceNodeIDs.push_back(_tetMesh->vertexID(node));
}

//////////////////////////////////////////////////////////////////////
// pull the mesh in a deterministic way for debugging
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::pull()
{
  VEC3F clicked(0.462357, 0.567332, 0.883899);
  VEC3F dragged(0.540074, -0.359055, 0.654308);
  VEC3F direction = dragged - clicked;

  // scale pulling by the current timestep
  //direction *= (Real)_totalSteps / 100.0;
  //direction *= 10.0;

  // cache the first one we see
  if (_pulledNode == NULL)
    _pulledNode = _tetMesh->closestSurfaceNode(clicked);
  _forceVectors.push_back(direction);
  _forceNodes.push_back(_pulledNode);
  _forceNodeIDs.push_back(_tetMesh->vertexID(_pulledNode));
}

//////////////////////////////////////////////////////////////////////
// Draw the clicked node
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::drawClickedNode()
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
void SUBSPACE_INTEGRATOR::drawForceVector()
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
// Print out a detailed timing breakdown
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::printTimingBreakdown()
{
  // create an inverse map so that it will sort by time
  map<double, string> inverseMap;
  map<string, double>::iterator forwardIter;
  for (forwardIter = _timingBreakdown.begin(); forwardIter != _timingBreakdown.end(); forwardIter++)
    inverseMap[forwardIter->second] = forwardIter->first;

  // print the map out backwards since it sorts from least to greatest
  cout << "SUBSPACE_INTEGRATOR TIMING BREAKDOWN: " << endl;
  cout << "===============================================================================" << endl;
  map<double,string>::reverse_iterator backwardIter;
  double totalSeen = 0.0;
  for (backwardIter = inverseMap.rbegin(); backwardIter != inverseMap.rend(); backwardIter++)
  {
    string name = (*backwardIter).second + string("               ");
    name = name.substr(0,15);

    cout << "[" << (*backwardIter).first / _totalTime * 100.0 << "%\t]: "
         << name.c_str() << "\t" << (*backwardIter).first / _totalSteps << "s / frame" << endl;
    totalSeen += (*backwardIter).first;
  }
  Real misc = (_totalTime - totalSeen) / _totalTime * 100.0;
  cout << "[" << misc << "%\t]: " << "Misc. " << endl;
  cout << "===============================================================================" << endl;
  cout << " Average CG iterations: " << _A.meanIterations() << endl;
  cout << " Current FPS: " << _totalSteps / _totalTime << endl;
  cout << "===============================================================================" << endl;
  if (misc < 0.0)
  {
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << " BREAKDOWN ADDS UP TO MORE THAN 100! TIMERS ARE OVERLAPPING! " << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// solve using the invertible finite element stiffness matrix
// from "Robust Quasisatic Finite Elements and Flesh Simulation"
// Teran, Sifakis, Irving, Fedkiw, SCA 2005
//
// There is a lot of magic that has to occur to do the matrix-vector
// multiply, so the whole CG solve is done explicitly here.
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::solveInvertibleStiffness(VECTOR& x, VECTOR& b, bool quasistatic)
{
  // get the temp CG arrays
  VECTOR& residual = _A.residual();
  VECTOR& direction = _A.direction();
  VECTOR& q = _A.q();
  Real eps = _A.eps();
  Real iterations = _A.maxIterations();

  // resize if necessary
  if (residual.size() != _rank)
  {
    residual.resizeAndWipe(_rank);
    direction.resizeAndWipe(_rank);
    q.resizeAndWipe(_rank);
  }
  
	// r = b - Ax
  residual = _A * x;
  residual = b - residual;

	// d = r
  direction.equals(residual);
  
	// deltaNew = transpose(r) * r
  Real deltaNew = residual ^ residual;
  
	// delta0 = deltaNew
  //Real delta0 = deltaNew;

	// While deltaNew > (eps^2) * delta0
	Real maxR = 2.0f * eps;
  int i = 0;
	while ((i < iterations) && (maxR > eps))
	{
    // hideous matrix-vector multiply
    invertibleStiffnessMultiply(direction, q);

    // if this is an implicit step, multiply the mass and damping in 
    // as well
    if (!quasistatic)
      q += _A * direction;
    
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
  
  // add the stats to the _A matrix so that it gets passed on to 
  // printTimingBreakdown
  _A.totalSolves()++;
  _A.totalIterations() += i;
}

//////////////////////////////////////////////////////////////////////
// Cache the deformation gradient diagonalizations
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::cacheKryslDiagonalizations()
{
  unsigned int totalTets = _tetMesh->totalTets();
  if (_Us.size() != totalTets)
  {
    _Us.resize(totalTets);
    _Vs.resize(totalTets);
    _Fhats.resize(totalTets);
    _stiffnesses.resize(totalTets);

    for ( unsigned int i = 0; i < totalTets; i++ )
    {
      _stiffnesses[ i ].resizeAndWipe( 9, 9 );
    }
  }

  // if the material isn't INVERTIBLE, you did something really bad
  vector<TET>& tets = _tetMesh->tets();
#if USING_OPENMP
#pragma omp parallel
#endif
  { 
#if USING_OPENMP
  int id  = omp_get_thread_num();
#pragma omp for  schedule(static)
#else
  const int id  = 0;
#endif
  for (unsigned int x = 0; x < totalTets; x++)
  {
    int materialIndex = tets[x].materialIndex();
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

#if 0
    MATRIX stiffness(9,9,stiffnessData);
    _stiffnesses[x] = stiffness;
#endif
    _stiffnesses[x].copyInplace( stiffnessData, 9, 9 );
  }
  } //OMP
}

//////////////////////////////////////////////////////////////////////
// Cache the deformation gradient diagonalizations, but only for the
// key tets.
//
// If checkTetState == true, then we will use the updateTetDeformation
// function to check whether or not a key tet needs to have its
// stiffness matrix contribution updated.  If not, we will not
// bother caching its stiffness density.
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::cacheDiagonalizations( bool checkTetState )
{
  unsigned int totalKeyTets = _tetMesh->totalKeyTets();
  if (_Us.size() != totalKeyTets)
  {
    _Us.resize(totalKeyTets);
    _Vs.resize(totalKeyTets);
    _Fhats.resize(totalKeyTets);
    _stiffnesses.resize(totalKeyTets);

    for ( unsigned int i = 0; i < totalKeyTets; i++ )
    {
      _stiffnesses[ i ].resizeAndWipe( 9, 9 );
    }
  }

  // if the material isn't INVERTIBLE, you did something really bad
  vector<TET>& tets = _tetMesh->tets();
  int*& keyTets = _tetMesh->keyTets();
  Real* precomputedF = _tetMesh->F().data();

  // This version is faster, but causes the online integrator to die
#if 0
  for ( int materialIndex = 0;
        materialIndex < _tetMesh->totalMaterials(); materialIndex++ )
  {
    for ( int id = 0; id < _tetMesh->totalCores(); id++ )
    {
      if ( _tetMesh->materialCopies()[id]
        && _tetMesh->materialCopies()[id][materialIndex] )
      {
        INVERTIBLE* material
          = (INVERTIBLE*)(_tetMesh->materialCopies()[id][materialIndex]);

        material->stiffnessDensityTime() = 0.0;
        material->materialStiffnessDensityTime() = 0.0;
      }
    }
  }

  //TIMER stiffnessTimer;

#undef USING_OPENMP

  // FIXME
  //INVERTIBLE::NUMINVERTED = 0;

#if USING_OPENMP  
#pragma omp parallel 
  {
	  const int num = omp_get_num_threads();
    const int id = omp_get_thread_num();
    unsigned int begin = (totalKeyTets / num) * id;
    unsigned int end   = (totalKeyTets / num) * (id+1);
    if (id == num-1)
      end = totalKeyTets;
    for (unsigned int x = begin; x < end; x++)
    {
      TET& tet = tets[keyTets[x]];
      int materialIndex = tet.materialIndex();
      INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materialCopies()[id][materialIndex]);
      
      // give the precomputed F as raw data to MATRIX3
      MATRIX3 F(&precomputedF[9 * x]); 
#if 0
      MATRIX3 U;
      MATRIX3 Fhat;
      MATRIX3 V;
      Real stiffnessData[81];
      material->diagonalizeF(F, U, Fhat, V);
      material->stiffnessDensity(U, Fhat, V, stiffnessData);

      _Us[x] = U;
      _Vs[x] = V;
      _Fhats[x] = Fhat;

#if 0
      MATRIX stiffness(9,9,stiffnessData);
      _stiffnesses[x] = stiffness;
#endif
      _stiffnesses[x].copyInplace( stiffnessData, 9, 9 );
#endif

      material->diagonalizeF( F, _Us[x], _Fhats[x], _Vs[x] );

      if ( checkTetState && !updateTetDeformation( x, _deformationTolerance,
                                                   F, _Us[x], _Fhats[x], _Vs[x] ) )
      {
        _updateKeyTetDeformation[ x ] = false;

        continue;
      }

      _updateKeyTetDeformation[ x ] = true;

      material->stiffnessDensity( _Us[x], _Fhats[x], _Vs[x],
                                  _stiffnesses[x].data() );
    }
  }

#else
  for (unsigned int x = 0; x < totalKeyTets; x++)
  {
    TET& tet = tets[keyTets[x]];
    int materialIndex = tet.materialIndex();
    INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materialCopies()[0][materialIndex]);
    
    // give the precomputed F as raw data to MATRIX3
    MATRIX3 F(&precomputedF[9 * x]); 

#if 0
    MATRIX3 U;
    MATRIX3 Fhat;
    MATRIX3 V;
    Real stiffnessData[81];
    material->diagonalizeF(F, U, Fhat, V);
    material->stiffnessDensity(U, Fhat, V, stiffnessData);

    _Us[x] = U;
    _Vs[x] = V;
    _Fhats[x] = Fhat;

#if 0
    MATRIX stiffness(9,9,stiffnessData);
    _stiffnesses[x] = stiffness;
#endif
    _stiffnesses[x].copyInplace( stiffnessData, 9, 9 );
#endif

    //TIMER diagonalizeTimer;
    material->diagonalizeF( F, _Us[x], _Fhats[x], _Vs[x] );
    //_timingBreakdown["Diagonalize F"] += diagonalizeTimer.timing();

    if ( checkTetState && !updateTetDeformation( x, _deformationTolerance,
                                                 F, _Us[x], _Fhats[x], _Vs[x] ) )
    {
      _updateKeyTetDeformation[ x ] = false;

      continue;
    }

    _updateKeyTetDeformation[ x ] = true;

    material->stiffnessDensity( _Us[x], _Fhats[x], _Vs[x],
                                _stiffnesses[x].data() );
  }
#endif

  //cout << INVERTIBLE::NUMINVERTED << " inverted elements" << endl;

  //_timingBreakdown["Stiffness Factorization"] += stiffnessTimer.timing();
//#define USING_OPENMP 1

  // Add material timing info
  for ( int materialIndex = 0;
        materialIndex < _tetMesh->totalMaterials(); materialIndex++ )
  {
    for ( int id = 0; id < _tetMesh->totalCores(); id++ )
    {
      if ( _tetMesh->materialCopies()[id]
        && _tetMesh->materialCopies()[id][materialIndex] )
      {
        INVERTIBLE* material
          = (INVERTIBLE*)(_tetMesh->materialCopies()[id][materialIndex]);

        //_timingBreakdown["Invertible Stiffness Density"] += material->stiffnessDensityTime();
        //_timingBreakdown["Material Stiffness Density"] += material->materialStiffnessDensityTime();
      }
    }
  }
#endif

  for (int x = 0; x < totalKeyTets; x++)
  {
    TET& tet = tets[keyTets[x]];
    int materialIndex = tet.materialIndex();
    INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materials()[materialIndex]);
    
    // give the precomputed F as raw data to MATRIX3
    MATRIX3 F(&precomputedF[9 * x]); 
    MATRIX3 U;
    MATRIX3 Fhat;
    MATRIX3 V;
    Real stiffnessData[81];
    TIMER svdTimer;
    material->diagonalizeF(F, U, Fhat, V);
    _timingBreakdown["SVD"] += svdTimer.timing();

    TIMER stiffnessTimer;
    material->stiffnessDensity(U, Fhat, V, stiffnessData);
    _timingBreakdown["Stiffness Factorization"] += stiffnessTimer.timing();

    TIMER diagonalTimer;
    _Us[x] = U;
    _Vs[x] = V;
    _Fhats[x] = Fhat;

    MATRIX stiffness(9,9,stiffnessData);
    _stiffnesses[x] = stiffness;
    _timingBreakdown["Misc. Diagonalization"] += diagonalTimer.timing();
  }
}

//////////////////////////////////////////////////////////////////////
// Do the matrix-vector multiply for solveInvertibleStiffness()
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::invertibleStiffnessMultiply(VECTOR& direction,
                                                      VECTOR& answer,
                                                      int verbose)
{
  // wipe any old answer
  answer.clear();
  
  // if the material isn't INVERTIBLE, you did something really bad
  int totalKeyTets = _tetMesh->totalKeyTets();
  vector<TET>& tets = _tetMesh->tets();
  int*& keyTets = _tetMesh->keyTets();
  VECTOR& keyWeights = _tetMesh->keyWeights();
  //Real* precomputedF = _tetMesh->F().data();

  for (int x = 0; x < totalKeyTets; x++)
  {
    TET& tet = tets[keyTets[x]];
    int materialIndex = tet.materialIndex();
    
    INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materials()[materialIndex]);

    // get cached stiffness matrix
    MATRIX& stiffnessMatrix = _stiffnesses[x];
    MATRIX3& U = _Us[x];
    MATRIX3& V = _Vs[x];

    // get the subbasis for the tet
    MATRIX subbasis = _tetMesh->tetSubBasis(&tet);

    // get PFPu from the tet - only depends on rest state,
    // so don't need to update tet state at all
    MATRIX pFpu(9,12);
    material->computePFPu(tet, pFpu);
    
    // unproject the direction using the tet subbasis
    // (not so bad, subbasis is only (12 x r))
    VECTOR deltaX = subbasis * direction;

    // compute deltaF
    VECTOR deltaF = pFpu * deltaX;

    // rotate deltaF
    MATRIX3 rotated = U.transpose() * TET::repackF(deltaF) * V;
    deltaF = TET::flattenF(rotated);

    VECTOR contraction = stiffnessMatrix * deltaF;
    MATRIX3 deltaP = TET::repackF(contraction);

    // rotate deltaP back
    deltaP = U * deltaP * V.transpose();

    VECTOR flatP = TET::flattenF(deltaP);
    VECTOR forces = pFpu ^ flatP;
   
    // project the basis using the tet subbasis
    VECTOR deltaG = subbasis ^ forces;

    // add to the answer 
    // (taking into account negation of stiffness matrix)
    answer -= keyWeights(x) * deltaG;
  }
}

//////////////////////////////////////////////////////////////////////
// debug the stiffness matrix-vector multiply
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::debugStiffness()
{
  int rank = _tetMesh->rank();
  VECTOR q(rank);
  VECTOR direction(rank);

  for (int x = 0; x < rank; x++)
    direction(x) = x;

  cout << " direction: " << direction << endl;
  
  // q = Ad
  q = _A * direction;
  cout << " gold standard q: " << q << endl;

  // try a naive multiply here, just to see what the stages should
  // look like
  int totalKeyTets = _tetMesh->totalKeyTets();
  VECTOR goldSum(rank);
  Real* precomputedF = _tetMesh->F().data();
  goldSum.clear();
  for (int x = 0; x < totalKeyTets; x++)
  {
    TET& tet = (_tetMesh->tets())[_tetMesh->keyTets()[x]];
    int materialIndex = tet.materialIndex();

    Real keyWeight = _tetMesh->keyWeights()(x);
    INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materials()[materialIndex]);
    MATRIX subbasis = _tetMesh->tetSubBasis(&tet);
    MATRIX PFPu(9,12);
    material->computePFPu(tet, PFPu);

    // compute the stiffness matrix for the first tet
    MATRIX3 F(&precomputedF[9 * x]); 
    VECTOR flatF = TET::flattenF(F);
    Real stiffnessDensity[81];
    material->stiffnessDensity(flatF.data(), stiffnessDensity);
    MATRIX stiffness(9,9,stiffnessDensity);

    // unproject to get the displacement vector of the tet
    VECTOR deltaX = subbasis * direction;

    // convert it to deformation space
    VECTOR deltaF = PFPu * deltaX;

    // multiply by stiffness matrix
    VECTOR deltaP = stiffness * deltaF;
    deltaP *= -1.0;
    
    // get out of deformation space
    VECTOR forces = PFPu ^ deltaP;

    // project final deltaG
    VECTOR deltaG = subbasis ^ forces;

    goldSum += keyWeight * deltaG;
   
    /*
    cout << " Gold standard tet " << x << " ------ " << endl;
    cout << " gold deltaX: " << deltaX << endl;
    cout << " gold PFPu:   " << PFPu << endl;
    cout << " gold F: " << F << endl;
    cout << " gold deltaF: " << deltaF << endl;
    cout << " gold stiffness: " << stiffness << endl;
    cout << " gold deltaP: " << deltaP << endl;
    cout << " gold forces: " << forces << endl;
    cout << " gold deltaG: " << deltaG << endl;
    cout << " gold keyWeight: " << keyWeight << endl;
    cout << " gold final for tet " << x << ": " << keyWeight * deltaG << endl;
    cout << " gold current sum: " << goldSum << endl;
    */
    cout << endl;
  }
  
  VECTOR debug(q.size());
  invertibleStiffnessMultiply(direction, debug, true);
  cout << " q returned by the multiply: " << debug << endl;
  cout << endl;
  exit(0);
}

//////////////////////////////////////////////////////////////////////
// Form the A and b so that they can be passed out for solution
// by the partitioned integrator -- implicit dynamics case
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::initializeImplicitStep()
{
  TIMER preamble;

  // accumulate external forces
  computeExternalForces();
  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " EXTERNAL FORCES SPECIFIED FROM FILE " << endl;

  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld = _velocity;
  _accelerationOld = _acceleration;

  _tetMesh->inertiaTensorOld() = _tetMesh->inertiaTensor();
  _tetMesh->inertiaTensorDtOld() = _tetMesh->inertiaTensorDt();
  _timingBreakdown["Preamble"] += preamble.timing();
}

//////////////////////////////////////////////////////////////////////
// Update the integrator and mesh state using and acceleration
// level update. Does not update the "olds" however -- this should
// have been done at the beginning of a time step.
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::updateState()
{
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);
  
	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], _accelerationOld);
  _acceleration.axpy(-_alpha[1], _velocityOld);
  _acceleration.axpy(-_alpha[0], _positionOld);
  _acceleration.axpy(_alpha[0], _position);

  VECTOR& q = _tetMesh->q();
  q.equals(_position);
}

//////////////////////////////////////////////////////////////////////
// Update the integrator and mesh state using and acceleration
// level update. Does not update the "olds" however -- this should
// have been done at the beginning of a time step.
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::updateStateUsingAcceleration()
{
  Real* alpha = _accelerationAlpha;

  _position.equals(_positionOld);
  _position.axpy(alpha[2], _velocityOld);
  _position.axpy(alpha[3], _accelerationOld);
  _position.axpy(alpha[4], _acceleration);

  //_position = _positionOld + alpha[2] * _velocityOld +
  //                           alpha[3] * _accelerationOld + 
  //                           alpha[4] * _acceleration;

  _velocity.equals(_velocityOld); 
  _velocity.axpy(alpha[0], _accelerationOld);
  _velocity.axpy(alpha[1], _acceleration);
  //_velocity = _velocityOld + alpha[0] * _accelerationOld +
  //                           alpha[1] * _acceleration;

  VECTOR& q = _tetMesh->q();
  q.equals(_position);

  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " DEACTIVATING VARIABLE MASS " << endl;
  _tetMesh->refreshSitBar();
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::updateStateUsingVelocity()
{
  Real* alpha = _velocityAlpha;

  _acceleration = alpha[0] * _velocity +
                  alpha[1] * _velocityOld +
                  alpha[2] * _accelerationOld;
  _position = _positionOld + alpha[3] * _velocityOld +
                             alpha[4] * _velocity + 
                             alpha[5] * _accelerationOld;

  VECTOR& q = _tetMesh->q();
  q.equals(_position);

  _tetMesh->refreshSitBar();
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an implicit step -- should be called 
// right after initializeImplicitStep
//
// Note that this does not load the mass matrix into A, or add
// the MA term to the residual!
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::generateImplicitAccelerationMatrices(map<string, double>& timingBreakdown, bool useKrysl)
{
  Real* alpha = _accelerationAlpha;
  MATRIX& M = _tetMesh->reducedMass();

  // compute the new deformation gradient F (needed by both R and K)
  //TIMER timerF;
  _tetMesh->SUBSPACE_TET_MESH::generateF();
  //timingBreakdown["Generate F"] += timerF.timing();

  // precompute diagonalizations --
  // this call is not timed out because timers are inside the function call
  //TIMER diagonalizeTimer;
  cacheDiagonalizations();
  //timingBreakdown["Diagonalize F"] += diagonalizeTimer.timing();

  // get the reduced internal forces
  //TIMER internalTimer;
  if (useKrysl)
  {
    finalizeImplicitAccelerationStep();
    _tetMesh->updateFullMesh();
    _tetMesh->generateKryslInternalForces();
  }
  else
    _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);

  VECTOR& R = _tetMesh->reducedInternalForce();
  //timingBreakdown["Internal Forces"] += internalTimer.timing();

  // get the reduced stiffness matrix -- do this even if we are 
  // doing the conjugate gradient solve, because the damping matrix
  // needs it
  if (useKrysl)
    _tetMesh->generateKryslStiffnessMatrix();
  else
    _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, timingBreakdown);
  MATRIX& K = _tetMesh->stiffness();
  
  // compute the damping matrix, but only if a new K is available
  //TIMER RHSTimer;
  // disabling dynamic damping as long as comparing to FULLSPACE_INTEGRATOR::generateImplicitMatrices
  /*
  _CDamping.clearingAxpy(_rayleighAlpha, M);
  _CDamping.axpy(_rayleighBeta, K);
  */

  //_temp = M * _acceleration;
  //_residual = _CDamping * _velocity;
  _CDamping.gemvInplace(1.0, _velocity, _residual, 0.0);

  _residual -= R;
  _residual -= _externalForces;

  // M acceleration is left out because the multibody solver adds in the MA term later
  //_residual += M * _acceleration;

  //timingBreakdown["RHS Assembly"] += RHSTimer.timing();

  // use the acceleration level Newmark consts
  //TIMER ATimer;
  //_A = alpha[4] * K;
  _A.clearingAxpy(alpha[4], K);
  _A.axpy(alpha[1], _CDamping);

  // M is left out because the multibody solver adds in the MA term later
  //_A.axpy(1.0, M);

  //timingBreakdown["A Assembly"] += ATimer.timing();

  /*
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " reduced acceleration residual: " << _tetMesh->U() * _residual << endl;

  VECTOR velocity = _alpha[5] * _accelerationOld + _alpha[4] * _velocityOld;
  velocity += -_alpha[3] * _positionOld + _alpha[3] * _position;
  VECTOR accel = _alpha[0] * _position - _alpha[0] * _positionOld;
  accel += -_alpha[1] * _velocityOld - _alpha[2] * _accelerationOld;
  VECTOR residual = _CDamping * velocity;
  residual += M * accel;
  residual -= R;
  residual -= _externalForces;

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  cout << " reduced position residual: " << _tetMesh->U() * _residual << endl;
  */
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an implicit step -- should be called 
// right after initializeImplicitStep
//
// this has been forked out into its own function in case we want
// to do more than one Newton step
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::generateImplicitMatrices(map<string, double>& timingBreakdown, bool useKrysl)
{
  TIMER updateMesh;
  // get the reduced mass matrix
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
  cacheDiagonalizations();

  // get the reduced internal forces
  TIMER internalTimer;
  //VECTOR& R = _tetMesh->generateInternalForces();
  if (useKrysl)
  {
    finalizeImplicitAccelerationStep();
    _tetMesh->updateFullMesh();
    _tetMesh->generateKryslInternalForces();
  }
  else
    _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);

  VECTOR& R = _tetMesh->reducedInternalForce();
  timingBreakdown["Internal Forces"] += internalTimer.timing();

  // get the reduced stiffness matrix -- do this even if we are 
  // doing the conjugate gradient solve, because the damping matrix
  // needs it
  TIMER stiffnessTimer;
  if (useKrysl)
    _tetMesh->generateKryslStiffnessMatrix();
  else
    _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, timingBreakdown);
  //MATRIX& K = _tetMesh->generateStiffnessMatrix();
  MATRIX& K = _tetMesh->stiffness();
  timingBreakdown["Generate Stiffness Matrix"] += stiffnessTimer.timing();
  
  TIMER RHSTimer;
  _CDamping.clearingAxpy(_rayleighAlpha, M);
  _CDamping.axpy(_rayleighBeta, K);

  // compute the LHS of the residual:
  // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
  _temp.clearingAxpy(-_alpha[2], _accelerationOld);
  _temp.axpy(-_alpha[1], _velocityOld);
  _temp.axpy(-_alpha[0], _positionOld);
  _temp.axpy(_alpha[0], _position);
  _temp = M * _temp;

  /*
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
  */
  // compute the RHS of the residual:
  // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
  _totalMeshForces.clearingAxpy(_alpha[5], _accelerationOld);
  _totalMeshForces.axpy(_alpha[4], _velocityOld);
  _totalMeshForces.axpy(-_alpha[3], _positionOld);
  _totalMeshForces.axpy(_alpha[3], _position);
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
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an implicit step -- should be called 
// right after initializeImplicitStep
//
// this has been forked out into its own function in case we want
// to do more than one Newton step
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::generateImplicitMatricesInvertible(map<string, double>& timingBreakdown)
{
  TIMER updateMesh;
  // get the reduced mass matrix
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
#if 1
  MATRIX& K = _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, timingBreakdown);
#else
  TIMER generateTimer;
  map<string,double> dummy;
  MATRIX& K = _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, dummy);
  timingBreakdown["Generate K"] += generateTimer.timing();
#endif
  
  TIMER RHSTimer;
  _CDamping.clearingAxpy(_rayleighAlpha, M);
  _CDamping.axpy(_rayleighBeta, K);

  // compute the LHS of the residual:
  // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
  _temp.clearingAxpy(-_alpha[2], _accelerationOld);
  _temp.axpy(-_alpha[1], _velocityOld);
  _temp.axpy(-_alpha[0], _positionOld);
  _temp.axpy(_alpha[0], _position);
  _temp = M * _temp;

  // compute the RHS of the residual:
  // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
  _totalMeshForces.clearingAxpy(_alpha[5], _accelerationOld);
  _totalMeshForces.axpy(_alpha[4], _velocityOld);
  _totalMeshForces.axpy(-_alpha[3], _positionOld);
  _totalMeshForces.axpy(_alpha[3], _position);
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
}

//////////////////////////////////////////////////////////////////////
// Generate multibody-specific, acceleration-based matrices
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::generateMultibodyMatrices()
{
  // get the reduced mass matrix
  MATRIX& M = _tetMesh->reducedMass();

  // compute the new deformation gradient F (needed by both R and K)
  _tetMesh->generateF();

  // precompute diagonalizations --
  // this call is not timed out because timers are inside the function call
  cacheDiagonalizations();

  // get the reduced internal forces
  VECTOR& R = _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);

  // get the reduced stiffness matrix -- do this even if we are 
  // doing the conjugate gradient solve, because the damping matrix
  // needs it
  map<string, double> dummy;
  MATRIX& K = _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, dummy);
  
  _CDamping.clearingAxpy(_rayleighAlpha, M);
  _CDamping.axpy(_rayleighBeta, K);

  _totalMeshForces.clear();
  _totalMeshForces += M * _acceleration;
  _totalMeshForces += _CDamping * _velocity;

  // assemble the mesh forces LHS + RHS + R
  _totalMeshForces -= R;

  // assemble full residual: LHS + RHS + R - F
  _residual = _totalMeshForces;
  _residual -= _externalForces;
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an implicit step -- should be called 
// right after initializeImplicitStep
//
// this has been forked out into its own function in case we want
// to do more than one Newton step
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::generateImplicitMatricesDebug()
{
  // get the reduced mass matrix
  MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // update tet mesh with the new q vector
  // an assignment is overly expensive -- the "=" forces an allocation
  TIMER updateMesh;
  q.equals(_position);
  _timingBreakdown["Update Mesh"] += updateMesh.timing();

  // compute the new deformation gradient F (needed by both R and K)
  _tetMesh->generateF();

  // precompute diagonalizations --
  // this call is not timed out because timers are inside the function call
  cacheDiagonalizations();

  // get the reduced internal forces
  TIMER internalTimer;
  //VECTOR& R = _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);
  //VECTOR& R = _tetMesh->generateInternalForces();
  VECTOR& R = _tetMesh->generateInternalForcesDebug();
  _timingBreakdown["Internal Forces"] += internalTimer.timing();

  // get the reduced stiffness matrix -- do this even if we are 
  // doing the conjugate gradient solve, because the damping matrix
  // needs it
  MATRIX& K = _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
  //MATRIX& K = _tetMesh->generateStiffnessMatrix();
  
  TIMER RHSTimer;
  _CDamping.clearingAxpy(_rayleighAlpha, M);
  _CDamping.axpy(_rayleighBeta, K);

  // compute the LHS of the residual:
  // M (alpha_1 (q_{i+1} - q_{i}) - alpha_2(q^{dot}_{i} - alpha_3(q^{dot dot}_{i})))
  _temp.clearingAxpy(-_alpha[2], _accelerationOld);
  _temp.axpy(-_alpha[1], _velocityOld);
  _temp.axpy(-_alpha[0], _positionOld);
  _temp.axpy(_alpha[0], _position);
  _temp = M * _temp;

  /*
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
  */
  // compute the RHS of the residual:
  // C (alpha_4 (q_{i+1} - q_{i}) + alpha_5 q^{dot}_{i} + alpha_6(q^{dot dot}_{i}))
  _totalMeshForces.clearingAxpy(_alpha[5], _accelerationOld);
  _totalMeshForces.axpy(_alpha[4], _velocityOld);
  _totalMeshForces.axpy(-_alpha[3], _positionOld);
  _totalMeshForces.axpy(_alpha[3], _position);
  _totalMeshForces = _CDamping * _totalMeshForces;

  // assemble the mesh forces LHS + RHS + R
  _totalMeshForces += _temp;
  _totalMeshForces -= R;

  // assemble full residual: LHS + RHS + R - F
  _residual = _totalMeshForces;
  _residual -= _externalForces;
  _timingBreakdown["RHS Assembly"] += RHSTimer.timing();

  // form A
  _A = K;
  _A.axpy(_alpha[0], M);
  _A.axpy(_alpha[3], _CDamping);
}

//////////////////////////////////////////////////////////////////////
// Form the A and b so that they can be passed out for solution
// by the partitioned integrator -- quasistatic case
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::initializeQuasistaticStep()
{
  TIMER preamble;

  // accumulate external forces
  computeExternalForces();
  
  // copy the current values into the old values
  _positionOld = _position;

  _timingBreakdown["Preamble"] += preamble.timing();
}

//////////////////////////////////////////////////////////////////////
// Generate the matrices for an quasistatic step -- should be called 
// right after initializeQuasistaticStep
//
// this has been forked out into its own function in case we want
// to do more than one Newton step
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::generateQuasistaticMatrices(map<string, double>& timingBreakdown)
{
  VECTOR& q = _tetMesh->q();
  //TIMER updateMesh;

  // update q vector with current tet mesh state
  // an assignment is overly expensive -- the "=" forces an allocation
  q.equals(_position);
  //timingBreakdown["Update Mesh"] += updateMesh.timing();
  
  // compute the new deformation gradient F (needed by both R and K)
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
// Update the mesh after the partitioned implicit step
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::finalizeImplicitStep()
{
  _velocity.clearingAxpy(_alpha[5], _accelerationOld);
  _velocity.axpy(_alpha[4], _velocityOld);
  _velocity.axpy(-_alpha[3], _positionOld);
  _velocity.axpy(_alpha[3], _position);
  
	// update acceleration
  _acceleration.clearingAxpy(-_alpha[2], _accelerationOld);
  _acceleration.axpy(-_alpha[1], _velocityOld);
  _acceleration.axpy(-_alpha[0], _positionOld);
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
void SUBSPACE_INTEGRATOR::finalizeImplicitAccelerationStep()
{
  _velocity = _velocityOld + 
              _accelerationAlpha[0] * _accelerationOld +
              _accelerationAlpha[1] * _acceleration;

  _position = _positionOld +
              _accelerationAlpha[2] * _velocityOld +
              _accelerationAlpha[3] * _accelerationOld +
              _accelerationAlpha[4] * _acceleration;

  _tetMesh->q() = _position;
  _tetMesh->qOld() = _positionOld;  
}

//////////////////////////////////////////////////////////////////////
// compute position and velocity using final velocity values
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::finalizeImplicitVelocityStep()
{
  _position = _positionOld +
              _velocityAlpha[3] * _velocityOld +
              _velocityAlpha[4] * _accelerationOld +
              _velocityAlpha[5] * _velocity;
  
  _acceleration = 
              _velocityAlpha[0] * _velocity +
              _velocityAlpha[1] * _velocityOld +
              _velocityAlpha[2] * _accelerationOld;

  _tetMesh->q() = _position;
  _tetMesh->qOld() = _positionOld;  
}

//////////////////////////////////////////////////////////////////////
// Update the mesh after the partitioned quasistatic step
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::finalizeQuasistaticStep()
{
  // update node positions
  _tetMesh->updateSurfaceMesh();
  _externalForces.clear();
  _totalSteps++;
}

//////////////////////////////////////////////////////////////////////
// update quantities after the tet mesh basis has been altered
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::updateBasis(MATRIX& oldBasis, MATRIX& newBasis)
{
  _rank = newBasis.cols();

  if (oldBasis.cols() == 0)
  {
    _position.resizeAndWipe(_rank);
    _velocity.resizeAndWipe(_rank);
    _acceleration.resizeAndWipe(_rank);
    _positionOld.resizeAndWipe(_rank);
    _velocityOld.resizeAndWipe(_rank);
    _accelerationOld.resizeAndWipe(_rank);
    _residual.resizeAndWipe(_rank);
    _temp.resizeAndWipe(_rank);
    _externalForces.resizeAndWipe(_rank);
    _gravity.resizeAndWipe(_rank);
  }
  else
  {
    MATRIX changeOfBasis = newBasis ^ oldBasis;

    _position = changeOfBasis * _position;
    _velocity = changeOfBasis * _velocity;
    _acceleration = changeOfBasis * _acceleration;
    _positionOld = changeOfBasis * _positionOld;
    _velocityOld = changeOfBasis * _velocityOld;
    _accelerationOld = changeOfBasis * _accelerationOld;
    _residual = changeOfBasis * _residual;
    _temp = changeOfBasis * _temp;
    _externalForces = changeOfBasis * _externalForces;
    _gravity = changeOfBasis * _gravity;
  }

  // precompute the gravity vector
  int nodes = _tetMesh->unconstrainedNodes();
  VECTOR fullGravity(nodes * 3);
  VEC3F gravityVector = _gravityMagnitude * _gravityDown;
  for (int x = 0; x < nodes; x++)
  {
    fullGravity(3 * x) = gravityVector[0];
    fullGravity(3 * x + 1) = gravityVector[1];
    fullGravity(3 * x + 2) = gravityVector[2];
  }
  fullGravity = _tetMesh->massMatrix() * fullGravity;

  _gravity = newBasis ^ fullGravity;
  
  _CDamping.resizeAndWipe(_rank, _rank);
  _A.resizeAndWipe(_rank, _rank);
}

//////////////////////////////////////////////////////////////////////
// update quantities after the cubature points have been altered
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::updateCubature(bool invertible)
{
  // backup the q values
  VECTOR qBackup = _tetMesh->q();
  
  // get the rest pose damping matrix in case we want to do just
  // linear Rayleigh damping
  _tetMesh->q().clear();
  if (_useKryslStiffness)
  {
    _tetMesh->updateFullMesh();
    _tetMesh->TET_MESH::generateF();
  }
  else
    _tetMesh->generateF();

  _CDamping.clearingAxpy(_rayleighAlpha, _tetMesh->reducedMass());

  // check that a cubature was in fact loaded by checking the size
  // of the stiffness matrix
  if (invertible)
  {
    if (_useKryslStiffness)
    {
      cacheKryslDiagonalizations();
      _tetMesh->generateKryslStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
    }
    else
    {
      cacheDiagonalizations();
      _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
    }
  }
  else
  {
    if (_useKryslStiffness)
      _tetMesh->generateKryslStiffnessMatrix();
    else
    {
      _tetMesh->generateStiffnessMatrix();
    }
  }
  MATRIX& stiffness = _tetMesh->stiffness();
  if (stiffness.rows() > 0)
    _CDamping.axpy(_rayleighBeta, stiffness);

  // restore the q values
  _tetMesh->q() = qBackup;
}

//////////////////////////////////////////////////////////////////////
// stomp everything related to the reduced model
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::resetBasis()
{
  _rank = 0;

  _position.resizeAndWipe(0);
  _velocity.resizeAndWipe(0);
  _acceleration.resizeAndWipe(0);
  _positionOld.resizeAndWipe(0);
  _velocityOld.resizeAndWipe(0);
  _accelerationOld.resizeAndWipe(0);
  _residual.resizeAndWipe(0);
  _temp.resizeAndWipe(0);
  _externalForces.resizeAndWipe(0);
  _gravity.resizeAndWipe(0);
  _Us.resize(0);
  _Vs.resize(0);
  _Fhats.resize(0);
  _stiffnesses.resize(0);

  _CDamping.resizeAndWipe(0,0);
  _A.resizeAndWipe(0,0);
}

//////////////////////////////////////////////////////////////////////
// change gravity and recompute
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::changeGravity(VEC3F gravityDown, Real gravityMagnitude)
{
  _gravityDown = gravityDown;
  _gravityMagnitude = gravityMagnitude;

  // compute the gravity vector
  int nodes = _tetMesh->unconstrainedNodes();
  VECTOR fullGravity(nodes * 3);
  VEC3F gravityVector = _gravityMagnitude * _gravityDown;
  for (int x = 0; x < nodes; x++)
  {
    fullGravity(3 * x) = gravityVector[0];
    fullGravity(3 * x + 1) = gravityVector[1];
    fullGravity(3 * x + 2) = gravityVector[2];
  }
  fullGravity = _tetMesh->massMatrix() * fullGravity;
  if (_tetMesh->U().rows() > 0)
    _gravity = _tetMesh->U() ^ fullGravity;
}

//////////////////////////////////////////////////////////////////////
// write out the state
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::writeState(string filename)
{
  cout << " Dumping state to file " << filename.c_str() << endl;
  FILE* file = fopen(filename.c_str(), "wb");

  _gravity.write(file);
  _position.write(file);
  _velocity.write(file);
  _acceleration.write(file);
  _positionOld.write(file);
  _velocityOld.write(file);
  _accelerationOld.write(file);
  _residual.write(file);
  _temp.write(file);
  _CDamping.write(file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read out the state
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::readState(string filename)
{
  cout << " Reading state from file " << filename.c_str() << endl;
  FILE* file = fopen(filename.c_str(), "rb");

  _gravity.read(file);
  _position.read(file);
  _velocity.read(file);
  _acceleration.read(file);
  _positionOld.read(file);
  _velocityOld.read(file);
  _accelerationOld.read(file);
  _residual.read(file);
  _temp.read(file);
  _CDamping.read(file);

  // resize external forces to the right length
  _externalForces.resizeAndWipe(_position.size());

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// reset everything to zero
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::reset()
{
  _externalForces *= 0.0;
  _positionOld *= 0.0;
  _velocityOld *= 0.0;
  _accelerationOld *= 0.0;
  _A *= 0.0;
  _position *= 0.0;
  _velocity *= 0.0;
  _acceleration *= 0.0;
  _temp *= 0.0;
  _residual *= 0.0;
}

//////////////////////////////////////////////////////////////////////
// compute velocity using the current _tetMesh q
//////////////////////////////////////////////////////////////////////
VECTOR SUBSPACE_INTEGRATOR::computeNewestReducedVelocity()
{
  VECTOR reducedVelocity(_rank);

  // get the most recent displacement
  VECTOR& qNewest = _tetMesh->q();

  reducedVelocity.clearingAxpy(_alpha[5], _accelerationOld);
  reducedVelocity.axpy(_alpha[4], _velocityOld);
  reducedVelocity.axpy(-_alpha[3], _positionOld);
  reducedVelocity.axpy(_alpha[3], qNewest);

  return reducedVelocity;
}

//////////////////////////////////////////////////////////////////////
// compute and return full coordinate velocity
//////////////////////////////////////////////////////////////////////
VECTOR SUBSPACE_INTEGRATOR::fullCoordinateVelocity()
{
  VECTOR reducedVelocity(_rank);
  MATRIX& U = _tetMesh->U();

  // get the most recent displacement
  VECTOR& qNewest = _tetMesh->q();

  reducedVelocity.clearingAxpy(_alpha[5], _accelerationOld);
  reducedVelocity.axpy(_alpha[4], _velocityOld);
  reducedVelocity.axpy(-_alpha[3], _positionOld);
  reducedVelocity.axpy(_alpha[3], qNewest);

  return U * reducedVelocity;
}

//////////////////////////////////////////////////////////////////////
// Checks the cached deformation state of the given key tet
// against the current deformation.  If the state has changed
// by a given threshold, then mark the tet as dirty, and
// update its deformation state.
//////////////////////////////////////////////////////////////////////
bool SUBSPACE_INTEGRATOR::updateTetDeformation(
                                        int i, Real threshold,
                                        MATRIX3 &F,
                                        MATRIX3 &U, MATRIX3 &Fhat, MATRIX3 &V )
{
  bool           updateDeformation = false;
  Real           oldDeformation, newDeformation;
  Real           diff;

  SUBSPACE_TET_MESH::DeformationDiagonalization &tetState
    = _tetMesh->cachedDeformation( i );

  // FIXME: This is probably the simplest thing we can do
  // here, though it probably makes sense to do something
  // that is more material aware (eg. look at changes in
  // strain energy, or maybe even internal forces, since
  // we need to get these anyways).
  oldDeformation = det( tetState._F );
  newDeformation = det( F );

#if 0
  printf( "threshold = %f\n", threshold );
  printf( "oldDeformation = %f\n", oldDeformation );
  printf( "newDeformation = %f\n\n", newDeformation );
#endif

  // Have to check for small values here
  if ( abs( oldDeformation ) < 1e-7 )
  {
    updateDeformation = abs( oldDeformation - newDeformation ) > threshold;
  }
  else
  {
    diff = abs( oldDeformation - newDeformation ) / abs( oldDeformation );

    updateDeformation = diff > threshold;
  }

  if ( updateDeformation )
  {
    tetState._F = F;
    tetState._U = U;
    tetState._Fhat = Fhat;
    tetState._V = V;
  }
#if 0
  else
  {
    cout << "NOT UPDATING DEFORMATION!!!" << endl;
    printf( "threshold = %f\n", threshold );
    printf( "oldDeformation = %f\n", oldDeformation );
    printf( "newDeformation = %f\n\n", newDeformation );

    cout << F << endl;
  }
#endif

  return updateDeformation;
}

//////////////////////////////////////////////////////////////////////
// add explicit collision forces
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::addExplicitCollisions(vector<pair<SURFACE*, int> >& collisions)
{
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

  // get the force for each point and add it to _externalForces
  map<pair<VEC3F*, SURFACE*>, bool>::iterator iter;
  for (iter = vertexHash.begin(); iter != vertexHash.end(); iter++)
  {
    VEC3F* vertex = iter->first.first;
    SURFACE* plane = iter->first.second;
    int vertexID = _tetMesh->vertexID(vertex);

    // retrieve the velocity for damping
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);
    VEC3F vertexVelocity = VEC3F(vertexSubBasis * _velocity);

    // retrieve the damped force
    VEC3F force = plane->force(*vertex, vertexVelocity);

    VECTOR projectedForce = vertexSubBasis ^ force.toVector();

    _externalForces += projectedForce;
  }
}

//////////////////////////////////////////////////////////////////////
// add implicit collision forces directly to the residual
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::addImplicitCollisionResiduals(vector<pair<SURFACE*, int> >& collisions)
{
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

  // get the force for each point and add it to the residual
  map<pair<VEC3F*, SURFACE*>, bool>::iterator iter;
  for (iter = vertexHash.begin(); iter != vertexHash.end(); iter++)
  {
    VEC3F* vertex = iter->first.first;
    SURFACE* surface = iter->first.second;
    int vertexID = _tetMesh->vertexID(vertex);

    // retrieve the velocity for damping
    MATRIX vertexSubBasis = _tetMesh->vertexSubBasis(vertex);
    VEC3F vertexVelocity = VEC3F(vertexSubBasis * _velocity);

    // don't use the tet mesh vertex position -- unproject from
    // the current integrator state
    VEC3F& restVertex= *(_tetMesh->restVertices(vertexID));
    VECTOR displace = vertexSubBasis * _position;
    VEC3F updatedPosition = restVertex + VEC3F(displace);

    // retrieve the damped force
    VEC3F force = surface->force(updatedPosition, vertexVelocity);

    // DEBUG
    //if (norm2(force) > 1e-6) cout << " vertex force: " << force << endl;

    VECTOR projectedForce = vertexSubBasis ^ force.toVector();

    _residual -= projectedForce;
  }
}

//////////////////////////////////////////////////////////////////////
// add implicit collision forces directly to the residual
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::addImplicitCollisionJacobians(vector<pair<SURFACE*, int> >& collisions)
{
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

  // get the force for each point and add it to the residual
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

    MATRIX springJacobian = surface->springJacobian(updatedPosition);
    MATRIX dampingJacobian = surface->dampingJacobian(updatedPosition, vertexVelocity);
    springJacobian *= _accelerationAlpha[4];
    dampingJacobian *= _accelerationAlpha[1];
    MATRIX jacobian = springJacobian - dampingJacobian;

    // DEBUG: see if the Jacobian looks correct
    //surface->verifySpringJacobian(*vertex);

    MATRIX reducedJacobian = vertexSubBasis ^ jacobian * vertexSubBasis;
    _A += reducedJacobian;
  }
}

//////////////////////////////////////////////////////////////////////
// compute the kinetic energy of the entire mesh
//////////////////////////////////////////////////////////////////////
Real SUBSPACE_INTEGRATOR::kineticEnergy()
{
  // unproject the velocity
  VECTOR velocity = _tetMesh->U() * _velocity;

  Real totalEnergy = 0.0;

  for (unsigned int x = 0; x < velocity.size() / 3; x++)
  {
    VEC3F vertexVelocity;
    vertexVelocity[0] = velocity[3 * x];
    vertexVelocity[1] = velocity[3 * x + 1];
    vertexVelocity[2] = velocity[3 * x + 2];

    Real mass = _tetMesh->mass(x);

    totalEnergy += mass * (vertexVelocity * vertexVelocity);
  }

  return totalEnergy;
}

//////////////////////////////////////////////////////////////////////
// Implicit Integration, using velocity as the primary variable
//////////////////////////////////////////////////////////////////////
int SUBSPACE_INTEGRATOR::stepImplicitVelocity()
{
  TIMER total;
  TIMER preamble;
  // get the reduced mass matrix
  MATRIX& M = _tetMesh->reducedMass();

  // get the state vector
  VECTOR& q = _tetMesh->q();

  // if there isn't a basis yet, do nothing
  if (q.size() == 0 || 
      (!_useKryslStiffness && _tetMesh->totalKeyTets() == 0))
    return 0;

  // accumulate external forces
  computeExternalForces();
  
  // copy the current values into the old values
  _positionOld = _position;
  _velocityOld =  _velocity;
  _accelerationOld = _acceleration;
  _tetMesh->qOld() = q;
  _timingBreakdown["Preamble"] += preamble.timing();

  // do Newton-Raphson
  Real eps = _solverEps;
  Real maxR = eps * 10;
  Real initialR = maxR;
  int step = 0;
  Real* alpha = _velocityAlpha;
 
  cout << "=========================" << endl;
  cout << " Velocity Timestep " << _totalSteps << endl;
  cout << "=========================" << endl;
  
  //while (step < _maxNewtonSteps && maxR > eps)
  while (step < _maxNewtonSteps && (maxR > initialR * eps && maxR > 1e-8))
  {
    // update position and acceleration as well
    _acceleration = alpha[0] * _velocity +
                    alpha[1] * _velocityOld +
                    alpha[2] * _accelerationOld;
    _position = _positionOld + alpha[3] * _velocityOld +
                               alpha[4] * _velocity + 
                               alpha[5] * _accelerationOld;

    // update tet mesh with the new q vector
    // an assignment is overly expensive -- the "=" forces an allocation
    TIMER updateMesh;
    q.equals(_position);
    _timingBreakdown["Update Mesh"] += updateMesh.timing();
    
    // compute the new deformation gradient F (needed by both R and K)
    _tetMesh->SUBSPACE_TET_MESH::generateF();

    // precompute diagonalizations --
    // this call is not timed out because timers are inside the function call
    cacheDiagonalizations();

    // get the reduced internal forces
    TIMER internalTimer;
    _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);

    VECTOR& R = _tetMesh->reducedInternalForce();
    _timingBreakdown["Internal Forces"] += internalTimer.timing();

    // get the reduced stiffness matrix -- do this even if we are 
    // doing the conjugate gradient solve, because the damping matrix
    // needs it
    _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, _timingBreakdown);
    MATRIX& K = _tetMesh->stiffness();
    
    // compute the damping matrix, but only if a new K is available
    TIMER RHSTimer;
    _CDamping.clearingAxpy(_rayleighAlpha, M);
    _CDamping.axpy(_rayleighBeta, K);

    _temp = M * _acceleration;
    _residual = _CDamping * _velocity;
    _residual += _temp;
    _residual -= R;
    _residual -= _externalForces;
    _timingBreakdown["RHS Assembly"] += RHSTimer.timing();
    maxR = _residual.norm2();
    if (step == 0) 
      initialR = maxR; 
    step++;

    //cout << " residual: " << _residual << endl;

    // use the acceleration level Newmark consts
    TIMER ATimer;
    _A = alpha[0] * M;
    _A.axpy(alpha[4], K);
    _A.axpy(1.0, _CDamping);
    _timingBreakdown["A Assembly"] += ATimer.timing();

    TIMER solverTimer;
    
    // solve with LU factorization
    //bool success = _A.factorLU();
    bool success = _A.factorCholesky();
    if (!success)
    {
      cout << " Choleky factor bombed! " << endl;
      exit(0);
    }
    cout << " residual " << step << ": " << maxR << endl;
    _A.solveCholesky(_residual);

    // update acceleration according to solution
    _velocity -= _residual;

    _timingBreakdown["Solver"] += solverTimer.timing();
  }

	// update velocity
  TIMER finalUpdatePosition;
  _tetMesh->q() = _position;
  _timingBreakdown["Final Update Position"] += finalUpdatePosition.timing();
  
  // update node positions
  TIMER finalUpdateForces;
  // pushing to main
  //_tetMesh->updateSurfaceMesh();
  _externalForces.clear();

  _timingBreakdown["Final Update Forces"] += finalUpdateForces.timing();
  _totalSteps++;
  _totalTime += total.timing();

  return step;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::generateImplicitVelocityMatrices(map<string,double>& timingBreakdown)
{
  Real* alpha = _velocityAlpha;
  MATRIX& M = _tetMesh->reducedMass();

  // compute the new deformation gradient F (needed by both R and K)
  //TIMER timerF;
  _tetMesh->SUBSPACE_TET_MESH::generateF();
  //timingBreakdown["Generate F"] += timerF.timing();

  // precompute diagonalizations --
  // this call is not timed out because timers are inside the function call
  //TIMER diagonalizeTimer;
  cacheDiagonalizations();
  //timingBreakdown["Diagonalize F"] += diagonalizeTimer.timing();

  // get the reduced internal forces
  _tetMesh->generateInternalForces(_Us, _Fhats, _Vs);

  VECTOR& R = _tetMesh->reducedInternalForce();

  // get the reduced stiffness matrix -- do this even if we are 
  // doing the conjugate gradient solve, because the damping matrix
  // needs it
  _tetMesh->generateStiffnessMatrix(_Us, _Vs, _stiffnesses, timingBreakdown);
  MATRIX& K = _tetMesh->stiffness();
  
  // compute the damping matrix, but only if a new K is available
  //TIMER RHSTimer;
  _CDamping.clearingAxpy(_rayleighAlpha, M);
  _CDamping.axpy(_rayleighBeta, K);

  _temp = M * _acceleration;
  _residual = _CDamping * _velocity;
  _residual += _temp;
  //_CDamping.gemvInplace(1.0, _velocity, _residual, 0.0);
  _residual -= R;
  _residual -= _externalForces;
  //timingBreakdown["RHS Assembly"] += RHSTimer.timing();

  // use the velocity level Newmark consts
  _A = alpha[0] * M;
  _A.axpy(alpha[4], K);
  _A.axpy(1.0, _CDamping);
}

//////////////////////////////////////////////////////////////////////
// get the velocity of a specific vertex
//////////////////////////////////////////////////////////////////////
VEC3F SUBSPACE_INTEGRATOR::vertexVelocity(int vertexID)
{
  assert(vertexID >= 0);
  if (_tetMesh->isConstrained(vertexID))
    return VEC3F();

  MATRIX subbasis = _tetMesh->vertexSubBasis(vertexID);
  VECTOR velocity = subbasis * _velocity;
  return VEC3F(velocity); 
}

//////////////////////////////////////////////////////////////////////
// get the old velocity of a specific vertex
//////////////////////////////////////////////////////////////////////
VEC3F SUBSPACE_INTEGRATOR::vertexVelocityOld(int vertexID)
{
  assert(vertexID >= 0);
  if (_tetMesh->isConstrained(vertexID))
    return VEC3F();

  MATRIX subbasis = _tetMesh->vertexSubBasis(vertexID);
  VECTOR velocity = subbasis * _velocityOld;
  return VEC3F(velocity); 
}

//////////////////////////////////////////////////////////////////////
// set the timestep size to something else (i.e. recompute the Newmark consts)
//////////////////////////////////////////////////////////////////////
void SUBSPACE_INTEGRATOR::setTimestep(Real dt)
{
  _dt = dt;

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
#define USING_BACKWARD_EULER 0
#if USING_BACKWARD_EULER
  // For fully implicit Euler, these have to be set explicitly. 
  // It is not possible to get fully implicit Euler using the Newmark consts.
  _accelerationAlpha[3] = 0.0;
  _accelerationAlpha[4] = _dt * _dt;
#else
  _accelerationAlpha[3] = _dt * _dt * (0.5 - _beta);
  _accelerationAlpha[4] = _dt * _dt * _beta;
#endif

  // alphas for velocity-level update
  // acceleration = alpha[0] * v_new + alpha[1] * v_old + alpha[2] * accel_old
  _velocityAlpha[0] = 1.0 / (_gamma * _dt);
  _velocityAlpha[1] = -1.0 / (_gamma * _dt);
  _velocityAlpha[2] = -(1.0 - _gamma) / _gamma;

  // position = position_old + alpha[3] * v_old + alpha[4] * v_new + alpha[5] * accel_old
  Real dtSqBeta = _dt * _dt * _beta;
  _velocityAlpha[3] = (_dt + dtSqBeta * _velocityAlpha[1]);
  _velocityAlpha[4] = dtSqBeta * _velocityAlpha[0];
  _velocityAlpha[5] = (_dt * _dt * 0.5) * (1.0 - 2.0 * _beta) + dtSqBeta * _velocityAlpha[2];
}
