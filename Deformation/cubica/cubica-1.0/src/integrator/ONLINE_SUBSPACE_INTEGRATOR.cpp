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
// ONLINE_SUBSPACE_INTEGRATOR.h: interface for the ONLINE_SUBSPACE_INTEGRATOR class.
//
//////////////////////////////////////////////////////////////////////

#include "ONLINE_SUBSPACE_INTEGRATOR.h"
#include <float.h>

#ifdef _WIN32
#define isnan _isnan
#endif

#define BEAM_EXAMPLE 1

//////////////////////////////////////////////////////////////////////
// Constructor for newmark integrator
//////////////////////////////////////////////////////////////////////
ONLINE_SUBSPACE_INTEGRATOR::ONLINE_SUBSPACE_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, Real dt, Real gravity, Real alpha, Real beta) :
  _tetMesh(tetMesh),
  _fullStep(true),
  _skippedFullSteps(0),
  _timesteps(0),
  _unreducedSteps(0),
  _reducedSteps(0),
  _consecutiveReducedSteps(0),
  _cubatureInitialized(false),
  _gravityMagnitude(gravity),
  _forceFullStep(false),
  _totalTime(0.0),
  _dt(dt),
  _twister(314159265),
  _lastTrueError(0.0),
  _skeleton(NULL),
  _errorProjectionRMS(0.0),
  _errorIntegrationRMS(0.0),
  _subspaceSteps(0),
  _maxSubspaceSteps(0),
  _previousStepWasFull(false),
  _renderBoth(true),
  _downdated(false),
  _skippedStep(false),
  _basisType(0),
  _fullActive(true),
  _maxConsecutiveSubsteps(0),
  _consecutiveSubsteps(0),
  _errorColumnPresent(false),
  _storeFullTrainingColumns(true)
  //_storeFullTrainingColumns(false)
{
  _rayleighAlpha = alpha;
  _rayleighBeta= beta;
  cout << " Damping: " << _rayleighAlpha << " " << _rayleighBeta << endl;
  
  _subspaceIntegrator = new SUBSPACE_INTEGRATOR(tetMesh, dt, _rayleighAlpha, _rayleighBeta, _gravityMagnitude);
  _integrator = new FULLSPACE_INTEGRATOR(tetMesh, dt, _rayleighAlpha, _rayleighBeta, _gravityMagnitude);
  _cubatureGenerator = new ONLINE_CUBATURE_GENERATOR(tetMesh);

  _solverEps = 1e-2;

  _subspaceIntegrator->solverEps() = _solverEps;
  _integrator->solverEps() = _solverEps;

  // The big control settings
  _trueErrorThreshold = 0.01; // position
  _projectionErrorSamples = 5;
  _snapshotDiscardThreshold = _trueErrorThreshold;
  _maxRank = 20;

  _downdateSize = _maxRank / 2;
}

ONLINE_SUBSPACE_INTEGRATOR::~ONLINE_SUBSPACE_INTEGRATOR()
{
  delete _subspaceIntegrator;
  delete _integrator;
  delete _cubatureGenerator;
}

//////////////////////////////////////////////////////////////////////
// Set the true error threshold, and couple it to the discard
// threshold for snapshots
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::setTrueErrorThreshold(Real trueErrorThreshold)
{
  _trueErrorThreshold = trueErrorThreshold;
}

//////////////////////////////////////////////////////////////////////
// Implicit integrator
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::stepImplicitFull(bool invertible, bool render)
{
  TIMER total;

  cout << "=========================================================" << endl;
  cout << " New timestep " << _timesteps << endl;
  cout << "=========================================================" << endl;
  _timesteps++;
 
  if (_skeleton)
    _skeleton->updateSkeleton();

  // compute the external forces
  _integrator->computeExternalForces();

  // step the full integrator
  cout << " Stepping full integrator ... ";
  if (invertible)
    _integrator->stepSparseImplicitInvertible(false);
  else
    _integrator->stepSparseImplicit(false);
  cout << "done. " << endl;
  if (render)
    forceRenderManDrawFull();

  // stomp the forces in both integrators when done
  // so that they don't accumulate in the integrator we didn't
  // use
  _integrator->clearForces();

  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// Remove the error direction from the basis
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::removeDirection(bool invertible)
{
  // strip the modal derivative from the basis
  int rows = _tetMesh->U().rows();
  int cols = _tetMesh->U().cols();
  MATRIX subBasis = _tetMesh->U().getSubmatrix(0, rows, 0, cols-1);

  // update integrators with new basis
  _subspaceIntegrator->updateBasis(_tetMesh->U(), subBasis);
  _tetMesh->updateBasis(subBasis);
  _tetMesh->updateCubature();
  _subspaceIntegrator->updateCubature(invertible);

  // erase all error direction info
  _Qraw.pop_back();
  _Rraw.pop_back();
  _forceSnapshots.pop_back();
  _Usnapshots.pop_back();

  if (_trainingColumns.size() > 0)
    _trainingColumns.pop_back();
}

//////////////////////////////////////////////////////////////////////
// force the next step to be unreduced
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::forceFullStep(bool invertible)
{
  if (_errorColumnPresent)
    removeDirection(invertible);
  _fullActive = true;
}

//////////////////////////////////////////////////////////////////////
// Store a force density snapshot
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::storeForceDensity()
{
  vector<TET>& tets = _tetMesh->tets();
  VECTOR densities(tets.size() * 12);

  for (unsigned int x = 0; x < tets.size(); x++)
  {
    TET& tet = tets[x];
    MATERIAL* material = _tetMesh->materials()[tet.materialIndex()];
    MATRIX3 F3 = tet.F();
    VECTOR F = TET::flattenF(F3);

    VECTOR forceDensity(9);
    Real* data = forceDensity.data();
    material->forceDensity(F.data(), data);
    MATRIX pFpu(9,12);

    // convert density to a force
    material->computePFPu(tet, pFpu);
    VECTOR densityVector = forceDensity * pFpu;

    // store the force
    data = densities.data();
    for (int y = 0; y < 12; y++)
      data[x * 12 + y] = densityVector[y];
  }

  _trainingColumns.push_back(densities);
}

//////////////////////////////////////////////////////////////////////
// Implicit integrator
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::beamStepImplicit(bool invertible, bool render, Real consecutiveAmp)
{
  TIMER total;
  static vector<Real> addedToBasis;
  static vector<Real> earlyExits;
  static int skippedSteps = 0;

  // store default results for the tracking vectors
  _skipped.push_back(0.0);
  addedToBasis.push_back(0.0);
  earlyExits.push_back(0.0);

  cout << "====================================" << endl;
  cout << " BEAM Timestep " << _timesteps << endl;
  cout << "====================================" << endl;
  //cout << " Error threshold currently: " << _trueErrorThreshold << endl;

  if (_fullActive)
  {
    cout << " Taking FULL step" << endl;
    _skipped[_timesteps] = 0.0;
  }
  else
  {
    cout << " Taking REDUCED step" << endl;
    cout << _consecutiveSubsteps << " of " << _maxConsecutiveSubsteps << " subspace steps" << endl;
    _skipped[_timesteps] = 1.0;
    skippedSteps++;
  }

  if (_fullActive)
  {
    if (invertible)
      _integrator->stepSparseImplicitInvertible();
    else
      _integrator->stepSparseImplicit();

    Real positionError = _trueErrorThreshold * 2.0;
    Real velocityError = _trueErrorThreshold * 2.0;
    if (_tetMesh->rank() > 0)
    {
      // take a subspace step
      if (invertible)
        _subspaceIntegrator->stepInvertibleImplicit();
      else
        _subspaceIntegrator->stepImplicit();

      // compute the current error
      VECTOR velocity = _tetMesh->U() * _subspaceIntegrator->velocity();
      VECTOR error = _integrator->velocity() - velocity;
      Real trueRMS = _integrator->velocity().rms();
      velocityError = error.rms();
      //cout << " Current velocity subspace error:  " << velocityError << endl;
      velocityError = error.rms() / trueRMS;
      //cout << " Current relative velocity subspace error:  " << velocityError << endl;
      _subspaceErrors.push_back(velocityError);

      VECTOR position = _tetMesh->U() * _subspaceIntegrator->position();
      error = _integrator->position() - position;
      trueRMS = _integrator->position().rms();
      positionError = error.rms();
      //cout << " Current position subspace error:  " << positionError << endl;
      //cout << " Current realtive position subspace error:  " << positionError / trueRMS << endl;

      //cout << " Current error threshold: " << _trueErrorThreshold << endl;
    }
    else
      _subspaceErrors.push_back(_errorFloor);

    // add results to the basis
    bool added = addSnapshot(_integrator->position(), invertible, true);
    if (added)
    {
      cout << " Added to basis" << endl;
      addedToBasis[_timesteps] = 2.0;

      // make sure it doesn't treat the previously computed error as
      // the real error if the basis was downdated
      if (_downdated)
        velocityError = _trueErrorThreshold * 2.0;
    }
    else
    {
      cout << " Did not add to the basis " << endl;
      addedToBasis[_timesteps] = 1.0;
    }
  
    // check if the basis is good enough
    if (velocityError < _trueErrorThreshold)
    {
      // compute the modal derivative
      VECTOR hessian = _tetMesh->hessianProduct(_integrator->velocity(), _integrator->velocity());

      // dynamic solve
      SPARSE_PETSC_MATRIX& petscSolver = _integrator->solver();
      Real mass = _tetMesh->mass(0);
      hessian *= -_dt * _dt * 1.0 / mass;

      // solve the system to get the final error mode
      petscSolver.finalize();
      VECTOR solution(hessian.size());
      bool hessianSuccess = petscSolver.solveCG(solution, hessian);
      if (!hessianSuccess)
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << " Hessian add failed! Reduced matrix could go indefinite!" << endl;
        cout << " Add more PCG iterations to Petsc!" << endl;
      }

      // orthogonalize the error mode
      MATRIX& U = _tetMesh->U();
      int rows = U.rows();
      int cols = U.cols();
      for (int x = 0; x < cols; x++)
      {
        VECTOR column(rows);
        column = U.getColumn(x);
        Real R = column * solution;
        solution.axpy(-R, column);
      }

      // how fast was the error growing last step?
      _maxConsecutiveSubsteps = consecutiveAmp * (_snapshotDiscardThreshold / _snapshotNormRMS);
      cout << " Computed max subspace steps: " << consecutiveAmp * (_snapshotDiscardThreshold / _snapshotNormRMS) << endl;
      cout << " Preamp max subspace steps: " << _snapshotDiscardThreshold / _snapshotNormRMS << endl;

      // add the modal derivative to the basis
      if (_maxConsecutiveSubsteps > 0)
      {
        cout << " Adding Error Hessian Column " << endl;
        Real discardThreshold = _snapshotDiscardThreshold;
        _snapshotDiscardThreshold = 1e-8;
        _errorColumnPresent = addSnapshot(solution, invertible, true, true);
        _fullActive = false;
        _consecutiveSubsteps = 0;
        _snapshotDiscardThreshold = discardThreshold;
      }
      else
      {
        if (_maxConsecutiveSubsteps <= 0)
          cout << " Error mode detected that state would drift out of subspace too fast" << endl;
        else
          cout << " Hessian solve bombed" << endl;
      }
    }
    copyFullToReduced();
  }
  // else subspace integrator only
  else
  {
    if (invertible)
      _subspaceIntegrator->stepInvertibleImplicit();
    else
      _subspaceIntegrator->stepImplicit();

    copyReducedToFull();

    // get the magnitude of the last vector in the basis
    Real currentError = _trueErrorThreshold * 0.5;
    if (_errorColumnPresent)
    {
      Real qError = _subspaceIntegrator->position()[_tetMesh->rank() - 1];
      Real fullRank = _tetMesh->TET_MESH::rank();
      currentError = sqrt(qError * qError / fullRank);
    }

    //cout << " Current RMS error: " << currentError << endl;
    //cout << " Current error threshold: " << _trueErrorThreshold << endl;

    // this error thresholding needs work
    _consecutiveSubsteps++;
    if (currentError > _trueErrorThreshold || _consecutiveSubsteps == _maxConsecutiveSubsteps)
    {
      // if one is present, remove the error direction from the basis
      if (_errorColumnPresent)
      {
        removeDirection(invertible);
        _errorColumnPresent = false;
      }

      _fullActive = true;

      if (_consecutiveSubsteps != _maxConsecutiveSubsteps)
        earlyExits[_timesteps] = 1.0;
    }

    // make sure to record that no subspace error was computed
    _subspaceErrors.push_back(0.0);
  }

  // store the current displacement in the history --
  // fullIntegrator keeps the gold standard position, so store that for now.
  // Perhaps come back and store the reduced coordinates instead later.
  if (_tetMesh->rank() > 0)
    _displacementHistory.push_back(_integrator->position());

  // keep the size of displacement history tractable
  if (_displacementHistory.size() > (unsigned int)_downdateSize)
    _displacementHistory.pop_front();

  // clean up and step the time
  _subspaceIntegrator->clearForces();
  _integrator->clearForces();

  _timesteps++;
  _totalTime += total.timing();

  int totalAdds = 0;
  cout << " Added to basis: " << endl;
  for (int i = 0; i < _timesteps; i++)
    if (addedToBasis[i] < 0.5)
      cout << "-";
    else if (addedToBasis[i] < 1.5)
      cout << "=";
    else
    {
      cout << "+";
      totalAdds++;
    }
  cout << endl;
  cout << " Total adds: " << totalAdds << endl;

  cout << " Skipped steps: " << endl;
  for (int i = 0; i < _timesteps; i++)
    if (_skipped[i] < 0.5)
      cout << "-";
    else
      cout << "+";
  cout << endl;

  cout << "Total skipped steps: " << skippedSteps << " of " << _timesteps << " (" << 100.0 * (Real)skippedSteps / (_timesteps+1) << "%)" << endl;
  cout << " Subspace early exits: " << endl;
  Real sum = 0.0;
  for (int i = 0; i < _timesteps; i++)
    if (earlyExits[i] < 0.5)
      cout << "-";
    else
    {
      cout << "+";
      sum += 1.0;
    }
  cout << endl;
  cout << "Total subspace early exits: " << sum << endl;
  cout << "Current mesh rank: " << _tetMesh->rank() << endl;
}

//////////////////////////////////////////////////////////////////////
// Quasistatic integrator
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::stepQuasistatic(bool invertible, bool render, Real consecutiveAmp)
{
  TIMER total;
  static vector<Real> addedToBasis;

  static Real tol = -1;
  static Real positionTol = _trueErrorThreshold;
  static Real throttle = 1;
  static int newTrainingData = 0;
  _cubatureGenerator->errorTolerance() = positionTol;

  // store default results for the tracking vectors
  _skipped.push_back(1.0);
  addedToBasis.push_back(0.0);
  throttle = consecutiveAmp;

  cout << "====================================" << endl;
  cout << " Online quasistatic timestep " << _timesteps << endl;
  cout << "====================================" << endl;
  //cout << " Error threshold currently: " << _trueErrorThreshold << endl;
  //cout << " Throttle set to: " << throttle << endl;

  // always step the subspace if there is one
  Real forceResidualNorm = tol * 2.0;
  if (_tetMesh->rank() > 0)
  {
    cout << " Provisionally stepping the subspace model" << endl;
    // do the step
    if (invertible)
      _subspaceIntegrator->stepInvertibleQuasistatic();
    else
      _subspaceIntegrator->stepQuasistatic();

    TIMER residualTimer;
    // compute the force residual (this needs to be redone for quasistatic)
    VECTOR residual = estimateQuasistaticResidual(invertible);
    forceResidualNorm = residual.rms();
    _timingBreakdown["Sparse residual error"] += residualTimer.timing();
  }
  //cout << " force residual: " << forceResidualNorm << endl;
  //cout << " tolerance: " << tol << endl;

  // if the error is unacceptable or no tol has ever been set
  if (forceResidualNorm > throttle * tol || tol < 0.0)
  {
    VECTOR position = _integrator->position();
    VECTOR positionOld = _integrator->positionOld();

    _skipped.back() = 0.0;

    // take a full step
    if (invertible)
      _integrator->stepSparseQuasistaticInvertible();
    else
      _integrator->stepSparseQuasistatic();

    // try to add to the basis
    bool added = addSnapshot(_integrator->position(), invertible, true);

    if (added)
    {
      addedToBasis.back() = 2.0;
      newTrainingData = 0;
    }
    else
      addedToBasis.back() = 1.0;

    // copy the old state into the subspace
    MATRIX& U = _tetMesh->U();
    _subspaceIntegrator->position() = U ^ position;
    _subspaceIntegrator->positionOld() = U ^ positionOld;

    // step the subspace integrator
    if (invertible)
      _subspaceIntegrator->stepInvertibleQuasistatic();
    else
      _subspaceIntegrator->stepQuasistatic();

    // compute the positional error
    TIMER errorCheck;
    VECTOR reducedPosition = U * _subspaceIntegrator->position();
    VECTOR diff = _integrator->position() - reducedPosition;
    Real errorRMS = diff.rms();

    //cout << " Positional error: " << errorRMS << endl;
    //cout << " Positional tolerance: " << positionTol << endl;

    // if the error is good
    if (errorRMS < positionTol)
    {
      // compute the force residual
      VECTOR residual = trueQuasistaticResidual(invertible); // TODO: SPARSIFY
      Real residualRMS = residual.rms();
      cout << " Error is good enough. " << endl;
      cout << " Current residual RMS: " << residualRMS << endl;

      if (residualRMS > tol)
        tol = residualRMS;

      cout << " Residual tolerance: " << tol << endl;
    }
    _timingBreakdown["Error check"] += errorCheck.timing();
  }
  else
  {
    TIMER copyTimer;
    copyReducedToFull();
    _skippedFullSteps++;
    _timingBreakdown["Copy reduced to full"] += copyTimer.timing();
  }

  // store the current displacement in the history --
  // fullIntegrator keeps the gold standard position, so store that for now.
  // Perhaps come back and store the reduced coordinates instead later.
  TIMER historyTimer;
  if (_tetMesh->rank() > 0)
    _displacementHistory.push_back(_integrator->position());

  // keep the size of displacement history tractable
  if (_displacementHistory.size() > (unsigned int)_downdateSize)
    _displacementHistory.pop_front();
  _timingBreakdown["Tracking history"] += historyTimer.timing();

  // clean up and step the time
  _subspaceIntegrator->clearForces();
  _integrator->clearForces();

  _timesteps++;
  _totalTime += total.timing();

  int totalAdds = 0;
  cout << " Added to basis: " << endl;
  for (int i = 0; i < _timesteps; i++)
    if (addedToBasis[i] < 0.5)
      cout << "-";
    else if (addedToBasis[i] < 1.5)
      cout << "=";
    else
    {
      cout << "+";
      totalAdds++;
    }
  cout << endl;
  cout << " Total adds: " << totalAdds << endl;

  cout << " Skipped steps: " << endl;
  for (int i = 0; i < _timesteps; i++)
    if (_skipped[i] < 0.5)
      cout << "-";
    else
      cout << "+";
  cout << endl;

  cout << "Total skipped steps: " << _skippedFullSteps << " of " << _timesteps << " (" << 100.0 * (Real)_skippedFullSteps/ (_timesteps+1) << "%)" << endl;
  cout << "Current mesh rank: " << _tetMesh->rank() << endl;
}

//////////////////////////////////////////////////////////////////////
// Implicit integrator
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::stepImplicit(bool invertible, bool render, Real consecutiveAmp)
{
  TIMER total;
  static vector<Real> addedToBasis;
  static vector<Real> earlyExits;
  static int skippedSteps = 0;

  // store default results for the tracking vectors
  _skipped.push_back(0.0);
  addedToBasis.push_back(0.0);
  earlyExits.push_back(0.0);

  cout << "====================================" << endl;
  cout << " Timestep " << _timesteps << endl;
  cout << "====================================" << endl;
  cout << " Error threshold currently: " << _trueErrorThreshold << endl;

  if (_fullActive)
  {
    cout << " Taking FULL step" << endl;
    _skipped[_timesteps] = 0.0;
  }
  else
  {
    cout << " Taking REDUCED step" << endl;
    cout << _consecutiveSubsteps << " of " << _maxConsecutiveSubsteps << " subspace steps" << endl;
    _skipped[_timesteps] = 1.0;
    skippedSteps++;
  }

  if (_fullActive)
  {
    if (invertible)
      _integrator->stepSparseImplicitInvertible();
    else
      _integrator->stepSparseImplicit();

    Real currentError = _trueErrorThreshold * 2.0;
    if (_tetMesh->rank() > 0)
    {
      // take a subspace step
      if (invertible)
        _subspaceIntegrator->stepInvertibleImplicit();
      else
        _subspaceIntegrator->stepImplicit();

      // compute the current error
      VECTOR position = _tetMesh->U() * _subspaceIntegrator->position();
      VECTOR error = _integrator->position() - position;
      currentError = error.rms();
      cout << " Current subspace error:  " << currentError << endl;
      cout << " Current error threshold: " << _trueErrorThreshold << endl;
      _subspaceErrors.push_back(currentError);

      VECTOR velocity = _tetMesh->U() * _subspaceIntegrator->velocity();
      VECTOR velocityDiff = _integrator->velocity() - velocity;
      Real velocityError = velocityDiff.rms();
      cout << " Current velocity RMS: " << _integrator->velocity().rms() << endl;
      cout << " Current velocity error: " << velocityError << endl;
      cout << " Current relative velocity error: " << velocityError / _integrator->velocity().rms() << endl;

      currentError = velocityError / _integrator->velocity().rms();
      cout << " Current relative position error: " << currentError << endl;
    }
    else
      _subspaceErrors.push_back(_errorFloor);

    // add results to the basis
    cout << " Trying a position basis! " << endl;
    bool added = addSnapshot(_integrator->position(), invertible, true);
    
    if (added)
    {
      cout << " Added to basis" << endl;
      addedToBasis[_timesteps] = 2.0;

      // make sure it doesn't treat the previously computed error as
      // the real error if the basis was downdated
      if (_downdated)
        currentError = _trueErrorThreshold * 2.0;
    }
    else
    {
      cout << " Did not add to the basis " << endl;
      addedToBasis[_timesteps] = 1.0;
    }
  
    // check if the basis is good enough
    if (currentError < _trueErrorThreshold)
    {
      // compute the modal derivative
      VECTOR hessian = _tetMesh->hessianProduct(_integrator->velocity(), _integrator->velocity());

      // dynamic solve
      SPARSE_PETSC_MATRIX& petscSolver = _integrator->solver();
      Real mass = _tetMesh->mass(0);
      hessian *= -_dt * _dt * 1.0 / mass;

      // solve the system to get the final error mode
      petscSolver.finalize();
      VECTOR solution(hessian.size());
      petscSolver.solveCG(solution, hessian);

      // orthogonalize the error mode
      MATRIX& U = _tetMesh->U();
      int rows = U.rows();
      int cols = U.cols();
      for (int x = 0; x < cols; x++)
      {
        VECTOR column(rows);
        column = U.getColumn(x);
        Real R = column * solution;
        solution.axpy(-R, column);
      }

      // how fast was the error growing last step?
      _maxConsecutiveSubsteps = consecutiveAmp * (_snapshotDiscardThreshold / _snapshotNormRMS);
      cout << " Computed max subspace steps: " << consecutiveAmp * (_snapshotDiscardThreshold / _snapshotNormRMS) << endl;
      cout << " Preamp max subspace steps: " << _snapshotDiscardThreshold / _snapshotNormRMS << endl;

      // add the modal derivative to the basis
      if (_maxConsecutiveSubsteps > 0)
      {
        cout << " Adding Error Hessian Column " << endl;
        Real discardThreshold = _snapshotDiscardThreshold;
        _snapshotDiscardThreshold = 1e-8;
        _errorColumnPresent = addSnapshot(solution, invertible, true, true);
        _fullActive = false;
        _consecutiveSubsteps = 0;
        _snapshotDiscardThreshold = discardThreshold;
      }
      else
        cout << " Error mode detected that state would drift out of subspace too fast" << endl;
    }
    copyFullToReduced();
  }
  // else subspace integrator only
  else
  {
    if (invertible)
      _subspaceIntegrator->stepInvertibleImplicit();
    else
      _subspaceIntegrator->stepImplicit();

    copyReducedToFull();

    // get the magnitude of the last vector in the basis
    Real currentError = _trueErrorThreshold * 0.5;
    if (_errorColumnPresent)
    {
      Real qError = _subspaceIntegrator->position()[_tetMesh->rank() - 1];
      Real fullRank = _tetMesh->TET_MESH::rank();
      currentError = sqrt(qError * qError / fullRank);
    }

    cout << " Current RMS error: " << currentError << endl;
    cout << " Current error threshold: " << _trueErrorThreshold << endl;

    // this error thresholding needs work
    _consecutiveSubsteps++;
    if (currentError > _trueErrorThreshold || _consecutiveSubsteps == _maxConsecutiveSubsteps)
    {
      // if one is present, remove the error direction from the basis
      if (_errorColumnPresent)
      {
        removeDirection(invertible);
        _errorColumnPresent = false;
      }

      _fullActive = true;

      if (_consecutiveSubsteps != _maxConsecutiveSubsteps)
        earlyExits[_timesteps] = 1.0;
    }

    // make sure to record that no subspace error was computed
    _subspaceErrors.push_back(0.0);
  }

  // store the current displacement in the history --
  // fullIntegrator keeps the gold standard position, so store that for now.
  // Perhaps come back and store the reduced coordinates instead later.
  if (_tetMesh->rank() > 0)
    _displacementHistory.push_back(_integrator->position());

  // keep the size of displacement history tractable
  if (_displacementHistory.size() > (unsigned int)_downdateSize)
    _displacementHistory.pop_front();

  // clean up and step the time
  _subspaceIntegrator->clearForces();
  _integrator->clearForces();

  _timesteps++;
  _totalTime += total.timing();

  int totalAdds = 0;
  cout << " Added to basis: " << endl;
  for (int i = 0; i < _timesteps; i++)
    if (addedToBasis[i] < 0.5)
      cout << "-";
    else if (addedToBasis[i] < 1.5)
      cout << "=";
    else
    {
      cout << "+";
      totalAdds++;
    }
  cout << endl;
  cout << " Total adds: " << totalAdds << endl;

  cout << " Skipped steps: " << endl;
  for (int i = 0; i < _timesteps; i++)
    if (_skipped[i] < 0.5)
      cout << "-";
    else
      cout << "+";
  cout << endl;

  cout << "Total skipped steps: " << skippedSteps << " of " << _timesteps << " (" << 100.0 * (Real)skippedSteps / (_timesteps+1) << "%)" << endl;
  cout << " Subspace early exits: " << endl;
  Real sum = 0.0;
  for (int i = 0; i < _timesteps; i++)
    if (earlyExits[i] < 0.5)
      cout << "-";
    else
    {
      cout << "+";
      sum += 1.0;
    }
  cout << endl;
  cout << "Total subspace early exits: " << sum << endl;
  cout << "Current mesh rank: " << _tetMesh->rank() << endl;
}

//////////////////////////////////////////////////////////////////////
// Blow away all the snapshots
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::resetBasis()
{
  _Qraw.resize(0);
  _Usnapshots.resize(0);
  _forceSnapshots.resize(0);
  _snapshotTimes.resize(0);
  _Rraw.resize(0);
  _Q.resizeAndWipe(0,0);
  _R.resizeAndWipe(0,0);

  _errorMap.resizeAndWipe(0);
  _errorMapCDF.resizeAndWipe(0);
  _errorColors.resizeAndWipe(0);

  /*
  _errorBoundLS.clear();
  _errorBoundFR.clear();
  _trueError.clear();
  */

  _tetMesh->resetBasis();
  _subspaceIntegrator->resetBasis();
  _cubatureGenerator->reset();
}

//////////////////////////////////////////////////////////////////////
// Add a snapshot to the QR factorization
//////////////////////////////////////////////////////////////////////
bool ONLINE_SUBSPACE_INTEGRATOR::addSnapshot(VECTOR& toAdd, bool invertible, bool verbose, bool errorColumn)
{
  // if the basis has grown too large, cut it in half
  if (_tetMesh->rank() >= _maxRank)
  {
    cout << " DOWNDATING BASIS ...." << endl;
    downdateBasis();

    // need to call here in order to compute the new damping matrix,
    // in case the snapshot is then discarded
    _tetMesh->updateCubature();
    _subspaceIntegrator->updateCubature(invertible);
    _downdated = true;
  }
  else
    _downdated = false;
 
  VECTOR snapshot(toAdd);

  // NEW: store the original norm of the snapshot for relative calcuation later
  Real originalNormRMS = snapshot.rms();

  // Update the raw data QR factorization
  TIMER qrTimer;
  VECTOR R(_Qraw.size() + 1);
  if (verbose) cout << " Q raw size: " << _Qraw.size() << endl;
  for (unsigned int x = 0; x < _Qraw.size(); x++)
  {
    // store the dot product
    R[x] = _Qraw[x] * snapshot;

    // project out this component
    snapshot.axpy(-R[x], _Qraw[x]);
  }

  // store the remaining magnitude and normalize
  Real norm2 = snapshot.norm2();
  Real normRMS = snapshot.rms();
  R[_Qraw.size()] = norm2;
  if (verbose) cout << " Snapshot RMS (R): " << normRMS << endl;
  if (verbose) cout << " Snapshot 2 (R): " << norm2 << endl;
  if (verbose) cout << " Discard threshold: " << _snapshotDiscardThreshold << endl;

  // hacky pass-thru
  _snapshotNormRMS = normRMS;
  _snapshotNorm2   = norm2;

  // check for divide by zero
  if (norm2 > 1e-8)
    snapshot *= 1.0 / norm2;
  // else this vector is really, really small, so make sure
  // it is discarded
  else
  {
    cout << " Snapshot is too small to normalize" << endl;
    normRMS = _snapshotDiscardThreshold / 10.0;
    norm2   = _snapshotDiscardThreshold / 10.0;
  }

  static int discarded = 0;
  bool used = true;

  // Is this snapshot too inconsequential to add?
  if (normRMS / originalNormRMS > _snapshotDiscardThreshold)
  {
    if (verbose) cout << " Adding snapshot. " << endl;
    if (verbose && _storeFullTrainingColumns) cout << " Using explicit training columns" << endl;
    TIMER snapshotTimer;

    // record the QR factorization
    _Qraw.push_back(snapshot);
    _Rraw.push_back(R);

    // backup the current mesh state
    VECTOR& backup = _tetMesh->x();
    if (!errorColumn)
    {
      // displace the mesh
      _tetMesh->x() = toAdd;

      // record the displacements
      _Usnapshots.push_back(_integrator->position());
    }
    else
    {
      // generate the error pose
      VECTOR newAdd = _integrator->position() + toAdd;
      _tetMesh->x() = newAdd;

      // record the error pose
      _Usnapshots.push_back(newAdd);
    }
    _tetMesh->TET_MESH::updateFullMesh();
    VECTOR& internal = _integrator->computeInternalForces(invertible);

    if (_storeFullTrainingColumns)
      storeForceDensity();

    _tetMesh->x() = backup;
    _tetMesh->TET_MESH::updateFullMesh();
    _forceSnapshots.push_back(internal);
    _snapshotTimes.push_back(_timesteps * _dt);

    _timingBreakdown["Storing Snapshots"] += snapshotTimer.timing();

    // record the displacments of the constrained nodes
    _constraintSnapshots.push_back(_tetMesh->constraintDisplacements());

    // build the new basis
    TIMER finalizeQR;
    buildQRMatrices();

    // update the basis in SUBSPACE_INTEGRATOR and SUBSPACE_TET_MESH
    _subspaceIntegrator->updateBasis(_tetMesh->U(), _Q);
    _tetMesh->updateBasis(_Q);
    _timingBreakdown["QR Finalize"] += finalizeQR.timing();

    // if this is the first basis vector, set q explicitly
    if (_tetMesh->rank() == 1)
    {
      _tetMesh->q() = _tetMesh->U() ^ _integrator->position();
      _tetMesh->qOld() = _tetMesh->q();
    }

    // if we're not doing Krysl-style stiffness matrices
    if (!_subspaceIntegrator->useKryslStiffness())
    {
      // This is instrumented separately
      trainCubature();

      // update the cubature-related quantities
      TIMER reinitTimer;
      _tetMesh->updateCubature();
      _subspaceIntegrator->updateCubature(invertible);

      _timingBreakdown["Reinit subspace"] += reinitTimer.timing();
    }
    // else still update cubature because it updates the damping matrix
    else
    {
      _tetMesh->updateCubature();
      _subspaceIntegrator->updateCubature(invertible);
    }
  }
  else
  {
    if (verbose) cout << " Discarding snapshot. " << endl;
    discarded++;
    used = false;
  }

  if (verbose)
  {
    cout << (Real)discarded / _timesteps * 100.0 << "\% discarded" << endl;
    cout << " rank is " << _Qraw.size() << endl;
  }

  return used;
}

//////////////////////////////////////////////////////////////////////
// Sanity check the QR factorization
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::checkQR()
{
  buildQRMatrices();

  MATRIX sanity = _Q * _R;
  cout << " sanity: " << sanity << endl;

  cout << " snapshots: " << endl;
  for (unsigned int x = 0; x < _Usnapshots.size(); x++)
    cout << _Usnapshots[x] << endl;
}

//////////////////////////////////////////////////////////////////////
// Build the full Q and R matrices
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::buildQRMatrices()
{
  if (_Qraw.size() == 0) return;
 
  // build Q
  _Q.resizeAndWipe(_Qraw[0].size(), _Qraw.size());
  for (unsigned int col = 0; col < _Qraw.size(); col++)
    for (int row = 0; row < _Qraw[0].size(); row++)
      _Q(row, col) = _Qraw[col](row);

  /*
  // build R
  _R.resizeAndWipe(_Rraw.size(), _Rraw.size());
  for (int col = 0; col < _Qraw.size(); col++)
    for (int row = 0; row < col + 1; row++)
      _R(row, col) = _Rraw[col][row];
      */
}

//////////////////////////////////////////////////////////////////////
// Copy full state to reduced
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::copyFullToReduced()
{
  if (_tetMesh->rank() == 0) return;
  
  MATRIX& UBasis = _tetMesh->U();

  // copy integrator variables
  VECTOR& fullPosition = _integrator->position();
  VECTOR& fullPositionOld = _integrator->positionOld();
  VECTOR& subPosition = _subspaceIntegrator->position();
  VECTOR& subPositionOld = _subspaceIntegrator->positionOld();

  VECTOR& fullVelocity = _integrator->velocity();
  VECTOR& fullVelocityOld = _integrator->velocityOld();
  VECTOR& subVelocity = _subspaceIntegrator->velocity();
  VECTOR& subVelocityOld = _subspaceIntegrator->velocityOld();

  VECTOR& fullAcceleration = _integrator->acceleration();
  VECTOR& fullAccelerationOld = _integrator->accelerationOld();
  VECTOR& subAcceleration = _subspaceIntegrator->acceleration();
  VECTOR& subAccelerationOld = _subspaceIntegrator->accelerationOld();
 
  subPosition = UBasis ^ fullPosition;
  subPositionOld = UBasis ^ fullPositionOld;
  subVelocity = UBasis ^ fullVelocity;
  subVelocityOld = UBasis ^ fullVelocityOld;
  subAcceleration = UBasis ^ fullAcceleration;
  subAccelerationOld = UBasis ^ fullAccelerationOld;

  // copy tet mesh variables
  VECTOR& q = _tetMesh->q();
  VECTOR& qOld = _tetMesh->qOld();
  q = UBasis ^ fullPosition;
  qOld = UBasis ^ fullPositionOld;
}

//////////////////////////////////////////////////////////////////////
// Copy reduced state to full
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::copyReducedToFull()
{
  if (_tetMesh->rank() == 0) return;

  MATRIX& UBasis = _tetMesh->U();

  // copy integrator variables
  VECTOR& fullPosition = _integrator->position();
  VECTOR& fullPositionOld = _integrator->positionOld();
  VECTOR& subPosition = _subspaceIntegrator->position();
  VECTOR& subPositionOld = _subspaceIntegrator->positionOld();

  VECTOR& fullVelocity = _integrator->velocity();
  VECTOR& fullVelocityOld = _integrator->velocityOld();
  VECTOR& subVelocity = _subspaceIntegrator->velocity();
  VECTOR& subVelocityOld = _subspaceIntegrator->velocityOld();

  VECTOR& fullAcceleration = _integrator->acceleration();
  VECTOR& fullAccelerationOld = _integrator->accelerationOld();
  VECTOR& subAcceleration = _subspaceIntegrator->acceleration();
  VECTOR& subAccelerationOld = _subspaceIntegrator->accelerationOld();
 
  fullPosition = UBasis * subPosition;
  fullPositionOld = UBasis * subPositionOld;
  fullVelocity = UBasis * subVelocity;
  fullVelocityOld = UBasis * subVelocityOld;
  fullAcceleration = UBasis * subAcceleration;
  fullAccelerationOld = UBasis * subAccelerationOld;

  // copy tet mesh variables
  VECTOR& q = _tetMesh->q();
  VECTOR& x = _tetMesh->x();
  x = UBasis * q;

  _tetMesh->updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
// flip solver between full and subspace
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::swapSolvers()
{
  if (_fullStep)
    copyFullToReduced();
  else
    copyReducedToFull();
    
  _fullStep = !_fullStep;
}

//////////////////////////////////////////////////////////////////////
// populate the training set of ONLINE_SUBSPACE_INTEGRATOR
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::populateTrainingSet()
{
  vector<VECTOR>& trainingForces = _cubatureGenerator->trainingForces();
  vector<VECTOR>& trainingQs = _cubatureGenerator->trainingQs();
  vector<VECTOR>& trainingConstraints = _cubatureGenerator->trainingConstraints();
  vector<Real>& snapshotTimes = _cubatureGenerator->snapshotTimes();
  vector<VECTOR>& trainingColumns = _cubatureGenerator->trainingColumns();

  // stomp any previous sets
  trainingForces.clear();
  trainingQs.clear();
  trainingConstraints.clear();
  trainingColumns.clear();
  snapshotTimes.clear();

  MATRIX& UBasis = _tetMesh->U();

  cout << " Training snapshots: " << _Usnapshots.size() << endl;
  for (unsigned int x = 0; x < _Usnapshots.size(); x++)
  {
    VECTOR q = UBasis ^ _Usnapshots[x];
    trainingQs.push_back(q);
  }
  for (unsigned int x = 0; x < _forceSnapshots.size(); x++)
  {
    VECTOR forces = UBasis ^ _forceSnapshots[x];
    trainingForces.push_back(forces);
  }
  for (unsigned int x = 0; x < _constraintSnapshots.size(); x++)
    trainingConstraints.push_back(_constraintSnapshots[x]);
  for (unsigned int x = 0; x < _snapshotTimes.size(); x++)
    snapshotTimes.push_back(_snapshotTimes[x]);
  for (unsigned int x = 0; x < _trainingColumns.size(); x++)
    trainingColumns.push_back(_trainingColumns[x]);
}

//////////////////////////////////////////////////////////////////////
// train a new cubature
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::trainCubature()
{
  // convert the snapshot Fs to the reduced space
  TIMER setTimer;
  populateTrainingSet();
  _timingBreakdown["Training Set Generation"] += setTimer.timing();

  // set cubature error tolerance to current error tol
  //_cubatureGenerator->errorTolerance() = _trueErrorThreshold;

  // call the cubature generator --
  // it always adds new tets unless the existing ones
  // are explicitly stomped
  //
  // this is timed internally
  _cubatureGenerator->generateCubature(_timingBreakdown);

  TIMER keyReadTimer;
  // copy the new cubature to SUBSPACE_TET_MESH
  int*& keyTets = _tetMesh->keyTets();
  VECTOR& keyWeights = _tetMesh->keyWeights();

  vector<TET*>& newKeyTets = _cubatureGenerator->keyTets();
  vector<Real>& newKeyWeights = _cubatureGenerator->keyWeights();

  // stomp the old key tets
  delete[] keyTets;
  keyWeights.clear();

  // allocate new key tets
  int size = newKeyTets.size();
  keyTets = new int[size];
  keyWeights.resizeAndWipe(size);

  // copy key tet info over
  for (int x = 0; x < size; x++)
  {
    keyTets[x] = _tetMesh->tetID(newKeyTets[x]);
    keyWeights(x) = newKeyWeights[x];
  }
  _timingBreakdown["Reading new key tets"] += keyReadTimer.timing();
 
  _cubatureInitialized = true;
}

//////////////////////////////////////////////////////////////////////
// Print out a detailed timing breakdown
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::printTimingBreakdown()
{
  // get the full integrator timings
  map<string, double>& fullTimings = _integrator->timingBreakdown();
  map<string, double>::iterator iter;
  double fullTotal = 0.0;
  for (iter = fullTimings.begin(); iter!= fullTimings.end(); iter++)
  {
    string name = string("(FULL) ") + iter->first;
    double timing = iter->second;
    _timingBreakdown[name] = timing;
    fullTotal += timing;
  }
  int fullSteps = _integrator->totalSteps();

  // get the subspace integrator timings
  map<string, double>& subspaceTimings = _subspaceIntegrator->timingBreakdown();
  double subspaceTotal = 0.0;
  for (iter = subspaceTimings.begin(); iter!= subspaceTimings.end(); iter++)
  {
    string name = string("(SUB) ") + iter->first;
    double timing = iter->second;
    _timingBreakdown[name] = timing;
    subspaceTotal += timing;
  }

  // create an inverse map so that it will sort by time
  map<double, string> inverseMap;
  map<string, double>::iterator forwardIter;
  for (forwardIter = _timingBreakdown.begin(); forwardIter != _timingBreakdown.end(); forwardIter++)
    inverseMap[forwardIter->second] = forwardIter->first;

  // print the map out backwards since it sorts from least to greatest
  cout << "TIMING BREAKDOWN: " << endl;
  cout << "===============================================================================" << endl;
  map<double,string>::reverse_iterator backwardIter;
  double totalSeen = 0.0;
  for (backwardIter = inverseMap.rbegin(); backwardIter != inverseMap.rend(); backwardIter++)
  {
    string name = (*backwardIter).second + string("                               ");
    name = name.substr(0,30);

    cout << "[" << (*backwardIter).first / _totalTime * 100.0 << "%\t]: "
         << name.c_str() << "\t" << (*backwardIter).first / _timesteps << "s / frame" << endl;
    totalSeen += (*backwardIter).first;
  }
  Real misc = (_totalTime - totalSeen) / _totalTime * 100.0;
  cout << "[" << misc << "%\t]: " << "Misc. " << endl;
  cout << "===============================================================================" << endl;
  cout << " Current FPS: " << _timesteps / _totalTime << endl;
  float fullPercentage = ((float)fullTotal / _totalTime) * 100.0;
  cout << " Full integrator " << fullTotal << "s (" << fullPercentage << "%) "<< endl;
  cout << " Subspace integrator " << _totalTime - fullTotal << "s (" << 100 - fullPercentage << "%)" << endl;
  cout << " Total running time: " << _totalTime << endl;
  cout << " Estimated time savings: " << (fullTotal / fullSteps) * _skippedFullSteps << " seconds" << endl;
  cout << " Skipped steps: " << _skippedFullSteps << " (" << ((float)_skippedFullSteps / _timesteps) * 100.0 << "%)" << endl;
  cout << "===============================================================================" << endl;

  if (misc < 0.0)
  {
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << " BREAKDOWN ADDS UP TO MORE THAN 100! TIMERS ARE OVERLAPPING! " << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::printSkippedSteps()
{
  // show when the full integrator fired
  cout << " Sync to full: " << endl;
  for (unsigned int x = 0; x < _syncToFull.size(); x++)
    if (_syncToFull[x] < 0.5)
      cout << "-";
    else
      cout << "+";
  cout << endl;
  cout << " Added to basis: " << endl;
  for (unsigned int x = 0; x < _addedToBasis.size(); x++)
    if (_addedToBasis[x] < 0.5)
      cout << "-";
    else
      cout << "+";
  cout << endl;

  cout << " Attempted a subspace step: " << endl;
  for (unsigned int x = 0; x < _subspaceTries.size(); x++)
    if (_subspaceTries[x] > 0)
      cout << "o";
    else if (_subspaceTries[x] < 0)
      cout << "-";
    else
      cout << "x";
  cout << endl;
  
  cout << " Max allowed steps: " << endl;
  for (unsigned int x = 0; x < _maxAllowedHistory.size(); x++)
    cout << _maxAllowedHistory[x] << " ";
  cout << endl;

  cout << " Acceptable projection error?: " << endl;
  for (unsigned int x = 0; x < _acceptableProjectionError.size(); x++)
    if (_acceptableProjectionError[x])
      cout << "+";
    else
      cout << "-";
  cout << endl;
}

//////////////////////////////////////////////////////////////////////
// Compare the full to the reduced results
//////////////////////////////////////////////////////////////////////
Real ONLINE_SUBSPACE_INTEGRATOR::computeTrueError()
{
  if (_tetMesh->rank() == 0) 
    return 0.0;
  
  // unproject reduced position
  MATRIX& UBasis = _tetMesh->U();
  VECTOR reduced = UBasis * _tetMesh->q();
  VECTOR& full = _tetMesh->x();

  // unproject reduced acceleration
  VECTOR reducedAccel = UBasis * _subspaceIntegrator->acceleration();
  VECTOR& fullAccel = _integrator->acceleration();
  VECTOR accelError = fullAccel - reducedAccel;
 
  // unproject reduced velocity
  VECTOR reducedVel = UBasis * _subspaceIntegrator->velocity();
  VECTOR& fullVel = _integrator->velocity();
  VECTOR velError = fullVel - reducedVel;

  // compute the error map
  _errorMap.resizeAndWipe(full.size());

  // use velocity
  _errorMap = reducedVel - fullVel;
  return _errorMap.normInf() * _dt;
}

//////////////////////////////////////////////////////////////////////
// Draw the one-ring of a vertex
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::drawOneRing()
{
  vector<VEC3F>& vertices = _tetMesh->vertices();
  //int index = vertices.size() * 0.2;
  int index = 1361;

  VEC3F* picked = &(vertices[index]);
  vector<TET*> tetRing;
  _tetMesh->oneRing(picked, tetRing);

  for (unsigned int x = 0; x < tetRing.size(); x++)
    //tetRing[x]->drawLines();
    tetRing[x]->drawFaces();
  
  glColor4f(1.0f, 0.0f, 0.0f, 10.0f);
  glBegin(GL_POINTS);
#ifdef SINGLE_PRECISION      
    glVertex3fv(*picked);
#else
    glVertex3dv(*picked);
#endif
  glEnd();
}

//////////////////////////////////////////////////////////////////////
// Draw the mesh colored with the error
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::drawErrorMap()
{
  // if there is no error map, just draw the mesh as usual
  if (_errorMap.size() == 0)
  {
    _tetMesh->drawSurfaceFaces();
    return;
  }
  
  vector<pair<int,int> >& surfaceFaces = _tetMesh->surfaceFaces();
  vector<TET>& tets = _tetMesh->tets();
  for (unsigned int x = 0; x < surfaceFaces.size(); x++)
  {
    Real color = 1.0 - _errorColors[x];
    glColor4f(1.0, color, color, 1.0);

    // draw it
    TRIANGLE triangle = tets[surfaceFaces[x].first].face(surfaceFaces[x].second);
    triangle.draw();
  }
}

//////////////////////////////////////////////////////////////////////
// draw the error points
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::drawErrorPoints()
{
  glDisable(GL_LIGHTING);
  
  glDisable(GL_DEPTH_TEST);
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < _errorPoints.size(); x++)
  {
    // color the points proportional to the error
    int index = _tetMesh->vertexID(_errorPoints[x]);

    if (_errorMap.size() > 0)
    {
      VEC3F error;
      error[0] = _errorMap[3 * index];
      error[1] = _errorMap[3 * index + 1];
      error[2] = _errorMap[3 * index + 2];
      Real color = norm(error) / _maxError;

      glColor4f(color, 0.0, 0.0, 1.0);
    }
#ifdef SINGLE_PRECISION      
    glVertex3fv(*(_errorPoints[x]));
#else
    glVertex3dv(*(_errorPoints[x]));
#endif
  }
  glEnd();
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
}

//////////////////////////////////////////////////////////////////////
// sample uniformly across the mesh
//////////////////////////////////////////////////////////////////////
int ONLINE_SUBSPACE_INTEGRATOR::sampleUniform(double random)
{
  int nodes = _tetMesh->unconstrainedNodes();
  return nodes * random;
}

//////////////////////////////////////////////////////////////////////
// downdate the basis by cutting it in half
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::downdateBasis()
{
  cout << endl;
  cout << "========================================" << endl;
  cout << " Max basis size exceeded. Downdating ..." << endl;
  cout << "========================================" << endl;
    
  // assemble the q history into a matrix so that we can perform an
  // SVD on it
  MATRIX qHistory(_tetMesh->rank(), _displacementHistory.size());
  list<VECTOR>::iterator iter;
  int col = 0;
  for (iter = _displacementHistory.begin(); iter != _displacementHistory.end(); iter++, col++)
  {
    // project down to q
    VECTOR q = _tetMesh->U() ^ (*iter);

    for (int x = 0; x < _tetMesh->rank(); x++) 
      qHistory(x, col) = q[x];
  }
  
  // perform SVD on the history matrix
  MATRIX U;
  VECTOR S;
  MATRIX VT;
  qHistory.SVD(U, S, VT);

  // rotate the Q basis
  MATRIX rotated = _Q * U;

  // stomp out orthogonalization drift
  rotated.orthogonalize();

  // stomp the old raw QR
  _Qraw.resize(0);
  _Rraw.resize(0);

  // store the columns on this new basis
  for (int x = 0; x < rotated.cols(); x++)
  {
    VECTOR Qcol(rotated.rows());
    for (int y = 0; y < rotated.rows(); y++)
      Qcol[y] = rotated(y,x);
    _Qraw.push_back(Qcol);
  }

  _Usnapshots.resize(0);
  _forceSnapshots.resize(0);
  _snapshotTimes.resize(0);
  if (_storeFullTrainingColumns)
    _trainingColumns.resize(0);

  _subspaceIntegrator->updateBasis(_tetMesh->U(), rotated);
  _tetMesh->updateBasis(rotated);
}

//////////////////////////////////////////////////////////////////////
// Cache the deformation gradient diagonalizations
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::cacheOneRingDiagonalizations(TET_MESH* oneRing)
{
  int totalTets = oneRing->totalTets();
  if (_Us.size() != (unsigned int)totalTets)
  {
    _Us.resize(totalTets);
    _Vs.resize(totalTets);
    _Fhats.resize(totalTets);
    _stiffnesses.resize(totalTets);
  }

  // if the material isn't INVERTIBLE, you did something really bad
  vector<TET>& tets = oneRing->tets();
  for (int x = 0; x < totalTets; x++)
  {
    int materialIndex = tets[x].materialIndex();
    INVERTIBLE* material = (INVERTIBLE*)(oneRing->materialCopies()[0][materialIndex]);
    
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

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR ONLINE_SUBSPACE_INTEGRATOR::estimateQuasistaticResidual(bool invertible)
{
  // backup tet mesh state
  VECTOR backup = _tetMesh->TET_MESH::x();

  // set tet mesh to reduced positions
  _tetMesh->TET_MESH::x() = _tetMesh->U() * _subspaceIntegrator->position();
  _tetMesh->updateFullMesh();

  int tries = _tetMesh->unconstrainedNodes() * 0.0025;
  VECTOR allResiduals(3 * tries);
  for (int x = 0; x < tries; x++)
  {
    int pick = sampleUniform(_twister.rand());
    VEC3F* pickedPoint = _tetMesh->vertices(pick);

    // construct the vertex's one ring
    map<VEC3F*, int> originalToCopy;
    TET_MESH oneRing(pickedPoint, _tetMesh, originalToCopy, _timingBreakdown);

    // move all the vertices to their current position
    map<VEC3F*, int>::iterator iter;
    bool first = true;
    vector<Real> allDisplacements;
    allDisplacements.resize(oneRing.vertices().size() * 3);
    for (iter = originalToCopy.begin(); iter != originalToCopy.end(); iter++)
    {
      VEC3F* copy = oneRing.vertices(iter->second);
      VEC3F* restCopy = oneRing.restVertices(iter->second);

      int vertexID = iter->second;
      
      // if it is constrained, move it directly
      if (_tetMesh->isConstrained(iter->first))
      {
        *(copy) = *(iter->first);
        allDisplacements[vertexID * 3] = 0;
        allDisplacements[vertexID * 3 + 1] = 0;
        allDisplacements[vertexID * 3 + 2] = 0;
        continue;
      }
      
      // move the vertices to the current position
      MATRIX vertexBasis = _tetMesh->vertexSubBasis(iter->first);
      VECTOR newDisplacement = vertexBasis * _subspaceIntegrator->position();

      (*copy)[0] = (*restCopy)[0] + newDisplacement[0];
      (*copy)[1] = (*restCopy)[1] + newDisplacement[1];
      (*copy)[2] = (*restCopy)[2] + newDisplacement[2];
      allDisplacements[vertexID * 3] = newDisplacement[0];
      allDisplacements[vertexID * 3 + 1] = newDisplacement[1];
      allDisplacements[vertexID * 3 + 2] = newDisplacement[2];
      if (first)
      {
        oneRing.x()[0] = newDisplacement[0];
        oneRing.x()[1] = newDisplacement[1];
        oneRing.x()[2] = newDisplacement[2];
        first = false;
      }
    }

    // copy position, velocity, and acceleration into the integrator for the 
    // point of interest
    MATRIX vertexBasis = _tetMesh->vertexSubBasis(pickedPoint);
    VECTOR position    = vertexBasis * _subspaceIntegrator->position();
    VECTOR positionOld = vertexBasis * _subspaceIntegrator->positionOld();

    // grab the external forces for this node
    VECTOR externalForces(3);
    externalForces[0] = _integrator->externalForces()[3 * pick];
    externalForces[1] = _integrator->externalForces()[3 * pick + 1];
    externalForces[2] = _integrator->externalForces()[3 * pick + 2];

    if (invertible)
    {
      cacheOneRingDiagonalizations(&oneRing);
      oneRing.TET_MESH::generateInternalForces(_Us, _Fhats, _Vs);
    }
    else
      oneRing.TET_MESH::generateInternalForces();

    // generate internal forces
    VECTOR& R = oneRing.TET_MESH::internalForce();

    // generate the stiffness matrix
    if (invertible)
      oneRing.TET_MESH::generateOneRingStiffness(_Us, _Vs, _stiffnesses);
    else
      oneRing.TET_MESH::generateSparseStiffnessMatrix();
    SPARSE_MATRIX& K = oneRing.TET_MESH::stiffnessMatrix();

    // compute the balance
    VECTOR displacements(allDisplacements);
    VECTOR product = K * displacements;

    allResiduals[3 * x] = product[0] - R[0];
    allResiduals[3 * x + 1] = product[1] - R[1];
    allResiduals[3 * x + 2] = product[2] - R[2];
  }

  // restore the tet mesh state
  _tetMesh->TET_MESH::x() = backup;
  _tetMesh->TET_MESH::updateFullMesh();

  return allResiduals;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR ONLINE_SUBSPACE_INTEGRATOR::trueQuasistaticResidual(bool invertible)
{
  // backup tet mesh state
  VECTOR backup = _tetMesh->TET_MESH::x();

  // set tet mesh to reduced positions
  _tetMesh->TET_MESH::x() = _tetMesh->U() * _subspaceIntegrator->position();
  _tetMesh->updateFullMesh();

  if (invertible)
  {
    _integrator->cacheDiagonalizations();
    _tetMesh->TET_MESH::generateInternalForces(
        _integrator->_Us,
        _integrator->_Fhats,
        _integrator->_Vs);
  }
  else
    _tetMesh->TET_MESH::generateInternalForces();

  // generate internal forces
  VECTOR& R = _tetMesh->TET_MESH::internalForce();

  // generate the stiffness matrix
  if (invertible)
  {
    _tetMesh->TET_MESH::generateSparseStiffnessMatrix(
        _integrator->_Us,
        _integrator->_Vs,
        _integrator->_stiffnesses);
  }
  else
    _tetMesh->TET_MESH::generateSparseStiffnessMatrix();
  SPARSE_MATRIX& K = _tetMesh->TET_MESH::stiffnessMatrix();

  // compute the balance
  VECTOR product = K * _tetMesh->TET_MESH::x();
  VECTOR diff = product - R;

  // restore the tet mesh state
  _tetMesh->TET_MESH::x() = backup;
  _tetMesh->TET_MESH::updateFullMesh();

  return diff;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR ONLINE_SUBSPACE_INTEGRATOR::trueProjectionErrorFR()
{
  // get the mass and damping matrices
  SPARSE_MATRIX& M = _tetMesh->TET_MESH::massMatrix();
  SPARSE_MATRIX& C = _integrator->dampingMatrix();

  // unproject the velocity and acceleration
  VECTOR acceleration = _tetMesh->U() * _subspaceIntegrator->acceleration();
  VECTOR velocity     = _tetMesh->U() * _subspaceIntegrator->velocity();

  // backup tet mesh state
  VECTOR backup = _tetMesh->TET_MESH::x();

  // set tet mesh to reduced positions
  _tetMesh->TET_MESH::x() = _tetMesh->U() * _subspaceIntegrator->position();
  _tetMesh->updateFullMesh();

  // generate internal forces
  VECTOR& R = _tetMesh->TET_MESH::generateInternalForces();

  VEC3F v0;
  v0[0] = _tetMesh->TET_MESH::x()[0];
  v0[1] = _tetMesh->TET_MESH::x()[1];
  v0[2] = _tetMesh->TET_MESH::x()[2];

  // compute the final force balance
  VECTOR residual = M * acceleration;
  residual += C * velocity;
  residual -= R;
  residual -= _integrator->externalForces();

  // restore the tet mesh state
  _tetMesh->TET_MESH::x() = backup;
  _tetMesh->TET_MESH::updateFullMesh();

  //return residual.rms();
  return residual;
}

//////////////////////////////////////////////////////////////////////
// force a GL draw of the reduced mesh
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::forceDrawReduced()
{
  if (_tetMesh->q().size() == 0) return;

  // only do something special if we just did a full step -- then the
  // full and reduced tet meshes would be different
  if (_fullStep)
    _tetMesh->SUBSPACE_TET_MESH::updateSurfaceMesh();

  _tetMesh->SUBSPACE_TET_MESH::drawSurfaceFaces();

  // restore old position
  if (_fullStep)
    _tetMesh->TET_MESH::updateFullMesh();
}

//////////////////////////////////////////////////////////////////////
// force a RenderMan draw of the reduced mesh
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::forceRenderManDrawReduced()
{
#if USING_RENDERMAN
  if (_tetMesh->q().size() == 0) return;

  // backup the tet mesh state
  VECTOR backup = _tetMesh->TET_MESH::x();

  // set the tet mesh state to whatever the subspace mesh thinks its is
  _tetMesh->x() = _tetMesh->U() * _subspaceIntegrator->position();
  if (_skeleton)
    _skeleton->addSkinningDisplacements();
  else
    //_tetMesh->TET_MESH::updateSurfaceMesh();
    _tetMesh->TET_MESH::updateFullMesh();

  // update and draw
  RtColor red;
  red[0] = 1.0;
  red[1] = 0.0;
  red[2] = 0.0;
  RiColor(red);
  //_tetMesh->drawSurfaceToRenderMan();
  _tetMesh->drawExhaustiveToRenderMan();
  if (_skeleton) _skeleton->undoSkinningDisplacements();

  // restore the tet mesh state to whatever it was before
  _tetMesh->x() = backup;
  _tetMesh->TET_MESH::updateFullMesh();
#endif
}

//////////////////////////////////////////////////////////////////////
// force a RenderMan draw of the full mesh
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::forceRenderManDrawFull()
{
#if USING_RENDERMAN
  // backup the tet mesh state
  VECTOR backup = _tetMesh->TET_MESH::x();

  // set the tet mesh state to whatever the fullspace mesh thinks its is
  _tetMesh->x() = _integrator->position();
  if (_skeleton)
    _skeleton->addSkinningDisplacements();
  else
    _tetMesh->TET_MESH::updateFullMesh();

  // update and draw
  RtColor white;
  white[0] = 1.0;
  white[1] = 1.0;
  white[2] = 1.0;
  RiColor(white);
  //_tetMesh->drawSurfaceToRenderMan();
  _tetMesh->drawExhaustiveToRenderMan();
  if (_skeleton) _skeleton->undoSkinningDisplacements();

  // restore the tet mesh state to whatever it was before
  _tetMesh->x() = backup;
  _tetMesh->TET_MESH::updateFullMesh();
#endif
}

//////////////////////////////////////////////////////////////////////
// Passthrus
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::click(VEC3F& point)
{
  _subspaceIntegrator->click(point);
  _integrator->click(point);
}
void ONLINE_SUBSPACE_INTEGRATOR::drag(VEC3F& point)
{
  _subspaceIntegrator->drag(point);
  _integrator->drag(point);
}
void ONLINE_SUBSPACE_INTEGRATOR::unclick()
{
  _subspaceIntegrator->unclick();
  _integrator->unclick();
}
void ONLINE_SUBSPACE_INTEGRATOR::drawClickedNode()
{
  // delegate to unreduced
  _integrator->drawClickedNode();
}
void ONLINE_SUBSPACE_INTEGRATOR::drawForceVector()
{
  _integrator->drawForceVector();
}
void ONLINE_SUBSPACE_INTEGRATOR::poke()
{
  _integrator->poke();
  _subspaceIntegrator->poke();
}
void ONLINE_SUBSPACE_INTEGRATOR::addGravity(VEC3F direction, Real magnitude)
{
  _integrator->changeGravity(direction, magnitude);
  _subspaceIntegrator->changeGravity(direction, magnitude);

  _integrator->addGravity(direction, magnitude);
  _subspaceIntegrator->addGravity();
}

//////////////////////////////////////////////////////////////////////
// write state to a binary file
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::writeState(string filename)
{
  cout << " Dumping state to file " << filename.c_str() << endl;
  FILE* file = fopen(filename.c_str(), "wb");
  int size;

  size = _Qraw.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (unsigned int x = 0; x < _Qraw.size(); x++)
    _Qraw[x].write(file);

  size = _Usnapshots.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (unsigned int x = 0; x < _Usnapshots.size(); x++)
    _Usnapshots[x].write(file);

  size = _forceSnapshots.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (unsigned int x = 0; x < _forceSnapshots.size(); x++)
    _forceSnapshots[x].write(file);

  size = _constraintSnapshots.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (unsigned int x = 0; x < _constraintSnapshots.size(); x++)
    _constraintSnapshots[x].write(file);

  _Q.write(file);
  _R.write(file);

  fwrite((void*)&_timesteps, sizeof(int), 1, file);
  fwrite((void*)&_unreducedSteps, sizeof(int), 1, file);
  fwrite((void*)&_reducedSteps, sizeof(int), 1, file);
  fwrite((void*)&_consecutiveReducedSteps, sizeof(int), 1, file);

  size = _displacementHistory.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  list<VECTOR>::iterator iter;
  for (iter = _displacementHistory.begin(); iter != _displacementHistory.end(); iter++)
    iter->write(file);

  _errorMap.write(file);
  _errorColors.write(file);
  _errorMapCDF.write(file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Write a real vector to a file
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::write(vector<Real>& vec, FILE* file)
{
  int size = vec.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for (unsigned int x = 0; x < vec.size(); x++)
  {
    Real val = vec[x];
    fwrite((void*)&val, sizeof(Real), 1, file);
  }
}

//////////////////////////////////////////////////////////////////////
// Read a real vector from a file
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::read(vector<Real>& vec, FILE* file)
{
  vec.clear();
  int size;
  fread((void*)&size, sizeof(int), 1, file);
  for (unsigned int x = 0; x < vec.size(); x++)
  {
    Real val;
    fread((void*)&val, sizeof(Real), 1, file);
    vec.push_back(val);
  }
}

//////////////////////////////////////////////////////////////////////
// read state from a binary file
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::readState(string filename)
{
  cout << " Reading state from file " << filename.c_str() << endl;
  FILE* file = fopen(filename.c_str(), "rb");
  int size;

  fread((void*)&size, sizeof(int), 1, file);
  _Qraw.resize(size);
  for (int x = 0; x < size; x++)
    _Qraw[x].read(file);

  fread((void*)&size, sizeof(int), 1, file);
  _Usnapshots.resize(size);
  for (unsigned int x = 0; x < _Usnapshots.size(); x++)
    _Usnapshots[x].read(file);

  fread((void*)&size, sizeof(int), 1, file);
  _forceSnapshots.resize(size);
  for (unsigned int x = 0; x < _forceSnapshots.size(); x++)
    _forceSnapshots[x].read(file);

  fread((void*)&size, sizeof(int), 1, file);
  _constraintSnapshots.resize(size);
  for (unsigned int x = 0; x < _constraintSnapshots.size(); x++)
    _constraintSnapshots[x].read(file);

  _Q.read(file);
  _R.read(file);

  fread((void*)&_timesteps, sizeof(int), 1, file);
  fread((void*)&_unreducedSteps, sizeof(int), 1, file);
  fread((void*)&_reducedSteps, sizeof(int), 1, file);
  fread((void*)&_consecutiveReducedSteps, sizeof(int), 1, file);

  fread((void*)&size, sizeof(int), 1, file);
  for (int x = 0; x < size; x++)
  {
    VECTOR history;
    history.read(file);
    _displacementHistory.push_back(history);
  }

  _errorMap.read(file);
  _errorColors.read(file);
  _errorMapCDF.read(file);

  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// Read state from a file
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::read(int frameNumber)
{
  char buffer[256];
  sprintf(buffer, "%i", frameNumber);

  string number = string(buffer);
  if (frameNumber < 10) number = std::string("0") + number;
  if (frameNumber < 100) number = std::string("0") + number;
  if (frameNumber < 1000) number = std::string("0") + number;

  string fullfile = string("./data/data.");
  fullfile += number;
  fullfile += string(".full.state");
  _integrator->readState(fullfile);

  string subfile = string("./data/data.");
  subfile += number;
  subfile += string(".sub.state");
  _subspaceIntegrator->readState(subfile);

  string meshfile = string("./data/data.");
  meshfile += number;
  meshfile += string(".mesh.state");
  _tetMesh->readState(meshfile);

  string onlinefile = string("./data/data.");
  onlinefile += number;
  onlinefile += string(".online.state");
  readState(onlinefile);
}

//////////////////////////////////////////////////////////////////////
// Write state to files
//////////////////////////////////////////////////////////////////////
void ONLINE_SUBSPACE_INTEGRATOR::write(int frameNumber)
{
  char buffer[256];
  sprintf(buffer, "%i", frameNumber);

  string number = string(buffer);
  if (frameNumber < 10) number = std::string("0") + number;
  if (frameNumber < 100) number = std::string("0") + number;
  if (frameNumber < 1000) number = std::string("0") + number;

  string fullfile = string("./data/data.");
  fullfile += number;
  fullfile += string(".full.state");
  _integrator->writeState(fullfile);

  string subfile = string("./data/data.");
  subfile += number;
  subfile += string(".sub.state");
  _subspaceIntegrator->writeState(subfile);

  string meshfile = string("./data/data.");
  meshfile += number;
  meshfile += string(".mesh.state");
  _tetMesh->writeState(meshfile);

  string onlinefile = string("./data/data.");
  onlinefile += number;
  onlinefile += string(".online.state");
  writeState(onlinefile);
}

//////////////////////////////////////////////////////////////////////
// orthogonalize a with respect to a basis -- returns the 2 norm of the
// final vector
//////////////////////////////////////////////////////////////////////
Real ONLINE_SUBSPACE_INTEGRATOR::orthogonalize(VECTOR& snapshot)
{
  // Update the raw data QR factorization
  TIMER qrTimer;
  VECTOR R(_Qraw.size() + 1);
  for (unsigned int x = 0; x < _Qraw.size(); x++)
  {
    // store the dot product
    R[x] = _Qraw[x] * snapshot;

    // project out this component
    snapshot.axpy(-R[x], _Qraw[x]);
  }

  return snapshot.norm2();
}

//////////////////////////////////////////////////////////////////////
// find the state vector with the maximum information
//////////////////////////////////////////////////////////////////////
VECTOR ONLINE_SUBSPACE_INTEGRATOR::maxOrthogonalizedVector()
{
  VECTOR position(_integrator->position());
  VECTOR velocity(_integrator->velocity());
  VECTOR acceleration(_integrator->acceleration());

  Real normPosition = orthogonalize(position);
  Real normVelocity = orthogonalize(velocity);
  Real normAcceleration = orthogonalize(acceleration);

  if (normPosition > normVelocity)
    if (normPosition > normAcceleration)
    {
      cout << " Position is the max orthogonal" << endl;
      return position;
    }
    else
    {
      cout << " Acceleration is the max orthogonal" << endl;
      return acceleration;
    }
  else
    if (normVelocity > normAcceleration)
    {
      cout << " Velocity is the max orthogonal" << endl;
      return velocity;
    }

  cout << " Acceleration is the max orthogonal" << endl;
  return acceleration;
}
