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

#ifndef ONLINE_SUBSPACE_INTEGRATOR_H
#define ONLINE_SUBSPACE_INTEGRATOR_H

#include <FULLSPACE_INTEGRATOR.h>
#include <SUBSPACE_INTEGRATOR.h>
#include <ONLINE_CUBATURE_GENERATOR.h>
#include <SKELETON.h>
#include <INVERTIBLE.h>

class ONLINE_SUBSPACE_INTEGRATOR {

public:
  ONLINE_SUBSPACE_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, Real dt = 1.0 / 60.0, Real gravity = 9.8, Real alpha = 0.001, Real beta = 0.001);
  ~ONLINE_SUBSPACE_INTEGRATOR();

  // integrator options
  void stepImplicit(bool invertible = false, bool render = true, Real consecutiveAmp = 1.0);
  void stepQuasistatic(bool invertible = false, bool render = true, Real consecutiveAmp = 1.0);
  void beamStepImplicit(bool invertible = false, bool render = true, Real consecutiveAmp = 1.0);
  void forceFullStep(bool invertible);

  // set integrator precision
  void setSolverEps(Real precision) { 
    _solverEps = precision;
    _integrator->solverEps() = precision; 
    _subspaceIntegrator->solverEps() = precision;
  };
  void setPCGDigits(int digits) { _integrator->PCGDigits() = digits; };
  void stepImplicitFull(bool invertible = false, bool render = true);
  
  // mouse support
  void click(VEC3F& point);
  void drag(VEC3F& point);
  void unclick();
 
  // drawing options
  void drawClickedNode();
  void drawForceVector();
  void forceDrawReduced();
  void forceDrawFull();
  void drawOneRing();
  void drawErrorMap();
  void drawErrorPoints();

  void forceRenderManDrawReduced();
  void forceRenderManDrawFull();

  // get/sets
  void setForceMultiplier(Real multiplier) {
    _subspaceIntegrator->forceMultiplier() = multiplier;
    _integrator->forceMultiplier() = multiplier;
    _forceMultiplier = multiplier;
  };

  // poke the mesh in a deterministic way for debugging
  void poke();

  // add gravity to both integrators
  void addGravity(VEC3F direction, Real magnitude);

  // flip solver between full and subspace
  void swapSolvers();

  // print out timing info
  void printTimingBreakdown();
  void printSkippedSteps();

  // accessors
  ONLINE_CUBATURE_GENERATOR* cubatureGenerator()  { return _cubatureGenerator; };
  FULLSPACE_INTEGRATOR* fullIntegrator()          { return _integrator; };
  SUBSPACE_INTEGRATOR* subspaceIntegrator()       { return _subspaceIntegrator; };
  SUBSPACE_TET_MESH* subspaceTetMesh()            { return _tetMesh; };
  int& bccRes()                                   { return _bccRes; };
  vector<VEC3F*>& fishnetPoints()                 { return _fishnetPoints; };
  bool& didFullStep()                             { return _fullStep; };
  SKELETON*& skeleton()                           { return _skeleton; };
  bool& forceFullStep()                           { return _forceFullStep; };
  Real& dt()                                      { return _dt; };
  void setTrueErrorThreshold(Real trueErrorThreshold);
  void setMaxNewtonSteps(int steps)               { _integrator->maxNewtonSteps() = steps; _subspaceIntegrator->maxNewtonSteps() = steps; };
  bool& renderBoth()                              { return _renderBoth; };
  string& logFileDirectory()                      { return _logFileDirectory; };
  Real& snapshotDiscardThreshold()                { return _snapshotDiscardThreshold; };
  bool& useKryslStiffness()                       { return _subspaceIntegrator->useKryslStiffness(); };
  list<VECTOR>& displacementHistory()             { return _displacementHistory; };
  int& downdateSize()                             { return _downdateSize; };
  bool& downdated()                               { return _downdated; };
  bool& skipped()                                 { return _skippedStep; };
  vector<Real>& skippedSteps()                    { return _skipped; };
  int& basisType()                                { return _basisType; };
  int& maxRank()                                  { return _maxRank; };
  Real& rayleighAlpha()                           { return _rayleighAlpha; };
  Real& rayleighBeta()                            { return _rayleighBeta; };
  Real& trueErrorThreshold()                      { return _trueErrorThreshold; };
  vector<VECTOR>& Qraw()                          { return _Qraw; };
  vector<VECTOR>& Rraw()                          { return _Rraw; };
  bool storeFullTrainingColumns()                 { return _storeFullTrainingColumns; };

  // read and write state
  void writeState(string filename);
  void readState(string filename);

  void write(int frameNumber);
  void read(int frameNumber);
  void addForce(VEC3F* point, VEC3F& forceVector) { _integrator->addForce(point, forceVector); };

  // blow away all the snapshots
  void resetBasis();

  // this is a hacky passthru for data
  Real _snapshotNorm2;
  Real _snapshotNormRMS;

protected:
  ///////////////////////////////////////////////////////////////////////////////////
  // These are the big control variables
  ///////////////////////////////////////////////////////////////////////////////////
  Real _snapshotDiscardThreshold;
  Real _trueErrorThreshold;
  int _projectionErrorSamples;

  SUBSPACE_INTEGRATOR* _subspaceIntegrator;
  FULLSPACE_INTEGRATOR* _integrator;
  SUBSPACE_TET_MESH* _tetMesh;
  ONLINE_CUBATURE_GENERATOR* _cubatureGenerator;

  // precision of both of the solvers
  Real _solverEps;
  
  // take a full integrator step?
  bool _fullStep;

  // number of full steps skipped
  int _skippedFullSteps;

  // BCC resolution of original grid
  int _bccRes;

  // QR factorization
  vector<VECTOR> _Qraw;
  vector<VECTOR> _Rraw;
  vector<VECTOR> _Usnapshots;
  vector<VECTOR> _forceSnapshots;
  vector<Real> _snapshotTimes;

  // constrained nodes snapshots
  vector<VECTOR> _constraintSnapshots;
  
  // full QR factorization
  MATRIX _Q;
  MATRIX _R;

  // timesteps so far
  int _timesteps;
  int _unreducedSteps;
  int _reducedSteps;
  int _consecutiveReducedSteps;

  bool _cubatureInitialized;

  // damping parameters for both integrators
  Real _rayleighAlpha;
  Real _rayleighBeta;
 
  // gravity constant
  Real _gravityMagnitude;

  // force multiplier
  Real _forceMultiplier;

  // force a full step
  bool _forceFullStep;

  // timings
  map<string, double> _timingBreakdown;
  double _totalTime;
  Real _dt;
  
  // RNG for picking point samples
  MERSENNETWISTER _twister;

  // maximum subspace rank allowed before a downdate
  int _maxRank;
 
  // history of previous displacements -- needed by downdateBasis
  list<VECTOR> _displacementHistory;

  // error map - diff between last full step and subspace step
  VECTOR _errorMap;
  VECTOR _errorColors;
  Real _maxError;

  // cumulative distribution function of the error map -- used
  // to pick which nodes to point sample
  VECTOR _errorMapCDF;
 
  // points used for point error esimation
  vector<VEC3F*> _errorPoints;

  vector<Real> _addedToBasis;
  vector<Real> _attempedAddToBasis;
  vector<Real> _syncToFull;
  vector<Real> _subspaceTries;

  vector<int> _maxAllowedHistory;
  vector<bool> _acceptableProjectionError;

  Real _lastTrainingProjectionError;
  Real _lastTrueError;
 
  // surface points for fishnet forces
  vector<VEC3F*> _fishnetPoints;

  // skinning skeleton
  SKELETON* _skeleton;

  // drop factor for snapshots
  Real _dropFactor;

  // log file directory
  string _logFileDirectory;

  // store a snapshot and add it to the QR factorization
  bool addSnapshot(VECTOR& toAdd, bool invertible, bool verbose = true, bool errorColumn = false);

  // downdate the basis
  void downdateBasis();

  // transfer full state to reduced and vice versa
  void copyFullToReduced();
  void copyReducedToFull();

  // build the full _Q and _R matrices from _Qraw and _Rraw
  void buildQRMatrices();

  // populate the training set of ONLINE_SUBSPACE_INTEGRATOR
  void populateTrainingSet();
  
  // train a new cubature
  void trainCubature();
  
  // sanity check the QR factorization
  void checkQR();

  // compare the full and the reduced solver results
  Real computeTrueError();

  // sanity check: compute the full FR error
  VECTOR trueProjectionErrorFR();
  VECTOR trueQuasistaticResidual(bool invertible);
  VECTOR estimateQuasistaticResidual(bool invertible);

  void cacheOneRingDiagonalizations(TET_MESH* oneRing);
  vector<MATRIX3> _Us;  // U from F diagonalization
  vector<MATRIX3> _Vs;  // V from F diagonalization
  vector<MATRIX3> _Fhats;  // Fhat from F diagonalization
  vector<MATRIX> _stiffnesses;  // diagonalized and clamped stiffness

  // sample the error CDF and return its corresponding vertex
  int sampleUniform(double random);
 
  // update the q history for a change of basis
  void changeHistoryBasis(MATRIX& oldBasis, MATRIX& newBasis);

  // write a Real vector to a file
  void write(vector<Real>& vec, FILE* file);

  // read a Real vector from a file
  void read(vector<Real>& vec, FILE* file);

  // variables for sparseImplicit
  Real _errorProjectionRMS;
  Real _errorIntegrationRMS;

  int _subspaceSteps;
  int _maxSubspaceSteps;
  int _downdateSize;
  bool _previousStepWasFull;

  // render both full and reduced models
  bool _renderBoth;

  // downdated basis this timestep
  bool _downdated;

  // last step was a skipped step
  bool _skippedStep;

  // basis type - 0 = position, 1 = velocity, 2 = acceleration
  int _basisType;

  // orthogonalize a with respect to a basis -- returns the 2 norm of the
  // final vector
  Real orthogonalize(VECTOR& snapshot);

  // find the state vector with the maximum information
  VECTOR maxOrthogonalizedVector();

  // stepImplicit vars
  bool _fullActive;
  int  _maxConsecutiveSubsteps;
  int  _consecutiveSubsteps;
  bool _errorColumnPresent;

  // remove the last direction from the basis
  //void removeErrorDirection(bool invertible);
  void removeDirection(bool invertible);

  Real _errorFloor;
  Real _errorCeil;
  vector<Real> _subspaceErrors;
  vector<Real> _skipped;
  vector<Real> _concaveErrors;

  // full training column support
  bool _storeFullTrainingColumns;
  vector<VECTOR> _trainingColumns;
  void storeForceDensity();
};

#endif
