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
// ONLINE_CUBATURE_GENERATOR.h: interface for the ONLINE_CUBATURE_GENERATOR class.
//
// The sparse eigensystem solve uses Slepc, which in turn uses PetSc.
// It may not be lightweight, but it is robust and well-documented.
//
//////////////////////////////////////////////////////////////////////

#ifndef ONLINE_CUBATURE_GENERATOR_H
#define ONLINE_CUBATURE_GENERATOR_H

#include <SETTINGS.h>
#include <SUBSPACE_TET_MESH.h>
#include <MERSENNETWISTER.h>
#include <vector>
#include <list>
#include <VECTOR.h>
#include <NNLS.h>
#include <TIMER.h>
#include <SKELETON.h>

//////////////////////////////////////////////////////////////////////
// Generate the training set for the cubature scheme and run the
// Non-Negative Least Squares (NNLS) on it as well
//////////////////////////////////////////////////////////////////////
class ONLINE_CUBATURE_GENERATOR {

public:
  ONLINE_CUBATURE_GENERATOR(SUBSPACE_TET_MESH* tetMesh);
  ~ONLINE_CUBATURE_GENERATOR() {};

  void generateCubature(map<string, double>& timings);
  void updateCubature(map<string, double>& timings);

  bool bailed()            { return _bailed; };
  int& maxKeyTets()        { return _maxKeyTets; };
  int& candidatesPerTry()  { return _candidatesPerTry; };
  double& errorTolerance() { return _errorTolerance; };
  void randSeed(int seed)  { _twister.seed(seed); };
  int trainingSetSize()    { return _trainingForces.size(); };
  vector<VECTOR>& trainingForces()      { return _trainingForces; };
  vector<VECTOR>& trainingQs()          { return _trainingQs; };
  vector<VECTOR>& trainingConstraints() { return _trainingConstraints; };
  vector<VECTOR>& trainingColumns()     { return _trainingColumns; };
  vector<Real>& snapshotTimes()         { return _snapshotTimes; };
  vector<TET*>& keyTets()               { return _keyTets; };
  vector<Real>& keyWeights()            { return _keyWeights; };
  vector<Real>& kronrodWeights()        { return _kronrodWeights; };
  vector<int>& totalTetsAfterUpdates()  { return _totalTetsAfterUpdates; };
  vector<int>& totalDuds()              { return _totalDuds; };
  vector<Real>& priorWeights()          { return _priorWeights; };
  SKELETON*& skeleton()                 { return _skeleton; };
  Real& relativeError()                 { return _relativeError; };
 
  void reset();

  void writeState(string filename);
  void readState(string filename);
  
protected:
  // mesh to generate cubature for
  SUBSPACE_TET_MESH* _tetMesh;

  // random number generator
  MERSENNETWISTER _twister;

  // training samples
  vector<VECTOR> _trainingForces;
  vector<VECTOR> _trainingQs;
  vector<VECTOR> _trainingConstraints;
  vector<Real> _snapshotTimes;
  vector<VECTOR> _trainingColumns;

  // inverse magnitudes of forces
  VECTOR _forceMagnitudes;

  // magnitude that the samples are amplified by
  Real _magnitude;

  // tolerance of the iterative eigenvalue solver
  Real _tolerance;
  
  // currently selected key tets
  vector<TET*> _keyTets;

  // key tet weights
  vector<Real> _keyWeights;
  
  // key tet Gauss-Kronrod weights for error estimation
  vector<Real> _kronrodWeights;
  
  // cubature generation params
  int _maxKeyTets;
  int _candidatesPerTry;
  double _errorTolerance;

  // number of key tets after each updates
  vector<int> _totalTetsAfterUpdates;

  // number of key tets with zero weight that were discarded
  vector<int> _totalDuds;

  // weights given to the previous tet points upon update
  vector<Real> _priorWeights;

  // relative error of last training run
  Real _relativeError;
  
  // pick a number of candidate tets
  vector<TET*> pickCandidates(int candidatesPerTry);
  vector<TET*> pickCandidatesAll();

  // find the best current tet from a list of candidates
  int findBestTet(vector<TET*>& candidates, VECTOR& residual);

  // generate the "g" training column for a tet
  VECTOR getTrainingColumn(TET* tet);

  // retrieve an explicitly stored training column
  VECTOR retrieveTrainingColumn(TET* tet);

  // filter out the dud tets and return a vector of the non-duds
  vector<int> filterKeyTets();

  // train the Gauss-Kronrod pairs
  void trainGaussKronrod(double* A, VECTOR& b, double* weightsNNLS, NNLS_SOLVER& nnls, vector<int>& validTets);

  SKELETON* _skeleton;

  // maximum number of needed NNLS iterations seen so far;
  int _nnlsMaxIter;

  // did the last training run bail?
  bool _bailed;

};

#endif
