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
// CUBATURE_GENERATOR.h: interface for the CUBATURE_GENERATOR class.
//
// The sparse eigensystem solve uses Slepc, which in turn uses PetSc.
// It may not be lightweight, but it is robust and well-documented.
//
//////////////////////////////////////////////////////////////////////

#ifndef CUBATURE_GENERATOR_H
#define CUBATURE_GENERATOR_H

#include <SETTINGS.h>
#include <SUBSPACE_TET_MESH.h>
#include <MERSENNETWISTER.h>
#include <vector>
#include <list>
#include <VECTOR.h>
#include <NNLS.h>
#include <slepceps.h>

//////////////////////////////////////////////////////////////////////
// Generate the training set for the cubature scheme and run the
// Non-Negative Least Squares (NNLS) on it as well
//////////////////////////////////////////////////////////////////////
class CUBATURE_GENERATOR {

public:
  CUBATURE_GENERATOR(SUBSPACE_TET_MESH* tetMesh, const char* filename = NULL);
  ~CUBATURE_GENERATOR() {};

  void generateTrainingSamples(int totalSamples = 100, Real magnitude = 1.0f);
  void generateCubature();

  // write and read training sets
  void write(const char* filename = NULL);
  void read(const char* filename = NULL);

  // write the generated cubature to a file
  void writeCubature(const char* filename = NULL);

  // compute the mesh basis using LMA
  void generateLMABasis(int totalModes = 20, bool finalize = true);

  int& maxKeyTets()        { return _maxKeyTets; };
  int& candidatesPerTry()  { return _candidatesPerTry; };
  double& errorTolerance() { return _errorTolerance; };
  void randSeed(int seed)  { _twister.seed(seed); };
  int trainingSetSize()    { return _trainingForces.size(); };
  int& maxIterationsNNLS()  { return _maxIterationsNNLS; };

  vector<VECTOR>& trainingQs() { return _trainingQs; };
  vector<VECTOR>& trainingForces() { return _trainingForces; };

  void pcaBasis(int maxRank);

protected:
  // mesh to generate cubature for
  SUBSPACE_TET_MESH* _tetMesh;

  // random number generator
  MERSENNETWISTER _twister;

  // training samples
  vector<VECTOR> _trainingForces;
  vector<VECTOR> _trainingQs;

  // inverse magnitudes of forces
  VECTOR _forceMagnitudes;

  string _filename;

  // magnitude that the samples are amplified by
  Real _magnitude;

  // tolerance of the iterative eigenvalue solver
  Real _tolerance;
  
  // currently selected key tets
  vector<TET*> _keyTets;

  // key tet weights
  vector<Real> _keyWeights;
  
  // cubature generation params
  int _maxKeyTets;
  int _candidatesPerTry;
  double _errorTolerance;

  int _maxIterationsNNLS;

  // pick a number of candidate tets
  vector<TET*> pickCandidates(int candidatesPerTry);

  // find the best current tet from a list of candidates
  int findBestTet(vector<TET*>& candidates, VECTOR& residual);

  // generate the "g" training column for a tet
  VECTOR getTrainingColumn(TET* tet);

  // Filter out the dud tets
  void filterKeyTets();

  // get the first 'rank' principal components of 'data'
  void pcaLapack(MATRIX& data, int rank, MATRIX& components, VECTOR& values, bool shift = false);
};

#endif
