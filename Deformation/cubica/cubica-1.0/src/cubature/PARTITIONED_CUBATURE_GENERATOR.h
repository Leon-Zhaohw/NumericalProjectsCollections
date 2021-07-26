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
// PARTITIONED_CUBATURE_GENERATOR.h: interface for the PARTITIONED_CUBATURE_GENERATOR class.
//
// The sparse eigensystem solve uses Slepc, which in turn uses PetSc.
// It may not be lightweight, but it is robust and well-documented.
//
//////////////////////////////////////////////////////////////////////

#ifndef PARTITIONED_CUBATURE_GENERATOR_H
#define PARTITIONED_CUBATURE_GENERATOR_H

#include <SETTINGS.h>
#include <PARTITIONED_SUBSPACE_TET_MESH.h>
#include <MERSENNETWISTER.h>
#include <vector>
#include <list>
#include <VECTOR.h>
#include <NNLS.h>
#include <slepceps.h>
#include <slepcsvd.h>

//////////////////////////////////////////////////////////////////////
// Generate the training set for the cubature scheme and run the
// Non-Negative Least Squares (NNLS) on it as well
//////////////////////////////////////////////////////////////////////
class PARTITIONED_CUBATURE_GENERATOR {

public:
  PARTITIONED_CUBATURE_GENERATOR(PARTITIONED_SUBSPACE_TET_MESH* tetMesh, const char* filename = NULL);
  ~PARTITIONED_CUBATURE_GENERATOR();

  void generateTrainingSamples(int totalSamples = 100, Real magnitude = 1.0f);
  void generateCubature();
  void generatePCACubature();

  // generate the cubature for a specific partition from PCA data
  void generatePCACubature(int partition);
  
  // generate the cubature for a specific partition from PCA data
  void generatePCACubatureDebug(int partition);

  // write and read training sets
  void write(const char* filename = NULL);
  void read(const char* filename = NULL);

  // write out the multiple cubatures
  void writePartitionedCubatures(const char* filename = NULL);
  void writePartitionedCubature(int partition, const char* filename = NULL);

  // compute the mesh basis using LMA
  void generateLMABasis(int totalModes = 20);

  int& maxKeyTets()        { return _maxKeyTets; };
  int& candidatesPerTry()  { return _candidatesPerTry; };
  double& errorTolerance() { return _errorTolerance; };
  void randSeed(int seed)  { _trainingTwister.seed(seed); };
  void randSeed(int partition, int seed)  { _twisters[partition].seed(seed); };
  int trainingSetSize()    { return _trainingForcesPartitioned[0].size(); };
  int& maxIterationsNNLS() { return _maxIterationsNNLS; };
 
  // distribute the subbases, making sure they are not rank deficient
  void distributeSubbases(vector<int> maxRanks);
  
  // distribute the subbases, making sure they are not rank deficient
  void distributeSubbases(MATRIX& modalDerivatives, VECTOR& modalEigenvalues, vector<int> maxRanks);

  // run PCA on all the partition subbases
  void pcaSubbases(vector<int> maxRanks);

  // run PCA on a specific subbases
  void pcaSubbasis(vector<int> maxRanks, int partition);

  // reinitialize for a different tet mesh
  void reinitialize(PARTITIONED_SUBSPACE_TET_MESH* tetMesh);

  vector<VECTOR>& trainingForcesPartitioned(int partition) { return _trainingForcesPartitioned[partition]; };
  vector<VECTOR>& trainingQsPartitioned(int partition) { return _trainingQsPartitioned[partition]; };
  vector<VECTOR>& trainingColumnsPartitioned(int partition) { return _trainingColumnsPartitioned[partition]; };
  vector<VECTOR>& trainingConstraintsPartitioned(int partition) { return _trainingConstraintsPartitioned[partition]; };
  
protected:
  // mesh to generate cubature for
  PARTITIONED_SUBSPACE_TET_MESH* _partitionedMesh;

  // random number generator
  MERSENNETWISTER* _twisters;
  MERSENNETWISTER _trainingTwister;

  // training samples for full mesh
  vector<VECTOR> _trainingForces;
  vector<VECTOR> _trainingQs;
 
  // inverse magnitudes of forces
  VECTOR _forceMagnitudes;

  string _filename;

  // magnitude that the samples are amplified by
  Real _magnitude;

  // tolerance of the iterative eigenvalue solver
  Real _tolerance;

  // PetSc solve matrix
  SPARSE_PETSC_MATRIX _A;

  // cubature generation params
  int _maxKeyTets;
  int _candidatesPerTry;
  double _errorTolerance;

  // pick a number of candidate tets
  //vector<TET*> pickCandidates(int candidatesPerTry);
  vector<TET*> pickCandidates(int partition, int candidatesPerTry);
  vector<TET*> pickCandidatesPCA(int partition, int candidatesPerTry);

  // find the best current tet from a list of candidates
  int findBestTet(int partition, vector<TET*>& candidates, VECTOR& residual);
  int findBestTetPCA(int partition, vector<TET*>& candidates, VECTOR& residual);

  // generate the "g" training column for a tet
  VECTOR getTrainingColumn(int partition, TET* tet);
  VECTOR getTrainingColumnPCA(int partition, TET* tet);

  // get the first 'rank' principal components of 'data'
  void pcaSlepc(MATRIX& data, int rank, MATRIX& components, VECTOR& values, bool shift = false);
  void pcaLapack(MATRIX& data, int rank, MATRIX& components, VECTOR& values, bool shift = false);

  // Filter out the dud tets
  void filterKeyTets(int partition);

  // retrieve an explicitly stored training column
  VECTOR retrieveTrainingColumn(int partition, TET* tet);

  ///////////////////////////////////////////////////////////////////////
  // Partitioning stuff starts here
  ///////////////////////////////////////////////////////////////////////
  
  // training samples for partitions
  vector<VECTOR>* _trainingForcesPartitioned;
  vector<VECTOR>* _trainingQsPartitioned;
  vector<VECTOR>* _trainingColumnsPartitioned;
  vector<VECTOR>* _trainingConstraintsPartitioned;
  vector<VECTOR> _forceMagnitudesPartitioned;

  // currently selected key tets
  vector<TET*>* _keyTets;

  // key tet weights
  vector<Real>* _keyWeights;
  
  // total number of partitions
  int _partitions;
 
  // maximum number of NNLS iterations 
  int _maxIterationsNNLS;

  // generate the cubature for a specific partition
  void generateCubature(int partition);
  
  // generate the training samples for a specific partition
  void generateTrainingSamples(int partition, int totalSamples, Real magnitude);

  // project out rotations of the rest pose
  void projectOutRotations(int partition);
};

#endif
