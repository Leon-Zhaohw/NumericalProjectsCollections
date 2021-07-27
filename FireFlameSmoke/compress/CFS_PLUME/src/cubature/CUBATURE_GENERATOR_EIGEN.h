/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2013 Theodore Kim
*/
// CUBATURE_GENERATOR_EIGEN.h: interface for the CUBATURE_GENERATOR_EIGEN class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CUBATURE_GENERATOR_EIGEN_H
#define CUBATURE_GENERATOR_EIGEN_H

#include "EIGEN.h"

#include <SETTINGS.h>
#include <SUBSPACE_FLUID_3D_EIGEN.h>
#include <MERSENNETWISTER.h>
#include <vector>
#include <list>
#include <VECTOR.h>
#include <NNLS.h>

//////////////////////////////////////////////////////////////////////
// Generate the training set for the cubature scheme and run the
// Non-Negative Least Squares (NNLS) on it as well
//////////////////////////////////////////////////////////////////////
class CUBATURE_GENERATOR_EIGEN {

public:
  CUBATURE_GENERATOR_EIGEN(SUBSPACE_FLUID_3D_EIGEN* fluid, const char* filename = NULL);
  ~CUBATURE_GENERATOR_EIGEN() {};

  // use importance sampling of some kind to generate the cubature
  void generateImportanceSampledCubature();

  void generateImportanceSampledCubatureMemory();

  int& candidatesPerTry()  { return _candidatesPerTry; };
  double& errorTolerance() { return _errorTolerance; };
  bool& grabAllCandidates() { return _grabAllCandidates; };
  int& importanceSamples() { return _importanceSamples; };

  // write the generated cubature to a file
  void writeCubature(const string& filename);

  const SUBSPACE_FLUID_3D_EIGEN* fluid() { return _fluid; };

  vector<VectorXd>& trainingPreadvection() { return _trainingPreadvection; };
  vector<VectorXd>& trainingPostadvection() { return _trainingPostadvection; };

  void allocate_postU_T(int rows, int cols) { _postU_T.resize(rows, cols); };
  void preprocess_prediffuse(const string& prediffusePath);

  double* postU_T_data() { return _postU_T.data(); };

  void set_prediffuseU_rows(int rows) { _prediffuseU_rows = rows; };
  void set_prediffuseU_cols(int cols) { _prediffuseU_cols = cols; };
  int prediffuseU_rows() const { return _prediffuseU_rows; };
  int prediffuseU_cols() const { return _prediffuseU_cols; };

  FILE* prediffuseFile() const { return _prediffuseFile; };

protected:
  // fluid volume to generate cubature for
  SUBSPACE_FLUID_3D_EIGEN* _fluid;

  // random number generator
  MERSENNETWISTER _twister;

  // training samples - unprojected!
  vector<VectorXd> _trainingPreadvection;
  vector<VectorXd> _trainingPostadvection;

  // cubature generation params
  int _maxKeyCells;
  int _candidatesPerTry;
  double _errorTolerance;

  // tolerance of the iterative eigenvalue solver
  Real _tolerance;

  int _maxIterationsNNLS;

  // inverse magnitudes of forces
  VectorXd _postadvectionMagnitudes;

  string _filename;

  // magnitude that the samples are amplified by
  Real _magnitude;
  
  // currently selected key cells
  vector<int> _keyCells;

  // key tet weights
  vector<Real> _keyWeights;

  // use a subset of candidates, or grab all the cells as candidates?
  bool _grabAllCandidates;

  // how many importance samples to take?
  int _importanceSamples;

  // prediffusion submatrix, transpose
  MatrixXd _postU_T;

  int _prediffuseU_rows;
  int _prediffuseU_cols;

  FILE* _prediffuseFile;

  // pick a number of candidate tets
  void pickImportanceSampledCandidatesOMP(const int pass, map<int, bool> alreadyUsed, const VectorXd& residual, vector<int>& candidates, int totalCandidates);
  void pickImportanceSampledCandidatesOMPMemory(const int pass, map<int, bool> alreadyUsed, const VectorXd& residual, vector<int>& candidates, int totalCandidates);

  // "index" is assumed to be in peeled coordinates, not full grid coordinates
  VectorXd getStagedTrainingColumn(int index);

  VectorXd getStagedTrainingColumnMemory(int index);

  // call Matlab to solve the non-negative least squares problem using the Bro and Jong
  // fast modification
  VectorXd matlabFNNLS(const vector<VectorXd>& columns, const VectorXd& b, bool verbose = false);
  VectorXd matlabFNNLS(const MatrixXd& A, const VectorXd& b, bool verbose = false);
  
  // refine the key cell set
  void cullCubatureZeros(VectorXd& weights, vector<int>& candidates);
};

#endif
