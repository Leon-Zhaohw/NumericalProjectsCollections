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
// SPARSE_PCG_MATRIX.h: interface for the SPARSE_PCG_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SPARSE_PCG_MATRIX_H
#define SPARSE_PCG_MATRIX_H

#include <SPARSE_MATRIX.h>
#include <PRECONDITIONER.h>

//////////////////////////////////////////////////////////////////////
// A sparse matrix class that can solve using Pardiso
//////////////////////////////////////////////////////////////////////
class SPARSE_PCG_MATRIX : public SPARSE_MATRIX {

public:
  SPARSE_PCG_MATRIX(Real eps = 1e-4, int iterations = 1000);
  SPARSE_PCG_MATRIX(MATRIX& matrix, Real eps = 1e-4, int iterations = 1000);
  SPARSE_PCG_MATRIX(int rows, int cols, Real eps = 1e-4, int iterations = 1000);
  ~SPARSE_PCG_MATRIX();

  // solve an SPD system using PCG
  void solveCG(VECTOR& x, VECTOR& b);
  //void solveICCG(VECTOR& x, VECTOR& b);
  void solvePCG(VECTOR& x, VECTOR& b, PRECONDITIONER* M);

  void setSparsity(SPARSE_MATRIX& A);

  int& maxIterations() { return _iterations; };
  Real& eps() { return _eps; };
  Real meanIterations() { return (Real)_totalIterations / _totalSolves; };
  Real meanResidual()   { return _totalResidual / _totalSolves; };
  VECTOR& residual()       { return _residual; };
  VECTOR& direction()      { return _direction; };
  VECTOR& q()              { return _q; };
  map<string, double>& timingBreakdown() { return _timingBreakdown; };
  double& totalTime() { return _totalTime; };
  void add(int row, int col, Real value) {
    pair<int,int> index;
    index.first = row;
    index.second = col;
    _matrix[index] = value;
  };

protected:
  // do allocations for symmetric matrix
  void initCG();
  void initPCG();

  // has a factorization been done before?
  bool _firstSolve;

  // store the size of the matrix that the solver expects
  int _matrixSize;

  // CG variables
  VECTOR _residual;
  VECTOR _direction;
  VECTOR _q;
  int _iterations;
  Real _eps;

  // PCG variables
  VECTOR _s;
 
  // keep track of some stats
  int _totalSolves;
  int _totalIterations;
  Real _totalResidual;

  // do a matrix-vector multiply that OpenMP can handle
  void parallelMultiply();
  //vector<pair<int,int> > _pairs;
  vector<int> _rowIndices;
  vector<int> _columnIndices;
  vector<Real> _entries;
  //vector<VECTOR> _qs;
  VECTOR* _qs;
  int _totalCores;

  // precache as much as possible before going into the main PCG loop
  void initParallelMultiply();

  map<string, double> _timingBreakdown;
  double _totalTime;
};

#endif
