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
// SPARSE_PETSC_MATRIX.h: interface for the SPARSE_PETSC_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef SPARSE_PETSC_MATRIX_H
#define SPARSE_PETSC_MATRIX_H

#include <SPARSE_MATRIX.h>
#include <BLOCK_MATRIX.h>

#if USING_PETSC
#define EXTERN extern
#include <petscksp.h>
#endif

//////////////////////////////////////////////////////////////////////
// A sparse matrix class that can solve using Petsc and Slepc
//////////////////////////////////////////////////////////////////////
class SPARSE_PETSC_MATRIX {

public:
  SPARSE_PETSC_MATRIX();
  SPARSE_PETSC_MATRIX(int rows, int cols);
  SPARSE_PETSC_MATRIX(MATRIX& A);
  ~SPARSE_PETSC_MATRIX();

  //void solvePCG(VECTOR& x, VECTOR& b);
  bool solveCG(VECTOR& x, VECTOR& b);

  void PCA(MATRIX& data, int rank, MATRIX& component, VECTOR& values);

  // perform PCA, but the matrix being passed in should already have
  // subtracted its mean!
  void meanNormalizedPCA(int rank, MATRIX& component, VECTOR& values);

  VECTOR smallestEigenvalues(int totalEigenvalues = -1) { return eigenvalues(true, totalEigenvalues); };
  VECTOR largestEigenvalues(int totalEigenvalues = -1) { return eigenvalues(false, totalEigenvalues); };
  VECTOR eigenvalues(bool smallest = true, int totalEigenvalues = -1);

  int& maxIterations() { return _iterations; };
  Real& eps() { return _eps; };

  map<string, double>& timingBreakdown() { return _timingBreakdown; };
  double& totalTime() { return _totalTime; };

  void setSparsity(SPARSE_MATRIX& A);
  void eraseRowColumn(int rowCol);

  void set(int row, int col, Real value) {
#if USING_PETSC
    // make sure it casts to double
    double val = value;
    MatSetValues(_APetsc,1,&row,1,&col,&val,INSERT_VALUES);
#endif
  };
  void add(int row, int col, Real value) {
#if USING_PETSC
    // make sure it casts to double
    double val = value;
    MatSetValues(_APetsc,1,&row,1,&col,&val,ADD_VALUES);
#endif
  };
  void clear();
  void scale(Real scalar) { 
#if USING_PETSC
    MatScale(_APetsc, scalar); 
#endif
  };
  void axpy(Real scalar, SPARSE_MATRIX& A);
  void finalize();

  int& rows() { return _rows; };
  int& cols() { return _cols; };
  const bool firstSolve() { return _firstSolve; };
  
  // assumes setSparsity has been called
  void equals(SPARSE_MATRIX& matrix);

  // assumes setSparsity has been called
  void equals(BLOCK_MATRIX& matrix);

  // assumes setSparsity has been called
  void equals(MATRIX& matrix);
#if USING_PETSC
  Mat petscData() { return _APetsc; };
#endif

  // these enums are ugly, but they get rid of the compiler warnings
  enum PRECONDITIONER { JACOBI, SOR, ICC, ILU, NONE };
  enum SOLVER { PCG, MINRES, GMRES };

  // select the preconditioner to use
  void useNone() { _whichPreconditioner = NONE; };
  void useJacobi() { _whichPreconditioner = JACOBI; };
  void useSOR()    { _whichPreconditioner = SOR; };
  void useICC()    { _whichPreconditioner = ICC; };
  void useILU()    { _whichPreconditioner = ILU; };

  // select the solver to use
  void usePCG()    { _whichSolver = PCG; };
  void useMINRES() { _whichSolver = MINRES; };
  void useGMRES()  { _whichSolver = GMRES; };

  // only initialize PetSc once
  static bool singleton;

protected:
  int _rows;
  int _cols;

  // has a factorization been done before?
  bool _firstSolve;

  // had the matrix been finalized before?
  bool _finalized;

  // store the size of the matrix that the solver expects
  int _matrixSize;

  // keep track of some stats
  int _totalSolves;
  int _totalIterations;
  Real _totalResidual;

  // error tolerances
  int _iterations;
  Real _eps;

  // timing info
  map<string, double> _timingBreakdown;
  double _totalTime;

#if USING_PETSC
  // PetSc variables
  KSP _ksp;
  Mat _APetsc;
  Mat _PPetsc;
  PC _preconditioner;
  Vec _XPetsc;
  Vec _BPetsc;
#endif

  // cache the columns and rows
  vector<int> _rowsIndices;
  vector<int> _colsIndices;
  vector<const Real*> _entries;

  // a copy of the sparsity structure
  map<pair<int,int>, bool> _sparsity;

  // which preconditioner to use
  PRECONDITIONER _whichPreconditioner;

  // which solver to use
  SOLVER _whichSolver;
};

#endif
