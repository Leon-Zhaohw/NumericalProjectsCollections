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
// PCG_MATRIX.h: interface for the PCG_MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PCG_MATRIX_H
#define PCG_MATRIX_H

#include <MATRIX.h>

//////////////////////////////////////////////////////////////////////
// A matrix that supports preconditioned conjugate gradient.
//
// This has been subclassed because there are significant temp arrays
// for the solver, and all matrix objects should not have to incur
// this memory overhead.
//////////////////////////////////////////////////////////////////////
class PCG_MATRIX : public MATRIX {

public:
  PCG_MATRIX();
  PCG_MATRIX(int rows, int cols);
  PCG_MATRIX(int rows, int cols, Real* data);
  PCG_MATRIX(const char* filename);
  PCG_MATRIX(const MATRIX& m);
  PCG_MATRIX(VECTOR& vec);
  PCG_MATRIX(MATRIX3& matrix3);
  virtual ~PCG_MATRIX();

  void solveCG(VECTOR& x, VECTOR& b);
  void solveICCG(VECTOR& x, VECTOR& b, Real dropTolerance = 1e-8);
  Real meanIterations()  { return (Real)_totalIterations / _totalSolves; };

  // accessors
  int& maxIterations()     { return _iterations; };
  int& maxIterationsSeen() { return _maxIterationsSeen; };
  Real& eps()              { return _eps; };
  VECTOR& residual()       { return _residual; };
  VECTOR& direction()      { return _direction; };
  VECTOR& q()              { return _q; };
  int& totalSolves()       { return _totalSolves; };
  int& totalIterations()   { return _totalIterations; };

protected:
  // CG variables
  VECTOR _residual;
  VECTOR _direction;
  VECTOR _q;

  // ICCG variables
  VECTOR _s;
  MATRIX _IC;
  
  bool _firstSolve;
  int _iterations;
  Real _eps;

  // keep track of some stats
  int _totalSolves;
  int _totalIterations;
  int _maxIterationsSeen;
  
  void initCG();
  void initICCG();

  // compute Incomplete Cholesky factorization
  void computeIC(Real dropTolerance);

  // apply the Incomplete Cholesky factorization
  void solveIC(VECTOR& x, VECTOR& b);
};

#endif

