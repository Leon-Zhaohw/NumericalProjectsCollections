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

#include "PCG_MATRIX.h"

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
PCG_MATRIX::PCG_MATRIX() :
  MATRIX(),
  _firstSolve(true),
  _iterations(100),
  _eps(1e-4),
  _totalSolves(0),
  _totalIterations(0),
  _maxIterationsSeen(0)
{
}

PCG_MATRIX::PCG_MATRIX(int rows, int cols) :
  MATRIX(rows, cols),
  _firstSolve(true),
  _iterations(100),
  _eps(1e-4),
  _totalSolves(0),
  _totalIterations(0),
  _maxIterationsSeen(0)
{
}

PCG_MATRIX::PCG_MATRIX(int rows, int cols, Real* data) :
  MATRIX(rows, cols, data),
  _firstSolve(true),
  _iterations(100),
  _eps(1e-4),
  _totalSolves(0),
  _totalIterations(0),
  _maxIterationsSeen(0)
{
}

PCG_MATRIX::PCG_MATRIX(const char* filename) :
  MATRIX(filename),
  _firstSolve(true),
  _iterations(100),
  _eps(1e-4),
  _totalSolves(0),
  _totalIterations(0),
  _maxIterationsSeen(0)
{
}

PCG_MATRIX::PCG_MATRIX(const MATRIX& m) :
  MATRIX(m),
  _firstSolve(true),
  _iterations(100),
  _eps(1e-4),
  _totalSolves(0),
  _totalIterations(0),
  _maxIterationsSeen(0)
{
}

PCG_MATRIX::PCG_MATRIX(VECTOR& vec) :
  MATRIX(vec),
  _firstSolve(true),
  _iterations(100),
  _eps(1e-4),
  _totalSolves(0),
  _totalIterations(0),
  _maxIterationsSeen(0)
{
}

PCG_MATRIX::PCG_MATRIX(MATRIX3& matrix3) :
  MATRIX(matrix3),
  _firstSolve(true),
  _iterations(100),
  _eps(1e-4),
  _totalSolves(0),
  _totalIterations(0),
  _maxIterationsSeen(0)
{
}

PCG_MATRIX::~PCG_MATRIX()
{
}

//////////////////////////////////////////////////////////////////////
// init CG vectors
//////////////////////////////////////////////////////////////////////
void PCG_MATRIX::initCG()
{
  int size = this->rows();
  _residual.resizeAndWipe(size);
  _direction.resizeAndWipe(size);
  _q.resizeAndWipe(size);
}

//////////////////////////////////////////////////////////////////////
// init ICCG vars
//////////////////////////////////////////////////////////////////////
void PCG_MATRIX::initICCG()
{
  int size = this->rows();
  _s.resizeAndWipe(size);
}

//////////////////////////////////////////////////////////////////////
// solve system with CG
//////////////////////////////////////////////////////////////////////
void PCG_MATRIX::solveCG(VECTOR& x, VECTOR& b)
{
  if (_firstSolve)
  {
    if (_cols != _rows)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Matrix is not square! Cannot solve!" << endl;
      return;
    }
    
    initCG();
    _firstSolve = false;
  }
	// r = b - Ax
  _residual = (*this) * x;
  _residual = b - _residual;

	// d = r
  _direction.copyInplace(_residual);
  
	// deltaNew = transpose(r) * r
  Real deltaNew = _residual ^ _residual;
  
	// delta0 = deltaNew
  //Real delta0 = deltaNew;

	// While deltaNew > (eps^2) * delta0
	Real maxR = 2.0f * _eps;
  int i = 0;
	while ((i < _iterations) && (maxR > _eps))
	{
		// q = Ad
    _q = (*this) * _direction;

		// alpha = deltaNew / (transpose(d) * q)
    Real alpha = _direction ^ _q;
		if (fabs(alpha) > 0.0)
      alpha = deltaNew / alpha;

		// x = x + alpha * d
    x.axpy(alpha, _direction);

		// r = r - alpha * q
    _residual.axpy(-alpha, _q);

		// deltaOld = deltaNew
		Real deltaOld = deltaNew;

		// deltaNew = transpose(r) * r
		deltaNew = _residual ^ _residual;

		// beta = deltaNew / deltaOld
		Real beta = deltaNew / deltaOld;

		// d = r + beta * d
    _direction *= beta;
    _direction += _residual;

    // maxR = max(r);
    maxR = _residual.maxValue();
    
		// i = i + 1
		i++;
  }
  _totalIterations += i;
  _totalSolves++;
  cout << " total iterations: " << i << endl;

  if (i > _maxIterationsSeen)
    _maxIterationsSeen = i;
}

//////////////////////////////////////////////////////////////////////
// solve system with ICCG
//////////////////////////////////////////////////////////////////////
void PCG_MATRIX::solveICCG(VECTOR& x, VECTOR& b, Real dropTolerance)
{
  if (_firstSolve)
  {
    if (_cols != _rows)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Matrix is not square! Cannot solve!" << endl;
      return;
    }
    
    initCG();
    initICCG();
    _firstSolve = false;
  }

  // compute IC factorization
  computeIC(dropTolerance);
  
	// r = b - Ax
  _residual = (*this) * x;
  _residual = b - _residual;

	// d = (M^-1) * r
  solveIC(_direction, _residual);
  
	// deltaNew = transpose(r) * d
  Real deltaNew = _residual ^ _direction;
  
	// delta0 = deltaNew
  //Real delta0 = deltaNew;

	// While deltaNew > (eps^2) * delta0
	Real maxR = 2.0f * _eps;
  int i = 0;
	while ((i < _iterations) && (maxR > _eps))
	{
		// q = Ad
    _q = (*this) * _direction;

		// alpha = deltaNew / (transpose(d) * q)
    Real alpha = _direction ^ _q;
		if (fabs(alpha) > 0.0)
      alpha = deltaNew / alpha;

		// x = x + alpha * d
    x.axpy(alpha, _direction);

		// r = r - alpha * q
    _residual.axpy(-alpha, _q);

		// s = (M^-1) * r
    solveIC(_s, _residual);

		// deltaOld = deltaNew
		Real deltaOld = deltaNew;

		// deltaNew = transpose(r) * s
		deltaNew = _residual ^ _s;

		// beta = deltaNew / deltaOld
		Real beta = deltaNew / deltaOld;

		// d = s + beta * d
    _direction *= beta;
    _direction += _s;

    // maxR = max(r);
    maxR = _residual.maxValue();
    
		// i = i + 1
		i++;
  }
  _totalIterations += i;
  _totalSolves++;

  if (i > _maxIterationsSeen)
    _maxIterationsSeen = i;
}

//////////////////////////////////////////////////////////////////////
// compute Incomplete Cholesky factorization
//////////////////////////////////////////////////////////////////////
void PCG_MATRIX::computeIC(Real dropTolerance)
{
  _IC = *this;

  // See Golub and Van Loan, 10.3.2
  int n = _rows;
  for (int k = 0; k < n; k++)
  {
    _IC(k,k) = sqrt(_IC(k,k));
    Real diagonal = 1.0 / _IC(k,k);

    for (int i = k + 1; i < n; i++)
      if (fabs(_IC(i,k)) > dropTolerance)
        _IC(i,k) *= diagonal;

    for (int j = k + 1; j < n; j++)
      for (int i = j; i < n; i++)
        if (fabs(_IC(i,j)) > dropTolerance)
          _IC(i,j) -= _IC(i,k) * _IC(j,k);
  }
}

//////////////////////////////////////////////////////////////////////
// apply the Incomplete Cholesky factorization
//////////////////////////////////////////////////////////////////////
void PCG_MATRIX::solveIC(VECTOR& x, VECTOR& b)
{
  // forward substitute
  for (int row = 0; row < _rows; row++)
  {
    Real sum = b(row);
    for (int col = 0; col < row; col++)
      sum -= _IC(row,col) * x(col);
    x(row) = sum / _IC(row, row);
  }
  
  // backward substitute
  for (int row = _rows - 1; row >= 0; row--)
  {
    Real sum = x(row);
    for (int col = _cols - 1; col > row; col--)
      sum -= _IC(col, row) * x(col);
    x(row) = sum / _IC(row, row);
  }
}
