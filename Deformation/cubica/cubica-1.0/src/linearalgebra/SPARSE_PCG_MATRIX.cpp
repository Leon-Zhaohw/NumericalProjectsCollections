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

#include "SPARSE_PCG_MATRIX.h"
#include <TIMER.h>

#ifdef USING_MKL
#include <mkl_pardiso.h>
#endif
#ifdef USING_OPENMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor for the sparse matrix
//////////////////////////////////////////////////////////////////////
SPARSE_PCG_MATRIX::SPARSE_PCG_MATRIX(int rows, int cols, Real eps, int iterations) :
  SPARSE_MATRIX(rows, cols),
  _iterations(iterations),
  _eps(eps),
  _firstSolve(true),
  _totalIterations(0),
  _totalSolves(0),
  _totalResidual(0.0),
  _totalTime(0.0),
  _qs(NULL),
  _totalCores(0)
{
}

SPARSE_PCG_MATRIX::SPARSE_PCG_MATRIX(Real eps, int iterations) :
  SPARSE_MATRIX(),
  _iterations(iterations),
  _eps(eps),
  _firstSolve(true),
  _totalIterations(0),
  _totalSolves(0),
  _totalResidual(0.0),
  _totalTime(0.0),
  _qs(NULL),
  _totalCores(0)
{
}

SPARSE_PCG_MATRIX::SPARSE_PCG_MATRIX(MATRIX& matrix, Real eps, int iterations) :
  SPARSE_MATRIX(matrix),
  _iterations(iterations),
  _eps(eps),
  _firstSolve(true),
  _totalIterations(0),
  _totalSolves(0),
  _totalResidual(0.0),
  _totalTime(0.0),
  _qs(NULL),
  _totalCores(0)
{
}

SPARSE_PCG_MATRIX::~SPARSE_PCG_MATRIX()
{
  if (_qs) delete[] _qs;
}

void SPARSE_PCG_MATRIX::setSparsity(SPARSE_MATRIX& A)
{
  // iterate through all the entries
  map<pair<int,int>, Real>::const_iterator i;
  const map<pair<int,int>, Real>& data = A.matrix();
  for (i = data.begin(); i != data.end(); i++)
  {
    const pair<int,int> index = i->first;
    const Real value = i->second;
    _matrix[index] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////
// solve an SPD system using CG
//////////////////////////////////////////////////////////////////////
void SPARSE_PCG_MATRIX::solveCG(VECTOR& x, VECTOR& b)
{
  TIMER total;
  if (_firstSolve)
  {
    initCG();
    _firstSolve = false;
  }

	// r = b - Ax
  //_residual = b - (*this) * x;
  _residual = (*this) * x;
  _residual = b - _residual;

	// d = r
  _direction.copyInplace(_residual);
  
	// deltaNew = transpose(r) * r
  Real deltaNew = _residual ^ _residual;
  
	// delta0 = deltaNew
  Real delta0 = deltaNew;

#ifdef USING_OPENMP
  // precache as much as possible before going into the main PCG loop
  initParallelMultiply();
#endif

	// While deltaNew > (eps^2) * delta0
	Real maxR = 2.0f * _eps;
  int i = 0;
#ifndef _WIN32  
  cout << " Residual: ";
#endif
	while ((i < _iterations) && (maxR > _eps))
	{
		// q = Ad
#ifdef USING_OPENMP
    parallelMultiply();
#else
    TIMER multiply;
    _q = (*this) * _direction;
    _timingBreakdown["Matrix multiply"] += multiply.timing();
#endif

		// alpha = deltaNew / (transpose(d) * q)
    TIMER alphaTimer;
    Real alpha = _direction ^ _q;
		if (fabs(alpha) > 0.0)
      alpha = deltaNew / alpha;
    _timingBreakdown["alpha compute"] += alphaTimer.timing();

		// x = x + alpha * d
    TIMER xTimer;
    x.axpy(alpha, _direction);
    _timingBreakdown["x update"] += xTimer.timing();

		// r = r - alpha * q
    TIMER residualTimer;
    _residual.axpy(-alpha, _q);
    _timingBreakdown["residual update"] += residualTimer.timing();

		// deltaOld = deltaNew
		Real deltaOld = deltaNew;

		// deltaNew = transpose(r) * r
    TIMER deltaTimer;
		deltaNew = _residual ^ _residual;
    _timingBreakdown["delta update"] += deltaTimer.timing();

		// beta = deltaNew / deltaOld
		Real beta = deltaNew / deltaOld;

		// d = r + beta * d
    TIMER directionTimer;
    _direction *= beta;
    _direction += _residual;
    _timingBreakdown["direction update"] += directionTimer.timing();

    // maxR = max(r);
    TIMER maxRTimer;
    maxR = _residual.maxValue();
    _timingBreakdown["maxR update"] += maxRTimer.timing();

		// i = i + 1
		i++;

#ifndef _WIN32  
    if (i % 10 == 0)
      cout << maxR << " ";
#endif
  }
#ifndef _WIN32  
  cout << endl;
#endif
  _totalIterations += i;
  _totalSolves++;
  _totalResidual += maxR;

  cout << " PCG iterations: " << i << endl;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// solve an SPD system using PCG
//////////////////////////////////////////////////////////////////////
void SPARSE_PCG_MATRIX::solvePCG(VECTOR& x, VECTOR& b, PRECONDITIONER* M)
{
  TIMER total;

  if (_firstSolve)
  {
    initCG();
    initPCG();
    _firstSolve = false;
  }

  TIMER Msetup;
  // initialize the preconditioner
  M->init();
  _timingBreakdown["Precond. Setup"] += Msetup.timing();

	// r = b - Ax
  //_residual = b - (*this) * x;
  _residual = (*this) * x;
  _residual = b - _residual;

  // d = (M^-1) * r
  M->solve(_direction, _residual);

	// deltaNew = transpose(r) * d
  Real deltaNew = _residual ^ _direction;
  
	// delta0 = deltaNew
  Real delta0 = deltaNew;

#ifdef USING_OPENMP
  // precache as much as possible before going into the main PCG loop
  initParallelMultiply();
#endif

	// While deltaNew > (eps^2) * delta0
	Real maxR = 2.0f * _eps;
  int i = 0;
  //cout << " Residual: ";
	while ((i < _iterations) && (maxR > _eps))
	{
#ifdef USING_OPENMP
    parallelMultiply();
#else
    TIMER multiply;
    _q = (*this) * _direction;
    _timingBreakdown["Matrix multiply"] += multiply.timing();
#endif

		// alpha = deltaNew / (transpose(d) * q)
    Real alpha = _direction ^ _q;
		if (fabs(alpha) > 0.0)
      alpha = deltaNew / alpha;
    else
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " PCG BREAKDOWN" << endl;
    }

		// x = x + alpha * d
    x.axpy(alpha, _direction);

		// r = r - alpha * q
    _residual.axpy(-alpha, _q);

		// s = (M^-1) * r
    TIMER Msolve;
    M->solve(_s, _residual);
    _timingBreakdown["Precond. solve"] += Msolve.timing();

		// deltaOld = deltaNew
		Real deltaOld = deltaNew;

		// deltaNew = transpose(r) * s
		deltaNew = _residual ^ _s;

		// beta = deltaNew / deltaOld
		Real beta = deltaNew / deltaOld;

		// d = r + beta * d
    _direction *= beta;
    _direction += _s;

    // maxR = max(r);
    maxR = _residual.maxValue();

		// i = i + 1
		i++;

    //if (i % 10 == 0)
    //  cout << maxR << " ";
  }
  cout << endl;

  _totalIterations += i;
  _totalSolves++;
  _totalResidual += maxR;
  cout << " total iterations: " << i << endl;
  //_residual = b - (*this) * x;
  _residual = (*this) * x;
  _residual = b - _residual;
  cout << " residual: " << _residual.normInf() << endl;
  _totalTime += total.timing();
}

//////////////////////////////////////////////////////////////////////
// initialize CG solver for first solve
//////////////////////////////////////////////////////////////////////
void SPARSE_PCG_MATRIX::initCG()
{
  int size = _rows;
  _residual.resizeAndWipe(size);
  _direction.resizeAndWipe(size);
  _q.resizeAndWipe(size);
}

//////////////////////////////////////////////////////////////////////
// initialize PCG solver for first solve
//////////////////////////////////////////////////////////////////////
void SPARSE_PCG_MATRIX::initPCG()
{
  int size = _rows;
  _s.resizeAndWipe(size);
}

//////////////////////////////////////////////////////////////////////
// precache as much as possible before going into the main PCG loop
//////////////////////////////////////////////////////////////////////
void SPARSE_PCG_MATRIX::initParallelMultiply()
{
  // if this is the first call, cache the pair and Real locations
  if (_entries.size() == 0)
  {
    TIMER cacheTimer;
    map<pair<int,int>, Real>::iterator iter;
    for (iter = _matrix.begin(); iter != _matrix.end(); iter++)
    {
      int row = iter->first.first;
      int col = iter->first.second;
      //_pairs.push_back(iter->first);
      _rowIndices.push_back(row);
      _columnIndices.push_back(col);
      _entries.push_back(iter->second);
    }
#ifdef USING_OPENMP
    //_qs.resize(omp_get_max_threads());
    _totalCores = omp_get_max_threads();
#else
    //_qs.resize(1);
    _totalCores = 1;
#endif
    _qs = new VECTOR[_totalCores];
    //for (int x = 0; x < _qs.size(); x++)
    for (int x = 0; x < _totalCores; x++)
      _qs[x].resizeAndWipe(_q.size());
    _timingBreakdown["Pair caching"] += cacheTimer.timing();
  }
  else
  {
    TIMER cacheTimer;
    int i = 0;
    map<pair<int,int>, Real>::iterator iter;
    for (iter = _matrix.begin(); iter != _matrix.end(); iter++, i++)
      _entries[i] = iter->second;
    _timingBreakdown["Entry caching"] += cacheTimer.timing();
  }
}

//////////////////////////////////////////////////////////////////////
// do a parallel multiply that OpenMP can handle
//////////////////////////////////////////////////////////////////////
void SPARSE_PCG_MATRIX::parallelMultiply()
{
  TIMER multiply;
#ifdef USING_OPENMP
#pragma omp parallel
#endif
  { 
#ifdef USING_OPENMP
    const int id  = omp_get_thread_num();
#else
    const int id  = 0;
#endif
    VECTOR& q = _qs[id];
    q.clear();
#pragma omp for  schedule(static)
    for (int x = 0; x < _entries.size(); x++)
    {
      //const pair<int,int> index = _pairs[x];
      //const Real value = *(_entries[x]);
      //int i = index.first;
      //int j = index.second;
      //_qs[id](i) += _direction[j] * value;
      //_qs[id](index.first) += _direction[index.second] * _entries[x];
      //_qs[id](_rowIndices[x]) += _direction[_columnIndices[x]] * _entries[x];
      q[_rowIndices[x]] += _direction[_columnIndices[x]] * _entries[x];
    }
  }
  _timingBreakdown["Matrix multiply"] += multiply.timing();

  // combine qs
  TIMER combineTimer;
  _q.clear();
  for (int x = 0; x < _totalCores; x++)
    _q += _qs[x];
  _timingBreakdown["Final combine"] += combineTimer.timing();
}

