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
 
 The Modified Incomplete Cholesky preconditioner was written by
 John Delaney
 */
#include "FLUID_3D_MIC.h"
#include <zlib.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FLUID_3D_MIC::FLUID_3D_MIC(int xRes, int yRes, int zRes, int amplify, unsigned int* boundaries) :
	_xRes(xRes), _yRes(yRes), _zRes(zRes), _res(0.)
{
	// set simulation consts
  _dt = 0.10; 

  _iterations = 1000;
  _buoyancy = 0.1f;
  // ADJ: modifying buoyancy here
  // _buoyancy = 0.05f;
  _wind = 0.05f;
  _heatDiffusion = 1e-3;
  _vorticityEps = 2.0;
	_totalTime = 0.0f;
	_totalSteps = 0;
	_res = VEC3I(_xRes,_yRes,_zRes);
	_maxRes = _res.maxElement();

  // scale the constants according to the refinement of the grid
	_dx = 1.0f / (Real)_maxRes;
	Real scaling = 64.0f / _maxRes;
	scaling = (scaling < 1.0f) ? 1.0f : scaling;
	_vorticityEps /= scaling;

  // precision of the CG solver
  _solverEps = 1e-10;

  _center = VEC3F(0,0,0);
  _lengths = VEC3F(_dx * _xRes, _dx * _yRes, _dx * _zRes);

  _density    = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _densityOld = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _heat    = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _heatOld = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _pressure = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _divergence = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _residual = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _direction = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _q = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);

  _dudt0 = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _dudt1 = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);

  _velocity = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _velocityOld = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);

  _force = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _vorticity = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);

	_precon = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_z      = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Adiag  = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aplusi = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aplusj = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aplusk = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aminui = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aminuj = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aminuk = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	
	// allocate arrays
	_totalCells   = _xRes * _yRes * _zRes;
	_slabSize = _xRes * _yRes;
	_xs = 1; _ys = _xRes; _zs = _slabSize;
	_obstacles    = new unsigned char[_totalCells];
	_dirichletObstacleField = new bool[_totalCells];
	_neumannObstacleField = new bool[_totalCells];

  _neumannVelocityField = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	
	for (int x = 0; x < _totalCells; x++)
  {
		_obstacles[x]    = false;
    _dirichletObstacleField[x] = false;
    _neumannObstacleField[x] = false;
  }
  
  // set up the boundary conditions
  if (boundaries != NULL)
  {
    _domainBcFront  = boundaries[0];
    _domainBcBack   = boundaries[1];
    _domainBcLeft   = boundaries[2];
    _domainBcRight  = boundaries[3];
    _domainBcTop    = boundaries[4];
    _domainBcBottom = boundaries[5];
  }
  else
  {
    _domainBcFront  = 1; 
    _domainBcBack   = 1;
    _domainBcLeft   = 1;
    _domainBcRight  = 1;
    _domainBcTop    = 0;
    _domainBcBottom = 0;
  }

	// set side obstacles
  int index;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      // front slab
      index = x + y * _xRes;
      if(_domainBcFront==1) _obstacles[index] = 1;

      // back slab
      index += _totalCells - _slabSize;
      if(_domainBcBack==1) _obstacles[index] = 1;
    }
  for (int z = 0; z < _zRes; z++)
    for (int x = 0; x < _xRes; x++)
    {
      // bottom slab
      index = x + z * _slabSize;
      if(_domainBcBottom==1) _obstacles[index] = 1;

      // top slab
      index += _slabSize - _xRes;
      if(_domainBcTop==1) _obstacles[index] = 1;
    }
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
    {
      // left slab
      index = y * _xRes + z * _slabSize;
      if(_domainBcLeft==1) _obstacles[index] = 1;

      // right slab
      index += _xRes - 1;
      if(_domainBcRight==1) _obstacles[index] = 1;
    }
}

FLUID_3D_MIC::FLUID_3D_MIC() :
  _xRes(-1), _yRes(-1), _zRes(-1), _totalCells(-1), 
  _obstacles(NULL), _dirichletObstacleField(NULL), _neumannObstacleField(NULL)
{
}

FLUID_3D_MIC::~FLUID_3D_MIC()
{
	if (_obstacles) delete[] _obstacles;
	if (_dirichletObstacleField) delete[] _dirichletObstacleField;
	if (_neumannObstacleField) delete[] _neumannObstacleField;
}

//////////////////////////////////////////////////////////////////////
// step simulation once with fanciness like vorticity confinement 
// and diffusion removed
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::step()
{
  Real goalTime = 0.1;
  Real currentTime = 0;

  // compute the CFL condition
  _dt = goalTime;

  // wipe forces
  _force.clear();

  // wipe boundaries
  _velocity.setZeroBorder();
  _density.setZeroBorder();

  // compute the forces
  addBuoyancy(_heat.data());
  _velocity.axpy(_dt, _force);
  _force.clear();

  _prevorticity = _velocity;

  addVorticity();
  _velocity.axpy(_dt, _force);

  // advect everything
  advectStam();

  _prediffusion = _velocity;

  if (_peeledDampingFull.rows() == 0)
    buildPeeledDampingMatrixFull();
  VectorXd after = _peeledDampingFull * _velocity.peelBoundary().flattenedEigen();
  _velocity.setWithPeeled(after);

  _preprojection = _velocity;

  // run the solvers
  project();

  currentTime += _dt;

	_totalTime += goalTime;
	_totalSteps++;
}

//////////////////////////////////////////////////////////////////////
// step simulation once with fanciness like vorticity confinement 
// and diffusion removed
// Has a moving fan obstacle
// Is *NOT* reordered!
// We do:
// 1) forces
// 2) advection
// 3) diffusion
// 4) boundary + pressure project (IOP)
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::stepWithMovingObstacle(BOX* box)
{
  Real goalTime = 0.1;
  Real currentTime = 0;

  // compute the CFL condition
  _dt = goalTime;

  // wipe forces
  _force.clear();

  // wipe boundaries
  _velocity.setZeroBorder();
  _density.setZeroBorder();

  // compute the forces
  _velocity.axpy(_dt, _force);
  _force.clear();

  _prevorticity = _velocity;

  addVorticity();
  _velocity.axpy(_dt, _force);

  // advect everything (this caches preadvection)
  advectStam();

  // cache prediffusion
  _prediffusion = _velocity;

  // if the matrix isn't built yet, build it
  if (_peeledDampingFull.rows() <= 0) {
    buildPeeledDampingMatrixFull();
  }

  VectorXd after = _peeledDampingFull * _velocity.peelBoundary().flattenedEigen();
  _velocity.setWithPeeled(after);

  // cache preprojection before doing IOP
  _preprojection = _velocity;

  // let's test IOP with the full matrix
 
  setPeeledSparseMovingIOP(box);
  setVelocityNeumann();
  VectorXd afterIOP = _neumannIOP * _velocityNeumann; 
  _velocity.setWithPeeled(afterIOP);
  
  // store the postIOP velocity (QUESTION: before or after the project?)
  _postIOP = _velocity;

  // project via Poisson
  project();

  // ADJ: commenting out the repetition of IOP for the time being.
  /*
  // repeat to do the iteration!
  setVelocityNeumann();
  afterIOP = _neumannIOP * _velocityNeumann;
  _velocity.setWithPeeled(afterIOP);

  project();
  */

  // we don't really need to cache this since it is the final one, but it doesn't hurt 
  _postIOPAndPressure = _velocity;

  // this corresponds to doing only 2 iterations of IOP, but that may be sufficient
  // for practical use

  currentTime += _dt;

	_totalTime += goalTime;
	_totalSteps++;
}
//////////////////////////////////////////////////////////////////////
// project into divergence free field
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::project()
{
  TIMER functionTimer(__FUNCTION__);
	int index, x, y, z;
	setObstacleBoundaries();

	// copy out the boundaries
	if(_domainBcLeft == 0)  
    _velocity.setNeumannX();
	else 
    _velocity.setZeroX();

	if(_domainBcTop == 0)
    _velocity.setNeumannY();
	else 
    _velocity.setZeroY();

	if(_domainBcFront == 0) 
    _velocity.setNeumannZ();
	else
    _velocity.setZeroZ();

	// calculate divergence
	index = _slabSize + _xRes + 1;
	for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
		for (y = 1; y < _yRes - 1; y++, index += 2)
			for (x = 1; x < _xRes - 1; x++, index++)
			{
				Real xright = _velocity[index + 1][0];
				Real xleft  = _velocity[index - 1][0];
				Real yup    = _velocity[index + _xRes][1];
				Real ydown  = _velocity[index - _xRes][1];
				Real ztop   = _velocity[index + _slabSize][2];
				Real zbottom = _velocity[index - _slabSize][2];

				if(_obstacles[index+1]) xright = - _velocity[index][0];
				if(_obstacles[index-1]) xleft  = - _velocity[index][0];
				if(_obstacles[index+_xRes]) yup    = - _velocity[index][1];
				if(_obstacles[index-_xRes]) ydown  = - _velocity[index][1];
				if(_obstacles[index+_slabSize]) ztop    = - _velocity[index][2];
				if(_obstacles[index-_slabSize]) zbottom = - _velocity[index][2];

				_divergence[index] = -_dx * 0.5f * (
						xright - xleft +
						yup - ydown + 
						ztop - zbottom );
				_pressure[index] = 0.0f;
			}
	_pressure.copyBorderAll();

  _cachedDivergence = _divergence;

	// solve Poisson equation
	solvePoisson(_pressure, _divergence, _obstacles, false);

  _cachedPressure = _pressure;

	// project out solution -- more matrix-friendly version
	Real invDx = 1.0f / _dx;
	index = _slabSize + _xRes + 1;
	for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
		for (y = 1; y < _yRes - 1; y++, index += 2)
			for (x = 1; x < _xRes - 1; x++, index++)
			{
        if (x == 1)
				  _velocity[index][0] -= (_pressure[index + 1] - _pressure[index]) * invDx;
        else if (x == _xRes - 2)
				  _velocity[index][0] -= (_pressure[index] - _pressure[index - 1]) * invDx;
        else
				  _velocity[index][0] -= 0.5f * (_pressure[index + 1]     - _pressure[index - 1])     * invDx;

        if (y == 1)
  				_velocity[index][1] -= (_pressure[index + _xRes]  - _pressure[index]) * invDx;
        else if (y == _yRes - 2)
  				_velocity[index][1] -= (_pressure[index]  - _pressure[index - _xRes]) * invDx;
        else
				  _velocity[index][1] -= 0.5f * (_pressure[index + _xRes]  - _pressure[index - _xRes]) * invDx;

        if (z == 1)
				  _velocity[index][2] -= (_pressure[index + _slabSize] - _pressure[index]) * invDx;
        else if (z == _zRes - 2)
				  _velocity[index][2] -= (_pressure[index] - _pressure[index - _slabSize]) * invDx;
        else
				  _velocity[index][2] -= 0.5f * (_pressure[index + _slabSize] - _pressure[index - _slabSize]) * invDx;
			}
}

//////////////////////////////////////////////////////////////////////
// calculate the obstacle directional types
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::setObstacleBoundaries()
{
	// cull degenerate obstacles , move to addObstacle?
	for (int z = 1, index = _slabSize + _xRes + 1; 
			z < _zRes - 1; z++, index += 2 * _xRes)
		for (int y = 1; y < _yRes - 1; y++, index += 2)
			for (int x = 1; x < _xRes - 1; x++, index++)
				if (_obstacles[index] != EMPTY)
				{
					const int top   = _obstacles[index + _slabSize];
					const int bottom= _obstacles[index - _slabSize];
					const int up    = _obstacles[index + _xRes];
					const int down  = _obstacles[index - _xRes];
					const int left  = _obstacles[index - 1];
					const int right = _obstacles[index + 1];

					int counter = 0;
					if (up)    counter++;
					if (down)  counter++;
					if (left)  counter++;
					if (right) counter++;
					if (top)  counter++;
					if (bottom) counter++;

					if (counter < 3)
						_obstacles[index] = EMPTY;
				}

	// tag remaining obstacle blocks
	for (int z = 1, index = _slabSize + _xRes + 1; 
			z < _zRes - 1; z++, index += 2 * _xRes)
		for (int y = 1; y < _yRes - 1; y++, index += 2)
			for (int x = 1; x < _xRes - 1; x++, index++)
		{
			// could do cascade of ifs, but they are a pain
			if (_obstacles[index] != EMPTY)
			{
				const int top   = _obstacles[index + _slabSize];
				const int bottom= _obstacles[index - _slabSize];
				const int up    = _obstacles[index + _xRes];
				const int down  = _obstacles[index - _xRes];
				const int left  = _obstacles[index - 1];
				const int right = _obstacles[index + 1];

        _velocity[index] = 0.0f;
				_pressure[index] = 0.0f;

				// average pressure neighbors
				Real pcnt = 0.;
				if (left && !right) {
					_pressure[index] += _pressure[index + 1];
					pcnt += 1.;
				}
				if (!left && right) {
					_pressure[index] += _pressure[index - 1];
					pcnt += 1.;
				}
				if (up && !down) {
					_pressure[index] += _pressure[index - _xRes];
					pcnt += 1.;
				}
				if (!up && down) {
					_pressure[index] += _pressure[index + _xRes];
					pcnt += 1.;
				}
				if (top && !bottom) {
					_pressure[index] += _pressure[index - _xRes];
					pcnt += 1.;
				}
				if (!top && bottom) {
					_pressure[index] += _pressure[index + _xRes];
					pcnt += 1.;
				}
				_pressure[index] /= pcnt; 

				// TODO? set correct velocity bc's
				// velocities are only set to zero right now
				// this means it's not a full no-slip boundary condition
				// but a "half-slip" - still looks ok right now
			}
		}
}


//////////////////////////////////////////////////////////////////////
// add buoyancy forces
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::addBuoyancy(Real *field)
{
  TIMER functionTimer(__FUNCTION__);
	int index = 0;

	Real beta = _buoyancy;
	if(beta==0.) return;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
			for (int x = 0; x < _xRes; x++, index++) 
        _force[index][1] += beta * field[index];
}


//////////////////////////////////////////////////////////////////////
// add 'wind' forces
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::addWind(Real* field)
{
 TIMER functionTimer(__FUNCTION__);
	int index = 0;

	Real beta = _wind;
	if(beta==0.) return;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
			for (int x = 0; x < _xRes; x++, index++) 
        _force[index][0] += beta * field[index];
}

//////////////////////////////////////////////////////////////////////
// add vorticity to the force field
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::addVorticity()
{
  TIMER functionTimer(__FUNCTION__);
	int x,y,z,index;
	if(_vorticityEps<=0.) return;

	// calculate vorticity
	Real gridSize = 0.5f / _dx;
	index = _slabSize + _xRes + 1;
	for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
		for (y = 1; y < _yRes - 1; y++, index += 2)
			for (x = 1; x < _xRes - 1; x++, index++)
			{
				int up    = _obstacles[index + _xRes] || y == _yRes - 2? index : index + _xRes;
				int down  = _obstacles[index - _xRes] || y == 1? index : index - _xRes;
				Real dy  = (up == index || down == index) ? 1.0f / _dx : gridSize;
				int out   = _obstacles[index + _slabSize] || z == _zRes - 2? index : index + _slabSize;
				int in    = _obstacles[index - _slabSize] || z == 1 ? index : index - _slabSize;
				Real dz  = (out == index || in == index) ? 1.0f / _dx : gridSize;
				int right = _obstacles[index + 1] || x == _xRes - 2 ? index : index + 1;
				int left  = _obstacles[index - 1] || x == 1 ? index : index - 1;
				Real dx  = (left == index || right == index) ? 1.0f / _dx : gridSize;

				_vorticity[index][0] = (_velocity[up][2] - _velocity[down][2]) * dy + (-_velocity[out][1] + _velocity[in][1]) * dz;
				_vorticity[index][1] = (_velocity[out][0] - _velocity[in][0]) * dz + (-_velocity[right][2] + _velocity[left][2]) * dx;
				_vorticity[index][2] = (_velocity[right][1] - _velocity[left][1]) * dx + (-_velocity[up][0] + _velocity[down][0])* dy;
			}

  FIELD_3D vorticity = _vorticity.magnitudeField();

	// calculate normalized vorticity vectors
	Real eps = _vorticityEps;
	index = _slabSize + _xRes + 1;
	for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
		for (y = 1; y < _yRes - 1; y++, index += 2)
			for (x = 1; x < _xRes - 1; x++, index++)
				if (!_obstacles[index])
				{
					Real N[3];

				  int up    = _obstacles[index + _xRes] || y == _yRes - 2? index : index + _xRes;
				  int down  = _obstacles[index - _xRes] || y == 1? index : index - _xRes;
					Real dy  = (up == index || down == index) ? 1.0f / _dx : gridSize;
				  int out   = _obstacles[index + _slabSize] || z == _zRes - 2? index : index + _slabSize;
				  int in    = _obstacles[index - _slabSize] || z == 1 ? index : index - _slabSize;
					Real dz  = (out == index || in == index) ? 1.0f / _dx : gridSize;
  				int right = _obstacles[index + 1] || x == _xRes - 2 ? index : index + 1;
	  			int left  = _obstacles[index - 1] || x == 1 ? index : index - 1;
					Real dx  = (right == index || left == index) ? 1.0f / _dx : gridSize;
					N[0] = (vorticity[right] - vorticity[left]) * dx; 
					N[1] = (vorticity[up] - vorticity[down]) * dy;
					N[2] = (vorticity[out] - vorticity[in]) * dz;

					Real magnitude = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
					if (magnitude > 0.0f)
					{
						magnitude = 1.0 / magnitude;
						N[0] *= magnitude;
						N[1] *= magnitude;
						N[2] *= magnitude;

						_force[index][0] += (N[1] * _vorticity[index][2] - N[2] * _vorticity[index][1]) * _dx * eps;
						_force[index][1] -= (N[0] * _vorticity[index][2] - N[2] * _vorticity[index][0]) * _dx * eps;
						_force[index][2] += (N[0] * _vorticity[index][1] - N[1] * _vorticity[index][0]) * _dx * eps;
					}
				}
}

//////////////////////////////////////////////////////////////////////
// Advect using the semi-Lagrangian method of Stam
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::advectStam()
{
  TIMER functionTimer(__FUNCTION__);
	VEC3I res = VEC3I(_xRes,_yRes,_zRes);

	if(_domainBcLeft == 0) 
    _velocity.copyBorderX();
	else 
    _velocity.setZeroX();

	if(_domainBcTop == 0) 
    _velocity.copyBorderY();
	else 
    _velocity.setZeroY();

	if(_domainBcFront == 0) 
    _velocity.copyBorderZ();
	else 
    _velocity.setZeroZ();

  _preadvection = _velocity;

  _velocity.swapPointers(_velocityOld);
	_density.swapPointers(_densityOld);
	_heat.swapPointers(_heatOld);

	const Real dt0 = _dt / _dx;
  VECTOR3_FIELD_3D::advect(dt0, _velocityOld, _densityOld, _density);
  VECTOR3_FIELD_3D::advect(dt0, _velocityOld, _heatOld, _heat);
  VECTOR3_FIELD_3D::advect(dt0, _velocityOld, _velocityOld, _velocity);

  _cachedQuadratic = _velocity - _velocityOld;
  _postadvection = _velocity;

  // this seems to fix things relative to explicit advection --
  // otherwise new values appear at the boundaries
  _velocity.setZeroBorder();

	if(_domainBcLeft == 0) 
    _velocity.copyBorderX();
	else 
    _velocity.setZeroX();

	if(_domainBcTop == 0)
    _velocity.copyBorderY();
  else
    _velocity.setZeroY();

	if(_domainBcFront == 0) 
    _velocity.copyBorderZ();
	else 
    _velocity.setZeroZ();

	_density.setZeroBorder();
	_heat.setZeroBorder();
}

//////////////////////////////////////////////////////////////////////
// add a test cube of density to the center
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::addSmokeColumn() 
{
	addSmokeTestCase(_density.data(), _res);
	addSmokeTestCase(_heat.data(), _res);
}

//////////////////////////////////////////////////////////////////////
// add a test portion of a sphere of density to the center
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::addSmokeSphere() 
{
	addSmokeSphereTestCase(_density.data(), _res);
	addSmokeSphereTestCase(_heat.data(), _res);
}

//////////////////////////////////////////////////////////////////////
// generic static version, so that it can be applied to the
// WTURBULENCE grid as well
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::addSmokeTestCase(Real* field, VEC3I res)
{
	const int slabSize = res[0]*res[1]; 
  int maxRes = res.maxElement();
	Real dx = 1.0f / (Real)maxRes;

	Real xTotal = dx * res[0];
	Real zTotal = dx * res[2];

  Real heighMin = 0.05;
  Real heighMax = 0.10;

  for (int z = 0; z < res[2]; z++)
    for (int y = (int)(heighMin*res[1]); y <= (int)(heighMax * res[1]); y++)
      for (int x = 0; x < res[0]; x++)
      {
        Real xLength = x * dx - xTotal * 0.4f;
        Real zLength = z * dx - zTotal * 0.5f;
        Real radius = sqrtf(xLength * xLength + zLength * zLength);

        if (radius < 0.075f * xTotal)
        {
          int index = x + y * res[0] + z * slabSize;
          field[index] = 1.0f;
        }
      }
}

//////////////////////////////////////////////////////////////////////
// setting the initial conditions for the spinning obstacle 
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::addSmokeSphereTestCase(Real* field, VEC3I res)
{
	const int slabSize = res[0]*res[1]; 

  Real r = 0.2;
  Real dr = 0.03;
  VEC3F center(0.5, 0.5, 0.5);
  Real heightMin1 = center[1] - r;
  Real heightMax1 = heightMin1 + 0.1; 
  Real heightMax2 = center[1] + r; 
  Real heightMin2 = heightMax2 - 0.1;

  VEC3F point(0.0, 0.0, 0.0);

  for (int z = 0; z < res[2]; z++) {
    for (int y = 0; y < res[1]; y++) {
      for (int x = 0; x < res[0]; x++) {
        if ( ( y >= heightMin1 * res[1] && y <= heightMax1 * res[1] ) ||
             ( y >= heightMin2 * res[1] && y <= heightMax2 * res[1] ) ) {

          point = this->cellCenter(x, y, z);
          point -= center;
          Real length = norm(point); 

          if (length <= r && length >= r - dr) {
            int index = x + y * res[0] + z * slabSize;
            field[index] = 1.0f;
          }
        }
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////
// set the fluid initial velocity at the obstacle
//////////////////////////////////////////////////////////////////////

void FLUID_3D_MIC::setInitialVelocity(BOX* box)
{
  TIMER functionTimer(__FUNCTION__);

  VEC3F r(0.0, 0.0, 0.0);
  VEC3F point(0.0, 0.0, 0.0);
  VEC3F* velocityPointer = box->get_velocity();

  for (int z = 1; z < _zRes - 1; z++) {
    for (int y = 1; y < _yRes - 1; y++) {
      for (int x = 1; x < _xRes - 1; x++) {
        point = this->cellCenter(x, y, z);
        if ( box->inside(point) ) {
          box->update_r(point, &r);
          box->update_linearVelocity(r);
          _velocity(x, y, z) = (*velocityPointer);
        }
      }
    }
  }        
}

//////////////////////////////////////////////////////////////////////
// Builds the Modified Incomplete Cholesky Preconditioner per Bridson
//
// if heat is true, solve the heat equation, otherwise solve a 
// generic Poisson equation
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::buildA(unsigned char* skip, bool heat)
{
  TIMER functionTimer(__FUNCTION__);
	int x, y, z, index;
	Real Aoff = (heat) ? _dt * _heatDiffusion / (_dx * _dx) : 1.0;
	index = _slabSize + _xRes + 1;
	for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
		for (y = 1; y < _yRes - 1; y++, index += 2)
			for (x = 1; x < _xRes - 1; x++, index++)
			{
				// if the cell is a variable
				if (!skip[index])
				{
					_Adiag[index] = (heat) ? 1.0f : 0.0f;
					// set the matrix to the Poisson stencil in order
					if (!skip[index + _xs]) _Adiag[index] += Aoff;
					if (!skip[index - _xs]) _Adiag[index] += Aoff;
					if (!skip[index + _ys]) _Adiag[index] += Aoff;
					if (!skip[index - _ys]) _Adiag[index] += Aoff;
					if (!skip[index + _zs]) _Adiag[index] += Aoff;
					if (!skip[index - _zs]) _Adiag[index] += Aoff;
					_Aplusi[index - _xs] = (skip[index - _xs] ? 0.0 : -Aoff); 
					_Aplusi[index] = (skip[index + _xs] ? 0.0 : -Aoff); 
					_Aplusj[index - _ys] = (skip[index - _ys] ? 0.0 : -Aoff); 
					_Aplusj[index] = (skip[index + _ys] ? 0.0 : -Aoff); 
					_Aplusk[index - _zs] = (skip[index - _zs] ? 0.0 : -Aoff); 
					_Aplusk[index] = (skip[index + _zs] ? 0.0 : -Aoff); 					
				}
			}
}

void FLUID_3D_MIC::multA(const FIELD_3D& xf, FIELD_3D& bf)
{
  TIMER functionTimer(__FUNCTION__);
	int x, y, z,
	index = _slabSize + _xRes + 1;

  const int xResTwo = 2 * _xRes;

	for (z = 1; z < _zRes - 1; z++, index += xResTwo)
		for (y = 1; y < _yRes - 1; y++, index += 2)
			for (x = 1; x < _xRes - 1; x++, index++)
			{
				bf[index] = _Adiag[index] * xf[index] + 
							xf[index - _xs] * _Aplusi[index - _xs] + 
							xf[index + _xs] * _Aplusi[index] +
							xf[index - _ys] * _Aplusj[index - _ys] + 
							xf[index + _ys] * _Aplusj[index] +
							xf[index - _zs] * _Aplusk[index - _zs] + 
							xf[index + _zs] * _Aplusk[index];
			}					
}

void FLUID_3D_MIC::buildMICPreconditioner()
{
  TIMER functionTimer(__FUNCTION__);
	int x, y, z, index;
	Real tau = 0.97; //tuning constant
	Real sig = 0.25; //safety constant

	index = _slabSize + _xRes + 1;
	for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
		for (y = 1; y < _yRes - 1; y++, index += 2)
			for (x = 1; x < _xRes - 1; x++, index++)
			{
				Real pmini2 = _precon[index-_xs]*_precon[index-_xs];
				Real pminj2 = _precon[index-_ys]*_precon[index-_ys];
				Real pmink2 = _precon[index-_zs]*_precon[index-_zs];
				
				Real e0 =
				_Aplusi[index - _xs]*_Aplusi[index - _xs]*pmini2 +
				_Aplusj[index - _ys]*_Aplusj[index - _ys]*pminj2 +
				_Aplusk[index - _zs]*_Aplusk[index - _zs]*pmink2;				
				
				Real e1 = tau*(_Aplusi[index - _xs] * (_Aplusj[index - _xs] + _Aplusk[index - _xs]) * pmini2 + 
							    _Aplusj[index - _ys] * (_Aplusi[index - _ys] + _Aplusk[index - _ys]) * pminj2 +
							    _Aplusk[index - _zs] * (_Aplusi[index - _zs] + _Aplusj[index - _zs]) * pmink2);
				
				Real e = _Adiag[index] - e0 - e1;
				if (e < sig*_Adiag[index]) {e = _Adiag[index];}     //if we fall under safety
				
				if (fabs(_Adiag[index]) < 0.000001) {_precon[index] = 0;} //if there is nothing in the cell
				else              {_precon[index] = 1.0/sqrt(e);}
			}	
}

//////////////////////////////////////////////////////////////////////
// Apply MIC preoconditioner for use in conjugate gradient solve
//
// 
// 
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::applyMICPreconditioner(FIELD_3D& zf, FIELD_3D& r, FIELD_3D& q)
{
  TIMER functionTimer(__FUNCTION__);
	int x, y, z, index;
	//int i = 0;
	
	//Lq=r 
	index = _slabSize + _xRes + 1;
	for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
		for (y = 1; y < _yRes - 1; y++, index += 2)
			for (x = 1; x < _xRes - 1; x++, index++)
			{
        int xMinus = index - _xs;  
        int yMinus = index - _ys;  
        int zMinus = index - _zs;  
				Real t = r[index] - 				
        				_Aplusi[xMinus] * _precon[xMinus] * q[xMinus] -
				        _Aplusj[yMinus] * _precon[yMinus] * q[yMinus] -
        				_Aplusk[zMinus] * _precon[zMinus] * q[zMinus];				
				q[index] = t*_precon[index];
			}

	//L^tz = q
	index = _xRes*_yRes*_zRes - (_slabSize + _xRes + 2);	
	for (z = _zRes-2; z > 0; z--, index -= 2 * _xRes)
		for (y = _yRes-2; y > 0; y--, index -= 2)
			for (x = _xRes-2; x > 0; x--, index--)
			{
        Real precon = _precon[index];
				Real t = q[index] - 
        				_Aplusi[index] * precon * zf[index + _xs] -
        				_Aplusj[index] * precon * zf[index + _ys] -
        				_Aplusk[index] * precon * zf[index + _zs];
				zf[index] = t * precon;
			}
}

void FLUID_3D_MIC::zeroBoundary(FIELD_3D& x, unsigned char* skip)
{
  TIMER functionTimer(__FUNCTION__);
	//int xi,yi,zi,
	//index = _slabSize + _xRes + 1;
	for (int i = 0; i < _xRes*_yRes*_zRes; i++) {
		x[i] = skip[i] ? 0.0 : x[i];
	}
}


//////////////////////////////////////////////////////////////////////
// solve the Poisson equation with CG
//
// if heat is true, solve the heat equation, otherwise solve a 
// generic Poisson equation
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::solvePoisson(FIELD_3D& x, FIELD_3D& b, unsigned char* skip, bool heat)
{
  TIMER functionTimer(__FUNCTION__);
	buildA(skip, heat);
	buildMICPreconditioner();
	
	multA(x,_residual);

  TIMER residualTimer("Build residual");
	int xi,yi,zi,
	index = _slabSize + _xRes + 1;	
	for (zi = 1; zi < _zRes - 1; zi++, index += 2 * _xRes)
		for (yi = 1; yi < _yRes - 1; yi++, index += 2)
			for (xi = 1; xi < _xRes - 1; xi++, index++)
				_residual[index] = b[index] - _residual[index];
  zeroBoundary(_residual, skip);
  residualTimer.stop();

	applyMICPreconditioner(_z,_residual, _q);
	
  _direction = _z;

  // deltaNew = transpose(r) * r
  Real deltaNew = _residual.dot(_z);

  _solverEps = 1e-10;

  // While deltaNew > (eps^2) * delta0
  const Real eps  = _solverEps;
  Real maxR = 2.0f * eps;
  //Real maxR = _residual.max();
  int i = 0;
	
  while ((i < _iterations) && (maxR > eps))
  {
    TIMER innerTimer("PCG inner loop");
    // q = Ad
	  multA(_direction, _z);	  
	  // alpha = deltaNew / (transpose(d) * q)
    TIMER dotTimer("PCG dot");
    Real alpha = _direction.dot(_z);
    dotTimer.stop();

    // if alpha broke down, bail  
    if (fabs(alpha) > 0.0f) 
      alpha = deltaNew / alpha;
    else
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " PCG alpha broke down! Bailing." << endl;
      break;
    }

    TIMER axpyTimer("PCG axpy");
    // x = x + alpha * d
    x.axpy(alpha, _direction);
    // r = r - alpha * z
    _residual.axpy(-alpha, _z);
    axpyTimer.stop();
	
	  applyMICPreconditioner(_z,_residual, _q);
	  zeroBoundary(_q, skip);
	 
    Real oldMaxR = maxR;

    _residual.setZeroBorder();
    maxR = _residual.max();  
	  // deltaOld = deltaNew
    Real deltaOld = deltaNew;
	
    // deltaNew = transpose(r) * r
    TIMER dotTimer2("PCG dot");
    deltaNew = _residual.dot(_z);
    dotTimer2.stop();

    // beta = deltaNew / deltaOld
    TIMER directionTimer("PCG direction");
    Real beta = deltaNew / deltaOld;
	  _direction *= beta;
	  _direction += _z;
	  zeroBoundary(_direction, skip);
    directionTimer.stop();
    // i = i + 1
    i++;

    // if we didn't make any relative progress, bail
    Real relative = fabs((maxR - oldMaxR) / oldMaxR);
    if (relative < 1e-7)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Bailing! Residual before iteration: " << oldMaxR << " after: " << maxR << " relative change: " << relative << endl;
      break;
    }

    //cout << " maxR: " << maxR << " alpha: " << alpha  << " beta: " << beta << endl;
  }
  cout << "   " << i << " iterations converged to inf norm: " << maxR << " two norm: " << _residual.flattened().norm2() << endl;
}

//////////////////////////////////////////////////////////////////////
// Dump state to a file
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::writeGz(string filename) const
{
  TIMER functionTimer(__FUNCTION__);

  // make sure it's named gz
  int size = filename.size();
  if (filename[size - 1] != 'z' || filename[size - 2] != 'g')
    filename = filename + string(".gz");

  gzFile file;
  file = gzopen(filename.c_str(), "wb1"); 
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FIELD_3D write failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  cout << " Writing file " << filename.c_str() << " ..."; flush(cout);

  //int _xRes, _yRes, _zRes, _maxRes;
  gzwrite(file, (void*)&_xRes, sizeof(int));
  gzwrite(file, (void*)&_yRes, sizeof(int));
  gzwrite(file, (void*)&_zRes, sizeof(int));
  gzwrite(file, (void*)&_maxRes, sizeof(int));
	//int _xs, _ys, _zs;
  gzwrite(file, (void*)&_xs, sizeof(int));
  gzwrite(file, (void*)&_ys, sizeof(int));
  gzwrite(file, (void*)&_zs, sizeof(int));

  gzwrite(file, (void*)&_totalCells, sizeof(int));
  gzwrite(file, (void*)&_slabSize, sizeof(int));
  gzwrite(file, (void*)&_totalSteps, sizeof(int));

  gzwrite(file, (void*)&_dx, sizeof(Real));
  gzwrite(file, (void*)&_totalTime, sizeof(Real));
  gzwrite(file, (void*)&_dt, sizeof(Real));
  gzwrite(file, (void*)&_buoyancy, sizeof(Real));
  gzwrite(file, (void*)&_vorticityEps, sizeof(Real));
  gzwrite(file, (void*)&_heatDiffusion, sizeof(Real));
  gzwrite(file, (void*)&_solverEps, sizeof(Real));

  _density.writeGz(file);
  _densityOld.writeGz(file);
  _heat.writeGz(file);
  _heatOld.writeGz(file);
  _pressure.writeGz(file);
  _divergence.writeGz(file);

  _velocity.writeGz(file);
  _velocityOld.writeGz(file);
  _force.writeGz(file);
  _vorticity.writeGz(file);

  _preprojection.writeGz(file);
  _prediffusion.writeGz(file);
  _preadvection.writeGz(file);
  _cachedQuadratic.writeGz(file);
  _postadvection.writeGz(file);
  _cachedPressure.writeGz(file);
  _cachedDivergence.writeGz(file);

  gzwrite(file, (void*)_obstacles, sizeof(unsigned char) * _totalCells);

  _center.writeGz(file);
  _lengths.writeGz(file);

  // write out IOP stuff
  int howMany = _cachedPressuresIOP.size();
  gzwrite(file, (void*)&howMany, sizeof(int));
  for (int x = 0; x < howMany; x++)
  {
    _cachedPressuresIOP[x].writeGz(file);
    _cachedDivergencesIOP[x].writeGz(file);
    _cachedVelocitiesIOP[x].writeGz(file);
  }

  // write out prevorticity
  _prevorticity.writeGz(file);
  _midMaccormack.writeGz(file);

  // write out IOP stuff
  _postIOP.writeGz(file);

  if (_neumannIOP.rows() > 0)
    _neumannIOP.writeGz(file);
  gzclose(file);

  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// Read state from a zipped file
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::readGz(string filename)
{
  TIMER functionTimer(__FUNCTION__);
  // make sure it's named gz
  int size = filename.size();
  if (filename[size - 1] != 'z' || filename[size - 2] != 'g')
    filename = filename + string(".gz");

  gzFile file;
  file = gzopen(filename.c_str(), "rb"); 
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " FLUID_3D_MIC read failed! " << endl;
    cout << " Could not open file " << filename.c_str() << endl;
    exit(0);
  }

  cout << " Reading file " << filename.c_str() << " ..."; flush(cout);

  //int _xRes, _yRes, _zRes, _maxRes;
  gzread(file, (void*)&_xRes, sizeof(int));
  gzread(file, (void*)&_yRes, sizeof(int));
  gzread(file, (void*)&_zRes, sizeof(int));
  gzread(file, (void*)&_maxRes, sizeof(int));

	//int _xs, _ys, _zs;
  gzread(file, (void*)&_xs, sizeof(int));
  gzread(file, (void*)&_ys, sizeof(int));
  gzread(file, (void*)&_zs, sizeof(int));

  gzread(file, (void*)&_totalCells, sizeof(int));
  gzread(file, (void*)&_slabSize, sizeof(int));
  gzread(file, (void*)&_totalSteps, sizeof(int));

  gzread(file, (void*)&_dx, sizeof(Real));
  gzread(file, (void*)&_totalTime, sizeof(Real));
  gzread(file, (void*)&_dt, sizeof(Real));
  gzread(file, (void*)&_buoyancy, sizeof(Real));
  gzread(file, (void*)&_vorticityEps, sizeof(Real));
  gzread(file, (void*)&_heatDiffusion, sizeof(Real));
  gzread(file, (void*)&_solverEps, sizeof(Real));

  _density.readGz(file);
  _densityOld.readGz(file);
  _heat.readGz(file);
  _heatOld.readGz(file);
  _pressure.readGz(file);
  _divergence.readGz(file);
	
  _velocity.readGz(file);
  _velocityOld.readGz(file);
  _force.readGz(file);
  _vorticity.readGz(file);

  _preprojection.readGz(file);
  _prediffusion.readGz(file);
  _preadvection.readGz(file);
  _cachedQuadratic.readGz(file);
  _postadvection.readGz(file);
  _cachedPressure.readGz(file);
  _cachedDivergence.readGz(file);

  if (_obstacles) delete[] _obstacles;
	_obstacles = new unsigned char[_totalCells];

  gzread(file, (void*)_obstacles, sizeof(unsigned char) * _totalCells);

  _center.readGz(file);
  _lengths.readGz(file);

  // read in IOP stuff
  _cachedPressuresIOP.clear();
  _cachedDivergencesIOP.clear();
  _cachedVelocitiesIOP.clear();

  int howMany;
  int success = gzread(file, (void*)&howMany, sizeof(int));

  // see if any IOP snapshots are in the file
  if (success)
  {
    //cout << " Found " << howMany << " IOP snapshots ... " << flush;
    for (int x = 0; x < howMany; x++)
    {
      FIELD_3D pressureIOP;
      FIELD_3D divergenceIOP;
      VECTOR3_FIELD_3D velocityIOP;

      pressureIOP.readGz(file);
      divergenceIOP.readGz(file);
      velocityIOP.readGz(file);

      _cachedPressuresIOP.push_back(pressureIOP);
      _cachedDivergencesIOP.push_back(divergenceIOP);
      _cachedVelocitiesIOP.push_back(velocityIOP);
    }
  }

  // check if some optional fields are available
  if (!gzeof(file))
  {
    _prevorticity.readGz(file);
    _midMaccormack.readGz(file);
    if (!gzeof(file))
    {
      _postIOP.readGz(file);
      if (!gzeof(file))
        _neumannIOP.readGz(file);
    }
  }
  else
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " GZ end of file was reached! " << endl;
  }

  gzclose(file);

  // resize all PCG objects
	_precon = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_z      = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Adiag  = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aplusi = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aplusj = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aplusk = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aminui = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aminuj = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
	_Aminuk = FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);

  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// build the peeled damping matrix
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::buildPeeledDampingMatrixFull()
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Building peeled damping matrix ... ";flush(cout);
	const Real w = 0.9;
  int xPeeled = _xRes - 2;
  int yPeeled = _yRes - 2;
  int zPeeled = _zRes - 2;
  int slabPeeled = xPeeled * yPeeled;
  int peeledDims = xPeeled * yPeeled * zPeeled;
  int totalCells = 3 * peeledDims;
  const Real wMinus = 1.0 - w;
  const Real sixth = (1./ 6.)* w;

  // initially everybody is identity
  _peeledDampingFull = SPARSE_MATRIX(totalCells, totalCells);
  _peeledDampingFull.setToIdentity();

  // populate the rest
  cout << " Setting entries " << endl;
  for (int z = 0; z < zPeeled; z++) 
    for (int y = 0; y < yPeeled; y++)
      for (int x = 0; x < xPeeled; x++) 
      {
        int index = 3 * (x + y * xPeeled + z * slabPeeled);
        int left  = 3 * ((x + 1) + y * xPeeled + z * slabPeeled);
        int right = 3 * ((x - 1) + y * xPeeled + z * slabPeeled);
        int up   = 3 * (x + (y + 1) * xPeeled + z * slabPeeled);
        int down = 3 * (x + (y - 1) * xPeeled + z * slabPeeled);
        int near = 3 * (x + y * xPeeled + (z + 1) * slabPeeled);
        int far  = 3 * (x + y * xPeeled + (z - 1) * slabPeeled);

        for (int i = 0; i < 3; i++)
        {
          _peeledDampingFull(index + i, index + i) = wMinus;

          if (x != xPeeled - 1)
            _peeledDampingFull(index + i, left + i) = sixth;

          if (x != 0)
            _peeledDampingFull(index + i, right + i) = sixth;

          if (y != yPeeled - 1)
            _peeledDampingFull(index + i, up + i)    = sixth;

          if (y != 0)
            _peeledDampingFull(index + i, down + i)  = sixth;

          if (z != zPeeled - 1)
            _peeledDampingFull(index + i, near + i)  = sixth;

          if (z != 0)
            _peeledDampingFull(index + i, far + i)   = sixth;
        }
      }
}


void FLUID_3D_MIC::buildPeeledSparseIOP(SPARSE_MATRIX& A, const VEC3I& center, double radius) 
{
  A = SPARSE_MATRIX(3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2),
                    3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2));
  A.setToIdentity();
  VEC3F centerCoords = cellCenter(center[0], center[1], center[2]);
  int index = 0;
  for (int z = 1; z < _zRes - 1; z++) {
    for (int y = 1; y < _yRes - 1; y++) {
      for (int x = 1; x < _xRes - 1; x++) {
        for (int i = 0; i < 3; i++, index++) {
          if ( norm2(cellCenter(x, y, z) - centerCoords) < radius * radius ) {
          A(index, index) = 0.0;
          }
        }
      }
    }
  }        
}

void FLUID_3D_MIC::appendStreams() const
{
  TIMER functionTimer(__FUNCTION__);
  assert(_snapshotPath.size() > 0);

  cout << " Appending to streams in path " << _snapshotPath.c_str() << " ... " << flush;
  vector<int> finalRows(6);
  vector<int> finalCols(6);

  // read in the initial file dims
  vector<string> filenames;
  filenames.push_back(_snapshotPath + string("velocity.final.matrix.dims"));
  filenames.push_back(_snapshotPath + string("velocity.preproject.matrix.dims"));
  filenames.push_back(_snapshotPath + string("velocity.prediffuse.matrix.dims"));
  filenames.push_back(_snapshotPath + string("velocity.preadvect.matrix.dims"));
  filenames.push_back(_snapshotPath + string("velocity.iop.matrix.dims"));
  filenames.push_back(_snapshotPath + string("pressure.matrix.dims"));
  for (int x = 0; x < filenames.size(); x++)
  {
    FILE* dimsFile = fopen(filenames[x].c_str(), "rb");
    if (dimsFile == NULL) continue;
    int rows, cols;

    // write dimensions
    fread((void*)&rows, sizeof(int), 1, dimsFile);
    fread((void*)&cols, sizeof(int), 1, dimsFile);
    fclose(dimsFile);

    finalRows[x] = rows;
    finalCols[x] = cols;
  }

  FILE* fileFinal;
  FILE* filePreproject;
  FILE* filePrediffuse;
  FILE* filePreadvect;
  FILE* fileIOP;
  FILE* fileNeumannIOP;
  FILE* filePressure;

  string velocityFinalFile = _snapshotPath + string("velocity.final.matrix.transpose");
  string velocityPreprojectFile = _snapshotPath + string("velocity.preproject.matrix.transpose");
  string velocityPrediffuseFile = _snapshotPath + string("velocity.prediffuse.matrix.transpose");
  string velocityPreadvectFile = _snapshotPath + string("velocity.preadvect.matrix.transpose");
  string velocityIOPFile = _snapshotPath + string("velocity.iop.matrix.transpose");
  string neumannIOPFile = _snapshotPath + string("neumann.iop.matrices");
  string pressureFile = _snapshotPath + string("pressure.matrix.transpose");

  // if there's a Neumann matrix at all, it's using IOP
  bool usingIOP = (_neumannIOP.rows() > 0);
  cout << "usingIOP deduced from neumann.rows: " << usingIOP << '\n';

  fileFinal      = fopen(velocityFinalFile.c_str(), "ab");
  filePreproject = fopen(velocityPreprojectFile.c_str(), "ab");
  filePrediffuse = fopen(velocityPrediffuseFile.c_str(), "ab");
  filePreadvect  = fopen(velocityPreadvectFile.c_str(), "ab");
  fileIOP        = fopen(velocityIOPFile.c_str(), "ab");
  if (usingIOP) 
    fileNeumannIOP  = fopen(neumannIOPFile.c_str(), "ab");
  filePressure   = fopen(pressureFile.c_str(), "ab");
    
  VECTOR velocity = _velocity.peelBoundary().flattened();
  int cols = velocity.size();
  assert(cols % 3 == 0);
  int scalarCols = cols / 3;

  // set the rows and cols, in case this is the first time
  // first 5 are vector fields, but the last one is a scalar field
  for (int x = 0; x < 5; x++)
    if (finalCols[x] == 0)
      finalCols[x] = cols;
  if (finalCols[5] == 0)
    finalCols[5] = scalarCols;


  if (velocity.norm2() > 1e-7)
  {
    fwrite((void*)(velocity.data()), sizeof(double), cols, fileFinal);
    finalRows[0]++;
  }

  velocity = _preprojection.peelBoundary().flattened();
  if (velocity.norm2() > 1e-7)
  {
    fwrite((void*)(velocity.data()), sizeof(double), cols, filePreproject);
    finalRows[1]++;
  }
  
  velocity = _prediffusion.peelBoundary().flattened();
  if (velocity.norm2() > 1e-7)
  {
    fwrite((void*)(velocity.data()), sizeof(double), cols, filePrediffuse);
    finalRows[2]++;
  }
  
  velocity = _preadvection.peelBoundary().flattened();
  if (velocity.norm2() > 1e-7)
  {
    fwrite((void*)(velocity.data()), sizeof(double), cols, filePreadvect);
    finalRows[3]++;
  }

  if (usingIOP)
  {
    velocity = _postIOP.peelBoundary().flattened();
    if (velocity.norm2() > 1e-7)
    {
      fwrite((void*)(velocity.data()), sizeof(double), cols, fileIOP);
      finalRows[4]++;
    }
  }
 
  VECTOR scalar = _pressure.peelBoundary().flattened();
  if (scalar.norm2() > 1e-7)
  {
    fwrite((void*)(scalar.data()), sizeof(double), scalarCols, filePressure);
    finalRows[5]++;
  }

  scalar = _divergence.peelBoundary().flattened();
  if (scalar.norm2() > 1e-7)
  {
    fwrite((void*)(scalar.data()), sizeof(double), scalarCols, filePressure);
    finalRows[5]++;
  }
  
  if (usingIOP) 
    _neumannIOP.write(fileNeumannIOP);

  fclose(fileFinal);
  fclose(filePreproject);
  fclose(filePrediffuse);
  fclose(filePreadvect);
  if (usingIOP)
    fclose(fileNeumannIOP);
  fclose(filePressure);

  // write out the dimensions of the matrix to a different file
  for (int x = 0; x < filenames.size(); x++)
  {
    FILE* dimsFile = fopen(filenames[x].c_str(), "wb");

    // write dimensions
    fwrite((void*)&finalRows[x], sizeof(int), 1, dimsFile);
    fwrite((void*)&finalCols[x], sizeof(int), 1, dimsFile);
    fclose(dimsFile);
  }
  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// Output directly to the files that will be sent to SVD. Uses IOP
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::appendStreamsIOP() const
{
  TIMER functionTimer(__FUNCTION__);
  assert(_snapshotPath.size() > 0);

  cout << " Appending to streams in path " << _snapshotPath.c_str() << " ... " << flush;
  vector<int> finalRows(6);
  vector<int> finalCols(6);

  // read in the initial file dims
  vector<string> filenames;
  filenames.push_back(_snapshotPath + string("velocity.final.matrix.dims"));
  filenames.push_back(_snapshotPath + string("velocity.preproject.matrix.dims"));
  filenames.push_back(_snapshotPath + string("velocity.prediffuse.matrix.dims"));
  filenames.push_back(_snapshotPath + string("velocity.preadvect.matrix.dims"));
  filenames.push_back(_snapshotPath + string("velocity.iop.matrix.dims"));
  filenames.push_back(_snapshotPath + string("pressure.matrix.dims"));
  for (int x = 0; x < filenames.size(); x++)
  {
    FILE* dimsFile = fopen(filenames[x].c_str(), "rb");
    if (dimsFile == NULL) continue;
    int rows, cols;

    // write dimensions
    fread((void*)&rows, sizeof(int), 1, dimsFile);
    fread((void*)&cols, sizeof(int), 1, dimsFile);
    fclose(dimsFile);

    finalRows[x] = rows;
    finalCols[x] = cols;
  }

  FILE* fileFinal;
  FILE* filePreproject;
  FILE* filePrediffuse;
  FILE* filePreadvect;
  FILE* fileIOP;
  FILE* fileNeumannIOP;
  FILE* filePressure;

  string velocityFinalFile = _snapshotPath + string("velocity.final.matrix.transpose");
  string velocityPreprojectFile = _snapshotPath + string("velocity.preproject.matrix.transpose");
  string velocityPrediffuseFile = _snapshotPath + string("velocity.prediffuse.matrix.transpose");
  string velocityPreadvectFile = _snapshotPath + string("velocity.preadvect.matrix.transpose");
  string velocityIOPFile = _snapshotPath + string("velocity.iop.matrix.transpose");
  string neumannIOPFile = _snapshotPath + string("neumann.iop.matrices");
  string pressureFile = _snapshotPath + string("pressure.matrix.transpose");

  bool usingIOP = 1; 

  fileFinal      = fopen(velocityFinalFile.c_str(), "ab");
  filePreproject = fopen(velocityPreprojectFile.c_str(), "ab");
  filePrediffuse = fopen(velocityPrediffuseFile.c_str(), "ab");
  filePreadvect  = fopen(velocityPreadvectFile.c_str(), "ab");
  fileIOP        = fopen(velocityIOPFile.c_str(), "ab");
  if (usingIOP)
    fileNeumannIOP = fopen(neumannIOPFile.c_str(), "ab");
  filePressure   = fopen(pressureFile.c_str(), "ab");
    
  VECTOR velocity = _velocity.peelBoundary().flattened();
  int cols = velocity.size();
  assert(cols % 3 == 0);
  int scalarCols = cols / 3;

  // set the rows and cols, in case this is the first time
  for (int x = 0; x < 5; x++)
    if (finalCols[x] == 0)
      finalCols[x] = cols;
  if (finalCols[5] == 0)
    finalCols[5] = scalarCols;

  if (velocity.norm2() > 1e-7)
  {
    fwrite((void*)(velocity.data()), sizeof(double), cols, fileFinal);
    finalRows[0]++;
    cout << " final rows is now: " << finalRows[0] << endl;
  }

  velocity = _preprojection.peelBoundary().flattened();
  cout << " preprojection norm2: " << velocity.norm2() << endl;
  
  // ADJ: the condition that the norm is greater than 1e-7 was skipping the first step in
  // preprojection, prediffusion, and preadvection, so it was removed.
  // QUESTION: This seems to break the SVD call, though.
  // Not sure why a row of zeros makes it not converge.

  // if (velocity.norm2() > 1e-7)
  {
    fwrite((void*)(velocity.data()), sizeof(double), cols, filePreproject);
    finalRows[1]++;
    cout << " preproject rows is now: " << finalRows[1] << endl;
  }
  
  velocity = _prediffusion.peelBoundary().flattened();
  cout << " prediffusion norm2: " << velocity.norm2() << endl;
  // if (velocity.norm2() > 1e-7)
  {
    fwrite((void*)(velocity.data()), sizeof(double), cols, filePrediffuse);
    finalRows[2]++;
    cout << " prediffuse rows is now: " << finalRows[2] << endl;
  }
  
  velocity = _preadvection.peelBoundary().flattened();
  cout << " preadvection norm2: " << velocity.norm2() << endl;
  // if (velocity.norm2() > 1e-7)
  {
    fwrite((void*)(velocity.data()), sizeof(double), cols, filePreadvect);
    finalRows[3]++;
    cout << " preadvect rows is now: " << finalRows[3] << endl;
  }
 
  if (usingIOP)
  { 
    velocity = _postIOP.peelBoundary().flattened();
    if (velocity.norm2() > 1e-7)
    {
      fwrite((void*)(velocity.data()), sizeof(double), cols, fileIOP);
      finalRows[4]++;
      cout << " IOP rows is now: " << finalRows[4] << endl;
    }
  }

  VECTOR scalar = _pressure.peelBoundary().flattened();
  if (scalar.norm2() > 1e-7)
  {
    fwrite((void*)(scalar.data()), sizeof(double), scalarCols, filePressure);
    finalRows[5]++;
  }

  scalar = _divergence.peelBoundary().flattened();
  if (scalar.norm2() > 1e-7)
  {
    fwrite((void*)(scalar.data()), sizeof(double), scalarCols, filePressure);
    finalRows[5]++;
    cout << " pressure rows is now: " << finalRows[5] << endl;
  }

  if (usingIOP)
    _neumannIOP.write(fileNeumannIOP);

  fclose(fileFinal);
  fclose(filePreproject);
  fclose(filePrediffuse);
  fclose(filePreadvect);
  fclose(fileIOP);
  if (usingIOP)
    fclose(fileNeumannIOP);
  fclose(filePressure);

  // write out the dimensions of the matrix to a different file
  for (int x = 0; x < filenames.size(); x++)
  {
    FILE* dimsFile = fopen(filenames[x].c_str(), "wb");

    // write dimensions
    fwrite((void*)&finalRows[x], sizeof(int), 1, dimsFile);
    fwrite((void*)&finalCols[x], sizeof(int), 1, dimsFile);
    fclose(dimsFile);
  }
  cout << " done. " << endl;
}

// set the neumannIOP obstacle matrix
void FLUID_3D_MIC::setPeeledSparseMovingIOP(BOX* box) 

{
  TIMER functionTimer(__FUNCTION__);
  // if it's the first time this function is called, resize _neumannIOP.
  // it has an extra homoegenous coordinate column
  if (_neumannIOP.cols() <= 0) {
     _neumannIOP = SPARSE_MATRIX(3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2),
                    1 + 3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2));
  }
  // setToIdentity still works on rectangular (non-square) matrices
  _neumannIOP.setToIdentity();
  int index = 0;
  int lastCol = _neumannIOP.cols() - 1;
  VEC3F r(0.0, 0.0, 0.0);
  VEC3F point(0.0, 0.0, 0.0);
  VEC3F* velocityPointer = box->get_velocity();

  // rotation matrix changes per step and must be updated
  // before calling box->inside
  // box->update_rotationMatrix();

  for (int z = 1; z < _zRes - 1; z++) {
    for (int y = 1; y < _yRes - 1; y++) {
      for (int x = 1; x < _xRes - 1; x++) {
        for (int i = 0; i < 3; i++, index++) {
          point = this->cellCenter(x, y, z);
          // if point is in the interior of the obstacle
          if ( box->inside(point) ) {
            _neumannIOP(index, index) = 0.0;
            // update the radial and linear velocities
            box->update_r(point, &r);
            box->update_linearVelocity(r);
            if (i == 0) {
              _neumannIOP(index, lastCol) = (*velocityPointer)[0];
            }
            else if (i == 1) {
              _neumannIOP(index, lastCol) = (*velocityPointer)[1];
            }
            else {
              _neumannIOP(index, lastCol) = (*velocityPointer)[2];
            }
          }
        }
      }
    }
  }        
}

// set the Neumann velocity from the velocity by appending a 1
void FLUID_3D_MIC::setVelocityNeumann()
{
  int peeledDims = 3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2);
  _velocityNeumann.resize(1 + peeledDims);
  _velocityNeumann.head(peeledDims) = _velocity.peelBoundary().flattenedEigen();
  _velocityNeumann[peeledDims] = 1.0;
}

// compute the spatial coordinates from integer indexing
VEC3F FLUID_3D_MIC::cellCenter(int x, int y, int z)
{
  TIMER functionTimer(__FUNCTION__);
  // assuming lengths are all 1.0
  VEC3F halfLengths(0.5, 0.5, 0.5); 

  // set it to the lower corner
  // **********************************************************************
  // ADJ: commenting this out for now. might break the static obstacle code
  // probably the variable _center was not (0.5, 0.5, 0.5) as expected?
  // VEC3F final = _center - halfLengths;
  // **********************************************************************

  VEC3F final(0.0, 0.0, 0.0);

  double dx = 1.0 / _xRes;
  double dy = 1.0 / _yRes;
  double dz = 1.0 / _zRes;
  // displace to the NNN corner
  final[0] += x * dx;
  final[1] += y * dy;
  final[2] += z * dz;

  // displace it to the cell center
  final[0] += dx * 0.5;
  final[1] += dy * 0.5;
  final[2] += dz * 0.5;

  return final;
}

