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
#include <fstream>

#include "FLUID_3D_MIC.h"
// #include <glog/logging.h>
// #include <zlib.h>

#include "util/write_density_pbrt.h"
#include "util/stringprintf.h"
#include "3D/drawer_3d.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FLUID_3D_MIC::FLUID_3D_MIC(int xRes, int yRes, int zRes, int amplify,const bool use_two_phase_smoke,
                           const double added_smoke_density, const std::string& buoyancy_dir,
                           const Real buoyancy,const double density_attenuate_factor,
                          const bool attenuate_smoke, const std::string density_folder, const double dt,
                           unsigned int* boundaries) :
	_xRes(xRes), _yRes(yRes), _zRes(zRes), _res(0.), use_two_phase_smoke_(use_two_phase_smoke),
	added_smoke_density_(added_smoke_density), _buoyancy(buoyancy),
	density_attenuate_factor_(density_attenuate_factor), attenuate_smoke_(attenuate_smoke), _dt(dt),
    density_folder_(density_folder)
{
	// set simulation consts
  // _dt = 0.60; 

  _iterations = 1000;
  _heatDiffusion = 1e-3;
  _vorticityEps = 1.0;
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

    SetSideObstacles();
    
    cell_center_[0] = static_cast<float>(xRes)*_dx*0.5; cell_center_[1] = static_cast<float>(yRes)*_dx*0.5;
    cell_center_[2] = static_cast<float>(zRes)*_dx*0.5;
    if (buoyancy_dir == "x") {
       buoyancy_dir_ = VEC3F(1, 0, 0);
    } else if (buoyancy_dir == "y") {
       buoyancy_dir_ = VEC3F(0, 1, 0);
    } else if (buoyancy_dir == "z") {
       buoyancy_dir_ = VEC3F(0, 0, 1);
    } else {
       std::cout << "FLUID_3D_MIC.cpp " << __LINE__ << " FATAL: " <<  "Invalid buoyancy_dir: " << buoyancy_dir << std::endl; exit(0);
    }
    if (use_two_phase_smoke_) {
      density_sum_ = FIELD_3D(_xRes, _yRes, _zRes);
    }
    // 
    temp1 = FIELD_3D(_xRes, _yRes, _zRes);
    temp2 = FIELD_3D(_xRes, _yRes, _zRes);
    vtemp1 = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes);
    vtemp2 = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes);
  // Particles.
  for (int i = 0; i < num_particles_; i++) {
    Particle3D particle;
    particles_.push_back(particle);
  }
    // Pointers.
  drawer_.reset(new Drawer3D());
  
  ReSeedParticles();
  render_initialized_ = false;
}

void FLUID_3D_MIC::SetSideObstacles() {
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
  Real goalTime = 0.06;
  Real currentTime = 0;

  if (use_obstacles_) {
    obstacles_wrapper_.get()->RasterToGrid(_obstacles);
    SetSideObstacles();
  }
  // compute the CFL condition

  // wipe forces
  _force.clear();

  // wipe boundaries
  _velocity.setZeroBorder();
  _density.setZeroBorder();

  // compute the forces
  // addBuoyancy(density().data());
  AddBuoyancy();
  
  _velocity.axpy(_dt, _force);
  _force.clear();

  _prevorticity = _velocity;

  addVorticity();
  _velocity.axpy(_dt, _force);

  // advect everything
  advectStam();
AdvectParticles();
  _prediffusion = _velocity;


  _preprojection = _velocity;

  // run the solvers
  project();

  currentTime += _dt;

	_totalTime += goalTime;
	_totalSteps++;
    
  if (attenuate_smoke_) {
    AttenuateSmoke();
  }
  
  // Move the obstacle.
  if (use_obstacles_ && move_obstacle_) {
    obstacles_wrapper_.get()->MoveObstacle(_dt);
  }
}

//////////////////////////////////////////////////////////////////////
// project into divergence free field
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::project()
{

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
                
				if(_obstacles[index+1]) {
                  VEC3F obs_velo(0,0,0);
                  if (use_obstacles_ && move_obstacle_) {
                    VEC3F pnt(static_cast<float>(x + 1)*_dx, static_cast<float>(y)*_dx,
                              static_cast<float>(z)*_dx);
                    obs_velo = obstacles_wrapper_.get()->getVelocityAt(pnt);
                  }
                  xright = - ( _velocity[index][0] - obs_velo[0]);
                }
				if(_obstacles[index-1]) {
                  VEC3F obs_velo(0,0,0);
                  if (use_obstacles_ && move_obstacle_) {
                    VEC3F pnt(static_cast<float>(x - 1)*_dx, static_cast<float>(y)*_dx,
                              static_cast<float>(z)*_dx);
                    obs_velo = obstacles_wrapper_.get()->getVelocityAt(pnt);
                  }
                  xleft  = - (_velocity[index][0] - obs_velo[0]);
                }
				if(_obstacles[index+_xRes]) {
                  VEC3F obs_velo(0,0,0);
                  if (use_obstacles_ && move_obstacle_) {
                    VEC3F pnt(static_cast<float>(x)*_dx, static_cast<float>(y + 1)*_dx,
                              static_cast<float>(z)*_dx);
                    obs_velo = obstacles_wrapper_.get()->getVelocityAt(pnt);
                  }
                  yup    = - (_velocity[index][1] - obs_velo[1]);
                }
                if(_obstacles[index-_xRes]) {
                  VEC3F obs_velo(0,0,0);
                  if (use_obstacles_ && move_obstacle_) {
                    VEC3F pnt(static_cast<float>(x)*_dx, static_cast<float>(y - 1)*_dx,
                              static_cast<float>(z)*_dx);
                    obs_velo = obstacles_wrapper_.get()->getVelocityAt(pnt);
                  }
                  ydown  = - (_velocity[index][1] - obs_velo[1]);
                }
				if(_obstacles[index+_slabSize]) {
                  VEC3F obs_velo(0,0,0);
                  if (use_obstacles_ && move_obstacle_) {
                    VEC3F pnt(static_cast<float>(x)*_dx, static_cast<float>(y)*_dx,
                              static_cast<float>(z + 1)*_dx);
                    obs_velo = obstacles_wrapper_.get()->getVelocityAt(pnt);
                  }
                  ztop    = - (_velocity[index][2] - obs_velo[2]);
                }
				if(_obstacles[index-_slabSize]) {
                  VEC3F obs_velo(0,0,0);
                  if (use_obstacles_ && move_obstacle_) {
                    VEC3F pnt(static_cast<float>(x)*_dx, static_cast<float>(y)*_dx,
                              static_cast<float>(z - 1)*_dx);
                    obs_velo = obstacles_wrapper_.get()->getVelocityAt(pnt);
                  }
                  zbottom = - (_velocity[index][2] - obs_velo[2]);
                }

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
              if (_obstacles[index]) {
                
                VEC3F obs_velo(0,0,0);
                if (use_obstacles_ && move_obstacle_) {
                  VEC3F pnt(static_cast<float>(x)*_dx, static_cast<float>(y)*_dx,
                            static_cast<float>(z)*_dx);
                  obs_velo = obstacles_wrapper_.get()->getVelocityAt(pnt);
                }
                _velocity[index] = obs_velo;
                
                continue;
              }
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
                
                VEC3F obs_velo(0,0,0);
                if (use_obstacles_ && move_obstacle_) {
                  VEC3F pnt(static_cast<float>(x)*_dx, static_cast<float>(y)*_dx,
                            static_cast<float>(z)*_dx);
                  obs_velo = obstacles_wrapper_.get()->getVelocityAt(pnt);
                }
                _velocity[index] = obs_velo;
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

	int index = 0;

	Real beta = _buoyancy;
	if(beta==0.) return;
	for (int z = 0; z < _zRes; z++)
		for (int y = 0; y < _yRes; y++)
			for (int x = 0; x < _xRes; x++, index++) 
        _force[index][1] += beta * field[index];
}

//////////////////////////////////////////////////////////////////////
// add vorticity to the force field
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::addVorticity()
{

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


	const Real dt0 = _dt / _dx;
  // VECTOR3_FIELD_3D::advect(dt0, _velocityOld, _densityOld, _density);
  VECTOR3_FIELD_3D::advectMacCormack(dt0, _velocityOld, _densityOld, _density, temp1, temp2);
  
  if (use_two_phase_smoke_) {
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].SwapPointer();
    }
    for (int i = 0; i < denWarppers_.size(); i++) {
      VECTOR3_FIELD_3D::advectMacCormack(dt0, _velocityOld, denWarppers_[i].density_old_,
                                               denWarppers_[i].density_, temp1, temp2);
    }
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].density_.setZeroBorder();
    }
  }  
  
   VECTOR3_FIELD_3D::advectMacCormack(dt0, _velocityOld, _velocityOld, _velocity, vtemp1, vtemp2);
// VECTOR3_FIELD_3D::advect(dt0, _velocityOld, _velocityOld, _velocity);
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
}

//////////////////////////////////////////////////////////////////////
// add a test cube of density to the center
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::addSmokeColumn() 
{
	addSmokeTestCase(_density.data(), _res);
}

//////////////////////////////////////////////////////////////////////
// generic static version, so that it can be applied to the
// WTURBULENCE grid as well
//////////////////////////////////////////////////////////////////////
void FLUID_3D_MIC::addSmokeTestCase(Real* field, VEC3I res)
{
  // VEC3I res(_xRes, _yRes, _zRes);
  const int slabSize = res[0]*res[1]; 
  int maxRes = res.maxElement();
  Real dx = 1.0f / (Real)maxRes;

  Real yTotal = dx * res[1];
  Real zTotal = dx * res[2];

  Real heighMin = 0.05;
  Real heighMax = 0.10;

  for (int z = 0; z < res[2]; z++)
    for (int y = 0; y < res[1]; y++)
      for ( int x = (int)(heighMin*res[0]); x <= (int)(heighMax * res[0]); x++)
      {
        Real yLength = y * dx - yTotal * 0.48f;
        Real zLength = z * dx - zTotal * 0.52f;
        Real radius = sqrtf(yLength * yLength + zLength * zLength);

        if (radius < 0.075f * yTotal)
        {
          int index = x + y * res[0] + z * slabSize;
          _density[index] = added_smoke_density_;
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

	buildA(skip, heat);
	buildMICPreconditioner();
	
	multA(x,_residual);


	int xi,yi,zi,
	index = _slabSize + _xRes + 1;	
	for (zi = 1; zi < _zRes - 1; zi++, index += 2 * _xRes)
		for (yi = 1; yi < _yRes - 1; yi++, index += 2)
			for (xi = 1; xi < _xRes - 1; xi++, index++)
				_residual[index] = b[index] - _residual[index];
  zeroBoundary(_residual, skip);

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

    // q = Ad
	  multA(_direction, _z);	  
	  // alpha = deltaNew / (transpose(d) * q)

    Real alpha = _direction.dot(_z);


    // if alpha broke down, bail  
    if (fabs(alpha) > 0.0f) 
      alpha = deltaNew / alpha;
    else
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " PCG alpha broke down! Bailing." << endl;
      break;
    }


    // x = x + alpha * d
    x.axpy(alpha, _direction);
    // r = r - alpha * z
    _residual.axpy(-alpha, _z);

	
	  applyMICPreconditioner(_z,_residual, _q);
	  zeroBoundary(_q, skip);
	 
    Real oldMaxR = maxR;

    _residual.setZeroBorder();
    maxR = _residual.max();  
	  // deltaOld = deltaNew
    Real deltaOld = deltaNew;
	
    // deltaNew = transpose(r) * r
 
    deltaNew = _residual.dot(_z);


    // beta = deltaNew / deltaOld

    Real beta = deltaNew / deltaOld;
	  _direction *= beta;
	  _direction += _z;
	  zeroBoundary(_direction, skip);

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

void FLUID_3D_MIC::OutputSmokeToFolder() {
  std::string fname = StringPrintf("%s%04d.pbrt", density_folder_.c_str(), _totalSteps);
  std::cout <<  "Output the smoke to file: " << fname << std::endl;
  if (!use_two_phase_smoke_) {
    WriteDensityPBRT(fname, _density, 1.0);
  } else {
    density_sum_ = _density; // + density1_;
    for (int i = 0; i < denWarppers_.size(); i++) {
      density_sum_ += denWarppers_[i].density_;
    }
    WriteDensityPBRT(fname, density_sum_, 1.0);
  }
 if (use_obstacles_) {
    obstacles_wrapper_.get()->WriteObstacleToPBRT(density_folder_.c_str(), _totalSteps);
  }
}

void FLUID_3D_MIC::AttenuateSmoke() {
  for (uint i = 0; i < _totalCells; i++) {
    _density[i] *= density_attenuate_factor_;
    if (_density[i] < 0.003) {
      _density[i] = 0;
    }
  }
  if (use_two_phase_smoke_) {
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].AttenuateSmoke(density_attenuate_factor_);
    }
  }
}

void FLUID_3D_MIC::DrawSmoke() {
  
  if (!use_two_phase_smoke_) {
    _density.draw(0.3);
  } else {
    density_sum_ = _density;
    for (int i = 0; i < denWarppers_.size(); i++) {
      density_sum_ += denWarppers_[i].density_;
    }
    density_sum_.draw(0.3);
  }
}

void FLUID_3D_MIC::AddBuoyancy() {

  VEC3F dir_ = buoyancy_dir_;
  for (uint i = 0; i < _totalCells; i++) {
    _force[i][0] += _buoyancy*_density[i]*dir_[0];
    _force[i][1] += _buoyancy*_density[i]*dir_[1];
    _force[i][2] += _buoyancy*_density[i]*dir_[2];
  }
  
  if (use_two_phase_smoke_) {
    // point upwards.
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].AddBuoyancyField(&_force);
    }
  }
}

#define EXPECT_STRING(str) \
in >> temp; \
if (temp != str) { \
  std::cout << "FLUID_3D_MIC.cpp " << __LINE__ << " FATAL: " <<  "Error: " << temp << std::endl; exit(0); \
}\
temp.clear();
void FLUID_3D_MIC::ParseSourceFromFile(const std::string& fname) {
  std::ifstream in(fname);
  if (!in.is_open()) {
     std::cout << "FLUID_3D_MIC.cpp " << __LINE__ << " FATAL: " <<  "Cannot open: " << fname << std::endl; exit(0);
  }
  do {
    VEC3F position, ssize, buoyancy_dir;
    std::string temp;
    EXPECT_STRING("position")
    in >> position[0]; in >> position[1]; in >> position[2];
    EXPECT_STRING("size")
    in >> ssize[0]; in >> ssize[1]; in >> ssize[2];
    EXPECT_STRING("buoyancy_direction")
    in >> buoyancy_dir[0]; in >> buoyancy_dir[1]; in >> buoyancy_dir[2];
    
    DensityWarpper denWarpper(VEC3I(_xRes, _yRes, _zRes),
                                          buoyancy_dir, _buoyancy, added_smoke_density_);
    denWarpper.InitialzeSourceSomke(position, ssize, _maxRes);
    denWarppers_.push_back(denWarpper);
  } while(!in.eof());
}

#undef EXPECT_STRING

void FLUID_3D_MIC::AddDensity(const int xpos, const int ypos, const int zpos, 
                              const int length, const int width, const int height,
                              FIELD_3D* field) {
  int xbegin = xpos - length / 2;
  int xend = xpos + length / 2;
  int ybegin = ypos - width / 2;
  int yend = ypos + width / 2;
  int zbegin = zpos - height / 2;
  int zend = zpos + height / 2;
  
  // clamp.
  xbegin = (xbegin < 0) ? 0 : xbegin;
  xend = (xend < 0) ? 0 : xend;
  ybegin = (ybegin < 0) ? 0 : ybegin;
  yend = (yend < 0) ? 0 : yend;
  zbegin = (zbegin < 0 ) ? 0 : zbegin;
  zend = (zend < 0 ) ? 0 : zend;
  
  xbegin = (xbegin > _xRes - 1) ? _xRes - 1 : xbegin;
  xend = (xend > _xRes - 1) ? _xRes - 1 : xend;
  ybegin = (ybegin > _yRes - 1) ? _yRes - 1 : ybegin;
  yend = (yend > _yRes - 1) ? _yRes - 1 : yend;
  zbegin = (zbegin > _zRes - 1) ? _zRes - 1 : zbegin;
  zend = (zend > _zRes - 1) ? _zRes - 1 : zend;
  
  uint idx_begin = xbegin + ybegin*_xRes + zbegin*_slabSize;
  
  for (int k = zbegin; k <= zend; k++) {
    for (int j = ybegin; j <= yend; j++) {
      for (int i = xbegin; i <= xend; i++) {
        int ind = idx_begin + i + j*_xRes + k*_slabSize;
        (*field)(i,j,k) = added_smoke_density_;
      }
    }
  }
}

void FLUID_3D_MIC::AddSmokeTestCase() {
  VEC3F pointTransformed = source_pos_flt_ - cell_center_;
  pointTransformed += cell_center_;
  souce_pos_[0] = static_cast<int>(pointTransformed[0]*_maxRes);
  souce_pos_[1] = static_cast<int>(pointTransformed[1]*_maxRes);
  souce_pos_[2] = static_cast<int>(pointTransformed[2]*_maxRes);
  
  AddDensity(souce_pos_[0], souce_pos_[1], souce_pos_[2],
             souce_size_[0], souce_size_[1], souce_size_[2], &_density);
  if (use_two_phase_smoke_) {
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].AddSmokeTestCase();
    }
  }
}

void FLUID_3D_MIC::InitializeObstacles(const ObstacleParams3D& param_) {
  std::cout <<  "Initialize the obstacle..." << std::endl;
  if (!param_.use_obstacles) {
    use_obstacles_ = false;
    return;
  } else  {
    // obstacle_.reset(new Obstacle2D(xRes_, yRes_, param_.obstacle_force_scale));
    // obstacle_.get()->InitializeUsingImage(param_.img_file_name);
    use_obstacles_ = true;
    handle_obstacle_implicit_ = param_.handle_obstacle_implicit;
    // force_amp_ = param_.obstacle_force_scale;
    move_obstacle_ = param_.move_obstacle;
  }
  obstacles_wrapper_.reset(new ObstacleWrapper3D(_dx, param_.obstacle_type,
                param_.obstacle_file, param_.obstacle_list, VEC3I(_xRes, _yRes, _zRes)));
  /*if (handle_obstacle_implicit_) {
    obs_proj_mat_.resize(basis_dim_, basis_dim_);
    std::cout <<  "start to compute the obstacle projection matrix..." << std::endl;
    timer_.Reset();
    obstacles_wrapper_.get()->ComputeObstacleProjectionMatrix(
      (*basis_set_.get()), VEC3I(xRes_, yRes_, zRes_),
      dx_, force_amp_, &obs_proj_mat_);
    // multiply force_amp_ and add I.
    obs_proj_mat_ *= force_amp_;
    for (int i = 0; i < basis_dim_; i++) {
      obs_proj_mat_(i,i) += 1.0;
    }
    std::cout << "Time spend to compute the obstacle matrix: "
        << timer_.ElapsedTimeInSeconds();
  }*/
}

void FLUID_3D_MIC::DrawObstacles() {
  if (use_obstacles_) {
    obstacles_wrapper_.get()->Draw();
  }
}

void FLUID_3D_MIC::DrawParticles(const double ptl_length) {
  drawer_.get()->DrawParticles(particles_, _dx, _dt, ptl_length);
}

void FLUID_3D_MIC::AdvectParticles() {
  const double dt0 = _dt / _dx;
  for (int i = 0; i < num_particles_; i++) {
    // Get the velocity at the particles.
    const VEC3& position = particles_[i].position;
    // This function assumen the velocity is stored on the vertices of the grid, while for laplacian fluid
    // case, the velocity is stored on the center of the grid.
    VEC3 p_v = _velocity.GetVelocity(position[0]-0.5, position[1]-0.5, position[2]-0.5);
    // Forward Eular.
    particles_[i].position[0] += p_v[0] * dt0;
    particles_[i].position[1] += p_v[1] * dt0;
    particles_[i].position[2] += p_v[2] * dt0;
    particles_[i].velocity = p_v;
    if (particles_[i].position[0] < 0. || particles_[i].position[0] > _xRes|| 
        particles_[i].position[1] < 0. || particles_[i].position[1] > _yRes ||
        particles_[i].position[2] < 0. || particles_[i].position[2] > _zRes) {
      particles_[i].position[0] = std::rand() % _xRes + 0.5;
      particles_[i].position[1] = std::rand() % _yRes + 0.5;
      particles_[i].position[2] = std::rand() % _zRes + 0.5;
      particles_[i].velocity = 0.;
    }
  }
}

void FLUID_3D_MIC::ReSeedParticles() {
  // Reseed the particles at random position.
  for (int i = 0; i < num_particles_; i++) {
    int px = std::rand() % _xRes;
    int py = std::rand() % _yRes;
    int pz = std::rand() % _zRes;
    particles_[i].position[0] = px + 0.5;
    particles_[i].position[1] = py + 0.5;
    particles_[i].position[2] = pz + 0.5;
    particles_[i].velocity = 0;
  }
}

void FLUID_3D_MIC::InitializeSourceSmoke(const VEC3F& pos, const VEC3F& ssize,
                                             const std::string& source_file) {
  source_pos_flt_[0] = pos[0];
  source_pos_flt_[1] = pos[1];
  source_pos_flt_[2] = pos[2];
  souce_size_[0] = static_cast<int>(ssize[0]*_maxRes);
  souce_size_[1] = static_cast<int>(ssize[1]*_maxRes);
  souce_size_[2] = static_cast<int>(ssize[2]*_maxRes);

  if (source_file.size() > 0) {
    ParseSourceFromFile(source_file);
  }
  if (use_two_phase_smoke_) {
    AddSmokeTestCase();
  }
}
