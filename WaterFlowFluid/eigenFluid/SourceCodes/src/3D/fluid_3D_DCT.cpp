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
 
 Copyright 2018 Qiaodong Cui (qiaodong@ucsb.edu)
 */

#include <fftw3.h>
#include <fstream>

#include "3D/fluid_3D_DCT.h"
// #include <glog/logging.h>

#include "util/write_density_pbrt.h"
#include "util/stringprintf.h"
#include "3D/dirichlet_basis_3d.h"
#include "3D/drawer_3d.h"
#include "util/timer.h"

FLUID_3D_DCT::FLUID_3D_DCT(int xRes, int yRes, int zRes, int amplify,const bool use_two_phase_smoke,
                           const double added_smoke_density, const std::string& buoyancy_dir,
                           const Real buoyancy,const double density_attenuate_factor,
                          const bool attenuate_smoke, const std::string density_folder, const double dt,
                           const int buoyancy_step,
                           unsigned int* boundaries) :
	_xRes(xRes), _yRes(yRes), _zRes(zRes), _res(0.), use_two_phase_smoke_(use_two_phase_smoke),
	added_smoke_density_(added_smoke_density), _buoyancy(buoyancy),
	density_attenuate_factor_(density_attenuate_factor), attenuate_smoke_(attenuate_smoke), _dt(dt),
	buoyancy_step_(buoyancy_step),
    density_folder_(density_folder)
{
	// set simulation consts
  // _dt = 0.60; 

  _iterations = 1000;
  _heatDiffusion = 1e-3;
  _vorticityEps = .0;
	_totalTime = 0.0f;
	_totalSteps = 0;
	_res = VEC3I(_xRes,_yRes,_zRes);
	_maxRes = _res.maxElement();
   invxRes_ = 1.0 / static_cast<double>(_xRes);
   invyRes_ = 1.0 / static_cast<double>(_yRes);
   invzRes_ = 1.0 / static_cast<double>(_zRes);
   
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

  _velocity = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _velocityOld = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);

  vxFreq_ = (Real*) malloc(sizeof(Real)*_xRes*_yRes*_zRes);
  memset(vxFreq_, 0x00, sizeof(Real)*_xRes*_yRes*_zRes);
  vyFreq_ = (Real*) malloc(sizeof(Real)*_xRes*_yRes*_zRes);
  memset(vyFreq_, 0x00, sizeof(Real)*_xRes*_yRes*_zRes);
  vzFreq_ = (Real*) malloc(sizeof(Real)*_xRes*_yRes*_zRes);
  memset(vzFreq_, 0x00, sizeof(Real)*_xRes*_yRes*_zRes);
  
  _force = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);
  _vorticity = VECTOR3_FIELD_3D(_xRes, _yRes, _zRes, _center, _lengths);

	
	// allocate arrays
	_totalCells   = _xRes * _yRes * _zRes;
    invTotalSize_ = 1.0 / static_cast<double>(_totalCells);
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

    cell_center_[0] = static_cast<float>(xRes)*_dx*0.5; cell_center_[1] = static_cast<float>(yRes)*_dx*0.5;
    cell_center_[2] = static_cast<float>(zRes)*_dx*0.5;
    if (buoyancy_dir == "x") {
       buoyancy_dir_ = VEC3F(1, 0, 0);
    } else if (buoyancy_dir == "y") {
       buoyancy_dir_ = VEC3F(0, 1, 0);
    } else if (buoyancy_dir == "z") {
       buoyancy_dir_ = VEC3F(0, 0, 1);
    } else {
       std::cout << "fluid_3D_DCT.cpp " << __LINE__ << " FATAL: " <<  "Invalid buoyancy_dir: " << buoyancy_dir << std::endl; exit(0);
    }
    if (use_two_phase_smoke_) {
      density_sum_ = FIELD_3D(_xRes, _yRes, _zRes);
    }
    // 
    temp1 = FIELD_3D(_xRes, _yRes, _zRes);
    temp2 = FIELD_3D(_xRes, _yRes, _zRes);
  // Particles.
  for (int i = 0; i < num_particles_; i++) {
    Particle3D particle;
    particles_.push_back(particle);
  }
    // Pointers.
  drawer_.reset(new Drawer3D());
  
  ReSeedParticles(); 
  render_initialized_ = false;
  InitializeFFTW();
  totalDCTtime_ = 0;
}

FLUID_3D_DCT::~FLUID_3D_DCT() {
	if (_obstacles) delete[] _obstacles;
	if (_dirichletObstacleField) delete[] _dirichletObstacleField;
	if (_neumannObstacleField) delete[] _neumannObstacleField;
    free(vxFreq_);
    free(vyFreq_);
    free(vzFreq_);
}

void FLUID_3D_DCT::step()
{
  Real goalTime = 0.06;
  Real currentTime = 0;

  if (use_obstacles_) {
    obstacles_wrapper_.get()->RasterToGrid(_obstacles);
 //   SetSideObstacles();
  }
  // compute the CFL condition

  // wipe forces
  _force.clear();

  // wipe boundaries
  _velocity.setZeroBorder();
  _density.setZeroBorder();

  // compute the forces
  // addBuoyancy(density().data());
    // Add buoyancy to the force field.
  if (_buoyancy > 0.0 && _totalSteps < buoyancy_step_) {
    AddBuoyancy();
  }
  _velocity.axpy(_dt, _force);
  _force.clear();


  addVorticity();
  _velocity.axpy(_dt, _force);

  // advect everything
  advectStam();
AdvectParticles();

  // run the solvers
  projectDCT();

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
// add a test cube of density to the center
//////////////////////////////////////////////////////////////////////
void FLUID_3D_DCT::addSmokeColumn() 
{
	addSmokeTestCase(_density.data(), _res);
}

//////////////////////////////////////////////////////////////////////
// generic static version, so that it can be applied to the
// WTURBULENCE grid as well
//////////////////////////////////////////////////////////////////////
void FLUID_3D_DCT::addSmokeTestCase(Real* field, VEC3I res)
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

void FLUID_3D_DCT::AddSmokeTestCase() {
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

void FLUID_3D_DCT::InitializeObstacles(const ObstacleParams3D& param_) {
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
}

void FLUID_3D_DCT::AddDensity(const int xpos, const int ypos, const int zpos, 
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
//////////////////////////////////////////////////////////////////////
// Advect using the semi-Lagrangian method of Stam
//////////////////////////////////////////////////////////////////////
void FLUID_3D_DCT::advectStam()
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
}

//////////////////////////////////////////////////////////////////////
// add vorticity to the force field
//////////////////////////////////////////////////////////////////////
void FLUID_3D_DCT::addVorticity()
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
// Forward transform velocity to the frequency.
void FLUID_3D_DCT::ForwardTransformToFrequency(const VECTOR3_FIELD_3D& v, Real* vxF, Real* vyF, Real* vzF) {
  // Put vx, vy, vz into vxF, vyF, vzF.
  for (int i = 0; i < _totalCells; i++) {
    vxF[i] = v[i][0];
    vyF[i] = v[i][1];
    vzF[i] = v[i][2];
  }
  
  fftw_plan plan_x_3D;
  // x: RODFT10, y: REDFT10 z: REDFT10.
  plan_x_3D = fftw_plan_r2r_3d(_zRes, _yRes, _xRes, vxF, vxF,
                               FFTW_REDFT10 , FFTW_REDFT10, FFTW_RODFT10, FFTW_MEASURE);
  fftw_execute(plan_x_3D);
  fftw_destroy_plan(plan_x_3D);
  
  fftw_plan plan_y_3D;
  // x: REDFT10, y: RODFT10 z: REDFT10.
  plan_y_3D = fftw_plan_r2r_3d(_zRes, _yRes, _xRes, vyF, vyF,
                               FFTW_REDFT10, FFTW_RODFT10, FFTW_REDFT10, FFTW_MEASURE);
  fftw_execute(plan_y_3D);
  fftw_destroy_plan(plan_y_3D);
  
  fftw_plan plan_z_3D;
  // x: REDFT10, y: REDFT10 z: RODFT10.
  plan_z_3D = fftw_plan_r2r_3d(_zRes, _yRes, _xRes, vzF, vzF,
                               FFTW_RODFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_MEASURE);
  fftw_execute(plan_z_3D);
  fftw_destroy_plan(plan_z_3D);
  // normalized.
  for (int i = 0; i < _totalCells; i++) {
    vxF[i] *= invTotalSize_*0.125;
    vyF[i] *= invTotalSize_*0.125;
    vzF[i] *= invTotalSize_*0.125;
  }  
}

void FLUID_3D_DCT::InverseTransformToVelocity(Real* vxF, Real* vyF, Real* vzF, VECTOR3_FIELD_3D* v) {
  /*for (int i = 0; i < _totalCells; i++) {
    vxF[i] *= 0.125;
    vyF[i] *= 0.125;
    vzF[i] *= 0.125;
  }*/
  fftw_plan plan_x_3D;
  plan_x_3D = fftw_plan_r2r_3d(_zRes, _yRes, _xRes, vxF, vxF,
                                FFTW_REDFT01, FFTW_REDFT01, FFTW_RODFT01, FFTW_MEASURE);
  fftw_execute(plan_x_3D);  
  fftw_destroy_plan(plan_x_3D);
  
  fftw_plan plan_y_3D;
  plan_y_3D = fftw_plan_r2r_3d(_zRes, _yRes, _xRes, vyF, vyF,
                               FFTW_REDFT01, FFTW_RODFT01, FFTW_REDFT01, FFTW_MEASURE);
  fftw_execute(plan_y_3D);  
  fftw_destroy_plan(plan_y_3D);
  
  fftw_plan plan_z_3D;
  plan_z_3D = fftw_plan_r2r_3d(_zRes, _yRes, _xRes, vzF, vzF,
                               FFTW_RODFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_MEASURE);
  fftw_execute(plan_z_3D);
  fftw_destroy_plan(plan_z_3D);
  
  // Put the value back into field.
  for (uint i = 0; i < _totalCells; i++) {
    (*v)[i][0] = vxF[i];
    (*v)[i][1] = vyF[i]; 
    (*v)[i][2] = vzF[i];
  }
}

void FLUID_3D_DCT::InitializeFFTW() {
  std::cout <<  "InitializeFFTW begin" << std::endl;
  for (uint i = 0; i < _totalCells; i++) {
    _velocity[i][0] = std::rand() / static_cast<double>(RAND_MAX);
    _velocity[i][1] = std::rand() / static_cast<double>(RAND_MAX);
    _velocity[i][2] = std::rand() / static_cast<double>(RAND_MAX);
  }
  ForwardTransformToFrequency(_velocity, vxFreq_, vyFreq_, vzFreq_);
  InverseTransformToVelocity(vxFreq_, vyFreq_, vzFreq_, &_velocity);
  std::cout <<  "InitializeFFTW end." << std::endl;
  _velocity.clear();
  memset(vxFreq_, 0x00, sizeof(Real)*_xRes*_yRes*_zRes);
  memset(vyFreq_, 0x00, sizeof(Real)*_xRes*_yRes*_zRes);
  memset(vzFreq_, 0x00, sizeof(Real)*_xRes*_yRes*_zRes);
}

// Use DCT to enforce the divergence free condition. Note only Dirichlet BC
// are supported.
void FLUID_3D_DCT::projectDCT() {
  Timer timer;
  timer.Reset();
  
  ForwardTransformToFrequency(_velocity, vxFreq_, vyFreq_, vzFreq_);
  totalDCTtime_ += timer.ElapsedTimeInSeconds();
  
  // Project the component that along the wavenumber away.
  for (int kz = 0; kz < _zRes; kz++) {
    for (int ky = 0; ky < _yRes; ky++) {
      for (int kx = 0; kx < _xRes; kx++) {
        const int index = kx + ky*_xRes + kz*_slabSize;
        
        // Wavenumber.
        Eigen::Vector3d K(static_cast<double>(kx)*invxRes_,
                          static_cast<double>(ky)*invyRes_,
                          static_cast<double>(kz)*invzRes_);
        Eigen::Vector3d W; // (vxFreq_[index], vyFreq_[index], vzFreq_[index]);
        if (kx == 0) W[0] = 0; else W[0] = vxFreq_[index - 1];
        if (ky == 0) W[1] = 0; else W[1] = vyFreq_[index - _xRes];
        if (kz == 0) W[2] = 0; else W[2] = vzFreq_[index - _slabSize];
        
        if (kz != 0 || ky!= 0 || kx != 0) {
          W = W - 1.0 / K.squaredNorm() * (K.dot(W))*K;
        }
        
        if (kx != 0) vxFreq_[index - 1] = W[0];
        if (ky != 0) vyFreq_[index - _xRes] = W[1];
        if (kz != 0) vzFreq_[index - _slabSize] = W[2];
      }
    }
  }
  // Inverse transform to velocity.
  timer.Reset();
  InverseTransformToVelocity(vxFreq_, vyFreq_, vzFreq_, &_velocity);
  totalDCTtime_ += timer.ElapsedTimeInSeconds();
}

void FLUID_3D_DCT::AttenuateSmoke() {
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

void FLUID_3D_DCT::AddBuoyancy() {
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

void FLUID_3D_DCT::AdvectParticles() {
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

void FLUID_3D_DCT::ReSeedParticles() {
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

void FLUID_3D_DCT::DrawSmoke() {

  if (!use_two_phase_smoke_) {
    _density.draw(1.0);
  } else {
    density_sum_ = _density;
    for (int i = 0; i < denWarppers_.size(); i++) {
      density_sum_ += denWarppers_[i].density_;
    }
    density_sum_.draw(1.0);
  }
}

void FLUID_3D_DCT::DrawObstacles() {
  if (use_obstacles_) {
    obstacles_wrapper_.get()->Draw();
  }
}

void FLUID_3D_DCT::OutputSmokeToFolder() {
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

void FLUID_3D_DCT::DrawParticles(const double ptl_length) {
  drawer_.get()->DrawParticles(particles_, _dx, _dt, ptl_length);
}

#define EXPECT_STRING(str) \
in >> temp; \
if (temp != str) { \
  std::cout << "fluid_3D_DCT.cpp " << __LINE__ << " FATAL: " <<  "Error: " << temp << std::endl; exit(0); \
}\
temp.clear();
void FLUID_3D_DCT::ParseSourceFromFile(const std::string& fname) {
  std::ifstream in(fname);
  if (!in.is_open()) {
     std::cout << "fluid_3D_DCT.cpp " << __LINE__ << " FATAL: " <<  "Cannot open: " << fname << std::endl; exit(0);
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

void FLUID_3D_DCT::InitializeSourceSmoke(const VEC3F& pos, const VEC3F& ssize,
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
