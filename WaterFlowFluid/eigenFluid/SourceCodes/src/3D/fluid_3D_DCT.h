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

#ifndef FLUID_3D_DCT_H
#define FLUID_3D_DCT_H

#include "Eigen"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>

#include "3D/obstacle_wrapper_3d.h"
#include "3D/Ted_obstacle.h"
#include "Alg/VEC3.h"
#include "3D/FIELD_3D.h"
#include "3D/VECTOR3_FIELD_3D.h"
// #include "MERSENNETWISTER.h"
// #include "SPARSE_MATRIX.h"

#include "3D/density_warpper.h"
#include "3D/particle_3d.h"
#include "3D/drawer_3d.h"

class FLUID_3D_DCT  
{
public:
   FLUID_3D_DCT(int xRes, int yRes, int zRes, int amplify, const bool use_two_phase_smoke, 
               const double added_smoke_density, const std::string& buoyancy_dir, const Real buoyancy,
               const double density_attenuate_factor_, const bool attenuate_smoke,
               const std::string density_folder, const double dt, const int buoyancy_step,
               unsigned int* boundaries = NULL);
  Real& vorticityEps() { return _vorticityEps; };
    
  virtual ~FLUID_3D_DCT();
  void addSmokeColumn();
  void addSmokeTestCase(Real* field, VEC3I res);

  void step();
  void ReSeedParticles();
  void InitializeFFTW();
  
  void AddDensity(const int xpos, const int ypos, const int zpos, 
                  const int length, const int width, const int height,
                  FIELD_3D* field);
  void AddSmokeTestCase();
  void InitializeSourceSmoke(const VEC3F& pos, const VEC3F& ssize,
                             const std::string& source_file);
  void DrawSmoke();
  void OutputSmokeToFolder();
  void InitializeObstacles(const ObstacleParams3D& param_);
  void DrawObstacles();
  void DrawParticles(const double ptl_length);
  float totalDCTtime_ ;
  
protected:
  // dimensions
  int _xRes, _yRes, _zRes, _maxRes;
  int _xs, _ys, _zs;
  Real invxRes_, invyRes_, invzRes_;
  VEC3I _res;
  int _totalCells;
  Real invTotalSize_;
  int _slabSize;
  Real _dx;
  Real _totalTime;
  int _totalSteps;
  int _totalImgDumps;
  int _totalVelDumps;

  // fields
  FIELD_3D _density;
  FIELD_3D _densityOld;
  FIELD_3D _pressure;
  FIELD_3D _divergence;
	
  VECTOR3_FIELD_3D _velocity;
  VECTOR3_FIELD_3D _velocityOld;
  
  Real* vxFreq_;
  Real* vyFreq_;
  Real* vzFreq_;
  void ForwardTransformToFrequency(const VECTOR3_FIELD_3D& v, Real* vxF, Real* vyF, Real* vzF);
  void InverseTransformToVelocity(Real* vxF, Real* vyF, Real* vzF, VECTOR3_FIELD_3D* v);
  
  VECTOR3_FIELD_3D _force;
  VECTOR3_FIELD_3D _vorticity;
  unsigned char*  _obstacles;
  bool* _dirichletObstacleField;

  // store the pre-projection velocity in case the reduced solver needs them
  VECTOR3_FIELD_3D _preprojection;
  
  // store the pre-diffusion velocity in case the reduced solver needs them
  VECTOR3_FIELD_3D _prediffusion;
  
  // store the pre-advection velocity in case the reduced solver needs them
  VECTOR3_FIELD_3D _preadvection;
  
  // cache the result of the Treuille-style quadratic in case the reduced solver needs it
  VECTOR3_FIELD_3D _cachedQuadratic;
  
  // store the post-advection velocity in case the reduced solver needs them
  VECTOR3_FIELD_3D _postadvection;

  // cache the velocity field from beofre vorticity confinement has been added
  VECTOR3_FIELD_3D _prevorticity;

  // cache the intermediate Maccormack field
  VECTOR3_FIELD_3D _midMaccormack;
  
  // cache the result of Neumann IOP
  VECTOR3_FIELD_3D _postIOP;

  FIELD_3D _cachedPressure;
  FIELD_3D _cachedDivergence;

  vector<FIELD_3D> _cachedPressuresIOP;
  vector<FIELD_3D> _cachedDivergencesIOP;
  vector<VECTOR3_FIELD_3D> _cachedVelocitiesIOP;

  VEC3F _center;
  VEC3F _lengths;

  int _iterations;

  // simulation constants
  const Real _dt;
  Real _dtOld;
  const Real _buoyancy;
  const int buoyancy_step_;
  
  
  Real _vorticityEps;
  Real _heatDiffusion;
  Real _solverEps;

  // array of Dirichlet obstacles
  vector<OBSTACLE*> _dirichletObstacles;
  
  // array of Neumann obstacles
  vector<OBSTACLE*> _neumannObstacles;

  VECTOR3_FIELD_3D _neumannVelocityField;
  bool* _neumannObstacleField;

  // file paths
  string _snapshotPath;

  // domain boundary conditions
  unsigned int _domainBcFront;
  unsigned int _domainBcBack;
  unsigned int _domainBcBottom;
  unsigned int _domainBcTop;
  unsigned int _domainBcLeft;
  unsigned int _domainBcRight;

  // procedurally generated velocity basis
  VECTOR3_FIELD_3D _proceduralVelocity;

  vector<MATRIX3> _glyphFrames;

  vector<VEC3F> _glyphPoints;
  vector<VEC3F> _glyphNormals;

  // timestepping functions
  void addVorticity();
  void addBuoyancy(Real *field);

  //preconditioner stuff
  void zeroBoundary(FIELD_3D& x, unsigned char* skip);
	
  void buildMICPreconditioner();
  void applyMICPreconditioner(FIELD_3D& zf, FIELD_3D& r, FIELD_3D& q);	
	
  void buildA(unsigned char* skip, bool heat);	
  void multA(const FIELD_3D& xf, FIELD_3D& bf);
  // solver stuff
  void projectDCT();
  
  // handle obstacle boundaries
  void setObstacleBoundaries();

  // build the damping matrix
  void buildPeeledDampingMatrixFull();

  // advect using first order semi-Lagrangian
  void advectStam();
  
  const bool use_two_phase_smoke_;
  VEC3I souce_pos_;
  VEC3I souce_size_;
  VEC3F cell_center_;
  VEC3F source_pos_flt_;
  const double added_smoke_density_;
  VEC3F buoyancy_dir_;
  void ParseSourceFromFile(const std::string& fname);
  void AddBuoyancy();

  std::vector<DensityWarpper> denWarppers_;
  // Temporary field for MacCormack.
  FIELD_3D temp1;
  FIELD_3D temp2;
  FIELD_3D density_sum_;
  void AttenuateSmoke();
  const double density_attenuate_factor_;
  const bool attenuate_smoke_;
  const std::string density_folder_;
  
  // Obstacles.
  std::unique_ptr<ObstacleWrapper3D> obstacles_wrapper_;
  bool use_obstacles_;
  bool handle_obstacle_implicit_;
  bool move_obstacle_;
  
  std::vector<Particle3D> particles_;
  const int num_particles_ = 0;
  void AdvectParticles();
  std::unique_ptr<Drawer3D> drawer_;
  // Render.
  bool render_initialized_;
};

#endif  // FLUID_3D_DCT_H
