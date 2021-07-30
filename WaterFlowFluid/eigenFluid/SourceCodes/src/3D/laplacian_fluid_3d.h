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

#ifndef LAPLACIAN_FLUID_3D_H
#define LAPLACIAN_FLUID_3D_H

#include "Eigen"
#include <memory>
#include <fstream>
#include <vector>
// #include <glog/logging.h>

#include "solver/integrator_2d.h"
#include "3D/drawer_3d.h"
#include "3D/FIELD_3D.h"
#include "3D/laplacian_basis_set_3d.h"
#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/particle_3d.h"
#include "setting.h"
#include "util/timer.h"
#include "3D/obstacle_wrapper_3d.h"
#include "Alg/MATRIX3.h"
#include "3D/density_warpper.h"

class LaplacianFluid3D {
public:
  LaplacianFluid3D(const int xRes, const int yRes, const int zRes, const int des_basis_dim, 
    const double dt, const double visc, double buoyancy, const std::string& buoyancy_dir,
    const double added_smoke_density,
    const int num_particles, const std::string& basis_type, const std::string& const_strategy,
    const Eigen::Vector3f lightpos,
    const std::string& basis_tensor_file, const int total_frame, const std::string& dens_folder,
    const std::string& solver_type, const bool attenuate_smoke, const double density_attenuate_factor,
    const bool use_MacCormack, const std::string& coefficient_file, const bool use_two_phase_smoke,
    const int buoyancy_step, const std::string& coefficient_file_in, const Real weightMultiplierC)
  :xRes_(xRes), yRes_(yRes), zRes_(zRes), des_basis_dim_(des_basis_dim), dt_(dt),visc_(visc),
  buoyancy_(buoyancy), added_smoke_density_(added_smoke_density),
  num_particles_(num_particles), basis_type_(basis_type), lightpos_(lightpos), 
  basis_tensor_file_(basis_tensor_file), total_frame_(total_frame),
  density_folder_(dens_folder), constant_init_strategy_(const_strategy),
  solver_type_(solver_type), attenuate_smoke_(attenuate_smoke),
  density_attenuate_factor_(density_attenuate_factor), use_MacCormack_(use_MacCormack),
  coefficient_file_(coefficient_file), use_two_phase_smoke_(use_two_phase_smoke), 
  buoyancy_step_(buoyancy_step), coefficient_file_in_(coefficient_file_in), weightMultiplierC_(weightMultiplierC){
    // cell length. The maximum edge length is 1.
     maxRes_ = std::max(std::max(xRes_, yRes_), zRes_);
     dx_ = 1.0 / static_cast<double>(maxRes_);
     cell_center_[0] = static_cast<float>(xRes)*dx_*0.5; cell_center_[1] = static_cast<float>(yRes)*dx_*0.5;
     cell_center_[2] = static_cast<float>(zRes)*dx_*0.5;
     
     Slabsize_ = xRes_*yRes_;
     Totalsize_ = Slabsize_*zRes_;
     basis_dim_ = 0;
     frame_simulated_ = 0;
     rest_frame_ = total_frame_;
     use_obstacles_ = false;
     move_obstacle_ = false;
     force_amp_ = 0;
     finished_ = false;
     render_initialized_ = false;
     buoyancy_dir_ = 0;
    
     if (buoyancy_dir == "x") {
       buoyancy_dir_ = VEC3F(1, 0, 0);
     } else if (buoyancy_dir == "y") {
       buoyancy_dir_ = VEC3F(0, 1, 0);
     } else if (buoyancy_dir == "z") {
       buoyancy_dir_ = VEC3F(0, 0, 1);
     } else {
       std::cout << "laplacian_fluid_3d.h " << __LINE__ << " FATAL: " <<  "Invalid buoyancy_dir: " << buoyancy_dir << std::endl; exit(0);
     }
     
     if (coefficient_file_in_.size() != 0) {
       coefficient_file_ = "";   // read coefficient in, disable write out of coefficient.
       coefficient_fptr_in_.reset(new std::ifstream(coefficient_file_in_));
       if (!coefficient_fptr_in_.get()->is_open()) {
         std::cout << "laplacian_fluid_3d.h " << __LINE__ << " FATAL: " <<  "Cannot open file: " << coefficient_file_in_ << std::endl; exit(0);
       }
     }
     
     if (coefficient_file_.size() != 0) {
       std::cout <<  "rrr" << std::endl;
       coefficient_file_out_.reset(new std::ofstream(coefficient_file_,ios::out | ios::binary | ios::trunc));
       if (!coefficient_file_out_.get()->is_open()) {
         std::cout << "laplacian_fluid_3d.h " << __LINE__ << " FATAL: " <<  "Cannot open file: " << coefficient_file_ << std::endl; exit(0);
       }
     }
     Initialize();
  };
  ~LaplacianFluid3D(){};
   void Step();
   void ReSeedParticles();
   void DrawParticles(const double ptl_length);
   void DrawCoefficients(const double multi_factor);
   void DrawSmoke(const float alpha_multipler, const int low_cut);
   void DrawObstacles();
   void AddSmokeTestCase();
   void AddSmokeTestCaseCylinder();
   void AddSmokeTestCaseSphere();
   
   void ClearDensity();
   void ResetCoeff();
   int GetBasiDim() const {return basis_dim_;}
   
   void InitializeObstacles(const ObstacleParams3D& param_);
   void InitializeSourceSmoke(const VEC3F& pos, const VEC3F& ssize,
                              const std::string& source_file, bool use_cylinder, bool addForceSoureSmoke);
   void ParseSourceFromFile(const std::string& fname);
   bool isFinished() const {return finished_;}
   void PrintDebugInfo();
   void OutputSmokeToFolder();
   // pass in the delta of mouse on screen, screen size is normalized to [-1, 1].
   void SetOmega(const float mouse_dx, const float mouse_dy,
                const float pos_x, const float pos_y);
   void SetObstacleVelo(const VEC3F& velo);
   // Add external force at the position of source smoke.
   void AddforceAtSourceSmoke();
   void readBasisCoefficients();
protected:
  void Initialize();
  // ATTENTION: We use FFTW_MEASURE flag for transformation.
  // FFTW will overwite the actuall data while planning, need to
  // do one bogus transformation first.
  void InitializeFFTW();
  // The resolution of the velocity field.
  const int xRes_;
  const int yRes_;
  const int zRes_;
  int Slabsize_;
  int Totalsize_;
  
  // The dimension (length) of the grid cell.
  double dx_;
  int maxRes_;
  
  // Viscosity.
  const double visc_;
  // Buoyancy.
  const double buoyancy_;
  VEC3F buoyancy_dir_;
  
  const double added_smoke_density_;
  // Timestep
  const double dt_;
  
  // The energy of current timestep.
  double current_energy_;
  // The dimention of coefficient of the basis field.
  const int des_basis_dim_;   // Desired number of basis.
  int basis_dim_;   // Actually allocated basis.
  
  // Number of wave number along each axis.
  // basis_dim_x_*basis_dim_y_*basis_dim_z_ <= basis_dim_.
  int basis_dim_x_;
  int basis_dim_y_;
  int basis_dim_z_;
  
  const std::string basis_type_; 

  Eigen::VectorXd coefficients_;
  Eigen::VectorXd coefficients_old_;
  Eigen::VectorXd basisWeights_;

  std::string coefficient_file_;
  std::unique_ptr<std::ofstream> coefficient_file_out_;
  const std::string coefficient_file_in_;
  std::unique_ptr<std::ifstream> coefficient_fptr_in_;
  
  std::vector<Adv_Tensor_Type> Adv_tensor_;
  std::unique_ptr<LaplacianBasisSet3D> basis_set_;
  const std::string constant_init_strategy_;
  
  std::string basis_tensor_file_;
  int basis_typeint_;
  Eigen::VectorXd eigenValues_;
  
  // Solver.
  std::unique_ptr<Integrator2D> integrator_;
  const std::string solver_type_;
  
  // Grids.
  // Uniform grid to store the velocity. ATTENTION: velocity is located in the center of the grid,
  // rather than on the vertices of the grid.
  VECTOR3_FIELD_3D velocity_;
  // Force field.
  VECTOR3_FIELD_3D force_;
  // Normal velocity.
  VECTOR3_FIELD_3D normal_velocity_;
  // Density field.
  FIELD_3D density_;
  FIELD_3D density_old_;
  FIELD_3D density_sum_;
  std::vector<DensityWarpper> denWarppers_;
  
  const bool use_two_phase_smoke_;
  
  // Temporary field for MacCormack.
  FIELD_3D temp1;
  FIELD_3D temp2;
  const std::string density_folder_;
  
  Timer timer_;
  int frame_simulated_;
  int rest_frame_;
  const int total_frame_;
  double contraction_time_;
  double solver_time_;
  double transformation_time_;
  double density_advection_time_;
  double obstacle_time_;
  
  std::vector<Particle3D> particles_;
  const int num_particles_;
  void AdvectParticles();
  void AdvectDensity();
  const bool use_MacCormack_;
  void DissipateEnergy();
  void AddDensity(const int xpos, const int ypos, const int zpos, 
                  const int length, const int width, const int height,
                  FIELD_3D* field);
  void AttenuateSmoke();
  
  const bool attenuate_smoke_;
  const double density_attenuate_factor_;
  const int buoyancy_step_;
  
  void AddBuoyancy();
  void AddExternalForce();
  bool InitializeBasisTensorFromFile();
  
  std::unique_ptr<Drawer3D> drawer_;
  // statistics.
  double maximum_condition_number_ = 0;
  
  // Render.
  const Eigen::Vector3f lightpos_;
  bool render_initialized_;
  
  // Obstacles.
  std::unique_ptr<ObstacleWrapper3D> obstacles_wrapper_;
  bool use_obstacles_;
  bool handle_obstacle_implicit_;
  bool move_obstacle_;
  
  // Projected object matrix.
  Eigen::MatrixXd obs_proj_mat_;
  // The spring constant for the obstacle force.
  double force_amp_;
  void ProjectOutObstacle();
  // Zero out things inside the obstacle.
  void ZeroOutInsideObstacle();
  VEC3I souce_pos_;
  VEC3I souce_size_;
  VEC3F cell_center_;
  VEC3F source_pos_flt_;
  
  // Weight coefficient in eq 22 of the paper.
  const Real weightMultiplierC_;
  
  // Quit and print some statistics.
  void Quit();
  bool finished_;
  bool addforceAtSourceSmoke_ = false;
  
};

#endif  // LAPLACIAN_FLUID_3D_H
