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

#include "Eigen"
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sys/stat.h>
 
#include "laplacian_fluid_3d.h"
#include "Alg/MATRIX3.h"
#include "3D/drawer_3d.h"
#include "3D/laplacian_basis_set_3d.h"
#include "3D/dirichlet_basis_set_3d.h"
#include "3D/one_neumann_basis_set_3d.h"
#include "3D/two_neumann_x_3d_basis_set.h"
#include "3D/four_neumann_basis_set_3d.h"
#include "3D/six_neumann_basis_set_3d.h"
#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/FIELD_3D.h"
#include "3D/particle_3d.h"
#include "3D/obstacle_wrapper_3d.h"
#include "3D/density_warpper.h"

#include "solver/integrator_2d.h"
#include "solver/trapezoidal.h"
#include "solver/obstacle_solver.h"
#include "util/util.h"
#include "util/read_write_tensor.h"
#include "util/stringprintf.h"
#include "util/write_density_pbrt.h"

void LaplacianFluid3D::Initialize() {
  // Initialize all members.
  // Basis.
  if (basis_type_ == "all_dirichlet") {
    basis_set_.reset(new DirichletBasisSet3D(des_basis_dim_, xRes_, yRes_, zRes_, 
                     constant_init_strategy_));
    basis_typeint_ = 0;
  } else if (basis_type_ == "one_neumann") {
    basis_set_.reset(new OneNeumannBasisSet3D(des_basis_dim_, xRes_, yRes_, zRes_, 
                     constant_init_strategy_));
    basis_typeint_ = 1;
  } else if (basis_type_ == "two_neumann_x") {
    basis_set_.reset(new TwoNeumannXBasisSet3D(des_basis_dim_, xRes_, yRes_, zRes_, 
                     constant_init_strategy_));
    basis_typeint_ = 2;
  } else if (basis_type_ == "four_neumann_xz") {
    basis_set_.reset(new FourNeumannBasisSet3D(des_basis_dim_, xRes_, yRes_, zRes_, 
                     constant_init_strategy_));
    basis_typeint_ = 4;
  } else if (basis_type_ == "six_neumann") {
    basis_set_.reset(new SixNeumannBasisSet3D(des_basis_dim_, xRes_, yRes_, zRes_, 
                     constant_init_strategy_));
    basis_typeint_ = 6;
  } else {
    std::cout << "laplacian_fluid_3d.cpp " << __LINE__ << " FATAL: " <<  "Unknown basis type: " << basis_type_ << std::endl; exit(0);
  }
  bool from_file_success = false;
  if (basis_tensor_file_.size() > 0) {
    from_file_success = InitializeBasisTensorFromFile();
  }
  if (!from_file_success) {
    // Allocate basis.
    basis_dim_ = basis_set_.get()->AllocateAllBasis();
    // Fill the advection tensor.
    timer_.Reset();
    basis_set_.get()->FillVariationalTensor(&Adv_tensor_);
    std::cout <<  "Time to compute the tensor: " << timer_.ElapsedTimeInSeconds() << std::endl;
  }
  
  if (basis_dim_ < 0) {std::cout << __LINE__ << " FATAL: "; exit(0);}
  velocity_ = VECTOR3_FIELD_3D(xRes_, yRes_, zRes_);
  force_ = VECTOR3_FIELD_3D(xRes_, yRes_, zRes_);
  normal_velocity_ = VECTOR3_FIELD_3D(xRes_, yRes_, zRes_);
  density_ = FIELD_3D(xRes_, yRes_, zRes_);
  density_old_ = FIELD_3D(xRes_, yRes_, zRes_);
  if (use_two_phase_smoke_) {
    density_sum_ = FIELD_3D(xRes_, yRes_, zRes_);
  }
  
  if (use_MacCormack_) {
    temp1 = FIELD_3D(xRes_, yRes_, zRes_);
    temp2 = FIELD_3D(xRes_, yRes_, zRes_);
  }
  coefficients_.resize(basis_dim_);
  coefficients_.setZero();
  coefficients_old_.resize(basis_dim_);
  coefficients_old_.setZero();
  basisWeights_.resize(basis_dim_);
  basisWeights_.setOnes();
  if (coefficient_file_in_.size() == 0) {
    basis_set_.get()->ComputeBasisWeights(&basisWeights_, weightMultiplierC_);
    basis_set_.get()->MultiplyweightTensor(Adv_tensor_, basisWeights_);
  }
  eigenValues_.resize(basis_dim_);
  eigenValues_.setZero();
  basis_set_.get()->ComputeEigenValues(&eigenValues_);
  // integrator_.
  if (solver_type_ == "trapezoidal") {
    integrator_.reset(new Trapezoidal());
  } else {
    std::cout << "laplacian_fluid_3d.cpp " << __LINE__ << " FATAL: " <<  "Unknown solver_type_: " << solver_type_ << std::endl; exit(0);
  }
  // Pointers.
  drawer_.reset(new Drawer3D());
  // Particles.
  for (int i = 0; i < num_particles_; i++) {
    Particle3D particle;
    particles_.push_back(particle);
  }
  ReSeedParticles();
  // ATTENTION: We use FFTW_MEASURE flag for transformation.
  // FFTW will overwite the actuall data while planning, need to
  // do one transformation first.
  InitializeFFTW();

}

bool LaplacianFluid3D::InitializeBasisTensorFromFile() {
  std::ifstream infile(basis_tensor_file_.c_str());
  if (!infile.is_open()) {
    std::cout << "laplacian_fluid_3d.cpp " << __LINE__ << " FATAL: " <<  "Cannot open tensor file: " << basis_tensor_file_ << std::endl; exit(0);
    return false;
  }
  basis_dim_ = basis_set_.get()->ReadBasis(infile);
  std::cout <<  "tensor dimension: " << basis_dim_ << std::endl;
  // Read tensor.
  if (coefficient_file_in_.size() == 0) {
    ReadTensor(infile, basis_dim_, basis_typeint_, &Adv_tensor_);
    std::cout <<  "read tensor completed." << std::endl;
  }
  return true;
}

void LaplacianFluid3D::InitializeObstacles(const ObstacleParams3D& param_) {
  std::cout <<  "Initialize the obstacle..." << std::endl;
  if (!param_.use_obstacles) {
    use_obstacles_ = false;
    return;
  } else  {
    // obstacle_.reset(new Obstacle2D(xRes_, yRes_, param_.obstacle_force_scale));
    // obstacle_.get()->InitializeUsingImage(param_.img_file_name);
    use_obstacles_ = true;
    handle_obstacle_implicit_ = param_.handle_obstacle_implicit;
    force_amp_ = param_.obstacle_force_scale;
    move_obstacle_ = param_.move_obstacle;
  }
  obstacles_wrapper_.reset(new ObstacleWrapper3D(dx_, param_.obstacle_type,
                param_.obstacle_file, param_.obstacle_list, VEC3I(xRes_, yRes_, zRes_)));
  if (handle_obstacle_implicit_) {
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
  }
}

void LaplacianFluid3D::Step() {

  //MATRIX3::rotation(VEC3F(0,1,0), 0.2*dt_*static_cast<float>(frame_simulated_));
  if (use_obstacles_) {
      if (obstacles_wrapper_.get()->GetTrigger()) {
      coefficients_.setZero(); coefficients_old_.setZero();
       density_.clear();
       density_old_.clear();
       AddSmokeTestCaseSphere();
    }
  }
  
  float DCT_time = 0;
  float obstacle_time = 0;
  float solver_time = 0;
  float contraction_time = 0;
  float advection_time = 0;

  current_energy_ = coefficients_.squaredNorm();
  
  // Swap the coefficients.
  for (int i = 0; i < coefficients_.size(); i++) {
    coefficients_old_[i] = coefficients_[i];
  }
  
  double condition_number = 1.0;
  if (coefficient_file_in_.size() == 0) {
    integrator_.get()->IntegrateForward(Adv_tensor_, coefficients_old_,
                                        dt_,
                                        &coefficients_, &condition_number,
                                        &contraction_time, &solver_time);
  } else {
    readBasisCoefficients();
  }
  if (coefficient_file_.size() != 0) {
    writeEigenDense_binary(*coefficient_file_out_.get(), coefficients_);
  }
  
  if (condition_number > maximum_condition_number_) {
    maximum_condition_number_ = condition_number;
  }
  
  // Disspate energy from viscosity.
  DissipateEnergy();
  timer_.Reset();
  // Project out of obstacle, explicit method.
  if (use_obstacles_ && !handle_obstacle_implicit_) {
    ProjectOutObstacle();
  }
  obstacle_time += timer_.ElapsedTimeInSeconds(); 
  
  // Add buoyancy to the force field.
  if (buoyancy_ > 0.0 && frame_simulated_ < buoyancy_step_ && !addforceAtSourceSmoke_) {
    AddBuoyancy();
  }
  
  if (addforceAtSourceSmoke_ && frame_simulated_ < buoyancy_step_) {
    AddforceAtSourceSmoke();
  }
  
  // Add external forces. (DCT is used here)
  timer_.Reset();
  AddExternalForce();
  DCT_time += timer_.ElapsedTimeInSeconds();
  
  timer_.Reset();
  // Implicit method for handling obstacles. TODO: Handle moving obstacles.
  if (use_obstacles_ && handle_obstacle_implicit_) {
    std::cout <<  "Implicit method used....." << std::endl;
    // Use coefficients_old_ as a temp buff.
    SolveObstacleImplicit(coefficients_, obs_proj_mat_, &coefficients_old_);
    coefficients_ = coefficients_old_;
  }
  obstacle_time += timer_.ElapsedTimeInSeconds();
  
  // Reconstruct the velocity field.
  timer_.Reset();

  basis_set_.get()->InverseTramsformToVelocity(coefficients_, &velocity_);
  // std::cout <<  "Div: " << velocity_.computeAbsDivergence() << std::endl;

  DCT_time += timer_.ElapsedTimeInSeconds();
  
  timer_.Reset();
  // Zero out things inside the obstacle.
  if (use_obstacles_) {
    ZeroOutInsideObstacle();
  }
  obstacle_time += timer_.ElapsedTimeInSeconds();
 
  timer_.Reset();
  // Advect density.
  AdvectDensity();
  if (attenuate_smoke_) {
    AttenuateSmoke();
  }
  // Move the obstacle.
  if (use_obstacles_ && move_obstacle_) {
    obstacles_wrapper_.get()->MoveObstacle(dt_);
   }
  // Advect particles.
  AdvectParticles();
  advection_time = timer_.ElapsedTimeInSeconds();

  frame_simulated_ ++;
  //std::fstream fs;
  //fs.open ("./enegy.txt", std::fstream::in | std::fstream::out | std::fstream::app);
  //fs << current_energy_ << "\n";
  //fs.close();
  
  std::cout <<  "Total energy: " << current_energy_ << std::endl;
  // std::cout <<  "Solver: " << " time: " << solver_time << " contraction_time: " << contraction_time << std::endl;
  // std::cout <<  "DCT time: " << DCT_time << std::endl;
  // std::cout <<  "obstacle time: " << obstacle_time << std::endl;
  // std::cout <<  "Advection time: " << advection_time << std::endl;
  std::cout <<  "Frame simulated: " << frame_simulated_ << std::endl;
  // std::cout <<  density_.max() << " " << denWarppers_[0].density_.max() << std::endl;
  contraction_time_ += contraction_time;
  solver_time_ += solver_time;
  density_advection_time_ += advection_time;  
  transformation_time_ += DCT_time;
  obstacle_time_ += obstacle_time;
  
  rest_frame_ --;
  if (rest_frame_ <=0) {
    Quit();
  }

}

void LaplacianFluid3D::DissipateEnergy() {
  for (int i = 0; i < basis_dim_; i++) {
    coefficients_[i] *= exp(-1.0 * eigenValues_[i] * dt_ * visc_);
  }
}

void LaplacianFluid3D::AddSmokeTestCase() {
  VEC3F pointTransformed = source_pos_flt_ - cell_center_;
  pointTransformed += cell_center_;
  
  souce_pos_[0] = static_cast<int>(pointTransformed[0]*maxRes_);
  souce_pos_[1] = static_cast<int>(pointTransformed[1]*maxRes_);
  souce_pos_[2] = static_cast<int>(pointTransformed[2]*maxRes_);
  
  AddDensity(souce_pos_[0], souce_pos_[1], souce_pos_[2],
             souce_size_[0], souce_size_[1], souce_size_[2], &density_);
  if (use_two_phase_smoke_) {
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].AddSmokeTestCase();
    }
  }
}

void LaplacianFluid3D::AddforceAtSourceSmoke() {
  VEC3F pointTransformed = source_pos_flt_ - cell_center_;
  pointTransformed += cell_center_;
  
  int xpos = static_cast<int>(pointTransformed[0]*maxRes_);
  int ypos = static_cast<int>(pointTransformed[1]*maxRes_);
  int zpos = static_cast<int>(pointTransformed[2]*maxRes_);
  int length = maxRes_*2; //souce_size_[0]*1.5;
  int width = souce_size_[1];
  int height = souce_size_[2];
  
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
  
  xbegin = (xbegin > xRes_ - 1) ? xRes_ - 1 : xbegin;
  xend = (xend > xRes_ - 1) ? xRes_ - 1 : xend;
  ybegin = (ybegin > yRes_ - 1) ? yRes_ - 1 : ybegin;
  yend = (yend > yRes_ - 1) ? yRes_ - 1 : yend;
  zbegin = (zbegin > zRes_ - 1) ? zRes_ - 1 : zbegin;
  zend = (zend > zRes_ - 1) ? zRes_ - 1 : zend;
  
  uint idx_begin = xbegin + ybegin*xRes_ + zbegin*Slabsize_;
  
  for (int k = zbegin; k <= zend; k++) {
    for (int j = ybegin; j <= yend; j++) {
      for (int i = xbegin; i <= xend; i++) {
        if (use_obstacles_) {
        // Prevent seeding smoke inside obstacles.       
        VEC3F pt((static_cast<float>(i) + 0.5)*dx_,
                (static_cast<float>(j) + 0.5)*dx_,
                (static_cast<float>(k) + 0.5)*dx_);
          if (! obstacles_wrapper_.get()->inside(pt)) {
            int ind = idx_begin + i + j*xRes_ + k*Slabsize_;
            force_(i,j,k)[0] += buoyancy_;
          }
        } else {
          int ind = idx_begin + i + j*xRes_ + k*Slabsize_;
          force_(i,j,k)[0] += buoyancy_;
        }
      }
    }
  }

  Eigen::VectorXd force_coef(basis_dim_); force_coef.setZero();
  basis_set_.get()->ForwardTransformtoFrequency(force_, &force_coef);
  for (int i = 0; i < basis_dim_; i++) {
    coefficients_[i] += force_coef[i] * dt_;
  }
  force_.clear();
}

void LaplacianFluid3D::AddSmokeTestCaseCylinder() {
  VEC3I res(xRes_, yRes_, zRes_);
  const int slabSize = res[0]*res[1]; 
  int maxRes = res.maxElement();
  Real dx = 1.0f / (Real)maxRes;

  Real yTotal = dx * res[1];
  Real zTotal = dx * res[2];

  // Add cylinder of smoke to target position. The center is source_pos_flt_
  // Oriented along x position.
  // Height of the smoke is in source_size_[0]
  
  Real heighMin = source_pos_flt_[0] - 0.5*souce_size_[0]*dx_;
  Real heighMax = source_pos_flt_[0] + 0.5*souce_size_[0]*dx_;
  

  for (int z = 0; z < res[2]; z++)
    for (int y = 0; y < res[1]; y++)
      for ( int x = (int)(heighMin*res[0]); x <= (int)(heighMax * res[0]); x++)
      {
        Real yLength = y * dx - source_pos_flt_[1];
        Real zLength = z * dx - source_pos_flt_[2];
        Real radius = sqrtf(yLength * yLength + zLength * zLength);

        if (radius < souce_size_[1] * yTotal*dx_)
        {
          int index = x + y * res[0] + z * slabSize;
          if (index >= 0 && index < Totalsize_) {
            density_[index] = added_smoke_density_;
          }
        }
      }
}

void LaplacianFluid3D::AddSmokeTestCaseSphere() {
  // Use center of box smoke,
#pragma omp parallel for
  for (int k = 0; k  < zRes_; k++) {
    for (int j = 0; j < yRes_; j++) {
      for (int i = 0; i < xRes_; i++) {
        const uint index = i + j*xRes_ + k*Slabsize_;
        const float xpos = static_cast<float>(i)*dx_;
        const float ypos = static_cast<float>(j)*dx_;
        const float zpos = static_cast<float>(k)*dx_;
        
        const float rad = (xpos - source_pos_flt_[0])*(xpos - source_pos_flt_[0]) +
                          (ypos - source_pos_flt_[1])*(ypos - source_pos_flt_[1]) + 
                          (zpos - source_pos_flt_[2])*(zpos - source_pos_flt_[2]);
        // take the width as rad ...
        if (rad < souce_size_[0]*dx_*souce_size_[0]*dx_) {
          density_[index] = added_smoke_density_;
        }
      }
    }
  }
}

void LaplacianFluid3D::ClearDensity() {
  density_.clear();
  density_old_.clear();
  if (use_two_phase_smoke_) {
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].density_.clear();
      denWarppers_[i].density_old_.clear();
    }
  }
}

void LaplacianFluid3D::ResetCoeff() {
  coefficients_.setZero();
  coefficients_old_.setZero();
}

void LaplacianFluid3D::AddDensity(const int xpos, const int ypos, const int zpos, 
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
  
  xbegin = (xbegin > xRes_ - 1) ? xRes_ - 1 : xbegin;
  xend = (xend > xRes_ - 1) ? xRes_ - 1 : xend;
  ybegin = (ybegin > yRes_ - 1) ? yRes_ - 1 : ybegin;
  yend = (yend > yRes_ - 1) ? yRes_ - 1 : yend;
  zbegin = (zbegin > zRes_ - 1) ? zRes_ - 1 : zbegin;
  zend = (zend > zRes_ - 1) ? zRes_ - 1 : zend;
  
  uint idx_begin = xbegin + ybegin*xRes_ + zbegin*Slabsize_;
  
  for (int k = zbegin; k <= zend; k++) {
    for (int j = ybegin; j <= yend; j++) {
      for (int i = xbegin; i <= xend; i++) {
        if (use_obstacles_) {
        // Prevent seeding smoke inside obstacles.       
        VEC3F pt((static_cast<float>(i) + 0.5)*dx_,
                (static_cast<float>(j) + 0.5)*dx_,
                (static_cast<float>(k) + 0.5)*dx_);
          if (! obstacles_wrapper_.get()->inside(pt)) {
            int ind = idx_begin + i + j*xRes_ + k*Slabsize_;
            (*field)(i,j,k) = added_smoke_density_;
          }
        } else {
          int ind = idx_begin + i + j*xRes_ + k*Slabsize_;
          (*field)(i,j,k) = added_smoke_density_;
        }
      }
    }
  }
}

void LaplacianFluid3D:: AdvectDensity() {
  const double dt0 = dt_ / dx_;
  density_.swapPointers(density_old_);
  if (!use_MacCormack_) {
  VECTOR3_FIELD_3D::advect(dt0, velocity_, density_old_, density_);
  } else {
    VECTOR3_FIELD_3D::advectMacCormack(dt0, velocity_, density_old_, density_, temp1, temp2);
  }
  density_.setZeroBorder();
  
  if (use_two_phase_smoke_) {
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].SwapPointer();
    }
    if (!use_MacCormack_) {
      for (int i = 0; i < denWarppers_.size(); i++) {
        VECTOR3_FIELD_3D::advect(dt0, velocity_, denWarppers_[i].density_old_, denWarppers_[i].density_);
      }
    }   else {
        for (int i = 0; i < denWarppers_.size(); i++) {
          std::cout <<  "using maccc" << std::endl;
          VECTOR3_FIELD_3D::advectMacCormack(dt0, velocity_,denWarppers_[i].density_old_,
                                               denWarppers_[i].density_, temp1, temp2);
      }
    }
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].density_.setZeroBorder();
    }
  }
}

void LaplacianFluid3D::AddExternalForce() {
  Eigen::VectorXd force_coef(basis_dim_); force_coef.setZero();
  basis_set_.get()->ForwardTransformtoFrequency(force_, &force_coef);
  for (int i = 0; i < basis_dim_; i++) {
    coefficients_[i] += force_coef[i] * dt_;
  }
  force_.clear();
}

void LaplacianFluid3D::AddBuoyancy() {

  VEC3F dir_ = buoyancy_dir_;
  for (uint i = 0; i < Totalsize_; i++) {
    force_[i][0] += buoyancy_*density_[i]*dir_[0];
    force_[i][1] += buoyancy_*density_[i]*dir_[1];
    force_[i][2] += buoyancy_*density_[i]*dir_[2];
  }
  
  if (use_two_phase_smoke_) {
    /*VEC3F dir_ = buoyancy_dir_;
    for (uint i = 0; i < Totalsize_; i++) {
      force_[i][0] += buoyancy_*density1_[i]*dir_[0];
      force_[i][1] += buoyancy_*density1_[i]*dir_[1];
      force_[i][2] += buoyancy_*density1_[i]*dir_[2];
    }*/
    // point upwards.
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].AddBuoyancyField(&force_);
    }
  }
}

void LaplacianFluid3D::ProjectOutObstacle() {
  obstacles_wrapper_.get()->CalculateNormalForce(velocity_, dt_, dx_, force_amp_, &force_);
}

void LaplacianFluid3D::ZeroOutInsideObstacle() {
  obstacles_wrapper_.get()->ZeroOutVelocityDensityObstalce(dx_, &velocity_, &density_);
  if (use_two_phase_smoke_) {
    for (int i = 0; i < denWarppers_.size(); i++) {
      obstacles_wrapper_.get()->ZeroOutVelocityDensityObstalce(dx_, &velocity_, &denWarppers_[i].density_);
    }
  }
}

#define EXPECT_STRING(str) \
in >> temp; \
if (temp != str) { \
  std::cout << "laplacian_fluid_3d.cpp " << __LINE__ << " FATAL: " <<  "Error: " << temp << std::endl; exit(0); \
}\
temp.clear();
void LaplacianFluid3D::ParseSourceFromFile(const std::string& fname) {
  std::ifstream in(fname);
  if (!in.is_open()) {
    std::cout << "laplacian_fluid_3d.cpp " << __LINE__ << " FATAL: " <<  "Cannot open: " << fname << std::endl; exit(0);
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
    
    DensityWarpper denWarpper(VEC3I(xRes_, yRes_, zRes_),
                                          buoyancy_dir, buoyancy_, added_smoke_density_);
    denWarpper.InitialzeSourceSomke(position, ssize, maxRes_);
    denWarppers_.push_back(denWarpper);
  } while(!in.eof());
}
#undef EXPECT_STRING

void LaplacianFluid3D::InitializeSourceSmoke(const VEC3F& pos, const VEC3F& ssize,
                                             const std::string& source_file, bool use_cylinder, bool addForceSoureSmoke) {
  source_pos_flt_[0] = pos[0];
  source_pos_flt_[1] = pos[1];
  source_pos_flt_[2] = pos[2];
  souce_size_[0] = static_cast<int>(ssize[0]*maxRes_);
  souce_size_[1] = static_cast<int>(ssize[1]*maxRes_);
  souce_size_[2] = static_cast<int>(ssize[2]*maxRes_);
  addforceAtSourceSmoke_ = addForceSoureSmoke;
  
  if (source_file.size() > 0) {
    ParseSourceFromFile(source_file);
  }
  if (use_two_phase_smoke_) {
    if (!use_cylinder) { 
      AddSmokeTestCase();
    } else {
      AddSmokeTestCaseCylinder();
    }
  }
}

void LaplacianFluid3D::Quit() {
  std::cout <<  "Total frame simulated: " << total_frame_ << std::endl;
  std::cout << "Average time for contraction per frame: " 
      << contraction_time_ / total_frame_;
  std::cout << "Average time for solver per frame: " 
      << solver_time_ / total_frame_;
  std::cout << "Average time for DCT and DST: "
      << transformation_time_ / total_frame_;
  std::cout << "Average time for advection of density and particles: "
      << density_advection_time_ / total_frame_;
  std::cout << "Average time for handling obstacles: "
      << obstacle_time_ / total_frame_;
      
  std::cout <<  "Maximum condition number: " << maximum_condition_number_ << std::endl;
  finished_ = true;
  if (coefficient_file_.size() != 0) {
    coefficient_file_out_.get()->close();
  }
}

void LaplacianFluid3D::OutputSmokeToFolder() {
  std::string fname = StringPrintf("%s%04d.pbrt", density_folder_.c_str(), frame_simulated_);
  std::cout <<  "Output the smoke to file: " << fname << std::endl;
  if (!use_two_phase_smoke_) {
    WriteDensityPBRT(fname, density_, 1.0);
  } else {
    density_sum_ = density_; // + density1_;
    for (int i = 0; i < denWarppers_.size(); i++) {
      density_sum_ += denWarppers_[i].density_;
    }
    WriteDensityPBRT(fname, density_sum_, 1.0);
  }
  if (use_obstacles_) {
    obstacles_wrapper_.get()->WriteObstacleToPBRT(density_folder_.c_str(), frame_simulated_);
  }
}

void LaplacianFluid3D::AdvectParticles() {
  const double dt0 = dt_ / dx_;
  for (int i = 0; i < num_particles_; i++) {
    // Get the velocity at the particles.
    const VEC3& position = particles_[i].position;
    // This function assumen the velocity is stored on the vertices of the grid, while for laplacian fluid
    // case, the velocity is stored on the center of the grid.
    VEC3 p_v = velocity_.GetVelocity(position[0]-0.5, position[1]-0.5, position[2]-0.5);
    // Forward Eular.
    particles_[i].position[0] += p_v[0] * dt0;
    particles_[i].position[1] += p_v[1] * dt0;
    particles_[i].position[2] += p_v[2] * dt0;
    particles_[i].velocity = p_v;
    if (particles_[i].position[0] < 0. || particles_[i].position[0] > xRes_ || 
        particles_[i].position[1] < 0. || particles_[i].position[1] > yRes_ ||
        particles_[i].position[2] < 0. || particles_[i].position[2] > zRes_) {
      particles_[i].position[0] = std::rand() % xRes_ + 0.5;
      particles_[i].position[1] = std::rand() % yRes_ + 0.5;
      particles_[i].position[2] = std::rand() % zRes_ + 0.5;
      particles_[i].velocity = 0.;
    }
  }
}

void LaplacianFluid3D::ReSeedParticles() {
  // Reseed the particles at random position.
  for (int i = 0; i < num_particles_; i++) {
    int px = std::rand() % xRes_;
    int py = std::rand() % yRes_;
    int pz = std::rand() % zRes_;
    particles_[i].position[0] = px + 0.5;
    particles_[i].position[1] = py + 0.5;
    particles_[i].position[2] = pz + 0.5;
    particles_[i].velocity = 0;
  }
}

void LaplacianFluid3D::DrawObstacles() {
  if (use_obstacles_) {
    obstacles_wrapper_.get()->Draw();
  }
}

void LaplacianFluid3D::DrawParticles(const double ptl_length){
  drawer_.get()->DrawParticles(particles_, dx_, dt_, ptl_length);
}

void LaplacianFluid3D::AttenuateSmoke() {
  for (uint i = 0; i < Totalsize_; i++) {
    density_[i] *= density_attenuate_factor_;
    if (density_[i] < 0.003) {
      density_[i] = 0;
    }
  }
  if (use_two_phase_smoke_) {
    for (int i = 0; i < denWarppers_.size(); i++) {
      denWarppers_[i].AttenuateSmoke(density_attenuate_factor_);
    }
  }
}

void LaplacianFluid3D::DrawSmoke(const float alpha_multipler, const int low_cut) {
  timer_.Reset();
  VEC3F rescaleRatio(static_cast<float>(xRes_) / maxRes_, static_cast<float>(yRes_) / maxRes_,
                     static_cast<float>(zRes_) / maxRes_);
  if (!render_initialized_) {

glEnable(GL_TEXTURE_3D);
  glDisable(GL_DEPTH_TEST);
  glCullFace(GL_FRONT);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  } render_initialized_;
  
  if (!use_two_phase_smoke_) {

    glPushMatrix();
    glScalef(rescaleRatio[0], rescaleRatio[1], rescaleRatio[2]);
    density_.draw(alpha_multipler);
    glPopMatrix();
  } else {
    density_sum_ = density_;
    for (int i = 0; i < denWarppers_.size(); i++) {
      density_sum_ += denWarppers_[i].density_;
    }
    
    glPushMatrix();
    glScalef(rescaleRatio[0], rescaleRatio[1], rescaleRatio[2]);
    density_sum_.draw(alpha_multipler);
    glPopMatrix();
  }
}

void LaplacianFluid3D::DrawCoefficients(const double multi_factor) {

  Eigen::VectorXd wEnergy = coefficients_.cwiseProduct(coefficients_);
  drawer_.get()->DrawCoefficients(wEnergy, dx_, multi_factor*20.0);
  //drawer_.get()->DrawSpectogram(wEnergy, multi_factor*100.0);  
}

// ATTENTION: We use FFTW_MEASURE flag for transformation.
// FFTW will overwite the actuall data while planning, need to
// do one bogus transformation first.
void LaplacianFluid3D::InitializeFFTW() {
  std::cout <<  "Intialize FFTW, because we use FFTW_MEASURE." << std::endl;
  for (int i = 0; i < basis_dim_; i++) {
    coefficients_[i] = std::rand() / static_cast<double>(RAND_MAX);
  }
  basis_set_.get()->InverseTramsformToVelocity(coefficients_, &velocity_);
  basis_set_.get()->ForwardTransformtoFrequency(velocity_, &coefficients_);
  coefficients_.setZero();
  velocity_.clear();
  std::cout <<  "FFTW init finished" << std::endl;
}

void LaplacianFluid3D::PrintDebugInfo() {
  basis_set_.get()->PrintDebugInfo(coefficients_);
  exit(0);
}

inline bool fileExists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

void LaplacianFluid3D::readBasisCoefficients() {
  // coefficient_fptr_in_.get()->read(reinterpret_cast<char*>(coefficients_.data()), sizeof(double)*coefficients_.size());
  readEigenDense_binary(*coefficient_fptr_in_.get(), coefficients_);
}

// pass in the delta of mouse on screen, screen size is normalized to [0, 1].
void LaplacianFluid3D::SetOmega(const float mouse_dx, const float mouse_dy,
                const float pos_x, const float pos_y) {
  obstacles_wrapper_.get()->SetOmega(mouse_dx, mouse_dy, pos_x, pos_y, dt_);
}

void LaplacianFluid3D::SetObstacleVelo(const VEC3F& velo) {
  obstacles_wrapper_.get()->SetObstacleVelo(velo);
}
