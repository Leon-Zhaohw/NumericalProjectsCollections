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

#ifndef LAPLACIAN_FLUID_H
#define LAPLACIAN_FLUID_H

#include "Eigen"
#include "Sparse"
// #include <glog/logging.h>
#include <vector>

#include "2D/drawer_2d.h"
#include "solver/integrator_2d.h"
#include "2D/laplacian_basis_2D.h"
#include "2D/FIELD2D.h"
#include "2D/particle_2d.h"
#include "util/timer.h"
#include "setting.h"
#include "2D/VFIELD2D.h"
#include "2D/obstacle_2d.h"

class LaplacianFluid2D{
public:
  LaplacianFluid2D(
    const int xRes, const int yRes, const int basis_dim_root,
    const double dt,
    const std::string& basis_type,
    const std::string& integrator_type,
    double buoyancy, double visc, 
    const int num_particles, const int total_frame, const double added_smoke_density,
    const std::string tensor_fname, const bool addmirrordens): 
    visc_(visc), buoyancy_(buoyancy), dt_(dt),
    xRes_(xRes), yRes_(yRes),
    basis_dim_(basis_dim_root*basis_dim_root),
    basis_dim_root_(basis_dim_root), basis_type_(basis_type),
    integrator_type_(integrator_type),num_particles_(num_particles),
    added_smoke_density_(added_smoke_density),
    total_frame_(total_frame),tensor_fname_(tensor_fname), add_mirror_density_(addmirrordens)
    {
      // TODO: modify the LaplacianBasis2D to support differemt xRes and yRes.
      if (xRes != yRes) {std::cout << "laplacian_fluid_2D.h " << __LINE__ << " FATAL: "  << "Only support square domain for now." << std::endl; exit(0);}
      dx_ = 1.0 / std::max(xRes_, yRes_);
      use_obstacles_ = false;
      Initialize();
    }
  ~LaplacianFluid2D(){}
  // Perform the forward step.
  void Step();
  void DrawDensity();
  void DrawVelocity();
  void DrawVort();
  void DrawParticles(const double ptl_length);
  void DrawCoefficients(const double multi_factor);
  void DrawObstacles();
  void AddSmokeTestCase(const int xpos, const int ypos, const int width,
                        const int height);
  
  void AddForceImpluse(const int x, const int y, const int width,
                       const int height, const double mag);
  
  void AddForce(const int xpos, const int ypos, const float fx, const float fy);
  void AddSmoke(const int xpos, const int ypos, const float amount);
  void MoveObstacles(const double delta_x, const double delta_y);
  void AddVelocityToObstacles(const double dvx, const double dvy);
  void ReSeedParticles();
  void Quit();
  void InitializeObstacles(const ObstacleParams& param_);
  void ComputeEnergyDerivative();
  void ComputeBoundaryIntegral();
  void InitVelocityFromDedalus(const std::string& Ufname, const std::string& Vfname,
                               const int dedalusXRes, const int dedalusYRes);
  
  /*
  // Dual step.
  void DualStep();
  void DualInitialize();*/
  bool addBounyancyOnce = false;
  bool quit_ = false;
private:
  void Initialize();
  void InitCheckBoardPattern();
  // Viscosity.
  const double visc_;
  // Buoyancy.
  const double buoyancy_;
  const double added_smoke_density_;
  // Timestep
  const double dt_;
  // The energy of current timestep.
  double current_energy_;
  // The resolution of the velocity field.
  const int xRes_;
  const int yRes_;
  double dx_;
  void ZeroOutUnWantedCoefficients();
  // The dimention of coefficient of the basis field.
  const int basis_dim_;
  const int basis_dim_root_;
  Eigen::VectorXd basis_coefficients_;
  Eigen::VectorXd basis_coefficients_old_;
  Eigen::VectorXd force_dw_;
  Eigen::VectorXd energy_derivative_;
  Eigen::VectorXd tempVector_;
  Eigen::MatrixXd eigLargeVectors_;
  
  // Basis field eigenvalues
  std::vector<double> eigs_;
  std::vector<double> eigs_inv_;
  std::vector<double> eigs_inv_root_;

  std::unique_ptr<LaplacianBasis2D> basis_;
  std::string basis_type_;
  const std::string tensor_fname_;
  std::vector<Adv_Tensor_Type> Adv_tensor_;
  // Projected object matrix.
  Eigen::MatrixXd obs_proj_mat_;
  
  // Vorticity transfer matrix.
  Eigen::MatrixXd TransferMatrix_;
  std::unique_ptr<Integrator2D> integrator_;
  std::string integrator_type_;
  
  // Class to draw stuff.
  std::unique_ptr<Drawer2D> drawer_;
  // Obstacle.
  std::unique_ptr<Obstacle2D> obstacle_;
  bool use_obstacles_;
  bool handle_obstacle_implicit_;
  bool is_neumann_;
  bool add_mirror_density_;
  
  // Velocity field.
  VFIELD2D velocity_;
  // Force field.
  VFIELD2D force_;
  // Normal velocity.
  VFIELD2D normal_velocity_;
  // Density field.
  FIELD2D density_;
  FIELD2D density_old_;

  FIELD2D temp1;
  FIELD2D temp2;
  FIELD2D vort_;

  // Particles.
  std::vector<Particle2D> particles_;
  Timer timer_;
  
  /*
  const bool use_2nd_order_;
  CTF::Tensor<double>* Q_;
  Eigen::VectorXd coefficients_dot_old_;
  Eigen::VectorXd coefficients_dot_;
  std::unique_ptr<IntegratorImplicit2ndOrder> integrator_2nd_order_;
  CTF::Tensor<double>* C_;
  */
  Eigen::VectorXd basisWeights_;
  const int num_particles_;
  const int total_frame_;
  int rest_frame_;
  int frame_simulated_;
  double integrator_time_;
  double transformation_time_;
  double density_advection_time_;
  double CalculateEnergy();
  // Compute the derivative of energy for each spectrum.
  void ComputeEnergyDerivativePerFreq();
  void SetEnergy(double desired_e);
  void DissipateEnergy();
  void AddExternalForce();
  void AddDensity(const int x, const int y, const int width,
                  const int height, FIELD2D* field);
  void AddBuoyancy(const FIELD2D& field);
  void AdvectDensity();
  void AdvectParticles();
  void ProjectOutObstacle();
  void EigenAnalysis();
  // Init a vortices in the middle.
  void InitVelocityVortices();
  double maximum_condition_number_;
};

#endif  // LAPLACIAN_FLUID_H
