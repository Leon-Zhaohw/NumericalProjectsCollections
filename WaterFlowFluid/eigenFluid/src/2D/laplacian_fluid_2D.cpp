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

#include <cstdlib>
#include <ctime>
#include "Eigen"
#include "Sparse"
#include <fstream>
// #include <glog/logging.h>
#include <math.h>

#include "2D/all_dirichlet_basis_2D.h"
#include "2D/drawer_2d.h"
#include "2D/FIELD2D.h"
#include "solver/integrator_2d.h"
#include "solver/integrator_RK4.h"
#include "solver/integrator_semi_implicit.h"
#include "2D/laplacian_basis_2D.h"
#include "2D/laplacian_fluid_2D.h"
#include "2D/particle_2d.h"
#include "2D/three_dirichlet_one_neumann.h"
#include "2D/obstacle_2d.h"
#include "util/timer.h"
// #include "2D/two_neumann_x.h"
#include "2D/two_neumann_x_2D.h"
#include "util/util.h"
#include "util/read_write_tensor.h"

#include "setting.h"
#include "solver/trapezoidal.h"
#include "Alg/VEC2.h"
#include "2D/VFIELD2D.h"

void LaplacianFluid2D::Initialize() {
  // Initialize the density and velocity field.
  velocity_ = VFIELD2D(xRes_, yRes_);
  force_    = VFIELD2D(xRes_, yRes_);
  normal_velocity_ = VFIELD2D(xRes_, yRes_);
  density_  = FIELD2D(xRes_, yRes_);
  density_old_ = FIELD2D(xRes_, yRes_);
  
  temp1 = FIELD2D(xRes_, yRes_);
  temp2 = FIELD2D(xRes_, yRes_);
  vort_ = FIELD2D(xRes_, yRes_);

  timer_ = Timer();
  // Initialize the coefficients.
  basis_coefficients_.resize(basis_dim_);
  basis_coefficients_.setZero();
  TransferMatrix_.resize(basis_dim_, basis_dim_);
  TransferMatrix_.setIdentity();
  eigLargeVectors_.resize(basis_dim_, basis_dim_);
  eigLargeVectors_.setZero();
  // GetNormalizedRandomVector(basis_dim_, &basis_coefficients_);
  basis_coefficients_old_.resize(basis_dim_);
  basis_coefficients_old_.setZero();
  force_dw_.resize(basis_dim_);
  force_dw_.setZero();
  eigs_.resize(basis_dim_, 0.0);
  eigs_inv_.resize(basis_dim_, 0.0);
  eigs_inv_root_.resize(basis_dim_, 0.0);

  energy_derivative_.resize(basis_dim_);
  energy_derivative_.setZero();
  tempVector_.resize(basis_dim_);
  tempVector_.setZero();

  // Initialize drawer.
  drawer_.reset(new Drawer2D());
  
  // Initialize the integrator.
  if (integrator_type_ == "RK4") {
    integrator_.reset(new IntegratorRK4());
  } else if (integrator_type_ == "semi_implicit"){
    integrator_.reset(new IntegratorSemiImplicit());
  } else if (integrator_type_ == "trapezoidal") {
    integrator_.reset(new Trapezoidal());
  } else {
    std::cout << "laplacian_fluid_2D.cpp " << __LINE__ << " FATAL: " <<  "Unknow integrator type: " << integrator_type_ << std::endl; exit(0);
  } 
  // Initialize the basis.
  int basis_type_int;
  if (basis_type_ == "all_dirichlet") {
    basis_type_int = 0;
    basis_.reset(new AllDirichletBasis2D(xRes_, basis_dim_root_));
    is_neumann_ = false;
  } else if (basis_type_ == "three_dirichlet_one_neumann") {
    basis_type_int = 1;
    basis_.reset(new ThreeDirichletOneNeumann(xRes_, basis_dim_root_));
    is_neumann_ = true;
  } else if (basis_type_ == "two_neumann_x") {
    basis_type_int = 2;
    basis_.reset(new TwoNeumannX2D(xRes_, basis_dim_root_));
    is_neumann_ = true;
  } else {
    std::cout << "laplacian_fluid_2D.cpp " << __LINE__ << " FATAL: " <<  "Unknow basis tpye: " << basis_type_ << std::endl; exit(0);
  }
  basisWeights_.resize(basis_dim_);
  basisWeights_.setOnes();
  
  // Fill the eigen values.
  basis_.get()->ComputeEigenValues(&eigs_, &eigs_inv_, &eigs_inv_root_);
  // Fill the tensor.
  if (tensor_fname_.size() == 0) {
    std::cout <<  "Fill the advection tensor..." << std::endl;
    basis_.get()->FillVariationalTensor(&Adv_tensor_);
  } else {
    std::cout <<  "Read tensor from file." << std::endl;
    ReadTensor(tensor_fname_, basis_dim_, basis_type_int, &Adv_tensor_);
    basis_.get()->ComputeBasisWeights(&basisWeights_);
    basis_.get()->MultiplyweightTensor(Adv_tensor_, basisWeights_);
  }

  for (int i = 0; i < num_particles_; i++) {
    int px = std::rand() % xRes_;
    int py = std::rand() % yRes_;
    Particle2D particle(px + 0.5, py + 0.5);
    particles_.push_back(particle);
  }
  rest_frame_ = total_frame_;
  integrator_time_ = 0.;
  transformation_time_ = 0.;
  density_advection_time_ = 0.;
  maximum_condition_number_ = 0.;
  frame_simulated_ = 0;


}

void LaplacianFluid2D::InitializeObstacles(const ObstacleParams& param_) {
  std::cout <<  "Initialize the obstacle..." << std::endl;
  if (!param_.use_obstacles) {
    use_obstacles_ = false;
    return;
  } else  {
    obstacle_.reset(new Obstacle2D(xRes_, yRes_, param_.obstacle_force_scale));
    // obstacle_.get()->InitializeUsingImage(param_.img_file_name);
    use_obstacles_ = true;
    handle_obstacle_implicit_ = param_.handle_obstacle_implicit;
  }
  obs_proj_mat_.resize(basis_dim_, basis_dim_);
  obstacle_.get()->ComputeObstacleProjectionMatrix((*basis_.get()),
  false, &obs_proj_mat_);
}

void LaplacianFluid2D::Step() {
  float DCT_time = 0;
  float obstacle_time = 0;
  float integrator_time = 0;
  float advection_time = 0;
  float contraction_time = 0;
  float solver_time = 0;
  
  current_energy_ = CalculateEnergy();

  // Swap the coefficients.
  for (int i = 0; i < basis_coefficients_.size(); i++) {
    basis_coefficients_old_[i] = basis_coefficients_[i];
  }

  double condition_number = 1.0;
  //ComputeEnergyDerivativePerFreq();
  integrator_.get()->IntegrateForward(Adv_tensor_, basis_coefficients_old_,
                                      dt_,
                                      &basis_coefficients_, &condition_number,
                                      &contraction_time, &solver_time);
  

  //ComputeEnergyDerivativePerFreq();


  if (basis_type_ == "two_neumann_x" || basis_type_ == "two_neumann_x_new") {
     // ZeroOutUnWantedCoefficients();
  }
  if (condition_number > maximum_condition_number_) {
    maximum_condition_number_ = condition_number;
  }
  integrator_time = contraction_time + solver_time;
 
  // Disspate energy from viscosity.
  DissipateEnergy();
  // Project out of obstacle, explicit method.
  if (use_obstacles_ && !handle_obstacle_implicit_) {
    ProjectOutObstacle();
  }
  // Add buoyancy by using DCT to transform to freq domain.
  timer_.Reset();
  if (addBounyancyOnce && frame_simulated_ == 0) {
    AddBuoyancy(density_);
  }
  if (!addBounyancyOnce) {
    AddBuoyancy(density_);
  }
  DCT_time += timer_.ElapsedTimeInSeconds();
  // Add external forces.
  AddExternalForce();
  if (basis_type_ == "two_neumann_x" || basis_type_ == "two_neumann_x_new") {
  //   ZeroOutUnWantedCoefficients();
  }
  timer_.Reset();
  // Implicit method for handling obstacles. TODO: Handle moving obstacles.
  if (use_obstacles_ && handle_obstacle_implicit_) {
    obstacle_.get()->SolveObstacleImplicit(basis_coefficients_,
    obs_proj_mat_, &basis_coefficients_old_);
    basis_coefficients_ = basis_coefficients_old_;
  }
  obstacle_time += timer_.ElapsedTimeInSeconds();
  
  // Reconstruct the velocity field.
  timer_.Reset();
  basis_.get()->InverseTramsformToVelocity(basis_coefficients_, &velocity_);
  DCT_time += timer_.ElapsedTimeInSeconds();
  
  // Handle obstacles.
  if (use_obstacles_) {
    obstacle_.get()->ZeroOutVelocityDensityObstalce(&velocity_, &density_);
  }
  timer_.Reset();
  // Advect density.
  AdvectDensity();
  // Move the obstacle.
  if (use_obstacles_) {
    obstacle_.get()->MoveForward(dt_/dx_);
  }
  // Advect particles.
  AdvectParticles();
  advection_time = timer_.ElapsedTimeInSeconds();
  rest_frame_ --;
  if (rest_frame_ <=0) {
    Quit();
  }
  // convert to vort field for visualization
  // basis_.get()->ComputeVorticity(basis_coefficients_, &vort_);
  
  std::cout <<  "Frame: " << frame_simulated_ << std::endl;
  std::cout <<  "Total energy: " << current_energy_ << std::endl;

  std::cout <<  "DCT time: " << DCT_time << std::endl;
  std::cout <<  "obstacle time: " << obstacle_time << std::endl;
  std::cout <<  "Advection time: " << advection_time << std::endl;
  
  integrator_time_ += integrator_time;
  density_advection_time_ += advection_time;  
  transformation_time_ += DCT_time;
  
  frame_simulated_ ++;
}

void LaplacianFluid2D::AdvectDensity() {
  const double dt0 = dt_ / dx_;
  density_.swapPointer(density_old_);
  // VFIELD2D::advect(dt0, velocity_, density_old_, density_);
  VFIELD2D::advectMacCormack(dt0, velocity_, density_old_, density_, temp1, temp2);
  temp1.clear();
  temp2.clear();
  density_.setZeroBorder();
}

void LaplacianFluid2D::AddSmokeTestCase(
    const int xpos, const int ypos, const int width, const int height) {
  AddDensity(xpos, ypos, width, height, &density_);  
}

void LaplacianFluid2D::AddDensity(const int x, const int y, const int width,
                                const int height, FIELD2D* field) {
  int startx = x - width / 2;
  int starty = y - height / 2;
  int endx = x + width / 2;
  int endy = y + height / 2;
  
  if (startx < 0) startx = 0;
  if (startx >= xRes_) startx = xRes_ - 1;
  if (starty < 0) starty = 0;
  if (starty >= yRes_) starty = yRes_ - 1;
  if (endx < 0) endx = 0;
  if (endx >= xRes_) endx = xRes_ - 1;
  if (endy < 0) endy = 0;
  if (endy >= yRes_) endy = yRes_ - 1;
  
  for(int j = starty; j <= endy; j++) {
    for (int i = startx; i <= endx; i ++) {
      int index = i + j*xRes_;
      (*field)[index] = 1.f;
    }
  }
}

void LaplacianFluid2D::AddForceImpluse(const int x, const int y, const int width,
                                       const int height, const double mag) {
  int startx = x - width / 2;
  int starty = y - height / 2;
  int endx = x + width / 2;
  int endy = y + height / 2;
  
  if (startx < 0) startx = 0;
  if (startx >= xRes_) startx = xRes_ - 1;
  if (starty < 0) starty = 0;
  if (starty >= yRes_) starty = yRes_ - 1;
  if (endx < 0) endx = 0;
  if (endx >= xRes_) endx = xRes_ - 1;
  if (endy < 0) endy = 0;
  if (endy >= yRes_) endy = yRes_ - 1;
  for(int j = starty; j <= endy; j++) {
    for (int i = startx; i <= endx; i ++) {
      int index = i + j*xRes_;
      force_[index][1] += mag / dt_;
    }
  }  
}

void LaplacianFluid2D::InitVelocityFromDedalus(const std::string& Ufname, const std::string& Vfname,
                                               const int dedalusXRes, const int dedalusYRes) {
  // VFIELD2D v(xRes_, yRes_);
  // ReadVelocityFromDedalus(dedalusYRes, dedalusYRes, Ufname, Ufname, v);
  // basis_.get()->ForwardTransformtoFrequency(v, &basis_coefficients_);
  
}

void LaplacianFluid2D::AddBuoyancy(const FIELD2D& field) {
  int index = 0;
  for (int y = 0; y < yRes_; y++) {
    for (int x = 0; x < xRes_ ; x++) 
    {
      index = x + y*xRes_;
      force_[index][1] += buoyancy_ * field[index];
    }
  }
}

void LaplacianFluid2D::ProjectOutObstacle() {
  obstacle_.get()->CalculateNormalForce(velocity_, dt_, &force_);
}

double LaplacianFluid2D::CalculateEnergy() {
  // Calculate current energy, sum of squares of coefficients
  // since laplacian eigenfunction basis is orthogonal.
  double energy=0.0;
  for (int i = 0; i < basis_dim_; i++)
#ifndef USE_ROOT_LAMBDA
    energy += eigs_inv_[i] * (basis_coefficients_[i] * basis_coefficients_[i]);
#else
    energy += (basis_coefficients_[i] * basis_coefficients_[i]);
#endif
  return energy;
}

void LaplacianFluid2D::SetEnergy(double desired_e) {
  double cur_e = CalculateEnergy();
  double fact = sqrt(desired_e)/sqrt(cur_e);
  for (int i = 0;i < basis_dim_; i++) 
    basis_coefficients_[i] *= fact;
}

void LaplacianFluid2D::DissipateEnergy() {
  for (int i = 0; i < basis_dim_; i++) {
#ifndef USE_ROOT_LAMBDA
    basis_coefficients_[i] *= exp(-1.0 * eigs_[i] * dt_ * visc_);
#else
    basis_coefficients_[i] *= exp(-1.0 * eigs_[i] * dt_ * visc_);
#endif
  }
}
// Convert the external force field into the coefficients, and then
// add the external forces.
void LaplacianFluid2D::AddExternalForce() {
  // Convert the force field to coefficients.
  Eigen::VectorXd force_coef(basis_dim_); force_coef.setZero();
  basis_.get()->ForwardTransformtoFrequency(force_, &force_coef);
  //basis_.get()->ProjectVelocituNumerical(force_, &force_coef);
  force_.clear();
  for (int i = 0; i < basis_dim_; i++) {
    force_dw_[i] += force_coef[i];
  }
  for (int i = 0; i < basis_dim_; i++) {
    basis_coefficients_[i] += force_dw_[i] * dt_;
    force_dw_[i] = 0.;
  }
}

void LaplacianFluid2D::DrawDensity() {
  drawer_.get()->DrawDensity(density_, xRes_, yRes_, dx_);
}

void LaplacianFluid2D::DrawVelocity() {
  drawer_.get()->DrawVelocity(velocity_, xRes_, yRes_, dx_);
}

void LaplacianFluid2D::DrawCoefficients(const double multi_factor) {

  Eigen::VectorXd energy_w = basis_coefficients_.cwiseProduct(basis_coefficients_);
  drawer_.get()->DrawCoefficients(energy_w,  dx_, multi_factor*20.0, 0);
//drawer_.get()->Draw2Coefficients(drawVec1, drawVec2, dx_, multi_factor*3.0, 0.3);
  
//drawer_.get()->DrawCoefficients(drawVec2, dx_, multi_factor*10.0, -0.3);

}

void LaplacianFluid2D::DrawParticles(const double ptl_length) {
  drawer_.get()->DrawParticles(particles_, dx_, dt_, ptl_length);
}

void LaplacianFluid2D::DrawObstacles() {
  if (use_obstacles_) {
    drawer_.get()->DrawObstacles(*obstacle_.get(), dx_);
  }
}

void LaplacianFluid2D::DrawVort() {
  drawer_.get()->DrawVort(vort_, xRes_, yRes_, dx_);
}

void LaplacianFluid2D::MoveObstacles(const double delta_x, const double delta_y) {
  if (use_obstacles_) {
    obstacle_.get()->MoveObstacles(delta_x, delta_y);
  } else {
    std::cout <<  "obstacles handling is disabled." << std::endl;
  } 
}

void LaplacianFluid2D::AddVelocityToObstacles(const double dvx, const double dvy) {
  if (use_obstacles_) {
  obstacle_.get()->AddVelocity(dvx, dvy);
  } else {
    std::cout <<  "obstacle handling is disabled." << std::endl;
  }
}

void LaplacianFluid2D::AddForce(const int xpos, const int ypos,
                              const float fx, const float fy) {
  const int index = xpos + ypos * xRes_;
  force_[index][0] += fx;
  force_[index][1] += fy;
}

void LaplacianFluid2D::AddSmoke(const int xpos, const int ypos, 
                              const float amount){
  const int index = xpos + ypos * xRes_;
  density_[index] += amount;
}

void LaplacianFluid2D::AdvectParticles() {
   const double dt0 = dt_ / dx_;
  for (int i = 0; i < num_particles_; i++) {
    // Get the velocity at the particles.
    const VEC2& position = particles_[i].position;
    VEC2 p_v = velocity_.GetVelocity(position[0], position[1]);
    // Forward Eular.
    particles_[i].position[0] += p_v[0] * dt0;
    particles_[i].position[1] += p_v[1] * dt0;
    particles_[i].velocity = p_v;

    if (particles_[i].position[0] < 0. || particles_[i].position[0] > xRes_ || 
      particles_[i].position[1] < 0. || particles_[i].position[1] > yRes_ ) {
      particles_[i].position[0] = std::rand() % xRes_ + 0.5;
      particles_[i].position[1] = std::rand() % yRes_ + 0.5;
      particles_[i].velocity = 0.;
    }
  }
}


void LaplacianFluid2D::ReSeedParticles() {
  // Reseed the particles at random position.
  for (int i = 0; i < num_particles_; i++) {
    int px = std::rand() % xRes_;
    int py = std::rand() % yRes_;
    particles_[i].position[0] = px + 0.5;
    particles_[i].position[1] = py + 0.5;
    particles_[i].velocity = 0.;
  }
}

void LaplacianFluid2D::Quit() {
  std::cout <<  "Total frame simulated: " << total_frame_ << std::endl;
  std::cout << "Average time for integrator per frame: " 
      << integrator_time_ / total_frame_;
  std::cout << "Average time for DCT and DST: "
      << transformation_time_ / total_frame_;
  std::cout << "Average time for advection of density and particles: "
      << density_advection_time_ / total_frame_;
  std::cout <<  "Maximum condition number: " << maximum_condition_number_ << std::endl;
  quit_ = true;
 // exit(0);
}

void LaplacianFluid2D::ZeroOutUnWantedCoefficients() {
  for (int i = 0; i < basis_dim_; i++) {
    int i1 = basis_.get()->LookUpBasis(i, 0) - 1;
    int i2 = basis_.get()->LookUpBasis(i, 1) - 1;
    if (i2 == 0 && i1 > 0) {
      basis_coefficients_[i] = 0.;
    }
  }
}

void LaplacianFluid2D::ComputeEnergyDerivative() {
  Eigen::MatrixXd C_W(basis_dim_, basis_dim_);
  C_W.setZero();
  for (int i = 0; i < basis_dim_; i++) {
    C_W.row(i) = basis_coefficients_.transpose() * Adv_tensor_[i];
  }
  double derivative = basis_coefficients_.transpose() * C_W * basis_coefficients_;
  std::cout <<  "Energy derivative: " << derivative << std::endl;
}

// Compute the energy derivative. Which is C*basis_coefficients_.
void LaplacianFluid2D::ComputeEnergyDerivativePerFreq() {

  Eigen::VectorXd v(basis_dim_);
  v.setZero();
  Eigen::VectorXd temp(basis_dim_);
  temp.setZero();
  for (int i = 0; i < basis_dim_; i++) {
    v[i] = basis_coefficients_old_.transpose() * Adv_tensor_[i] * basis_coefficients_;
    // energy_derivative_[i] = basis_coefficients_[i]*v[i];
    temp[i] = basis_coefficients_old_[i]*v[i];
  }
  energy_derivative_ = temp;
  //for (int i = 0; i < basis_dim_; i++) {
  //  energy_derivative_[i] = basis_coefficients_[i]*basis_coefficients_[i] - basis_coefficients_old_[i]*basis_coefficients_old_[i];
  //}
  double el = 0, ed = 0;
  for (int i = basis_dim_ - 51; i < basis_dim_; i++) {
    el += basis_coefficients_old_[i]*basis_coefficients_old_[i];
    ed += temp[i];
  }
  // std::cout <<  basis_coefficients_old_[1660] << std::endl;
  std::fstream f("test.txt", std::ios::out | std::ios::app);
  f << el << " " << ed << "\n";
  f.close();
  std::cout <<  "low freq energy: " << el << " " << ed << std::endl;
}

void LaplacianFluid2D::InitCheckBoardPattern() {
  // Init a checkboard pattern smoke.
  double checkSize = 0.1;
  double delta = 0.15;
  
  int iCheckSize = static_cast<int>(checkSize*xRes_);
  int iDelta = static_cast<int>(delta*xRes_);
  
  for (int y = 0; y < yRes_; y++) {
    for (int x = 0; x < xRes_; x++) {
      int ix = x / iCheckSize;
      int iy = y / iCheckSize;
      const int index = x + y*xRes_;
      if ((ix + iy) % 2 == 0) {
        // Fill the block with smoke.
        density_[index] = 1.0;
      } else {
        density_[index] = 0.0;
      }
    }
  }
  /*
  for (int y = 0; y < yRes_; y++) {
    for (int x = 0; x < xRes_; x++) {
      int ix = x / iCheckSize;
      int iy = (y + iDelta) / iCheckSize;
      const int index = x + y*xRes_;
      if (iy % 2 == 0) {
        // Fill the block with smoke.
        density_[index] = 1.0;
      } else {
        density_[index] = 0.0;
      }
    }
  }*/
}

void LaplacianFluid2D::InitVelocityVortices() {
  // domain [0,1]^2
  double theta = 0.02, cx = 0.4, cy = 0.5, mag = 0.04;
  // stream function Phi = exp(-[(x - cx)^2 + (y - cy)^2] / theta);
  // double* phi = new double[xRes_*yRes_];
  
  for (int y = 0; y < yRes_; y++) {
    for (int x = 0; x < xRes_; x++) {
      double xx = (static_cast<double>(x) + 0.5)*dx_;
      double yy = (static_cast<double>(y) + 0.5)*dx_;
      double nt = (- (xx - cx)*(xx - cx) - (yy - cy)*(yy - cy)) / theta;
      double phi = exp(nt);
      velocity_[x + y*xRes_][0] = - mag * phi / theta * 2.0*(yy - cy);
      velocity_[x + y*xRes_][1] = mag * phi / theta * 2.0*(xx - cx);
    }
  }
  cx = 1.0 - cx;
  // cy = 1.0 - cy;
  // second vortices.
  for (int y = 0; y < yRes_; y++) {
    for (int x = 0; x < xRes_; x++) {
      double xx = (static_cast<double>(x) + 0.5)*dx_;
      double yy = (static_cast<double>(y) + 0.5)*dx_;
      double nt = (- (xx - cx)*(xx - cx) - (yy - cy)*(yy - cy)) / theta;
      double phi = exp(nt);
      velocity_[x + y*xRes_][0] += - mag * phi / theta * 2.0*(yy - cy);
      velocity_[x + y*xRes_][1] += mag * phi / theta * 2.0*(xx - cx);
    }
  }
  // Init the coefficients vector using this.
  basis_.get()->ForwardTransformtoFrequency(velocity_, &basis_coefficients_);
  velocity_.clear();
}

void LaplacianFluid2D::EigenAnalysis() {
  Eigen::MatrixXd P(basis_dim_, basis_dim_);
  P.setZero();
  #pragma omp parallel for
  for (int i = 0; i < basis_dim_; i++) {
    P.row(i).noalias() = basis_coefficients_.transpose() * Adv_tensor_[i]; 
  }
  
  Eigen::EigenSolver<Eigen::MatrixXd> eigs(P + P.transpose());
  Eigen::VectorXd Re(basis_dim_);
  int index = 0; double largest_eig = -10000;
  int j = 0; int k = 0;
  Eigen::MatrixXd curUMatrix(basis_dim_, basis_dim_);
  curUMatrix.setZero();
  Eigen::MatrixXd negUMatrix(basis_dim_, basis_dim_);
  negUMatrix.setZero();
  for (int i = 0; i < basis_dim_; i++) {
     Re[i] = eigs.eigenvalues()[i].real();
     if (Re[i] > largest_eig) {
       index = i;
       largest_eig = Re[i];
     }
     if (Re[i] > 0) {
       curUMatrix.col(j).noalias() = eigs.eigenvectors().col(i).real();
       j++;
     }
     if (Re[i] < 0) {
       negUMatrix.col(k).noalias() = eigs.eigenvectors().col(i).real();
       k++;
     }
  }
  Eigen::VectorXcd a = eigs.eigenvectors().col(largest_eig);
  for (int i = 0; i < basis_dim_; i++) {
    tempVector_[i] = a[i].real();
  }
  std::cout <<  "larger than zero eig vectors: " << j << std::endl;
  std::cout << "Increase part: " << (basis_coefficients_.transpose()*curUMatrix).squaredNorm() / basis_coefficients_.squaredNorm() << " Decrease part: " <<
  (basis_coefficients_.transpose()*negUMatrix).squaredNorm() / basis_coefficients_.squaredNorm() ;
  double ecur = basis_coefficients_.squaredNorm();
  double eprev = basis_coefficients_old_.squaredNorm(); 
  std::cout <<  "Energy increase rate: " << (ecur - eprev) / ecur << std::endl;
  
  eigLargeVectors_ = curUMatrix;
}

void LaplacianFluid2D::ComputeBoundaryIntegral() {
  // Compute the integral at x = pi.
  double integral = 0.;
  for (int i = 0; i < yRes_; i++) {
    int index = i * xRes_ + xRes_ - 1;
    double vx = velocity_[index][0];
    integral += vx*vx*vx;
  }
  std::cout << "The total energy flow at boundary: " << integral;
}


