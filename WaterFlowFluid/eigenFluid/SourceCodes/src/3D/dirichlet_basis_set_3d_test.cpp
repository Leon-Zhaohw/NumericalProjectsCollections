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
#include "Sparse"
#include <fftw3.h>
#include <omp.h>
#include <math.h>
#include <vector>

#include "3D/dirichlet_basis_3d.h"
#include "3D/dirichlet_basis_set_3d.h"
#include "util/timer.h"
#include "util/util.h"
#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/FIELD_3D.h"
#include "setting.h"

namespace {
 
void AllocateBasisTest(std::vector<DirichletBasis3D>* all_basis_) {
  //Test purpose
  DirichletBasis3D basis(10,9,15,-12,-5,11);
  all_basis_->push_back(basis);
  
  DirichletBasis3D basis2(3,2,5,-3,2,1);
  all_basis_->push_back(basis2);
  
  DirichletBasis3D basis3(0,1,1,0,-1,1);
  all_basis_->push_back(basis3);
  
  DirichletBasis3D basis1(1,0,1,-1,0,1);
  all_basis_->push_back(basis1);
  
  DirichletBasis3D basis4(1,1,0,-1,1,0);
  all_basis_->push_back(basis4);
}
 
void TransformationTest() {
  const int basis_dim = 5;
  const int xRes = 200, yRes = 133, zRes = 100;
  int total = xRes*yRes*zRes;
  
  VECTOR3_FIELD_3D vf1(xRes, yRes, zRes);
  VECTOR3_FIELD_3D vf2(xRes, yRes, zRes);
  DirichletBasisSet3D basis_set(basis_dim, xRes, yRes, zRes, "principle_x");
  
  std::vector<DirichletBasis3D> all_basis_;
  AllocateBasisTest(&all_basis_);
  basis_set.AllocateAllBasis(all_basis_);
  
  Timer timer;
  
  Eigen::VectorXd coef(basis_dim);
  coef.setZero();
  coef[0] = 0.1;
  coef[1] = -1;
  coef[2] = -1.3;
  coef[3] = 1.5;
  coef[4] = 2.67;
  
  Eigen::VectorXd coef1(basis_dim);
  coef1.setZero();
  
  timer.Reset();
  basis_set.FillBasisNumerical(coef, &vf1);
  std::cout <<  "Time for neumerical: " << timer.ElapsedTimeInSeconds() << std::endl;
  
  basis_set.InverseTramsformToVelocity(coef, &vf2);
  timer.Reset();
  basis_set.InverseTramsformToVelocity(coef, &vf2);
  std::cout <<  "Timer for DCT: " << timer.ElapsedTimeInSeconds() << std::endl;
  basis_set.ForwardTransformtoFrequency(vf1, &coef1);
  timer.Reset();
  basis_set.ForwardTransformtoFrequency(vf1, &coef1);
  std::cout <<  "Timer for DCT: " << timer.ElapsedTimeInSeconds() << std::endl;
   
  Real diff_x = 0, diff_y = 0, diff_z = 0;
  Real vf1_nsquared = 0;
  Real vf2_nsquared = 0;
  for (int i = 0; i < total; i++) {
    Real temp = vf1[i][0] - vf2[i][0];
    diff_x += temp*temp;
    diff_y += (vf1[i][1] - vf2[i][1])*(vf1[i][1] - vf2[i][1]);
    diff_z += (vf1[i][2] - vf2[i][2])*(vf1[i][2] - vf2[i][2]);
  }
  std::cout <<  "X Velocity diff squared : " << diff_x / static_cast<double>(total) << std::endl;
  std::cout <<  "Y Velocity diff squared : " << diff_y / static_cast<double>(total) << std::endl;
  std::cout <<  "Z Velocity diff squared : " << diff_z / static_cast<double>(total) << std::endl;
  for (int i = 0; i < basis_dim; i++) {
    std::cout <<  coef1[i] << std::endl;
  }
  std::cout <<  "Vnorm: " << vf2.twoNormSqared()*PI_CUBE/(xRes*yRes*zRes) << std::endl;
  std::cout <<  "CoefNorm: " << coef.squaredNorm() << std::endl;
  // std::cout <<  "Vf1 norm squared: " << vf1_nsquared << "  Vf2 norm squared: " << vf2_nsquared << std::endl;
  // std::cout <<  "Maxv: " << vf2.maxAbsScalar() << std::endl;
}

void TensorTest () {
  const int basis_dim = 400;
  const int xRes = 200, yRes = 266, zRes = 200;
  int total = xRes*yRes*zRes;
  
  DirichletBasisSet3D basis_set(basis_dim, xRes, yRes, zRes, "principle_x");
  
  basis_set.AllocateAllBasis();
  
  std::vector<Adv_Tensor_Type> C;
  basis_set.FillVariationalTensor(&C);
  basis_set.VerifyAntisymmetric(C);
}

// Test the time for numerical computing the bases.
void TimingTest () {
  const int basis_dim = 300;
  const int xRes = 220, yRes = 220, zRes = 220;
  int total = xRes*yRes*zRes;
  
  DirichletBasisSet3D basis_set(basis_dim, xRes, yRes, zRes, "principle_x");

  int num_basis = basis_set.AllocateAllBasis();
  std::cout <<  "Total number of basis: " << num_basis << std::endl;
  Eigen::VectorXd coef(num_basis);
  for (int i = 0; i < num_basis; i++) {
    coef(i) = static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
  }
  
  VECTOR3_FIELD_3D vf(xRes, yRes, zRes);
  
  Timer timer;
  timer.Reset();
  basis_set.FillBasisNumerical(coef, &vf);
  std::cout <<  "Elaps: " << timer.ElapsedTimeInSeconds() << std::endl;
}

// This will eat tones of memory. Test the time for cached basis.
void CachedTest() {
  const int basis_dim = 200;
  const int xRes = 220, yRes = 220, zRes = 220;
  int64_t num_velo = xRes*yRes*zRes*3;
  
  Eigen::MatrixXd bases(num_velo, basis_dim); bases.setRandom();
  Eigen::VectorXd v(num_velo);
  Eigen::VectorXd coef(basis_dim); coef.setRandom();
  
  Timer timer;
  timer.Reset();
  for (int i = 0; i < 10; i++) {
  v = bases*coef;
  }
  std::cout <<  "Elaps: " << timer.ElapsedTimeInSeconds() / 10.0 << std::endl;
}

}  // namespace

int main(int argc, char ** argv) {
//  google::ParseCommandLineFlags(&argc, &argv, true);
//  google::InitGoogleLogging(argv[0]);
  fftw_init_threads();
  fftw_plan_with_nthreads(8);
  
  TransformationTest();
  TensorTest();
  //TimingTest();
  // CachedTest();
  
  std::cout <<  "Finished." << std::endl;
  
  return 0;
}
