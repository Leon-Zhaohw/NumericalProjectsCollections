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

#include "3D/four_neumann_basis_3d.h"
#include "3D/four_neumann_basis_set_3d.h"
#include "util/timer.h"
#include "util/util.h"
#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/FIELD_3D.h"
#include "setting.h"

namespace {
  
void AllocateBasisTest(std::vector<FourNeumannBasis3D>* all_basis_) {
  // -k1*A + k2*B - k3*C = 0. 
  // k1_ == 0 && k3_ == 0 zero basis.
  // k1_ == 0 // vy, vz == 0
  // k3_ == 0 // vx, vy == 0
  //Test purpose
  FourNeumannBasis3D basis(1,4,0,0,0,1);
  all_basis_->push_back(basis);
  
  FourNeumannBasis3D basis1(4,0,1,1,0,-4);
  all_basis_->push_back(basis1);
  
  FourNeumannBasis3D basis2(0,8,9,1,0,0);
  all_basis_->push_back(basis2);
  
  FourNeumannBasis3D basis3(12,14,9,-5,-3,2);
  all_basis_->push_back(basis3);
  
  FourNeumannBasis3D basis4(23,5,15,10,-8,-18);
  all_basis_->push_back(basis4);
  
  FourNeumannBasis3D basis5(7,1,0,0,0,1);
  all_basis_->push_back(basis5);
  
  FourNeumannBasis3D basis6(0,0,1,1,0,0);
  all_basis_->push_back(basis6);
  
  FourNeumannBasis3D basis7(0,0,0,1,0,0);
  all_basis_->push_back(basis7);
  FourNeumannBasis3D basis8(0,0,0,0,0,1);
  all_basis_->push_back(basis8);
}

void NormalizationTest1() {
  const int basis_dim = 9;
  const int xRes = 200, yRes = 266, zRes = 200;
  int total = xRes*yRes*zRes;
  
  VECTOR3_FIELD_3D vf1(xRes, yRes, zRes);
  VECTOR3_FIELD_3D vf2(xRes, yRes, zRes);
  FourNeumannBasisSet3D basis_set(basis_dim, xRes, yRes, zRes, "principle_x");
  Timer timer;
  
  std::vector<FourNeumannBasis3D> all_basis_;
  AllocateBasisTest(&all_basis_);
  basis_set.AllocateAllBasis(all_basis_);
  
  Eigen::VectorXd coef(basis_dim);
  coef.setZero();
  coef[0] = -1.5;
  coef[1] = 0.5;
  coef[2] = 0.3;
  coef[3] = -1;
  coef[4] = 2;
  coef[5] = -0.4;
  coef[6] = 1.2;
  coef[7] = 0;
  coef[8] = 0;
  Eigen::VectorXd coef1(basis_dim);
  Eigen::VectorXd coef2(basis_dim);
  coef1.setZero();
  coef2.setZero();
  
  timer.Reset();
  basis_set.FillBasisNumerical(coef, &vf1);
  // ATTENTION: FFTW_MEASURE will overwrite the data while planning, need to do one transform ahead.
  
  basis_set.InverseTramsformToVelocity(coef, &vf2);
  basis_set.InverseTramsformToVelocity(coef, &vf2);
  
  basis_set.ForwardTransformtoFrequency(vf1, &coef2);
  basis_set.ForwardTransformtoFrequency(vf1, &coef2);
  
  std::cout <<  "Time for neumerical: " << timer.ElapsedTimeInSeconds() << std::endl;
  
  std::cout <<  "Vnorm: " << vf1.twoNormSqared()*PI_CUBE/(xRes*yRes*zRes) << std::endl;
  std::cout <<  "CoefNorm: " << coef.squaredNorm() << std::endl;
  for (int i = 0; i < coef2.size(); i++) {
    std::cout <<  coef2[i] << std::endl;
  }
  Real diff_x = 0, diff_y = 0, diff_z = 0;
  Real vf1_nsquared = 0;
  Real vf2_nsquared = 0;
  
  for (int i = 0; i < total; i++) {
    Real temp = vf1[i][0] - vf2[i][0];
    diff_x += temp*temp;
    diff_y += (vf1[i][1] - vf2[i][1])*(vf1[i][1] - vf2[i][1]);
    diff_z += (vf1[i][2] - vf2[i][2])*(vf1[i][2] - vf2[i][2]);
    //std::cout << "s :  " << vf1[i][0] << " s: " << vf2[i][] 
  }
  
  std::cout <<  "X Velocity diff squared : " << diff_x << std::endl;
  std::cout <<  "Y Velocity diff squared : " << diff_y << std::endl;
  std::cout <<  "Z Velocity diff squared : " << diff_z << std::endl;
  
}

void NormalizationTest() {
  const int basis_dim_des = 1024;
  const int xRes = 128, yRes = 64, zRes = 64;
  
  FourNeumannBasisSet3D basis_set(basis_dim_des, xRes, yRes, zRes, "principle_x");
  
  int basis_dim = basis_set.AllocateAllBasis();
  Eigen::VectorXd coef(basis_dim);
  for (int i = 0; i < basis_dim; i++) {
    coef(i) = static_cast<double> (std::rand()) / RAND_MAX;
  }
  VECTOR3_FIELD_3D vf(xRes, yRes, zRes);
  basis_set.InverseTramsformToVelocity(coef, &vf);
  basis_set.InverseTramsformToVelocity(coef, &vf);
  
  std::cout <<  "Vnorm: " << vf.twoNormSqared()*PI_CUBE/(xRes*yRes*zRes) << std::endl;
  std::cout <<  "CoefNorm: " << coef.squaredNorm() << std::endl;
}

void TensorTest () {
  const int basis_dim = 400;
  const int xRes = 200, yRes = 200, zRes = 200;
  int total = xRes*yRes*zRes;
  
  FourNeumannBasisSet3D basis_set(basis_dim, xRes, yRes, zRes, "principle_x");
  
  basis_set.AllocateAllBasis();
  
  std::vector<Adv_Tensor_Type> C;
  basis_set.FillVariationalTensor(&C);
  basis_set.VerifyAntisymmetric(C);
}

}  // namespace

int main(int argc, char ** argv) {
//  google::ParseCommandLineFlags(&argc, &argv, true);
//  google::InitGoogleLogging(argv[0]);
  fftw_init_threads();
  fftw_plan_with_nthreads(2);
  
  // NormalizationTest1();
  // NormalizationTest();
  TensorTest();
  
  std::cout <<  "Finished." << std::endl;
  return 0;
}
