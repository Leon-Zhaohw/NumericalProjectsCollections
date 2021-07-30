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
#include <gtest/gtest.h>
// #include <glog/logging.h>
#include <math.h>
#include <vector>

#include "2D/all_dirichlet_basis_2D.h"
#include "util/timer.h"
#include "setting.h"
#include "2D/VFIELD2D.h"

int main() {
  const int basis_dim_root = 20;
  const int basis_dim = basis_dim_root * basis_dim_root;
  const int xRes = 512;
  Timer timer_;
  VFIELD2D* basis_v = new VFIELD2D(xRes, xRes);
  VFIELD2D* basis_f = new VFIELD2D(xRes, xRes);
  AllDirichletBasis2D diri(xRes, basis_dim_root);
  // Test transformations.
  Eigen::VectorXd basis_coefficients(basis_dim);
  basis_coefficients.setZero();
  for (int i = 0; i < 7; i++) {
    basis_coefficients[i] = 1.0;
  }
  for (int i = 7; i < basis_dim; i++) {
    basis_coefficients[i] = -0.5;
  }
  
  Eigen::VectorXd inv_basis_coefficients(basis_dim);
  inv_basis_coefficients.setZero();
  timer_.Reset();
  diri.FillBasisNumerical(basis_coefficients, basis_v);
  std::cout <<  "time11: " << timer_.ElapsedTimeInSeconds() << std::endl;
  
  diri.ForwardTransformtoFrequency(*basis_v, &inv_basis_coefficients);
  
  basis_f = new VFIELD2D(xRes, xRes);
  
  timer_.Reset();
  diri.InverseTramsformToVelocity(basis_coefficients, basis_f);
  std::cout <<  "time22: " << timer_.ElapsedTimeInSeconds() << std::endl;
  
  float diff = 0.0;
  int index = 0;
  for (int j = 0;j < xRes; j++) {
    for (int i = 0; i < xRes; i++) {
      float temp = ((*basis_f)[index][1] - (*basis_v)[index][1]);
      diff += temp * temp;
      temp = ((*basis_f)[index][0] - (*basis_v)[index][0]);
      diff += temp * temp;
      index ++;
    }
  }
  std::cout <<  "Basis root mean squared difference: " << sqrtf(diff) / xRes << std::endl;
  diff = 0;
  for (int i = 0; i < basis_dim; i++) {
    float temp = basis_coefficients[i] - inv_basis_coefficients[i];
    diff += (temp) * (temp);
  }
  std::cout <<  "Frequency root mean squared difference:" << sqrtf(diff) / basis_dim_root << std::endl;
  
  // Test the tensor.
  std::vector<Adv_Tensor_Type> C;
  diri.FillVariationalTensor(&C);
  diri.VerifyAntisymmetric(C);
  std::cout <<  "VerifyAntisymmetric finished." << std::endl;
  basis_coefficients.setZero();
  basis_coefficients[0] = 1.0;
  basis_v->clear();
  diri.FillBasisNumerical(basis_coefficients, basis_v);
  Eigen::VectorXd proj_coef(basis_dim);
  proj_coef.setZero();
  diri.ProjectVelocityNumerical(*basis_v, &proj_coef);
  std::cout <<  "Factor: " << proj_coef(0) << std::endl;
  return 0;
}
