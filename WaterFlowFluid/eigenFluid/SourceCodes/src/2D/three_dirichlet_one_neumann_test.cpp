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

#include "2D/three_dirichlet_one_neumann.h"
#include "util/timer.h"
#include "2D/VFIELD2D.h"

int main() {
  const int basis_dim_root = 20;
  const int basis_dim = basis_dim_root * basis_dim_root;
  const int xRes = 128;
  
  Timer timer_;
  //fluid = new FLUID2D(xRes, yRes);
  VFIELD2D* basis_v = new VFIELD2D(xRes, xRes);
  VFIELD2D* basis_f = new VFIELD2D(xRes, xRes);
  ThreeDirichletOneNeumann oneNeu(xRes, basis_dim_root);
  Eigen::VectorXd basis_coefficients(basis_dim);
  basis_coefficients.setZero();
  for (int i = 0; i < basis_dim / 2; i++) {
    basis_coefficients[i] = 1.0;
  }
  for (int i = basis_dim / 2; i < basis_dim; i++) {
    basis_coefficients[i] = -0.5;
  }
  Eigen::VectorXd inv_basis_coefficients(basis_dim);
  inv_basis_coefficients.setZero();
  timer_.Reset();
  oneNeu.FillBasisNumerical(basis_coefficients, basis_v);
  std::cout <<  "time11: " << timer_.ElapsedTimeInSeconds() << std::endl;
  
  oneNeu.ForwardTransformtoFrequency(*basis_v, &inv_basis_coefficients);
 
  timer_.Reset();
  oneNeu.InverseTramsformToVelocity(basis_coefficients, basis_f);
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
  oneNeu.FillVariationalTensor(&C);
  oneNeu.VerifyAntisymmetric(C);
  std::cout <<  "VerifyAntisymmetric finished." << std::endl;
  return 0;
}
