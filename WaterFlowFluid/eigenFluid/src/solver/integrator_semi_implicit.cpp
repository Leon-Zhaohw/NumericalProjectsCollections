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
#include <vector>

#include "integrator_semi_implicit.h"
#include "setting.h"

void IntegratorSemiImplicit::IntegrateForward(
  const std::vector<Adv_Tensor_Type>& Adv_tensor,
  const Eigen::VectorXd& coefficients_old,
  const double dt,
  Eigen::VectorXd* coefficients_new,
  double* condition_number_,
  float* contraction_time,
  float* solver_time) {
  
  const int num_basis = coefficients_old.size();
  // Allocate the contracted tensor. Generally, it should be dense.
  Eigen::MatrixXd w_C(num_basis, num_basis);
  w_C.setIdentity();
  Eigen::VectorXd w(num_basis);
  w.setZero();
  w = Eigen::VectorXd::Map(coefficients_old.data(), num_basis);
  for (int i = 0; i < num_basis; i++) {
    w_C.row(i) -= w.transpose() * Adv_tensor[i] * dt; 
  }
  // Estimate the condition number.
  //Eigen::JacobiSVD<Eigen::MatrixXd> svd(w_C);
  //(*condition_number_) = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
  
  Eigen::VectorXd x = w_C.colPivHouseholderQr().solve(w);
  for (int i = 0; i < num_basis; i++) {
    (*coefficients_new)[i] = x(i);
  }
}
