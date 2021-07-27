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
// #include <glog/logging.h>
#include <vector>
#include <iostream>

#include "integrator_RK4.h"
#include "setting.h"
// Solve the equation dw_k/dt = wT C_k w
void IntegratorRK4::IntegrateForward(
    const std::vector<Adv_Tensor_Type>& Adv_tensor,
    const Eigen::VectorXd& coefficients_old,
    const double dt,
    Eigen::VectorXd* coefficients_new,
    double* condition_number_ ,
    float* contraction_time,
    float* solver_time) {
    /* RK4: higher order explicit integrator */
  const int num_coefficient = coefficients_old.size();
  if (num_coefficient != coefficients_new->size()) {std::cout << "integrator_RK4.cpp " << __LINE__ << " FATAL: " << "The dimention of the coefficients must agree."  << std::endl; exit(0);}
  if (num_coefficient != Adv_tensor[0].cols()) {std::cout << "integrator_RK4.cpp " << __LINE__ << " FATAL: " <<  "The dimention of the coefficients must be same as the tensor" << std::endl; exit(0);}
  if (num_coefficient != Adv_tensor.size()) {std::cout << "integrator_RK4.cpp " << __LINE__ << " FATAL: " << "The dimention of the coefficients must be same as the tensor" << std::endl; exit(0);}
  if (num_coefficient != Adv_tensor[0].rows()) {std::cout << "integrator_RK4.cpp " << __LINE__ << " FATAL: "  << "The dimention of the coefficients must be same as the tensor" << std::endl; exit(0);}
  std::vector<Eigen::VectorXd> qn;
  std::vector<Eigen::VectorXd> dwt;
  for(int i = 0; i < 4; i++) {
    Eigen::VectorXd q(num_coefficient);
    q.setZero();
    qn.push_back(q);
    Eigen::VectorXd dw(num_coefficient);
    dw.setZero();
    dwt.push_back(dw);
  }
  qn[0] = Eigen::VectorXd::Map(coefficients_old.data(), num_coefficient);
  Eigen::VectorXd dw(num_coefficient);
  // dwt[0] = RHS.
  // qn[1] = w_t + RHS * 0.5dt.
  for (int k = 0; k < num_coefficient; k++) {
    dwt[0](k) = qn[0].transpose() * Adv_tensor[k] * qn[0]; // k1
    qn[1](k) = qn[0](k) + 0.5*dt*dwt[0](1);
  } 
  for (int k = 0; k < num_coefficient; k++) {
    dwt[1](k) = qn[1].transpose() * Adv_tensor[k] * qn[1];  // k2
    qn[2](k) = qn[0](k) + 0.5*dt * dwt[1](k);	
  }
  for (int k = 0; k < num_coefficient; k++) {
    dwt[2](k) = qn[2].transpose() * Adv_tensor[k] * qn[2];  // k3
    qn[3](k) = qn[0](k) + dt * dwt[2](k);	
  }
  for (int k = 0; k < num_coefficient; k++) {
    dwt[3](k) = qn[3].transpose() * Adv_tensor[k] * qn[3]; // k4
    dw(k) = (dwt[0](k) + 2*dwt[1](k) + 2*dwt[2](k) + dwt[3](k))/6.0;
  }
  for (int k = 0; k < num_coefficient; k++) {
    (*coefficients_new)[k] = coefficients_old[k] + dw(k) * dt;
  }
}
