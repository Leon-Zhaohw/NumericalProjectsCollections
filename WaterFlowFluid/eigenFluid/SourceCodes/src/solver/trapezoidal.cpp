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
#include <omp.h>
#include <vector>

#include "solver/integrator_2d.h"
#include "solver/trapezoidal.h"
#include "util/timer.h"

#include "Eigenvalues"

void Trapezoidal::IntegrateForward(const std::vector<Adv_Tensor_Type>& Adv_tensor,
                        const Eigen::VectorXd & coefficients_old,
                        const double dt,
                        Eigen::VectorXd* coefficients_new,
                        double* condition_number_,
                        float* contraction_time,
                        float* solver_time) {
  Timer timer;
  const int num_basis = coefficients_old.size();
  // Allocate the contracted tensor. Generally, it should be dense.
  Eigen::MatrixXd w_C(num_basis, num_basis);
  w_C.setIdentity();
  Eigen::MatrixXd P(num_basis, num_basis);
  P.setZero();
  Eigen::VectorXd w(num_basis);
  w.setZero();
  w = Eigen::VectorXd::Map(coefficients_old.data(), num_basis);
  timer.Reset();
  #pragma omp parallel for
  for (int i = 0; i < num_basis; i++) {
    P.row(i).noalias() = w.transpose() * Adv_tensor[i]; 
  }
  // P = TransferMatrix.transpose() * P;
  *contraction_time = timer.ElapsedTimeInSeconds();
  
  w_C = w_C - P*(0.5*dt);
  w  = w + P*(0.5*dt*w);
  
  Eigen::DiagonalPreconditioner<Eigen::MatrixXd::RealScalar> precond(w_C);
  int maxIters = 200;
  double tol_error = 1e-10;
  timer.Reset();
  conjugate_gradient_NE(w_C, w_C.transpose()*w, *coefficients_new, precond, maxIters, tol_error);
  *solver_time = timer.ElapsedTimeInSeconds();
  // std::cout <<  "maxIter: " << maxIters << " tolErr: " << tol_error << std::endl;
  
}
