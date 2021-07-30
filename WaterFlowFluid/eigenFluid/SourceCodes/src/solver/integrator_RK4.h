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

#ifndef INTEGRATOR_RK4_H
#define INTEGRATOR_RK4_H

#include "Eigen"
#include "Sparse"
#include <vector>

#include "integrator_2d.h"
// Runge Kutta 4 explicit integrator.
class IntegratorRK4 : public Integrator2D {
public:
  IntegratorRK4() : Integrator2D() {};
  ~IntegratorRK4(){};
  void IntegrateForward(const std::vector<Adv_Tensor_Type>& Adv_tensor,
                        const Eigen::VectorXd & coefficients_old,
                        const double dt,
                        Eigen::VectorXd* coefficients_new,
                        double* condition_number_,
                        float* contraction_time,
                        float* solver_time) override;
protected:
  
};

#endif  //INTEGRATOR_RK4_H
