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

#ifndef INTEGRATOR_SEMI_IMPLICIT_H
#define INTEGRATOR_SEMI_IMPLICIT_H

#include "Eigen"
#include "Sparse"
#include <vector>

#include "integrator_2d.h"
#include "setting.h"

// We call this integrator 'semi' implicit. The complete implicit formulation
// contains quadratic terms of unknows in the right hand side of the equation
// and linear terms on left hand side of the equation. By replace part of the
// unknow with the known coefficients, we reduced the quadratic terms to linear.

class IntegratorSemiImplicit : public Integrator2D {
public:
  IntegratorSemiImplicit() : Integrator2D() {};
  ~IntegratorSemiImplicit(){};
  void IntegrateForward(const std::vector<Adv_Tensor_Type>& Adv_tensor,
                        const Eigen::VectorXd& coefficients_old,
                        const double dt,
                        Eigen::VectorXd* coefficients_new,
                        double* condition_number_,
                        float* contraction_time,
                        float* solver_time) override;
private:
};

#endif  // INTEGRATOR_SEMI_IMPLICIT_H
