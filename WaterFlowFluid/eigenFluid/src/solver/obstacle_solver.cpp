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
#include "obstacle_solver.h"

void SolveObstacleImplicit(const Eigen::VectorXd& omega_in,
        const Eigen::MatrixXd& proj_matrix, Eigen::VectorXd* omega_out) {
  
  Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
  solver.compute(proj_matrix);
  (*omega_out) = solver.solve(omega_in);
}
