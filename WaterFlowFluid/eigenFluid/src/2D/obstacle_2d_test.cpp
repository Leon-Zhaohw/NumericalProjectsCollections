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
#include <math.h>
// #include <glog/logging.h>

#include "2D/all_dirichlet_basis_2D.h"
#include "2D/laplacian_basis_2D.h"
#include "2D/obstacle_2d.h"
#include "util/timer.h"

int main() {
  const int basis_dim_root = 20;
  const int basis_dim = basis_dim_root * basis_dim_root;
  const int xRes = 512;
  AllDirichletBasis2D diri(xRes, basis_dim_root);
  // Create an empty obstacle_2d object.
  Obstacle2D obstacle(xRes, xRes, 20);
  Eigen::MatrixXd obs_proj_matrix;
  obs_proj_matrix.resize(basis_dim, basis_dim);
  obs_proj_matrix.setRandom();
  obstacle.ComputeObstacleProjectionMatrix(diri, false, &obs_proj_matrix);
  for (int j = 0; j < basis_dim; j++) {
    for (int i = 0; i < basis_dim; i++) {
      if (obs_proj_matrix(i,j) != 0 && i != j) {
        std::cout << "obstacle_2d_test.cpp " << __LINE__ << " FATAL: " <<  "Should be 0." << std::endl; exit(0);
      }
    }
  }
  Timer timer;
  /*obs_proj_matrix.setRandom();
  timer.Reset();
  obstacle.ComputeObstacleProjectionMatrix(diri, 20.0, true, &obs_proj_matrix);
  std::cout <<  "Full rank time: " << timer.ElapsedTimeInSeconds() << std::endl;
  for (int j = 0; j < basis_dim; j++) {
    for (int i = 0; i < basis_dim; i++) {
      if (i == j) {
        std::cout <<  "Diag entry: " << obs_proj_matrix(i,j) << std::endl;
      } else {
        if (abs(obs_proj_matrix(i,j)) > 1e-7) {
          std::cout <<  "Should be zero."  << std::endl;
        }
      }
    }
  }*/
  obs_proj_matrix.setZero();
  obstacle.InitializeUsingImage("obstacles/boxes.jpg");
  std::cout <<  "Measuring time ... " << std::endl; 
  timer.Reset();
  obstacle.ComputeObstacleProjectionMatrix(diri, false, &obs_proj_matrix);
  std::cout <<  "Actually time " <<timer.ElapsedTimeInSeconds() << std::endl;
  
  return 0;
}
