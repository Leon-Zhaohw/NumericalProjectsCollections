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

#ifndef OBSTACLE_2D_H
#define OBSTACLE_2D_H

#include "Eigen"
#include <cstring>
#include <stdlib.h>
#include <vector>

#include "2D/laplacian_basis_2D.h"
#include "2D/FIELD2D.h"
#include "2D/VFIELD2D.h"

enum OBSTACLE_FLAGS {
  EMPTY = 0, 
  MARCHED = 2, 
  RETIRED = 4 
};  

struct ObstacleParams{
  ObstacleParams() {
    img_file_name.clear();
    use_obstacles = false;
    obstacle_force_scale = 20.;
  }
  std::string img_file_name;
  bool use_obstacles;
  double obstacle_force_scale;
  bool handle_obstacle_implicit;
};

// Store the coordinate of the particular position where the speed need to 
// be handled because of obstacles (Along a narrow band around obstacles).
struct ProjObstacleIdx {
  ProjObstacleIdx(const int x_, const int y_, const int mode_):
      x(x_),y(y_),mode(mode_){}
  int x;
  int y;
  // mode = 0, vx need to be handled. mode = 1, vy need to be handled. mode = 2
  // both need to be handled. mode = 3, 45 degree / obstacles.
  // mode = 4, 45 degree \ obstacles.
  int mode;
};

class Obstacle2D {
public:
  Obstacle2D(int xRes, int yRes, const double amp_):xRes_(xRes), yRes_(yRes)
  , force_amp_(amp_) {
    obstacles = (char*) malloc(sizeof(char)*xRes*yRes);
    std::memset(obstacles, 0x00, sizeof(char)*xRes*yRes);
    dx_ = 0.;
    dy_ = 0.;
    vx_ = 0.;
    vy_ = 0.;
  };
  void InitializeObstacles();
  // Calculate the normal component of a given velocity field along the 
  // border of the obstacles. This function won't clear the content of
  // normal_component. Instead it will add the addtion value to it.
  void CalculateNormalForce(
      const VFIELD2D& velocity, const double dt, VFIELD2D* normal_component);
  // Zero out the velocity inside the obstacle of a given velocity field. 
  void ZeroOutVelocityDensityObstalce(VFIELD2D* velocity, FIELD2D* density);
  ~Obstacle2D(){};
  void InitializeAsBox(const int topleft_x, const int topleft_y,
                       const int rightbottom_x, const int rightbottom_y);
  // Use an image to initialize the obstacle.
  // void InitializeUsingImage(const std::string& img_name);
  inline bool WithinRange(const int x, const int y);
  // Check whether a position has obstacle. The displacement of the obstacle
  // will be considered.
  const char Access(const int pos_x, const int pos_y) const;
  void MoveObstacles(const double delta_x, const double delta_y) {
    dx_ += delta_x;
    dy_ += delta_y;
  };
  void AddVelocity(const double vx, const double vy) {
    vx_ += vx;
    vy_ += vy;
  }
  void MoveForward(const double dt) {
    dx_ += vx_ * dt;
    dy_ += vy_ * dt;
  }
  void ComputeObstacleProjectionMatrix(const LaplacianBasis2D& basis, bool test,
                                       Eigen::MatrixXd* proj_matrix);
  const int xRes_;
  const int yRes_;
  void GenerateTestProjIndices(std::vector<ProjObstacleIdx>* proj_indices);
  void SolveObstacleImplicit(const Eigen::VectorXd& omega_in,
                             const Eigen::MatrixXd& proj_matrix,
                             Eigen::VectorXd* omega_out);
protected:
  void ExtractProjObstacleIndex(std::vector<ProjObstacleIdx>* proj_indices);
  double vx_;
  double vy_;
  double dx_;
  double dy_;
  const double force_amp_;
  char* obstacles;
};

#endif
