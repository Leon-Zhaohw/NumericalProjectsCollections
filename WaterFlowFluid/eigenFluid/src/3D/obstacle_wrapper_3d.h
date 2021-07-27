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

#ifndef OBSTACLE_WRAPPER_3D_H
#define OBSTACLE_WRAPPER_3D_H

#include <vector>
#include <set>

#include "3D/obstacle_3d.h"
#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/laplacian_basis_set_3d.h"

class ObstacleWrapper3D{
public:
  ObstacleWrapper3D(const double dx, const std::string& obs_type, const std::string& obs_file,
                    const std::string& obstacle_list, const VEC3I& gridRes):
  dx_(dx), obstacle_type_(obs_type), obstacle_file_(obs_file), obstacle_list_(obstacle_list),
  gridRes_(gridRes){
    // InitializeObstacles();
    ReadObstaclesFromList();
    total_cells_ = gridRes_[0]*gridRes_[1]*gridRes_[2];
    check_ = (char*) malloc(sizeof(char)*total_cells_);
    memset(check_, 0x00, sizeof(char)*total_cells_);
  };
  ~ObstacleWrapper3D(){
    for (int i = 0; i < obstacles_.size(); i++) {
      delete obstacles_[i];
    }
    obstacles_.clear();
    for (std::set<OBSTACLE3D*>::iterator it = random_boxs_.begin(); it != random_boxs_.end(); it++) {
      delete *it;
    }
    random_boxs_.clear();
    
    free(check_);
  };
  void Draw() const;
  // Return the index of the obstacle this one is in.
  // ATTENTION: Assume the obstacles do not overlap, if they overlap, only
  // return the first one.
  bool inside(const VEC3F& point) const;
  void CalculateNormalForce(
      const VECTOR3_FIELD_3D& velocity, const double dt, const double cell_len, const double force_amp,
      VECTOR3_FIELD_3D* normal_component);
  // Zero out the velocity inside the obstacle of a given velocity field. 
  void ZeroOutVelocityDensityObstalce(const double cell_len, VECTOR3_FIELD_3D* velocity, FIELD_3D* density);
  
  // Extract the index of the cell grid which need to be handled for implicit obstacle handling.
  void ExtractProjObstacleIndex(const VEC3I& dim_, const double cell_len, 
                                std::vector<ObstacleIdx3D>* proj_indices);
  // Get the projection matrix for implict obstacle handling.
  // When the projection indices convering all the cells in the volume, the result projection
  // Matrix is identity.
  void ComputeObstacleProjectionMatrix(const LaplacianBasisSet3D& basis,
                                       const VEC3I& dim_, const double cell_len,
                                       const double force_amp,
                                       Eigen::MatrixXd* proj_matrix);
  
  void RasterToGrid(unsigned char* obstacle_filed);
  
  void InitializeObstacles();
  void ReadObstaclesFromList();
  
  void MoveObstacle(const double dt) {
    for (int i = 0; i < obstacles_.size(); i++) {
      obstacles_[i]->MoveObstacle(dt);
    }
    for (std::set<OBSTACLE3D*>::iterator it = random_boxs_.begin(); it != random_boxs_.end(); it++) {
      (*it)->MoveObstacle(dt);
    }
  }
  void WriteObstacleToPBRT(const std::string& foldername, const int frame_idx);
  bool GetTrigger();
  
  VEC3F getVelocityAt(const VEC3F& point) const;
  // pass in the delta of mouse on screen, screen size is normalized to [0, 1].
  void SetOmega(const float mouse_dx, const float mouse_dy,
                const float pos_x, const float pos_y, const double dt);
  void SetObstacleVelo(const VEC3F& velo);
  
protected:
  std::vector<OBSTACLE3D*> obstacles_;
  std::set<OBSTACLE3D*> random_boxs_;
  OBSTACLE3D* GenerateRandomBox();
  void GenerateRandomBoxes();
  
  char* check_;
  int total_cells_;
  // cell dimension,
  const double dx_;
  // Which kind of obstacle want...
  const std::string obstacle_type_;
  const std::string obstacle_file_;
  const std::string obstacle_list_;
  
  const VEC3I gridRes_;
};

#endif  // OBSTACLE_WRAPPER_3D_H
