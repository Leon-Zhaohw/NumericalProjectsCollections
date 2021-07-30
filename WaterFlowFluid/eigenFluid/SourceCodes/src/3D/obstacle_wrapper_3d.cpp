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

#include <fstream>
#include <omp.h>

#include "3D/obstacle_wrapper_3d.h"
#include "3D/obstacle_box_3d.h"
#include "3D/obstacle_sphere_3d.h"
#include "3D/obstacle_cylinder_3d.h"
#include "3D/obstacle_cylinder_static_3d.h"
#include "3D/obstacle_cylinder_interact.h"

#include "util/stringprintf.h"
#include "util/util.h"

OBSTACLE3D* ObstacleWrapper3D::GenerateRandomBox() {
  //const float r1 = GenerateUnifomRnd();
  //const float r2 = GenerateUnifomRnd();
  const float r3 = GenerateUnifomRnd();
  const VEC3F center(0.1, 0.5*0.8+0.1, 0.5*0.8+0.1);
  // const VEC3F velo(,0,0);
  const VEC3F velo(0.1 , 0, 0);
  const VEC3F lenght(0.1,0.3,0.3);
  const VEC3F color(0.1, 0.2, 0.8);
  
  OBSTACLE3D* obs= new ObstacleBox3D(center, lenght, 90, 90, 0, VEC3F(0,0,1), gridRes_, dx_);
  // std::cout <<  obs << std::endl;
  obs->SetVelocity(velo);
  obs->set_color(color);
  return obs;
}


void ObstacleWrapper3D::Draw() const {
  for (int i = 0; i < obstacles_.size(); i++) {
    obstacles_[i]->draw();
  }
  for (std::set<OBSTACLE3D*>::iterator it = random_boxs_.begin(); it != random_boxs_.end(); it++) {
    (*it)->draw();
  }
}

bool ObstacleWrapper3D::inside(const VEC3F& point) const {
  const int num_obstacle = obstacles_.size();
  
  for (int i = 0; i < num_obstacle; i++) {
    if (obstacles_[i]->inside(point)) {
      return true;
    }
  }
  for (std::set<OBSTACLE3D*>::iterator it = random_boxs_.begin(); it != random_boxs_.end(); it++) {
    if ((*it)->inside(point))
      return true;
  }
  // Not in all the obstacles.
  return false;
}

void ObstacleWrapper3D::InitializeObstacles() {
  if (obstacle_type_ == "box") {
    // A box
    VEC3F cen(0.4,0.24,0.26);
    VEC3F obslen(0.25, 0.25, 0.25);
    float angle = 0.0f; VEC3F axis(0,1,0);
    ObstacleBox3D *box = NULL;
    if (obstacle_file_.size() == 0) {
      box = new ObstacleBox3D(cen, obslen , 90,90, angle, axis, gridRes_, dx_);
    } else {
      box = new ObstacleBox3D(obstacle_file_, gridRes_, dx_);
    }
    obstacles_.push_back(box);
  } else if (obstacle_type_ == "sphere") {
     // A sphere.
     VEC3F center(0.5,0.25,0.25); float radius = 0.1f;
     ObstacleSphere3D *sphere = NULL;
     if (obstacle_file_.size() == 0) {
       sphere = new ObstacleSphere3D(center, radius, gridRes_, dx_);
     } else {
       sphere = new ObstacleSphere3D(obstacle_file_, gridRes_, dx_);
     }
     obstacles_.push_back(sphere);
  } else if (obstacle_type_ == "cylinder") {
    // A cylinder.
    float radius = 0.2f; float height = 0.1f; VEC3F center(0.2,0.24,0.26);
    ObstacleCylinder3D *cylinder = new ObstacleCylinder3D(center, radius, height,
                                                          5, true, gridRes_, dx_);
    // cylinder->SetVelocity(VEC3F( -0.2,0.0,0.0));
    obstacles_.push_back(cylinder);
  } else if (obstacle_type_ == "cylinder_static") {
    float radius = 0.1f; float height = 0.5f; VEC3F center(0.5,0.25,0.25);
    float angle = 0.0f; VEC3F axis(0,1,0);
    ObstacleCylinderStatic3D* cylinder_static = NULL;
    if (obstacle_file_.size() == 0) {
      cylinder_static = new ObstacleCylinderStatic3D(center, radius,
                                                     height,angle, axis, gridRes_, dx_);
    } else {
      cylinder_static = new ObstacleCylinderStatic3D(obstacle_file_, gridRes_, dx_);
    }
    obstacles_.push_back(cylinder_static);
  } else {
    std::cout << "obstacle_wrapper_3d.cpp " << __LINE__ << " ERROR: " <<  "Unknow obstacle type." << std::endl;
  }
}

void ObstacleWrapper3D::ReadObstaclesFromList() {
  std::ifstream in(obstacle_list_);
  if (!in.is_open()) {
    std::cout << "obstacle_wrapper_3d.cpp " << __LINE__ << " FATAL: " <<  "Cannot open: " << obstacle_list_ << std::endl; exit(0);
  }
  do {
    std::string obs_type;
    std::string filename;
    in >> obs_type; in >> filename;
    if (obs_type == "box") {
      ObstacleBox3D *box = new ObstacleBox3D(filename, gridRes_, dx_);
      obstacles_.push_back(box);
    } else if (obs_type == "sphere") {
      ObstacleSphere3D *sphere = new ObstacleSphere3D(filename, gridRes_, dx_);
      obstacles_.push_back(sphere);
    } else if (obs_type == "cylinder_static") {
      ObstacleCylinderStatic3D* cylinder_static = new ObstacleCylinderStatic3D(filename, gridRes_, dx_);
      obstacles_.push_back(cylinder_static);
    } else if (obs_type == "cylinder") {
      ObstacleCylinder3D* cylinder = new ObstacleCylinder3D(filename, gridRes_, dx_);
      obstacles_.push_back(cylinder);
    } else if (obs_type == "cylinder_interact") {
      ObstacleCylinderInterac3D* cylinder_interact = new ObstacleCylinderInterac3D(filename, gridRes_, dx_);
      obstacles_.push_back(cylinder_interact);
    } else {
      std::cout << "obstacle_wrapper_3d.cpp " << __LINE__ << " " <<  "Unknow type: " << obs_type << "Ingnored." << std::endl;
    }
  } while (!in.eof());
}

void ObstacleWrapper3D::CalculateNormalForce(
      const VECTOR3_FIELD_3D& velocity, const double dt, const double cell_len, const double force_amp,
      VECTOR3_FIELD_3D* normal_component) {
  memset(check_, 0x00, sizeof(char)*total_cells_);
  for (int i = 0; i < obstacles_.size(); i++) {
    obstacles_[i]->CalculateNormalForce(velocity, check_, dt, cell_len, force_amp, normal_component);
  }
  for (std::set<OBSTACLE3D*>::iterator it = random_boxs_.begin(); it != random_boxs_.end(); it++) {
    (*it)->CalculateNormalForce(velocity, check_, dt, cell_len, force_amp, normal_component);
  }
}

// Zero out the velocity inside the obstacle of a given velocity field. 
void ObstacleWrapper3D::ZeroOutVelocityDensityObstalce(const double cell_len,
                                    VECTOR3_FIELD_3D* velocity, FIELD_3D* density) {
  for (int i = 0; i < obstacles_.size(); i++) {
    obstacles_[i]->ZeroOutVelocityDensityObstalce(cell_len, velocity, density);
  }
  for (std::set<OBSTACLE3D*>::iterator it = random_boxs_.begin(); it != random_boxs_.end(); it++) {
    (*it)->ZeroOutVelocityDensityObstalce(cell_len, velocity, density);
  }
}

void ObstacleWrapper3D::RasterToGrid(unsigned char* obstacle_filed) {
  // obstacle_filed is preallocated.
  memset(obstacle_filed, 0x00, sizeof(char)*total_cells_);
  for (int i = 0; i < obstacles_.size(); i++) {
    obstacles_[i]->RasterToGrid(obstacle_filed);
  }
}

void ObstacleWrapper3D::WriteObstacleToPBRT(const std::string& foldername, const int frame_idx) {
 // std::string fname = StringPrintf("%sObs%04d", foldername.c_str(), frame_idx);
 // for (int i = 0; i < obstacles_.size(); i++) {
 //   obstacles_[i]->OutPutGeometryToPBRT(fname);
 // }
  
  // Write to a single file instead.
  std::string fname = StringPrintf("%sObs%04d.pbrt", foldername.c_str(), frame_idx);
  std::ofstream out(fname);
  if (!out.is_open()) {
    std::cout << "obstacle_wrapper_3d.cpp " << __LINE__ << " FATAL: " <<  "Cannot open: " << fname << std::endl; exit(0);
  }
  for (int i = 0; i < obstacles_.size(); i++) {
    obstacles_[i]->OutPutGeometryToPBRT(out);
  }
  for (std::set<OBSTACLE3D*>::iterator it = random_boxs_.begin(); it != random_boxs_.end(); it++) {
    (*it)->OutPutGeometryToPBRT(out);
  }
  out.close();
}

bool ObstacleWrapper3D::GetTrigger() {
  for (int i = 0; i < obstacles_.size(); i++) {
    // obstacles_[i]->OutPutGeometryToPBRT(out);
    if (obstacles_[i]->GetTrigger()) {
      return true;
    }
  }
  
  return false;
}

void ObstacleWrapper3D::ComputeObstacleProjectionMatrix(const LaplacianBasisSet3D& basis,
                                     const VEC3I& dim_, const double cell_len,
                                     const double force_amp,
                                     Eigen::MatrixXd* proj_matrix) {
  const int xRes = dim_[0];
  const int yRes = dim_[1];
  const int zRes = dim_[2];
  const int slab = xRes*yRes;
  
  const double normFactor = PI_CUBE/(xRes*yRes*zRes);
  
  const int basis_dim = basis.GetBasisDim();
  proj_matrix->setZero();
  if (proj_matrix->cols() != basis_dim) {std::cout << "obstacle_wrapper_3d.cpp " << __LINE__ << " FATAL: "  << "The dimention must agree." << std::endl; exit(0);}
  if (proj_matrix->rows() != basis_dim) {std::cout << "obstacle_wrapper_3d.cpp " << __LINE__ << " FATAL: "  << "The dimention must agree." << std::endl; exit(0);}
  std::vector<ObstacleIdx3D> indices;
  ExtractProjObstacleIndex(dim_, cell_len, &indices);
  const int indices_count = indices.size();
  
  std::cout <<  "Total number of basis: " << indices_count << std::endl;
  std::cout << "Percentage filled: " << 
      indices_count / (static_cast<float> (xRes*yRes*zRes));
  // Iterate over the entries of the matrix.
  #pragma omp parallel for
  for (int j = 0; j < basis_dim; j++) {
    for (int i = 0; i < basis_dim; i++) {
      double value  = 0.;
      // Iterate over all the position where the normal need to be handled.
      for (int k = 0; k < indices_count; k++) {
        const double vi = basis.ComputeBasisAt(indices[k].idx, i, indices[k].mode);
        const double vj = basis.ComputeBasisAt(indices[k].idx, j, indices[k].mode);
        value += vi*vj;
      }
      (*proj_matrix)(i,j) = value * normFactor;
    }
  }
}

void ObstacleWrapper3D::SetOmega(const float mouse_dx, const float mouse_dy,
                                 const float pos_x, const float pos_y, const double dt) {
  for (int i = 0; i < obstacles_.size(); i++) {
    obstacles_[i]->setOmega(mouse_dx, mouse_dy, pos_x, pos_y, dt);
  }
}

void ObstacleWrapper3D::SetObstacleVelo(const VEC3F& velo) {
  for (int i = 0; i < obstacles_.size(); i++) {
    obstacles_[i]->SetVelocity(velo);
  }
}

VEC3F ObstacleWrapper3D::getVelocityAt(const VEC3F& point) const {
  VEC3F velo_sum(0,0,0);
  // If the point hit inside of multiple obstacles....
  for (int i = 0; i < obstacles_.size(); i++) {
    velo_sum += obstacles_[i]->getVelocityAt(point);
  }
  return velo_sum;
}

void ObstacleWrapper3D::ExtractProjObstacleIndex(const VEC3I& dim_, const double cell_len, 
                              std::vector<ObstacleIdx3D>* proj_indices) {
   proj_indices->clear();
  const int xRes = dim_[0];
  const int yRes = dim_[1];
  const int zRes = dim_[2];
  const int slab = xRes*yRes;

  for (int k = 0; k < zRes; k++) {
    for (int j = 0; j < yRes; j++) {
      for (int i = 0; i < xRes; i++) {
        const int index = i + j*xRes + k*slab;
        
        VEC3F pt((static_cast<float>(i) + 0.5)*cell_len,
                (static_cast<float>(j) + 0.5)*cell_len,
                (static_cast<float>(k) + 0.5)*cell_len);
        // Cell is inside the obstacle.
        if (this->inside(pt)) continue;
        int  left = 0, right = 0, top = 0, bottom = 0, front = 0, back = 0;
        // left
        pt[0] -= cell_len;
        if (this->inside(pt)) left ++;
        pt[0] += cell_len;
        // right
        pt[0] += cell_len;
        if (this->inside(pt)) right ++;
        pt[0] -= cell_len;
        // top
        pt[1] += cell_len;
        if (this->inside(pt)) top++;
        pt[1] -= cell_len;
        // bottom
        pt[1] -= cell_len;
        if (this->inside(pt)) bottom++;
        pt[1] += cell_len;
        // front
        pt[2] += cell_len;
        if (this->inside(pt)) front++;
        pt[2] -= cell_len;
        // back
        pt[2] -= cell_len;
        if (this->inside(pt)) back++;
        pt[2] += cell_len;
        int total = top + bottom + left + right + front + back;
        // No obstacle nearby.
        if (total == 0) continue;
        // Left side. x component.
        if (left != 0) {
          proj_indices->push_back(ObstacleIdx3D(i,j,k,0));
        }
        // right side. x component.
        if (right != 0) {
          proj_indices->push_back(ObstacleIdx3D(i,j,k,0));
        }
        // top, y component.
        if (top != 0) {
          proj_indices->push_back(ObstacleIdx3D(i,j,k,1));
        }
        // bottom, y component.
        if (bottom != 0) {
          proj_indices->push_back(ObstacleIdx3D(i,j,k,1));
        }
        // front, z component.
        if (front != 0) {
          proj_indices->push_back(ObstacleIdx3D(i,j,k,2));
        }
        // back, z component.
        if (back != 0) {
          proj_indices->push_back(ObstacleIdx3D(i,j,k,2));
        }
      }
    }
  }
}
