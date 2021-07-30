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

// #include <CImg.h>
#if defined(Success)
#undef Success
#endif

#include "Eigen"
// #include <glog/logging.h>
#include <math.h>
#include <omp.h>
#include <vector>

#include "2D/laplacian_basis_2D.h"
#include "2D/obstacle_2d.h"

// This is the explicit method.
// We don't need to take care the border of the box. As the basis already
// solved for them. The new value will be added to normal_component.
void Obstacle2D::CalculateNormalForce(const VFIELD2D& velocity,
                                      const double dt,
                                      VFIELD2D* normal_component) {
  if(velocity.getxRes() != normal_component->getxRes()) { std::cout << "FATAL: " << __LINE__ << "The dimention"
      << " must be equal" << std::endl; exit(0);}
  if(velocity.getyRes() != normal_component->getyRes()) { std::cout << "FATAL: " << __LINE__ << "The dimention"
      << " must be equal" << std::endl; exit(0);}
  const float amp = force_amp_ / dt;
  for (int j = 0; j < yRes_; j++) {
    for (int i = 0; i < xRes_; i++){
      int index = i + j * xRes_;
      if (Access(i, j)) continue;
      int up_box = 0;
      int down_box = 0;
      int left_box = 0;
      int right_box = 0;
      // Upper box.
      if (WithinRange(i, j - 1)) {
        if (Access(i, j - 1)) up_box = 1;
      }
      if (WithinRange(i, j + 1)) {
        if (Access(i, j + 1)) down_box = 1;
      }
      if (WithinRange(i - 1, j)) {
        if (Access(i - 1, j)) left_box = 1;
      }
      if (WithinRange(i + 1, j)) {
        if (Access(i + 1, j)) right_box = 1;
      }
      int total = up_box + down_box + left_box + right_box;
      // No obstacles around.
      if (total == 0) continue;
      // Obstacles equal or more than three, velocity should be zero.
      if (total >= 3) {
        (*normal_component)[index][1] -= amp*(velocity[index][1] - vy_);
        (*normal_component)[index][0] -= amp*(velocity[index][0] - vx_);
      }
      // One obstacles.
      if (total == 1) {
        if (up_box || down_box) {
        (*normal_component)[index][1] -= amp*(velocity[index][1] - vy_);
//          (*normal_component)[index][0] = 0.;
        }
        if (left_box || right_box) {
          (*normal_component)[index][0] -= amp*(velocity[index][0] - vx_);
//          (*normal_component)[index][1] = 0.;
        }
      }
      // Two obstacles.
      if (total == 2) {
        if (up_box && down_box) {
          (*normal_component)[index][1] -= amp*(velocity[index][1] - vy_);
//          (*normal_component)[index][0] = 0.;
        }
        if (left_box && right_box) {
          (*normal_component)[index][0] -= amp*(velocity[index][0] - vx_);
//          (*normal_component)[index][1] = 0.;
        }
        if ((up_box && right_box) || (down_box && left_box)) {
          (*normal_component)[index][0] -= (velocity[index][0] - velocity[index][1])*0.5*amp;
          (*normal_component)[index][1] -= (velocity[index][1] - velocity[index][0])*0.5*amp;
        }
        if ((up_box && left_box) || (down_box && right_box)) {
          (*normal_component)[index][0] -= (velocity[index][0] + velocity[index][1])*0.5*amp;
          (*normal_component)[index][1] -= (velocity[index][0] + velocity[index][1])*0.5*amp;
        }
      }
    }
  }
}

// Extract the index where the velocity around the obstacles need to be handled.
void Obstacle2D::ExtractProjObstacleIndex(std::vector<ProjObstacleIdx>* proj_indices) {
  proj_indices->clear();
  for (int j = 0; j < yRes_; j++) {
    for (int i = 0; i < xRes_; i++){
      if (Access(i, j)) continue;
      int up_box = 0;
      int down_box = 0;
      int left_box = 0;
      int right_box = 0;
      // Upper box.
      if (WithinRange(i, j - 1)) {
        if (Access(i, j - 1)) up_box = 1;
      }
      if (WithinRange(i, j + 1)) {
        if (Access(i, j + 1)) down_box = 1;
      }
      if (WithinRange(i - 1, j)) {
        if (Access(i - 1, j)) left_box = 1;
      }
      if (WithinRange(i + 1, j)) {
        if (Access(i + 1, j)) right_box = 1;
      }
      int total = up_box + down_box + left_box + right_box;
      // No obstacles around.
      if (total == 0) continue;
      // Obstacles equal or more than three, velocity should be zero.
      if (total >= 3) {
        ProjObstacleIdx proj_idx(i,j,2);
        proj_indices->push_back(proj_idx);
      }
      // One obstacles.
      if (total == 1) {
        if (up_box || down_box) {
          // Handle vy.
          ProjObstacleIdx proj_idx(i,j,1);
          proj_indices->push_back(proj_idx);
        }
        if (left_box || right_box) {
          // Handle vx.
          ProjObstacleIdx proj_idx(i,j,0);
          proj_indices->push_back(proj_idx);
        }
      }
      // Two obstacles.
      if (total == 2) {
        if (up_box && down_box) {
          // Handle y.
          ProjObstacleIdx proj_idx(i,j,1);
          proj_indices->push_back(proj_idx);
        }
        if (left_box && right_box) {
          // Handle x.
          ProjObstacleIdx proj_idx(i,j,0);
          proj_indices->push_back(proj_idx);
        }
        // \ mode.
        if ((up_box && right_box) || (down_box && left_box)) {
          ProjObstacleIdx proj_idx(i,j,4);
          proj_indices->push_back(proj_idx);
        }
        // / mode.
        if ((up_box && left_box) || (down_box && right_box)) {
          ProjObstacleIdx proj_idx(i,j,3);
          proj_indices->push_back(proj_idx);
        }
      }
    }
  }
}

void Obstacle2D::GenerateTestProjIndices(std::vector<ProjObstacleIdx>* proj_indices) {
  proj_indices->clear();
  for (int j = 0; j < yRes_; j++) {
    for (int i = 0; i < xRes_; i++){
      ProjObstacleIdx proj_idx(i,j,2);
      proj_indices->push_back(proj_idx);
    }
  }
}

// This is the implicit method. TODO: Optimize this function. Currently is
// too slow for computing on the fly.
void Obstacle2D::ComputeObstacleProjectionMatrix(
    const LaplacianBasis2D& basis, bool test,
    Eigen::MatrixXd* proj_matrix) {
  const double inv_xRes = 1.0 / xRes_;
  const double inv_yRes = 1.0 / yRes_;
  const int num_basis = basis.GetBasisNum();
  proj_matrix->setZero();
  if (proj_matrix->cols() != num_basis) {std::cout << "obstacle_2d.cpp " << __LINE__ << " FATAL: "  << "The dimention must agree." << std::endl; exit(0);}
  if (proj_matrix->rows() != num_basis) {std::cout << "obstacle_2d.cpp " << __LINE__ << " FATAL: "  << "The dimention must agree." << std::endl; exit(0);}
  std::vector<ProjObstacleIdx> indices;
  if (!test)
    ExtractProjObstacleIndex(&indices);
  else 
    GenerateTestProjIndices(&indices);
  const int indices_count = indices.size();
  std::cout <<  "Total number of basis: " << indices.size() << std::endl;
  std::cout << "Percentage filled: " << 
      indices.size() / (static_cast<float> (xRes_) * yRes_);
  // Iterate over the entries of the matrix.
#pragma omp parallel for
  for (int j = 0; j < num_basis; j++) {
    for (int i = 0; i < num_basis; i++) {
      double value  = 0.;
      // Iterate over all the position where the normal need to be handled.
      for (int k = 0; k < indices_count; k++) {
        // Vx. or both.
        if (indices[k].mode == 0 || indices[k].mode == 2) {
          const double vi = basis.ComputeBasisAt(indices[k].x, indices[k].y,
              i, 0);
          const double vj = basis.ComputeBasisAt(indices[k].x, indices[k].y,
              j, 0);
          value += vi*vj;
        }
        // Vy. or both.
        if (indices[k].mode == 1 || indices[k].mode == 2) {
          const double vi = basis.ComputeBasisAt(indices[k].x, indices[k].y,
              i, 1);
          const double vj = basis.ComputeBasisAt(indices[k].x, indices[k].y,
              j, 1);
          value += vi*vj;
        }
      }
      // TODO: handle the / and \ case.
      // Normalization.
      (*proj_matrix)(i,j) = value * inv_xRes * inv_yRes * 4 * force_amp_;
    }
  }
  // Add the indetity.
  for (int i = 0; i < num_basis; i++) {
    (*proj_matrix)(i,i) += 1.0;
  }
}

// Solve a equantion where proj_matrix*omega_out = omega_in;
void Obstacle2D::SolveObstacleImplicit(const Eigen::VectorXd& omega_in,
    const Eigen::MatrixXd& proj_matrix, Eigen::VectorXd* omega_out){
  Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
  solver.compute(proj_matrix);
  (*omega_out) = solver.solve(omega_in);
}

// We don't need to take care the border of the box. As the basis already
// solved for them.
void Obstacle2D::ZeroOutVelocityDensityObstalce(VFIELD2D* velocity, FIELD2D* density) {
  if (velocity->getxRes() != xRes_) {std::cout << "obstacle_2d.cpp " << __LINE__ << " FATAL: "  << "The dimention must equal" << std::endl; exit(0);}
  if (velocity->getyRes() != yRes_) {std::cout << "obstacle_2d.cpp " << __LINE__ << " FATAL: "  << "The dimention must equal" << std::endl; exit(0);}
  if (density->getxRes() != xRes_) {std::cout << "obstacle_2d.cpp " << __LINE__ << " FATAL: "  << "The dimention must equal" << std::endl; exit(0);}
  if (density->getyRes() != yRes_) {std::cout << "obstacle_2d.cpp " << __LINE__ << " FATAL: "  << "The dimention must equal" << std::endl; exit(0);}
  
  int index = 0;
  for (int j = 0; j < yRes_; j++) {
    for (int i = 0; i < xRes_; i++) {
      if(Access(i, j)) {
        (*velocity)[index] = 0;
        (*density)[index] = 0.;
      }
      index ++;
    }
  }
}

void Obstacle2D::InitializeAsBox(const int topleft_x, const int topleft_y,
                                 const int rightbottom_x, const int rightbottom_y) {
  if (topleft_x < 0 || topleft_x > (xRes_ - 1)|| rightbottom_x < 0
      || rightbottom_x > (xRes_ - 1)) {
    std::cout << "obstacle_2d.cpp " << __LINE__ << " FATAL: " <<  "Invalid range x" << std::endl; exit(0);
  }
  if (topleft_y < 0 || topleft_y > (yRes_ - 1)|| rightbottom_y < 0 
      || rightbottom_y > (yRes_ - 1)) {
    std::cout << "obstacle_2d.cpp " << __LINE__ << " FATAL: " <<  "Invalid range y" << std::endl; exit(0);
  }
  if (topleft_x > rightbottom_x || topleft_y > rightbottom_y) {
    std::cout << "obstacle_2d.cpp " << __LINE__ << " FATAL: " <<  "Invalid range" << std::endl; exit(0);
  }
  for (int j = topleft_y; j <= rightbottom_y; j++) {
    for (int i = topleft_x; i <= rightbottom_x; i++) {
      int index = i + j*xRes_;
      obstacles[index] = 1;
    }
  }
}

/*void Obstacle2D::InitializeUsingImage(const std::string& img_name) {
  cimg_library::CImg<unsigned char> image(img_name.c_str());
  if(xRes_!= image.width()) { std::cout << "The dimention of the image need to"
      << " match the domain." << std::endl; exit(0);};
  if(yRes_ != image.height()) {std::cout << "The dimention of the image need to"
      << " match the domain." << std::endl; exit(0);};
  int index = 0;
  for (int j = 0; j < yRes_; j++) {
    for (int i = 0; i < xRes_; i++) {
      // The colored positions are obstacles.
      if (image(i,j,0) > 100 || image(i,j,1) > 100 || image(i,j,2) > 100) {
        obstacles[index] = 1;
      }
      index ++;
    }
  }
  image.clear();
}*/

const char Obstacle2D::Access(const int pos_x, const int pos_y) const {
  const int offset_dx = pos_x - static_cast<int>(dx_);
  const int offset_dy = pos_y - static_cast<int>(dy_);
  if (offset_dx >= 0 && offset_dx < xRes_ && offset_dy >= 0 && offset_dy < yRes_) {
    return obstacles[offset_dx + offset_dy * xRes_];
  } else {
    return 0;
  }
}

inline bool Obstacle2D::WithinRange(const int x, const int y) {
  if (x >= 0 && x < xRes_ && y >= 0 && y < yRes_) {
    return true;
  } else {
    return false;
  }
}
