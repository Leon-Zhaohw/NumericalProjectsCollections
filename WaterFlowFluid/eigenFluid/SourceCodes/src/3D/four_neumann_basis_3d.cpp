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

#include <math.h>

#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/four_neumann_basis_3d.h"
#include "util/util.h"

void FourNeumannBasis3D::DiscretizeAdd(const double coef, VECTOR3_FIELD_3D* vfield) {
  if (coef == 0) {
    return;
  }
  const int xRes = vfield->xRes();
  const int yRes = vfield->yRes();
  const int zRes = vfield->zRes();
  
  const double dx_ = M_PI / static_cast<double>(xRes);
  const double dy_ = M_PI / static_cast<double>(yRes);
  const double dz_ = M_PI / static_cast<double>(zRes);
  int index = 0;
  Real vx,vy,vz;
  for (int k = 0; k < zRes; k++) {
    for (int j = 0; j < yRes; j++) {
      for (int i = 0; i < xRes; i++) {
        vx = 0; vy = 0; vz = 0;
        if (A_ != 0 && k3_ != 0) {
          vx = invNorm_*A_*cos((i+0.5)*k1_*dx_)*cos((j+0.5)*k2_*dy_)*sin((k+0.5)*k3_*dz_);
        }
        if (k2_ != 0 && k1_ != 0 && k3_ != 0 && B_ != 0) {
          vy = invNorm_*B_*sin((i+0.5)*k1_*dx_)*sin((j+0.5)*k2_*dy_)*sin((k+0.5)*k3_*dz_);
        }
        if (k1_ != 0 && C_ != 0) {
          vz = invNorm_*C_*sin((i+0.5)*k1_*dx_)*cos((j+0.5)*k2_*dy_)*cos((k+0.5)*k3_*dz_);
        }
        (*vfield)[index][0] += coef*vx;
        (*vfield)[index][1] += coef*vy;
        (*vfield)[index][2] += coef*vz;
        index ++;
      }
    }
  }
}

double FourNeumannBasis3D::ComputeBasisAt(const int xpos, const int ypos,
                        const int zpos, const int mode,
                        const Eigen::Vector3d& dxyz) const {
   // vx...
  if (mode == 0) {
    return invNorm_*A_*cos((xpos+0.5)*k1_*dxyz[0])*cos((ypos+0.5)*k2_*dxyz[1])
                      *sin((zpos+0.5)*k3_*dxyz[2]);
                      
  } else if (mode == 1) {   // vy...
    return invNorm_*B_*sin((xpos+0.5)*k1_*dxyz[0])*sin((ypos+0.5)*k2_*dxyz[1])
                      *sin((zpos+0.5)*k3_*dxyz[1]);
  } else {    // vz...
      return invNorm_*C_*sin((xpos+0.5)*k1_*dxyz[0])*cos((ypos+0.5)*k2_*dxyz[1])
                        *cos((zpos+0.5)*k3_*dxyz[2]);
  }
  
  return 0;
}

void FourNeumannBasis3D::Normalize() {
  int num_zeros = 0;
  if (k1_ == 0) num_zeros ++;
  if (k2_ == 0) num_zeros ++;
  if (k3_ == 0) num_zeros ++;
  if (num_zeros == 0) {
    double S_ = A_*A_ + B_*B_ + C_*C_;
    if (std::abs(S_) > 1e-15) {
      invNorm_ = 1.0 / (sqrt(S_*PI_CUBE*0.125));
    }
  }
  if (num_zeros == 1) { 
    double S_ = 0;
    if (k1_ == 0) {  // vy, vz == 0
      S_ = A_*A_;
    } else if (k2_ == 0) { // vy == 0
      S_ = A_*A_ + C_*C_;
    } else if (k3_ == 0) {  // vx, vy == 0
      S_ = C_*C_;
    }
    if (std::abs(S_) > 1e-15) {
      invNorm_ = 1.0 / (sqrt(S_*PI_CUBE*0.25));
    }
  }
  if (num_zeros == 2) {
    if (k1_ == 0 && k2_ == 0) { // vx != 0
      double S_ = A_*A_;
      if (std::abs(S_) > 1e-15) {
        invNorm_ = 1.0 / (sqrt(S_*PI_CUBE*0.5));
      }
    }
    if (k1_ == 0 && k3_ == 0) { // vx = 0, vy = 0, vz = 0
      invNorm_ = 0;
    }
    if (k2_ == 0 && k3_ == 0) { // vz != 0
      double S_ = C_*C_;
      if (std::abs(S_) > 1e-15) {
        invNorm_ = 1.0 / (sqrt(S_*PI_CUBE*0.5));
      }
    }
  }
  if (num_zeros == 3) {  // zero mode.
    invNorm_ = 0;
  }
}
