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

#include "2D/two_neumann_basis_2D.h"
#include "2D/VFIELD2D.h"
#include "util/util.h"

void TwoNeumannBasis2D::Normalize() {
  int num_zeros = 0;
  if (k1_ == 0) num_zeros ++;
  if (k2_ == 0) num_zeros ++;
  
  if (num_zeros == 0) {
    double S_ = A_*A_ + B_*B_;
    if (std::abs(S_) > 1e-15) {
      invNorm_ = 1.0 / (pow(S_*PI_SQUARE*0.25, 0.5));
    }
  }
  
  if (num_zeros == 1) {
    double S_ = 0;
    if (k1_ == 0) {
      S_ = A_*A_;
    } else if (k2_ == 0) {
      S_ = A_*A_;
    }
    
    if (std::abs(S_) > 1e-15) {
      invNorm_ = 1.0 / (pow(S_*PI_SQUARE*0.5, 0.5));
    }
  }
  
  if (num_zeros == 2) {
    if (std::abs(A_) > 1e-15) {
      invNorm_ = 1.0 / (A_*sqrt(PI_SQUARE));
    }
  }
}

double TwoNeumannBasis2D::ComputeBasisAt(const int xpos, const int ypos,
                                         const int mode, const VEC2& dxy) const {
  // vx...
  if (mode == 0) {
    return invNorm_*A_*cos((xpos+0.5)*k1_*dxy[0])*cos((ypos+0.5)*k2_*dxy[1]);
  } else {  // vy...
    return invNorm_*B_*sin((xpos+0.5)*k1_*dxy[0])*sin((ypos+0.5)*k2_*dxy[1]);
  }
}
                        
void TwoNeumannBasis2D::DiscretizeAdd(const double coef, VFIELD2D* vfield) {
  if (coef == 0) {
    return;
  }
  const int xRes = vfield->getxRes();
  const int yRes = vfield->getyRes();
  
  const double dx_ = M_PI / static_cast<double>(xRes);
  const double dy_ = M_PI / static_cast<double>(yRes);
  int index = 0;
  Real vx,vy;
  for (int j = 0; j < yRes; j++) {
    for (int i = 0; i < xRes; i++) {
      vx = 0; vy = 0;
      if (A_ != 0) {
        vx = invNorm_*A_*cos((i+0.5)*k1_*dx_)*cos((j+0.5)*k2_*dy_);
      }
      if (k2_ != 0 && k1_ != 0 && B_ != 0) {
        vy = invNorm_*B_*sin((i+0.5)*k1_*dx_)*sin((j+0.5)*k2_*dy_);
      }
      (*vfield)[index][0] += coef*vx;
      (*vfield)[index][1] += coef*vy;
      index ++;
    }
  }
}
