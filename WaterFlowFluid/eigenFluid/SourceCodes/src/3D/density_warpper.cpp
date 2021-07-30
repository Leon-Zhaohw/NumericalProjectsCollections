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

#include "Alg/VEC3.h"
#include "3D/FIELD_3D.h"
#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/density_warpper.h"



void DensityWarpper::AddBuoyancyField(VECTOR3_FIELD_3D* force) {
  if (! initialized_) {std::cout << "density_warpper.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  
  for (uint i = 0; i < totalsize_; i++) {
    (*force)[i][0] += buoyancy_*density_[i]*buoyancy_dir_[0];
    (*force)[i][1] += buoyancy_*density_[i]*buoyancy_dir_[1];
    (*force)[i][2] += buoyancy_*density_[i]*buoyancy_dir_[2];
  }
}

void  DensityWarpper::AddSmokeTestCase() {
  if (! initialized_) {std::cout << "density_warpper.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  VEC3F pointTransformed = source_pos_flt_ - cell_center_;
  pointTransformed += cell_center_;
  
  souce_pos_[0] = static_cast<int>(pointTransformed[0]*maxRes_);
  souce_pos_[1] = static_cast<int>(pointTransformed[1]*maxRes_);
  souce_pos_[2] = static_cast<int>(pointTransformed[2]*maxRes_);
  AddDensity(souce_pos_[0], souce_pos_[1], souce_pos_[2],
             souce_size_[0], souce_size_[1], souce_size_[2], &density_);
}

void DensityWarpper::AddDensity(const int xpos, const int ypos, const int zpos, 
                const int length, const int width, const int height,
                FIELD_3D* field) {
  int xbegin = xpos - length / 2;
  int xend = xpos + length / 2;
  int ybegin = ypos - width / 2;
  int yend = ypos + width / 2;
  int zbegin = zpos - height / 2;
  int zend = zpos + height / 2;
  
  // clamp.
  xbegin = (xbegin < 0) ? 0 : xbegin;
  xend = (xend < 0) ? 0 : xend;
  ybegin = (ybegin < 0) ? 0 : ybegin;
  yend = (yend < 0) ? 0 : yend;
  zbegin = (zbegin < 0 ) ? 0 : zbegin;
  zend = (zend < 0 ) ? 0 : zend;
  
  xbegin = (xbegin > Res_[0] - 1) ? Res_[0] - 1 : xbegin;
  xend = (xend > Res_[0] - 1) ? Res_[0] - 1 : xend;
  ybegin = (ybegin > Res_[1] - 1) ? Res_[1] - 1 : ybegin;
  yend = (yend > Res_[1] - 1) ? Res_[1] - 1 : yend;
  zbegin = (zbegin > Res_[2] - 1) ? Res_[2] - 1 : zbegin;
  zend = (zend > Res_[2] - 1) ? Res_[2] - 1 : zend;
  const int Slabsize_ = Res_[0]*Res_[1];
  uint idx_begin = xbegin + ybegin*Res_[0] + zbegin*Slabsize_;
  
  for (int k = zbegin; k <= zend; k++) {
    for (int j = ybegin; j <= yend; j++) {
      for (int i = xbegin; i <= xend; i++) {
        int ind = idx_begin + i + j*Res_[0] + k*Slabsize_;
        (*field)(i,j,k) = addedDensity_;
      }
    }
  }
}
