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

#ifndef DENSITY_WARPPER_H
#define DENSITY_WARPPER_H

#include "Eigen"
#include "Alg/VEC3.h"
#include "3D/FIELD_3D.h"
#include "3D/VECTOR3_FIELD_3D.h"

class DensityWarpper {
  
public:
  DensityWarpper(const VEC3I& res, const VEC3F buoyancy_dir, const float buoyancy,
    const float addedDensity):Res_(res),
  buoyancy_dir_(buoyancy_dir), buoyancy_(buoyancy) , addedDensity_(addedDensity)
  {
    density_ = FIELD_3D(Res_[0], Res_[1], Res_[2]);
    density_old_ = FIELD_3D(Res_[0], Res_[1], Res_[2]);
    initialized_ = false;
    totalsize_ = Res_[0]*Res_[1]*Res_[2];
    maxRes_ = std::max(std::max(Res_[0], Res_[1]), Res_[2]);
  };
  ~DensityWarpper() {};
  void AddBuoyancyField(VECTOR3_FIELD_3D* force);
  void InitialzeSourceSomke(const VEC3F&  pos, const VEC3F& ssize, const int MaxRes) {
    source_pos_flt_[0] = pos[0];
    source_pos_flt_[1] = pos[1];
    source_pos_flt_[2] = pos[2];
    souce_size_[0] = static_cast<int>(ssize[0]*MaxRes);
    souce_size_[1] = static_cast<int>(ssize[1]*MaxRes);
    souce_size_[2] = static_cast<int>(ssize[2]*MaxRes);
    initialized_ = true;
  }
  void AddSmokeTestCase();
  void AddDensity(const int xpos, const int ypos, const int zpos, 
                const int length, const int width, const int height,
                FIELD_3D* field);
  void SwapPointer() {
    density_.swapPointers(density_old_);
  }
  void AttenuateSmoke(const float factor) {
    for (uint i = 0; i < totalsize_; i++) {
      density_[i] *= factor;
      if (density_[i] < 0.003) {
        density_[i] = 0;
      }
    }
  }
  FIELD_3D density_;
  FIELD_3D density_old_;
private:
  int maxRes_;
  bool initialized_;
  const VEC3F buoyancy_dir_;
  const VEC3I Res_;
  uint totalsize_;
 
  VEC3I souce_pos_;
  VEC3I souce_size_;
  VEC3F cell_center_;
  VEC3F source_pos_flt_;
  const float buoyancy_;
  const float addedDensity_;
  
};

#endif  // DENSITY_WARPPER_H
