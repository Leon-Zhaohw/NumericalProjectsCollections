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

#ifndef OBSTACLE_CYLINDER_INTERACT_H
#define OBSTACLE_CYLINDER_INTERACT_H

#include "3D/obstacle_3d.h"
#include "Alg/VEC3.h"

class ObstacleCylinderInterac3D : public OBSTACLE3D {
public:
  
  ObstacleCylinderInterac3D(const std::string& fname, const VEC3I& gridRes, 
                            const double dx):OBSTACLE3D(gridRes, dx) {
    GetFromFile(fname);
    Init();
  }
  // Calculate the normal component of a given velocity field along the 
  // border of the obstacles. This function won't clear the content of
  // normal_component. Instead it will add the addtion value to it.
  void CalculateNormalForce(
      const VECTOR3_FIELD_3D& velocity, char * check, const double dt, const double cell_len, const double force_amp,
      VECTOR3_FIELD_3D* normal_component) const override;
       // draw to OpenGL
  void draw() const override;
  // is the passed in point inside the cylinder ?
  bool inside(const VEC3F& point) const override;
    // convert to degrees for OpenGL
  double thetaDegrees() const {
                                return _theta * 360.0 / (2.0 * M_PI); } 
  // Do not consider rotation for now.
  void MoveObstacle(const double dt) override;
  
  void SetVelocity(const VEC3F& velo) override {
    velocity_ += velo;
  }
  VEC3F ComputeRotateVelocity(const VEC3F& point) const;
  void ZeroOutVelocityDensityObstalce(const double cell_len,
                                      VECTOR3_FIELD_3D* velocity, FIELD_3D* density) override;
  void RasterToGrid(unsigned char* obstacle_filed) override;
  void OutPutGeometryToPBRT(const std::string &fname) override;
  void OutPutGeometryToPBRT(std::ofstream& out) override;
  void GetFromFile(const std::string& fname);
  
  void update_rotationMatrix();
  bool GetTrigger() override {
    return false;
  }
  VEC3F getVelocityAt(const VEC3F& point) const override;
  void spin(double dt_);
  void setOmega(const float mouse_dx, const float mouse_dy,
                const float pos_x, const float pos_y, const double dt) override;
protected:
  void Init();
  void ComputeAABB(const double dx, const VEC3I& gridRes);
  double time_acc;
  
  VEC3F center_;
  VEC3F original_center_;
  float radius_;
  float radius_squared_;
  float height_;
  float half_height_;
  
  VEC3F velocity_;
  // VEC3F debug_point_;
  MATRIX3 rotation_;
  MATRIX3 inv_rotation_;
  
  double _theta;
  // Angular.
  double omega_;
  VEC3F rotationAxis_;
  bool is_rotating_;
  
  int list_index = 0;
  std::string current_axis__str_;
  int cur_axis_int_;
  bool trigger_;
  bool use_trigger_;
};

#endif  // OBSTACLE_CYLINDER_INTERACT_H
