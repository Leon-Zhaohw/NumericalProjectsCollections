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

#ifndef OBSTACLE_CYLINDER_3D_H
#define OBSTACLE_CYLINDER_3D_H

#include "3D/obstacle_3d.h"
#include "Alg/VEC3.h"

class ObstacleCylinder3D : public OBSTACLE3D {
public:
  ObstacleCylinder3D(const VEC3F& center, const float radius, const float height, const double period,
                     const bool rotating, const VEC3I& gridRes, const double dx):
  center_(center), radius_(radius), height_(height), _period(period), is_rotating_(rotating),
  OBSTACLE3D(gridRes, dx) { 
    Init();
  }
  
  ObstacleCylinder3D(const std::string& fname, const VEC3I& gridRes, const double dx):OBSTACLE3D(gridRes, dx) {
    GetFromFile(fname);
    Init();
  }
  
  ~ObstacleCylinder3D(){};
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
  bool inside_band(const VEC3F& point) const;
  double compute_band_weight(const VEC3F& point) const;
  
    // convert to degrees for OpenGL
  double thetaDegrees() const {
                                return _theta * 360.0 / (2.0 * M_PI); } 
  // Do not consider rotation for now.
  void MoveObstacle(const double dt) override;
  
  void SetVelocity(const VEC3F& velo) override {
    velocity_ = velo;
  }
  VEC3F ComputeRotateVelocity(const VEC3F& point) const;
  void ZeroOutVelocityDensityObstalce(const double cell_len,
                                      VECTOR3_FIELD_3D* velocity, FIELD_3D* density) override;
  void RasterToGrid(unsigned char* obstacle_filed) override;
  void OutPutGeometryToPBRT(const std::string &fname) override;
  void OutPutGeometryToPBRT(std::ofstream& out) override;
  void GetFromFile(const std::string& fname);
  void GetListOfPositions(const std::string& fname);
  
  void spin(double dt_);
  void update_rotationMatrix();
  bool GetTrigger() override {
    return trigger_ && use_trigger_;
  }
  VEC3F getVelocityAt(const VEC3F& point) const override;
protected:
  void Init();
  void ComputeAABB(const double dx, const VEC3I& gridRes);
  double time_acc;
  
  VEC3F center_;
  VEC3F original_center_;
  float radius_;
  float radius_squared_;
  float radius_band_squared_;
  
  float height_;
  float half_height_;
  float half_height_band_;
  
  VEC3F velocity_;
  // VEC3F debug_point_;
  MATRIX3 rotation_;
  MATRIX3 inv_rotation_;
  
  double _theta;
  double _period;
  VEC3F rotationAxis_;
  bool is_rotating_;
  bool use_sine_velo_;
  std::string reset_pos_name_;
  std::vector<VEC3F> list_center_;
  std::vector<VEC3F> list_velocity_;
  std::vector<std::string> list_axis_;
  std::vector<double> list_period_;
  
  int list_index = 0;
  std::string current_axis__str_;
  int cur_axis_int_;
  bool trigger_;
  bool use_trigger_;
  
  const float band_width_ = 0.07;
};

#endif  // OBSTACLE_CYLINDER_3D_H
