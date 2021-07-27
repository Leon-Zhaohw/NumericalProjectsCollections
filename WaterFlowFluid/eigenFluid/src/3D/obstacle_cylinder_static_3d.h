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

#ifndef OBSTACLE_CYLINDER_STATIC_3D_H
#define OBSTACLE_CYLINDER_STATIC_3D_H

#include "3D/obstacle_3d.h"
#include "Alg/VEC3.h"

class ObstacleCylinderStatic3D : public OBSTACLE3D {
public:
  ObstacleCylinderStatic3D(const VEC3F& center, const float radius, const float height, const float angle, 
    const VEC3F& axis, const VEC3I& gridRes, const double dx):
  center_(center), radius_(radius), height_(height), angle_degree_(angle), rotationAxis_(axis),
  OBSTACLE3D(gridRes, dx) {
    Init();
  }
  ObstacleCylinderStatic3D(const std::string& fname, const VEC3I& gridRes, const double dx):
  OBSTACLE3D(gridRes, dx){
    GetFromFile(fname);
    Init();
  }
  ~ObstacleCylinderStatic3D(){};
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
  void MoveObstacle(const double dt) override {
  }
  void SetVelocity(const VEC3F& velo) override {
    velocity_ = velo;
  }
  void ZeroOutVelocityDensityObstalce(const double cell_len,
                                      VECTOR3_FIELD_3D* velocity, FIELD_3D* density) override;
  void RasterToGrid(unsigned char* obstacle_filed) override;
  void OutPutGeometryToPBRT(const std::string &fname) override;
  void OutPutGeometryToPBRT(std::ofstream& out) override;
  
  void GetFromFile(const std::string& fname);
  void ComputeAABB(const double dx, const VEC3I& gridRes);
  bool GetTrigger() override {
    return false;
  }
  VEC3F getVelocityAt(const VEC3F& point) const override;
protected:
  void Init();
  double time_acc;
  
  VEC3F center_;
  VEC3F original_center_;
  float radius_;
  float radius_squared_;
  float height_;
  float half_height_;
  float angle_degree_;
  
  VEC3F velocity_;
  MATRIX3 rotation_;
  MATRIX3 inv_rotation_;
  
  double _theta;
  VEC3F rotationAxis_;
};

#endif  // OBSTACLE_CYLINDER_STATIC_3D_H

