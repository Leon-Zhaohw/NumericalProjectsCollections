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

 Copyright 2013 Theodore Kim
*/
#ifndef OBSTACLEBOX_3D_H
#define OBSTACLEBOX_3D_H

#include "Eigen"
#include <iostream>
#include <vector>

#include "3D/laplacian_basis_set_3d.h"
#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/obstacle_3d.h"
#include "Alg/VEC3.h"
#include "Alg/MATRIX3.h"


class ObstacleBox3D : public OBSTACLE3D
{
public:
  ObstacleBox3D(const VEC3F& center, const VEC3F& lengths, double period, double translationPeriod,
                const float angle, const VEC3F& axis,
                const VEC3I& gridRes, const double dx);
  ObstacleBox3D(const std::string& fname, const VEC3I& gridRes, const double dx);
  ObstacleBox3D();
  ~ObstacleBox3D();

  // setters
  void set_center(const VEC3F& center) { _center = center; };

  // set the x-, y-, and z-dimensions of the box
  void set_lengths(const VEC3F& lengths) { _lengths = lengths; };
  

  void set_halfLengths() { _halfLengths = 0.5 * _lengths; };

  // time step should be in lock-step with FLUID_3D_MIC
  void set_dt(double dt) { _dt = dt; };

  void initialize_rotationMatrix();
  void update_rotationMatrix();
  
  // is the passed in point inside the box?
  bool inside(const VEC3F& point) const override;

  // convert to degrees for OpenGL
  double thetaDegrees() const { return _theta * 360 / (2 * M_PI); } 

  // getters
  VEC3F* get_velocity() { return &(_velocity); }
 
  // a *clockwise* rotation matrix
  MATRIX3* get_rotation() { return &(_rotation); }

  // spin the box around a fixed axis of rotation
  void spin();

  // draw to OpenGL
  void draw() const override;

  // Calculate the normal component of a given velocity field along the 
  // border of the obstacles. This function won't clear the content of
  // normal_component. Instead it will add the addtion value to it.
  void CalculateNormalForce(
      const VECTOR3_FIELD_3D& velocity, char * check, const double dt, const double cell_len, const double force_amp,
      VECTOR3_FIELD_3D* normal_component) const override;
  
  // Do not consider rotation for now....
  void MoveObstacle(const double dt) override {
     _center += dt * _velocity ;
    
    time_acc_ += dt;
   // _velocity[0] = 0.3;// *sin(time_acc_*2);
    ComputeAABB(this->dx_, this->gridRes_);
  }
  void SetVelocity(const VEC3F& velo) override {
    _velocity = velo;
  }
  void ZeroOutVelocityDensityObstalce(const double cell_len,
                                      VECTOR3_FIELD_3D* velocity, FIELD_3D* density) override;
  void RasterToGrid(unsigned char* obstacle_filed) override;
  void OutPutGeometryToPBRT(const std::string &fname) override;
  void OutPutGeometryToPBRT(std::ofstream& out) override;
  
  void GetFromFile(const std::string& fname);
  bool OutsideBoundary(const VEC3F& GridDim) {
    return !BoxIntersect(_center, _lengths, 0.5*GridDim, GridDim);
  }
  bool GetTrigger() override {
    return false;
  }
  VEC3F getVelocityAt(const VEC3F& point) const override;
  
protected:
  void Init();
  void ComputeAABB(const double dx, const VEC3I& gridRes);
  VEC3F _center;
  VEC3F _original_center;
  VEC3F _lengths;
  VEC3F _halfLengths;
  VEC3F _rotationAxis;

  VEC3F _angularVelocity;
  VEC3F _velocity;
  VEC3F _displacement;
  VEC3F _displacementVelocity;
  MATRIX3 _rotation;
  MATRIX3 static_rotation_;
  float static_angle_;
  VEC3F static_axis_;
  
  double _theta;
  double _period;

  double _phi;
  double _translationPeriod;

  double _currentTime;
  double _dt;
  double time_acc_;
};

#endif  // OBSTACLEBOX_3D_H
