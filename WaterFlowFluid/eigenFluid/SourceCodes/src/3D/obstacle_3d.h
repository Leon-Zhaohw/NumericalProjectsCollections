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

#ifndef OBSTACLE3D_3D_H
#define OBSTACLE3D_3D_H

#include "Eigen"
#include <iostream>
#include <fstream>

#include "3D/VECTOR3_FIELD_3D.h"
#include "Alg/VEC3.h"
#include "Alg/MATRIX3.h"
#include "FIELD_3D.h"

using namespace std;


struct ObstacleParams3D{
  ObstacleParams3D() {
    use_obstacles = false;
    obstacle_force_scale = 20.;
  }
  bool use_obstacles;
  double obstacle_force_scale;
  bool handle_obstacle_implicit;
  bool move_obstacle;
  std::string obstacle_type;
  std::string obstacle_file;
  std::string obstacle_list;
};

struct ObstacleIdx3D {
  ObstacleIdx3D(const int xpos, const int ypos,
                const int zpos, const int mode_):idx(VEC3I(xpos, ypos, zpos)), mode(mode_){
  }
  // x, y, z...
  VEC3I idx;
  int mode;
};

class OBSTACLE3D  
{
public:
  OBSTACLE3D(const VEC3I& gridRes, const double dx):gridRes_(gridRes),dx_(dx) {
    // Default color.
    color_ = VEC3F(0.1, 0.2, 0.8);
  };
  virtual ~OBSTACLE3D() {};

  bool inside(float x, float y, float z) {return inside(VEC3F(x,y,z));}
  virtual void draw() const {
    std::cout << "obstacle_3d.h " << __LINE__ << " FATAL: " <<  "Base class called" << std::endl; exit(0);
  }
  virtual bool inside(const VEC3F& point) const {
    std::cout << "obstacle_3d.h " << __LINE__ << " FATAL: " <<  "Base class called." << std::endl; exit(0);
  }
  // virtual Real distance(const VEC3F& point) = 0;

  // Eigen::Vector3f& translation()     { return Translation_; };
  // Eigen::Matrix3f& rotation()      { return Rotation_; };
  // Eigen::Vector3f& velocity()        { return V_; };

  // virtual void scale(const Real alpha) = 0;
  // virtual void boundingBox(VEC3F& mins, VEC3F& maxs) const = 0;
  virtual void CalculateNormalForce(
      const VECTOR3_FIELD_3D& velocity, char * check, const double dt, const double cell_len, const double force_amp,
      VECTOR3_FIELD_3D* normal_component) const {
        std::cout << "obstacle_3d.h " << __LINE__ << " FATAL: " <<  "Base class called." << std::endl; exit(0);
      }
  
  virtual void MoveObstacle(const double dt) {
    std::cout << "obstacle_3d.h " << __LINE__ << " FATAL: " <<  "Base class called." << std::endl; exit(0);
  }
  
  virtual void SetVelocity(const VEC3F& velo) {
    std::cout << "obstacle_3d.h " << __LINE__ << " FATAL: " <<  "Base class called." << std::endl; exit(0);
  }
  
  virtual
  void ZeroOutVelocityDensityObstalce(const double cell_len,
                                      VECTOR3_FIELD_3D* velocity, FIELD_3D* density) {
    std::cout << "obstacle_3d.h " << __LINE__ << " FATAL: " <<  "Base class called." << std::endl; exit(0);
  }
  
  virtual void RasterToGrid(unsigned char* obstacle_filed) {
    std::cout << "obstacle_3d.h " << __LINE__ << " FATAL: " <<  "Base method called." << std::endl; exit(0);
  }
  
  virtual void OutPutGeometryToPBRT(const std::string &fname) {
    std::cout << "obstacle_3d.h " << __LINE__ << " FATAL: " <<  "Base class called." << std::endl; exit(0);
  }
  
  virtual void OutPutGeometryToPBRT(std::ofstream& out) {
    std::cout << "obstacle_3d.h " << __LINE__ << " FATAL: " <<  "Base class called." << std::endl; exit(0);
  }
  void set_color(const VEC3F& color) {color_ = color;};
  virtual bool GetTrigger() {
    std::cout << "Based called" << std::endl;
    return false;
  }
  virtual VEC3F getVelocityAt(const VEC3F& point) const {
    std::cout << "obstacle_3d.h " << __LINE__ << " FATAL: " <<  "Base called." << std::endl; exit(0);
  }
  virtual void setOmega(const float mouse_dx, const float mouse_dy,
                        const float pos_x, const float pos_y, const double dt) {
    return;
  }
protected:
  const VEC3I gridRes_;
  // Axis aligned boundingBox
  VEC3I AABBStart_;
  VEC3I AABBEnd_;
  const double dx_;
  VEC3F color_;
  
};


#endif  // OBSTACLE3D_3D_H
