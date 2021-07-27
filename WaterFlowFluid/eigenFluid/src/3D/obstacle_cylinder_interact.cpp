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

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <float.h>

#include <fstream>
#include "3D/obstacle_cylinder_interact.h"
#include "util/util.h"
#include "util/stringprintf.h"

void ObstacleCylinderInterac3D::Init() {
  // rotation_ = MATRIX3::I();
  rotationAxis_ = VEC3I(0,0,1);
  current_axis__str_ = "z";
  cur_axis_int_ = 2;
  rotation_ =  MATRIX3::rotation(VEC3(0,1,0), M_PI/2);
  inv_rotation_ = MATRIX3::rotation(VEC3(0,1,0), -M_PI/2);
  radius_squared_ = radius_*radius_;
  half_height_ = height_*0.5;
  time_acc = 0;
  ComputeAABB(this->dx_, this->gridRes_);
  trigger_ = true;
  // debug_point_ = VEC3F(0.2,0.34,0.26);
}

// Calculate the normal component of a given velocity field along the 
// border of the obstacles. This function won't clear the content of
// normal_component. Instead it will add the addtion value to it.
void ObstacleCylinderInterac3D::CalculateNormalForce(
    const VECTOR3_FIELD_3D& velocity, char * check, const double dt, const double cell_len, const double force_amp,
    VECTOR3_FIELD_3D* normal_component) const {
  int xRes = velocity.xRes();
  int yRes = velocity.yRes();
  int zRes = velocity.zRes();
  int slab = velocity.slabSize();
  
  if (xRes != normal_component->xRes()) {std::cout << "obstacle_cylinder_interact.cpp " << __LINE__ << " FATAL: "  << "The dimention" << " must be equal" << std::endl; exit(0);}
  if (yRes != normal_component->yRes()) {std::cout << "obstacle_cylinder_interact.cpp " << __LINE__ << " FATAL: "  << "The dimention" << " must be equal" << std::endl; exit(0);}
  if (zRes != normal_component->zRes()) {std::cout << "obstacle_cylinder_interact.cpp " << __LINE__ << " FATAL: "  << "The dimention" << " must be equal" << std::endl; exit(0);}
  
  const double amp = force_amp / dt;
#pragma omp parallel for
  for (int k = AABBStart_[2]; k < AABBEnd_[2]; k++) {
    for (int j = AABBStart_[1]; j < AABBEnd_[1]; j++) {
      for (int i = AABBStart_[0]; i < AABBEnd_[0]; i++) {
        const int index = i + j*xRes + k*slab;
        
        VEC3F pt((static_cast<float>(i) + 0.5)*cell_len,
                (static_cast<float>(j) + 0.5)*cell_len,
                (static_cast<float>(k) + 0.5)*cell_len);
        // skip the cells that have already been checked by other obstacles.
        if (check[index]) continue;
        // Cell is inside the obstacle.
        if (inside(pt)) {
          VEC3F angular_velo(0,0,0);
          if (is_rotating_) {
            angular_velo = ComputeRotateVelocity(pt);
          }
          // Superimpose the velocity.
          angular_velo += velocity_;
          (*normal_component)[index][0] -= amp*(velocity[index][0] - angular_velo[0]);
          (*normal_component)[index][1] -= amp*(velocity[index][1] - angular_velo[1]);
          (*normal_component)[index][2] -= amp*(velocity[index][2] - angular_velo[2]);
          check[index] = 1;
          continue;
        }
        int  left = 0, right = 0, top = 0, bottom = 0, front = 0, back = 0;
        // left
        pt[0] -= cell_len;
        if (inside(pt)) left ++;
        pt[0] += cell_len;
        // right
        pt[0] += cell_len;
        if (inside(pt)) right ++;
        pt[0] -= cell_len;
        // top
        pt[1] += cell_len;
        if (inside(pt)) top++;
        pt[1] -= cell_len;
        // bottom
        pt[1] -= cell_len;
        if (inside(pt)) bottom++;
        pt[1] += cell_len;
        // front
        pt[2] += cell_len;
        if (inside(pt)) front++;
        pt[2] -= cell_len;
        // back
        pt[2] -= cell_len;
        if (inside(pt)) back++;
        pt[2] += cell_len;
        int total = top + bottom + left + right + front + back;
        // No obstacle nearby.
        if (total == 0) continue;
        else check[index] = 1;
        VEC3F angular_velo(0,0,0);
        if (is_rotating_) {
          angular_velo = ComputeRotateVelocity(pt);
        }
        angular_velo += velocity_;
        
        // Left side. x component.
        if (left != 0) {
          (*normal_component)[index][0] -= amp*(velocity[index][0] - angular_velo[0]);
        }
        // right side. x component.
        if (right != 0) {
          (*normal_component)[index][0] -= amp*(velocity[index][0] - angular_velo[0]);
        }
        // top, y component.
        if (top != 0) {
          (*normal_component)[index][1] -= amp*(velocity[index][1] - angular_velo[1]);
        }
        // bottom, y component.
        if (bottom != 0) {
          (*normal_component)[index][1] -= amp*(velocity[index][1] - angular_velo[1]);
        }
        // front, z component.
        if (front != 0) {
          (*normal_component)[index][2] -= amp*(velocity[index][2] - angular_velo[2]);
        }
        // back, z component.
        if (back != 0) {
          (*normal_component)[index][2] -= amp*(velocity[index][2] - angular_velo[2]);
        }
      }
    }
  }
}

VEC3F ObstacleCylinderInterac3D::getVelocityAt(const VEC3F& point) const {
  if (this->inside(point)) {
    VEC3F angular_velo(0,0,0);
    if (is_rotating_) {
      angular_velo = ComputeRotateVelocity(point);
    }
    // Superimpose the velocity.
    angular_velo += velocity_;
    return angular_velo;
  } else {
    return VEC3F(0,0,0);
  }
}

void ObstacleCylinderInterac3D::ZeroOutVelocityDensityObstalce(const double cell_len,
                                      VECTOR3_FIELD_3D* velocity, FIELD_3D* density) {
  const int xRes = velocity->xRes();
  const int yRes = velocity->yRes();
  const int zRes = velocity->zRes();
  const int Slab = velocity->slabSize();
#pragma omp parallel for
  for (int k = AABBStart_[2]; k < AABBEnd_[2]; k++) {
    for (int j = AABBStart_[1]; j < AABBEnd_[1]; j++) {
      for (int i = AABBStart_[0]; i < AABBEnd_[0]; i++) {
        VEC3F pt((static_cast<float>(i) + 0.5)*cell_len,
                 (static_cast<float>(j) + 0.5)*cell_len,
                 (static_cast<float>(k) + 0.5)*cell_len);
        if (this->inside(pt)) {
          const int idx = i + j*xRes + k*Slab;
          
          VEC3F angular_velo(0,0,0);
          if (is_rotating_) {
            angular_velo = ComputeRotateVelocity(pt);
          }
          // Superimpose the velocity.
          angular_velo += velocity_;
          (*velocity)[idx] = angular_velo;
          (*density)[idx] = 0;
        }
      }
    }
  }
}

void ObstacleCylinderInterac3D::RasterToGrid(unsigned char* obstacle_filed) {
  const int xRes = gridRes_[0];
  const int yRes = gridRes_[1];
  const int zRes = gridRes_[2];
  const int Slab = xRes*yRes;
  #pragma omp parallel for
  for (int k = AABBStart_[2]; k < AABBEnd_[2]; k++) {
    for (int j = AABBStart_[1]; j < AABBEnd_[1]; j++) {
      for (int i = AABBStart_[0]; i < AABBEnd_[0]; i++) {
        VEC3F pt((static_cast<float>(i) + 0.5)*dx_,
                 (static_cast<float>(j) + 0.5)*dx_,
                 (static_cast<float>(k) + 0.5)*dx_);
        if (this->inside(pt)) {
          const int idx = i + j*xRes + k*Slab;
          obstacle_filed[idx] = 1;
        }
      }
    }
  }
}

void ObstacleCylinderInterac3D::MoveObstacle(const double dt)  {
    center_ += dt * velocity_ ;
    time_acc += dt;

    if (is_rotating_) {
      spin(dt);
    }
    ComputeAABB(this->dx_, this->gridRes_);
    // if move out side..
    
    if (center_[0] < 0 || center_[0] > 1.0) {
      velocity_ = 0;
    }
}

void ObstacleCylinderInterac3D::spin(double dt_) {
  double d_theta = omega_ * dt_;
  _theta += d_theta;
  // printf("Theta is updated to: %f\n", _theta);
  // work modulo 2pi
  if (std::abs(_theta) >= 2 * M_PI) { _theta = 0.0; };
  this->update_rotationMatrix();
}

VEC3F ObstacleCylinderInterac3D::ComputeRotateVelocity(const VEC3F& point) const {
  // The cylinder is rotating along z axis in GLOBAL coordinate system.
  // angular velocity.
  float omega = omega_;
  VEC3F pointTransformed = point - center_;
  // pointTransformed = rotation_ * pointTransformed;
  if (cur_axis_int_ == 2) {  // z 
    const float xd = pointTransformed[0];
    const float yd = pointTransformed[1];
    const float r = sqrt(xd*xd + yd*yd);
  
    if (std::abs(r) < 1e-7) {
      return VEC3F(0,0,0);
    }
    float v_scalar = r*omega;
    // the velocity is in GLOBAL coordinate system.
    return VEC3F(-v_scalar*yd/r, v_scalar*xd/r, 0);
  } else if (cur_axis_int_ == 1) { // y
    const float xd = pointTransformed[0];
    const float zd = pointTransformed[2];
    const float r = sqrt(xd*xd + zd*zd);
  
    if (std::abs(r) < 1e-7) {
      return VEC3F(0,0,0);
    }
    float v_scalar = r*omega;
    // the velocity is in GLOBAL coordinate system.
    return VEC3F(v_scalar*zd/r, 0, -v_scalar*xd/r);
    
  } else {  // x
      const float yd = pointTransformed[1];
      const float zd = pointTransformed[2];
      const float r = sqrt(yd*yd + zd*zd);
      if (std::abs(r) < 1e-7) {
        return VEC3F(0,0,0);
      }
      float v_scalar = r*omega;
      // the velocity is in GLOBAL coordinate system.
      return VEC3F(0, -v_scalar*zd/r, v_scalar*yd/r);
  }
}

void ObstacleCylinderInterac3D::update_rotationMatrix() {
  // we must rotate *clockwise*, hence the minus sign
  rotation_ = MATRIX3::rotation(VEC3(0,1,0), M_PI/2)*MATRIX3::rotation(rotationAxis_, -1 * _theta);
  // inv_rotation_ is inverse of rotation_.
  inv_rotation_ = MATRIX3::rotation(rotationAxis_, _theta)*MATRIX3::rotation(VEC3(0,1,0), -M_PI/2);
}

#define EXPECT_STRING(str) \
in >> temp; \
if (temp != str) { \
  std::cout << "obstacle_cylinder_interact.cpp " << __LINE__ << " FATAL: " <<  "Error: " << temp << std::endl; exit(0); \
}\
temp.clear();
void ObstacleCylinderInterac3D::GetFromFile(const std::string& fname) {
  std::ifstream in(fname);
  if (!in.is_open()) {
    std::cout << "obstacle_cylinder_interact.cpp " << __LINE__ << " FATAL: " <<  "Cannot open: " << fname << std::endl; exit(0);
  }
  std::string temp;
  EXPECT_STRING("cylinder_interact")
  EXPECT_STRING("center")
  in >> center_[0]; in >> center_[1]; in >> center_[2];
  EXPECT_STRING("height")
  in >> height_;
  EXPECT_STRING("radius")
  in >> radius_;
  EXPECT_STRING("rotating")
  in >> is_rotating_;
  EXPECT_STRING("color")
  in >> color_[0]; in >> color_[1]; in >> color_[2];
  EXPECT_STRING("velocity")
  in >> velocity_[0]; in >> velocity_[1]; in >> velocity_[2];
  in.close();
}

#undef EXPECT_STRING

void ObstacleCylinderInterac3D::OutPutGeometryToPBRT(std::ofstream& outfile) {
  if (!outfile.is_open()) {
    std::cout << "obstacle_cylinder_interact.cpp " << __LINE__ << " FATAL: " <<  "Cannot open file." << std::endl; exit(0);
  }
  outfile << "AttributeBegin \n TransformBegin\n";
  outfile << "Translate " << center_[0] << " " << center_[1] << " " << center_[2] << "\n";
  outfile << "Rotate " << thetaDegrees() << " " << rotationAxis_[0] << " " << rotationAxis_[1]
      << " " << rotationAxis_[2] << "\n";
  // The axis is along x.
  outfile << "Rotate 90 0 1 0\n";
  // Texture and color of the obstacle.
    outfile << "Material \"matte\" \"color Kd\" [ " << color_[0] << " " <<  color_[1]
      << " " <<  color_[2] << " ]\n";
      
  outfile << "TransformBegin\n";
  outfile << "Translate 0 0 " << -half_height_;
  outfile << "Shape \"cylinder\" \"float radius\" " << radius_  << " \"float zmin\" 0 \"float zmax\" "
      << height_ << "\n";
  outfile << "Shape \"disk\" \"float radius\" " << radius_ << "\n";
  outfile << "Translate 0 0 " << height_ << "\n";
  outfile << "Shape \"disk\" \"float radius\" " << radius_ << "\n";
  outfile << "TransformEnd \n TransformEnd \n AttributeEnd\n";
}

void ObstacleCylinderInterac3D::OutPutGeometryToPBRT(const std::string &fname) {
  const std::string m_fname = StringPrintf("%scylinder.pbrt", fname.c_str());
  std::ofstream outfile(m_fname.c_str());
  OutPutGeometryToPBRT(outfile);
  outfile.close();
}

// draw to OpenGL    
void ObstacleCylinderInterac3D::draw() const {
   glPushMatrix();
    
    glTranslatef(center_[0], center_[1], center_[2]);
    // OpenGL uses degrees, so we need to convert
    glRotatef(thetaDegrees(), rotationAxis_[0], rotationAxis_[1], rotationAxis_[2]);
    glRotatef(90, 0,1,0);
    glColor4f(0.1, 0.7, 0.8, 0.6);
    GLUquadric *quad;
    quad = gluNewQuadric();
    glTranslatef(0,0, -half_height_);
    gluCylinder(quad, radius_, radius_, height_, 20,20);
    gluDisk(quad, 0, radius_, 20,20);
    glTranslatef(0,0,height_);
    gluDisk(quad, 0, radius_, 20,20);
  glPopMatrix();
}

void ObstacleCylinderInterac3D::ComputeAABB(const double dx, const VEC3I& gridRes) {
    // Get corners in local frame.
  std::vector<VEC3F> corners;
  corners.push_back(VEC3F(-radius_, -radius_, -half_height_));
  corners.push_back(VEC3F(-radius_, -radius_, half_height_));
  corners.push_back(VEC3F(radius_, -radius_, -half_height_));
  corners.push_back(VEC3F(radius_, -radius_, half_height_));
  corners.push_back(VEC3F(-radius_, radius_, -half_height_));
  corners.push_back(VEC3F(-radius_, radius_, half_height_));
  corners.push_back(VEC3F(radius_, radius_, -half_height_));
  corners.push_back(VEC3F(radius_, radius_, half_height_));
  // Transform into global frame.
  for (int i = 0; i < corners.size(); i++) {
    VEC3F transfomed = inv_rotation_ * corners[i];
    transfomed += center_;
    corners[i] = transfomed;
  }
  // Get min max for x,y,z.
  VEC3F minC(FLT_MAX, FLT_MAX, FLT_MAX);
  VEC3F maxC(FLT_MIN, FLT_MIN, FLT_MIN);
  for (int i = 0; i < corners.size(); i++) {
    for (int id = 0; id < 3; id++) {
      if (corners[i][id] < minC[id]) {
        minC[id] = corners[i][id];
      }
      if (corners[i][id] > maxC[id]) {
        maxC[id] = corners[i][id];
      }
    }
  }
    // Get grid indices on grid.
  AABBStart_[0] = static_cast<int>(minC[0]/dx);
  AABBStart_[1] = static_cast<int>(minC[1]/dx);
  AABBStart_[2] = static_cast<int>(minC[2]/dx);
  AABBEnd_[0] = static_cast<int>(maxC[0]/dx);
  AABBEnd_[1] = static_cast<int>(maxC[1]/dx);
  AABBEnd_[2] = static_cast<int>(maxC[2]/dx);
  // Clamp.
  for (int i = 0; i < 3; i++) {
    if (AABBStart_[i] < 0) AABBStart_[i] = 0;
    if (AABBEnd_[i] < 0) AABBEnd_[i] = 0;
    if (AABBStart_[i] >= gridRes[i]) AABBStart_[i] = gridRes[i] - 1;
    if (AABBEnd_[i] >= gridRes[i]) AABBEnd_[i] = gridRes[i] - 1;
  }
}

// pass in the delta of mouse on screen, screen size is normalized to [-1, 1].
void ObstacleCylinderInterac3D::setOmega(const float mouse_dx, const float mouse_dy,
                                         const float pos_x, const float pos_y, const double dt) {
  // Convert to world coord. [0, 1]
  float cur_pos_x = (pos_x + 1.0)*0.5*static_cast<float>(gridRes_[0])*dx_;
  float cur_pos_y = (pos_y + 1.0)*0.5*static_cast<float>(gridRes_[1])*dx_;
  // std::cout <<  "dx dy: " << mouse_dx << " " << mouse_dy << std::endl;
  // std::cout <<  "curx, cury: " << cur_pos_x << " " << cur_pos_y << std::endl;
  cur_pos_x -= center_[0];
  cur_pos_y -= center_[1];
  float r = sqrt(cur_pos_x*cur_pos_x + cur_pos_y*cur_pos_y);
  // rotate clockwise or anticlockwise.
  float cross_product = mouse_dx * pos_y - mouse_dy*pos_x;
  if (cross_product > 0) cross_product = -1.0; else cross_product = 1.0;
  
  if (r > 1e-6) {
    // Get normaized.
    cur_pos_x /= r; cur_pos_y /= r;
    // get dx and dy that are vertical to [cur_pos_x, cur_pos_y]
    float dotproduct = mouse_dx * cur_pos_x + mouse_dy * cur_pos_y;
    float vdx = mouse_dx - dotproduct*cur_pos_x;
    float vdy = mouse_dy - dotproduct*cur_pos_y;
    // Get magnitute of v;
    omega_ += cross_product*sqrtf(vdx*vdx + vdy*vdy) / r / dt *0.05;
    if (std::abs(omega_) > 10.0) {
      omega_ = getsgn(omega_)*10.0;
    }
  }
}

// is the passed in point inside the cylinder?
bool ObstacleCylinderInterac3D::inside(const VEC3F& point) const {
  // translate so that the center of rotation is the origin
  VEC3F pointTransformed = point - center_;
  // rotate *clockwise* back to the original position
  pointTransformed = rotation_ * pointTransformed;
  
  float disxy = pointTransformed[0]*pointTransformed[0] +
                 pointTransformed[1]*pointTransformed[1];
  return (std::abs(pointTransformed[2]) <= half_height_)  && 
         disxy <= radius_squared_;
}

