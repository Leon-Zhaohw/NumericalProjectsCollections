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

#include <fstream>
#include <float.h>

#include "util/stringprintf.h"
#include "3D/obstacle_sphere_3d.h"
#include "3D/VECTOR3_FIELD_3D.h"

void ObstacleSphere3D::Init() {
  radius_squared_ = radius_*radius_;
  // velocity_ = VEC3F(0.0,0.0,0.0);
  ComputeAABB(this->dx_, this->gridRes_);
}

void ObstacleSphere3D::CalculateNormalForce(
    const VECTOR3_FIELD_3D& velocity, char * check, const double dt, const double cell_len, const double force_amp,
    VECTOR3_FIELD_3D* normal_component) const {
  int xRes = velocity.xRes();
  int yRes = velocity.yRes();
  int zRes = velocity.zRes();
  int slab = velocity.slabSize();
  
  if(xRes != normal_component->xRes()) { std::cout << "The dimention" << " must be equal" << std::endl; exit(0);};
  if(yRes != normal_component->yRes()) { std::cout << "The dimention " << " must be equal " << std::endl; exit(0);};
  if(zRes != normal_component->zRes()) { std::cout  << "The dimention"  << " must be equal" << std::endl; exit(0);}
  
  const double amp = force_amp / dt;
#pragma omp parallel for
  for (int k = AABBStart_[2]; k < AABBEnd_[2]; k++) {
    for (int j = AABBStart_[1]; j < AABBEnd_[1]; j++) {
      for (int i = AABBStart_[0]; i < AABBEnd_[0]; i++) {
        const int index = i + j*xRes + k*slab;
        // skip the cells that have already been checked by other obstacles.
        if (check[index]) continue;
        
        VEC3F pt((static_cast<float>(i) + 0.5)*cell_len,
                (static_cast<float>(j) + 0.5)*cell_len,
                (static_cast<float>(k) + 0.5)*cell_len);
        // Cell is inside the obstacle.
        if (inside(pt)) {
          (*normal_component)[index][0] -= amp*(velocity[index][0] - velocity_[0]);
          (*normal_component)[index][1] -= amp*(velocity[index][1] - velocity_[1]);
          (*normal_component)[index][2] -= amp*(velocity[index][2] - velocity_[2]);
          check[index] = 1;
          continue;
        }
        int  left = 0, right = 0, top = 0, bottom = 0, front = 0, back = 0;
        // left
        pt[0] -= cell_len;
        if (this->inside(pt)) left ++;
        pt[0] += cell_len;
        // right
        pt[0] += cell_len;
        if (this->inside(pt)) right ++;
        pt[0] -= cell_len;
        // top
        pt[1] += cell_len;
        if (this->inside(pt)) top++;
        pt[1] -= cell_len;
        // bottom
        pt[1] -= cell_len;
        if (this->inside(pt)) bottom++;
        pt[1] += cell_len;
        // front
        pt[2] += cell_len;
        if (this->inside(pt)) front++;
        pt[2] -= cell_len;
        // back
        pt[2] -= cell_len;
        if (this->inside(pt)) back++;
        pt[2] += cell_len;
        int total = top + bottom + left + right + front + back;
        // No obstacle nearby.
        if (total == 0) continue;
        else check[index] = 1;
        // Left side. x component.
        if (left != 0) {
          (*normal_component)[index][0] -= amp*(velocity[index][0] - velocity_[0]);
        }
        // right side. x component.
        if (right != 0) {
          (*normal_component)[index][0] -= amp*(velocity[index][0] - velocity_[0]);
        }
        // top, y component.
        if (top != 0) {
          (*normal_component)[index][1] -= amp*(velocity[index][1] - velocity_[1]);
        }
        // bottom, y component.
        if (bottom != 0) {
          (*normal_component)[index][1] -= amp*(velocity[index][1] - velocity_[1]);
        }
        // front, z component.
        if (front != 0) {
          (*normal_component)[index][2] -= amp*(velocity[index][2] - velocity_[2]);
        }
        // back, z component.
        if (back != 0) {
          (*normal_component)[index][2] -= amp*(velocity[index][2] - velocity_[2]);
        }
      }
    }
  }
}

VEC3F ObstacleSphere3D::getVelocityAt(const VEC3F& point) const {
  if (!this->inside(point)) {
    return velocity_;
  } else {
    return VEC3F(0,0,0);
  } 
}

void ObstacleSphere3D::ComputeAABB(const double dx, const VEC3I& gridRes) {
  // Get corners in local frame.
  std::vector<VEC3F> corners;
  corners.push_back(VEC3F(-radius_, -radius_, -radius_));
  corners.push_back(VEC3F(-radius_, -radius_, radius_));
  corners.push_back(VEC3F(radius_, -radius_, -radius_));
  corners.push_back(VEC3F(radius_, -radius_, radius_));
  corners.push_back(VEC3F(-radius_, radius_, -radius_));
  corners.push_back(VEC3F(-radius_, radius_, radius_));
  corners.push_back(VEC3F(radius_, radius_, -radius_));
  corners.push_back(VEC3F(radius_, radius_, radius_));
  // Transform into global frame.
  for (int i = 0; i < corners.size(); i++) {
    VEC3F transfomed = corners[i];
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

void ObstacleSphere3D::ZeroOutVelocityDensityObstalce(const double cell_len,
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
         (*velocity)[idx] = velocity_;
         //(*density)[idx] = 0;
        }
      }
    }
  }
}

void ObstacleSphere3D::RasterToGrid(unsigned char* obstacle_filed) {
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

#define EXPECT_STRING(str) \
in >> temp; \
if (temp != str) { \
  std::cout << "obstacle_sphere_3d.cpp " << __LINE__ << " FATAL: " <<  "Error: " << temp << std::endl; exit(0); \
}\
temp.clear();

void ObstacleSphere3D::GetFromFile(const std::string& fname) {
  
  std::ifstream in(fname);
  if (!in.is_open()) {
    std::cout << "obstacle_sphere_3d.cpp " << __LINE__ << " FATAL: " <<  "Cannot open: " << fname << std::endl; exit(0);
  }
  std::string temp;
  EXPECT_STRING("sphere")
  EXPECT_STRING("center")
  in >> center_[0]; in >> center_[1]; in >> center_[2];
  EXPECT_STRING("radius")
  in >> radius_;
  EXPECT_STRING("color")
  in >> color_[0]; in >> color_[1]; in >> color_[2];
  EXPECT_STRING("velocity");
  
  in >> velocity_[0]; in >> velocity_[1]; in >> velocity_[2];
}

#undef EXPECT_STRING
// draw to OpenGL
void ObstacleSphere3D::draw() const {
  
  glPushMatrix();
  glTranslatef(center_[0], center_[1], center_[2]);
  // glPushMatrix();
  // glScalef((AABBEnd_[0] - AABBStart_[0])*dx_, (AABBEnd_[1] - AABBStart_[1])*dx_, 
  //         (AABBEnd_[2] - AABBStart_[2])*dx_);
  // glutSolidCube(1.0);
  // glPopMatrix();
  // OpenGL uses degrees, so we need to convert
  // glRotatef(this->thetaDegrees(), _rotationAxis[0], _rotationAxis[1], _rotationAxis[2]);
  glScalef(radius_, radius_, radius_);
  glColor4f(0.1, 0.7, 0.8, 0.6);
  GLUquadric *quad;
  quad = gluNewQuadric();
  gluSphere(quad,1, 20, 20);
  glPopMatrix();  
}

void ObstacleSphere3D::OutPutGeometryToPBRT(std::ofstream& outfile) {
  if (!outfile.is_open()) {
    std::cout << "obstacle_sphere_3d.cpp " << __LINE__ << " FATAL: " <<  "Cannot open tensor file." << std::endl; exit(0);
  }
  outfile << "AttributeBegin \n TransformBegin\n";
  outfile << "Translate " << center_[0] << " " << center_[1] << " " << center_[2] << "\n";
  // Texture and color of the obstacle.
  outfile << "Material \"matte\" \"color Kd\" [ " << color_[0] << " " <<  color_[1]
      << " " <<  color_[2] << " ]\n";
  outfile << "Shape \"sphere\" \"float radius\" " << radius_ ;
  outfile <<" TransformEnd \n AttributeEnd\n";
}

void ObstacleSphere3D::OutPutGeometryToPBRT(const std::string &fname) {
  const std::string m_fname = StringPrintf("%ssphere.pbrt", fname.c_str());
  std::ofstream outfile(m_fname.c_str());
  OutPutGeometryToPBRT(outfile);
  outfile.close();
}

// is the passed in point inside the sphere?
bool ObstacleSphere3D::inside(const VEC3F& point) const {
  // translate so that the center of rotation is the origin
  VEC3F pointTransformed = point - center_;
  // Distance.
  float dist = pointTransformed[0]*pointTransformed[0] + 
               pointTransformed[1]*pointTransformed[1] + 
               pointTransformed[2]*pointTransformed[2];
  if (dist < radius_squared_) {
    return true;
  } else {
    return false;
  }
}
