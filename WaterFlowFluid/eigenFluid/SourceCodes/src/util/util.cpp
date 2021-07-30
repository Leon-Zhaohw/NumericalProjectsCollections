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

#include "Eigen"
#include "util/util.h"

void GetNormalizedRandomVector(const int size, Eigen::VectorXd* vect) {
  (*vect) = Eigen::VectorXd::Random(size);
  vect->normalize();
}

void CompareTensor(const std::vector<Adv_Tensor_Type>& tensor1, 
                   const std::vector<Adv_Tensor_Type>& tensor2) {
  const int basis_dim = tensor1.size();
  if (tensor2.size() != basis_dim) {
    std::cout <<  "Tensor dimention is different, do not compare." << std::endl;
    return;
  }
  std::cout <<  "Begin compare tensor" << std::endl;
  double diff = 0;
  int num_diff = 0;
  for (int k = 0; k < basis_dim; k++) {
    for (int j = 0; j < basis_dim; j++) {
      for (int i = 0; i < basis_dim; i++) {
        double val1 = AccessMatrix(tensor1[k], i, j);
        double val2 = AccessMatrix(tensor2[k], i, j);
        double abs_dif = std::abs(val1 - val2);
        if (abs_dif > 1e-7) {
          num_diff ++;
            std::cout << "Error entry: i: " << i << " j: " 
               << j << " k: " << k << " value1: " << val1 << " value2: " << val2;
        }
        diff += abs_dif;
      }
    }
  }
  std::cout << "Total difference in tensor: " 
      << diff / basis_dim / basis_dim / basis_dim << " number different: " << num_diff;
}

bool PointInBox(const VEC3F& point, const VEC3F& box_center, const VEC3F& length_) {
  VEC3F pointTransformed = point - box_center;
  return ( fabs(pointTransformed[0]) <= length_[0]*0.5  &&
           fabs(pointTransformed[1]) <= length_[1]*0.5  &&
           fabs(pointTransformed[2]) <= length_[2]*0.5 );
}

// Whether two box intersect each other.
bool BoxIntersect(const VEC3F& center1_, const VEC3F& length1_,
                  const VEC3& center2_, const VEC3F& length2_) {
  
  VEC3F low1 = center1_ - length1_*0.5;
  VEC3F high1 = center1_ + length1_*0.5;
  VEC3F low2 = center2_ - length2_*0.5;
  VEC3F high2 = center2_ + length2_*0.5;
  // x
  if (low1[0] > high2[0]) {
    return false;
  }
  if (high1[0] < low2[0]) {
    return false;
  }
  // y
  if (low1[1] > high2[1]) {
    return false;
  }
  if (high1[1] < low2[1]) {
    return false;
  }
  // z
  if (low1[2] > high2[2]) {
    return false;
  }
  if (high1[2] < low2[2]) {
    return false;
  }
  return true;
}

void DrawCube(const float boxXLength, const float boxYLength, const float boxZLength)
{
  glPushMatrix(); 
  glTranslatef(boxXLength,boxYLength,boxZLength);
  glBegin(GL_LINE_STRIP);
  glColor3f(1.f,1.f,1.f);
  glVertex3f(-boxXLength,-boxYLength,-boxZLength);
  glVertex3f(boxXLength,-boxYLength,-boxZLength);
  glVertex3f(boxXLength,boxYLength,-boxZLength);
  glVertex3f(-boxXLength, boxYLength,-boxZLength);
  glVertex3f(-boxXLength,-boxYLength,-boxZLength);
  
  glVertex3f(-boxXLength,-boxYLength,boxZLength);
  glVertex3f(boxXLength,-boxYLength,boxZLength);
  glVertex3f(boxXLength,boxYLength,boxZLength);
  glVertex3f(-boxXLength,boxYLength,boxZLength);
  glVertex3f(-boxXLength,-boxYLength,boxZLength);
  glVertex3f(-boxXLength,boxYLength,boxZLength);
  glVertex3f(-boxXLength,boxYLength,-boxZLength);
  glVertex3f(boxXLength,boxYLength,-boxZLength);
  glVertex3f(boxXLength,boxYLength,boxZLength);
  glVertex3f(boxXLength,-boxYLength,boxZLength);
  glVertex3f(boxXLength,-boxYLength,-boxZLength);
  glEnd();
  glTranslatef(-0.1f, -0.1f, -0.1f);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
    // x axis is red
    glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(10.0f, 0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);

    // y axis is green 
    glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 10.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    
    // z axis is blue
    glColor4f(0.0f, 0.0f, 10.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 0.0f, 10.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
  glEnd();
  glLineWidth(1.0f);
  glPopMatrix();
}
