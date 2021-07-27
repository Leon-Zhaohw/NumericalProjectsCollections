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
#include <float.h>
#include <GL/glut.h>
#include <math.h>
#include <omp.h>
#include <vector>

#include "3D/drawer_3d.h"
#include "3D/particle_3d.h"
#include "Alg/VEC3.h"

void Drawer3D::DrawParticles(const std::vector<Particle3D>& particles,
                             const double dx_, const double dt_, const double ptl_length) {
  // get the color for every points.
  std::vector<Eigen::Vector4f> colors;
  const int ptls_size = particles.size();
  colors.reserve(ptls_size);
  
#pragma omp parallel for
  for (int i = 0; i < ptls_size; i++) {
    VEC3F velocity = particles[i].velocity;
    float velomag = sqrtf(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
    Eigen::Vector4f cols;
    cols[3] = velomag*5.0;
    getColorAtValue(velomag*0.8, cols[0], cols[1], cols[2]);
    colors[i] = cols;
  }
  float x,y,z;
  // glColor3f(1.0,1.0,1.0);
  // dx_ is the length of the cell.
  const float dt0 = dt_ / dx_;
  glPointSize(4.5);
  glBegin(GL_POINTS);
  for (int i = 0; i < ptls_size; i++) {
    x = (particles[i].position[0] ) * dx_;
    y = (particles[i].position[1] ) * dx_;
    z = (particles[i].position[2] ) * dx_;
    
    glColor4f(colors[i][0], colors[i][1], colors[i][2], colors[i][3]);
    // glColor3f(0,0,0);
    glVertex3f(x, y, z);

    //glBegin(GL_LINES);
    //glVertex3f(x, y, z);
    //glVertex3f(x + particles[i].velocity[0] *  dt0 * ptl_length,
    //           y + particles[i].velocity[1] * dt0 * ptl_length,
    //           z + particles[i].velocity[2] * dt0 * ptl_length);
    //glEnd();
  }
  glEnd();
}

void Drawer3D::subSampleVector(const Eigen::VectorXd& vecIn, Eigen::VectorXd& vecOut) {
  if (vecIn.size() <= 500) {
    vecOut.resize(vecIn.size());
    vecOut = Eigen::VectorXd::Map(vecIn.data(), vecIn.size());
    return;
  } else {
    vecOut.resize(500);
    vecOut.setZero();
    const int inSize = vecIn.size();
    float fIdx = 0;
    const float deltaF = 1.0 / 500.0;
    for (int i = 0; i < 500; i++) {
      int idxBegin = static_cast<int>(fIdx * inSize);
      int idxEnd = static_cast<int>((fIdx + deltaF)*inSize);
      double aveAbsVal = 0;
      for (int j = idxBegin; j<= idxEnd; j++)
        aveAbsVal += vecIn[j];
      aveAbsVal /= static_cast<double>(idxEnd - idxBegin);
      vecOut[i] = aveAbsVal;
      fIdx += deltaF;
    }
  }
}

void Drawer3D::DrawSpectogram(const Eigen::VectorXd& coefficients,
                              const double multi_factor) {
  Eigen::VectorXd sampledCoeff;
  subSampleVector(coefficients, sampledCoeff);
  
  float start_x = 0.00;
  float start_y = 0.00;
  
  spectogramCoeff_.push_back(sampledCoeff);
  const float deltaY = 1.0 / spectogramCoeff_[0].size();
  const float deltaX = 1.0 / 600.0;
  for (int i = 0; i < spectogramCoeff_.size(); i++) {
    const Eigen::VectorXd& curCoeff = spectogramCoeff_[i];
    float invMax = 0;
    if (curCoeff.maxCoeff() > 1e-12) {
      invMax = 1.0 / curCoeff.maxCoeff();
    }
    
    start_y = 0.01;
    
    for (int j = 0; j < curCoeff.size(); j++) {
      // Draw one spec.
      float color[3] = {1.0,1.0,1.0};
      getColorAtValue(sqrt((std::abs(curCoeff[j]))*invMax), color[0], color[1], color[2]);
      //getColorAtValue(static_cast<float>(j)*0.001, color[0], color[1], color[2]);
      glColor3f(color[0], color[1], color[2]);
      glBegin(GL_QUADS);
      glVertex2d(start_x, start_y);
      glVertex2d(start_x, start_y + deltaY);
      glVertex2d(start_x + deltaX, start_y + deltaY);
      glVertex2d(start_x + deltaX, start_y);
      glEnd();
      start_y += deltaY;
    }
    start_x += deltaX;
  }  
}

// Visualize the coefficients as a bar graph.
void Drawer3D::DrawCoefficients(const Eigen::VectorXd& coefficients,
                        const double dx_, const double multi_factor) {
  const int num_coefficients = coefficients.size();
  const double step_size = 1.0 / num_coefficients;
  const double start_y = 0.5;
  double start_x = 0;
  // Plot a quad per coefficients.
  glColor3f(0.0, 0.0, 0.0);
  for (int i = 0; i < num_coefficients; i++) {
    double height = coefficients[i] * multi_factor;
    glBegin(GL_QUADS);
    glVertex2d(start_x, start_y);
    glVertex2d(start_x, start_y + height);
    glVertex2d(start_x + step_size, start_y + height);
    glVertex2d(start_x + step_size, start_y);
    glEnd();
    start_x += step_size;
  }
}
