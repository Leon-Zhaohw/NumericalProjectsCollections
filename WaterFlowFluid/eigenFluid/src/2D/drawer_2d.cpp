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
// #include <glog/logging.h>
#include <GL/glut.h>
#include <vector>

#include "2D/drawer_2d.h"
#include "2D/FIELD2D.h"

void Drawer2D::DrawDensity(const FIELD2D& density, const int xRes, const int yRes,
                           const double dx_) {
  float x, y, h, d00, d01, d10, d11;
  h = dx_;
  glBegin ( GL_QUADS );
  for (int j = 0; j < yRes - 1; j++) {
    y = (j - 0.5f) * h;
      for (int i = 0; i < xRes - 1; i++) {
        x = (i - 0.5f) * h;
        int index = i + j * xRes;
        d00 = 1.0 - density[index];
        d01 = 1.0 - density[index + xRes];
        d10 = 1.0 - density[index + 1];
        d11 = 1.0 - density[index + 1 + xRes];

        glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
        glColor3f ( d10, d10, d10 ); glVertex2f ( x + h, y );
        glColor3f ( d11, d11, d11 ); glVertex2f ( x + h, y + h );
        glColor3f ( d01, d01, d01 ); glVertex2f ( x, y + h );
      }
  }
  glEnd ();
}

void Drawer2D::DrawVort(const FIELD2D& vort, const int xRes, const int yRes,
                const double dx_) {
  double maxabs = -1.0;
  for (uint i = 0; i < xRes*yRes; i++) {
    double absval = std::abs(vort[i]);
    if (absval > maxabs) {
      maxabs = absval;
    }
  }
  
  float x, y, h, d00, d01, d10, d11;
  h = dx_;
  glBegin ( GL_QUADS );
  for (int j = 0; j < yRes - 1; j++) {
    y = (j - 0.5f) * h;
      for (int i = 0; i < xRes - 1; i++) {
        x = (i - 0.5f) * h;
        int index = i + j * xRes;
        d00 = std::abs(vort[index]) / maxabs;
        d01 = std::abs(vort[index + xRes]) / maxabs;
        d10 = std::abs(vort[index + 1]) / maxabs;
        d11 = std::abs(vort[index + 1 + xRes]) / maxabs;
        Eigen::Vector3f cols0, cols1, cols2, cols3;
        
        getColorAtValue(d00, cols0[0], cols0[1], cols0[2]);
        getColorAtValue(d01, cols1[0], cols1[1], cols1[2]);
        getColorAtValue(d10, cols2[0], cols2[1], cols2[2]);
        getColorAtValue(d11, cols3[0], cols3[1], cols3[2]);
        
        glColor3f ( cols0[0], cols0[1], cols0[2] ); glVertex2f ( x, y );
        glColor3f ( cols2[0], cols2[1], cols2[2] ); glVertex2f ( x + h, y );
        glColor3f ( cols3[0], cols3[1], cols3[2] ); glVertex2f ( x + h, y + h );
        glColor3f ( cols1[0], cols1[1], cols1[2] ); glVertex2f ( x, y + h );
      }
  }
  glEnd ();
  
}
  
void Drawer2D::DrawVelocity(const VFIELD2D& velocity, const int xRes,
                            const int yRes, const double dx_) {
  float x, y, h;
  h = dx_;
  glColor3f ( 1.0f, 1.0f, 1.0f );
  glLineWidth ( 1.0f );
  glBegin ( GL_LINES );
  for (int j = 0; j < yRes - 1; j++) {
    y = (j - 0.5f) * h;
    for (int i = 0 ;i < xRes - 1; i++) {
      x = (i - 0.5f) * h;
      int index = i + j * xRes;
      glVertex2f ( x, y );
      glVertex2f ( x + velocity[index][0], y + velocity[index][1]);
    }
  }
  glEnd ();
}

void Drawer2D::DrawParticles(const std::vector<Particle2D>& particles,
                             const double dx_, const double dt_,
                             const double ptl_length) {
    // get the color for every points.
  std::vector<Eigen::Vector3f> colors;
  const int ptls_size = particles.size();
  // colors.reserve(ptls_size);
  
  for (int i = 0; i < ptls_size; i++) {
    VEC2 velocity = particles[i].velocity;
    float velomag = sqrtf(velocity[0]*velocity[0] + velocity[1]*velocity[1]);
    Eigen::Vector3f cols;
    getColorAtValue(velomag*7.0, cols[0], cols[1], cols[2]);
    colors.push_back(cols);
  }

  float x,y;
  const float dt0 = dt_ / dx_;
  for (int i = 0; i < ptls_size; i++) {
    x = (particles[i].position[0] - 0.5) * dx_;
    y = (particles[i].position[1] - 0.5) * dx_;
    float lx = particles[i].velocity[0] *  dt0 * ptl_length;
    float ly = particles[i].velocity[1] * dt0 * ptl_length;
    // glColor3f(colors[i][0], colors[i][1], colors[i][2]);
    glColor3f(0,0,0);
    glPointSize(1.5);
    glBegin(GL_POINTS);
    glVertex2f(x, y);
    glEnd();
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glVertex2f(x + 0.5*lx,y + 0.5*ly);
    glVertex2f(x - 0.5*lx,
               y - 0.5*ly);
    glEnd();
  }
}

void Drawer2D::DrawObstacles(const Obstacle2D& obstacle, const double dx_) {
  glColor3f(0.0,1.0,1.0);
  const int xRes = obstacle.xRes_;
  const int yRes = obstacle.yRes_;
  glBegin(GL_QUADS);
  for (int j = 0; j < yRes; j++) {
    for (int i = 0; i < xRes; i++) {
      if (obstacle.Access(i, j)) {
        glVertex2f(i*dx_, j*dx_);
        glVertex2f((i+1)*dx_, j*dx_);
        glVertex2f((i+1)*dx_, (j+1)*dx_);
        glVertex2f(i*dx_, (j+1)*dx_);
      }
    }
  }
  glEnd();
}

void Drawer2D::Draw2Coefficients(const Eigen::VectorXd& coefficients, const Eigen::VectorXd& col,
                                 const double dx_, const double multi_factor, const double offset) {
  const int num_coefficients = coefficients.size();
  const double step_size = 1.0 / num_coefficients;
  const double start_y = 0.5;
  double start_x = 0;
  // Plot a quad per coefficients.
  // glColor3f(1.0, 1.0, 1.0);
  glPushMatrix();
  glTranslatef(0.0,offset,0.0);
  for (int i = 0; i < num_coefficients; i++) {
    double height = coefficients[i] * multi_factor;
    glBegin(GL_QUADS);
    if (col[i] < 0) {
      glColor3f(1.0, 0.0, 0.0);
    } else {
      glColor3f(0.0, 1.0, 0.0);
    }
    glVertex2d(start_x, start_y);
    glVertex2d(start_x, start_y + height);
    glVertex2d(start_x + step_size, start_y + height);
    glVertex2d(start_x + step_size, start_y);
    glEnd();
    start_x += step_size;
  }
  glPopMatrix();
}

void Drawer2D::DrawCoefficients(const Eigen::VectorXd& coefficients,
                                const double dx_, const double multi_factor, const double offset) {
  const int num_coefficients = coefficients.size();
  const double step_size = 1.0 / num_coefficients;
  const double start_y = 0.5;
  double start_x = 0;
  // Plot a quad per coefficients.
  glColor3f(1.0, 1.0, 1.0);
  glPushMatrix();
  glTranslatef(0.0,offset,0.0);
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
  glPopMatrix();
}
