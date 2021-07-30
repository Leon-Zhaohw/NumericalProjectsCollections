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

#ifndef FLUID2D_H
#define FLUID2D_H

#include <math.h>
#include "2D/FIELD2D.h"
#include "2D/VFIELD2D.h"
#include "setting.h"

#include "2D/particle_2d.h"
#include "2D/drawer_2d.h"
#define GLUT_DISABLE_ATEXIT_HACK

#ifdef APPLE
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

struct StamBoundarySetting {
  StamBoundarySetting (const bool top, const bool bottom, const bool left,
    const bool right):wall_top(top), wall_bottom(bottom),wall_left(left),
    wall_right(right){}
    
  bool wall_bottom;
  bool wall_top;
  bool wall_left;
  bool wall_right;
};

class FLUID2D{
public:
  FLUID2D();
  FLUID2D(int _xRes, int _yRes, const int num_particles, const double dt_, 
          const double buoyancy, const double added_smoke_density,
          const StamBoundarySetting& BoundarySetting);
  void step();
  
  void addSmokeTest(const int xpos, const int ypos, const int width, const int height);
  void drawDensity();
  void drawVelocity();
  void drawVelocityX();
  void drawVelocityY();
  void ClearDens();
  void ClearVelo();
  void DrawObstal();
  void resetObsvelo();
  FIELD2D* get_density(){return &density;}
  VFIELD2D* get_velocity(){return &velocity;}
  void ReSeedParticles();
  void DrawParticles(const double ptl_length);
  bool addforce;
  void MouseAddForce(int xpos, int ypos,  const float fx, const float fy);
  bool quit_ = false;
  
private:

  Real _dx;
  const double _dt;
  void advect_stam();
  void addBuoyancy(Real* field);
  void project();
  void solvePoission(FIELD2D& x, FIELD2D& b,unsigned char* obs);
  void addField(FIELD2D&, int x, int y, const int width, const int height);
  void addVoriticity();
  void setObstaclBoundary();
  void zeroBoundary(FIELD2D& x, unsigned char* obs);
  void MultPossisonStencil(FIELD2D& in, FIELD2D& out);
  void AdvectParticles();
  void ForwardTransformToFrequency(const VFIELD2D& v, Real* uFreq, Real* vFreq);
  void InverseTransformToVelocity(Real* uFreq, Real* vFreq, VFIELD2D* v);
  void projectDCT();
  void computeStreamFunction(Real* uFreq, FIELD2D* stream);
  void InitCheckBoardPattern();
  unsigned char* obstacle;;
  // Class to draw stuff.
  std::unique_ptr<Drawer2D> drawer_;
  int xRes;
  int yRes;
  int total;
  FIELD2D density;
  FIELD2D density_old;
  FIELD2D temp1;
  FIELD2D temp2;
  
  FIELD2D pressure;
  FIELD2D divergence;
  FIELD2D streamF;
  
  VFIELD2D velocity;
  VFIELD2D velocity_old;
  VFIELD2D force;

  FIELD2D residual;
  FIELD2D direction;
  FIELD2D _q;

  FIELD2D voriticity;

  bool wall_right_;
  bool wall_left_;
  bool wall_top_;
  bool wall_bottom_;
  int Maxiter;
  // simulation constants
  //Real _dt;
  Real _dtOld;
  Real* vxFreq_;
  Real* vyFreq_;
  
  const double buoyancy_;
  const double added_smoke_density_;
  Real _vorticityEps;
  Real _heatDiffusion;
  Real _solverEps;
    
  const int num_particles_;
  std::vector<Particle2D> particles_;
};

#endif
