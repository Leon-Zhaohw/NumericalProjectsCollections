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
#include <fftw3.h>
#include <fstream>

#include "2D/FLUID2D.h"
#include "2D/laplacian_basis_2D.h"
#include "2D/three_dirichlet_one_neumann.h"
#include "util/util.h"
#include "util/stringprintf.h"

#define MAX_ITERACTION 1000


FLUID2D::FLUID2D(int _xres, int _yres, const int num_particles, const double dt_,
                 const double buoyancy, const double added_smoke_density,
                 const StamBoundarySetting& BoundarySetting):xRes(_xres), yRes(_yres),
                 num_particles_(num_particles), _dt(dt_), buoyancy_(buoyancy), 
                 added_smoke_density_(added_smoke_density) {  
                  const int num_basis_root = 20;
  const int num_basis = num_basis_root * num_basis_root;
  Eigen::VectorXd coeff;
  coeff.resize(num_basis); coeff.setZero();
  coeff[0] = 0.3;
  coeff[1] = 0.3;
  // ThreeDirichletOneNeumann baa(_xres, num_basis_root);
  // baa.FillAdvectionTensor(&Adv_tensor_);
  _heatDiffusion = 1e-3;
  _vorticityEps = 2.0;
  _solverEps = 1e-10;

  density = FIELD2D(xRes,yRes);
  density_old = FIELD2D(xRes,yRes);
  temp1 = FIELD2D(xRes, yRes);
  temp2 = FIELD2D(xRes, yRes);
  
  pressure = FIELD2D(xRes, yRes);
  streamF  = FIELD2D(xRes, yRes);
  
  divergence = FIELD2D(xRes, yRes);

  velocity = VFIELD2D(xRes, yRes);
  velocity_old = VFIELD2D(xRes, yRes);
    
  force = VFIELD2D(xRes,yRes);
  residual = FIELD2D(xRes,yRes);
  direction = FIELD2D(xRes,yRes);
  _q        = FIELD2D(xRes,yRes);

  voriticity = FIELD2D(xRes, yRes);
  vxFreq_ = new Real[xRes*yRes];
  memset(vxFreq_, 0x00, sizeof(Real)*xRes*yRes);
  vyFreq_ = new Real[xRes*yRes];
  memset(vyFreq_, 0x00, sizeof(Real)*xRes*yRes);
  
  int maxres = fmax(xRes, yRes);
  _dx = 1.0/(Real)maxres;

  obstacle = new unsigned char[xRes*yRes];
  for(int i=0;i<xRes*yRes;i++)
    obstacle[i] = EMPTY;
  Maxiter = MAX_ITERACTION;
  total = xRes*yRes;
    
  wall_right_ = BoundarySetting.wall_right;
  wall_left_ = BoundarySetting.wall_left;
  wall_top_ = BoundarySetting.wall_top;
  wall_bottom_ = BoundarySetting.wall_bottom;
  
  // Initialize particles.
  for (int i = 0; i < num_particles_; i++) {
    int px = std::rand() % xRes;
    int py = std::rand() % yRes;
    Particle2D particle(px + 0.5, py + 0.5);
    particles_.push_back(particle);
  }
  // Initialize drawer.
  drawer_.reset(new Drawer2D());
  
  // set side obstacles
  if (wall_right_) {
    for (int i = 0; i < yRes; i++) {
      int index = i * xRes + xRes - 1;
      obstacle[index] = 1;
    }
  }
  if (wall_left_) {
    for (int i = 0; i < yRes; i++) {
      int index = i * xRes;
      obstacle[index] = 1;
    }
  }
  if (wall_top_) {
    const int begin_idx = xRes*(yRes - 1);
    for (int i = 0; i < xRes; i++) {
      obstacle[begin_idx + i] = 1;
    }
  }
  if (wall_bottom_) {
    for (int i = 0; i < xRes; i++) {
      obstacle[i] = 1;
    }
  }
  addforce = true;
  
  //InitCheckBoardPattern();
}

void FLUID2D::step() {
  force.clear();
  // wipe boundaries
  velocity.set_zero_border();
  density.setZeroBorder();
  
 // density.CheckNegative();
  if (addforce) {
    // compute the forces
    addBuoyancy(density.getdata());
    velocity.axpy(_dt, force);
  }
  //addVoriticity();
  //velocity.axpy(_dt, force);
  force.clear();
  // advect things
  
  advect_stam();
  AdvectParticles();
  // run the solvers
  // project();
  projectDCT();
  // currentTime += _dt;
  // 
  // _totalTime += goalTime;
  // _totalSteps++;
}

//advect every thing
void FLUID2D::advect_stam() {
  if (wall_left_) {
    velocity.SetZeroLeft();
  } else {
    velocity.CopyBorderLeft();
  }
  
  if (wall_right_) { 
    velocity.SetZeroRight();
  } else {
    velocity.CopyBorderRight();
  }
  
  if (wall_top_) {
    velocity.SetZeroTop();
  } else {
    velocity.CopyBorderTop();
  }
  
  if (wall_bottom_) {
    velocity.SetZeroBottom();
  } else {
    velocity.CopyBorderBottom();
  }

  const Real dt0 = _dt / _dx;
  
  velocity.swapPointer(velocity_old);
  density.swapPointer(density_old);

  VFIELD2D::advect(dt0, velocity_old, density_old, density);
  VFIELD2D::advect(dt0, velocity_old, velocity_old, velocity);

  velocity.set_zero_border();
  density.setZeroBorder();
   
  if (wall_left_) {
    velocity.SetZeroLeft();
  } else {
    velocity.CopyBorderLeft();
  }
  
  if (wall_right_) { 
    velocity.SetZeroRight();
  } else {
    velocity.CopyBorderRight();
  }
  
  if (wall_top_) {
    velocity.SetZeroTop();
  } else {
    velocity.CopyBorderTop();
  }
  
  if (wall_bottom_) {
    velocity.SetZeroBottom();
  } else {
    velocity.CopyBorderBottom();
  }
}

//add buoyancey to x direction
void FLUID2D::addBuoyancy(Real* field) {
  Real beta = buoyancy_;
  int index = 0;
  if(beta==0.) return;
  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++, index++) {
      force[index][0] += beta * field[index];
    }
}

void FLUID2D::ForwardTransformToFrequency(const VFIELD2D& v, Real* uFreq, Real* vFreq) {
  for (int i = 0; i < xRes*yRes; i++) {
    uFreq[i] = v[i][0];
    vFreq[i] = v[i][1];
  }
  double invTotalSize_ = 1.0 / static_cast<double>(xRes*yRes);
  fftw_plan plan_x_2D;
  // x: RODFT10, y: REDFT10
  plan_x_2D = fftw_plan_r2r_2d(yRes, xRes, uFreq, uFreq,
                               FFTW_REDFT10, FFTW_RODFT10, FFTW_ESTIMATE);
  fftw_execute(plan_x_2D);
  fftw_destroy_plan(plan_x_2D);
  
  fftw_plan plan_y_2D;
  // x: REDFT10, y: RODFT10
  plan_y_2D = fftw_plan_r2r_2d(yRes, xRes, vFreq, vFreq,
                               FFTW_RODFT10, FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(plan_y_2D);
  fftw_destroy_plan(plan_y_2D);
  // normalized.
  for (int i = 0; i < xRes*yRes; i++) {
    uFreq[i] *= invTotalSize_*0.25;
    vFreq[i] *= invTotalSize_*0.25;
  }
}

void FLUID2D::InverseTransformToVelocity(Real* uFreq, Real* vFreq, VFIELD2D* v) {
  fftw_plan plan_x_2D;
  plan_x_2D = fftw_plan_r2r_2d(yRes, xRes, uFreq, uFreq,
                               FFTW_REDFT01, FFTW_RODFT01, FFTW_ESTIMATE);
  fftw_execute(plan_x_2D);  
  fftw_destroy_plan(plan_x_2D);
  
  fftw_plan plan_y_2D;
  plan_y_2D = fftw_plan_r2r_2d(yRes, xRes, vFreq, vFreq,
                               FFTW_RODFT01, FFTW_REDFT01, FFTW_ESTIMATE);
  fftw_execute(plan_y_2D);  
  fftw_destroy_plan(plan_y_2D);
  
  for (int i = 0; i < xRes*yRes; i++) {
    (*v)[i][0] = uFreq[i];
    (*v)[i][1] = vFreq[i];
  }
}

void FLUID2D::computeStreamFunction(Real* uFreq, FIELD2D* stream) {
  // compute the stream function using u, assume the velocity is div free,
  // so kx*u(kx,ky) ky*v(kx,ky) = 0
  // u DST on x, DCT on y.
  Real* phi = new double[xRes*yRes];
  memset(phi, 0x00, sizeof(double)*xRes*yRes);
  for (int ky = 0; ky < yRes; ky++) {
    for (int kx = 0; kx < xRes; kx++) {
      const int index = kx + ky*xRes;
      double uhat = 0;
      if (kx == 0) uhat = 0; else uhat = uFreq[index - 1];
      if (ky != 0 && kx != 0) {
        phi[index - 1 - xRes] = 1.0 / static_cast<double>(ky) * uhat;
      } else {
        phi = 0;
      }
    }
  }
  fftw_plan plan_x_2D;
  plan_x_2D = fftw_plan_r2r_2d(yRes, xRes, phi, phi,
                               FFTW_RODFT01, FFTW_RODFT01, FFTW_ESTIMATE);
  fftw_execute(plan_x_2D);  
  fftw_destroy_plan(plan_x_2D);
  
  for (int i = 0; i < xRes*yRes; i++) {
    (*stream)[i] = phi[i];
  }
  delete [] phi;
}

void FLUID2D::projectDCT() {
  const double invxRes_ = 1.0 / static_cast<double>(xRes);
  const double invyRes_ = 1.0 / static_cast<double>(yRes);
  
  ForwardTransformToFrequency(velocity, vxFreq_, vyFreq_);
  
  for (int ky = 0; ky < yRes; ky++) {
    for (int kx = 0; kx < xRes; kx++) {
      const int index = kx + ky*xRes;
      
      // Wavenumber.
      Eigen::Vector2d K(static_cast<double>(kx)*invxRes_,
                        static_cast<double>(ky)*invyRes_);
      
      Eigen::Vector2d W; // (vxFreq_[index], vyFreq_[index], vzFreq_[index]);
      if (kx == 0) W[0] = 0; else W[0] = vxFreq_[index - 1];
      if (ky == 0) W[1] = 0; else W[1] = vyFreq_[index - xRes];
      
      if (ky!= 0 || kx != 0) {
        W = W - 1.0 / K.squaredNorm() * (K.dot(W))*K;
      }
      
      if (kx != 0) vxFreq_[index - 1] = W[0];
      if (ky != 0) vyFreq_[index - xRes] = W[1];
    }
  }
  InverseTransformToVelocity(vxFreq_, vyFreq_, &velocity);
}

void FLUID2D::project() {
  setObstaclBoundary();
    
  if (wall_left_) {
    velocity.SetZeroLeft();
  } else {
    velocity.SetNeumannLeft();
  }
  
  if (wall_right_) { 
    velocity.SetZeroRight();
  } else {
    velocity.SetNeumannRight();
  }
  
  if (wall_top_) {
    velocity.SetZeroTop();
  } else {
    velocity.SetNeumannTop();
  }
  
  if (wall_bottom_) {
    velocity.SetZeroBottom();
  } else {
    velocity.SetNeumannBottom();
  }

  //div

  // calculate divergence
  int index = xRes + 1;
  int x,y;
  for (y = 1; y < yRes - 1; y++, index += 2) {
    for (x = 1; x < xRes - 1; x++, index++) {
      
      Real xright = velocity[index + 1][0];
      Real xleft  = velocity[index - 1][0];
      Real yup    = velocity[index + xRes][1];
      Real ydown  = velocity[index - xRes][1];
      if(obstacle[index+1]) xright = - velocity[index][0];
      if(obstacle[index-1]) xleft  = - velocity[index][0];
      if(obstacle[index+xRes]) yup    = - velocity[index][1];
      if(obstacle[index-xRes]) ydown  = - velocity[index][1];

      divergence[index] = -_dx * 0.5f * (
        xright - xleft +
        yup - ydown);
        pressure[index] = 0.0f;
    }
  }
  pressure.copyBorderAll();
  solvePoission(pressure, divergence, obstacle);

  Real invDx = 1.0/_dx;
  index = xRes + 1;
  for (y = 1; y < yRes - 1; y++, index += 2) {	
    for (x = 1; x < xRes - 1; x++, index++) {
      //int index = x + y * (xRes);
      if (x == 1)
        velocity[index][0] -= (pressure[index + 1] - pressure[index]) * invDx;
      else if (x == xRes - 2)
        velocity[index][0] -= (pressure[index] - pressure[index - 1]) * invDx;
      else
        velocity[index][0] -= 0.5 * (pressure[index + 1] - pressure[index - 1]) * invDx;
      if (y == 1)
        velocity[index][1] -= (pressure[index + xRes]  - pressure[index]) * invDx;
      else if (y == yRes - 2)
        velocity[index][1] -= (pressure[index]  - pressure[index - xRes]) * invDx;
      else
        velocity[index][1] -= 0.5f * (pressure[index + xRes]  - pressure[index - xRes]) * invDx;
    }
  }
}

void FLUID2D::solvePoission(FIELD2D& _x, FIELD2D& b, unsigned char* obs) {
//residual
// for(int y=1;y<yRes-1;y++) {
//   for(int x=1;x<xRes-1;x++) {
//     int index = x + y*xRes;
//     //std::cout<<residual.getxRes()<<" a "<<residual.getyRes()<<std::endl;
//     residual[index] = b[index] - (4*_x[index] - _x[index - 1] - _x[index + 1] - _x[index + xRes] - _x[index - xRes]);
//     direction[index] = residual[index]; //d  = r
//    }
// }
  MultPossisonStencil(_x, residual);
  int index = xRes + 1;
  for(int y=1;y<yRes-1;y++, index += 2) {
    for(int x=1;x<xRes-1;x++, index ++) {
      residual[index] = b[index] - residual[index];
      direction[index] = residual[index]; //d  = r
    }
  }
  zeroBoundary(residual,obs);
  //direction = residual;

  Real DeltaNew = 0.f;
  //deltaNew = r^t*r
  DeltaNew = residual.dot(residual);
  Real delta0=DeltaNew;
  //While deltaNew ? (eps^2)*delta0
  const Real eps=_solverEps;
  Real maxR = 2.0f*eps;
  int i=0;
  while ((i<Maxiter) && (DeltaNew>eps)) {
//  for(int y=1;y<yRes-1;y++) {
//    for(int x=1;x<xRes-1;x++) {
//       int index = x + y*xRes;
//      _q[index] = 4*direction[index] - direction[index - 1] - direction[index + 1] - direction[index + xRes] - direction[index - xRes];
//    }
//  }
    MultPossisonStencil(direction, _q);
    zeroBoundary(_q,obs);
    float alpha=0.0f;
    alpha = direction.dot(_q);

    if(fabs(alpha)>0.0f)
      alpha=DeltaNew/alpha;
    //x=x+alpha*direction
    _x.axpy(alpha,direction);
    //r=r-alpha*q r+= -alpha*q
    residual.axpy(-alpha, _q);

    //deltaOld=deltaNew
    float deltaOld=DeltaNew;
    //deltaNew = transpose(r)*r
    DeltaNew=0.0f;
    DeltaNew = residual.dot(residual);
    //beta=deltaNew/deltaOld
    float beta=DeltaNew/deltaOld;

    direction *= beta;
    direction += residual;

    i++;
  }
  // std::cout<<i<<std::endl;
}

void FLUID2D::addSmokeTest(const int xpos, const int ypos, const int width, const int height) {
  addField(density, xpos, ypos, width, height);
}

void FLUID2D::addField(FIELD2D& field, const int x, const int y, const int width, const int height) {
  int startx = x - width / 2;
  int starty = y - height / 2;
  int endx = x + width / 2;
  int endy = y + height / 2;
  
  if (startx < 0) startx = 0;
  if (startx >= xRes) startx = xRes - 1;
  if (starty < 0) starty = 0;
  if (starty >= yRes) starty = yRes - 1;
  if (endx < 0) endx = 0;
  if (endx >= xRes) endx = xRes - 1;
  if (endy < 0) endy = 0;
  if (endy >= yRes) endy = yRes - 1;
  
  for(int j = starty; j <= endy; j++) {
    for (int i = startx; i <= endx; i ++) {
      int index = i + j*xRes;
      field[index] = 1.f;
    }
  }
}

void FLUID2D::setObstaclBoundary() {
  
  int index = xRes + 1;
  int x,y;
  for(y = 1; y<yRes - 1;y++, index += 2) {
    for(x = 1; x<xRes - 1;x++, index++) {
      if(obstacle[index] != EMPTY) {
        const int up    = obstacle[index + xRes];
        const int down  = obstacle[index - xRes];
        const int left  = obstacle[index - 1];
        const int right = obstacle[index + 1];

        velocity[index] = 0.;
        pressure[index] = 0.;

        Real pcnt = 0.;
        if (left && !right) {
          pressure[index] += pressure[index + 1];
          pcnt += 1.;
        }
        if (!left && right) {
          pressure[index] += pressure[index - 1];
          pcnt += 1.;
        }
        if (up && !down) {
          pressure[index] += pressure[index - xRes];
          pcnt += 1.;
        }
        if (!up && down) {
          pressure[index] += pressure[index + xRes];
          pcnt += 1.;
        }
        pressure[index] /= pcnt; 
      }
    }
  }
}

void FLUID2D::zeroBoundary(FIELD2D& x, unsigned char* obs) {
  #pragma omp parallel for
  for(int i=0;i<total;i++)
    x[i] = obs[i] ? 0.0 : x[i];
}

void FLUID2D::MultPossisonStencil(FIELD2D& in, FIELD2D& out) {
  //int index = xRes + 1;
  #pragma omp parallel for
  for(int y=1;y<yRes-1;y++) {
    for(int x=1;x<xRes-1;x++) {
      int index = x + y*xRes;
      if(!obstacle[index]) {
        Real diag = 0;
        Real dleft = 0;
        Real drigh = 0;
        Real dupp  = 0;
        Real ddown = 0;
        if(!obstacle[index-1]) {diag++;dleft++;}
        if(!obstacle[index+1]) {diag++;drigh++;}
        if(!obstacle[index+xRes]) {diag++;dupp++;}
        if(!obstacle[index-xRes]) {diag++;ddown++;}
        //Real result = diag*in[index];
        out[index] = diag*in[index] - dleft*in[index - 1] - drigh*in[index + 1] - dupp*in[index + xRes] - ddown*in[index - xRes];
      }
      else out[index] = 0.;
    }
  }
}

void FLUID2D::addVoriticity() {
  //calculate voriticity
  int index, x, y;
  index = xRes + 1;
  Real gridSize = 0.5f / _dx;
  for (y = 1; y < yRes - 1; y++, index += 2) {	
    for (x = 1; x < xRes - 1; x++, index++) {
      int up    = obstacle[index + xRes] || y == yRes - 2? index : index + xRes;
      int down  = obstacle[index - xRes] || y == 1? index : index - xRes;
      Real dy  = (up == index || down == index) ? 1.0f / _dx : gridSize;
      int right = obstacle[index + 1] || x == xRes - 2 ? index : index + 1;
      int left  = obstacle[index - 1] || x == 1 ? index : index - 1;
      Real dx  = (left == index || right == index) ? 1.0f / _dx : gridSize;

      voriticity[index] = (velocity[right][1] - velocity[left][1]) * dx + (- velocity[up][0] + velocity[down][0])* dy;
    }
  }
  index = xRes + 1;
  for (y = 1; y < yRes - 1; y++, index += 2) {	
    for (x = 1; x < xRes - 1; x++, index++) {
      Real N[2];
      int up    = obstacle[index + xRes] || y == yRes - 2? index : index + xRes;
      int down  = obstacle[index - xRes] || y == 1? index : index - xRes;
      Real dy  = (up == index || down == index) ? 1.0f / _dx : gridSize;
      int right = obstacle[index + 1] || x == xRes - 2 ? index : index + 1;
      int left  = obstacle[index - 1] || x == 1 ? index : index - 1;
      Real dx  = (left == index || right == index) ? 1.0f / _dx : gridSize;
      N[0] = (voriticity[right] - voriticity[left]) * dx; 
      N[1] = (voriticity[up] - voriticity[down]) * dy;

      Real Magnitude = sqrt(N[0]*N[0] + N[1]*N[1]);
      if (Magnitude > 0.0) {
        Magnitude = 1.0 / Magnitude;
        N[0] *= Magnitude;
        N[1] *= Magnitude;

        force[index][0] += (N[1] * voriticity[index] - 0.) * _dx * _vorticityEps;
        force[index][1] -= (N[0] * voriticity[index] - 0.) * _dx * _vorticityEps;
      }
    }
  }
}

void FLUID2D::ClearDens() {
  density.clear();
  //std::cout<<"sa "<<std::endl;
}

void FLUID2D::ClearVelo() {
  velocity.clear();
}

void FLUID2D::drawDensity() {
  float x, y, h, d00, d01, d10, d11;
  h = _dx;
  glPushMatrix();
  // glTranslatef(-1.0,0.0,0);
  glBegin ( GL_QUADS );

  for (int j=0 ; j<yRes-1; j++ ) {
    y = (j-0.5f)*h;
    for (int i=0 ; i<xRes-1; i++ ) {
      x = (i-0.5f)*h;

      int index = i + j*xRes;
      d00 = 1.0 - density[index];
      d01 = 1.0 - density[index + xRes];
      d10 = 1.0 - density[index + 1];
      d11 = 1.0 - density[index + 1 + xRes];

      glColor3f ( d00, d00, d00 );
      glVertex2f ( x, y );
      glColor3f ( d10, d10, d10 ); glVertex2f ( x+h, y );
      glColor3f ( d11, d11, d11 ); glVertex2f ( x+h, y+h );
      glColor3f ( d01, d01, d01 ); glVertex2f ( x, y+h );
    }
  }
 
  
  glEnd ();
  glPopMatrix();
}

void FLUID2D::DrawObstal() {
  float x, y, h, d00, d01, d10, d11;
  h = _dx;

  glBegin ( GL_QUADS );

  for (int j = 0 ; j < yRes; j++ ) {
    y = j*h;
    for (int i = 0 ; i < xRes; i++ ) {
      x = i*h;

      int index = i + j*xRes;
      if (obstacle[index]) {
        glColor3f ( 1.0, 0, 0 );
        glVertex2f ( x, y );
        glVertex2f ( x+h, y );
        glVertex2f ( x+h, y+h );
        glVertex2f ( x, y+h );
      }
    }
  }
  glEnd ();
}

void FLUID2D::drawVelocity() {
  int i, j;
  float x, y, h;

  h = _dx;

  glColor3f ( 1.0f, 1.0f, 1.0f );
  glLineWidth ( 1.0f );

  glBegin ( GL_LINES );

  for (int j=0 ; j<yRes-1; j++ ) {
  y = (j-0.5f)*h;
   for (int i=0 ; i<xRes-1; i++ ) {
      x = (i-0.5f)*h;
      int index = i + j*xRes;
      if(obstacle[index] == EMPTY) {
        glVertex2f ( x, y );
        glVertex2f ( x+velocity[index][0], y+velocity[index][1]);
      }
    }
  }
  glEnd ();
}

const float AMPLIFACATION = 2.5;

void FLUID2D::drawVelocityX() {
  float x, y, h, d00, d01, d10, d11;
  h = _dx;

  glBegin ( GL_QUADS );

  for (int j=0 ; j<yRes-1; j++ ) {
    y = (j-0.5f)*h;
    for (int i=0 ; i<xRes-1; i++ ) {
      x = (i-0.5f)*h;

      int index = i + j*xRes;
      d00 = velocity[index][0]*AMPLIFACATION;
      d01 = velocity[index + xRes][0]*AMPLIFACATION;
      d10 = velocity[index + 1][0]*AMPLIFACATION;
      d11 = velocity[index + 1 + xRes][0]*AMPLIFACATION;

      glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
      glColor3f ( d10, d10, d10 ); glVertex2f ( x+h, y );
      glColor3f ( d11, d11, d11 ); glVertex2f ( x+h, y+h );
      glColor3f ( d01, d01, d01 ); glVertex2f ( x, y+h );
    }
  }
  glEnd ();
}

void FLUID2D::drawVelocityY() {
  float x, y, h, d00, d01, d10, d11;
  h = _dx;

  glBegin ( GL_QUADS );

  for (int j=0 ; j<yRes-1; j++ ) {
    y = (j-0.5f)*h;
    for (int i=0 ; i<xRes-1; i++ ) {
      x = (i-0.5f)*h;

      int index = i + j*xRes;
      d00 = velocity[index][1]*AMPLIFACATION;
      d01 = velocity[index + xRes][1]*AMPLIFACATION;
      d10 = velocity[index + 1][1]*AMPLIFACATION;
      d11 = velocity[index + 1 + xRes][1]*AMPLIFACATION;

      glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
      glColor3f ( d10, d10, d10 ); glVertex2f ( x+h, y );
      glColor3f ( d11, d11, d11 ); glVertex2f ( x+h, y+h );
      glColor3f ( d01, d01, d01 ); glVertex2f ( x, y+h );
    }
  }
  glEnd ();
}

void FLUID2D::InitCheckBoardPattern() {
  // Init a checkboard pattern smoke.
  double checkSize = 0.1;
  double delta = 0.15;
  
  int iCheckSize = static_cast<int>(checkSize*xRes);
  int iDelta = static_cast<int>(delta*xRes);
  
  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      int ix = (x + iDelta) / iCheckSize;
      int iy = (y + iDelta) / iCheckSize;
      const int index = x + y*xRes;
      if (ix % 2 == 0) {
        // Fill the block with smoke.
        density[index] = 1.0;
      } else {
        density[index] = 0.0;
      }
    }
  }
}

void FLUID2D::AdvectParticles() {
   const Real dt0 = _dt / _dx;
  for (int i = 0; i < num_particles_; i++) {
    // Get the velocity at the particles.
    const VEC2& position = particles_[i].position;
    VEC2 p_v = velocity.GetVelocity(position[0], position[1]);
    // Forward Eular.
    particles_[i].position[0] += p_v[0] * dt0;
    particles_[i].position[1] += p_v[1] * dt0;
    particles_[i].velocity = p_v;
    if (particles_[i].position[0] < 0. || particles_[i].position[0] > xRes || 
      particles_[i].position[1] < 0. || particles_[i].position[1] > yRes ) {
      particles_[i].position[0] = std::rand() % xRes + 0.5;
      particles_[i].position[1] = std::rand() % yRes + 0.5;
      particles_[i].velocity = 0.;
    }
  }
}

void FLUID2D::ReSeedParticles() {
  // Reseed the particles at random position.
  for (int i = 0; i < num_particles_; i++) {
    int px = std::rand() % xRes;
    int py = std::rand() % yRes;
    particles_[i].position[0] = px + 0.5;
    particles_[i].position[1] = py + 0.5;
    particles_[i].velocity = 0.;
  }
}

void FLUID2D::MouseAddForce(int xpos, int ypos, const float fx, const float fy) {
  std::cout <<  xpos << " sa " << ypos << std::endl;
}

void FLUID2D::DrawParticles(const double ptl_length) {
  drawer_.get()->DrawParticles(particles_, _dx, _dt, ptl_length);
}
