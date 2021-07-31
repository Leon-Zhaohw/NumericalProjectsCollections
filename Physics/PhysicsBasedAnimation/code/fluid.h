
#ifndef FLUID_H
#define FLUID_H

#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "array3D.h"

struct SimulationParameters {
  double dt, total_time, density, flip_ratio, h;
  int nx, ny, nz;
  Eigen::Vector3d lc, uc;
  std::string output_fname;
};

struct Particle {
  Eigen::Vector3d pos, vel;
};

class StaggeredGrid {
public:
  StaggeredGrid(const SimulationParameters &params);
  
  Eigen::Vector3d interpVel(const Eigen::Vector3d &x) const;
  double interpxface(const Eigen::Vector3d &x, const Array3D<double> &f) const;
  double interpyface(const Eigen::Vector3d &x, const Array3D<double> &f) const;
  double interpzface(const Eigen::Vector3d &x, const Array3D<double> &f) const;
  
  void advectParticles(double dt, std::vector<Particle> &particles);
  void particlesToGrid(const std::vector<Particle> &particles);
  void gridToParticles(double flip_ration, std::vector<Particle> &particles);

  void setBoundaryVelocity();
  void computePressure();
  
  void gravity(double dt);

  static const unsigned char AIR = 1;
  static const unsigned char LIQUID = 2;
  static const unsigned char OBSTACLE = 4;
  
private:
  Array3D<double> u,v,w;
  Array3D<unsigned char> cellLabels;
  unsigned int nx, ny, nz, nynz;
  double h,halfh;
  Eigen::Vector3d lc, uc;
  
  void splatParticleValue(Array3D<double> &u, Array3D<double> &fu,
      double w0, double w1, double w2,
      int i, int j, int k, double v);
  void applyA(Array3D<double> &x, Array3D<double> &b);
  inline void clipToGrid(Eigen::Vector3d &x) const;
  inline bool boundary(unsigned int i, unsigned int j, unsigned int k) const;
  Array3D<double> nu,nv,nw,fu,fv,fw;
  Array3D<double> p, r, d, q;
  Array3D<unsigned short> laplacian;
};


// I/O

void readParticles(const char *fname, std::vector<Particle> &particles);
void writeParticles(const char *fname, const std::vector<Particle> &particles);
void readInputFile(const char *fname, SimulationParameters &params, std::vector<Particle> &particles);

  
#endif
