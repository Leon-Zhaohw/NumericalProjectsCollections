#include "fluid.h"

static const double FEPS = 1e-6;

inline double trilinearInterp(const Array3D<double> &x, int i, int j, int k,
    double w0, double w1, double w2) {
  return (1.0-w0)*(1.0-w1)*(1.0-w2)*x(i,j,k)+
    (1.0-w0)*(1.0-w1)*(w2)*x(i,j,k+1)+
    (1.0-w0)*(w1)*(1.0-w2)*x(i,j+1,k)+
    (1.0-w0)*(w1)*(w2)*x(i,j+1,k+1)+
    (w0)*(1.0-w1)*(1.0-w2)*x(i+1,j,k)+
    (w0)*(1.0-w1)*(w2)*x(i+1,j,k+1)+
    (w0)*(w1)*(1.0-w2)*x(i+1,j+1,k)+
    (w0)*(w1)*(w2)*x(i+1,j+1,k+1);
}

inline void floor(const Eigen::Vector3d &x, double h, int &i, int &j, int &k) {
  i = (int)floor(x[0] / h);
  j = (int)floor(x[1] / h);
  k = (int)floor(x[2] / h);
}

inline void weights(const Eigen::Vector3d &x, double h, int &i, int &j, int &k, double &w0, double &w1, double &w2) {
  floor(x, h, i, j, k);
  w0 = (x[0] - i*h)/h;
  w1 = (x[1] - j*h)/h;
  w2 = (x[2] - k*h)/h;
}

double StaggeredGrid::interpxface(const Eigen::Vector3d &x, const Array3D<double> &f) const {
  Eigen::Vector3d y(x-lc);
  int i,j,k;
  double w0, w1, w2;
  y[1] -= halfh;
  y[2] -= halfh;
  weights(y, h, i, j, k, w0, w1, w2);

  return trilinearInterp(f, i, j, k, w0, w1, w2);
}

double StaggeredGrid::interpyface(const Eigen::Vector3d &x, const Array3D<double> &f) const {
  Eigen::Vector3d y(x-lc);
  int i,j,k;
  double w0, w1, w2;
  y[0] -= halfh;
  y[2] -= halfh;
  weights(y, h, i, j, k, w0, w1, w2);
  
  return trilinearInterp(f, i, j, k, w0, w1, w2);
}

double StaggeredGrid::interpzface(const Eigen::Vector3d &x, const Array3D<double> &f) const {
  Eigen::Vector3d y(x-lc);
  int i,j,k;
  double w0, w1, w2;
  y[0] -= halfh;
  y[1] -= halfh;
  weights(y, h, i, j, k, w0, w1, w2);
  
  return trilinearInterp(f, i, j, k, w0, w1, w2);
}

Eigen::Vector3d StaggeredGrid::interpVel(const Eigen::Vector3d &x) const {
  return Eigen::Vector3d(interpxface(x,u), interpyface(x,v), interpzface(x,w));
}

void StaggeredGrid::advectParticles(double dt, std::vector<Particle> &particles) {
  for (std::vector<Particle>::iterator p = particles.begin(); p!=particles.end(); p++) {
    p->pos += dt*interpVel(p->pos);
    clipToGrid(p->pos);
  }
}

void StaggeredGrid::gridToParticles(double flip_ratio, std::vector<Particle> &particles) {
  for (std::vector<Particle>::iterator p = particles.begin(); p!=particles.end(); p++) {
    Eigen::Vector3d oVel(interpxface(p->pos,fu), interpyface(p->pos,fv), interpzface(p->pos,fw));
    Eigen::Vector3d nVel(interpVel(p->pos));
    p->vel = flip_ratio * (p->vel-oVel) + nVel;
  }
}

inline void contribute(double w, double v, double &u, double &fu) {
  fu += w;
  u += w*v;
}

void StaggeredGrid::splatParticleValue(Array3D<double> &u, Array3D<double> &fu,
    double w0, double w1, double w2, int i, int j, int k, double v) {
  contribute ((1.0-w0)*(1.0-w1)*(1.0-w2), v, u(i,j,k), fu(i,j,k));
  contribute ((1.0-w0)*(1.0-w1)*(w2), v, u(i,j,k+1), fu(i,j,k+1));
  contribute ((1.0-w0)*(w1)*(1.0-w2), v, u(i,j+1,k), fu(i,j+1,k));
  contribute ((1.0-w0)*(w1)*(w2), v, u(i,j+1,k+1), fu(i,j+1,k+1));
  contribute ((w0)*(1.0-w1)*(1.0-w2), v, u(i+1,j,k), fu(i+1,j,k));
  contribute ((w0)*(1.0-w1)*(w2), v, u(i+1,j,k+1), fu(i+1,j,k+1));
  contribute ((w0)*(w1)*(1.0-w2), v, u(i+1,j+1,k), fu(i+1,j+1,k));
  contribute ((w0)*(w1)*(w2), v, u(i+1,j+1,k+1), fu(i+1,j+1,k+1));
}

void StaggeredGrid::particlesToGrid(const std::vector<Particle> &particles) {
  u = 0.0;
  fu = 0.0;
  v = 0.0;
  fv = 0.0;
  w = 0.0;
  fw = 0.0;

  for (unsigned int i=0; i<nx; i++) {
    for (unsigned int j=0; j<ny; j++) {
      for (unsigned int k=0; k<nz; k++) {
	if (boundary(i,j,k)) {
	  cellLabels(i,j,k) = OBSTACLE;
	  continue;
	}
	cellLabels(i,j,k) = AIR;
      }
    }
  }
  
  for (std::vector<Particle>::const_iterator p = particles.begin(); p!=particles.end(); p++) {
    Eigen::Vector3d y(p->pos-lc);
    int i,j,k;
    double w0, w1, w2;

    floor(y, h, i, j, k);
    cellLabels(i,j,k) = LIQUID;
		
    y[1] -= halfh;
    y[2] -= halfh;
    weights(y, h, i, j, k, w0, w1, w2);
    splatParticleValue(u, fu, w0, w1, w2, i, j, k, p->vel[0]);
    
    y[0] -= halfh;
    y[1] += halfh;
    weights(y, h, i, j, k, w0, w1, w2);
    splatParticleValue(v, fv, w0, w1, w2, i, j, k, p->vel[1]);
	
    y[1] -= halfh;
    y[2] += halfh;
    weights(y, h, i, j, k, w0, w1, w2);
    splatParticleValue(w, fw, w0, w1, w2, i, j, k, p->vel[2]);
  }
  
  for (unsigned int i=0; i<nx+1; i++) {
    for (unsigned int j=0; j<ny; j++) {
      for (unsigned int k=0; k<nz; k++) {
	if (i == 0 || i == 1 || i == nx || i == nx-1 || fu(i,j,k) < FEPS) {
	  u(i,j,k) = 0.0;
	} else {
	  u(i,j,k) /= fu(i,j,k);
	}
	fu(i,j,k) = u(i,j,k);
      }
    }
  }
  for (unsigned int i=0; i<nx; i++) {
    for (unsigned int j=0; j<ny+1; j++) {
      for (unsigned int k=0; k<nz; k++) {
	if (j == 0 || j == 1 || j == ny || j == ny-1 || fv(i,j,k) < FEPS) {
	  v(i,j,k) = 0.0;
	} else {
	  v(i,j,k) /= fv(i,j,k);
	}
	fv(i,j,k) = v(i,j,k);
      }
    }
  }
  for (unsigned int i=0; i<nx; i++) {
    for (unsigned int j=0; j<ny; j++) {
      for (unsigned int k=0; k<nz+1; k++) {
	if (k == 0 || k == 1 || k == nz || k == nz-1 || fw(i,j,k) < FEPS) {
	  w(i,j,k) = 0.0;
	} else {
	  w(i,j,k) /= fw(i,j,k);
	}
	fw(i,j,k) = w(i,j,k);
      }
    }
  }
  
  setBoundaryVelocity();
}

void StaggeredGrid::computePressure() {
  double sigma_new, sigma, alpha, beta;
  double tol=1e-6;
  int iter=0, max_iter=1000;
  int n=nx*ny*nz;
  
  // Cache neighbor data to make applying the Laplacian easier/faster.
  laplacian = 0;
  for (unsigned int i=1; i<nx-1; i++) {
    for (unsigned int j=1; j<ny-1; j++) {
      for (unsigned int k=1; k<nz-1; k++) {
	if (cellLabels(i,j,k) != LIQUID) continue;
	if (cellLabels(i-1,j,k) & (~OBSTACLE)) {
	  laplacian(i,j,k) += 1;
	  if (cellLabels(i-1,j,k) & LIQUID) laplacian(i,j,k) |= 8;
	}
	if (cellLabels(i,j-1,k) & (~OBSTACLE)) {
	  laplacian(i,j,k) += 1;
	    if (cellLabels(i,j-1,k) & LIQUID) laplacian(i,j,k) |= 16;
	}
	if (cellLabels(i,j,k-1) & (~OBSTACLE)) {
	  laplacian(i,j,k) += 1;
	  if (cellLabels(i,j,k-1) & LIQUID) laplacian(i,j,k) |= 32;
	}
	if (cellLabels(i+1,j,k) & (~OBSTACLE)) {
	  laplacian(i,j,k) += 1;
	  if (cellLabels(i+1,j,k) & LIQUID) laplacian(i,j,k) |= 64;
	}
	if (cellLabels(i,j+1,k) & (~OBSTACLE)) {
	  laplacian(i,j,k) += 1;
	    if (cellLabels(i,j+1,k) & LIQUID) laplacian(i,j,k) |= 128;
	}
	if (cellLabels(i,j,k+1) & (~OBSTACLE)) {
	  laplacian(i,j,k) += 1;
	  if (cellLabels(i,j,k+1) & LIQUID) laplacian(i,j,k) |= 256;
	}
      }
    }
  }

  // compute residual (weighted divergence)
  for (unsigned int i=1; i<nx-1; i++) {
    for (unsigned int j=1; j<ny-1; j++) {
      for (unsigned int k=1; k<nz-1; k++) {
	if (cellLabels(i,j,k) & (AIR | OBSTACLE)) {
	  r(i,j,k) = 0.0;
	} else {
	  r(i,j,k) = -(u(i+1,j,k) - u(i,j,k) + v(i,j+1,k) - v(i,j,k) +
	      w(i,j,k+1) - w(i,j,k));
	}
      }
    }
  }

  // conjugate gradients
  p = 0.0;
  d = r;
  
  sigma = dot(r,r);
  tol *= sigma;
  
  while (iter++ < max_iter && sigma > tol) {
    applyA(d,q);
    alpha = sigma / dot(d,q);
    daxpy(alpha, d, p);
    daxpy(-alpha, q, r);
    double sigma_old = sigma;
    sigma = dot(r,r);
    beta = sigma / sigma_old;
    q = r;
    daxpy(beta, d, q);
    d = q;
  }

  // subtract pressure gradient
  for (unsigned int i=1; i<nx-1; i++) {
    for (unsigned int j=1; j<ny-1; j++) {
      for (unsigned int k=1; k<nz-1; k++) {
	if (cellLabels(i,j,k) & (OBSTACLE)) continue;
	if (cellLabels(i-1,j,k) & (~OBSTACLE)) u(i,j,k) -= (p(i,j,k) - p(i-1,j,k));
	if (cellLabels(i,j-1,k) & (~OBSTACLE)) v(i,j,k) -= (p(i,j,k) - p(i,j-1,k));
	if (cellLabels(i,j,k-1) & (~OBSTACLE)) w(i,j,k) -= (p(i,j,k) - p(i,j,k-1));
      }
    }
  }
}

void StaggeredGrid::applyA(Array3D<double> &s, Array3D<double> &z) {
  for (unsigned int i=1; i<nx-1; i++) {
    for (unsigned int j=1; j<ny-1; j++) {
      for (unsigned int k=1; k<nz-1; k++) {
	unsigned short &l = laplacian(i,j,k);
	if (!l) {z(i,j,k) = 0.0; continue;}
	z(i,j,k) = (l & 7)*s(i,j,k) -
	  ((l & 8) ? s(i-1,j,k) : 0) -
	  ((l & 16) ? s(i,j-1,k) : 0) -
	  ((l & 32) ? s(i,j,k-1) : 0) -
	  ((l & 64) ? s(i+1,j,k) : 0) -
	  ((l & 128) ? s(i,j+1,k) : 0) -
	  ((l & 256) ? s(i,j,k+1) : 0);
      }
    }
  }
}

void StaggeredGrid::gravity(double dt) {
  for (unsigned int i=0; i<w.nx(); i++) {
    for (unsigned int j=0; j<w.ny(); j++) {
      for (unsigned int k=0; k<w.nz(); k++) {
	w(i,j,k) -= dt*9.8;
      }
    }
  }
}

inline void StaggeredGrid::clipToGrid(Eigen::Vector3d &x) const {
  static const double FEPS = 1e-4;
  if (x[0] <= lc[0]+h+FEPS) x[0] = lc[0]+h+FEPS;
  if (x[0] >= uc[0]-h-FEPS) x[0] = uc[0]-h-FEPS;
  if (x[1] <= lc[1]+h+FEPS) x[1] = lc[1]+h+FEPS;
  if (x[1] >= uc[1]-h-FEPS) x[1] = uc[1]-h-FEPS;
  if (x[2] <= lc[2]+h+FEPS) x[2] = lc[2]+h+FEPS;
  if (x[2] >= uc[2]-h-FEPS) x[2] = uc[2]-h-FEPS;
};

inline bool StaggeredGrid::boundary(unsigned int i, unsigned int j, unsigned int k) const {
  if (i == 0 || i == nx-1 || j == 0 || j == ny-1 || k == 0 || k == nz-1)
    return true;
  return false;
};

void StaggeredGrid::setBoundaryVelocity() {
  unsigned int i,j,k;

  // set boundary conditions
  for (j=0; j<ny; j++) {
    for (k=0; k<nz; k++) {
      u(0,j,k) = u(1,j,k) = 0.0;
      v(0,j,k) = v(1,j,k);
      w(0,j,k) = w(1,j,k);
      u(nx,j,k) = u(nx-1,j,k) = 0.0;
      v(nx-1,j,k) = v(nx-2,j,k);
      w(nx-1,j,k) = w(nx-2,j,k);
    }
  }
  for (i=0; i<nx; i++) {
    for (k=0; k<nz; k++) {
      u(i,0,k) = u(i,1,k);
      v(i,0,k) = v(i,1,k) = 0.0;
      w(i,0,k) = w(i,1,k);
      u(i,ny-1,k) = u(i,ny-2,k);
      v(i,ny,k) = v(i,ny-1,k) = 0.0;
      w(i,ny-1,k) = w(i,ny-2,k);
    }
  }
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      u(i,j,0) = u(i,j,1);
      v(i,j,0) = v(i,j,1);
      w(i,j,0) = w(i,j,1) = 0.0;
      u(i,j,nz-1) = u(i,j,nz-2);
      v(i,j,nz-1) = v(i,j,nz-2);
      w(i,j,nz) = w(i,j,nz-1) = 0.0;
    }
  }
}

StaggeredGrid::StaggeredGrid(const SimulationParameters &params){
  nx = params.nx;
  ny = params.ny;
  nz = params.nz;
  nynz = ny*nz;
  h = params.h;
  halfh = h/2.0;
  lc = params.lc;
  uc = params.uc;
  u.allocate(nx+1, ny, nz);
  v.allocate(nx, ny+1, nz);
  w.allocate(nx, ny, nz+1);
  p.allocate(nx, ny, nz);
  nu.allocate(nx+1, ny, nz);
  nv.allocate(nx, ny+1, nz);
  nw.allocate(nx, ny, nz+1);
  fu.allocate(nx+1, ny, nz);
  fv.allocate(nx, ny+1, nz);
  fw.allocate(nx, ny, nz+1);
  cellLabels.allocate(nx, ny, nz);
  laplacian.allocate(nx, ny, nz);
  r.allocate(nx, ny, nz);
  d.allocate(nx, ny, nz);
  q.allocate(nx, ny, nz);
}



int main(int argc, char *argv[]) {
  char fname[80];
  SimulationParameters params;
  double time = 0;
  int frame = 0;
  double frame_time = -1.0;
  std::vector<Particle> particles;

  readInputFile(argv[1], params, particles);

  StaggeredGrid grid(params);

  grid.particlesToGrid(particles);

  while (time < params.total_time) {
    if (frame_time < 0) {
      sprintf(fname, params.output_fname.c_str(), frame);
      writeParticles(fname, particles);
      frame_time = 1.0/30.0-0.0001;
      frame++;
    }
    
    grid.advectParticles(params.dt, particles);
    grid.particlesToGrid(particles);

    grid.gravity(params.dt);
    grid.setBoundaryVelocity();

    grid.computePressure();
    
    grid.gridToParticles(params.flip_ratio, particles);
    
    time += params.dt;
    frame_time -= params.dt;
  }
}















/////////////////////////////////////////
// I/O
/////////////////////////////////////////
#include "json/json.h"
#include <fstream>

void readParticles(const char *fname, std::vector<Particle> &particles) {
  std::ifstream in(fname, std::ios::in);
  int nparts;
  in>>nparts;
  particles.resize(nparts);
  for (std::vector<Particle>::iterator p = particles.begin(); p != particles.end(); p++) {
    in>>p->pos[0]>>p->pos[1]>>p->pos[2]>>p->vel[0]>>p->vel[1]>>p->vel[2];
  }
  in.close();
  std::cout<<"inputfile " << fname <<" read"<<std::endl;
}

void writeParticles(const char *fname, const std::vector<Particle> &particles) {
  std::ofstream out(fname, std::ios::out);
  out<<particles.size()<<std::endl;
  for (std::vector<Particle>::const_iterator p = particles.begin(); p != particles.end(); p++) {
    out<<p->pos[0]<<" "<<p->pos[1]<<" "<<p->pos[2]<<" "<<p->vel[0]<<" "<<p->vel[1]<<" "<<p->vel[2]<<" "<<std::endl;
  }
  out.close();
  std::cout<<"outputfile " << fname <<" written"<<std::endl;
}

void readInputFile(const char *fname, 
    SimulationParameters &params,
    std::vector<Particle> &particles) {
  
  std::ifstream in(fname, std::ios::in);

  Json::Value root;
  Json::Reader jReader;

  if(!jReader.parse(in, root)){
	std::cout << "couldn't read input file: " << fname << '\n'
		  << jReader.getFormattedErrorMessages() << std::endl;
	exit(1);
  }

  params.dt = root.get("dt", 1.0/300.0).asDouble();
  params.total_time = root.get("total_time", 1.0).asDouble();
  params.density = root.get("density", 1.0/300.0).asDouble();
  params.flip_ratio = root.get("flip_ratio", 0.95).asDouble();

  params.nx = root["res"][0].asInt();
  params.ny = root["res"][1].asInt();
  params.nz = root["res"][2].asInt();
  params.h = root["h"].asDouble();
  params.lc << root["lc"][0].asDouble(), root["lc"][1].asDouble(), root["lc"][2].asDouble();
  params.uc = params.lc + Eigen::Vector3d(params.h*params.nx, params.h*params.ny, params.h*params.nz);

  params.output_fname = root.get("output_fname", std::string("output.%04d.obj")).asString();

  readParticles(root["particles"].asString().c_str(), particles);
}

