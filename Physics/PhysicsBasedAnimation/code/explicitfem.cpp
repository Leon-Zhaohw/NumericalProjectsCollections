#include "explicitfem.h"

static const double FEPS = 1e-6;

void polar(const Eigen::Matrix3d &F, Eigen::Matrix3d &Q) {
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();
    Q = U * V.transpose();
}

void initialize(Object &obj, double density) {
  for (std::vector<Particle>::iterator p = obj.particles.begin(); p != obj.particles.end(); p++) {
    p->mass = 0.0;
  }

  for (std::vector<Element>::iterator e = obj.elements.begin(); e != obj.elements.end(); e++) {
    Eigen::Vector3d &x0 = obj.particles[(*e)[0]].pos;
    Eigen::Vector3d &x1 = obj.particles[(*e)[1]].pos;
    Eigen::Vector3d &x2 = obj.particles[(*e)[2]].pos;
    Eigen::Vector3d &x3 = obj.particles[(*e)[3]].pos;

    e->normals[0] = (x2-x1).cross(x3-x1);
    e->normals[1] = (x3-x0).cross(x2-x0);
    e->normals[2] = (x1-x0).cross(x3-x0);
    e->normals[3] = (x2-x0).cross(x1-x0);

    e->basis << x1-x0, x2-x0, x3-x0;
    double det = e->basis.determinant();
    double mass = params.density * det / 24.0;
    e->basis = e->basis.inverse().eval();

    obj.particles[(*e)[0]].mass += mass;
    obj.particles[(*e)[1]].mass += mass;
    obj.particles[(*e)[2]].mass += mass;
    obj.particles[(*e)[3]].mass += mass;
  }
}

int main(int argc, char *argv[]) {
  char fname[80];
  SimulationParameters params;
  std::vector<Object> objects;
  readInputFile(argv[1], params, objects);

  double time = 0;
  int frame = 0;
  double frameTime = -1.0;

  for (std::vector<Object>::iterator obj = objects.begin(); obj != objects.end(); obj++) {
    initialize(*obj, params.density);
  }

  while (time < params.total_time) {
    // write frame
    if (frameTime < 0) {
      for (unsigned int o=0; o<objects.size(); o++) {
	sprintf(fname, params.output_fname.c_str(), o, frame);
	writeObj(fname, objects[o].particles, objects[o].triangles);
	frameTime = 1.0/30.0-0.0001;
	frame++;
      }
    }

    // takes a timestep
    for (std::vector<Object>::iterator obj = objects.begin(); obj != objects.end(); obj++) {
      // initialize force accumulator, apply gravity
      for (std::vector<Particle>::iterator p = obj->particles.begin(); p != obj->particles.end(); p++) {
	p->frc = Eigen::Vector3d::Zero();
	p->frc[2] = -9.8*p->mass;
      }

      // comupte and accumulate spring forces
      for (std::vector<Element>::const_iterator e = obj->elements.begin(); e != obj->elements.end(); e++) {
	Eigen::Matrix3d X, F, Q, Ftilde, strain, stress;
	Eigen::Matrix3d Xdot, Fdot, Fdottilde, strainrate, stressrate;

	Eigen::Vector3d &x0 = obj->particles[(*e)[0]].pos;
	Eigen::Vector3d &x1 = obj->particles[(*e)[1]].pos;
	Eigen::Vector3d &x2 = obj->particles[(*e)[2]].pos;
	Eigen::Vector3d &x3 = obj->particles[(*e)[3]].pos;

	Eigen::Vector3d &v0 = obj->particles[(*e)[0]].vel;
	Eigen::Vector3d &v1 = obj->particles[(*e)[1]].vel;
	Eigen::Vector3d &v2 = obj->particles[(*e)[2]].vel;
	Eigen::Vector3d &v3 = obj->particles[(*e)[3]].vel;

	X << x1-x0, x2-x0, x3-x0;
	F = X * e->basis;
	polar(F, Q);
	Ftilde = Q.transpose() * F;
	strain = 0.5 * (Ftilde + Ftilde.transpose()) - Eigen::Matrix3d::Identity();
	stress = params.lambda * strain.trace() * Eigen::Matrix3d::Identity() + 2 * params.mu * strain;
	
	Xdot << v1-v0, v2-v0, v3-v0;
	Fdot = Xdot * e->basis;
	Fdottilde = Q.transpose() * Fdot;
	strainrate = 0.5 * (Ftildedot + Ftildedot.transpose()) - Eigen::Matrix3d::Identity();
	stressrate = params.damp*(params.lambda * strainrate.trace() * Eigen::Matrix3d::Identity() + 2 * params.mu * strainrate);
	
	obj->particles[(*e)[0]].frc += Q * (stress+stressrate) * e->normals[0] / 6.0;
	obj->particles[(*e)[1]].frc += Q * (stress+stressrate) * e->normals[1] / 6.0;
	obj->particles[(*e)[2]].frc += Q * (stress+stressrate) * e->normals[2] / 6.0;
	obj->particles[(*e)[3]].frc += Q * (stress+stressrate) * e->normals[3] / 6.0;
      }

      for (std::vector<Particle>::iterator p = obj->particles.begin(); p != obj->particles.end(); p++) {
	// time integration
	p->vel += params.dt*(p->frc / p->mass);
	p->pos += params.dt*(p->vel);

	// ground collision
	if (p->pos[2] < 0.0) {
	  p->pos[2] = 0.0;
	  p->vel[2] = 0.0;
	}
      }
    }
    
    time += params.dt;
    frameTime -= params.dt;
  }
}




















/////////////////////////////////////////
// I/O
/////////////////////////////////////////

#include "json/json.h"
#include <fstream>

bool readInputFile(const char *fname, 
    SimulationParameters &params,
    std::vector<Object> &objects) {
  
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
  params.density = root.get("density", 1e4).asDouble();
  params.lambda = root.get("lambda", 1e4).asDouble();
  params.mu = root.get("mu", 1e4).asDouble();
  params.damp = root.get("damping", 1e-2).asDouble();
  params.output_fname = root.get("output_fname", std::string("output-%02d.%04d.obj")).asString();
  
  Json::Value objectsIn = root["objects"];
  objects.resize(objectsIn.size());
  for (unsigned int i=0; i<objectsIn.size(); i++) {
    readObject((objectsIn[i]["filename"]).asString().c_str(), objects[i]);
  }
  return true;
}

bool readObject(const char *fname, Object &object) {
  char ch;
  Particle p;
  Element e;
  Tri t;
  

  std::ifstream in(fname, std::ios::in);
  while (in>>ch) {
    if (ch == 'p') {
      in>>p.pos[0]>>p.pos[1]>>p.pos[2]>>p.vel[0]>>p.vel[1]>>p.vel[2];
      object.particles.push_back(p);
      continue;
    }
    if (ch == 'e') {
      in>>e[0]>>e[1]>>e[2]>>e[3];
      object.elements.push_back(e);
      continue;
    }
    if (ch == 't') {
      in>>t[0]>>t[1]>>t[2];
      object.triangles.push_back(t);
      continue;
    }
  }
  in.close();
  std::cout<<"inputfile " << fname <<" read"<<std::endl;
  return true;
}

void writeObj(char *fname, const std::vector<Particle> &meshPts, const std::vector<Tri> &triangles) {
  std::cout<<"writing "<<fname<<std::endl;
  std::ofstream out;
  std::vector<Particle>::const_iterator p;
  std::vector<Tri>::const_iterator t;
  
  out.open(fname);
  
  for (p=meshPts.begin(); p!=meshPts.end(); p++) 
    out<<"v "<<p->pos[0]<<" "<<p->pos[1]<<" "<<p->pos[2]<<std::endl;
  
  for (t=triangles.begin(); t!=triangles.end(); t++) 
    out<<"f "<<(*t)[0]+1<<" "<<(*t)[1]+1<<" "<<(*t)[2]+1<<std::endl;
  
  out.close();
}


