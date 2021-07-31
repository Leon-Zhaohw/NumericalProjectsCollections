#include "springmass.h"

static const double FEPS = 1e-6;

int main(int argc, char *argv[]) {
  char fname[80];
  SimulationParameters params;
  std::vector<Object> objects;
  readInputFile(argv[1], params, objects);

  double time = 0;
  int frame = 0;
  double frameTime = -1.0;

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
      for (std::vector<Spring>::const_iterator s = obj->springs.begin(); s != obj->springs.end(); s++) {
	Eigen::Vector3d d = obj->particles[s->j].pos - obj->particles[s->i].pos;
	double l = d.norm();
	Eigen::Vector3d v = obj->particles[s->j].vel - obj->particles[s->i].vel;
	Eigen::Vector3d frc = (params.k_s*((l / s->r) - 1.0) + params.k_d*(v.dot(d)/(l*s->r))) * (d/l);
	obj->particles[s->i].frc += frc;
	obj->particles[s->j].frc -= frc;
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
  params.k_s = root.get("stiffness", 1.0/300.0).asDouble();
  params.k_d = root.get("damping", 1.0/300.0).asDouble();
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
  Spring s;
  Tri t;
  

  std::ifstream in(fname, std::ios::in);
  while (in>>ch) {
    if (ch == 'p') {
      in>>p.mass>>p.pos[0]>>p.pos[1]>>p.pos[2]>>p.vel[0]>>p.vel[1]>>p.vel[2];
      object.particles.push_back(p);
      continue;
    }
    if (ch == 's') {
      in>>s.i>>s.j>>s.r;
      object.springs.push_back(s);
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
  std::cout<<"read "<<object.particles.size()<<" vertices and "<<object.triangles.size()<<" triangles."<<std::endl;
  return true;
}

void writeObj(char *fname, const std::vector<Particle> &meshPts, const std::vector<Tri> &triangles) {
  std::cout<<"writing "<<fname<<std::endl;
  std::ofstream out;
  std::vector<Particle>::const_iterator p;
  std::vector<Tri>::const_iterator t;
  
  out.open(fname);

  for (p=meshPts.begin(); p!=meshPts.end(); p++) {
    out<<"v "<<p->pos[0]<<" "<<p->pos[1]<<" "<<p->pos[2]<<std::endl;
  }
  for (t=triangles.begin(); t!=triangles.end(); t++) {
    out<<"f "<<(*t)[0]+1<<" "<<(*t)[1]+1<<" "<<(*t)[2]+1<<std::endl;
  }
  
  out.close();
}


