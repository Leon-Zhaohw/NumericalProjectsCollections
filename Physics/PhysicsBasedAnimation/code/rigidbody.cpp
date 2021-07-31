#include "rigidbody.h"

static const double FEPS = 1e-6;

Eigen::Matrix3d star(const Eigen::Vector3d &x) {
  Eigen::Matrix3d s;
  s << 0, -x[2], x[1],
    x[2], 0, -x[0],
    -x[1], x[0], 0;
  return s;
}

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
	writeObj(fname, objects[o]);
	frameTime = 1.0/30.0-0.0001;
	frame++;
      }
    }

    // takes a timestep
    for (std::vector<Object>::iterator obj = objects.begin(); obj != objects.end(); obj++) {
      // apply gravity
      obj->rb.frc = Eigen::Vector3d::Zero();
      obj->rb.trq = Eigen::Vector3d::Zero();
      obj->rb.frc[2] -= 9.8 * obj->rb.mass;

      // compute Inv, omega
      Eigen::Vector3d vel = obj->rb.P / obj->rb.mass;
      Eigen::Matrix3d R = obj->rb.q.toRotationMatrix();
      Eigen::Matrix3d Iinv = R * obj->rb.I0.inverse() * R.transpose();
      Eigen::Vector3d omega = Iinv * obj->rb.L;

      // time integration
      obj->rb.P += params.dt * obj->rb.frc;
      obj->rb.L += params.dt * obj->rb.trq;
      obj->rb.pos += params.dt * (obj->rb.P / obj->rb.mass);
      Eigen::Quaterniond update(0.0, omega[0], omega[1], omega[2]);
      update *= obj->rb.q;
      obj->rb.q.w() += 0.5*params.dt*update.w();
      obj->rb.q.x() += 0.5*params.dt*update.x();
      obj->rb.q.y() += 0.5*params.dt*update.y();
      obj->rb.q.z() += 0.5*params.dt*update.z();
      obj->rb.q.normalize();

      vel = obj->rb.P / obj->rb.mass;
      R = obj->rb.q.toRotationMatrix();
      Iinv = R * obj->rb.I0.inverse() * R.transpose();
      omega = Iinv * obj->rb.L;

      // ground plane / collision normal
      Eigen::Vector3d n(0.0, 0.0, 1.0);

      for (unsigned int i=0; i<obj->particles.size(); i++) {
	Eigen::Vector3d pt = obj->rb.pos + (R * obj->particles[i].pos);
	if (pt[2] < 0.0) {
	  std::cout<<"collision "<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<std::endl;
	  Eigen::Vector3d r = pt - obj->rb.pos; // lever arm
	  Eigen::Vector3d v = vel + omega.cross(r);  // velocity of the point
	  double vrel = v[2];
	  if (vrel > 0) continue;
	  double j = (-(1.0+params.cor)*vrel) / ((1.0/obj->rb.mass) + n.dot((Iinv * (r.cross(n))).cross(r)));
	  Eigen::Vector3d impulse = j * n;
	  obj->rb.P += impulse;
	  obj->rb.L += r.cross(impulse);

	  // update velocity
	  vel = obj->rb.P / obj->rb.mass;
	  omega = Iinv * obj->rb.L;
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
  params.cor = root.get("cor", 0.75).asDouble();
  params.output_fname = root.get("output_fname", std::string("output-%02d.%04d.obj")).asString();
  
  Json::Value objectsIn = root["objects"];
  objects.resize(objectsIn.size());
  for (unsigned int i=0; i<objectsIn.size(); i++) {
    readObject((objectsIn[i]["filename"]).asString().c_str(), objects[i]);
    objects[i].rb.mass = objectsIn[i].get("mass", 1.0).asDouble();
    objects[i].rb.I0 << objectsIn[i]["inertia"][0].asDouble(), objectsIn[i]["inertia"][1].asDouble(), objectsIn[i]["inertia"][2].asDouble(),
      objectsIn[i]["inertia"][3].asDouble(), objectsIn[i]["inertia"][4].asDouble(), objectsIn[i]["inertia"][5].asDouble(),
      objectsIn[i]["inertia"][6].asDouble(), objectsIn[i]["inertia"][7].asDouble(), objectsIn[i]["inertia"][8].asDouble();
    objects[i].rb.pos << objectsIn[i]["position"][0].asDouble(), objectsIn[i]["position"][1].asDouble(), objectsIn[i]["position"][2].asDouble();
    objects[i].rb.P = Eigen::Vector3d::Zero();
    //objects[i].rb.R = Eigen::Matrix3d::Identity();
    objects[i].rb.q.setIdentity();
    objects[i].rb.L = Eigen::Vector3d::Zero();
    Json::Value lmomentum = objectsIn[i]["lmomentum"];
    if (!lmomentum.isNull()) {
      objects[i].rb.P << lmomentum[0].asDouble(), lmomentum[1].asDouble(), lmomentum[2].asDouble();
    }
    Json::Value amomentum = objectsIn[i]["amomentum"];
    if (!amomentum.isNull()) {
      objects[i].rb.L << amomentum[0].asDouble(), amomentum[1].asDouble(), amomentum[2].asDouble();
    }
  }
  return true;
}

bool readObject(const char *fname, Object &object) {
  char ch;
  Particle p;
  Tri t;
  
  std::ifstream in(fname, std::ios::in);
  while (in>>ch) {
    if (ch == 'p') {
      in>>p.pos[0]>>p.pos[1]>>p.pos[2];
      object.particles.push_back(p);
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

void writeObj(char *fname, const Object &obj) {
  std::cout<<"writing "<<fname<<std::endl;
  std::ofstream out;
  std::vector<Particle>::const_iterator p;
  std::vector<Tri>::const_iterator t;
  Eigen::Matrix3d R = obj.rb.q.toRotationMatrix();
  
  out.open(fname);
  
  for (p=obj.particles.begin(); p!=obj.particles.end(); p++) {
    Eigen::Vector3d pt = obj.rb.pos + (R * (p->pos));
    out<<"v "<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<std::endl;
  }
  
  for (t=obj.triangles.begin(); t!=obj.triangles.end(); t++) 
    out<<"f "<<(*t)[0]+1<<" "<<(*t)[1]+1<<" "<<(*t)[2]+1<<std::endl;
  
  out.close();
}


