#ifndef rigidbody_h
#define rigidbody_h
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>

struct Particle {
  Eigen::Vector3d pos;
};

struct RigidBody {
  double mass;
  Eigen::Matrix3d I0;
  Eigen::Vector3d pos, P, frc;
  Eigen::Quaterniond q;
  Eigen::Vector3d L, trq;
};

class Tri {
public:
  int indices[3];
  inline int &operator[](const unsigned int &i) { return indices[i];};
  inline int operator[](const unsigned int &i) const { return indices[i];};
};

struct SimulationParameters {
  double dt, total_time, cor;
  std::string output_fname;
};

struct Object {
  RigidBody rb;
  std::vector<Particle> particles;
  std::vector<Tri> triangles;
};


// I/O
#include "json/json.h"
#include <fstream>

bool readObject(const char *fname, Object &object);
bool readInputFile(const char *fname, SimulationParameters &params, std::vector<Object> &objects);
void writeObj(char *fname, const Object &obj);

#endif
