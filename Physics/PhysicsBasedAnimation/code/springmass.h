#ifndef springmass_h
#define springmass_h
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>

struct Particle {
  double mass;
  Eigen::Vector3d pos, vel, frc;
};

struct Spring {
  double r; // rest length
  int i, j; 
};

class Tri {
public:
  int indices[3];
  inline int &operator[](const unsigned int &i) { return indices[i];};
  inline int operator[](const unsigned int &i) const { return indices[i];};
};

struct SimulationParameters {
  double dt, total_time, k_s, k_d;
  std::string output_fname;
};

struct Object {
  std::vector<Particle> particles;
  std::vector<Spring> springs;
  std::vector<Tri> triangles;
};


// I/O
#include "json/json.h"
#include <fstream>

bool readObject(const char *fname, Object &object);
bool readInputFile(const char *fname, SimulationParameters &params, std::vector<Object> &objects);
void writeObj(char *fname, const std::vector<Particle> &meshPts, const std::vector<Tri> &triangles);

#endif
