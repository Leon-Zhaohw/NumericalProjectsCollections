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

#ifndef PARTICLE_3D_H
#define PARTICLE_3D_H

#include "Eigen"

#include "Alg/VEC3.h"

class Particle3D {
public:
  Particle3D(const float px, const float py, const float pz){
    position[0] = px;
    position[1] = py;
    position[2] = pz;
  }
  Particle3D() {
    position = 0;
  }
  ~Particle3D(){}
  VEC3 position;
  VEC3 velocity;
};

#endif
