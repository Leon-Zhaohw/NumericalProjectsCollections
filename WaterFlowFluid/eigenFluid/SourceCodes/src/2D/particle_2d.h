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

#ifndef PARTICLE_2D_H
#define PARTICLE_2D_H

#include "Alg/VEC2.h"

struct Particle2D {
  Particle2D(const float px, const float py){
    position[0] = px;
    position[1] = py;
  }
  ~Particle2D(){}
  VEC2 position;
  VEC2 velocity;
};

#endif  // PARTICLE_2D_H
