/* 
 * File: SquareObstacle.cpp
 *
 * Copyright (C) 2019  Camille Schreck
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "SquareObstacle.hpp"
#include "settings.hpp"
#include <iostream>

using namespace settings;

void SquareObstacle::setBoundaries(uint w) {
  std::vector<InputPoint*> boundaries;
  std::vector<VEC2C> normals;
  FLOAT wl = wave_lenghts[w];
  FLOAT step_th = step_sampling_*wl/4.0;
  uint pts_per_side = size/step_th;
  FLOAT step = size/(FLOAT)pts_per_side;
  
  std::vector<VEC2> corners(4);
  corners[0] = VEC2(pos(0) - size/2.0, pos(1) - size/2.0);
  corners[1] = VEC2(pos(0) - size/2.0, pos(1) + size/2.0); 
  corners[2] = VEC2(pos(0) + size/2.0, pos(1) + size/2.0);
  corners[3] = VEC2(pos(0) + size/2.0, pos(1) - size/2.0);

  for (uint i = 0; i < 4; ++i) {
    InputPoint* ipc = new InputPoint(128, 0.05);
    ipc->setPos(corners[i]);
    VEC2 d = corners[i] - pos;
    d.normalize();
    VEC2C nc(COMPLEX(d(0), 0), COMPLEX(d(1), 0));
    normals.push_back(nc);
    boundaries.push_back(ipc);
    VEC2 dir = (corners[(i+1)%4] - corners[i])/size;
    for (uint j = 1; j < pts_per_side; ++j) {
      InputPoint* ip = new InputPoint(128, dt_);
      ip->setPos(corners[i] + j*step*dir);
      boundaries.push_back(ip);
      VEC2C n(COMPLEX(dir(0), 0), COMPLEX(dir(1), 0));
      normals.push_back(n);
    }
    
  }
  boundaries_l[w] = boundaries;
  normals_l[w] = normals;
}

void SquareObstacle::setEquivalentSources(uint w) {
  std::vector<EquivalentSource*> sources;
  FLOAT wl = wave_lenghts[w];
  
  FLOAT s = size - offset_;
  FLOAT step_th = step_sampling_*wl;
  uint pts_per_side = s/step_th;
  FLOAT step = s/(FLOAT)pts_per_side;
  
  std::vector<VEC2> corners(4);
  corners[0] = VEC2(pos(0) - s/2.0, pos(1) - s/2.0);
  corners[1] = VEC2(pos(0) - s/2.0, pos(1) + s/2.0); 
  corners[2] = VEC2(pos(0) + s/2.0, pos(1) + s/2.0);
  corners[3] = VEC2(pos(0) + s/2.0, pos(1) - s/2.0);

  EquivalentSource* es;
  
  for (uint i = 0; i < 4; ++i) {
    es = new EquivalentSource(wl, ampli_steps[w]);
    es->setPos(corners[i]);;
    sources.push_back(es);

    VEC2 dir = (corners[(i+1)%4] - corners[i])/s;
    for (uint j = 1; j <= pts_per_side; ++j) {
      es = new EquivalentSource(wl, ampli_steps[w]);
      es->setPos(corners[i] + j*step*dir);
      sources.push_back(es);
    }
  }
  sources_l[w] = sources;
}

SquareObstacle::SquareObstacle(FLOAT s): size(s), Obstacle() {}
  
SquareObstacle::SquareObstacle(VEC2 p, FLOAT s) : size(s), Obstacle(p) {}

SquareObstacle::~SquareObstacle() {}

FLOAT SquareObstacle::getGrid(int i, int j) const {
  FLOAT x = i*cell_size_obs - pos(0);
  FLOAT y = j*cell_size_obs - pos(1);
  FLOAT out = 0;
  if (fabs(x) < size/2.0 && fabs(y) < size/2.0) {
    out = 1;
  }
  return out;
}
