/* 
 * File: WavyObstacle.cpp
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

#include "WavyObstacle.hpp"
#include <iostream>
#include "settings.hpp"
#include "error.hpp"

using namespace settings;

void WavyObstacle::setBoundaries(uint w) {
  std::vector<InputPoint*> boundaries;
  std::vector<VEC2C> normals;
  FLOAT wl = wave_lenghts[w];

  FLOAT sample_rate_boundary = step_sampling_*wl/4.0;

  FLOAT size_max = 30;
  FLOAT x_shift = size_max - 2*ampli;
  FLOAT y_cur = 0;

  while (y_cur < size_max + sample_rate_boundary) {
    InputPoint* ip = new InputPoint(128, dt_);
    VEC2 n;
    if (pattern == 0) {
      VEC2 p_cur(x_shift + ampli*cos(freq*2*M_PI*y_cur), y_cur);
      ip->setPos(p_cur);

      n = VEC2(-ampli*sin(y_cur), 1);
      n.normalize();
    } else {
      int ns = floor(y_cur*freq);
      FLOAT weight = y_cur*freq - (FLOAT)ns;
      FLOAT p_prev, p_next;
      if (ns % 2 == 0) {
	p_prev = x_shift-ampli;
	p_next = x_shift+ampli;
	n = VEC2(-1.0/freq, -2*ampli);
      } else {
      	p_prev = x_shift+ampli;
      	p_next = x_shift-ampli;
      	n = VEC2(-1.0/freq, 2*ampli);
      }
      n.normalize();
      VEC2 p_cur((1-weight)*p_prev + weight*p_next, y_cur);
      ip->setPos(p_cur);

    }
    
    y_cur += sample_rate_boundary;
    normals.push_back(n);
    boundaries.push_back(ip);
  }
     
  boundaries_l[w] = boundaries;
  normals_l[w] = normals;
}

void WavyObstacle::setEquivalentSources(uint w) {
  std::vector<EquivalentSource*> sources;

  FLOAT wl = wave_lenghts[w];
  FLOAT sample_rate_boundary = step_sampling_*wl;

  EquivalentSource* es;

  FLOAT size_max = 30;
  FLOAT x_shift = size_max - 2*ampli + offset_;
  FLOAT y_cur = 0;
  while (y_cur < size_max + sample_rate_boundary) {
    VEC2 p_cur(0, 0);
    if (pattern == 0) {
      p_cur = VEC2(x_shift + ampli*cos(freq*2*M_PI*y_cur), y_cur);
      es = new EquivalentSource(wl, ampli_steps[w]);
      es->setPos(p_cur);
    } else {
      int ns = floor(y_cur*freq);
      FLOAT weight = y_cur*freq - (FLOAT)ns;
      FLOAT p_prev, p_next;
      if (ns % 2 == 0) {
       	p_prev = x_shift-ampli;
       	p_next = x_shift+ampli;
      } else {
       	p_prev = x_shift+ampli;
       	p_next = x_shift-ampli;
      }
      p_cur = VEC2((1.0-weight)*p_prev + weight*p_next, y_cur);
      //     INFO("weight "<< p_cur);
      es = new EquivalentSource(wl, ampli_steps[w]);
      es->setPos(p_cur);
    }
    sources.push_back(es);
    y_cur += sample_rate_boundary;
  }

#ifndef TIME_DOMAIN
  sources_l[w] = sources;
#endif
 
}

WavyObstacle::WavyObstacle(uint p, FLOAT f, FLOAT a): Obstacle(VEC2(0, 0)) {
  pattern = p;
  freq = f;
  ampli = a;
}

WavyObstacle::~WavyObstacle() {}

FLOAT WavyObstacle::getGrid(int i, int j) const {
  FLOAT size_max = 30;
  FLOAT out = 0;
  FLOAT x_shift = size_max - 2*ampli;
  VEC2 v(i*cell_size_obs, j*cell_size_obs);
  if (pattern == 0) {
    if (v(0) >= x_shift + ampli*cos(freq*2*M_PI*v(1))) {
      out = 1;
    }
  } else {
    int ns = floor(v(1)*freq);
    FLOAT weight = v(1)*freq - (FLOAT)ns;
    FLOAT p_prev, p_next;
    if (ns % 2 == 0) {
      p_prev = x_shift-ampli;
      p_next = x_shift+ampli;
    } else {
      p_prev = x_shift+ampli;
      p_next = x_shift-ampli;
    }
    if (v(0) >= (1-weight)*p_prev + weight*p_next) {
      out = 1;
    }
  }
  return out;
}
