/* 
 * File: EquivalentSource.cpp
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

#include "EquivalentSource.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include <iostream>
#include "settings.hpp"
#include "error.hpp"

using namespace definitions;
using namespace settings;

inline COMPLEX Hankel(FLOAT x) {
  return boost::math::cyl_bessel_j<FLOAT, FLOAT>(0, x) -
    std::complex<FLOAT>(0, 1)*boost::math::cyl_neumann<FLOAT, FLOAT>(0, x);
  //return sqrtf((FLOAT)2.0/((FLOAT)M_PI*x))*exp(i_*(x - (FLOAT)(M_PI/4.0)));
}

inline COMPLEX derHankel(FLOAT x) {
  // return 0.5f*(boost::math::cyl_bessel_j<FLOAT, FLOAT>(-1, x) -
  // 	       std::complex<FLOAT>(0, 1)*boost::math::cyl_neumann<FLOAT, FLOAT>(-1, x) -
  // 	       boost::math::cyl_bessel_j<FLOAT, FLOAT>(1, x) +
  // 	       std::complex<FLOAT>(0, 1)*boost::math::cyl_neumann<FLOAT, FLOAT>(1, x));
  return sqrtf((FLOAT)2.0/((FLOAT)M_PI*x))*(-1.0f/(2.0f*x) + i_)*exp(i_*(x - (FLOAT)(M_PI/4.0)));
}

EquivalentSource::EquivalentSource(): Wave() {
  wave_length = 1;
  reset();
}

EquivalentSource::EquivalentSource(EquivalentSource* es): Wave() {
  wave_length = es->wave_length;
  ampli_step = es->ampli_step;
  reset();
  pos = es->pos;
}

EquivalentSource::EquivalentSource(FLOAT wl, int as): Wave() {
  wave_length = wl;
  ampli_step = as;
  reset();
}


void EquivalentSource::reset() {
  wave_number = 2*M_PI/wave_length;
  angular_freq = angular_vel(wave_number);
  pos = VEC2(0, 0);
  ampli = 1;

  amplis = std::vector<COMPLEX>(size_tmp);
  for (uint i = 0; i < size_tmp; ++i) {
#ifdef INTERACTIVE_
    amplis[i] = 0;
#else
    amplis[i] = 1;
#endif
  }
  is_active = false;
  FLOAT size_grid = std::max(n_rows_, n_cols_)*cell_size_*sqrt(2);
  FLOAT speed = velocity(wave_number);
#ifdef PROJECTED_GRID
  FLOAT rmax = log(0.02)/(-damping_*wave_number*wave_number);
  FLOAT time_lim = rmax/speed;
#else 
  FLOAT time_lim = size_grid/speed;
#endif
  lim =  2*time_lim / (ampli_step*dt_) + 1;
  inactive_time = 0;
}

void EquivalentSource::update(FLOAT dt) {}

FLOAT EquivalentSource::height(FLOAT x, FLOAT y, FLOAT time) const {
  COMPLEX ampli_cur = ampli;
  FLOAT rx = x - pos(0);
  FLOAT ry = y - pos(1);
  FLOAT r  = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
  FLOAT damp = damping(r, wave_number);
  FLOAT out = 0;
  if (damp > 0.02) {
#ifdef INTERACTIVE_
    FLOAT ret = r/velocity(wave_number);
    int t = floor((time-ret)/((FLOAT)ampli_step*dt_));
    if (time-ret < 0) {
      --t;
    }
    
    // if (t >= size_tmp-1) {
    //   ampli_cur = amplis[size_tmp-1];
    // } else
    if (t < 0) {
      ampli_cur = amplis[0];
    } else {
      FLOAT w = interpolation(time-ret, t, ampli_step*dt_);
      ampli_cur = w * amplis[t%size_tmp] + (1-w) * amplis[(t+1)%size_tmp];
    }
#endif
    if (r != 0) {
      out = damp*real(addWaves(wave_number*r)*ampli_cur);
    }
  }
  return out;
}

FLOAT EquivalentSource::height(VEC2 p, FLOAT time) const {
  return height(p(0), p(1), time);
}

  
COMPLEX EquivalentSource::heightc(VEC2 p, FLOAT time) const {
  return heightc(p(0), p(1), time);
}


COMPLEX EquivalentSource::heightc(FLOAT x, FLOAT y, FLOAT time) const {
  COMPLEX ampli_cur = ampli;
  FLOAT rx = x - pos(0);
  FLOAT ry = y - pos(1);
  FLOAT r  = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
  FLOAT damp = damping(r, wave_number);
  COMPLEX out(0, 0);
  if (damp > 0.02) {
#ifdef INTERACTIVE_
    FLOAT ret = r/velocity(wave_number);
    int t = floor((time-ret)/((FLOAT)ampli_step*dt_));
    if (time-ret < 0) {
      --t;
    }
  
    // if (t >= size_tmp-1) {
    //   ampli_cur = amplis[size_tmp-1];
    // } else
    if (t < 0) {
      ampli_cur = amplis[0];
    } else {
      FLOAT w = interpolation(time-ret, t, ampli_step*dt_);
      ampli_cur = w * amplis[t%size_tmp] + (1-w) * amplis[(t+1)%size_tmp];
    }
#endif
    if (r != 0) {
      out = damp*addWaves(wave_number*r)*ampli_cur;
    }
  }
  return out;
  
}


COMPLEX EquivalentSource::heightc(VEC2 p) const {
  return heightc(p(0), p(1));
}

COMPLEX EquivalentSource::heightc(FLOAT x, FLOAT y) const {
  FLOAT rx = x - pos(0);
  FLOAT ry = y - pos(1);
  FLOAT r  = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
  FLOAT damp = damping(r, wave_number);
  COMPLEX out(0, 0);

  if (damp > 0.02) {
    if (r != 0) {
      out = damp*addWaves(wave_number*r);
    }
  }
  return out;
}



VEC2C EquivalentSource::gradHeightc(FLOAT x, FLOAT y, FLOAT time) const {
  COMPLEX ampli_cur = ampli;
  FLOAT rx = x - pos(0);
  FLOAT ry = y - pos(1);
  FLOAT r  = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
  COMPLEX out_x(0, 0);
  COMPLEX out_y(0, 0);
  FLOAT damp = damping(r, wave_number);
  if (damp > 0.02) {
#ifdef INTERACTIVE_
    FLOAT ret = r/velocity(wave_number);
    int t = floor((time-ret)/((FLOAT)ampli_step*dt_));
    if (time-ret < 0) {
      --t;
    }
  
    // if (t >= size_tmp-1) {
    //   ampli_cur = amplis[size_tmp-1];
    // } else
    if (t < 0) {
      ampli_cur = amplis[0];
    } else {
      FLOAT w = interpolation(time-ret, t, ampli_step*dt_);
      ampli_cur = w * amplis[t%size_tmp] + (1-w) * amplis[(t+1)%size_tmp];
    }
#endif
    if (r != 0) {
      FLOAT der_damp = 0;//exp(-damping*pow(wave_number, 2)*r);
      FLOAT cos_phi = rx/r;
      FLOAT sin_phi = ry/r;
      out_x = cos_phi*damp*(-i_/(FLOAT)4.0*wave_number*derHankel(wave_number*r)*ampli_cur) -
	sin_phi/r*der_damp*(-i_/(FLOAT)4.0*Hankel(wave_number*r)*ampli_cur);
      out_y = sin_phi*damp*(-i_/(FLOAT)4.0*wave_number*derHankel(wave_number*r)*ampli_cur) +
	cos_phi/r*der_damp*(-i_/(FLOAT)4.0*Hankel(wave_number*r)*ampli_cur);
    }
  }
  return VEC2C(out_x, out_y);
}

VEC2C EquivalentSource::gradHeightc(FLOAT x, FLOAT y) const {
  FLOAT rx = x - pos(0);
  FLOAT ry = y - pos(1);
  FLOAT r  = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
  COMPLEX out_x(0, 0);
  COMPLEX out_y(0, 0);
  FLOAT damp = damping(r, wave_number);
  if (damp > 0.02) {
    if (r != 0) {
      FLOAT der_damp = 0;//damping(r, wave_number);//exp(-damping*pow(wave_number, 2)*r);
      FLOAT cos_phi = rx/r;
      FLOAT sin_phi = ry/r;
      out_x = cos_phi*damp* (-i_/(FLOAT)4.0*wave_number*derHankel(wave_number*r)) -
	sin_phi/r*der_damp*(-i_/(FLOAT)4.0*Hankel(wave_number*r));
      out_y = sin_phi*damp* (-i_/(FLOAT)4.0*wave_number*derHankel(wave_number*r)) +
	cos_phi/r*der_damp*(-i_/(FLOAT)4.0*Hankel(wave_number*r));
    }
  }
  return VEC2C(out_x, out_y);
}

VEC2C EquivalentSource::gradHeightc(VEC2 p, FLOAT time) const {
  return gradHeightc(p(0), p(1), time);
}


VEC2C EquivalentSource::gradHeightc(VEC2 p) const {
  return gradHeightc(p(0), p(1));
}


void EquivalentSource::setPos(VEC2 p) {
  pos = p;
}

void EquivalentSource::setPos(FLOAT x, FLOAT y) {
  pos = VEC2(x, y);
}

void EquivalentSource::setAmplitude(COMPLEX a) {
  ampli = a;
  amplis[0] = a;
}

void EquivalentSource::setAmplitude(COMPLEX a, int t) {
  //assert(t < size_tmp);
  if (t <= 0) {
    //   amplis[0] = a;
  } else {
    if (t%size_tmp == 0) {
      for (int i = 0; i < size_tmp/2; ++i) {
	amplis[i] = 0;
      }
    }
    if (t%size_tmp == size_tmp/2) {
      for (int i = size_tmp/2; i < size_tmp; ++i) {
	amplis[i] = 0;
      }
    }
    amplis[t%size_tmp] = a;
  }
  ampli = a;
}

void EquivalentSource::setAmplitude(COMPLEX a, uint t1, uint t2) {
  ampli = a;
  assert(t1 < size_tmp && t2 < size_tmp);
  for (uint i = t1; i <= t2; ++i) {
    amplis[i] = a;
  }
}



VEC2 EquivalentSource::getPos() const {
  return pos;
}

FLOAT EquivalentSource::getDistance(FLOAT x, FLOAT y) const {
  VEC2 v = pos - VEC2(x, y);
  return v.norm();
}

FLOAT EquivalentSource::getDistance(VEC2 p) const {
  VEC2 v = pos - p;
  return v.norm();
}


void EquivalentSource::writeAmplis() const {
  for (uint i = 0; i < 100; ++i) {
    std::cout<<amplis[i]<<" ";
  }
  std::cout<<"\n";
}


COMPLEX EquivalentSource::getAmpli(int i) const {
  COMPLEX a(0, 0);
  if (i > 0) {
    i = i%size_tmp; 
    a = amplis[i];
  } else {
    a = amplis[0];
  }
  return a;
}

void EquivalentSource::setActive(int t) const {
  COMPLEX a(0, 0);
  if (t > 0) {
    t = t%size_tmp;
    a = amplis[t];
  }
  if (norm(a) > 1e-10) {
    is_active = true;
    inactive_time = 0;
  } else {
    ++inactive_time;
    if (inactive_time > lim) {
      is_active = false;
 
    }
  }
}

COMPLEX EquivalentSource::getAmpli() const {
  return ampli;
}

uint EquivalentSource::getIndex() const {
  return 0;
}

void EquivalentSource::move(VEC2 trans) {
  pos += trans;
}
