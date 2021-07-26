/* 
 * File: MovingEquivalentSource.cpp
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


#include "MovingEquivalentSource.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include "settings.hpp"
#include "error.hpp"
#include "definitions.hpp"

using namespace settings;
using namespace definitions;

inline COMPLEX Hankel(FLOAT x) {
  return boost::math::cyl_bessel_j<FLOAT, FLOAT>(0, x) -
    std::complex<FLOAT>(0, 1)*boost::math::cyl_neumann<FLOAT, FLOAT>(0, x);
  ///return sqrtf((FLOAT)2.0/((FLOAT)M_PI*x))*exp(i_*(x - (FLOAT)(M_PI/4.0)));
}

inline COMPLEX derHankel(FLOAT x) {
  return 0.5f*(boost::math::cyl_bessel_j<FLOAT, FLOAT>(-1, x) -
	       std::complex<FLOAT>(0, 1)*boost::math::cyl_neumann<FLOAT, FLOAT>(-1, x) -
	       boost::math::cyl_bessel_j<FLOAT, FLOAT>(1, x) +
	       std::complex<FLOAT>(0, 1)*boost::math::cyl_neumann<FLOAT, FLOAT>(1, x));
  //  return sqrtf((FLOAT)2.0/((FLOAT)M_PI*x))*(-1.0f/(2.0f*x) + i_)*exp(i_*(x - (FLOAT)(M_PI/4.0)));
}


MovingEquivalentSource::MovingEquivalentSource(): EquivalentSource() {}

MovingEquivalentSource::MovingEquivalentSource(FLOAT wl, int as): EquivalentSource(wl, as) {
  positions = std::vector<VEC2>(size_tmp);
}

VEC2 MovingEquivalentSource::getPos(int i) const {
  VEC2 a(0, 0);
  if (i > 0) {
    a = positions[i%size_tmp];
  } else {
    a = positions[0];
  }
  return a;
}


void MovingEquivalentSource::setPos(VEC2 p, int t) {
  //  TEST(t < size_tmp);
  if (t >= 0) {
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
    positions[t%size_tmp] = p;
  }
}

void MovingEquivalentSource::setPos(VEC2 p, int t1, int t2) {
  TEST(t1 < size_tmp);
  TEST(t2 < size_tmp);
  for (int i = t1; i <= t2; ++i) {
    positions[i] = p;
  }
}

void MovingEquivalentSource::setPos(FLOAT x, FLOAT y, int t) {
  VEC2 p(x, y);
  setPos(p, t);
}
void MovingEquivalentSource::setPos(FLOAT x, FLOAT y, int t1, int t2) {
  VEC2 p(x, y);
  setPos(p, t1, t2);
}

FLOAT MovingEquivalentSource::height(FLOAT x, FLOAT y, FLOAT time) const{
  COMPLEX out = heightc(x, y, time);
  return real(out);
}
FLOAT MovingEquivalentSource::height(VEC2 p, FLOAT time) const {
  return height(p(0), p(1), time);
}

COMPLEX MovingEquivalentSource::heightc(FLOAT x, FLOAT y, FLOAT time) const {
  COMPLEX ampli_cur = ampli;
  COMPLEX out = 0;
  int t_cur = floor((time)/((FLOAT)ampli_step*dt_));
  for (int past_t = std::min(t_cur-lim, 0); past_t <= t_cur; ++past_t) {
    FLOAT rx = x - positions[past_t%size_tmp](0);
    FLOAT ry = y - positions[past_t%size_tmp](1);
    FLOAT r  = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
    if (r != 0) {

      FLOAT damp = damping(r, wave_number);
      if (damp > 0.02) {
	  
	FLOAT ret = r/velocity(wave_number);
	int t = floor((time-ret)/((FLOAT)ampli_step*dt_));
	if (time-ret > 0) {
	   if (t >= past_t - 1 || t <= past_t +1) {
	     int pt = past_t;
	     //  for (int pt = past_t - 1; pt <= past_t + 1; ++pt) {
	    FLOAT w = interpolation(time-ret, pt, ampli_step*dt_);
	    ampli_cur = w * amplis[pt%size_tmp];
	    out += damp*addWaves(wave_number*r)*ampli_cur;
	  }
	}
      }
    }
  }
  return out;
}
COMPLEX  MovingEquivalentSource::heightc(VEC2 p, FLOAT time) const {
  return heightc(p(0), p(1), time);
}
  
VEC2C MovingEquivalentSource::gradHeightc(FLOAT x, FLOAT y, FLOAT time) const {
  COMPLEX ampli_cur = ampli;
  COMPLEX out_x = 0, out_y = 0;
  int t_cur = floor((time)/((FLOAT)ampli_step*dt_));
  for (int past_t = std::min(t_cur-lim, 0); past_t <= t_cur; ++past_t) {
    FLOAT rx = x - positions[past_t%size_tmp](0);
    FLOAT ry = y - positions[past_t%size_tmp](1);
    FLOAT r  = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
    if (r != 0) {

      FLOAT damp = damping(r, wave_number);
      if (damp > 0.02) {
	FLOAT der_damp = 0;
	FLOAT cos_phi = rx/r;
	FLOAT sin_phi = ry/r;
	FLOAT ret = r/velocity(wave_number);
	int t = floor((time-ret)/((FLOAT)ampli_step*dt_));
	if (time-ret > 0) {
	  if (t >= past_t - 1 || t <= past_t +1) {
	    FLOAT w = interpolation(time-ret, past_t, ampli_step*dt_);
	    ampli_cur = w * amplis[past_t%size_tmp];
	    out_x += cos_phi*damp*(-i_/(FLOAT)4.0*wave_number*derHankel(wave_number*r)*ampli_cur) -
	      sin_phi/r*der_damp*(-i_/(FLOAT)4.0*Hankel(wave_number*r)*ampli_cur);
	    out_y += sin_phi*damp* (-i_/(FLOAT)4.0*wave_number*derHankel(wave_number*r)*ampli_cur) +
	      cos_phi/r*der_damp*(-i_/(FLOAT)4.0*Hankel(wave_number*r)*ampli_cur);
	  }
	}
      }
    }
  }
  return  VEC2C(out_x, out_y);
}
VEC2C MovingEquivalentSource::gradHeightc(VEC2 p, FLOAT time) const {
  return gradHeightc(p(0), p(1), time);
}
  
void MovingEquivalentSource::move(VEC2 trans, uint i) {
  pos += trans;
  positions[i] = pos;
}
  
