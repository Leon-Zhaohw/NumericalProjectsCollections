/* 
 * File: LinearWave.cpp
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

#include "LinearWave.hpp"
#include <iostream>
#include "settings.hpp"

using namespace definitions;
using namespace settings;

LinearWave::LinearWave(): Wave() {
  wave_length = 1;
  reset();
}
  
LinearWave::LinearWave(FLOAT wl): Wave() {
  wave_length = wl;
  reset();
}

void LinearWave::reset() {
  wave_number = 2*M_PI/wave_length;
  angular_freq = angular_vel(wave_number);
  dir = VEC2(sqrt(2), sqrt(2));
  dir.normalize();
  ampli = 0.01;
}
  
void LinearWave::update(FLOAT dt) {

}

FLOAT LinearWave::height(FLOAT x, FLOAT y, FLOAT time) const {
  FLOAT kx = wave_number*VEC2(x, y).dot(dir);
#ifdef INTERACTIVE_
  if (time > x/velocity(wave_number)) {
    return real(ampli*exp(i_*(kx - angular_freq*time)));
  } else {
    return 0;
  }
#else
  return real(ampli*exp(i_*(kx - angular_freq*time))*phase_exp);
#endif
  
}

FLOAT LinearWave::height(VEC2 p, FLOAT time) const {
  FLOAT kx = wave_number*p.dot(dir);
#ifdef INTERACTIVE_
  if (time > p(0)/velocity(wave_number)) {
    return real(ampli*exp(i_*(kx - angular_freq*time)*phase_exp));
  } else {
    return 0;
  }
#else
  return real(ampli*exp(i_*(kx - angular_freq*time)));
#endif
}

COMPLEX LinearWave::heightc(FLOAT x, FLOAT y, FLOAT time) const {
  FLOAT kx = wave_number*VEC2(x, y).dot(dir);
#ifdef INTERACTIVE_
  if (time > x/velocity(wave_number)) {
    return ampli*exp(i_*(kx))*phase_exp;
  } else {
    return 0;
  }
#else
  return ampli*exp(i_*(kx));
#endif
}

COMPLEX LinearWave::heightc(VEC2 p, FLOAT time) const {
  FLOAT kx = wave_number*p.dot(dir);
#ifdef INTERACTIVE_
  if (time > p(0)/velocity(wave_number)) {
    return ampli*exp(i_*(kx))*phase_exp;
  } else {
    return 0;
  }
#else
  return ampli*exp(i_*(kx));
#endif
}

VEC2C LinearWave::gradHeightc(FLOAT x, FLOAT y, FLOAT time) const {
  FLOAT kdotx = wave_number*VEC2(x, y).dot(dir);
  FLOAT kx = wave_number*dir(0);
  FLOAT ky = wave_number*dir(1);
#ifdef INTERACTIVE_
  if (time > x/velocity(wave_number)) {
    return VEC2C(ampli*i_*kx*exp(i_*(kdotx - angular_freq*time)),
		 ampli*i_*ky*exp(i_*(kdotx - angular_freq*time)))*phase_exp;
  } else {
    return VEC2C(COMPLEX(0, 0), COMPLEX(0, 0));
  }
#else
  return VEC2C(ampli*i_*kx*exp(i_*(kdotx - angular_freq*time)),
	       ampli*i_*ky*exp(i_*(kdotx - angular_freq*time)));
#endif
}

VEC2C LinearWave::gradHeightc(VEC2 p, FLOAT time) const {
  FLOAT kdotx = wave_number*p.dot(dir);
  FLOAT kx = wave_number*dir(0);
  FLOAT ky = wave_number*dir(1);
#ifdef INTERACTIVE_
  if (time > p(0)/velocity(wave_number)) {
    return VEC2C(ampli*i_*kx*exp(i_*(kdotx - angular_freq*time)),
		 ampli*i_*ky*exp(i_*(kdotx - angular_freq*time)))*phase_exp;
  } else {
    return VEC2C(COMPLEX(0, 0), COMPLEX(0, 0));
  }
#else
  return VEC2C(ampli*i_*kx*exp(i_*(kdotx - angular_freq*time)),
	       ampli*i_*ky*exp(i_*(kdotx - angular_freq*time)));
#endif
}


void LinearWave::setDir(VEC2 d) {
  dir = d;
  dir.normalize();
}

void LinearWave::setDir(FLOAT x, FLOAT y) {
  dir = VEC2(x, y);
  dir.normalize();
}

void LinearWave::setAmplitude(COMPLEX a) {
  ampli = 0.1f*a;
}

FLOAT LinearWave::getDistance(FLOAT x, FLOAT y) const {
  return 0;
}

COMPLEX LinearWave::getAmpli(int i) const {
  return ampli;
}

COMPLEX LinearWave::getAmpli() const {
  return ampli;
}

uint LinearWave::getIndex() const {
  return 1;
}
VEC2 LinearWave::getPos() const {
  return dir;
}
