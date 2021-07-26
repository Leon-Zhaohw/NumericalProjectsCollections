/* 
 * File: LinearWave.hpp
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

#ifndef LINEARWAVE_HPP
#define LINEARWAVE_HPP

#include "Wave.hpp"

class LinearWave : public Wave {
protected: 
  FLOAT wave_length;
  FLOAT wave_number;
  FLOAT angular_freq;
  VEC2 dir;
  COMPLEX ampli;

public:
  LinearWave();
  LinearWave(FLOAT wl);

  void reset();
  
  void update(FLOAT dt);
  FLOAT height(FLOAT x, FLOAT y, FLOAT time) const;
  FLOAT height(VEC2 p, FLOAT time) const;
  COMPLEX heightc(FLOAT x, FLOAT y, FLOAT time) const;
  COMPLEX heightc(VEC2 p, FLOAT time) const;
  VEC2C gradHeightc(VEC2 p, FLOAT time) const;
  VEC2C gradHeightc(FLOAT x, FLOAT y, FLOAT time) const;
  
  void setDir(VEC2 d);
  void setDir(FLOAT x, FLOAT y);
  void setAmplitude(COMPLEX a);

  FLOAT getDistance(FLOAT x, FLOAT y) const;
  COMPLEX getAmpli(int i) const;
  COMPLEX getAmpli() const;
  uint getIndex() const;
  VEC2 getPos() const;
};


#endif
