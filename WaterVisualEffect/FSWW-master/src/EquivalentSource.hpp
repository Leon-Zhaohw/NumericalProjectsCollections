/* 
 * File: EquivalentSource.hpp
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

#ifndef EQUIVALENTSOURCE_HPP
#define EQUIVALENTSOURCE_HPP

#include "Wave.hpp"
#include <vector>


inline COMPLEX Hankel(FLOAT x);
inline COMPLEX derHankel(FLOAT x);

class EquivalentSource: public Wave {

protected:
  FLOAT wave_length;
  FLOAT wave_number;
  FLOAT angular_freq;
  int ampli_step;
  VEC2 pos;
  COMPLEX ampli;

  std::vector<COMPLEX> amplis;

  int lim;
  mutable int inactive_time;
  
public:
  EquivalentSource();
  EquivalentSource(EquivalentSource* es);
  EquivalentSource(FLOAT wl, int as);
  void reset();
  
  void update(FLOAT dt);
  virtual FLOAT height(FLOAT x, FLOAT y, FLOAT time) const;
  virtual FLOAT height(VEC2 p, FLOAT time) const;
  virtual COMPLEX heightc(FLOAT x, FLOAT y, FLOAT time) const;
  virtual COMPLEX heightc(VEC2 p, FLOAT time) const;
  COMPLEX heightc(VEC2 p) const;
  COMPLEX heightc(FLOAT x, FLOAT y) const;
  

  virtual VEC2C gradHeightc(FLOAT x, FLOAT y, FLOAT time) const;
  virtual VEC2C gradHeightc(VEC2 p, FLOAT time) const;
  VEC2C gradHeightc(FLOAT x, FLOAT y) const;
  VEC2C gradHeightc(VEC2 p) const;
 
  void setPos(VEC2 p);
  void setPos(FLOAT x, FLOAT y);

  void setAmplitude(COMPLEX a);
  void setAmplitude(COMPLEX a, int t);
  void setAmplitude(COMPLEX a, uint t1, uint t2);

  VEC2 getPos() const;
  FLOAT getDistance(FLOAT x, FLOAT y) const;
  FLOAT getDistance(VEC2 p) const;
  uint getIndex() const;
  
  void writeAmplis() const;

  COMPLEX getAmpli(int i) const;
  COMPLEX getAmpli() const;
  void setActive(int t) const;

  void move(VEC2 trans);
};

#endif
