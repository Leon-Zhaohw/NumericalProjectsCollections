/* 
 * File: Wave.hpp
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

#ifndef WAVE_HPP
#define WAVE_HPP

#include "definitions.hpp"

class Wave {
protected:
  mutable bool is_active;
  //  mutable bool changed_active;
  COMPLEX phase_exp;
  FLOAT phase;

public:
  Wave();
  
  virtual void update(FLOAT dt) = 0;
  virtual FLOAT height(FLOAT x, FLOAT y, FLOAT time) const = 0;
  virtual FLOAT height(VEC2 p, FLOAT time) const = 0;
  virtual COMPLEX heightc(FLOAT x, FLOAT y, FLOAT time) const = 0;
  virtual COMPLEX heightc(VEC2 p, FLOAT time) const = 0;
  virtual VEC2C gradHeightc(FLOAT x, FLOAT y, FLOAT time) const = 0;
  virtual VEC2C gradHeightc(VEC2 p, FLOAT time) const = 0;

  virtual FLOAT getDistance(FLOAT x, FLOAT y) const = 0;
  virtual COMPLEX getAmpli(int i) const = 0;
  virtual COMPLEX getAmpli() const = 0;
  virtual uint getIndex() const = 0;
  virtual VEC2 getPos() const = 0;

  bool isActive() const;
  virtual void setActive(int t) const;

  void setPhase(FLOAT p);
  COMPLEX getPhaseExp() const;
  FLOAT getPhase() const;
};


#endif
