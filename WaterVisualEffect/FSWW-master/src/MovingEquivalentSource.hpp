/* 
 * File: MovingEquivalentSource.hpp
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

#ifndef MOVINGEQUIVALENTSOURCE_HPP
#define MOVINGEQUIVALENTSOURCE_HPP

#include "EquivalentSource.hpp"

class MovingEquivalentSource: public EquivalentSource {

protected:
  std::vector<VEC2> positions;

public:
  MovingEquivalentSource();
  MovingEquivalentSource(FLOAT wl, int as);
  
  VEC2 getPos(int i) const;

  void setPos(VEC2 p, int t);
  void setPos(VEC2 p, int t1, int t2);
  void setPos(FLOAT x, FLOAT y, int t);
  void setPos(FLOAT x, FLOAT y, int t1, int t2);

  FLOAT height(FLOAT x, FLOAT y, FLOAT time) const;
  FLOAT height(VEC2 p, FLOAT time) const;
  COMPLEX heightc(FLOAT x, FLOAT y, FLOAT time) const;
  COMPLEX heightc(VEC2 p, FLOAT time) const;
  
  VEC2C gradHeightc(FLOAT x, FLOAT y, FLOAT time) const;
  VEC2C gradHeightc(VEC2 p, FLOAT time) const;

  void move(VEC2 trans, uint i);
};


#endif


