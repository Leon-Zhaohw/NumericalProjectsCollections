/* 
 * File: CircularObstacle.hpp
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

#ifndef CIRCULAROBSTACLE_HPP
#define CIRCULAROBSTACLE_HPP

#include "Obstacle.hpp"

class CircularObstacle: public Obstacle {

protected:
  FLOAT radius;

  void setBoundaries(uint w);
  void setEquivalentSources(uint w);

public:
  CircularObstacle(FLOAT r);
  CircularObstacle(VEC2 p, FLOAT r);
  ~CircularObstacle();

  FLOAT getGrid(int i, int j) const;
};

#endif
