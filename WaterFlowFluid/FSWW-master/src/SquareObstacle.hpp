/* 
 * File: SquareObstacle.hpp
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

#ifndef SQUAREOBSTACLE_HPP
#define SQUAREOBSTACLE_HPP

#include "Obstacle.hpp"

class SquareObstacle: public Obstacle {

protected:
  FLOAT size;

  void setBoundaries(uint w);
  void setEquivalentSources(uint w);

public:
  SquareObstacle(FLOAT s);
  SquareObstacle(VEC2 p, FLOAT s);
  ~SquareObstacle();

  FLOAT getGrid(int i, int j) const;
};

#endif
