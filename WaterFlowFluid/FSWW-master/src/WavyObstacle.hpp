/* 
 * File: WavyObstacle.hpp
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

#ifndef WAVYOBSTACLE_HPP
#define WAVYOBSTACLE_HPP

#include "Obstacle.hpp"

class WavyObstacle: public Obstacle {
  uint pattern; // 0: cosine, 1:triangle
  FLOAT freq;
  FLOAT ampli;
  
  void setBoundaries(uint w);
  void setEquivalentSources(uint w);
  
public:
  WavyObstacle(uint p, FLOAT f, FLOAT a);
  ~WavyObstacle();

  FLOAT getGrid(int i, int j) const;

};

#endif
