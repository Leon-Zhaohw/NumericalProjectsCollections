/* 
 * File: Grid.hpp
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

#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include "definitions.hpp"
#include <SDL2/SDL_image.h>

class Grid {

private:
  uint n_rows, n_cols, n_nodes;
  FLOAT cell_size;
  std::vector<FLOAT> nodes;

  inline uint index(uint i, uint j) const;
  inline uint row(uint ind) const ;
  inline uint col(uint ind) const;
  
public:
  Grid();
  Grid(uint n_rows, uint n_cols, FLOAT cs);
  ~Grid();
  
  bool isEmpty() const;
  
  uint getNbRows() const;
  uint getNbCols() const;
  FLOAT getCellSize() const;
  void setCellSize(FLOAT cs);
  uint getHeight() const;
  uint getWidth() const;
  VEC2 toWorld(uint i, uint j) const;

  FLOAT operator()(int i, int j) const;
  FLOAT &operator()(int i, int j);
  FLOAT value(FLOAT x, FLOAT y) const;
  FLOAT interpolatedValue(FLOAT x, FLOAT y) const;

  void reset(FLOAT val);
  void setValues(const SDL_Surface *texture);

  void getPixels(unsigned char *pixels) const;
  std::ostream& exportObj(std::ostream& os) const;
  
  Grid& operator=(const Grid& g);
  Grid& operator+=(const Grid& g);
  friend std::ostream& operator<<(std::ostream& os, const Grid& g);
  friend std::istream& operator>>(std::istream& is, Grid& F);
};

#endif
