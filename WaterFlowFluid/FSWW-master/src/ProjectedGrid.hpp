/* 
 * File: ProjectedGrid.hpp
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

#ifndef PROJECTED_GRID_HPP
#define PROJECTED_GRID_HPP

#include "settings.hpp"
#include <vector>

class ProjectedGrid {
private:
  uint n_rows, n_cols, n_nodes;
  std::vector<VEC3> viewer_pos;
  std::vector<VEC2> displacement;
  std::vector<FLOAT> sizes;
  FLOAT min_size;

  inline uint index(uint i, uint j) const;
  inline uint row(uint ind) const ;
  inline uint col(uint ind) const;

public:
  ProjectedGrid();
  ProjectedGrid(uint nr, uint nc);
  ~ProjectedGrid();

  VEC3 operator()(uint i, uint j) const;
  //VEC3 &operator()(uint i, uint j);
  VEC3 operator()(uint i) const;
  //VEC3 &operator()(uint i);
  VEC2 getPosWorld(int i, int j) const;
  VEC2 getPosWorld(int i) const;
  
  void setPosOnThePlane(uint i, uint j, FLOAT x, FLOAT y);
  void setHeight(uint i, uint j, FLOAT z);
  void setDisplacement(uint i, uint j, VEC2 d);
  void setDisplacement(uint i, uint j, FLOAT x, FLOAT y);
  
  void setSizes();
  FLOAT minSize();
  FLOAT size(uint i);

  std::ostream& exportObj(std::ostream& os) const;
  
  friend std::ostream& operator<<(std::ostream& os, const ProjectedGrid& g);
};

#endif
