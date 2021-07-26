/* 
 * File: Grid.cpp
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

#include "Grid.hpp"
#include <iostream>
#include "settings.hpp"
#include "error.hpp"

using namespace settings;

inline uint Grid::index(uint i, uint j) const {
  return n_cols*i + j;
}
  
inline uint Grid::row(uint ind) const {
  return ind/n_cols;
}

inline uint Grid::col(uint ind) const {
  return ind - (row(ind)*n_cols);
}


Grid::Grid() {
  n_rows = 0;
  n_cols = 0;
  cell_size = 0;
}

Grid::Grid(uint rows, uint cols, FLOAT cs) {
  n_rows = rows;
  n_cols = cols;
  n_nodes = n_rows*n_cols;
  cell_size = cs;

  nodes = std::vector<FLOAT>(n_nodes);
  reset(0);
}

Grid::~Grid() {}

bool Grid::isEmpty() const {
  return n_nodes == 0;
}

uint Grid::getNbRows() const {
  return n_rows;
}

uint Grid::getNbCols() const {
  return n_cols;
}

FLOAT Grid::getCellSize() const {
  return cell_size;
}

void Grid::setCellSize(FLOAT cs) {
  cell_size = cs;
}

uint Grid::getHeight() const {
  return n_cols*cell_size;
}
uint Grid::getWidth() const {
  return n_rows*cell_size;
}

VEC2 Grid::toWorld(uint i, uint j) const {
  VEC2 out(i*cell_size, j*cell_size);
  return out;
}

FLOAT Grid::operator()(int i, int j) const {
  if (i < 0 || i >= n_rows || j < 0 || j >= n_cols) {
    return 0;
  }
  uint ind = index(i, j);
  return nodes[ind];
}
FLOAT & Grid::operator()(int i, int j) {
  if (i < 0 || i >= n_rows || j < 0 || j >= n_cols) {
    nodes[0] = 0;
    return nodes[0];
  }
  uint ind = index(i, j);
  return nodes[ind];
}

FLOAT Grid::interpolatedValue(FLOAT x, FLOAT y) const {
  ERROR(false, "TODO: Grid::interpolatedValue", "");
}

void Grid::reset(FLOAT val) {
#pragma omp for
  for (uint i = 0; i < n_nodes; ++i) {
    nodes[i] = val;
  }
}
  
void Grid::setValues(const SDL_Surface * texture) {
  uint w = texture->w;
  uint h = texture->h;
  FLOAT pix_per_cell_x = (FLOAT)h/(FLOAT)n_rows;
  FLOAT pix_per_cell_y = (FLOAT)w/(FLOAT)n_cols;
  FLOAT byte_per_cell_x = texture->format->BytesPerPixel * pix_per_cell_x;
  FLOAT byte_per_cell_y = texture->format->BytesPerPixel * pix_per_cell_y;
  unsigned char* pixels = (unsigned char*) texture->pixels;
  
#pragma omp for
  for (uint i = 0; i < n_rows; ++i) {
    uint x = i*byte_per_cell_x/texture->format->BytesPerPixel;//
    x *= texture->format->BytesPerPixel;
    for (uint j = 0; j < n_cols; ++j) {
      uint y = j*byte_per_cell_y/texture->format->BytesPerPixel;
      y *= texture->format->BytesPerPixel;
      
      nodes[i*n_cols+j] = (uint)pixels[(uint)(x*w+y)]/256.0;
      if (nodes[i*n_cols+j] < 0.1) {
	nodes[i*n_cols+j] = 0;
      } else {
	nodes[i*n_cols+j] = 1;
      } 
    }
  }
}

std::ostream& Grid::exportObj(std::ostream& os) const {
  INFO("export Grid");
  os << "# Grid\n";
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) {
      VEC2 posv = VEC2(2*i*cell_size/scale_ - 1, 2*(FLOAT)j*cell_size/scale_ - 1);
      FLOAT h = nodes[index(i, j)];
      if (h == 0 && cell_size == cell_size_obs) { 
       	h = -1;
      } 
      os <<"v "<<posv(0)<<" "<<posv(1)<<" "<<h<<"\n";
    }
  }
  for (int i = 0; i < n_rows-1; i++) {
    for (int j = 0; j < n_cols-1; j++) {
      int idx = j + i * n_cols + 1;
      int J   = 1;
      int I   = n_cols;
      os << "f "<<idx<<" "<<idx + I<<" "<<idx + J<<"\n";
      os << "f "<<idx + J<<" "<<idx + I<<" "<<idx + I + J<<"\n";
    }
  }
  return os;
}


Grid& Grid::operator=(const Grid& g) {
  n_rows = g.n_rows;
  n_cols = g.n_cols;
  n_nodes = n_rows*n_cols;
  cell_size = g.cell_size;

  nodes = std::vector<FLOAT>(n_nodes);
  for (uint i = 0; i < n_nodes; ++i) {
    nodes[i] = g.nodes[i];
  }
  return *this;
}

std::ostream& operator<<(std::ostream& os, const Grid& g) {
  os << g.n_rows<<" "<<g.n_cols<<" ";
  for (uint i = 0; i < g.n_nodes; ++i) {
    os << g.nodes[i]<<" ";
  }
  return os;
}
std::istream& operator >> (std::istream& is, Grid& g) {
  is >> g.n_rows >> g.n_cols;
  g.n_nodes = g.n_rows*g.n_cols;
  g.nodes = std::vector<FLOAT>(g.n_nodes);
  for (uint i = 0; i < g.n_nodes; ++i) {
    is >> g.nodes[i];
  }
  return is;
}
