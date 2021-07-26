/* 
 * File: TextureObstacle.cpp
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

#include "TextureObstacle.hpp"
#include <iostream>
#include "settings.hpp"
#include "error.hpp"

using namespace settings;

void TextureObstacle::setBoundaries(uint w) {
  FLOAT wl = wave_lenghts[w];
#ifdef INTERACTIVE_
  if (wl < 20*cell_size || cell_size*grid.getNbRows() < wl) {
    return;
  }
#endif  
  std::vector<InputPoint*> boundaries;
  std::vector<VEC2C> normals;
  Grid tmp = grid;
  Grid tmp_prev = grid;

  float os = step_sampling_*wl/cell_size/2.0;
  if (wl >= 1) {
    os /= 2;
  }
  int offset = floor(os)+1;
  uint nb_ip = 0;
  for (int i = 1; i < grid.getNbRows()-1; ++i) {
    for (int j = 1; j < grid.getNbCols()-1; ++j) {
      if (tmp(i, j) != 0 &&
	  (tmp(i-1, j) == 0 || tmp(i, j-1) == 0 || tmp(i+1, j) == 0 || tmp(i, j+1) == 0)) {

	bool no_neigh = true;
	for (int k = -offset; k <= offset; ++k) {
	  for (int h = -offset; h <= offset; ++h) {
	    if (tmp(i+k, j+h) == 0.5) {
	      no_neigh = false;
	      break;
	    }
	  }
	}

	if (no_neigh) {
	  InputPoint* ip = new InputPoint(128, dt_);
	  ip->setPos(((FLOAT)i-(FLOAT)n_rows/2.0)*cell_size + pos(0),
		     ((FLOAT)j-(FLOAT)n_cols/2.0)*cell_size + pos(1));
	  VEC2 d(0, 0);
	  for (int k = -1; k <= 1; ++k) {
	    for (int h = -1; h<= 1; ++h) {
	      if (tmp(i+k, j+h) != 0) {
		d += VEC2(-k, -h);
	      }
	    }
	  }
	  ++nb_ip;
	  d.normalize();
	  VEC2C n(COMPLEX(d(0), 0), COMPLEX(d(1), 0));
	  normals.push_back(n);
	  boundaries.push_back(ip);
	  tmp(i, j) = 0.5;
	}
	
      }
    }
  }

  if (boundaries.size() > 10) {
    boundaries_l[w] = boundaries;
    normals_l[w] = normals;
  }
}


void TextureObstacle::setEquivalentSources(uint w) {
  FLOAT wl = wave_lenghts[w];
#ifdef INTERACTIVE_
  if (wl >= 20*cell_size || cell_size*grid.getNbRows() < wl) {
    return;
  }
#endif
  std::vector<EquivalentSource*> sources;
  
  FLOAT k = 2*M_PI/wl;
  FLOAT vel = velocity(k);
  uint k_max = wl*offset_/cell_size;
  if (k_max < 1) {
    k_max = 1;
  }
  if (k_max > 25) {
    k_max = 25;
  }

  Grid tmp_prev = grid;
  Grid tmp;
  for (uint k = 0; k < k_max; ++k) {
    tmp = tmp_prev;
    for (uint i = 1; i < grid.getNbRows()-1; ++i) {
      for (uint j = 1; j < grid.getNbCols()-1; ++j) {
	if (tmp_prev(i, j) != 0 &&
	    (tmp_prev(i-1, j) == 0 || tmp_prev(i, j-1) == 0 ||
	     tmp_prev(i+1, j) == 0 || tmp_prev(i, j+1) == 0)) {
	  tmp(i, j) = 0;
	}
      }
    }
    tmp_prev = tmp;
  }
  EquivalentSource* es;
  uint nb_es = 0;
  for (int i = 1; i < grid.getNbRows()-1; ++i) {
    for (int j = 1; j < grid.getNbCols()-1; ++j) {
      if (tmp(i, j) != 0 &&
	  (tmp(i-1, j) == 0 || tmp(i, j-1) == 0 || tmp(i+1, j) == 0 || tmp(i, j+1) == 0)) {
	bool no_neigh = true;

	float os = step_sampling_*wl/cell_size/2.0;
	if (wl >= 1) {
	  os /= 2;
	}
	int offset = floor(os)+1;
	for (int k = -offset; k <= offset; ++k) {
	  for (int h = -offset; h <= offset; ++h) {
	    if (tmp(i+k, j+h) == 0.5) {
	      no_neigh = false;
	      break;
	    }
	  }
	}
	
	if (no_neigh) {
	  if (movable_) {
	    es = new MovingEquivalentSource(wl,ampli_steps[w]);
	    ((MovingEquivalentSource*)es)->setPos((FLOAT)(i-n_rows/2)*cell_size+pos(0),
						  (FLOAT)(j-n_cols/2)*cell_size+pos(1), 0);
	  } else {
	    es = new EquivalentSource(wl,ampli_steps[w]);
	    es->setPos((FLOAT)(i-n_rows/2)*cell_size+pos(0),
		       (FLOAT)(j-n_cols/2)*cell_size+pos(1));
	  }
	  ++nb_es;
	  sources.push_back(es);
	  tmp(i, j) = 0.5;
	}
      }
    }
  }
  
  sources_l[w] = sources;
}

TextureObstacle::TextureObstacle(): Obstacle() {}
  
TextureObstacle::TextureObstacle(std::string file, uint nr, uint nc, FLOAT cs): Obstacle() {
  n_rows = nr;
  n_cols = nc;
  cell_size = cs;

  height_field = IMG_Load(file.c_str());
  if(height_field == 0) {
    std::cout << "Erreur : " << SDL_GetError() << std::endl;
    std::exit(1);
  }
  grid = Grid(n_rows, n_cols, cs);
  grid.setCellSize(cs);
      
  grid.setValues(height_field);
   
  pos = VEC2(0, 0);
  file_texture = file;
}

TextureObstacle::~TextureObstacle() {
  delete height_field;
}

FLOAT TextureObstacle::getGrid(int i, int j) const {
  VEC2 pw = gridObs2world(i, j);
  pw -= pos;
  Vector2i pi = Vector2i(pw(0)/cell_size + n_rows/2, pw(1)/cell_size + n_cols/2);
  return grid(pi(0), pi(1));
}

void TextureObstacle::setCellSize(FLOAT cs) {
  cell_size = cs;
  grid.setCellSize(cs);
}

std::ofstream &TextureObstacle::exportMitsuba(std::ofstream &file) const {
  FLOAT scalex = 1.01*cell_size * n_rows/scale_;
  FLOAT scaley = 1.01*cell_size * n_cols/scale_;
  FLOAT scalez = cell_size * 4;
  VEC2 pv = world2viewer(pos);
  file<<"<shape type=\"heightfield\">\n";
  file<<"<string name=\"filename\" value=\""<<file_texture<<"\"/>\n";
  file<<"<float name=\"scale\" value=\""<<scalez<<"\"/>\n";
  file<<"<transform name=\"toWorld\">\n";
  file<<"<rotate x=\"0\" y=\"0\" z=\"1\" angle=\"-90\"/>\n";
  file<<"<scale x=\""<<scalex<<"\" y=\""<<-scaley<<"\" z=\"1\"/>\n";
  file<<"<translate x=\""<<pv(0)<<"\" y=\""<<pv(1)<<"\" z=\""<<-3.0*scalez/4.0<<"\"/>\n";
  file<<" </transform>\n";
  file<<"<bsdf type=\"diffuse\">\n";
  file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  file<<"</bsdf>\n";
  file<<"</shape>\n";
}
