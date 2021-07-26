/* 
 * File: setting.hpp
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

#ifndef SETTINGS_HPP
#define SETTINGS_HPP

//#define INTERACTIVE_
//#define USE_CUDA
//#define PROJECTED_GRID // Note: only implemented with the cuda version
//#define PLOT_RESULT //Note: only for one frequency at the time, without projective grid

#include "definitions.hpp"
#include <vector>
#include <string>

#define NTHREADS_ 7

namespace settings {

  extern FLOAT dt_; //second

  // Grid spec
  extern uint n_rows_, n_cols_;
  extern FLOAT cell_size_; //meter
  extern FLOAT scale_; 
  // if the grid is projected, cell_size is irrelevant, so we use scale_ instead to go
  // from the world coordinates to the viewer coordinates

  extern uint n_rows_obs, n_cols_obs;
  extern FLOAT cell_size_obs;
 
  extern FLOAT gravity_;
  extern FLOAT h_;
  
  extern FLOAT damping_;
  FLOAT damping(FLOAT r, FLOAT k);

  // dispersion relation
  FLOAT angular_vel(FLOAT k);
  FLOAT velocity(FLOAT k);

  extern FLOAT offset_;
  extern FLOAT step_sampling_;

  extern FLOAT init_wl_;
  extern FLOAT height_ampli_;


  VEC2 grid2viewer(int i, int j);
  VEC2 gridObs2viewer(int i, int j);

  Vector2i world2grid(VEC2 p);
  Vector2i world2gridObs(VEC2 p);

  VEC2 grid2world(Vector2i p);
  VEC2 grid2world(int i, int j);
  VEC2 gridObs2world(Vector2i p);
  VEC2 gridObs2world(int i, int j);
  
  VEC2 world2viewer(FLOAT x, FLOAT y);
  VEC2 world2viewer(VEC2 p);

  VEC2 viewer2world(VEC2 p);
  VEC2 viewer2world(FLOAT x, FLOAT y);


  FLOAT interpolation(FLOAT t, int l, FLOAT dt);

  COMPLEX addWaves(FLOAT x0);
  void createTabs();

  extern std::vector<COMPLEX> hankel_tab;
  
  extern int size_tmp;

  extern bool neumann;

  extern FLOAT damping_source_;


  //profil buffer
  extern uint nb_profil;
  extern FLOAT step_profil;

  extern bool use_gersner;
  extern bool random_;

  extern bool approx_inter;

  extern std::string str_obstacle_grid;

  extern std::vector<int> ampli_steps;
};

#endif
