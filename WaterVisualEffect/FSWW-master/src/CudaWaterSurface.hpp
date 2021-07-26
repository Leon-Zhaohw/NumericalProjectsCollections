/* 
 * File: CudaWaterSurface.hpp
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

#ifndef CUDAWATERSURFACE_HPP
#define CUDAWATERSURFACE_HPP

#include "definitions.hpp"
#include "settings.hpp"
#include <vector>
#include <list>
//#include "Wave.hpp"

class CudaWaterSurface {

public:
  std::vector<FLOAT *>amplitudes; //for each freq : nb_sources*size_tmp
  std::vector<FLOAT *>indexes; //for each freq : nb_sources
  std::vector<FLOAT *>positions; //for each freq : 2*nb_sources
  std::vector<bool *>is_active; //for each freq : nb_sources
  std::vector<uint> nb_sources;
  std::vector<FLOAT> wave_lenghts;
  std::vector<uint> nb_sources_input;

  FLOAT* positions_grid; //2D positions of the node of the adaptive grid
  FLOAT* sizes;

  FLOAT* obs_grid;
  FLOAT* hankel_i_tab;
  FLOAT* hankel_r_tab;
  uint nb_wl;

  void init(uint nw);
  void allocMem(uint w, uint ns, uint na);

  void initM(uint nw);
  void allocMemM(uint w, uint ns, uint na);
  
 #ifndef INTERACTIVE_
  void setTimeHeight(int time);
   std::vector<FLOAT *>heights;
  FLOAT * time_heights;
  std::vector<FLOAT *> displacement;
  FLOAT* time_displacement; //2D positions of the node of the grid
 #else
  FLOAT * heights;
  FLOAT* displacement; //2D positions of the node of the grid
  FLOAT* buffer;

  std::vector<FLOAT *>amplitudes_m; //for each freq : nb_sources*size_tmp
  std::vector<FLOAT *>indexes_m; //for each freq : nb_sources
  std::vector<FLOAT *>positions_m; //for each freq : 2*nb_sources
  std::vector<bool *>is_active_m; //for each freq : nb_sources
  std::vector<uint> nb_sources_m;
  std::vector<uint> nb_sources_input_m;
#endif


  CudaWaterSurface();
  ~CudaWaterSurface();
  void clear();

  void setHeight(int time);
  void createTabs();
};

#endif
