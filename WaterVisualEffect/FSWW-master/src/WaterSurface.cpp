/* 
 * File: WaterSurface.cpp
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

#include "WaterSurface.hpp"
#include "Viewer.hpp"
#include "EquivalentSource.hpp"
#include "LinearWave.hpp"
#include "CircularObstacle.hpp"
#include "SquareObstacle.hpp"
#include "settings.hpp"
#include "ui_parameters.hpp"
#include "TextureObstacle.hpp"
#include "WavyObstacle.hpp"
#include "error.hpp"
#include "Times.hpp"
#include <fstream>
#include <boost/math/special_functions/bessel.hpp>


using namespace definitions;
using namespace settings;
using namespace ui_parameters;

WaterSurface::WaterSurface() {
  omp_set_num_threads(NTHREADS_);
  Eigen::setNbThreads(NTHREADS_);
  
  import_ = false;
  export_ = false;
  export_mit = false;
  load_conf = false;

  stop_time = 1e4;

  data_file = "";
}

WaterSurface::~WaterSurface() {
  clear();
}

void WaterSurface::clear() {
  if (t_solve.joinable()) {
    t_solve.join();
  }

  std::list<Wave*>::iterator it;
  std::vector<std::list<Wave*> >::iterator itwf;
  for (itwf= waves.begin(); itwf != waves.end(); ++itwf) {
    for (it = (*itwf).begin(); it != (*itwf).end(); ++it) {
      delete (*it);
    }
  }
  std::list<MovingEquivalentSource*>::iterator itms;
  std::vector<std::list<MovingEquivalentSource*> >::iterator itmsf;
  for (itmsf= moving_waves.begin(); itmsf != moving_waves.end(); ++itmsf) {
    for (itms = (*itmsf).begin(); itms != (*itmsf).end(); ++itms) {
      delete (*itms);
    }
  }
  std::list<Obstacle*>::iterator ito;
  for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
    delete (*ito);
  }
  obstacles.clear();
  waves.clear();
  moving_waves.clear();
  inter_src.clear();
  ampli_steps.clear();
#ifdef USE_CUDA
  cuda_surface.clear();
#endif
}

void WaterSurface::reset() {
  clear();
  createTabs();
  
  srand (std::time(NULL));
  export_step = 1;
  
  step_wl = init_wl_;
  min_wl = 1;
  nb_wl = 1;
  max_wl = 0;

  setLists();
  
  if (load_conf) {
    importConfig(conf_file);
  }
  
#ifdef PROJECTED_GRID
  proj_grid = ProjectedGrid(n_rows_, n_cols_);
#endif
  u = Grid(n_rows_, n_cols_, cell_size_);
  u_re = Grid(n_rows_, n_cols_, cell_size_);
  u_im = Grid(n_rows_, n_cols_, cell_size_);
  b = Grid(n_rows_obs, n_cols_obs, cell_size_obs);
  dispx = Grid(n_rows_, n_cols_, cell_size_);
  dispy = Grid(n_rows_, n_cols_, cell_size_);
  
  if (!import_) {
    INFO("Wavelength: "<<nb_wl<<" wavelenths between "<<min_wl<<" and "<<max_wl);
  

#ifdef INTERACTIVE_
    // pre-defined used for interaction with the user,
    // so that the memory can be created for the gpu computations
    for (uint wl = 0; wl < nb_wl; ++wl) {
      inter_src[wl] = std::vector<EquivalentSource*>(100);
      for (uint i = 0; i < 100; ++i) {
    	Wave *w = new EquivalentSource(wave_lenghts[wl], ampli_steps[wl]);
    	((EquivalentSource*)w)->setAmplitude(0, 0, size_tmp-1);
    	waves[wl].push_back(w);
    	inter_src[wl][i] = ((EquivalentSource*)w);
      }
    }
 

    /* boats */
    // addMovingSource(VEC2(8, 2), VEC2(2, 0), 0, 30);
    // addMovingSource(VEC2(11, 5), VEC2(1, 0), 20, 80);
    // addMovingSource(VEC2(15, 7.0), VEC2(0.5, 0), 50, 170);
#endif

    if (point_source1_) {
      addEqSource(1, 15, 0, height_ampli_);
    }
#ifdef INTERACTIVE_
    if (point_source2_) {
      for (uint wl = 0; wl < nb_wl; ++wl) {
	MovingEquivalentSource*w = new MovingEquivalentSource(wave_lenghts[wl], ampli_steps[wl]);
	for (int i = 1; i < 200; ++i) {
	  w->setPos(10, 10+0.05*i, i);
	  w->setAmplitude(height_ampli_, i);
	}
	moving_waves[wl].push_back(w);
      }
    }
#else 
    if (point_source2_) {
      addEqSource(5, 5, 0, height_ampli_);
    }
#endif
    if (point_source3_) {
      addEqSource(15, 15, 0, height_ampli_);
    }

    if (linear_wave1_) {
      addLinearWave(1.0, 0, 0, height_ampli_);
    }
    if (linear_wave2_) {
      addLinearWave(1.0, 1.0, 0, height_ampli_);
    }
    if (linear_wave3_) {
      addLinearWave(0.0, 1.0, 0, height_ampli_);
    }
    if (user_def_source_) {
      VEC2 pos = viewer2world(user_def_pos_);
      addEqSource(pos(0), pos(1), 0, height_ampli_);
    }
  } else {
#ifndef INTERACTIVE_
    importSurfaceFreq(import_file);
#endif
    
  }
  if (circle_) {
    CircularObstacle *ob = new CircularObstacle(VEC2(15, 15), 1);
    obstacles.push_back(ob);
  } if (square_) {
    SquareObstacle *ob = new SquareObstacle(VEC2(15, 15), 2);
    obstacles.push_back(ob);
  } if (harbour_) {
    TextureObstacle *ob = new TextureObstacle("./Textures/harbour.png", 600, 600, 0.05);
    ob->setPos(15, 15);
    obstacles.push_back(ob);
  } if (line_) {
    TextureObstacle *ob = new TextureObstacle("./Textures/line_hd.png", 300, 300, 0.05);
    ob->setPos(15, 15);
    obstacles.push_back(ob);
  } if (island_) {
    TextureObstacle *ob = new TextureObstacle("./Textures/island.png", 4800, 4800, 0.006125);
    ob->setPos(10, 10);
    obstacles.push_back(ob);
  } if (test_) {
    TextureObstacle *ob = new TextureObstacle("./Textures/test.png", 1200, 1200, 0.025);
    ob->setPos(15, 15);
    obstacles.push_back(ob);
  } if (wavy_) {
    WavyObstacle *ob = new WavyObstacle(1, 1, 2);
    ob->setPos(15, 15);
    obstacles.push_back(ob);
  } 
  time = 0;

  if (!import_) {
    if (!obstacles.empty()) {
      std::list<Obstacle*>::iterator ito;
      for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
	(*ito)->reset(waves,  wave_lenghts);
      }
#ifdef INTERACTIVE_
    }
#else 
    for (uint w = 0; w < nb_wl; ++w) {
      for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
	(*ito)->setAmpliSources(waves[w], w);
      }
    }
  }
#endif
#ifndef USE_CUDA
  setAmpli();
#else
  setCuda();
#ifdef INTERACTIVE_
  setCudaM();
  cuda_surface.setHeight(0);
  updateCuda();
  updateCudaM();
#else
  cuda_surface.setHeight(waves.size());
  updateCuda();
#endif
    
#ifndef INTERACTIVE_
#ifdef PLOT_RESULT
  for (uint w = 0; w < nb_wl; ++w) {
    for (int i = 0; i < n_rows_ - 1; i++) {
      for (int j = 0; j < n_cols_ - 1; j++) {
	ampli_in_re[w](i, j) = cuda_surface.heights[w][2*(i*n_cols_ + j)];
	ampli_in_im[w](i, j) =  cuda_surface.heights[w][2*(i*n_cols_ + j)+1];
	ampli_scattered_re[w](i, j) = cuda_surface.heights[w][2*n_cols_*n_rows_ + 2*(i*n_cols_ + j)];
	ampli_scattered_im[w](i, j) = cuda_surface.heights[w][2*n_cols_*n_rows_ + 2*(i*n_cols_ + j)+1];
      }
    }
  }
#endif
#endif
#endif
}
setObstacleGrid();
#ifndef INTERACTIVE_
if (export_) {
  exportSurfaceFreq(export_file);
 }
#endif
#ifdef PLOT_RESULT
std::stringstream ss; 
ss <<data_file<<"ampli.txt";
std::string str(ss.str());
std::stringstream ss2; 
ss2 <<data_file<<"ampli_in.txt";
std::string str2(ss2.str());
std::stringstream ss3; 
ss3 <<data_file<<"ampli_sc.txt";
std::string str3(ss3.str());
std::stringstream ss4; 
ss4 <<data_file<<"phase.txt";
std::string str4(ss4.str());
std::stringstream ss5; 
ss5 <<data_file<<"ampli_dir.txt";
std::string str5(ss5.str());
std::stringstream ss6; 
ss6 <<data_file<<"analytic.txt";
std::string str6(ss6.str());

exportAmplitude(str);
exportAmplitudeIn(str2);
exportAmplitudeScattered(str3);
exportPhase(str4);
exportDirectivityAmplitude(str5, 5);
exportDirectivityAnalytic(str6, 1, 5);
#endif

INFO("oooone update aaaand...");
update();
INFO("DONE!   "<<nb_wl);
}

void WaterSurface::setLists() {
  std::list<Wave*>::iterator it;
  wave_lenghts = std::vector<FLOAT>();
  ampli_steps = std::vector<int>();

  FLOAT wl = min_wl;
  nb_wl = 1;
  wave_lenghts.push_back(wl);
  while (wl < max_wl) {
    if (random_) {
      wl *= (1-0.3*(0.5-(FLOAT)rand()/(FLOAT)(RAND_MAX)))*step_wl;
    } else {
      wl *= step_wl;
    }
    wave_lenghts.push_back(wl);
    ++nb_wl;
  }
   
  ampli_steps = std::vector<int>(nb_wl);
  wl = wave_lenghts[0];
  FLOAT period = 0.5*wl/velocity(2*M_PI/wl);
  int d_period = period/(dt_);
  if (d_period == 0) {
    ampli_steps[0] = 1;
  } else {
    ampli_steps[0] = d_period;
  }
  for (uint w = 1; w < nb_wl; ++w) {
    wl = wave_lenghts[w];
    period = 0.5*wl/velocity(2*M_PI/wl);
    d_period = floor(period/(dt_*ampli_steps[0]))+1;
    if (d_period == 0) {
      ampli_steps[w] = 1;
    } else {
      ampli_steps[w] = d_period*ampli_steps[0];
    }
  }
  
  waves = std::vector<std::list<Wave*> >(nb_wl);
  moving_waves = std::vector<std::list<MovingEquivalentSource*> >(nb_wl);
#ifdef INTERACTIVE_
  inter_src = std::vector<std::vector<EquivalentSource*> >(nb_wl);
  index_inter_src = std::vector<int>(nb_wl);
  for (uint wl = 0; wl < nb_wl; ++wl) {
    index_inter_src[wl] = 0;
  }
#endif
  ampli_in_re = std::vector<Grid>(nb_wl);
  ampli_in_im = std::vector<Grid>(nb_wl);
  ampli_scattered_re = std::vector<Grid>(nb_wl);
  ampli_scattered_im = std::vector<Grid>(nb_wl);
  for (uint i = 0; i < nb_wl; ++i) {
    ampli_in_re[i] = Grid(n_rows_, n_cols_, cell_size_);
    ampli_in_im[i] = Grid(n_rows_, n_cols_, cell_size_);
    ampli_scattered_re[i] = Grid(n_rows_, n_cols_, cell_size_);
    ampli_scattered_im[i] = Grid(n_rows_, n_cols_, cell_size_);
  }
}

void WaterSurface::setAmpli() {
  std::list<Wave*>::iterator it;

  for (uint w = 0; w < nb_wl; ++w) {
    ampli_in_re[w].reset(0);
    ampli_in_im[w].reset(0);
    ampli_scattered_re[w].reset(0);
    ampli_scattered_im[w].reset(0);
    std::list<Wave*> waves_wl = waves[w];
    waves_wl.insert(waves_wl.end(), moving_waves[w].begin(), moving_waves[w].end());
    for (it = waves_wl.begin(); it != waves_wl.end(); ++it) {
#pragma omp parallel for
      for (int i = 0; i < n_rows_ - 1; i++) {
	FLOAT x = cell_size_*i;
	for (int j = 0; j < n_cols_ - 1; j++) {
	  FLOAT y = cell_size_*j;
	  ampli_in_re[w](i, j) += real((*it)->heightc(x, y, 0));
	  ampli_in_im[w](i, j) += imag((*it)->heightc(x, y, 0));
	}
      }
    }
    
    std::list<Obstacle*>::iterator ito;
    for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
#pragma omp parallel for
      for (int i = 0; i < n_rows_ - 1; i++) {
	FLOAT x = cell_size_*i;
	for (int j = 0; j < n_cols_ - 1; j++) {
	  FLOAT y = cell_size_*j;
  	  ampli_scattered_re[w](i, j) += real((*ito)->heightc_wl(x, y, w));
 	  ampli_scattered_im[w](i, j) += imag((*ito)->heightc_wl(x, y, w));
 	}
      }
    }
  }
}

void WaterSurface::setAmpli(FLOAT t) {
  std::list<Wave*>::iterator it;
  for (uint w = 0; w < nb_wl; ++w) {
    ampli_in_re[w].reset(0);
    ampli_in_im[w].reset(0);
    ampli_scattered_re[w].reset(0);
    ampli_scattered_im[w].reset(0);

    std::list<Wave*> waves_wl = waves[w];
    waves_wl.insert(waves_wl.end(), moving_waves[w].begin(), moving_waves[w].end());
    for (it = waves_wl.begin(); it != waves_wl.end(); ++it) {
#pragma omp parallel for
      for (int i = 0; i < n_rows_ - 1; i++) {
	FLOAT x = cell_size_*i;
	for (int j = 0; j < n_cols_ - 1; j++) {
	  FLOAT y = cell_size_*j;
	  ampli_in_re[w](i, j) += real((*it)->heightc(x, y, t));
	  ampli_in_im[w](i, j) += imag((*it)->heightc(x, y, t));
	}
      }
    }
    
    std::list<Obstacle*>::iterator ito;
    for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
#pragma omp parallel for
      for (int i = 0; i < n_rows_ - 1; i++) {
	FLOAT x = cell_size_*i;
	for (int j = 0; j < n_cols_ - 1; j++) {
	  FLOAT y = cell_size_*j;
	  ampli_scattered_re[w](i, j) += real((*ito)->heightc_wl(x, y, w, t));
	  ampli_scattered_im[w](i, j) += imag((*ito)->heightc_wl(x, y, w, t));
 	}
      }
    }
  }
}
 
void WaterSurface::setObstacleGrid() {
  std::list<Obstacle*>::iterator ito;
#pragma omp for
  for (int i = 0; i < n_rows_obs; i++) {
    for (int j = 0; j < n_cols_obs; j++) {
      b(i, j) = 0;
    }
  }
#pragma omp for
  for (int i = 0; i < n_rows_obs - 1; i++) {
    for (int j = 0; j < n_cols_obs - 1; j++) {
      for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {

	if ((*ito)->getGrid(i, j) != 0) {
	  b(i, j) = 0.02;
	  for (int k = -1; k <= 1; ++k) {
	    for (int h = -1; h <= 1; ++h) {
	      if ((*ito)->getGrid(i+k, j+h) == 0) {
		b(i, j) = 0.02;
		break;
	      }
	    }
	  }

	} else {
	  bool border = false;
	  for (int k = -1; k <= 1; ++k) {
	    for (int h = -1; h <= 1; ++h) {
	      if ((*ito)->getGrid(i+k, j+h) != 0) {
		b(i, j) = 0.01;
		break;
	      }
	    }
	  }

	}
	
      }
    }
  }

}

FLOAT WaterSurface::height(uint i, uint j) const {
  return u(i, j);
}

FLOAT WaterSurface::height_obs(uint i, uint j) const {
  if (b(i, j) == 0) {
    return -1;
  }
  return b(i, j);
}

void WaterSurface::addEqSource(FLOAT x, FLOAT y, FLOAT wl, COMPLEX ampli) {
  if (wl == 0) {
    wl = wave_lenghts[nb_wl-1];
  }
  for (uint ind = 0; ind < nb_wl && wave_lenghts[ind]<=wl; ++ind) {
    Wave *w = new EquivalentSource(wave_lenghts[ind], ampli_steps[ind]);
    ((EquivalentSource*)w)->setPos(x, y);
    ((EquivalentSource*)w)->setAmplitude(ampli, 1, size_tmp-1);
    waves[ind].push_back(w);
  }
}

void WaterSurface::addLinearWave(FLOAT x, FLOAT y, FLOAT wl, COMPLEX ampli) {
  FLOAT n = sqrt(x*x + y*y);
  FLOAT s = 0;
  if (wl == 0) {
    wl = wave_lenghts[nb_wl-1];
  }
  for (uint index = 0; index < nb_wl && wave_lenghts[index]<=wl; ++index) {
    s += wave_lenghts[index];
  }
  for (uint index = 0; index < nb_wl && wave_lenghts[index]<=wl; ++index) {
    FLOAT eps_x = 0;
    FLOAT eps_y = 0;
    if (random_) {
      eps_x = 0.25*n*(FLOAT)rand()/(FLOAT)(RAND_MAX);
      eps_y = 0.25*n*(FLOAT)rand()/(FLOAT)(RAND_MAX);
    }
    Wave *w = new LinearWave(wave_lenghts[index]);
    ((LinearWave*)w)->setDir(x + eps_x, y + eps_y);
    ((LinearWave*)w)->setAmplitude(wave_lenghts[index]*ampli/(FLOAT)s);
    waves[index].push_back(w);
  }
}


void WaterSurface::applyForce(FLOAT x, FLOAT y, FLOAT radius) {
  VEC2 pos = viewer2world(x, y);
  if (!inter_src.empty()) {
    for (uint i = 0; i < nb_wl; ++i) {
      inter_src[i][index_inter_src[i]]->setPos(pos);
      if (time >= ampli_steps[i]) {
	
	inter_src[i][index_inter_src[i]]->setAmplitude(wave_lenghts[i]*height_ampli_/(FLOAT)nb_wl, time/ampli_steps[i], time/ampli_steps[i]+2);
      } else {
	inter_src[i][index_inter_src[i]]->setAmplitude(wave_lenghts[i]*height_ampli_/(FLOAT)nb_wl, 1, 1);
      }
      ++index_inter_src[i];
      if (index_inter_src[i] > inter_src[i].size()-1) {
	index_inter_src[i] = 0;
      }
    }
  }
}

void WaterSurface::addDrop(FLOAT x, FLOAT y, FLOAT size) {
  VEC2 pos = viewer2world(x, y);
  addDropW(pos, size, time);
}

int WaterSurface::getTDrop(FLOAT size, int t) {
  int t_drop = (floor(t/ampli_steps[0]) + 1)*ampli_steps[0];
  for (int i = nb_wl - 1; i >= 0; --i) {
    FLOAT wl = wave_lenghts[i];
    if (wl <= 1.01*size) {
      t_drop = (floor(t/ampli_steps[i]) + 1)*ampli_steps[i];
      break;
    }
  }
  return t_drop;
}

void WaterSurface::addDropW(VEC2 pos, FLOAT size, int t) {
#ifdef INTERACTIVE_
  FLOAT size_r = ((FLOAT)rand()/(FLOAT)(RAND_MAX) + 0.5)*size;
  int t_drop = (floor(t/ampli_steps[0]) + 1)*ampli_steps[0];
  for (int i = nb_wl - 1; i >= 0; --i) {
    FLOAT wl = wave_lenghts[i];
    if (wl <= 1.01*size_r) {
      t_drop = (floor(t/ampli_steps[i]) + 1)*ampli_steps[i];
      break;
    }
  }
  for (int i = 0; i < nb_wl; ++i) {
    FLOAT wl = wave_lenghts[i];
    if (!inter_src[i].empty()) {
      if (wl <= 1.01*size_r) {
	FLOAT period = wl/velocity(2*M_PI/wl);
	int d_period = 2*period/(dt_*ampli_steps[i]);
	inter_src[i][index_inter_src[i]]->setPos(pos);
	if (t_drop >= ampli_steps[i]) {
	  inter_src[i][index_inter_src[i]]->setAmplitude(wl*height_ampli_, t_drop/ampli_steps[i], t_drop/ampli_steps[i]+d_period);
	} else {
	  inter_src[i][index_inter_src[i]]->setAmplitude(wl*height_ampli_, 1, 1+d_period);
	}
	++index_inter_src[i];
	if (index_inter_src[i] > inter_src[i].size()-1) {
	  index_inter_src[i] = 0;
	}
      }
    }
  }
#endif
}

void updateSourcesAmplis(WaterSurface *ws) {
  uint t_cur = ws->time; 
  Times::TIMES_UP->next_loop();
  Times::TIMES_UP->tick(Times::solve_time_);
  FLOAT t = t_cur*dt_;

#ifdef INTERACTIVE_
  if(!ws->obstacles.empty()) {
    std::list<Obstacle*>::iterator ito;
    std::list<Obstacle*>::iterator ito2;
    uint n_w = 0;
    for (int w = ws->nb_wl - 1; w >= 0; --w) {
      if ((!ws->waves[w].empty() || !ws->moving_waves[w].empty()) && t_cur%ampli_steps[w] == 0) {
	FLOAT wl = ws->wave_lenghts[w];
	FLOAT k = 2*M_PI/wl;
	if (wl/2.0 > ws->proj_grid.minSize()) {
	  uint n_o = 0;
	  for (ito = ws->obstacles.begin(); ito != ws->obstacles.end(); ++ito) {
	    VEC2 po = (*ito)->getPos();
	    std::list<Wave*> all_waves = ws->waves[w];
	    all_waves.insert(all_waves.end(), ws->moving_waves[w].begin(), ws->moving_waves[w].end());
	    for (ito2 = ws->obstacles.begin(); ito2 != ws->obstacles.end(); ++ito2) {
	      if (ito != ito2) {
		FLOAT d = (po - (*ito2)->getPos()).norm();
		if (damping(d, k) > 0.01) {
		  (*ito2)->getSources(all_waves, w);
		}
	      }
	    }
	    (*ito)->setAmpliSources(all_waves, w, t_cur/ampli_steps[w]+1);
	    ++n_o;
	  }
	}
	++n_w;
      }
    }

  }


  
# endif
  std::list<Obstacle*>::const_iterator ito;
  std::list<Wave*>::const_iterator it;
  std::list<MovingEquivalentSource*>::const_iterator itm;
  for (int w =0; w < ws->nb_wl; ++w) {
    if (t_cur%ampli_steps[w] == 0) {
      std::list<Wave*> waves_wl = ws->waves[w];
      std::list<MovingEquivalentSource*> moving_waves_wl = ws->moving_waves[w];
      for (ito = ws->obstacles.begin(); ito != ws->obstacles.end(); ++ito) {
  	(*ito)->getSources(waves_wl, w);
	(*ito)->getSourcesM(moving_waves_wl, w);
      }
      for (it = waves_wl.begin(); it != waves_wl.end(); ++it) {
  	(*it)->setActive(t_cur/ampli_steps[w]);
      }
      for (itm = moving_waves_wl.begin(); itm != moving_waves_wl.end(); ++itm) {
	((EquivalentSource*)(*itm))->setActive(t_cur/ampli_steps[w]);
      }

    }
  }
  Times::TIMES_UP->tock(Times::solve_time_);
}



void WaterSurface::update() {
#ifdef INTERACTIVE_
  if (import_) {
    std::stringstream ss;
    ss <<import_file<<time<<".grid";
    std::string str(ss.str());
    importSurfaceTime(str);
  } else {

#endif
    //    INFO("TIME: "<<time);
    if (time%ampli_steps[0] == 0) {
      /* we solve the amplitude of the sources every <ampli_steps> frames */
      /* this is done on another thread, so we can keep displaying the animation */
      if (t_solve.joinable()) {
       	t_solve.join();
      }
      t_solve = std::thread(updateSourcesAmplis, this);
      //       updateSourcesAmplis(this);
    }
  
    Times::TIMES->tick(Times::sum_up_time_);    
    updateHeight();
    Times::TIMES->tock(Times::sum_up_time_);   
#ifdef INTERACTIVE_
  }
  if (export_) {
    if (time%export_step == 0) {
      std::stringstream ss; 
      ss <<export_file<<time/export_step<<".grid";
      std::string str(ss.str());
      exportSurfaceTime(str);
    }

  }

#endif

  ++time;
  
  if (time > stop_time) {
    exit(0);
  }
}



void WaterSurface::updateHeight() {
  u.reset(0.0);
  u_re.reset(0.0);
  u_im.reset(0.0);
  FLOAT t = time*dt_;
#ifndef USE_CUDA
  setAmpli(t);
  for (uint f = 0; f < nb_wl; ++f) {
    FLOAT k =  2*M_PI/wave_lenghts[f];
    FLOAT omega = angular_vel(k);
    if (show_in_field) {
#pragma omp parallel for
      for (int i = 0; i < n_rows_ - 1; i++) {
	for (int j = 0; j < n_cols_ - 1; j++) {
	  u(i, j) += real(COMPLEX(ampli_in_re[f](i, j), ampli_in_im[f](i, j))*exp(-omega*t*i_));
	}
      }
    }
       
    std::list<Obstacle*>::iterator ito;
       
    if (show_scattered_field) {
#pragma omp parallel for
      for (int i = 0; i < n_rows_ - 1; i++) {
	for (int j = 0; j < n_cols_ - 1; j++) {
	  u(i, j) += real(COMPLEX(ampli_scattered_re[f](i, j), ampli_scattered_im[f](i, j))*exp(-omega*t*i_));
	}
      }
    }
  }

#else
#ifdef INTERACTIVE_
  if (time%ampli_steps[0] == 0) {
    updateCuda();
    updateCudaM();
  }
  cuda_surface.setHeight(time);
 
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n_rows_ - 1; i++) {
    for (int j = 0; j < n_cols_ - 1; j++) {
#ifdef PROJECTED_GRID
      proj_grid.setHeight(i, j, cuda_surface.heights[i*n_cols_ + j]);
      IS_DEF(cuda_surface.heights[i*n_cols_ + j]);
#else
      u(i, j) = cuda_surface.heights[i*n_cols_ + j];
#ifdef PLOT_RESULT
      ampli_scattered_re[0](i, j) = cuda_surface.heights[n_cols_*n_rows_ + i*n_cols_ + j];
#endif
#endif
    }
  }
#else
  updateCuda();
  //  updateCudaM();
#pragma omp parallel for
  for (int i = 0; i < n_rows_ - 1; i++) {
    for (int j = 0; j < n_cols_ - 1; j++) {
#ifdef PROJECTED_GRID
      proj_grid.setHeight(i, j, cuda_surface.time_heights[i*n_cols_ + j]);
      proj_grid.setDisplacement(i, j, -cuda_surface.time_displacement[2*(i*n_cols_ + j)],
				-cuda_surface.time_displacement[2*(i*n_cols_ + j)+1]);
#else
      u(i, j) = cuda_surface.time_heights[i*n_cols_ + j];
      dispx(i, j) = -cuda_surface.time_displacement[2*(i*n_cols_ + j)];
      dispy(i, j) = -cuda_surface.time_displacement[2*(i*n_cols_ + j)+1];
#endif
    }
  }

#endif
#endif
}



 
void WaterSurface::getObstacleBoundary(std::list<VEC2> &boundaries) const {
  std::list<Obstacle*>::const_iterator ito;
  for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
    (*ito)->getBoundariesPos(boundaries);
  }
}


void WaterSurface::exportAmplitude(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a =	sqrt(pow(ampli_in_re[0](i, j) + ampli_scattered_re[0](i, j), 2) +
		     pow(ampli_in_im[0](i, j) + ampli_scattered_im[0](i, j), 2));
      if (b(i, j) != 0 || a < 0.001) {
	a = 0.0;
      } else {
	a =20*log10(a/0.001);
      }
      out_file<<i<<" "<<j<<" "<<a<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}

void WaterSurface::exportAmplitudeRe(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = ampli_in_re[0](i, j) + ampli_scattered_re[0](i, j);
      if (b(i, j) == 1) {
	a = 0;
      }
      out_file<<i<<" "<<j<<" "<<a<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}

void WaterSurface::exportAmplitudeIm(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = ampli_in_im[0](i, j) + ampli_scattered_im[0](i, j);
      if (b(i, j) == 1) {
	a = 0;
      }
      out_file<<i<<" "<<j<<" "<<a<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}

void WaterSurface::exportPhase(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = sqrt(pow(ampli_in_re[0](i, j) + ampli_scattered_re[0](i, j), 2) +  pow(ampli_in_im[0](i, j) + ampli_scattered_im[0](i, j), 2));
      if (b(i, j) == 1) {
	a = 0;
      }
      out_file<<i<<" "<<j<<" "<<acos((ampli_in_re[0](i, j) + ampli_scattered_re[0](i, j))/a)<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}

void WaterSurface::exportAmplitudeIn(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = sqrt(pow(ampli_in_re[0](i, j), 2) +  pow(ampli_in_im[0](i, j), 2));
      if (b(i, j) != 0 || a < 0.0001) {
	a = 0.0;
      } else {
	a =20*log10(a/0.0001);
      }
      out_file<<i<<" "<<j<<" "<<a<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}

void WaterSurface::exportAmplitudeInRe(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = ampli_in_re[0](i, j);
      if (b(i, j) == 1) {
	a = 0;
      }
      out_file<<i<<" "<<j<<" "<<a<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}

void WaterSurface::exportAmplitudeInIm(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = ampli_in_im[0](i, j);
      if (b(i, j) == 1) {
	a = 0;
      }
      out_file<<i<<" "<<j<<" "<<a<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}

void WaterSurface::exportPhaseIn(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = sqrt(pow(ampli_in_re[0](i, j), 2) +  pow(ampli_in_im[0](i, j), 2));
      if (b(i, j) == 1) {
	a = 0;
      }
      out_file<<i<<" "<<j<<" "<<acos(ampli_in_re[0](i, j)/a)<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}

void WaterSurface::exportAmplitudeScattered(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = sqrt(pow(ampli_scattered_re[0](i, j), 2) +  pow(ampli_scattered_im[0](i, j), 2));
      if (b(i, j) != 0 || a < 0.0001) {
	a = 0.0;
      } else {
	a =20*log10(a/0.0001);
      }
      out_file<<i<<" "<<j<<" "<<a<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();

}

void WaterSurface::exportAmplitudeScatteredRe(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = ampli_scattered_re[0](i, j);
      if (b(i, j) == 1) {
	a = 0;
      }
      out_file<<i<<" "<<j<<" "<<a<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}

void WaterSurface::exportAmplitudeScatteredIm(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = ampli_scattered_im[0](i, j);
      if (b(i, j) == 1) {
	a = 0;
      }
      out_file<<i<<" "<<j<<" "<<a<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}


void WaterSurface::exportPhaseScattered(std::string file) const {
  std::ofstream  out_file;
  out_file.open(file);

  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = sqrt(pow(ampli_scattered_re[0](i, j), 2) +  pow(ampli_scattered_im[0](i, j), 2));
      if (b(i, j) == 1) {
	a = 0;
      }
      out_file<<i<<" "<<j<<" "<<acos(ampli_scattered_re[0](i, j)/a)<<"\n";
    }
    out_file<<"\n";
  }
  out_file.close();
}


void WaterSurface::exportDirectivityAmplitude(std::string file, FLOAT dist) const {
  INFO("exporting "<<file);
  std::ofstream  out_file;
  out_file.open(file);
  uint nb = 1000;
  FLOAT step = 2*M_PI/(float)nb;
  VEC2 dir(1, 0);
  MAT2 rotation;
  rotation <<
    cos(step), -sin(step),
    sin(step), cos(step);
  for (uint i = 0; i <= nb; ++i) {
    VEC2 p = VEC2(15, 15) + dist*dir;
    Vector2i n = world2grid(p);
    FLOAT a = sqrt(pow(ampli_in_re[0](n(0), n(1)) + ampli_scattered_re[0](n(0), n(1)), 2) +
		   pow(ampli_in_im[0](n(0), n(1)) + ampli_scattered_im[0](n(0), n(1)), 2));
    // a =20*log10(a/0.001);
    // if (a < 0) {
    //   a = 0;
    // }
    out_file<<i*step<<" "<<a<<"\n";
    dir = rotation*dir;
  }
  out_file.close();
}

void WaterSurface::exportDirectivityAnalytic(std::string file, FLOAT r_obs, FLOAT r) const {
  INFO("Exporting "<<file);
  std::ofstream  out_file;
  out_file.open(file);
  uint nb = 1000;
  FLOAT step = 2*M_PI/(float)nb;
  FLOAT theta = 0;
  VEC2 dir(1, 0);
  MAT2 rotation;
  rotation <<
    cos(step), -sin(step),
    sin(step), cos(step);
  FLOAT k = 2*M_PI/wave_lenghts[0];
  FLOAT err = 0;
  for (uint i = 0; i <= nb; ++i) {
    COMPLEX a = 0;
    for (int m = -45; m <= 45; ++m) {
      FLOAT Jmkr = boost::math::cyl_bessel_j<FLOAT, FLOAT>(m, k*r);
      FLOAT Ymkr = boost::math::cyl_neumann<FLOAT, FLOAT>(m, k*r);
      if (neumann) {
	FLOAT Jmka_der = 0.5*(boost::math::cyl_bessel_j<FLOAT, FLOAT>(m-1, k*r_obs) -
			      boost::math::cyl_bessel_j<FLOAT, FLOAT>(m+1, k*r_obs));
	FLOAT Ymka_der = 0.5*(boost::math::cyl_neumann<FLOAT, FLOAT>(m-1, k*r_obs) -
			      boost::math::cyl_neumann<FLOAT, FLOAT>(m+1, k*r_obs));
	a += pow(i_, m)*
	  (Jmkr*Ymka_der - Jmka_der*Ymkr)/(Jmka_der+i_*Ymka_der)*exp(i_*(FLOAT)m*theta);
      } else {
	FLOAT Jmka = boost::math::cyl_bessel_j<FLOAT, FLOAT>(m, k*r_obs);
	FLOAT Ymka = boost::math::cyl_neumann<FLOAT, FLOAT>(m, k*r_obs);
	
	a += pow(i_, m)*
	  (Jmkr*Ymka - Jmka*Ymkr)/(Jmka+i_*Ymka)*exp(i_*(FLOAT)m*theta);
      }
    }
    a *= i_;
    FLOAT aa = 0.1*height_ampli_*sqrtf(pow(real(a), 2) + pow(imag(a), 2));
    // aa =20*log10(aa/0.001);
    // if (aa < 0) {
    //   aa = 0;
    // }
    out_file<<theta<<" "<<aa<<"\n";
    theta += step;
    
    /* compute error compared between analytical sol and simu */
    // VEC2 p = VEC2(15, 15) + r*dir;
    // Vector2i n = world2grid(p);
    
    // FLOAT a_sim = sqrt(pow(ampli_in_re[0](n(0), n(1)) + ampli_scattered_re[0](n(0), n(1)), 2) +
    // 		       pow(ampli_in_im[0](n(0), n(1)) + ampli_scattered_im[0](n(0), n(1)), 2));
    // // a =20*log10(a/0.001);
    // // if (a < 0) {
    // //   a = 0;
    // // }
    // dir = rotation*dir;
    // err += fabs(a_sim - aa)/aa;
    // INFO("err "<<a_sim<<" "<<aa<<" "<<fabs(a_sim - aa)<<" "<<i);
  }
  // err /= nb;
  // INFO("MEAN ERR"<<err);
  out_file.close();
}


void WaterSurface::exportSurfaceFreq(std::string file) const {
  std::cout<<"Export freq "<<file<<std::endl;
  std::ofstream os(file.c_str());

  ERROR(os.good(), "Cannot open file "<<file, "");

  os<<nb_wl<<"\n";

  os<<"f ";
  for (uint i = 0; i < nb_wl; ++i) {
    os<<wave_lenghts[i]<<" ";
  }
  os<<"\n";
  os<<"gir\n";
  for (uint i = 0; i < nb_wl; ++i) {
    os<<ampli_in_re[i]<<"\n";
  }
  os<<"gii\n";
  for (uint i = 0; i < nb_wl; ++i) {
    os<<ampli_in_im[i]<<"\n";
  }
  os<<"gsr\n";
  for (uint i = 0; i < nb_wl; ++i) {
    os<<ampli_scattered_re[i]<<"\n";
  }
  os<<"gsi\n";
  for (uint i = 0; i < nb_wl; ++i) {
    os<<ampli_scattered_im[i]<<"\n";
  }
  os.close();
}
 
void WaterSurface::importSurfaceFreq(std::string file) {
  std::cout<<"Import freq "<<file<<std::endl;
  std::ifstream is(file.c_str());

  ERROR(is.good(), "Cannot open file "<<file, "");
  
  is>>nb_wl;
  wave_lenghts = std::vector<FLOAT>(nb_wl);

  ampli_in_re = std::vector<Grid>(nb_wl);
  ampli_in_im = std::vector<Grid>(nb_wl);
  ampli_scattered_re = std::vector<Grid>(nb_wl);
  ampli_scattered_im = std::vector<Grid>(nb_wl);
  for (uint i = 0; i < nb_wl; ++i) {
    ampli_in_re[i] = Grid(n_rows_, n_cols_, cell_size_);
    ampli_in_im[i] = Grid(n_rows_, n_cols_, cell_size_);
    ampli_scattered_re[i] = Grid(n_rows_, n_cols_, cell_size_);
    ampli_scattered_im[i] = Grid(n_rows_, n_cols_, cell_size_);
  }
  std::string line;
  while (getline(is, line)) {
    if (line.substr(0,2) == "f ") {
      std::istringstream s(line.substr(2));
      for (uint i = 0; i < nb_wl; ++i) {
	s>>wave_lenghts[i];
      }
    } else if (line.substr(0,3) == "gir") {
      for (uint i = 0; i < nb_wl; ++i) {
	is>>ampli_in_re[i];
      }
      getline(is, line);
    } else if (line.substr(0,3) == "gii") {
      for (uint i = 0; i < nb_wl; ++i) {
	is>>ampli_in_im[i];
      }
      getline(is, line);
    } else if (line.substr(0,3) == "gsr") {
      for (uint i = 0; i < nb_wl; ++i) {
	is>>ampli_scattered_re[i];
      }
      getline(is, line);
    } else if (line.substr(0,3) == "gsi") {
      for (uint i = 0; i < nb_wl; ++i) {
	is>>ampli_scattered_im[i];
      }
      getline(is, line);
    } else if (line[0] == '#' || line == "") {
      INFO("COMMENT "<<line);
    } else {
      ERROR(false, "Imported file possibly corrupted", line.substr(0,100));
    }
   
  }
  is.close();
}

void WaterSurface::exportSurfaceTime(std::string file) const {
  std::cout<<"Export time "<<file<<std::endl;
  std::ofstream os(file.c_str());
  ERROR(os.good(), "Cannot open file "<<file, "");
  os<<u;
  os.close();
}
 
void WaterSurface::importSurfaceTime(std::string file) {
  std::ifstream is(file.c_str());
  ERROR(is.good(), "Cannot open file "<<file, "");
  std::cout<<"Import time "<<file<<std::endl;
  is>>u;
  is.close();
}

void WaterSurface::exportSurfaceObj(std::string file) const {
  std::cout<<"Exporting "<<file<<std::endl;
  std::ofstream os(file.c_str());
  ERROR(os.good(), "Cannot open file "<<file, "");
#ifdef PROJECTED_GRID
  proj_grid.exportObj(os);
#else
  u.exportObj(os);
#endif
  os.close();
}

void  WaterSurface::exportObstacleGrid(std::string file) const {
  std::cout<<"Exporting "<<file<<std::endl;
  std::ofstream os(file.c_str());
  ERROR(os.good(), "Cannot open file "<<file, "");
  b.exportObj(os);
  os.close();
}

std::ofstream& WaterSurface::exportObstacleMitsuba(std::ofstream &file) const {
  std::list<Obstacle*>::const_iterator ito;
  for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
    (*ito)->exportMitsuba(file);
  }
  return file;
}

void WaterSurface::importConfig(std::string file) {
  std::ifstream is(file);
  ERROR(is.good(), "Cannot open file "<<file, "");
  std::string line;
  while (getline(is, line)) {
    if (line.substr(0,1) == "#") {
      //comment
    } else if (line.substr(0,11) == "<time_step>") {
      std::istringstream s(line.substr(11));
      s >> dt_ ;
    } else if (line.substr(0,12) == "<ampli_step>") {
      INFO("Configuration file: line \""<<line<<"\" ignored");
    } else if (line.substr(0,9) == "<gravity>") {
      std::istringstream s(line.substr(9));
      s >> gravity_ ;
    } else if (line.substr(0,9) == "<damping>") {
      std::istringstream s(line.substr(9));
      s >> damping_ ;
    } else if (line.substr(0,7) == "<ampli>") {
      std::istringstream s(line.substr(7));
      s >> height_ampli_ ;
    } else if (line.substr(0,16) == "<damping_source>") {
      std::istringstream s(line.substr(16));
      s >> damping_source_ ;
    } else if (line.substr(0,8) == "<approx>") {
      approx_inter = true;
    } else if (line.substr(0,11) == "<no_approx>") {
      approx_inter = false;
    } else if (line.substr(0,8) == "<offset>") {
      std::istringstream s(line.substr(8));
      s >> offset_ ;
    } else if (line.substr(0,10) == "<sampling>") {
      std::istringstream s(line.substr(10));
      s >> step_sampling_;
    } else if (line.substr(0,9) == "<neumann>") {
      neumann = true;
    } else if (line.substr(0,11) == "<dirichlet>") {
      neumann = false;
    } else if (line.substr(0,6) == "<grid>") { //GRID DEF 
      getline(is, line);
      while (line.substr(0,7) != "</grid>") {
	if (line.substr(0,6) == "<size>") {
	  std::istringstream s(line.substr(6));
	  s >> n_rows_ >> n_cols_;
	} else if (line.substr(0,11) == "<cell_size>") {
	  std::istringstream s(line.substr(11));
	  s >> cell_size_;
	} else {
	  ERROR(false, "Invalid configuration file (grid)\""<<file, line);
	}
	scale_ = cell_size_ * n_rows_;
	getline(is, line);
      }
    } else if (line.substr(0,10) == "<grid_obs>") { //GRID DEF 
      getline(is, line);
      while (line.substr(0,11) != "</grid_obs>") {
	if (line.substr(0,6) == "<size>") {
	  std::istringstream s(line.substr(6));
	  s >> n_rows_obs >> n_cols_obs;
	} else if (line.substr(0,11) == "<cell_size>") {
	  std::istringstream s(line.substr(11));
	  s >> cell_size_obs;
	} else {
	  ERROR(false, "Invalid configuration file (grid_obs)\""<<file, line);
	}
	getline(is, line);
      }
    } else if (line.substr(0,14) == "<wave_lenghts>") { //WAVE LENGHT RANGE 
      getline(is, line);
      while (line.substr(0,15) != "</wave_lenghts>") {
	if (line.substr(0,5) == "<min>") {
	  std::istringstream s(line.substr(5));
	  s >> min_wl;
	} else if (line.substr(0,5) == "<max>") {
	  std::istringstream s(line.substr(5));
	  s >> max_wl;
	} else if (line.substr(0,8) == "<number>") {
	  std::istringstream s(line.substr(8));
	  s >> nb_wl;
	  step_wl = (max_wl - min_wl)/(FLOAT)(nb_wl);
	} else if (line.substr(0,6) == "<step>") {
	  std::istringstream s(line.substr(6));
	  s >> step_wl;
	  nb_wl = (max_wl - min_wl)/(FLOAT)(step_wl);
	} else {
	  ERROR(false, "Invalid configuration file (wave lenghts)\""<<file, line);
	}
	getline(is, line);
      }
      setLists();
	 	 
    } else if (line.substr(0,6) == "<wave>") { //def of a wave
      getline(is, line);
      FLOAT wl = 0;
      uint ind = 0;
      FLOAT ampli = height_ampli_;
      while (line.substr(0,7) != "</wave>") {
	if (line.substr(0,13) == "<wave_lenght>") {
	  std::istringstream s(line.substr(11));
	  s >> wl;
	} else if (line.substr(0,11) == "<amplitude>") {
	  std::istringstream s(line.substr(11));
	  s >> ampli;
	} else if (line.substr(0,7) == "<ampli>") {
	  std::istringstream s(line.substr(7));
	  s >> ampli;
	} else if (line.substr(0,6) == "<line>") {
	  getline(is, line);
	  VEC2 dir;
	  while (line.substr(0,7) != "</line>") {
	    if  (line.substr(0,5) == "<dir>") {
	      std::istringstream s(line.substr(6));
	      s >> dir(0) >> dir(1);
	    } else {
	      ERROR(false, "Invalid configuration file (line wave)\""<<file, line);
	    }
	    getline(is, line);
	  }
	  addLinearWave(dir(0), dir(1), wl, ampli);
	     
	} else if (line.substr(0,8) == "<source>") {
	  getline(is, line);
	  VEC2 pos;
	  while (line.substr(0,9) != "</source>") {
	    if  (line.substr(0,5) == "<pos>") {
	      std::istringstream s(line.substr(6));
	      s >> pos(0) >> pos(1);
	    } else {
	      ERROR(false, "Invalid configuration file (source wave)\""<<file, line);
	    }
	    getline(is, line);
	  }
	  addEqSource(pos(0), pos(1), wl, ampli);
	} else if (line.substr(0,10) == "<line_all>") {
	  getline(is, line);
	  VEC2 dir;
	  while (line.substr(0,11) != "</line_all>") {
	    if  (line.substr(0,5) == "<dir>") {
	      std::istringstream s(line.substr(6));
	      s >> dir(0) >> dir(1);
	      dir.normalize();
	    } else {
	      ERROR(false, "Invalid configuration file (line all wave)\""<<file, line);
	    }
	    getline(is, line);
	  }
	  addLinearWave(dir(0), dir(1), 0, ampli);
	} else if (line.substr(0,12) == "<source_all>") {
	  getline(is, line);
	  VEC2 pos;
	  while (line.substr(0,13) != "</source_all>") {
	    if  (line.substr(0,5) == "<pos>") {
	      std::istringstream s(line.substr(6));
	      s >> pos(0) >> pos(1);
	    } else {
	      ERROR(false, "Invalid configuration file (spurce all wave)\""<<file, line);
	    }
	    getline(is, line);
	  }
	  for (uint i = 0; i < nb_wl; ++i) {
	    EquivalentSource *w = new EquivalentSource(wave_lenghts[i], ampli_steps[i]);
	    w->setPos(pos);
	    w->setAmplitude(ampli);
	    waves[i].push_back(w);
	  }
	} else {
	  ERROR(false, "Invalid configuration file (wave)\""<<file, line);
	}
	getline(is, line);
      }

    } else if (line.substr(0,10) == "<obstacle>") { //def of a obstacle
      getline(is, line);
      VEC2 pos(15, 15);
      while (line.substr(0,11) != "</obstacle>") {
	if (line.substr(0,8) == "<circle>") {
	  getline(is, line);
	  VEC2 pos(0, 0);
	  FLOAT radius = 1;
	  while (line.substr(0,9) != "</circle>") {
	    if  (line.substr(0,5) == "<pos>") {
	      std::istringstream s(line.substr(5));
	      s >> pos(0) >> pos(1);
	    } else if (line.substr(0,8) == "<radius>") {
	      std::istringstream s(line.substr(8));
	      s >> radius;
	    } else {
	      ERROR(false, "Invalid configuration file (circle obstacle)\""<<file, line);
	    }
	    getline(is, line);
	  }
	  CircularObstacle *ob = new CircularObstacle(pos, radius);
	  obstacles.push_back(ob);
	} else  if (line.substr(0,8) == "<square>") {
	  getline(is, line);
	  VEC2 pos(0, 0);
	  FLOAT size = 1;
	  while (line.substr(0,9) != "</square>") {
	    if  (line.substr(0,5) == "<pos>") {
	      std::istringstream s(line.substr(5));
	      s >> pos(0) >> pos(1);
	    } else if (line.substr(0,6) == "<size>") {
	      std::istringstream s(line.substr(6));
	      s >> size;
	    } else {
	      ERROR(false, "Invalid configuration file (dquare obstacle)\""<<file, line);
	    }
	    getline(is, line);
	  }
	  SquareObstacle *ob = new SquareObstacle(pos, size);
	  ob->setPos(pos);
	  obstacles.push_back(ob);
	} else  if (line.substr(0,9) == "<texture>") {
	  getline(is, line);
	  std::string name;
	  uint nc = n_cols_;
	  uint nr = n_rows_;
	  FLOAT cs = cell_size_;
	  while (line.substr(0,10) != "</texture>") {
	    if  (line.substr(0,6) == "<file>") {
	      file = line.substr(7);
	    } else if (line.substr(0,6) == "<size>") {
	      std::istringstream s(line.substr(6));
	      s >> nr >> nc;
	    } else if (line.substr(0,11) == "<cell_size>") {
	      std::istringstream s(line.substr(11));
	      s >> cs;
	    } else {
	      ERROR(false, "Invalid configuration file (texture obstacle)\""<<file, line);
	    }
	    getline(is, line);
	  }
	  TextureObstacle *ob = new TextureObstacle(file, nr, nc, cs);
	  ob->setPos(pos);
	  obstacles.push_back(ob);
	} else  if (line.substr(0,6) == "<wavy>") {
	  getline(is, line);
	  uint pattern = 0;
	  FLOAT freq = 1;
	  FLOAT ampli = 0.5;
	  while (line.substr(0,7) != "</wavy>") {
	    if  (line.substr(0,9) == "<pattern>") {
	      std::istringstream s(line.substr(9));
	      s >> pattern;
	    } else if (line.substr(0,11) == "<frequence>") {
	      std::istringstream s(line.substr(11));
	      s >> freq;
	    } else if (line.substr(0,7) == "<ampli>") {
	      std::istringstream s(line.substr(7));
	      s >> ampli;
	    } else {
	      ERROR(false, "Invalid configuration file (wavy obstacle)\"", line);
	    }
	    getline(is, line);
	  }
	  WavyObstacle *ob = new WavyObstacle(pattern, freq, ampli);
	  ob->setPos(pos);
	  obstacles.push_back(ob);
	} else  if (line.substr(0,5) == "<pos>") {
	  std::istringstream s(line.substr(5));
	  s >> pos(0) >> pos(1);
	} else {
	  ERROR(false, "Invalid configuration file (obstacle)\""<<file, line);
	}
	getline(is, line);
      }
    } else {
      ERROR(false, "Invalid configuration file \""<<file, line);
    }
  }
  is.close();
}



void WaterSurface::setImport(std::string file) {
  import_ = true;
  import_file = file;
}

void WaterSurface::setExport(std::string file) {
  export_ = true;
  export_file = file;
}

void WaterSurface::setExportMitsuba(std::string file) {
  export_mit = true;
  export_mit_file = file;
}

void WaterSurface::setImportConf(std::string file) {
  load_conf = true;
  conf_file = file;
}

void WaterSurface::setExportStep(uint es) {
  export_step = es;
}

void WaterSurface::setStopTime(uint end) {
  stop_time = end;
}

void WaterSurface::setData(std::string file) {
  data_file = file;
}

void WaterSurface::drawHeighField(std::string file) {
  std::ofstream  out_file;
  out_file.open(file);
  INFO("Exporting "<<file);
  out_file<<0<<" "<<0<<" "<<0.25*height_ampli_<<"\n";
  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      FLOAT a = u(i, j);
      if (b(i*(cell_size_/cell_size_obs), j*(cell_size_/cell_size_obs)) != 0/* || a > 0.01*/) {
	a = -0.25*height_ampli_;
      }
      if (a < -0.25*height_ampli_) {
	a = -0.25*height_ampli_;
      }
      if (a > 0.25*height_ampli_) {
	a = 0.25*height_ampli_;
      }
      if (i != 0 && j != 0) {
	out_file<<i<<" "<<j<<" "<<a<<"\n";
      }
    }
    out_file<<"\n";
  }
  out_file.close();
}

#ifdef USE_CUDA
#ifdef INTERACTIVE_

void WaterSurface::setCuda() {
  cuda_surface.init(nb_wl);
  
  std::list<Obstacle*>::const_iterator ito;
  std::list<Wave*>::const_iterator it;
  for (uint w = 0; w < nb_wl; ++w) {
    std::list<Wave*> waves_wl;
    std::list<Wave*> waves_wl_input = waves[w];
    for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
      (*ito)->getSources(waves_wl, w);
    }
    uint nb_sources = waves_wl.size();
    uint nb_sources_input = waves_wl_input.size();
    cuda_surface.allocMem(w, nb_sources, nb_sources_input);
    uint s = 0;
    cuda_surface.wave_lenghts[w] = wave_lenghts[w];

    for (it = waves_wl_input.begin(); it != waves_wl_input.end(); ++it, ++s) {
      COMPLEX a = (*it)->getAmpli(0);
      cuda_surface.amplitudes[w][2*s*size_tmp] = real(a);
      cuda_surface.amplitudes[w][2*s*size_tmp+1] = imag(a);
      cuda_surface.indexes[w][s] = (*it)->getIndex();
      VEC2 pos = (*it)->getPos();
      cuda_surface.positions[w][2*s] = pos(0);
      cuda_surface.positions[w][2*s+1] = pos(1);
      cuda_surface.is_active[w][s] = (*it)->isActive();
    }
    for (it = waves_wl.begin(); it != waves_wl.end(); ++it, ++s) {
      COMPLEX a = (*it)->getAmpli(0);
      cuda_surface.amplitudes[w][2*s*size_tmp] = real(a);
      cuda_surface.amplitudes[w][2*s*size_tmp+1] = imag(a);
      cuda_surface.indexes[w][s] = (*it)->getIndex();
      VEC2 pos = (*it)->getPos();
      cuda_surface.positions[w][2*s] = pos(0);
      cuda_surface.positions[w][2*s+1] = pos(1);

      cuda_surface.is_active[w][s] = (*it)->isActive();
    }
  }
#ifdef PROJECTED_GRID
  proj_grid.setSizes();
  for (uint i = 0; i < n_rows_*n_cols_; ++i) {
    VEC3 p = proj_grid(i);
    cuda_surface.positions_grid[2*i] = p(0);
    cuda_surface.positions_grid[2*i+1] = p(1);
    cuda_surface.sizes[i] = proj_grid.size(i);
  }
#endif
}

void WaterSurface::updateCuda() {
  std::list<Obstacle*>::const_iterator ito;
  std::list<Wave*>::const_iterator it;
  for (int w = nb_wl-1; w >= 0; --w) {
    if (wave_lenghts[w]/2.0 <= proj_grid.minSize()) {
      break;
    } else {
      std::list<Wave*> waves_wl;
      std::list<Wave*> waves_wl_input = waves[w];
      uint nb_in_wave = waves_wl.size();
      for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
	(*ito)->getSources(waves_wl, w);
      }
      uint s = 0;
      uint as = 0;
     
      for (it = waves_wl_input.begin(); it != waves_wl_input.end(); ++it, ++s) {
	COMPLEX a(0, 0);
	if (time/ampli_steps[w] > 0) {
	  int ta = (time/ampli_steps[w])%size_tmp+1;
	  a = (*it)->getAmpli(time/ampli_steps[w]);
	  if (real(a) != 0 || imag(a) != 0) {
	    cuda_surface.amplitudes[w][2*(s*size_tmp+ta)] = real(a);
	    cuda_surface.amplitudes[w][2*(s*size_tmp+ta)+1] = imag(a);
	  }
	}
	cuda_surface.is_active[w][s] = (*it)->isActive();
	cuda_surface.positions[w][2*s] = (*it)->getPos()[0];
	cuda_surface.positions[w][2*s+1] = (*it)->getPos()[1];
      }
      if (wave_lenghts[w]/2.0 > proj_grid.minSize()) {
	for (it = waves_wl.begin(); it != waves_wl.end(); ++it, ++s) {
	  if (time/ampli_steps[w] > 0) {
	    int ta = (time/ampli_steps[w])%size_tmp+1;
	    COMPLEX a = (*it)->getAmpli(time/ampli_steps[w]);
	    if (real(a) != 0 || imag(a) != 0) {
	      cuda_surface.amplitudes[w][2*(s*size_tmp+ta)] = real(a);
	      cuda_surface.amplitudes[w][2*(s*size_tmp+ta)+1] = imag(a);
	    }
	  }
	  cuda_surface.is_active[w][s] = (*it)->isActive();
	}
      }
    }

  }
}


void WaterSurface::setCudaM() {
  cuda_surface.initM(nb_wl);
  std::list<Obstacle*>::const_iterator ito;
  std::list<MovingEquivalentSource*>::const_iterator it;
  for (uint w = 0; w < nb_wl; ++w) {
    std::list<MovingEquivalentSource*> waves_wl;
    std::list<MovingEquivalentSource*> waves_wl_input = moving_waves[w];
    for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
      (*ito)->getSourcesM(waves_wl, w);
    }
    uint nb_sources = waves_wl.size();
    uint nb_sources_input = waves_wl_input.size();
    cuda_surface.allocMemM(w, nb_sources, nb_sources_input);
    uint s = 0;

    for (it = waves_wl_input.begin(); it != waves_wl_input.end(); ++it, ++s) {
      COMPLEX a = (*it)->getAmpli(0);
      cuda_surface.amplitudes_m[w][2*s*size_tmp] = real(a);
      cuda_surface.amplitudes_m[w][2*s*size_tmp+1] = imag(a);
      cuda_surface.indexes_m[w][s] = (*it)->getIndex();
      VEC2 pos = (*it)->getPos(0);
      cuda_surface.positions_m[w][2*s*size_tmp] = pos(0);
      cuda_surface.positions_m[w][2*s*size_tmp+1] = pos(1);
      cuda_surface.is_active_m[w][s] = (*it)->isActive();
    }
    for (it = waves_wl.begin(); it != waves_wl.end(); ++it, ++s) {
      COMPLEX a = (*it)->getAmpli(0);
      cuda_surface.amplitudes_m[w][2*s*size_tmp] = real(a);
      cuda_surface.amplitudes_m[w][2*s*size_tmp+1] = imag(a);
      cuda_surface.indexes_m[w][s] = (*it)->getIndex();
      VEC2 pos = (*it)->getPos(0);
      cuda_surface.positions_m[w][2*s*size_tmp] = pos(0);
      cuda_surface.positions_m[w][2*s*size_tmp+1] = pos(1);

      cuda_surface.is_active_m[w][s] = (*it)->isActive();
    }
  }

}

void WaterSurface::updateCudaM() {
  std::list<Obstacle*>::const_iterator ito;
  std::list<MovingEquivalentSource*>::const_iterator it;
  for (int w = nb_wl-1; w >= 0; --w) {
    if (wave_lenghts[w]/2.0 <= proj_grid.minSize()) {
      break;
    } else {
      std::list<MovingEquivalentSource*> waves_wl;
      std::list<MovingEquivalentSource*> waves_wl_input = moving_waves[w];
      uint nb_in_wave = waves_wl.size();
      for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
	(*ito)->getSourcesM(waves_wl, w);
      }
      uint s = 0;
      uint as = 0;
      for (it = waves_wl_input.begin(); it != waves_wl_input.end(); ++it, ++s) {
	COMPLEX a(0, 0);
	if (time/ampli_steps[w] > 0) {
	  int ta = (time/ampli_steps[w])%size_tmp+1;
	  a = (*it)->getAmpli(time/ampli_steps[w]);
	    VEC2 pos = (*it)->getPos(time/ampli_steps[w]-1);
	    if (real(a) != 0 || imag(a) != 0) {
	      cuda_surface.amplitudes_m[w][2*(s*size_tmp+ta)] = real(a);
	      cuda_surface.amplitudes_m[w][2*(s*size_tmp+ta)+1] = imag(a);
	    }
	    cuda_surface.positions_m[w][2*(s*size_tmp+ta)] = pos(0);
	    cuda_surface.positions_m[w][2*(s*size_tmp+ta)+1] =pos(1);	    
	  }
	  cuda_surface.is_active_m[w][s] = (*it)->isActive();
	}
	if (wave_lenghts[w]/2.0 > proj_grid.minSize()) {
	  for (it = waves_wl.begin(); it != waves_wl.end(); ++it, ++s) {
	    if (time/ampli_steps[w] > 0) {
	      int ta = (time/ampli_steps[w])%size_tmp+1;
	      COMPLEX a = (*it)->getAmpli(time/ampli_steps[w]-1);
	      VEC2 pos = (*it)->getPos(time/ampli_steps[w]);
	      if (real(a) != 0 || imag(a) != 0) {
		cuda_surface.amplitudes_m[w][2*(s*size_tmp+ta)] = real(a);
		cuda_surface.amplitudes_m[w][2*(s*size_tmp+ta)+1] = imag(a);
	      }
	      cuda_surface.positions_m[w][2*(s*size_tmp+ta)] = pos(0);
	      cuda_surface.positions_m[w][2*(s*size_tmp+ta)+1] =pos(1);
	    }
	    cuda_surface.is_active_m[w][s] = (*it)->isActive();
	  }
	}
      }

  }
}

#else  // not interactive
void WaterSurface::setCuda() {
  cuda_surface.init(nb_wl);
  
  std::list<Obstacle*>::const_iterator ito;
  std::list<Wave*>::const_iterator it;
  for (uint w = 0; w < nb_wl; ++w) {
    std::list<Wave*> waves_wl;
    std::list<Wave*> waves_wl_input = waves[w];
    for (ito = obstacles.begin(); ito != obstacles.end(); ++ito) {
      (*ito)->getSources(waves_wl, w);
    }
    uint nb_sources = waves_wl.size();
    uint nb_sources_input = waves_wl_input.size();
    cuda_surface.allocMem(w, nb_sources, nb_sources_input);
    uint s = 0;
    cuda_surface.wave_lenghts[w] = wave_lenghts[w];
    for (it = waves_wl_input.begin(); it != waves_wl_input.end(); ++it, ++s) {
      COMPLEX a = (*it)->getAmpli();
      cuda_surface.amplitudes[w][2*s] = real(a);
      cuda_surface.amplitudes[w][2*s+1] = imag(a);
      cuda_surface.indexes[w][s] = (*it)->getIndex();
      VEC2 pos = (*it)->getPos();
      cuda_surface.positions[w][2*s] = pos(0);
      cuda_surface.positions[w][2*s+1] = pos(1);
    }
    for (it = waves_wl.begin(); it != waves_wl.end(); ++it, ++s) {
      COMPLEX a = (*it)->getAmpli();
      cuda_surface.amplitudes[w][2*s] = real(a);
      cuda_surface.amplitudes[w][2*s+1] = imag(a);
      cuda_surface.indexes[w][s] = (*it)->getIndex();
      VEC2 pos = (*it)->getPos();
      cuda_surface.positions[w][2*s] = pos(0);
      cuda_surface.positions[w][2*s+1] = pos(1);
    }

  }
#ifdef PROJECTED_GRID
  proj_grid.setSizes();
  for (uint i = 0; i < n_rows_*n_cols_; ++i) {
    VEC3 p = proj_grid(i);
    cuda_surface.positions_grid[2*i] = p(0);
    cuda_surface.positions_grid[2*i+1] = p(1);
    cuda_surface.sizes[i] = proj_grid.size(i);
  }
#endif
}

void WaterSurface::updateCuda() {
  cuda_surface.setTimeHeight(time);
}

#endif

#ifdef PROJECTED_GRID

void WaterSurface::updatePosGridCuda() {
  proj_grid.setSizes();
#pragma omp parallel for
  for (uint i = 0; i < n_rows_*n_cols_; ++i) {
    VEC3 p = proj_grid(i);
    cuda_surface.positions_grid[2*i] = p(0);
    cuda_surface.positions_grid[2*i+1] = p(1);
    cuda_surface.sizes[i] = proj_grid.size(i);
  }

#ifdef INTERACTIVE_
  cuda_surface.setHeight(0);
#else
  cuda_surface.setHeight(waves.size());
#ifdef PLOT_RESULT
#pragma omp parallel for collapse(3)
  for (uint w = 0; w < nb_wl; ++w) {
    for (int i = 0; i < n_rows_ - 1; i++) {
      for (int j = 0; j < n_cols_ - 1; j++) {
	ampli_in_re[w](i, j) = cuda_surface.heights[w][2*(i*n_cols_ + j)];
	ampli_in_im[w](i, j) =  cuda_surface.heights[w][2*(i*n_cols_ + j)+1];
	ampli_scattered_re[w](i, j) = cuda_surface.heights[w][2*n_cols_*n_rows_ + 2*(i*n_cols_ + j)];
	ampli_scattered_im[w](i, j) = cuda_surface.heights[w][2*n_cols_*n_rows_ + 2*(i*n_cols_ + j)+1];
      }
    }
  }
#endif

#endif
  updateHeight();
}

void WaterSurface::setProjGrid(uint i, uint j, FLOAT x, FLOAT y) {
  proj_grid.setPosOnThePlane(i, j, x, y);
}

VEC3 WaterSurface::getPosProjGrid(uint i, uint j) const {
  VEC3 p = proj_grid(i, j);
  VEC2 wc = proj_grid.getPosWorld(i, j);
  Vector2i oc = world2gridObs(wc);
  VEC3 disp(0, 0, 0);
  if (use_gersner) {
#ifdef INTERACTIVE_
    disp(0) = -cuda_surface.displacement[2*(i*n_cols_ + j)];
    disp(1) = -cuda_surface.displacement[2*(i*n_cols_ + j)+1];
#else
    disp(0) = -cuda_surface.time_displacement[2*(i*n_cols_ + j)];
    disp(1) = -cuda_surface.time_displacement[2*(i*n_cols_ + j)+1];
#endif
  }
  disp(2) = 0;
  return p + disp;
}

VEC3 WaterSurface::getPosProjGrid(uint i) const {
  VEC3 p = proj_grid(i);
  VEC2 wc = proj_grid.getPosWorld(i/n_cols_, i - (int)(i/n_cols_)*n_cols_);
  VEC3 disp(0, 0, 0);
  if (use_gersner) {
#ifdef INTERACTIVE_
    disp(0) = -cuda_surface.displacement[2*i];
    disp(1) = -cuda_surface.displacement[2*i];
#else
    disp(0) = -cuda_surface.time_displacement[2*i];
    disp(1) = -cuda_surface.time_displacement[2*i+1];
#endif
  }
  disp(2) = 0;
  Vector2i oc = world2gridObs(wc);
  return p + disp;
}

void WaterSurface::setTargetLookAt(VEC2 target) {
  target_lookat = viewer2world(target);
}

#endif

#endif

VEC3 WaterSurface::getPosGrid(uint i, uint j) const {
  VEC2 pv = grid2viewer(i, j);
#ifdef USE_CUDA
  VEC2 wc = grid2world(i, j);
  VEC2 disp(0, 0);
  if (use_gersner) {
#ifdef INTERACTIVE_
    disp(0) = -cuda_surface.displacement[2*(i*n_cols_ + j)];
    disp(1) = -cuda_surface.displacement[2*(i*n_cols_ + j)+1];
#else
    disp(0) = -cuda_surface.time_displacement[2*(i*n_cols_ + j)];
    disp(1) = -cuda_surface.time_displacement[2*(i*n_cols_ + j)+1];
#endif
  }
  Vector2i oc = world2gridObs(wc);
  VEC3 p(pv(0) + disp(0), pv(1) + disp(1), u(i, j));
  if (b(oc(0), oc(1)) != 0 && p(2) > b(oc(0), oc(1))) {
    p(2) = b(oc(0), oc(1))- 0.01;
    disp = VEC2(0, 0);
  }
  return p;
#else
  return VEC3(pv(0), pv(1), u(i, j));
#endif
}

VEC3 WaterSurface::getPosGrid(uint i) const {
  VEC2 pv = grid2viewer(i/n_cols_, i - (int)(i/n_cols_)*n_cols_);
#ifdef USE_CUDA
  VEC2 wc = grid2world(i/n_cols_, i - (int)(i/n_cols_)*n_cols_);
  VEC2 disp(0, 0);
  if (use_gersner) {
#ifdef INTERACTIVE_
    disp(0) = -cuda_surface.displacement[2*i];
    disp(1) = -cuda_surface.displacement[2*i+1];
#else
    disp(0) = -cuda_surface.time_displacement[2*i];
    disp(1) = -cuda_surface.time_displacement[2*i+1];
#endif
  }
  Vector2i oc = world2gridObs(wc);
  VEC3 p(pv(0) + disp(0), pv(1) + disp(1), u(i/n_cols_, i - (int)(i/n_cols_)*n_cols_));
  if (b(oc(0), oc(1)) != 0 && p(2) > b(oc(0), oc(1))) {
    p(2) = b(oc(0), oc(1))- 0.01;
    disp = VEC2(0, 0);
  }
  return p;
#else
  return VEC3(pv(0), pv(1), u(i/n_cols_, i - (int)(i/n_cols_)*n_cols_));
#endif
}



FLOAT WaterSurface::minWL() const {
  return wave_lenghts[0];
}
FLOAT WaterSurface::maxWL() const {
  return wave_lenghts[nb_wl-1];
}


void WaterSurface::addMovingSource(VEC2 pos0, VEC2 speed, FLOAT t0, FLOAT t1) {

  FLOAT sp = speed.norm();
  FLOAT maink = gravity_/(4.0*pow(sp*cos(35.3*M_PI/180.0f), 2));
  FLOAT maxk = gravity_/(4.0*pow(sp, 2));
  FLOAT mainwl = 2.0*M_PI/maink;
  FLOAT maxwl = 2.0*M_PI/maxk;

  int start = t0/dt_;
  int stop = t1/dt_;
  int t_start = (floor(start/ampli_steps[0]) + 1)*ampli_steps[0];
  int t_stop = (floor(stop/ampli_steps[0]) + 1)*ampli_steps[0];
  for (int w = nb_wl - 1; w >= 0; --w) {
    FLOAT wl = wave_lenghts[w];
    if (wl <= 1.01*maxwl) {
      t_start = (floor(start/ampli_steps[w]) + 1)*ampli_steps[w];
      t_stop = (floor(stop/ampli_steps[w]) + 1)*ampli_steps[w];
      break;
    }
  }
  FLOAT s = 0;
  // for (uint w = 0; w < nb_wl; ++w) {
  //   FLOAT wl = wave_lenghts[w];
  //   if (wl <= 1.01*size && wl > 0.05*size) {
  // 	s += sqrt(wl);
  //   }
  // }

  for (int w = 0; w < nb_wl; ++w) {
    FLOAT wl = wave_lenghts[w];
    if (wl <= 1.01*maxwl && wl > 0.5*mainwl) {
      FLOAT period = wl/velocity(2*M_PI/wl);
      int d_period = period/(dt_*ampli_steps[w]);
      int nb_steps = (t_stop - t_start)/ampli_steps[w];
      VEC2 step2d = dt_*ampli_steps[w]*speed;
      VEC2 pos_cur = pos0;
      int t_cur = t_start/ampli_steps[w];
      for (uint i = 0; i < nb_steps; ++i) {
	Wave *ms = new EquivalentSource(wl, ampli_steps[w]);
	((EquivalentSource*)ms)->setPos(pos_cur);
	((EquivalentSource*)ms)->setAmplitude(sqrt(mainwl)*height_ampli_, t_cur, t_cur+3);
	waves[w].push_back(ms);
	pos_cur += step2d;
	t_cur += 1;
      }
    }
  }
}
