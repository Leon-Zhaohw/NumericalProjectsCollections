/* 
 * File: Obstacle.cpp
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

#include "Obstacle.hpp"

#include <iostream>
#include "WaterSurface.hpp"
#include "settings.hpp"
#include "error.hpp"

using namespace Eigen;
using namespace definitions;
using namespace settings;

void Obstacle::setAmpliSources(const std::list<Wave*> &waves, uint index, int t) {
  if (boundaries_l[index].size() != 0 && sources_l[index].size() != 0) {
    VectorXcf c(sources_l[index].size());
    VectorXcf p_in( boundaries_l[index].size());
    std::list<Wave*>::const_iterator it;
    bool to_solve = false;
    for (uint i = 0; i < boundaries_l[index].size(); ++i) {
      VEC2C n = normals_l[index][i];
      p_in(i) = 0;
      VEC2 pb = boundaries_l[index][i]->getPos();
      for (it = waves.begin(); it != waves.end(); ++it) {
	if ((*it)->isActive()) {
	  if (neumann) {
	    p_in(i) -= n.dot((*it)->gradHeightc(pb, t*dt_*ampli_steps[index]));
	  } else {
	    p_in(i) -= (*it)->heightc(pb,t*dt_*ampli_steps[index]);
	  }
	}
      }

      for (uint j = 0; j < sources_l[index].size(); ++j) {
	if (sources_l[index][j]->isActive()) {
	  FLOAT d = sources_l[index][j]->getDistance(pb(0), pb(1));
	  FLOAT ret = d/velocity(2*M_PI/wave_lenghts[index]);
	  FLOAT q = interpolation(-ret, 0, dt_*ampli_steps[index]);
	  if (q == 0) {
	    if (neumann) {
	      p_in(i) -= n.dot(sources_l[index][j]->gradHeightc(pb,t*dt_*ampli_steps[index]));
	    } else {
	      p_in(i) -= 0.8f*sources_l[index][j]->heightc(pb, t*dt_*ampli_steps[index]);
	    }
	  }
	  if (q != 0 && !approx_inter) {
	    if (neumann) {
	      p_in(i) -= n.dot(sources_l[index][j]->gradHeightc(pb,t*dt_*ampli_steps[index]));
	    } else {
	      p_in(i) -= sources_l[index][j]->heightc(pb, t*dt_*ampli_steps[index]);
	    }
	  }
	}
      }
    }
    for (uint i = 0; i < boundaries_l[index].size() && !to_solve; ++i) {
      to_solve = to_solve || (norm(p_in(i)) > 1e-8*wave_lenghts[index]);
    }
    if (to_solve) {
      c = svd_l[index].solve(p_in);
      for (uint i = 0; i < sources_l[index].size(); ++i) {
	if (norm(c(i)) < 1e-8) {
	  sources_l[index][i]->setAmplitude(0, t);
	} else {
	  c(i) = (1.0f - damping_source_)*sources_l[index][i]->getAmpli(t-1) + damping_source_*c(i);
	  IS_DEF(real(c(i)));
	  IS_DEF(imag(c(i)));
	  sources_l[index][i]->setAmplitude(c(i), t);
	}
      }
    }
  }
}


void Obstacle::setTransfeMatrix(uint index, int t) {
  INFO("Set Transfer matrix:  "<<index<<"  (wavelength "<<wave_lenghts[index]<<", "<<sources_l[index].size()<<" sources)");
  BDCSVD<MatrixXcf> svd;
  FLOAT k = 2*M_PI/wave_lenghts[index];
  std::vector<EquivalentSource*> &sources = sources_l[index];
  std::vector<InputPoint*> &boundaries = boundaries_l[index];
  if (boundaries.size() != 0) {
    std::vector<VEC2C> &normals = normals_l[index];
    MatrixXcf T = MatrixXcf(boundaries.size(), sources.size());
    if (t <= 1) {
      uint nb_min = sources.size();
      for (uint j = 0; j < sources.size(); ++j) {
	uint nb = 0;
	for (uint i = 0; i < boundaries.size(); ++i) {
	  VEC2C n = normals[i];
	  VEC2 pb = boundaries[i]->getPos();
	  FLOAT d = sources[j]->getDistance(pb);
	  FLOAT ret = d/velocity(k);
	  FLOAT q = interpolation(dt_*ampli_steps[index] -ret, 1, dt_*ampli_steps[index]);
	  if (q != 0) {
	    if (approx_inter) {
	      q = 1;
	    }
	    if (neumann) {
	      T(i, j) = q*n.dot(sources[j]->gradHeightc(pb));
	    } else {
	      T(i, j) = q*sources[j]->heightc(pb);
	    }
	    ++nb;
	  } else {
	    T(i, j) = 0;
	  }
	}
	if (nb < nb_min) {
	  nb_min = nb;
	}
      }
      INFO("--- nb xb influencing source "<<nb_min);
      INFO("SVD..."<<index<<"  (wavelength "<<wave_lenghts[index]<<")");
      svd = BDCSVD<MatrixXcf>(T,ComputeThinU | ComputeThinV);
      svd_l[index] = svd;
      INFO("SVD...done:  "<<index<<"  (wavelength "<<wave_lenghts[index]<<")");
    }
  }
}

// frequency domain solve
void Obstacle::setAmpliSources(const std::list<Wave*> &waves, uint index) {
  if (boundaries_l[index].size() != 0) {
    BDCSVD<MatrixXcf> svd = svd_l[index];
    std::vector<EquivalentSource*> &sources = sources_l[index];
    std::vector<InputPoint*> &boundaries = boundaries_l[index];
    if (boundaries.size() != 0) {
      VectorXcf c(sources.size());
      VectorXcf p_in(boundaries.size());

      // fill p_in
      std::list<Wave*>::const_iterator it;
      if (neumann) {
	std::vector<VEC2C> &normals = normals_l[index];
	for (uint i = 0; i < boundaries.size(); ++i) {
	  p_in(i) = 0;
	  VEC2 pos = boundaries[i]->getPos();
	  VEC2C n = normals[i];
	  for (it = waves.begin(); it != waves.end(); ++it) {
	    p_in(i) -= n.dot((*it)->gradHeightc(boundaries[i]->getPos(), 0));
	  }
	}
      } else {
	for (uint i = 0; i < boundaries.size(); ++i) {
	  p_in(i) = 0;
	  VEC2 pos = boundaries[i]->getPos();
	  for (it = waves.begin(); it != waves.end(); ++it) {
	    p_in(i) -= (*it)->heightc(pos, 0);
	  }
	}
      }
  
      // solve
      c = svd.solve(p_in);

      /* error of the least square */
      // VectorXcf err = T*c - p_in;
      // FLOAT max = 0;
      // FLOAT mean = 0;
      // FLOAT mean_in = 0;
      // FLOAT e_max = 0;
      // FLOAT in_max = 0;
      // for (uint i = 0; i < boundaries.size(); ++i) {
      //   float e = sqrt(norm(err(i)));
      //   float in = sqrt(norm(p_in(i)));
      //   mean += e;
      //   mean_in += in;
      //   if (e > max) {
      //     max = e;
      //     e_max = e;
      //     in_max = in;
      //   }
      // }
      // mean /= (FLOAT)boundaries.size();
      // mean_in /= (FLOAT)boundaries.size();
      // INFO("err\n"<<err);
      // INFO("err max "<<max<<" "<<e_max<<" "<<in_max);
      // INFO("err mean "<<mean<<" "<<mean_in);

      // update ampli sources
      for (uint i = 0; i < sources.size(); ++i) {
	sources[i]->setAmplitude(c(i));
      }
    }
  }
}



void Obstacle::setTransfeMatrix(uint index) {
  INFO("Set Transfer matrix:  "<<index<<"  (wavelength "<<wave_lenghts[index]<<", "<<sources_l[index].size()<<" sources)");
  if (boundaries_l[index].size() != 0) {
    std::vector<EquivalentSource*> &sources = sources_l[index];
    std::vector<InputPoint*> &boundaries = boundaries_l[index];
  
    MatrixXcf T = MatrixXcf(boundaries.size(), sources.size());
    if (neumann) {
      std::vector<VEC2C> &normals = normals_l[index];
      for (uint i = 0; i < boundaries.size(); ++i) {
	VEC2C n = normals[i];
	for (uint j = 0; j < sources.size(); ++j) {
	  T(i, j) = n.dot(sources[j]->gradHeightc(boundaries[i]->getPos()));
	}
      }
    } else {
      for (uint i = 0; i < boundaries.size(); ++i) {
	for (uint j = 0; j < sources.size(); ++j) {
	  T(i, j) = sources[j]->heightc(boundaries[i]->getPos());
	}
      }    
    }
    INFO("SVD..."<<index<<"  (wavelength "<<wave_lenghts[index]<<")");
    svd_l[index] = BDCSVD<MatrixXcf>(T,ComputeThinU | ComputeThinV);;
    INFO("SVD...done:  "<<index<<"  (wavelength "<<wave_lenghts[index]<<")");
  }
}


Obstacle::Obstacle() {
  pos = VEC2(0, 0);
  nb_wl = 0;
  movable_ = false;
}

Obstacle::Obstacle(VEC2 p) {
  pos = p;
  nb_wl = 0;
  movable_ = false;
}

void Obstacle::reset(const std::vector<std::list<Wave*> > &waves_l, const std::vector<FLOAT> &wl) {
  nb_wl = wl.size();

  wave_lenghts = wl;
  boundaries_l = std::vector<std::vector<InputPoint*> >(nb_wl);
  normals_l = std::vector<std::vector<VEC2C> >(nb_wl);
  sources_l = std::vector<std::vector<EquivalentSource*> >(nb_wl);
 svd_l =  std::vector<Eigen::BDCSVD<MatrixXcf> >(nb_wl);

 uint nm = 0;
 uint j = 0;
 uint nb_threads = NTHREADS_;
 while (nm < nb_wl) {
   nm = std::min((j+1)*nb_threads, nb_wl);
#pragma omp parallel for
   for (uint i = j*nb_threads; i < nm; ++i) {
     setBoundaries(i);
     setEquivalentSources(i);
     if (!sources_l[i].empty() && !boundaries_l[i].empty()) {
#ifdef INTERACTIVE_
       setTransfeMatrix(i, 0);
       setAmpliSources(waves_l[i], i, 0);
#else
       setTransfeMatrix(i);
       setAmpliSources(waves_l[i], i);
#endif
    }
  }
  ++j;
 }

 
// #pragma omp parallel for
//   for (uint i = 0; i < nb_wl; ++i) {
//     setBoundaries(i);
//     setEquivalentSources(i);
//     if (!sources_l[i].empty() && !boundaries_l[i].empty()) {
// #ifdef INTERACTIVE_
//       setTransfeMatrix(i, 0);
//       setAmpliSources(waves_l[i], i, 0);
// #else
//       setTransfeMatrix(i);
//       setAmpliSources(waves_l[i], i);
// #endif
//     }
//   }
}


Obstacle::~Obstacle() {
  for (uint w = 0; w < nb_wl; ++w) {
    for (uint i = 0; i < boundaries_l[w].size(); ++i) {
      delete boundaries_l[w][i];
    }
    for (uint i = 0; i < sources_l[w].size(); ++i) {
      delete sources_l[w][i];
    }
  }
}

void Obstacle::setPos(VEC2 p) {
  pos = p;
}

void Obstacle::setPos(FLOAT x, FLOAT y) {
  pos = VEC2(x, y);
}

VEC2 Obstacle::getPos() {
  return pos;
}

void Obstacle::update(const WaterSurface *surface) {
  for (uint w = 0; w < nb_wl; ++w) {
    for (uint i = 0; i < boundaries_l[w].size(); ++i) {
      VEC2 pos = boundaries_l[w][i]->getPos();
      FLOAT height = surface->height(pos[0], pos[1]);
      boundaries_l[w][i]->update(height);
    }
  }
}

void Obstacle::getBoundariesPos(std::list<VEC2> &boundaryPos) const {
  for (uint w = 0; w < nb_wl; ++w) {
    for (uint i = 0; i < boundaries_l[w].size(); ++i) {
      boundaryPos.push_back(world2viewer(boundaries_l[w][i]->getPos()));
    }
    for (uint i = 0; i < sources_l[w].size(); ++i) {
      boundaryPos.push_back(world2viewer(sources_l[w][i]->getPos()));
    }
  }


}



FLOAT Obstacle::height(FLOAT x, FLOAT y, FLOAT time) const{
  FLOAT a = 0;
  std::vector<std::vector<EquivalentSource*> >::const_iterator it;
  for (it = sources_l.begin(); it != sources_l.end(); ++it) {
    std::vector<EquivalentSource*> sources = *it;
    for (uint i = 0; i < sources.size(); ++i) {
      a += sources[i]->height(x, y, time);
    }
  }
  return a;
}

COMPLEX Obstacle::heightc(FLOAT x, FLOAT y, FLOAT time) const{
  COMPLEX a = 0;
  std::vector<std::vector<EquivalentSource*> >::const_iterator it;
  for (it = sources_l.begin(); it != sources_l.end(); ++it) {
    std::vector<EquivalentSource*> sources = *it;
    for (uint i = 0; i < sources.size(); ++i) {
      a += sources[i]->heightc(x, y, time);
    }
  }
  return a;
}

COMPLEX Obstacle::heightc_wl(FLOAT x, FLOAT y, uint w) const{
  COMPLEX a = 0;
  std::vector<EquivalentSource*> sources = sources_l[w];
  for (uint i = 0; i < sources.size(); ++i) {
    a += sources[i]->heightc(x, y, 0);
  }
  return a;
}


COMPLEX Obstacle::heightc_wl(FLOAT x, FLOAT y, uint w, FLOAT time) const{
  COMPLEX a = 0;
  
  std::vector<EquivalentSource*> sources = sources_l[w];
  for (uint i = 0; i < sources.size(); ++i) {
    a += sources[i]->heightc(x, y, time);
  }
  return a;
}


void Obstacle::getSources(std::list<Wave*> &waves, uint w) {
  if (!movable_) {
    for (uint i = 0; i < sources_l[w].size(); ++i) {
      waves.push_back(sources_l[w][i]);
    }
  } 
}

void Obstacle::getSourcesM(std::list<MovingEquivalentSource*> &waves, uint w) {
  if (movable_) {
    for (uint i = 0; i < sources_l[w].size(); ++i) {
      waves.push_back((MovingEquivalentSource*)sources_l[w][i]);
    }
  }
}


std::vector<InputPoint*> &Obstacle::getBoundaries(uint w) {
  return boundaries_l[w];
}

std::vector<VEC2C> &Obstacle::getNormales(uint w) {
  return normals_l[w];
}

std::vector<EquivalentSource*> &Obstacle::getSources(uint w) {
  return sources_l[w];
}
		    

int Obstacle::getBoundariesSize(uint w) const {
  return boundaries_l[w].size();
}

int Obstacle::getSourcesSize(uint w) const {
  return sources_l[w].size();
}


std::ofstream &Obstacle::exportMitsuba(std::ofstream &file) const {
  file<<"<shape type=\"obj\">\n";
  file<<"<string name=\"filename\" value=\""<<str_obstacle_grid<<"\"/>\n";
  file<<"<bsdf type=\"diffuse\">\n";
  file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  file<<"</bsdf>\n";
  file<<"</shape>\n";
  return file;
}

void Obstacle::move(VEC2 dir, FLOAT speed, int time) {
  VEC2 step2d = speed*dt_*dir;
  pos += step2d;
  for (uint w = 0; w < nb_wl; ++w) {
    for (uint i = 0; i < sources_l[w].size(); ++i) {
      if (time%(4*ampli_steps[w])==0) {
	((MovingEquivalentSource*)sources_l[w][i])->move(4*ampli_steps[w]*step2d,time/ampli_steps[w]);
      } else {
	((MovingEquivalentSource*)sources_l[w][i])->move(VEC2(0, 0),time/ampli_steps[w]);
      }
    }      
    for (uint i = 0; i < boundaries_l[w].size(); ++i) {
      boundaries_l[w][i]->move(step2d);
    }
  }

     
}


void Obstacle::setMovable(bool mv) {
  movable_ = mv;
}
bool Obstacle::isMovable() const {
  return movable_;
}
