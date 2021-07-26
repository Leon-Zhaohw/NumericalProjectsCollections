/* 
 * File: CudaSurface.cu
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

#include "CudaWaterSurface.hpp"
#include "settings.hpp"
#include "definitions.hpp"
#include <iostream>
#include "error.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include <cuda_fp16.h>
#include "ui_parameters.hpp"

using namespace settings;
using namespace definitions;
    
__global__
void addHeight(FLOAT *heights, FLOAT *displacement, FLOAT *amplitudes, FLOAT *indexes,
	       FLOAT *positions, bool *is_active, FLOAT *positions_grid, FLOAT *sizes,
	       uint nb_sources, uint nb_sources_input,
	       FLOAT k, FLOAT omega, FLOAT vel,
	       FLOAT time, FLOAT dt, uint nb_rows, uint nb_cols,
	       int size_ampli, int ampli_step, FLOAT damping, FLOAT scale,
	       bool show_input, bool show_scattered,
	       FLOAT *hankel_r_tab, FLOAT *hankel_i_tab ) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (uint i = index; i < nb_rows*nb_cols; i += stride) {
#ifdef PROJECTED_GRID
    FLOAT wl = 2*M_PI/k;
    if (sizes[i] >= wl/2.0) {
      return;
    }
      // position of the projected grid are recorded with the viewer coordinates
      // we need to set them in world coordinates
      FLOAT x = (1+positions_grid[2*i])*scale/2.0;
      FLOAT y = (1+positions_grid[2*i+1])*scale/2.0;
  
      // FLOAT dispxr = 0, dispyr = 0;
      // FLOAT dispxi = 0, dispyi = 0;
      FLOAT mod = 1;
      if (sizes[i] > wl/4.0) {
	mod = 1 - (sizes[i] - wl/4.0)/(wl/4.0);
      }
#else
    FLOAT x = scale/(FLOAT)nb_rows*(i/nb_cols);
    FLOAT y = scale/(FLOAT)nb_cols*(i%nb_cols);
    FLOAT mod = 1;
#endif
    FLOAT hr = 0, hi = 0;
      uint start = 0, end = nb_sources + nb_sources_input;
      if (!show_input) {
	start = nb_sources_input;
      }
      if (!show_scattered) {
	end = nb_sources_input;
      }

      for (uint s = start; s < end; ++s) {
      while (!is_active[s] && s < end) {
	++s;
      }
      if (s >= end) {
	break;
      }
      if (indexes[s] == 1) {
	if (time > x/vel) {
	  FLOAT dx = positions[2*s];
	  FLOAT dy = positions[2*s + 1];
	  FLOAT kx = k*(x*dx + y*dy);
	  FLOAT ar = amplitudes[2*s*size_ampli];
	  FLOAT ai = amplitudes[2*s*size_ampli + 1];
	  hr += ar*cos(kx) - ai*sin(kx);
	  hi += ai*cos(kx) + ar*sin(kx);
	  // dispxr -= hi*dx;
	  // dispxi += hr*dx;
	  // dispyr -= hi*dy;
	  // dispyi += hr*dy;
	}
      } else {
	FLOAT rx = x - positions[2*s];
	FLOAT ry = y - positions[2*s+1];
	  
	FLOAT r = sqrt((float)(rx*rx + ry*ry));
	// FLOAT kx = rx/r;
	// FLOAT ky = ry/r;
	FLOAT damp = exp(-damping*k*k*r);
	if (damp > 0.02) { 
	  FLOAT ret = r/vel;
	  int l = floor((time-ret)/((FLOAT)ampli_step*dt));
	  if (time-ret < 0) {
	    --l;
	  }

	  FLOAT ar = 0, ai = 0;
	  // if (l > size_ampli-1) {
	  //   // ar = amplitudes[2*(s*size_ampli + size_ampli-1)];
	  //   // ai = amplitudes[2*(s*size_ampli + size_ampli-1)+1];
	  // } else
	    if (l < 0) {
	    ar = amplitudes[2*s*size_ampli];
	    ai = amplitudes[2*s*size_ampli+1];
	  
	  } else {
	    FLOAT t = time-ret;
	    FLOAT w = 0;
	    FLOAT tl = l*ampli_step*dt;
	    FLOAT tl_prev = (l-1)*ampli_step*dt, tl_next = (l+1)*ampli_step*dt;
	    if (t < tl_prev || t > tl_next) {
	      w = 0;
	    } else if (t < tl) {
	      w = (t - tl_prev)/(dt*(FLOAT)ampli_step);
	    } else {
	      w = (tl_next - t)/(dt*(FLOAT)ampli_step);
	    }
	    //	    w = 1;
	    ar = w * amplitudes[2*(s*size_ampli+l%size_ampli)] +
	      (1-w)*amplitudes[2*(s*size_ampli + (l+1)%size_ampli)];
	    ai = w * amplitudes[2*(s*size_ampli+l%size_ampli)+1] +
	      (1-w)*amplitudes[2*(s*size_ampli + (l+1)%size_ampli)+1];
	  }
	  if (r > 0.001) {
	    // FLOAT han_r = sqrt((float)(2.0/((FLOAT)M_PI*k*r)))*cos(k*r - (FLOAT)M_PI/4.0);
	    // FLOAT han_i = sqrt((float)(2.0/((FLOAT)M_PI*k*r)))*sin(k*r - (FLOAT)M_PI/4.0);

	    uint ind = floor(k*r/0.025);
	     FLOAT coef = k*r/0.025 - ind;
	    if (ind >= 99999) {
	      ind = 0;
	      coef = 0;
	    }
#ifdef PROJECTED_GRID
	    if (s < nb_sources_input && ind < 10) {
	      ind = 10;
	      coef = 0;
	    }
#endif
	    FLOAT han_r = //hankel_r_tab[ind];
		    (1-coef)*hankel_r_tab[ind] + coef*hankel_r_tab[ind+1];
	    FLOAT han_i = //hankel_i_tab[ind];
		    (1-coef)*hankel_i_tab[ind] + coef*hankel_i_tab[ind+1];
	    FLOAT tmpr = han_r*ar - han_i*ai;
	    FLOAT tmpi = han_r*ai + han_i*ar;
	    tmpr *= damp;
	    tmpi *= damp;
	    hr += tmpr;
	    hi += tmpi;
	    // dispxr += tmpi*kx;
	    // dispxi -= tmpr*kx;
	    // dispyr += tmpi*ky;
	    // dispyi -= tmpr*ky;
	  }
	}
	//	}
      }
	
    }
    hr *= mod;
    hi *= mod;
    heights[i] += hr*cos(-omega*time) - hi*sin(-omega*time);
#ifdef PLOT_RESULT
    heights[nb_cols*nb_rows + i] = sqrt(hr*hr +  hi*hi);
#endif
    // displacement[2*i] += dispxr*cos(-omega*time) - dispxi*sin(-omega*time);
    // displacement[2*i+1] += dispyr*cos(-omega*time) - dispyi*sin(-omega*time);
    
  }
}

__global__

void addHeightM(FLOAT *heights, FLOAT *displacement, FLOAT *amplitudes, FLOAT *indexes,
		FLOAT *positions, 
		bool *is_active, FLOAT *positions_grid, FLOAT *sizes,
		uint nb_sources, uint nb_sources_input,
		FLOAT k, FLOAT omega, FLOAT vel,
		FLOAT time, FLOAT dt, uint nb_rows, uint nb_cols,
		int size_ampli, int ampli_step, FLOAT damping, FLOAT scale,
		bool show_input, bool show_scattered,
		FLOAT *hankel_r_tab, FLOAT *hankel_i_tab ) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (uint i = index; i < nb_rows*nb_cols; i += stride) {
    FLOAT wl = 2*M_PI/k;
#ifdef PROJECTED_GRID
    if (sizes[i] >= wl/2.0) {
      return;
    }
    FLOAT x = (1+positions_grid[2*i])*scale/2.0;
    FLOAT y = (1+positions_grid[2*i+1])*scale/2.0;
    // FLOAT dispxr = 0, dispyr = 0;
    // FLOAT dispxi = 0, dispyi = 0;
      FLOAT mod = 1;
      if (sizes[i] > wl/4.0) {
	mod = 1 - (sizes[i] - wl/4.0)/(wl/4.0);
      }
#else
      FLOAT x = scale/(FLOAT)nb_rows*(i/nb_cols);
      FLOAT y = scale/(FLOAT)nb_cols*(i%nb_cols);
      FLOAT mod = 1;
#endif
      FLOAT hr = 0, hi = 0;
      uint start = 0, end =  nb_sources + nb_sources_input;
      if (!show_input) {
	start = nb_sources_input;
      }
      if (!show_scattered) {
	end = nb_sources_input;
      }
      for (uint s = start; s < end; ++s) {
	
      while (!is_active[s] && s < end) {
	++s;
      }
      if (s >= end) {
	break;
      }
      int t_cur = floor((time)/((FLOAT)ampli_step*dt));
      FLOAT rmax = log(0.02)/(-damping*k*k);
      FLOAT time_lim = rmax/vel;
      int lim = 2*time_lim / (ampli_step*dt) + 1;
      int tb = t_cur-lim;
      if (tb < 0) {
	tb = 0;
      }
      for (int past_t = tb; past_t <= t_cur; ++past_t) {
	FLOAT rx = x - positions[2*(s*size_ampli + past_t)];
	FLOAT ry = y - positions[2*(s*size_ampli + past_t)+1];
	FLOAT r = sqrt((float)(rx*rx + ry*ry));
	FLOAT damp = exp(-damping*k*k*r);
	if (damp > 0.02) { 
	  FLOAT ret = r/vel;
	  int l = floor((time-ret)/((FLOAT)ampli_step*dt));
	  if (time-ret > 0) {
 	  
	    FLOAT ar = 0, ai = 0;
	     if (l >= past_t - 1 || l <= past_t +1) {
	       int pt = past_t;
	       //for (int pt = past_t - 1; pt <= past_t + 1; ++pt) {
	      FLOAT t = time-ret;
	      FLOAT w = 0;
	      FLOAT tl = pt*ampli_step*dt;
	      FLOAT tl_prev = (pt-1)*ampli_step*dt, tl_next = (pt+1)*ampli_step*dt;
	      if (t < tl_prev || t > tl_next) {
		w = 0;
	      } else if (t < tl) {
		w = (t - tl_prev)/(dt*(FLOAT)ampli_step);
	      } else {
		w = (tl_next - t)/(dt*(FLOAT)ampli_step);
	      }
	      ar = w * amplitudes[2*(s*size_ampli + pt)];
	      ai = w * amplitudes[2*(s*size_ampli + pt) + 1];
	    
	      uint ind = floor(k*r/0.025);
	      FLOAT coef = k*r/0.025 - ind;
	      if (ind >= 99999) {
		ind = 0;
		coef = 0;
	      }
#ifdef PROJECTED_GRID
	      if (s < nb_sources_input && ind < 10) {
		ind = 10;
		coef = 0;
	      }
#endif
	      FLOAT han_r = (1-coef)*hankel_r_tab[ind] + coef*hankel_r_tab[ind+1];
	      //hankel_r_tab[ind];
	      FLOAT han_i = (1-coef)*hankel_i_tab[ind] + coef*hankel_i_tab[ind+1];
	      //hankel_i_tab[ind];
	      FLOAT tmpr = han_r*ar - han_i*ai;
	      FLOAT tmpi = han_r*ai + han_i*ar;
	      tmpr *= damp;
	      tmpi *= damp;
	      hr += tmpr;
	      hi += tmpi;
   
	    }
	  }
	}
      }
    }
    hr *= mod;
    hi *= mod;
    heights[i] += hr*cos(-omega*time) - hi*sin(-omega*time);
#ifdef PLOT_RESULT
    heights[nb_cols*nb_rows + i] = sqrt(hr*hr +  hi*hi);
#endif
    // displacement[2*i] += dispxr*cos(-omega*time) - dispxi*sin(-omega*time);
    // displacement[2*i+1] += dispyr*cos(-omega*time) - dispyi*sin(-omega*time);
  }
}
      

__global__
void addHeight0(FLOAT *heights, FLOAT *displacement, FLOAT *amplitudes,
		      FLOAT *indexes, FLOAT *positions,
		      FLOAT *positions_grid, FLOAT *sizes,
		      uint nb_sources, uint nb_sources_input, FLOAT k, FLOAT omega, FLOAT vel,
		      uint nb_rows, uint nb_cols, FLOAT damping,
		      FLOAT scale, bool show_input, bool show_scattered,
		      FLOAT *hankel_r_tab, FLOAT *hankel_i_tab) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (uint i = index; i < nb_rows*nb_cols; i += stride) {
    FLOAT wl = 2.0*M_PI/k;
    FLOAT mod = 1;
#ifdef PROJECTED_GRID
    if (sizes[i] >= wl/2.0) {
      return;
    }
    if (sizes[i] > wl/4.0) {
      mod = 1 - (sizes[i] - wl/4.0)/(wl/4.0);
    }
    FLOAT x = (1+positions_grid[2*i])*scale/2.0;
    FLOAT y = (1+positions_grid[2*i+1])*scale/2.0;
#else
    FLOAT x = scale/(FLOAT)nb_rows*(i/nb_cols);
    FLOAT y = scale/(FLOAT)nb_cols*(i%nb_cols);
 #endif
     uint start = 0, end = nb_sources + nb_sources_input;
      if (!show_input) {
	start = nb_sources_input;
      }
      if (!show_scattered) {
	end = nb_sources + nb_sources_input;
      }
      for (uint s = start; s < end; ++s) {
      FLOAT hr = 0, hi = 0;
      FLOAT dispxr = 0, dispyr = 0;
      FLOAT dispxi = 0, dispyi = 0;
      if (indexes[s] == 1) {
	FLOAT dx = positions[2*s];
	FLOAT dy = positions[2*s + 1];
	FLOAT kx = k*(x*dx + y*dy);
	FLOAT ar = amplitudes[2*s];
	FLOAT ai = amplitudes[2*s + 1];
	hr = ar*cos(kx) - ai*sin(kx);
	hi = ai*cos(kx) + ar*sin(kx);
	hr *= mod;
	hi *= mod;
	dispxr = hi*dx;
	dispxi = -hr*dx;
	dispyr = hi*dy;
	dispyi = -hr*dy;
      } else {
	FLOAT rx = x - positions[2*s];
	FLOAT ry = y - positions[2*s+1];
      
	FLOAT r = sqrt((float)(rx*rx + ry*ry));
	FLOAT kx = rx/r;
	FLOAT ky = ry/r;
	FLOAT ar = amplitudes[2*s], ai = amplitudes[2*s+1];
     
	FLOAT damp = exp(-damping*k*k*r);
	//	  if (r > 0.0001) {
	// FLOAT han_r = sqrt((float)(2.0/((FLOAT)M_PI*k*r)))*cos(k*r - (FLOAT)M_PI/4.0);
	// FLOAT han_i = sqrt((float)(2.0/((FLOAT)M_PI*k*r)))*sin(k*r - (FLOAT)M_PI/4.0);

	uint ind = floor(k*r/0.025);
	FLOAT coef = k*r/0.025 - ind;
	if (ind >= 999999) {
	  ind = 0;
	  coef = 0;
	}
	FLOAT han_r = (1-coef)*hankel_r_tab[ind] + coef*hankel_r_tab[ind+1];
	FLOAT han_i = (1-coef)*hankel_i_tab[ind] + coef*hankel_i_tab[ind+1];

	hr = han_r*ar - han_i*ai;
	hi = han_r*ai + han_i*ar;
	if (k*sqrt(ai*ai + ar*ar) < 1) {
	  dispxr = hi*kx;
	  dispxi = -hr*kx;
	  dispyr = hi*ky;
	  dispyi = -hr*ky;
	}
	hr *= damp;
	hi *= damp;
	hr *= mod;
	hi *= mod;
	//  }
	 
	 
      }
      	   
#ifdef PLOT_RESULT
      if (s < nb_sources_input) {
	heights[2*i] += hr;
	heights[2*i+1] += hi;
      } else {
	heights[2*nb_rows*nb_cols + 2*i] += hr;
	heights[2*nb_rows*nb_cols + 2*i+1] += hi;
      }
      displacement[4*i] = 0;
      displacement[4*i+1] = 0;
      displacement[4*i+2] = 0;
      displacement[4*i+3] = 0;
#else
      heights[2*i] += hr;
      heights[2*i+1] += hi;
      displacement[4*i] += 0.5*dispxr;
      displacement[4*i+1] += 0.5*dispxi;
      displacement[4*i+2] += 0.5*dispyr;
      displacement[4*i+3] += 0.5*dispyi;
#endif
    }
  }
}

 
__global__
void updateHeight(FLOAT * heights, FLOAT* time_height,
		  FLOAT * displacement, FLOAT* time_displacement,
		  FLOAT omega, uint n_nodes, FLOAT time) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (uint i = index; i < n_nodes; i += stride) {
    FLOAT hr = heights[2*i];
    FLOAT hi = heights[2*i+1];
    FLOAT dispxr = displacement[4*i];
    FLOAT dispxi = displacement[4*i+1];
    FLOAT dispyr = displacement[4*i+2];
    FLOAT dispyi = displacement[4*i+3];
#ifdef PLOT_RESULT
    hr += heights[2*n_nodes + 2*i];
    hi += heights[2*n_nodes + 2*i+1];
#endif
    time_height[i] += hr*cos(-omega*time) - hi*sin(-omega*time);
    time_displacement[2*i] += dispxr*cos(-omega*time) - dispxi*sin(-omega*time);
    time_displacement[2*i+1] += dispyr*cos(-omega*time) - dispyi*sin(-omega*time);
  }
}


CudaWaterSurface::CudaWaterSurface() {
  nb_wl = 0;
}

CudaWaterSurface::~CudaWaterSurface() {}

void CudaWaterSurface::clear() {
#ifdef INTERACTIVE_
  cudaFree(heights);
#endif
  for (uint w = 0; w < nb_wl; ++w) {
    cudaFree(amplitudes[w]);
    cudaFree(indexes[w]);
    cudaFree(positions[w]);
#ifdef INTERACTIVE_
    cudaFree(is_active[w]);
#else
    cudaFree(heights[w]);
#endif
  }
#ifdef PROJECTED_GRID
  cudaFree(positions_grid);
  cudaFree(sizes);
#endif
  cudaFree(hankel_r_tab);
  cudaFree(hankel_i_tab);
  INFO("CLEAR CUDA");
}

 
#ifndef INTERACTIVE_

void CudaWaterSurface::init(uint nw) {
  nb_wl = nw;
  amplitudes = std::vector<FLOAT*>(nb_wl);
  indexes = std::vector<FLOAT*>(nb_wl);
  positions = std::vector<FLOAT*>(nb_wl);
  nb_sources = std::vector<uint>(nb_wl);
  nb_sources_input = std::vector<uint>(nb_wl);
  wave_lenghts = std::vector<FLOAT>(nb_wl);
  heights = std::vector<FLOAT*>(nb_wl);
  displacement = std::vector<FLOAT*>(nb_wl);

  cudaMallocManaged(&time_heights, n_cols_*n_rows_*sizeof(FLOAT));
  cudaMallocManaged(&time_displacement, 2*n_cols_*n_rows_*sizeof(FLOAT));
#ifdef PROJECTED_GRID
  cudaMallocManaged(&positions_grid, 2*n_cols_*n_rows_*sizeof(FLOAT));
  cudaMallocManaged(&sizes, n_cols_*n_rows_*sizeof(FLOAT));
#else
  positions_grid = NULL;
  sizes = NULL;
#endif
  createTabs();
}
 

void CudaWaterSurface::allocMem(uint wl, uint ns, uint na) {
  cudaMallocManaged(&displacement[wl], 4*n_cols_*n_rows_*sizeof(FLOAT));
#ifdef PLOT_RESULT
  cudaMallocManaged(&heights[wl], 4*n_cols_*n_rows_*sizeof(FLOAT));
  for (uint i = 0; i < 4*n_rows_*n_cols_; ++i) {
#else
  cudaMallocManaged(&heights[wl], 2*n_cols_*n_rows_*sizeof(FLOAT));
  for (uint i = 0; i < 2*n_rows_*n_cols_; ++i) {
#endif
    heights[wl][i] = 0;
		    
  }
  for (uint i = 0; i < 4*n_rows_*n_cols_; ++i) {
    displacement[wl][i] = 0;
  }

  cudaMallocManaged(&amplitudes[wl], 2*(ns+na)*sizeof(FLOAT));
  cudaMallocManaged(&indexes[wl], (ns+na)*sizeof(FLOAT));
  cudaMallocManaged(&positions[wl], 2*(ns+na)*sizeof(FLOAT));
  nb_sources[wl] = ns;
  nb_sources_input[wl] = na;
}

void CudaWaterSurface::setHeight(int nb_in_waves) {
  
  for (uint w = 0; w < nb_wl; ++w) {
    
#ifdef PLOT_RESULT
    uint nb_f = 4*n_rows_*n_cols_;
#else
    uint nb_f = 2*n_rows_*n_cols_;
#endif  
    for (uint i = 0; i < nb_f; ++i) {
      heights[w][i] = 0;
    }
    for (uint i = 0; i < 4*n_rows_*n_cols_; ++i) {
      displacement[w][i] = 0;
    }
    FLOAT wl = wave_lenghts[w];
    FLOAT k = 2*M_PI/wl;
    FLOAT omega = angular_vel(k);
    FLOAT v = velocity(k);
    int blockSize = 256;
    int numBlocks = ( n_rows_*n_cols_ + blockSize - 1) / blockSize;

    uint ns = nb_sources[w];
    uint na = nb_sources_input[w];
    addHeight0<<<numBlocks,blockSize>>>
      (heights[w], displacement[w], amplitudes[w], indexes[w], positions[w], 
       positions_grid, sizes, ns, na,
       k, omega, v, n_rows_, n_cols_,
       damping_, scale_,
       ui_parameters::show_in_field , ui_parameters::show_scattered_field,
       hankel_r_tab, hankel_i_tab);
    cudaDeviceSynchronize();
  }
  
}
 
void CudaWaterSurface::setTimeHeight(int time) {
  for (uint i = 0; i < n_rows_*n_cols_; ++i) {
    time_heights[i] = 0;
  }
  for (uint i = 0; i < 4*n_rows_*n_cols_; ++i) {
    time_displacement[i] = 0;
  }

  for (uint w = 0; w < nb_wl; ++w) {
    FLOAT wl = wave_lenghts[w];
    FLOAT k = 2*M_PI/wl;
    FLOAT omega = angular_vel(k);
    int blockSize = 256;
    int numBlocks = ( n_rows_*n_cols_ + blockSize - 1) / blockSize;
      
    updateHeight<<<numBlocks,blockSize>>>(heights[w], time_heights,
					  displacement[w], time_displacement,
					  omega, n_cols_*n_rows_, time*dt_);
    cudaDeviceSynchronize();
  }
}

#else

void CudaWaterSurface::init(uint nw) {
  nb_wl = nw;
  amplitudes = std::vector<FLOAT*>(nb_wl);
  indexes = std::vector<FLOAT*>(nb_wl);
  positions = std::vector<FLOAT*>(nb_wl);
  nb_sources = std::vector<uint>(nb_wl);
  nb_sources_input = std::vector<uint>(nb_wl);
  wave_lenghts = std::vector<FLOAT>(nb_wl);
  is_active = std::vector<bool *>(nb_wl);

#ifdef PLOT_RESULT
  uint h_size = 2*n_cols_*n_rows_;
#else
  uint h_size = n_cols_*n_rows_;
#endif
 cudaMallocManaged(&heights, h_size*sizeof(FLOAT));
 for (uint i = 0; i < h_size; ++i) {
   heights[i] = 0;
 }
 cudaMallocManaged(&displacement, 2*n_cols_*n_rows_*sizeof(FLOAT));
 for (uint i = 0; i < 2*n_rows_*n_cols_; ++i) {
   displacement[i] = 0;
 }
#ifdef PROJECTED_GRID
 cudaMallocManaged(&positions_grid, 2*n_cols_*n_rows_*sizeof(FLOAT));
 cudaMallocManaged(&sizes, n_cols_*n_rows_*sizeof(FLOAT));
#else
  positions_grid = NULL;
  sizes = NULL;
#endif
  createTabs();
  cudaMallocManaged(&buffer, nb_profil*sizeof(FLOAT));
}

void CudaWaterSurface::allocMem(uint wl, uint ns, uint na) {
  cudaError_t alloc_ok = cudaMallocManaged(&amplitudes[wl], 2*(ns+na)*size_tmp*sizeof(FLOAT));
  cudaMallocManaged(&indexes[wl], (ns+na)*sizeof(FLOAT));
  cudaMallocManaged(&is_active[wl], (ns+na)*sizeof(bool));
  cudaMallocManaged(&positions[wl], 2*(ns+na)*sizeof(FLOAT));
  nb_sources[wl] = ns;
  nb_sources_input[wl] = na;
}


void CudaWaterSurface::setHeight(int time) {
  for (uint i = 0; i < n_rows_*n_cols_; ++i) {
    heights[i] = 0;
  }
  for (uint i = 0; i < 2*n_rows_*n_cols_; ++i) {
    displacement[i] = 0;
  }
	
  for (uint w = 0; w < nb_wl; ++w) {
    FLOAT wl = wave_lenghts[w];
    FLOAT k = 2*M_PI/wl;
    FLOAT omega = angular_vel(k);
    int blockSize = 64;
    int numBlocks = ( n_rows_*n_cols_ + blockSize - 1) / blockSize;
    FLOAT v = velocity(k);

    addHeight<<<numBlocks,blockSize>>>
      (heights, displacement, amplitudes[w], indexes[w],
       positions[w],  is_active[w],
       positions_grid, sizes, 
       nb_sources[w], nb_sources_input[w],
       k, omega, v, time*dt_, dt_, n_rows_, n_cols_,
       size_tmp, ampli_steps[w], damping_, scale_,
       ui_parameters::show_in_field, ui_parameters::show_scattered_field,
       hankel_r_tab, hankel_i_tab);
     cudaDeviceSynchronize();
    addHeightM<<<numBlocks,blockSize>>>
      (heights, displacement, amplitudes_m[w], indexes_m[w],
       positions_m[w], is_active_m[w],
       positions_grid, sizes,
       nb_sources_m[w], nb_sources_input_m[w],
       k, omega, v, time*dt_, dt_, n_rows_, n_cols_,
       size_tmp, ampli_steps[w], damping_, scale_,
       ui_parameters::show_in_field, ui_parameters::show_scattered_field,
       hankel_r_tab, hankel_i_tab);
    cudaDeviceSynchronize();
  }

}


void CudaWaterSurface::initM(uint nw) {
  amplitudes_m = std::vector<FLOAT*>(nb_wl);
  indexes_m = std::vector<FLOAT*>(nb_wl);
  positions_m = std::vector<FLOAT*>(nb_wl);
  nb_sources_m = std::vector<uint>(nb_wl);
  nb_sources_input_m = std::vector<uint>(nb_wl);
  is_active_m = std::vector<bool *>(nb_wl);
}

void CudaWaterSurface::allocMemM(uint wl, uint ns, uint na) {
  cudaError_t alloc_ok = cudaMallocManaged(&amplitudes_m[wl], 2*(ns+na)*size_tmp*sizeof(FLOAT));
  cudaMallocManaged(&indexes_m[wl], (ns+na)*sizeof(FLOAT));
  cudaMallocManaged(&is_active_m[wl], (ns+na)*sizeof(bool));
  cudaMallocManaged(&positions_m[wl], 2*(ns+na)*size_tmp*sizeof(FLOAT));
  nb_sources_m[wl] = ns;
  nb_sources_input_m[wl] = na;
}


      

      
#endif


void CudaWaterSurface::createTabs() {
  cudaMallocManaged(&hankel_r_tab, nb_profil*sizeof(FLOAT));
  cudaMallocManaged(&hankel_i_tab, nb_profil*sizeof(FLOAT));

  for (uint i = 0; i < settings::nb_profil; ++i) {
    COMPLEX h = settings::hankel_tab[i];
    hankel_r_tab[i] = real(h);
    hankel_i_tab[i] = imag(h);
  }
  for (uint i = 0; i < 1; ++i) {
    hankel_r_tab[i] = hankel_r_tab[1];
    hankel_i_tab[i] = hankel_i_tab[1];
  }
}
			
