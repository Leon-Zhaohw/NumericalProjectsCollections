/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2018 Qiaodong Cui (qiaodong@ucsb.edu)
 */

#include <fftw3.h>
#include <omp.h>
#include "util/block_3d_dct.h"
#include <iostream>
/*
 *           int rank, const int *n,
			 int howmany,
			 R *in, const int *inembed,
			 int istride, int idist,
			 R *out, const int *onembed,
			 int ostride, int odist,
			 const X(r2r_kind) * kind, unsigned flags)
 */
namespace {

void CountNonZero(Real* x, const int xsize) {
  int nonzeros = 0;
  for (int i = 0; i < xsize; i++) {
    if (x[i] != 0) {
      nonzeros ++;
    }
  }
  std::cout <<  "Number of non-zeros: " << nonzeros << std::endl;
}
  
}  // namespace

void block_DCT_F2V(const int xRes, const int yRes, const int zRes,
               const int xN, const int yN, const int zN, Real* zF, Real* yF, Real* xF,
               fftw_r2r_kind kindX, fftw_r2r_kind kindY, fftw_r2r_kind kindZ) {
  
  // First along z axis.
  int xSlab, ySlab, zSlab;
  int n[] = {zRes};
  int idist = zRes;
  int odist = zRes;
  int istride = 1; int ostride = 1;
  int *inembed = n, *onembed = n;
  int howmany = xN*yN;
  
  fftw_plan plan_z = fftw_plan_many_r2r(1, n, howmany, zF, inembed, istride, idist,
                                        zF, onembed, ostride, odist, &kindZ, FFTW_MEASURE);
  fftw_execute(plan_z);
  ySlab = zRes*xN; 
  zSlab = yRes*xN;
  // Then along y axis.
#pragma omp parallel for
  for (int z = 0; z < zRes; z++) {
    for (int x = 0; x < xN; x++) {
      for (int y = 0; y < yN; y++) {
        const int idxZ = z + x*zRes + y*ySlab;
        const int idxY = y + x*yRes + z*zSlab;
        yF[idxY] = zF[idxZ];
      }
    }
  }
  n[0] = yRes;
  idist = yRes; odist = yRes;
  howmany = xN*zRes;
  
  fftw_plan plan_y = fftw_plan_many_r2r(1, n, howmany, yF, inembed, istride, idist,
                                        yF, onembed, ostride, odist, &kindY, FFTW_MEASURE);
  fftw_execute(plan_y);
  
  // Final along x axis.
  xSlab = xRes*yRes;
  
#pragma omp parallel for
  for (int z = 0; z < zRes; z++) {
    for (int y = 0; y < yRes; y++) {
      for (int x = 0; x < xN; x++) {
        const int idxX = x + y*xRes + z*xSlab;
        const int idxY = y + x*yRes + z*zSlab;
        xF[idxX] = yF[idxY];
      }
    }
  }
  
  n[0] = xRes;
  idist = xRes; odist = xRes;
  howmany = yRes*zRes;
  
  fftw_plan plan_x = fftw_plan_many_r2r(1, n, howmany, xF, inembed, istride, idist,
                                        xF, onembed, ostride, odist, &kindX, FFTW_MEASURE);
  fftw_execute(plan_x);
}

// Velocity transform to Frequency.
void block_DCT_V2F(const int xRes, const int yRes, const int zRes,
               const int xN, const int yN, const int zN, Real* xF, Real* yF, Real* zF,
               fftw_r2r_kind kindX, fftw_r2r_kind kindY, fftw_r2r_kind kindZ) {
  // First along X direction. size(xF) = xRes*yRes*zRes;
  // X is passed in, the memory stride along x is 1.
  int xSlab, ySlab, zSlab;
  int n[] = {xRes};
  int idist = xRes;
  int odist = xRes;
  int istride = 1; int ostride = 1;
  int *inembed = n, *onembed = n;
  int howmany = yRes*zRes;
  
  fftw_plan plan_x = fftw_plan_many_r2r(1, n, howmany, xF, inembed, istride, idist, 
                                        xF, onembed, ostride, odist, &kindX, FFTW_MEASURE);
  fftw_execute(plan_x);
  
  xSlab = xRes*yRes;
  ySlab = yRes*xN;
  // Then along y direction. size(yF) = xN*yRes*zRes.
#pragma omp parallel for
  for (int z = 0; z < zRes; z++) {
    for (int y = 0; y < yRes; y++) {
      for (int x = 0; x < xN; x++) {
        const int idxX = x + y*xRes + z*xSlab;
        const int idxY = y + x*yRes + z*ySlab; 
        yF[idxY] = xF[idxX];
      }
    }
  }
  
  n[0] = yRes;
  idist = yRes; odist = yRes;
  howmany = xN*zRes;
  
  fftw_plan plan_y = fftw_plan_many_r2r(1, n, howmany, yF, inembed, istride, idist,
                                        yF, onembed, ostride, odist, &kindY, FFTW_MEASURE);
  fftw_execute(plan_y);
  
  zSlab = xN*zRes;
  // Finally along z direction. size(zF) = xN*yN*zRes.
#pragma omp parallel for
  for (int z = 0; z < zRes; z++) {
    for (int y = 0; y < yN; y++) {
      for (int x = 0; x < xN; x++) {
        const int idxZ = z + x*zRes + y*zSlab;
        const int idxY = y + x*yRes + z*ySlab;
        zF[idxZ] = yF[idxY];
      }
    }
  }
  
  n[0] = zRes;
  idist = zRes; odist = zRes;
  howmany = xN*yN;
  
  fftw_plan plan_z = fftw_plan_many_r2r(1, n, howmany, zF, inembed, istride, idist,
                                        zF, onembed, ostride, odist, &kindZ, FFTW_MEASURE);
  
  fftw_execute(plan_z);
  // The result is along z direction.
}

/*
void block_DCT(const int xRes, const int yRes, const int zRes,
               const int xN, const int yN, const int zN, Real* xF, Real* yF, Real* zF,
               fftw_r2r_kind kindX, fftw_r2r_kind kindY, fftw_r2r_kind kindZ) {
  
  int xSlab, ySlab, zSlab;
  int n[] = {xRes};
  int idist = xRes;
  int odist = xRes;
  int istride = 1; int ostride = 1;
  int *inembed = n, *onembed = n;
  int howmany = yN*zN;
  
  fftw_plan plan_x = fftw_plan_many_r2r(1, n, howmany, xF, inembed, istride, idist,
                                        xF, onembed, ostride, odist, &kindX, FFTW_MEASURE);
  fftw_execute(plan_x);
  
  
  xSlab = xRes*yN;
  ySlab = xRes*yRes;
  
#pragma omp parallel for
  for (int z = 0; z < zN; z++) {
    for (int x = 0; x < xRes; x++) {
      for (int y = 0; y < yN; y++) {
        const int idxX = x + y*xRes + z*xSlab;
        const int idxY = y + x*yRes + z*ySlab;
        yF[idxY] = xF[idxX];
      }
    }
  }

  
  n[0] = yRes;
  idist = yRes; odist = yRes;
  howmany = zN*xRes;
  fftw_plan plan_y = fftw_plan_many_r2r(1, n, howmany, yF, inembed, istride, idist,
                                        yF, onembed, ostride, odist, &kindY, FFTW_MEASURE);
  fftw_execute(plan_y);
  
  zSlab = xRes*zRes;  
#pragma omp parallel for
  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      for (int z = 0; z < zN; z++) {
        const int idxZ = z + x*zRes + y*zSlab;
        const int idxY = y + x*yRes + z*ySlab;
        zF[idxZ] = yF[idxY];
      }
    }
  }

  n[0] = zRes;
  idist = zRes; odist = zRes;
  howmany = xRes*yRes;
  fftw_plan plan_z = fftw_plan_many_r2r(1, n, howmany, zF, inembed, istride, idist, zF, onembed,
    ostride, odist, &kindZ, FFTW_MEASURE);
  fftw_execute(plan_z);
}
*/
