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

#ifndef BLOCK_3D_DCT_H
#define BLOCK_3D_DCT_H

#include <fftw3.h>

#include "setting.h"
// Frequency transform to velocity.
void block_DCT_F2V(const int xRes, const int yRes, const int zRes,
               const int xN, const int yN, const int zN, Real* tempZ_, Real* tempY_, Real* tempX_,
               fftw_r2r_kind kindX, fftw_r2r_kind kindY, fftw_r2r_kind kindZ);

// Velocity transform to Frequency.
void block_DCT_V2F(const int xRes, const int yRes, const int zRes,
               const int xN, const int yN, const int zN, Real* tempX_, Real* tempY_, Real* tempZ_,
               fftw_r2r_kind kindX, fftw_r2r_kind kindY, fftw_r2r_kind kindZ);
#endif  // BLOCK_3D_DCT_H
