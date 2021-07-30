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

#ifndef SETTING_H
#define SETTING_H

#include "Eigen"
#include "Sparse"

#ifndef Real
#define Real double
#endif

//#if Real==float 
// #define SINGLE_PRECISION true
//#else 
// #define SINGLE_PRECISION false
//#endif

#ifndef Adv_Tensor_Sparse
#define Adv_Tensor_Sparse
#endif

typedef Eigen::SparseMatrix<double> Adv_Tensor_Type;
// typedef Eigen::MatrixXd Adv_Tensor_Type;

#ifndef USE_ROOT_LAMBDA
#define USE_ROOT_LAMBDA
#endif

//#define APPLE

//#ifndef max
//#define max(a,b) (((a) > (b)) ? (a) : (b))
//#endif

//#ifndef min
//#define min(a,b) (((a) < (b)) ? (a) : (b))
//#endif

#endif
