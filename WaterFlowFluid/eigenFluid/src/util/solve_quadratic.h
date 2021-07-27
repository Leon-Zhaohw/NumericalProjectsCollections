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

#ifndef SOLVE_QUADRATIC_H
#define SOLVE_QUADRATIC_H

// #include <glog/logging.h>
#include <math.h>

// Solve Ax^2 + bx + c = 0
void SolveQuadratic(const double A, const double b, const double c,
                    double* root_plus, double*root_minus);

#endif  // SOLVE_QUADRATIC_H
