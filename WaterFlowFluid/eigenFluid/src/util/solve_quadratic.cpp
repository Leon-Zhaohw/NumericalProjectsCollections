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

#include "util/solve_quadratic.h"
#include <iostream>
// Solve Ax^2 + bx + c = 0
void SolveQuadratic(const double A, const double b, const double c,
                    double* root_plus, double*root_minus) {
  const double Delta = b*b - 4*A*c;
  if (Delta < 0) {
    std::cout <<  "Imaginary roots." << std::endl;
    *root_plus = 0;
    *root_minus = 0;
    return;
  } else {
    *root_plus = 0.5*(-b + sqrt(Delta)) / A;
    *root_minus = 0.5*(-b - sqrt(Delta)) / A;
  }
}
