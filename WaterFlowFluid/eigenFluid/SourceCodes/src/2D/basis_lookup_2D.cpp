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

#include <stdio.h>
#include "2D/basis_lookup_2D.h"

int BasisLookUp::LookUpBasis(int i, int mode) {
 if (i > basis_num_ || i < 0) {
   printf("Invalid range.");
   return -1;
 }
 // Look up X.
 if (mode == 0) {
   return i % basis_num_root_ + 1;
 }
 // Look up Y.
 else if (mode == 1) {
   return i / basis_num_root_ + 1;
 }
 else {
   printf("Invalid mode.");
   return -1;
 }
}

int BasisLookUp::InverseLookUp(int i1, int i2) {
  if (i1 <= 0 || i2 <= 0 || i1 > basis_num_root_ || i2 > basis_num_root_) {
    return -1;
  }
  return (i1 - 1) + (i2 - 1)*basis_num_root_;
}
