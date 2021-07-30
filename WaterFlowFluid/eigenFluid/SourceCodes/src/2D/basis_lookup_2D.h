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

#ifndef BASIS_LOOKUP_H
#define BASIS_LOOKUP_H

// The index of basis arranged in 1D begin with 0.
// The index of basis arranged in pair of integers begin with 1.
class BasisLookUp {
public:
  BasisLookUp(int basis_num_root) : basis_num_root_(basis_num_root)
  { basis_num_ = basis_num_root_ * basis_num_root;};
  ~BasisLookUp(){};
  // Put int the index of the coieffcients of the basis arranged in 1D,
  // return the pair of integers. For example, 1D basis index [0-63],
  // return pair of integers, (i1, i2), i1 [1-8], i2 [1-8]. X axis associated
  // with first integer, Y axis associated with second integer. 
  int LookUpBasis(int i, int mode);
  // Do a inverse lookup of the basis. Give a pair of integers, return the
  // index of the basis arranged in 1D.
  int InverseLookUp(int i1, int i2);
private:
  int basis_num_root_;
  int basis_num_;
};

#endif
