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

#include <cmath>
#include <float.h>
#include "Eigen"
#include <vector>

#include "3D/laplacian_basis_set_3d.h"
#include "setting.h"
#include "util/util.h"

void LaplacianBasisSet3D::MultiplyweightTensor(std::vector<Adv_Tensor_Type>& C, const Eigen::VectorXd& weights) {
  if (C.size() != basis_dim_) {std::cout << "laplacian_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  if (weights.size() != basis_dim_) {std::cout << "laplacian_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  Eigen::MatrixXd WW(basis_dim_, basis_dim_);
  WW = weights * weights.transpose();
  
  for (int k = 0; k < basis_dim_; k++) {
    Adv_Tensor_Type Cknew = weights[k] * C[k].cwiseProduct(WW);
    Cknew.makeCompressed();
    C[k] = Cknew;
  }
}

void LaplacianBasisSet3D::VerifyAntisymmetric(const std::vector<Adv_Tensor_Type>& C) {
  int num_non_zero = 0;
  double max_entries = FLT_MIN;
  for (int k = 0; k < basis_dim_; k++) {
    for (int j = 0; j < basis_dim_; j++) {
      for (int i = 0; i < basis_dim_; i++) {
        double Cijk = AccessMatrix(C[k],i,j);
        if (Cijk != 0) {
          num_non_zero ++;
          double Cikj = AccessMatrix(C[j],i,k);
          if ((Cijk + Cikj) > 1e-10) {
            std::cout <<  "Non symeetric entries found." << std::endl;
          }
          if (std::abs(Cijk) > max_entries) {
            max_entries = std::abs(Cijk);
          }
        }
      }
    }
  }
  std::cout << "Number of non-zeros in tensor : " << num_non_zero
      << "  Maximum abs tensor entry value: " << max_entries;
}
