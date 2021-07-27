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
// #include <glog/logging.h>
#include "Eigen"
#include <vector>

#include "2D/laplacian_basis_2D.h"
#include "setting.h"
#include "util/util.h"

void LaplacianBasis2D::MultiplyweightTensor(std::vector<Adv_Tensor_Type>& C, const Eigen::VectorXd& weights) {
  const int basis_dim_= root_num_basis_*root_num_basis_;
  
  if (C.size() != basis_dim_) {std::cout << "laplacian_basis_2D.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  if (weights.size() != basis_dim_) {std::cout << "laplacian_basis_2D.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  Eigen::MatrixXd WW(basis_dim_, basis_dim_);
  WW = weights * weights.transpose();
  
  for (int k = 0; k < basis_dim_; k++) {
    Adv_Tensor_Type Cknew = weights[k] * C[k].cwiseProduct(WW);
    Cknew.makeCompressed();
    C[k] = Cknew;
  }
}

void LaplacianBasis2D::ElimiteDiagElements(std::vector<Adv_Tensor_Type>& C) {
  const int num_basis_ = C.size();
  for (int k = 0; k < num_basis_; k++) {
    for (int j = 0; j < num_basis_; j++) {
      C[k].coeffRef(j,j) = 0;
    }
    C[k].makeCompressed();
  }
}

void LaplacianBasis2D::VerifyAntisymmetric(const std::vector<Adv_Tensor_Type>& C) {
  const int num_basis_ = root_num_basis_*root_num_basis_;
  double max_abs_entries = FLT_MIN;
  int num_non_zero = 0;
  
  for (int k = 0; k < num_basis_; k++) {
    for (int j = 0; j < num_basis_; j++) {
      for (int i = 0; i < num_basis_; i++) {
        double Cijk = AccessMatrix(C[k],i,j);
        if (Cijk != 0) {
          num_non_zero ++;
          double Cikj = AccessMatrix(C[j],i,k);
          if ((Cijk + Cikj) > 1e-10) {
            std::cout <<  "Non symeetric entries found." << std::endl;
          }
          if (std::abs(Cijk) > max_abs_entries) {
            max_abs_entries = std::abs(Cijk);
          }
        }
      }
    }
  }
  std::cout << "Number of non-zero elements in tensor: " << num_non_zero
      << " Maximum abs entry: " << max_abs_entries;
}
