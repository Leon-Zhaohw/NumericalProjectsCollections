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

#ifndef DIRICHLET_BASIS_3D_H
#define DIRICHLET_BASIS_3D_H

#include "Eigen"
// #include <glog/logging.h>

#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/laplacian_basis_3d.h"

class DirichletBasis3D :public LaplacianBasis3D {
public:
  DirichletBasis3D(const int k1, const int k2, const int k3,
                   const double A, const double B, const double C):
                   LaplacianBasis3D(k1,k2,k3,A,B,C) {
                     // if ((A_*k1_ + B_*k2_ + C_*k3_ ) != 0) {std::cout << "dirichlet_basis_3d.h " << __LINE__ << " FATAL: "  << "Basis is not divergence free." << std::endl; exit(0);}
                     if (std::abs(A_*k1_ + B_*k2_ + C_*k3_) > 1e-13) {
                       std::cout << "dirichlet_basis_3d.h " << __LINE__ << " FATAL: " <<  "Basis is not divergence free." << std::endl; exit(0);
                     }
                     Normalize();
                   }
  ~ DirichletBasis3D(){}
  void DiscretizeAdd(const double coef, VECTOR3_FIELD_3D* vfield) override;
  void Normalize() override;
  double GetLambdaSquared() const override {
    return static_cast<double>(k1_*k1_ + k2_*k2_ + k3_*k3_);
  }
  double ComputeBasisAt(const int xpos, const int ypos,
                        const int zpos, const int mode,
                        const Eigen::Vector3d& dxyz) const override;
protected:
};

#endif  // DIRICHLET_BASIS_3D_H
