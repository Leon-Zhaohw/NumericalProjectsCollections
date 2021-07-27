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

#ifndef TWO_NEUMANN_X_3D_BASIS_H
#define TWO_NEUMANN_X_3D_BASIS_H

#include <cmath>
// #include <glog/logging.h>

#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/laplacian_basis_3d.h"

class TwoNeumannXBasis3D :public LaplacianBasis3D {
public:
  TwoNeumannXBasis3D(const int k1, const int k2, const int k3,
                    const double A, const double B, const double C):
                    LaplacianBasis3D(k1,k2,k3,A,B,C) {
    // Check for divergence free.
    if (std::abs(-A_*k1_ + B_*k2_ + C_*k3_) > 1e-13 && k1_ != 0) {
      std::cout << "two_neumann_x_3d_basis.h " << __LINE__ << " FATAL: " <<  "Divergence is non-zero." << std::abs(-A_*k1_ + B_*k2_ + C_*k3_)   << std::endl; exit(0);
    }
    // if (k1_ == 0), vy = 0, vz = 0, always divergence free.
    Normalize();
  }
  ~TwoNeumannXBasis3D(){}
  void DiscretizeAdd(const double coef, VECTOR3_FIELD_3D* vfield) override;
  void Normalize() override;
  double GetLambdaSquared() const override {
    return k1_*k1_ + k2_*k2_ + k3_*k3_;
  }
  double ComputeBasisAt(const int xpos, const int ypos,
                        const int zpos, const int mode,
                        const Eigen::Vector3d& dxyz) const override;
protected:
};

#endif  // TWO_NEUMANN_X_3D_BASIS_H
