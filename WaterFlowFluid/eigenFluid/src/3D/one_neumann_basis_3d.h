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

#ifndef ONE_NEUMANN_BASIS_3D_H
#define ONE_NEUMANN_BASIS_3D_H

#include "Eigen"
#include <cmath>
// #include <glog/logging.h>

#include "3D/VECTOR3_FIELD_3D.h"
#include "3D/laplacian_basis_3d.h"

// For Neuman basis, the actual wavenumber for x is integer - 0.5., here k1_ equals
// To actuall wavenumber + 0.5.
// The x compoenent can be zero.
class OneNeumannBasis3D :public LaplacianBasis3D {
public:
  OneNeumannBasis3D(const int k1, const int k2, const int k3,
                    const double A, const double B, const double C):
                    LaplacianBasis3D(k1,k2,k3,A,B,C) {
    // Check the divergence.
    // x component is zero.
    if (k1_ == 0) {
      if (std::abs(k2_*B + k3_*C_) > 1e-14) {
        std::cout << "one_neumann_basis_3d.h " << __LINE__ << " FATAL: " <<  "Divergence is non-zero." << std::endl; exit(0);
      }
    } else {  // x component is non-zero.
      if (std::abs((k1_ - 0.5)*A + k2_*B_ + k3_*C_ ) > 1e-14) {
        std::cout << "one_neumann_basis_3d.h " << __LINE__ << " FATAL: " <<  "Divergence is non-zero." << std::endl; exit(0);
      }
    }
    Normalize();
  }
  ~OneNeumannBasis3D(){}
 
  void DiscretizeAdd(const double coef, VECTOR3_FIELD_3D* vfield) override;
  void Normalize() override;
  double GetLambdaSquared() const override {
    if (k1_ == 0) {
      return static_cast<double>(k2_*k2_ + k3_*k3_);
    } else {
      return ((k1_ - 0.5)*(k1_ - 0.5) + k2_*k2_ + k3_*k3_);
    }
  }
  double ComputeBasisAt(const int xpos, const int ypos,
                        const int zpos, const int mode,
                        const Eigen::Vector3d& dxyz) const override;
protected:
};

#endif
