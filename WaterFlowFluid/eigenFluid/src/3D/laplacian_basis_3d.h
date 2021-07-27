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

#ifndef LAPLACIAN_BASIS_3D_H
#define LAPLACIAN_BASIS_3D_H

#include "Eigen"
#include "3D/VECTOR3_FIELD_3D.h"
#include "Alg/VEC3.h"

class LaplacianBasis3D{
public:
  LaplacianBasis3D(const int k1, const int k2, const int k3, 
                   const double A, const double B, const double C):
                   k1_(k1),k2_(k2),k3_(k3),A_(A),B_(B),C_(C) {
    invNorm_ = 0; // Default zero.
    // Initialize DCT norm.
    int num_zeros = 0;
    if (k1 == 0) num_zeros ++;
    if (k2 == 0) num_zeros ++;
    if (k3 == 0) num_zeros ++;
    
    if (num_zeros == 0) {
      DCTNorm_ = 0.125;
    } else if (num_zeros == 1) {
      DCTNorm_ = 0.25;
    } else if (num_zeros == 2) {
      DCTNorm_ = 0.5;
    } else {
      DCTNorm_ = 1.0;
    }
  }
  ~LaplacianBasis3D(){}
  virtual void DiscretizeAdd(const double coef, VECTOR3_FIELD_3D* vfield) = 0;
  virtual void Normalize() = 0;
  virtual double ComputeBasisAt(const int xpos, const int ypos,
                                const int zpos, const int mode,
                                const Eigen::Vector3d& dxyz) const = 0;
  VEC3I GetWaveNum() const {
    return VEC3I(k1_, k2_, k3_); 
  }
  // Get the const coeff for vx.
  double GetVxCoef() const {
    return invNorm_*A_;
  }
  // Get the const coeff for vy.
  double GetVyCoef() const {
    return invNorm_*B_;
  }
  // Get the const coeff for vz.
  double GetVzCoef() const {
    return invNorm_*C_;
  }
  double GetDCTNorm() const {
    return DCTNorm_;
  }
  double GetInvNorm() const {
    return invNorm_;
  }
  double GetA () const {return A_;}
  double GetB () const {return B_;}
  double GetC () const {return C_;}
  virtual double GetLambdaSquared() const = 0;

protected:
  // Wavenumber for this basis.
  const int k1_;
  const int k2_;
  const int k3_;
  // Coefficient.
  const double A_;
  const double B_;
  const double C_;
  // Normalization factor.
  double invNorm_;
  double DCTNorm_;
};

#endif  // LAPLACIAN_BASIS_3D_H
