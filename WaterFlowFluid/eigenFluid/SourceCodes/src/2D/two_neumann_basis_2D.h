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

#ifndef TWO_NEUMANN_BASIS_2D_H
#define TWO_NEUMANN_BASIS_2D_H

#include "Eigen"
// #include <glog/logging.h>
#include <math.h>

#include "2D/VFIELD2D.h"
#include "Alg/VEC2.h"

class TwoNeumannBasis2D {
public:
  TwoNeumannBasis2D(const int k1, const int k2, const double A, const double B):
      k1_(k1), k2_(k2), A_(A), B_(B) {
    invNorm_ = 0; // Default zero.
    // Check for divergence free.
    if (std::abs(-A_*k1_ + B_*k2_ ) > 1e-13 && k1_ != 0) {
      std::cout << "two_neumann_basis_2D.h " << __LINE__ << " FATAL: " <<  "Divergence is non-zero." << std::abs(-A_*k1_ + B_*k2_)   << std::endl; exit(0);
    }
    // normalize.
    Normalize();

    // Initialize DCT norm.
    int num_zeros = 0;
    if (k1 == 0) num_zeros ++;
    if (k2 == 0) num_zeros ++;
    if (num_zeros == 0) {
      DCTNorm_ = 0.25;
    } else if (num_zeros == 1) {
      DCTNorm_ = 0.5;
    } else {
      DCTNorm_ = 1.0;
    }
    weight_ = 1.0;
  }
  
  TwoNeumannBasis2D():k1_(0), k2_(0), A_(0), B_(0) {
    invNorm_ = 0.0;
  }
  
  ~TwoNeumannBasis2D(){}
  void Normalize();

  // Get the const coeff for vx.
  double GetVxCoef() const {
    return invNorm_*A_;
  }
  // Get the const coeff for vy.
  double GetVyCoef() const {
    return invNorm_*B_;
  }
  
  VEC2I GetWaveNum() const {
    return VEC2I(k1_, k2_); 
  }
  double getInvNorm() const {
    return invNorm_;
  }
  // The domain is [0,1]^2. This actually returns the *negative* of actual
  // Laplacian eigen values.
  double GetLambdaSquared() const {
    return k1_*k1_ + k2_*k2_;
  }
  void DiscretizeAdd(const double coef, VFIELD2D* vfield);
  double ComputeBasisAt(const int xpos, const int ypos,
                        const int mode, const VEC2& dxy) const;
  double GetDCTNorm() const {
    return DCTNorm_;
  }
  
  double GetA () const {return A_;}
  double GetB () const {return B_;}
  double GetInvNorm() const {
    return invNorm_;
  }
  double getWeight() const {
    return weight_;
  }
protected:
  // Wavenumber for this basis.
  int k1_;
  int k2_;
  // Coefficient.
  double A_;
  double B_;
  // Normalization factor.
  double invNorm_;
  double DCTNorm_;
  
  double weight_;
};

#endif  // TWO_NEUMANN_BASIS_2D_H
