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

#include "Eigen"
#include <fftw3.h>
// #include <glog/logging.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "2D/two_neumann_x_2D.h"
#include "2D/two_neumann_basis_2D.h"
#include "3D/trig_integral_3d.h"

#include "setting.h"
#include "util/util.h"
// Two neumann wall lies on the same axis of X. And two Dirichlet wall on y axis.
// Vx = root_inv_lambda*k2*cos(k1x)cos(k2y).
// Vy = root_inv_lambda*k1*sin(k1x)sin(k2y).
// k1 can goes to zero. Also, when both k1 and k2 are zero, the mode is constant
// velocity field along x axis.
// Default basis allocation strategy.

void TwoNeumannX2D::ComputeBasisWeights(Eigen::VectorXd* weights) {
  
  if (weights->size() != basis_dim_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  weights->setZero();
  
  for (int i = 0; i < basis_dim_; i++) {
    const VEC2I I = all_basis_[i].GetWaveNum();
    const int wsum = I[0]*I[0] + I[1]*I[1];
    (*weights)[i] = 1.0 - 0.9*wsum / 12484.0;
  }
  
  std::fstream out("weight.txt", std::ios::out);
  for (int i = 1; i < weights->size(); i++) {
    out << (*weights)[i] << "\n";
  }
  out.close();
}

int TwoNeumannX2D::AllocateAllBasis() {
  // need to allocate total basis_dim_ bases.
  // [0, root_basis_dim_ - 1] x [0, root_basis_dim_ - 1]
  std::vector<VEC2I> waveNs;
  for (int k2 = 0; k2 < root_basis_dim_; k2++) {
    for (int k1 = 0; k1 < root_basis_dim_; k1++) {
      VEC2I v(k1,k2);
      waveNs.push_back(v);
    }
  }
  
  std::sort(waveNs.begin(), waveNs.end(), [](const VEC2I& a, const VEC2I& b){
    return (a[0]*a[0]+a[1]*a[1]) < (b[0]*b[0]+b[1]*b[1]);
  });

  for (int i = 0; i < waveNs.size(); i++) {
    const int k1 = waveNs[i][0];
    const int k2 = waveNs[i][1];
    
    if (k1 == 0 && k2 == 0) { // zero mode.
      TwoNeumannBasis2D basis(0,0,1.0, 0.0);
      all_basis_.push_back(basis);
    } else if (k2 == 0 && k1 != 0) { // zero basis.
      TwoNeumannBasis2D basis(k1,k2,0.0,0.0);
      all_basis_.push_back(basis);
    } else if (k1 == 0 && k2 != 0) { // cos(k2y) for vx
      TwoNeumannBasis2D basis(k1,k2,1.0,0.0);
      all_basis_.push_back(basis);
    } else {
      TwoNeumannBasis2D basis(k1,k2,k2,k1);
      all_basis_.push_back(basis);
    }
  }
  if (all_basis_.size() != basis_dim_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  
  basis_allocated_ = true;
  return all_basis_.size();
}

void TwoNeumannX2D::FillBasisNumerical(const Eigen::VectorXd& coefficients,
                                             VFIELD2D* field) {
  if (coefficients.size() != basis_dim_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << " The number of coefficients" << " does not equal to number of basis." << std::endl; exit(0);}
  if (! basis_allocated_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  for (int i = 0; i < basis_dim_; i++) {
     all_basis_[i].DiscretizeAdd(coefficients[i], field);
  }
}

double TwoNeumannX2D::ComputeBasisAt(const int x_pos, const int y_pos,
                                     const int band, const int mode) const {
  if (x_pos < 0 || x_pos >= N_ || y_pos < 0 || y_pos >= N_) {
     std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: " <<  "Index out of range" << std::endl; exit(0);
  }
  if (band < 0 || band >= basis_dim_) {
    std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: " <<  "basis index is out of range: " << band << std::endl; exit(0);
  }

  return all_basis_[band].ComputeBasisAt(x_pos, y_pos, mode, dxy_);
}

void TwoNeumannX2D::InverseTramsformToVelocity(
      const Eigen::VectorXd& coefficients, VFIELD2D* field) {
  if (! basis_allocated_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  if (coefficients.size() != basis_dim_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  
  field->clear();
  ClearTemp();
  // Iterate over all the coefficients and put it into the field.
  for (int i = 0; i < basis_dim_; i ++) {
    if (coefficients[i] == 0) {
      continue;
    }
    VEC2I K_ = all_basis_[i].GetWaveNum();
    double C_x = all_basis_[i].GetVxCoef()*all_basis_[i].GetDCTNorm();
    
    // Process x field. x: REDFT01, y: REDFT01
    temp_[K_[0] + K_[1] * N_] = coefficients[i]*C_x;
  }
  
  fftw_plan plan_x_2D;
  // Process x field. x: REDFT01, y: REDFT01
  plan_x_2D = fftw_plan_r2r_2d(N_, N_, temp_, temp_, FFTW_REDFT01,
                                FFTW_REDFT01, FFTW_ESTIMATE);
  fftw_execute(plan_x_2D);  
  fftw_destroy_plan(plan_x_2D);
  // Put the value back into field.
  for (uint i = 0; i < N_*N_; i++) {
    (*field)[i][0] = temp_[i];
  }
  
  // Process vy.
  ClearTemp();
  // Iterate over all the coefficients and put it into the field.
  for (int i = 0; i < coefficients.size(); i ++) {
    if (coefficients[i] == 0) {
      continue;
    }
    
    // Get wavenumber.
    VEC2I K_ = all_basis_[i].GetWaveNum();
    double C_y = all_basis_[i].GetVyCoef()*all_basis_[i].GetDCTNorm();
    if (K_[0] < 1 || K_[1] < 1) {  // vy = 0.
      continue;
    }
    
    // Process y field. x: RODFT01, y: RODFT01
    temp_[K_[0] - 1 + (K_[1] - 1) * N_] = coefficients[i] * C_y;
  }
  // Do the transformation.
  fftw_plan plan_y_2D;
  // Inverse sine transformation on y, inverse sine transformation on x.
  plan_y_2D = fftw_plan_r2r_2d(N_, N_, temp_, temp_,
                                FFTW_RODFT01, FFTW_RODFT01,FFTW_ESTIMATE);
  fftw_execute(plan_y_2D);  
  fftw_destroy_plan(plan_y_2D);
  // Put the value back.
  for (uint i = 0; i < N_*N_; i++) {
    (*field)[i][1] = temp_[i];
  }
}

void TwoNeumannX2D::ForwardTransformtoFrequency(
    const VFIELD2D& field, Eigen::VectorXd* coefficients) {
  if (! basis_allocated_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  if (coefficients->size() != basis_dim_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  
  coefficients->setZero();
  ClearTemp();
  // Copy value to freq.
  for (int i = 0; i < N_*N_; i++) {
    temp_[i] = field[i][0];
  }
  
  // Do a transformation inplace on vx.
  fftw_plan plan_x_2D;

  // cosine transformation on X, cosine transformation on y.
  plan_x_2D = fftw_plan_r2r_2d(N_, N_, temp_, temp_, FFTW_REDFT10,
                                FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(plan_x_2D);
  fftw_destroy_plan(plan_x_2D);
  // Copy the basis back.
  for (int i = 0; i < basis_dim_; i++) {
    VEC2I K_ = all_basis_[i].GetWaveNum();
    double C_x = all_basis_[i].GetVxCoef();
    (*coefficients)[i] += temp_[K_[0] + K_[1]*N_]*inv_N_*inv_N_*0.25*PI_SQUARE*C_x;
  }
  
  // vy...
  ClearTemp();
  // Copy Y value into Temp.
  for (int i = 0; i < N_*N_; i++) {
    temp_[i] = field[i][1];
  }
  
  // Do a transformation inplace on vy.
  fftw_plan plan_y_2D;
  // Sine transformation on X, sine transformation on y.
  plan_y_2D = fftw_plan_r2r_2d(N_, N_, temp_, temp_,FFTW_RODFT10,
                                FFTW_RODFT10, FFTW_ESTIMATE);
  fftw_execute(plan_y_2D);
  fftw_destroy_plan(plan_y_2D);
  // Copy the basis back.
  for (int i = 0; i < basis_dim_; i++) {
    VEC2I K_ = all_basis_[i].GetWaveNum();
    double C_y = all_basis_[i].GetVyCoef();
    if (K_[0] < 1 || K_[1] < 1) {  // vy is zero.
      continue;
    }
    (*coefficients)[i] += temp_[K_[0] - 1 + (K_[1] - 1)*N_]*inv_N_*inv_N_*0.25*PI_SQUARE*C_y;
  }
}

void TwoNeumannX2D::ComputeVorticity(const Eigen::VectorXd& coefficients, FIELD2D* vort) {
  if (! basis_allocated_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  if (coefficients.size() != basis_dim_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  
  vort->clear();
   
  ClearTemp();
  // Iterate over all the coefficients and put it into the field.
  for (int i = 0; i < basis_dim_; i ++) {
    if (coefficients[i] == 0) {
      continue;
    }
    VEC2I K_ = all_basis_[i].GetWaveNum();
    // BK1
    double C_x = all_basis_[i].GetVyCoef()*all_basis_[i].GetDCTNorm()*static_cast<double>(K_[0]);
    // AK2
    double C_y = all_basis_[i].GetVxCoef()*all_basis_[i].GetDCTNorm()*static_cast<double>(K_[1]);
    if (K_[1] < 1) {
      continue;
    }
    
    // Process x field. x: REDFT01, y: RODFT01
    temp_[K_[0] + (K_[1] - 1) * N_] = coefficients[i]*(C_x + C_y);
  }
  
  fftw_plan plan_x_2D;
  // Process x field. x: REDFT01, y: RODFT01
  plan_x_2D = fftw_plan_r2r_2d(N_, N_, temp_, temp_, FFTW_RODFT01,
                                FFTW_REDFT01, FFTW_ESTIMATE);
  fftw_execute(plan_x_2D);  
  fftw_destroy_plan(plan_x_2D);
  // Put the value back into field.
  for (uint i = 0; i < N_*N_; i++) {
    (*vort)[i] = temp_[i];
  }
  
  ClearTemp();
}

void TwoNeumannX2D::FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) {
  if (! basis_allocated_) {std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  if (Adv_tensor != nullptr ){ Adv_tensor->clear();}
  Adv_tensor->reserve(basis_dim_);
  // Fill the tensor with zeros.
  for (int k = 0; k < basis_dim_; k++) {
    Adv_Tensor_Type Ck(basis_dim_, basis_dim_);
    Ck.setZero();
    Adv_tensor->emplace_back(Ck);
  }
  const int five_percent = basis_dim_ / 20;
  
  Eigen::VectorXd weight(basis_dim_);
  for (int i = 0; i < basis_dim_; i++) {
    weight[i] = all_basis_[i].getWeight();
  }
  // normalize weight vector.
  //double scale = sqrt(basis_dim_ / weight.squaredNorm());
  //weight *= scale;
  
  #pragma omp parallel for
  for (int k = 0; k < basis_dim_; k++) {
    if (k % five_percent == 0) {
      std::cout <<  "% 5 " << "Tensor computed." << std::endl;
    }

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;

    for (int j = 0; j < basis_dim_; j++) {
      for (int i = 0; i < basis_dim_; i++) {
        const TwoNeumannBasis2D& basis_i = all_basis_[i];
        const TwoNeumannBasis2D& basis_j = all_basis_[j];
        const TwoNeumannBasis2D& basis_k = all_basis_[k];
        const double Pz = ComputeZIntegral(basis_i, basis_j, basis_k);
        double Cijk = 0.25*(Pz);
        if (Cijk == 0) {
          continue;
        }
        double InvNormijk = basis_i.GetInvNorm()*basis_j.GetInvNorm()*basis_k.GetInvNorm();
        double tensorVal = Cijk*InvNormijk;
        if (tensorVal != 0) {
          tripletList.push_back(T(i,j, tensorVal));
        }
      }
    }
    (*Adv_tensor)[k].setFromTriplets(tripletList.begin(), tripletList.end());
  }
  // Make makeCompressed.
  for (int k = 0; k < basis_dim_; k++) {    
    (*Adv_tensor)[k].makeCompressed();
  }
}

double TwoNeumannX2D::ComputeZIntegral(const TwoNeumannBasis2D& basis_i,
                                       const TwoNeumannBasis2D& basis_j, const TwoNeumannBasis2D& basis_k) {
  const VEC2I I_ = basis_i.GetWaveNum();
  const VEC2I J_ = basis_j.GetWaveNum();
  const VEC2I K_ = basis_k.GetWaveNum();
  
  const double Ai = basis_i.GetA();
  const double Bi = basis_i.GetB();
  
  const double AkBj = basis_k.GetA()*basis_j.GetB();
  const double AjBk = basis_j.GetA()*basis_k.GetB();
  
  // (B_i i_1 + A_i i_2)
  const double CoefL = Bi*I_[0] + Ai*I_[1];
  if (CoefL == 0) return 0;
  // CoefR 
  double CoefR = 0;
  //++
  CoefR += (AkBj - AjBk)*ComputeIntegralCSSS(I_[0], I_[1], J_[0] + K_[0], J_[1] + K_[1]);
  //+-
  CoefR += (AkBj + AjBk)*ComputeIntegralCSSS(I_[0], I_[1], J_[0] + K_[0], J_[1] - K_[1]);
  //-+
  CoefR += (AkBj + AjBk)*ComputeIntegralCSSS(I_[0], I_[1], J_[0] - K_[0], J_[1] + K_[1]);
  //--
  CoefR += (AkBj - AjBk)*ComputeIntegralCSSS(I_[0], I_[1], J_[0] - K_[0], J_[1] - K_[1]);
  
  return CoefL*CoefR;
}

void TwoNeumannX2D::ComputeEigenValues(std::vector<double>* eigs, std::vector<double>* eigs_inv,
                          std::vector<double>* eigs_inv_root) {
  for (int i = 0; i < basis_dim_; i++) {
    double k1 = LookUpBasis(i,0) - 1;
    double k2 = LookUpBasis(i,1) - 1;
    (*eigs)[i] = (k1*k1 + k2*k2);
    if (k1 == 0 && k2 == 0) {
      (*eigs_inv)[i] = 0.;
      (*eigs_inv_root)[i] = 0.;
    }
    (*eigs_inv)[i] = 1.0/(k1*k1+k2*k2);
    (*eigs_inv_root)[i] = 1.0/sqrt(k1*k1+k2*k2);
  }
}  

void TwoNeumannX2D::FillAdvectionTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) {
  std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: " <<  "depredated." << std::endl; exit(0);
}

void TwoNeumannX2D::EnforceAntisymmetric(std::vector<Adv_Tensor_Type>* Adv_tensor) {
  std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: " <<  "Not implemented" << std::endl; exit(0);
}

void TwoNeumannX2D::ProjectVelocityNumerical(const VFIELD2D& field,
                                           Eigen::VectorXd* coefficients) {
   std::cout << "two_neumann_x_2D.cpp " << __LINE__ << " FATAL: " <<  "Not implemented" << std::endl; exit(0);
}
