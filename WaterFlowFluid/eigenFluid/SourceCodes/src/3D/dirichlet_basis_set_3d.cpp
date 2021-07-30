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

#include <algorithm>
#include "Eigen"
#include <fftw3.h>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <vector>

#include "3D/dirichlet_basis_set_3d.h"
#include "3D/dirichlet_basis_3d.h"
#include "3D/trig_integral_3d.h"
#include "3D/VECTOR3_FIELD_3D.h"
#include "Alg/VEC3.h"
#include "Alg/VEC2.h"
#include "util/timer.h"
#include "util/block_3d_dct.h"
#include "util/util.h"

void DirichletBasisSet3D::FillBasisNumerical(
      const Eigen::VectorXd& coefficients, VECTOR3_FIELD_3D* field) {
  if (coefficients.size() != basis_dim_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << " The number of coefficients" << " does not equal to number of basis." << std::endl; exit(0);}
  if (! basis_allocated_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
#pragma omp parallel for
  for (int i = 0; i < basis_dim_; i++) {
      if (i % 100 == 0) {
          std::cout <<  i << std::endl;
      }
    all_basis_[i].DiscretizeAdd(coefficients[i], field);
  }
}

void DirichletBasisSet3D::ComputeBasisWeights(Eigen::VectorXd* weights, const Real weightMultiplierC) {
 
  if (weights->size() != basis_dim_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  weights->setZero();
  
  std::vector<VEC2> WavenumAmp;
  
  for (int i = 0; i < basis_dim_; i++) {
    const VEC3I I = all_basis_[i].GetWaveNum();
    const int wsum = I[0]*I[0] + I[1]*I[1] + I[2]*I[2];
    (*weights)[i] = 1.0  + weightMultiplierC*(I[0]*I[0] + I[1]*I[1] + I[2]*I[2]);
     WavenumAmp.push_back(VEC2(wsum, (*weights)[i]));
  }
  std::sort(WavenumAmp.begin(), WavenumAmp.end(), [](VEC2 a, VEC2 b) {
    return a[0] < b [0];   
  });
  
  std::fstream out("weight.txt", std::ios::out);
  for (int i = 1; i < basis_dim_; i++) {
    out << WavenumAmp[i][0] << " " << WavenumAmp[i][1] << "\n";
  }
  out.close();
}

double DirichletBasisSet3D::ComputeBasisAt(const VEC3I& pos, 
                                           const int basis_idx, const int mode) const {
  if (basis_idx < 0 || basis_idx >= basis_dim_) {
    std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: " <<  "basis index is out of range: " << basis_idx << std::endl; exit(0);
  }
  return all_basis_[basis_idx].ComputeBasisAt(pos[0], pos[1], pos[2], mode, dxyz_);
}

void DirichletBasisSet3D::AllocateAllBasis(const std::vector<DirichletBasis3D>& basis) {
  for (int i = 0; i < basis.size(); i++) {
    all_basis_.push_back(basis[i]);
  }
  basis_dim_ = all_basis_.size();
  basis_allocated_ = true;
}

// Default basis allocation stragegy.

int DirichletBasisSet3D::AllocateAllBasis() {
  // Allocate more basis along longer dimension.
  const float Ydx = static_cast<float> (yRes_) / xRes_;
  const float Zdx = static_cast<float> (zRes_) / xRes_;
  int allow_freq = des_basis_dim_;
  
  float n_ = pow(static_cast<float>(allow_freq)/(Ydx*Zdx), 1.0 / 3.0);
  const int nx = static_cast<int>(n_);
  const int ny = static_cast<int>(n_*Ydx);
  const int nz = static_cast<int>(n_*Zdx);
  std::cout <<  "xRes: " << xRes_ << " Allocated basis: " << nx << std::endl;
  std::cout <<  "yRes: " << yRes_ << " Allocated basis: " << ny << std::endl;
  std::cout <<  "zRes: " << zRes_ << " Allocated basis: " << nz << std::endl;
  
  std::vector<VEC3I> waveNs;
  for (int k3 = 0; k3 < nz; k3 ++) {
    for (int k2 = 0; k2 < ny; k2 ++) {
      for (int k1 = 0; k1 < nx; k1 ++) {
        VEC3I V(k1,k2,k3);
        waveNs.push_back(V);
      }
    }
  }

  for (int i = 0; i < waveNs.size(); i ++) {
    const VEC3I & K_ = waveNs[i];
     
    Eigen::Vector3d coef = GenerateCoefficients(K_[0], K_[1], K_[2]);
    if (coef.norm() == 0) {
      continue;
    }
    DirichletBasis3D basis(K_[0], K_[1], K_[2], coef(0), coef(1), coef(2));
    all_basis_.push_back(basis);
  }
  
  basis_dim_ = all_basis_.size();
  basis_allocated_ = true;
  
  return all_basis_.size();
}

Eigen::Vector3d DirichletBasisSet3D::GenerateCoefficients
    (const int k1, const int k2, const int k3) {
  double a = 0, b = 0, c = 0;
  int num_zero = 0;
  if (k1 == 0) num_zero ++;
  if (k2 == 0) num_zero ++;
  if (k3 == 0) num_zero ++;
  if (num_zero == 3 || num_zero == 2) {
    return Eigen::Vector3d(0,0,0);
  }
  if (num_zero == 1) {
    if (k1 == 0) {
      return Eigen::Vector3d(0, k3, -k2);
      // return Eigen::Vector3d(0,0,0);
    }
    if (k2 == 0) {
      return Eigen::Vector3d(-k3, 0, k1);
      //return Eigen::Vector3d(0,0,0);
    }
    if (k3 == 0) {
      return Eigen::Vector3d(k2, -k1, 0);
      //return Eigen::Vector3d(0,0,0);
    }
  }
  if (num_zero == 0) {
    switch (int_type_const_strategy_) {
      case 0:  // principle_x
        return Eigen::Vector3d(-(k2*k2+k3*k3), k1*k2, k1*k3);
        break;
      case 1:  // principle_y
        return Eigen::Vector3d(k1*k2, -(k1*k1 + k3*k3), k2*k3);
        break;
      case 2: // principle_z
        return Eigen::Vector3d(k1*k3, k2*k3, -(k1*k1+k2*k2));
        break;
      case 3: // random
          {
            double r1 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) + 1.0;
            double r2 = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) + 1.0;
      
            const int which = rand() % 3;
            if (which == 0) { // rand on y,z
              const double A = (-r1*k2 - r2*k3) / static_cast<double>(k1);
              return Eigen::Vector3d(A, r1, r2);
            } else if (which == 1) { // rand on x, z
              const double B = (-r1*k1 - r2*k3) / static_cast<double>(k2);
              return Eigen::Vector3d(r1, B, r2);
            } else if (which == 2) {  // rand on x,y
              const double C = (-r1*k1 - r2*k2) / static_cast<double>(k3);
              return Eigen::Vector3d(r1, r2, C);
            } else {
                std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: " <<  "WTF.." << std::endl; exit(0);
            }
          }
        break;
      case 4:   // uniform.
          {
            const double A = -(k2*k2+k3*k3) + k1*k2 + k1*k3;
            const double B = k1*k2 -(k1*k1 + k3*k3) + k2*k3;
            const double C = k1*k3 +  k2*k3 -(k1*k1+k2*k2);
            return Eigen::Vector3d(A, B, C);
          }
      break;
      default:
        std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: " <<  "...." << std::endl; exit(0);
    }
  }
}

int DirichletBasisSet3D::ReadBasis(std::ifstream& infile) {
  
  if (! !basis_allocated_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << "Basis already allocated." << std::endl; exit(0);}
  
  if (!infile.is_open()) {
    std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: " <<  "Cannot read from stream." << std::endl; exit(0);
  }
  int read_type = -1;
  int read_dim = -1;
 
  // Read in the basis dimension and type.
  infile.read(reinterpret_cast<char*>(&read_dim), sizeof(int));
  infile.read(reinterpret_cast<char*>(&read_type), sizeof(int));
  
  if (read_type != 0) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << "expect dirichlet type." << std::endl; exit(0);}
  // Read and allocate the basis.
  for (int i = 0; i < read_dim; i++) {
    VEC3I K_;
    double A,B,C;
    infile.read(reinterpret_cast<char*>(&K_[0]), sizeof(int));
    infile.read(reinterpret_cast<char*>(&K_[1]), sizeof(int));
    infile.read(reinterpret_cast<char*>(&K_[2]), sizeof(int));
    infile.read(reinterpret_cast<char*>(&A), sizeof(double));
    infile.read(reinterpret_cast<char*>(&B), sizeof(double));
    infile.read(reinterpret_cast<char*>(&C), sizeof(double));
    DirichletBasis3D basis (K_[0], K_[1], K_[2], A, B, C);
    all_basis_.push_back(basis);
  }
  
  basis_dim_ = all_basis_.size();
  basis_allocated_ = true;
  return all_basis_.size();
}

void DirichletBasisSet3D::WriteBasis(std::ofstream& out) {
  if (!out.is_open()) {
    std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: " <<  "Cannot write to file stream" << std::endl; exit(0);
  }
  int basis_type_ = 0;
  if (basis_dim_ != all_basis_.size()) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: " << std::endl; exit(0);}
  // first write out the basis dimension.
  out.write(reinterpret_cast<const char *>(&basis_dim_), sizeof(int));
  // Then write out the tensor type. encoded as int.
  out.write(reinterpret_cast<const char *>(&basis_type_), sizeof(int));
  // Write the wavenumber and coefficients for each basis.
  for (int i = 0; i < basis_dim_; i++) {
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double A_ = all_basis_[i].GetA();
    double B_ = all_basis_[i].GetB();
    double C_ = all_basis_[i].GetC();
    out.write(reinterpret_cast<const char *>(&K_[0]), sizeof(int));
    out.write(reinterpret_cast<const char *>(&K_[1]), sizeof(int));
    out.write(reinterpret_cast<const char *>(&K_[2]), sizeof(int));
    out.write(reinterpret_cast<const char *>(&A_), sizeof(double));
    out.write(reinterpret_cast<const char *>(&B_), sizeof(double));
    out.write(reinterpret_cast<const char *>(&C_), sizeof(double));
  }
}

void DirichletBasisSet3D::ComputeEigenValues(Eigen::VectorXd* eigenValue) {
  if (! basis_allocated_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  for (int i = 0; i < basis_dim_; i++) {
    VEC3I K_ = all_basis_[i].GetWaveNum();
    (*eigenValue)(i) = K_[0]*K_[0] + K_[1]*K_[1] + K_[2]*K_[2];
  }
}

void DirichletBasisSet3D::FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) {
  if (! basis_allocated_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  if (Adv_tensor != nullptr ){ Adv_tensor->clear();}
  Adv_tensor->reserve(basis_dim_);
  // Fill the tensor with zeros.
  for (int k = 0; k < basis_dim_; k++) {
    Adv_Tensor_Type Ck(basis_dim_, basis_dim_);
    Ck.setZero();
    Adv_tensor->emplace_back(Ck);
  }
  const int five_percent = basis_dim_ / 20;
#pragma omp parallel for
  for (int k = 0; k < basis_dim_; k++) {
    if (k % five_percent == 0) {
      std::cout <<  "% 5 " << "Tensor computed." << std::endl;
    }
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    
    for (int j = 0; j < basis_dim_; j++) {
      for (int i = 0; i < basis_dim_; i++) {
        //
        const DirichletBasis3D& basis_i = all_basis_[i];
        const DirichletBasis3D& basis_j = all_basis_[j];
        const DirichletBasis3D& basis_k = all_basis_[k];
        // Cijk = 0.125*(P_x + P_y + P_z) * invNorm_i *invNorm_j * invNorm_k
        const double Px = ComputeXIntegral(basis_i, basis_j, basis_k);
        const double Py = ComputeYIntegral(basis_i, basis_j, basis_k);
        const double Pz = ComputeZIntegral(basis_i, basis_j, basis_k);
        
        double Cijk = 0.125*(Px + Py + Pz);
        if (Cijk == 0) {
          continue;
        }
        double InvNormijk = basis_i.GetInvNorm()*basis_j.GetInvNorm()*basis_k.GetInvNorm();
        // IncrementMatrixEntry(&(*Adv_tensor)[k],i, j, Cijk*InvNormijk);
        tripletList.push_back(T(i,j, Cijk*InvNormijk));
      }
    }
     (*Adv_tensor)[k].setFromTriplets(tripletList.begin(), tripletList.end());
  }
  // Make makeCompressed.
  for (int k = 0; k < basis_dim_; k++) {    
    (*Adv_tensor)[k].makeCompressed();
  }
}

// P_x
double DirichletBasisSet3D::ComputeXIntegral(const DirichletBasis3D& basis_i, const DirichletBasis3D& basis_j,
                                             const DirichletBasis3D& basis_k) {
  const VEC3I I_ = basis_i.GetWaveNum();
  const VEC3I J_ = basis_j.GetWaveNum();
  const VEC3I K_ = basis_k.GetWaveNum();
  
  const double Bi = basis_i.GetB();
  const double Ci = basis_i.GetC();

  const double BkCj = basis_k.GetB() * basis_j.GetC();
  const double BjCk = basis_j.GetB() * basis_k.GetC();
  // CoefL = Bi*i_3 - Ci*i_2
  const double CoefL = (Bi*I_[2] - Ci*I_[1]);
  if (CoefL == 0) return 0;
  
  // CoefR = \sum_{i = 0}^7 c_i^x \int_{\Omega} \xi^6(i_1, i_2, i_3) \xi^6(...) d\Omega
  double CoefR = 0;
  // +++
  CoefR +=  (BkCj - BjCk)*ComputeIntegralCSS(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] + K_[1], J_[2] + K_[2]);
  // ++-
  CoefR +=  (BkCj + BjCk)*ComputeIntegralCSS(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] + K_[1], J_[2] - K_[2]);
  // +-+
  CoefR += -(BkCj + BjCk)*ComputeIntegralCSS(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] - K_[1], J_[2] + K_[2]);
  // +--
  CoefR += (-BkCj + BjCk)*ComputeIntegralCSS(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] - K_[1], J_[2] - K_[2]);
  // -++
  CoefR +=  (BkCj - BjCk)*ComputeIntegralCSS(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] + K_[1], J_[2] + K_[2]);
  // -+-
  CoefR +=  (BkCj + BjCk)*ComputeIntegralCSS(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] + K_[1], J_[2] - K_[2]);
  // --+
  CoefR += -(BkCj + BjCk)*ComputeIntegralCSS(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] - K_[1], J_[2] + K_[2]);
  // ---
  CoefR += (-BkCj + BjCk)*ComputeIntegralCSS(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] - K_[1], J_[2] - K_[2]);
  
  return CoefL*CoefR;
}

// P_y
double DirichletBasisSet3D::ComputeYIntegral(const DirichletBasis3D& basis_i, const DirichletBasis3D& basis_j,
                                             const DirichletBasis3D& basis_k) {
  const VEC3I I_ = basis_i.GetWaveNum();
  const VEC3I J_ = basis_j.GetWaveNum();
  const VEC3I K_ = basis_k.GetWaveNum();
  
  const double Ci = basis_i.GetC();
  const double Ai = basis_i.GetA();
  
  const double AkCj = basis_k.GetA() * basis_j.GetC();
  const double AjCk = basis_j.GetA() * basis_k.GetC();
  // CoefL = Ci*i_1 - Ai*i_3
  const double CoefL = (Ci*I_[0] - Ai*I_[2]);
  if (CoefL == 0) return 0;
  // CoefR = \sum_{i = 0}^7 c_i^y \int_{\Omega} \xi^5(i_1, i_2, i_3) \xi^5(...) d\Omega
  double CoefR = 0;
  // +++
  CoefR += (-AkCj + AjCk)*ComputeIntegralSCS(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] + K_[1], J_[2] + K_[2]);
  // ++-
  CoefR += (-AkCj - AjCk)*ComputeIntegralSCS(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] + K_[1], J_[2] - K_[2]);
  // +-+
  CoefR += (-AkCj + AjCk)*ComputeIntegralSCS(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] - K_[1], J_[2] + K_[2]);
  // +--
  CoefR += (-AkCj - AjCk)*ComputeIntegralSCS(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] - K_[1], J_[2] - K_[2]);
  // -++
  CoefR +=  (AkCj + AjCk)*ComputeIntegralSCS(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] + K_[1], J_[2] + K_[2]);
  // -+-
  CoefR +=  (AkCj - AjCk)*ComputeIntegralSCS(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] + K_[1], J_[2] - K_[2]);
  // --+
  CoefR +=  (AkCj + AjCk)*ComputeIntegralSCS(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] - K_[1], J_[2] + K_[2]);
  // ---
  CoefR +=  (AkCj - AjCk)*ComputeIntegralSCS(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] - K_[1], J_[2] - K_[2]);
  
  return CoefL*CoefR;
}

// P_z
double DirichletBasisSet3D::ComputeZIntegral(const DirichletBasis3D& basis_i, const DirichletBasis3D& basis_j,
                                             const DirichletBasis3D& basis_k) {
  const VEC3I I_ = basis_i.GetWaveNum();
  const VEC3I J_ = basis_j.GetWaveNum();
  const VEC3I K_ = basis_k.GetWaveNum();
  
  const double Ai = basis_i.GetA();
  const double Bi = basis_i.GetB();
  
  const double AkBj = basis_k.GetA()*basis_j.GetB();
  const double AjBk = basis_j.GetA()*basis_k.GetB();
  // CoefL = Ai*i_2 - Bi*i_1
  const double CoefL = (Ai*I_[1] - Bi*I_[0]);
  if (CoefL == 0) return 0;
  // CoefR = \sum_{i = 0}^7 c_i^z \int_{\Omega} \xi^4(i_1, i_2, i_3) \xi^4(...) d\Omega
  double CoefR = 0;
  // +++
  CoefR +=  (AkBj - AjBk)*ComputeIntegralSSC(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] + K_[1], J_[2] + K_[2]);
  // ++-
  CoefR +=  (AkBj - AjBk)*ComputeIntegralSSC(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] + K_[1], J_[2] - K_[2]);
  // +-+
  CoefR +=  (AkBj + AjBk)*ComputeIntegralSSC(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] - K_[1], J_[2] + K_[2]);
  // +--
  CoefR +=  (AkBj + AjBk)*ComputeIntegralSSC(I_[0], I_[1], I_[2], J_[0] + K_[0], J_[1] - K_[1], J_[2] - K_[2]);
  // -++
  CoefR += (-AkBj - AjBk)*ComputeIntegralSSC(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] + K_[1], J_[2] + K_[2]);
  // -+-
  CoefR += (-AkBj - AjBk)*ComputeIntegralSSC(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] + K_[1], J_[2] - K_[2]);
  // --+
  CoefR += (-AkBj + AjBk)*ComputeIntegralSSC(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] - K_[1], J_[2] + K_[2]);
  // ---
  CoefR += (-AkBj + AjBk)*ComputeIntegralSSC(I_[0], I_[1], I_[2], J_[0] - K_[0], J_[1] - K_[1], J_[2] - K_[2]);
  
  return CoefL*CoefR;
}

// Input a vector field, transform and output the basis coefficients.
void DirichletBasisSet3D::ForwardTransformtoFrequency(
  const VECTOR3_FIELD_3D& field, Eigen::VectorXd* coefficients) {
#ifdef TRANSFORM_FAST
  ForwardTransformtoFrequencyFast(field, coefficients);
#else
  if (! basis_allocated_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  if (coefficients->size() != basis_dim_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  coefficients->setZero();
  ClearTemp();
  // vx...
  // Put the value into field.
  for (uint i = 0; i < TotalSize_; i++) {
     temp_[i] = field[i][0];
  }

  fftw_plan plan_x_3D;
  // x: RODFT10, y: REDFT10 z: REDFT10.
  plan_x_3D = fftw_plan_r2r_3d(zRes_, yRes_, xRes_, temp_, temp_,
                               FFTW_REDFT10 , FFTW_REDFT10, FFTW_RODFT10, FFTW_MEASURE);
  fftw_execute(plan_x_3D);
  fftw_destroy_plan(plan_x_3D);
  for (int i = 0; i < basis_dim_; i++) {
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_x = all_basis_[i].GetVxCoef();
    if (K_[0] < 1) {  // vx = 0., set coef as zero.
      continue;
    }
    (*coefficients)[i] += temp_[K_[0] -1 + K_[1]*xRes_ + K_[2]*SlabSize_]*invTotalSize_*PI_CUBE*0.125*C_x;
  }
  // vy...
  ClearTemp();
  // Put the value into field.
  for (uint i = 0; i < TotalSize_; i++) {
     temp_[i] = field[i][1];
  }
  fftw_plan plan_y_3D;
  // x: REDFT10, y: RODFT10 z: REDFT10.
  plan_y_3D = fftw_plan_r2r_3d(zRes_, yRes_, xRes_, temp_, temp_,
                               FFTW_REDFT10, FFTW_RODFT10, FFTW_REDFT10, FFTW_MEASURE);
  fftw_execute(plan_y_3D);
  fftw_destroy_plan(plan_y_3D);
  for (int i = 0; i < basis_dim_; i++) {
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_y = all_basis_[i].GetVyCoef();
    if (K_[1] < 1) {  // vy = 0., set coef as zero.
      continue;
    }
    (*coefficients)[i] += temp_[K_[0] + (K_[1] - 1)*xRes_ + K_[2]*SlabSize_]*invTotalSize_*PI_CUBE*0.125*C_y;
  }
  // vz...
  ClearTemp();
  // Put the value into field.
  for (uint i = 0; i < TotalSize_; i++) {
     temp_[i] = field[i][2];
  }
  fftw_plan plan_z_3D;
  // x: REDFT10, y: REDFT10 z: RODFT10.
  plan_z_3D = fftw_plan_r2r_3d(zRes_, yRes_, xRes_, temp_, temp_,
                               FFTW_RODFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_MEASURE);
  fftw_execute(plan_z_3D);
  fftw_destroy_plan(plan_z_3D);
  for (int i = 0; i < basis_dim_; i++) {
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_z = all_basis_[i].GetVzCoef();
    if (K_[2] < 1) {  // vz = 0., set coef as zero.
      continue;
    }
    (*coefficients)[i] += temp_[K_[0] + K_[1]*xRes_ + (K_[2] - 1)*SlabSize_]*invTotalSize_*PI_CUBE*0.125*C_z;
  }
#endif
}

void DirichletBasisSet3D::InverseTramsformToVelocity(
    const Eigen::VectorXd& coefficients, VECTOR3_FIELD_3D* field){
#ifdef TRANSFORM_FAST
  InverseTramsformToVelocityFast(coefficients, field);
#else
  if (! basis_allocated_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  if (coefficients.size() != all_basis_.size()) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: " << "The basis size and " << "coefficients size must be equal." << std::endl; exit(0);}
  field->clear();
  ClearTemp();
  // vx..
  // Iterate over all the coefficients and put it into the field.
  for (int i = 0; i < basis_dim_; i++) {
    if (coefficients[i] == 0) {
      continue;
    }
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_x = all_basis_[i].GetVxCoef()*all_basis_[i].GetDCTNorm();
    if (K_[0] < 1) {  // vx = 0.
      continue;
    }
    // Process x field. x: RODFT01, y: REDFT01 z: REDFT01. 
    temp_[K_[0] -1 + K_[1]*xRes_ + K_[2]*SlabSize_] = coefficients[i] * C_x;
  }
  fftw_plan plan_x_3D;
  // !!! DONT ASSUME ANYTHING, FFTW IS ROW MAJOR !!! 
  plan_x_3D = fftw_plan_r2r_3d(zRes_, yRes_, xRes_, temp_, temp_,
                                FFTW_REDFT01, FFTW_REDFT01, FFTW_RODFT01, FFTW_MEASURE);
  fftw_execute(plan_x_3D);  
  fftw_destroy_plan(plan_x_3D);
  // Put the value back into field.
  for (uint i = 0; i < TotalSize_; i++) {
     (*field)[i][0] = temp_[i];
  }
  // vy...
  ClearTemp();
  for (int i = 0; i < basis_dim_; i++) {
    if (coefficients[i] == 0) {
      continue;
    }
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_y = all_basis_[i].GetVyCoef()*all_basis_[i].GetDCTNorm();
    if (K_[1] < 1) {  // vy = 0.
      continue;
    }
    // Process y field. x: REDFT01, y: RODFT01 z: REDFT01.
    temp_[K_[0] + (K_[1] - 1)*xRes_ + K_[2]*SlabSize_] = coefficients[i] * C_y;
  }
  fftw_plan plan_y_3D;
  plan_y_3D = fftw_plan_r2r_3d(zRes_, yRes_, xRes_, temp_, temp_,
                               FFTW_REDFT01, FFTW_RODFT01, FFTW_REDFT01, FFTW_MEASURE);
  fftw_execute(plan_y_3D);  
  fftw_destroy_plan(plan_y_3D);
  // Put the value back into field.
  for (uint i = 0; i < TotalSize_; i++) {
     (*field)[i][1] = temp_[i];
  }
  // vz...
  ClearTemp();
  for (int i = 0; i < basis_dim_; i++) {
    if (coefficients[i] == 0) {
      continue;
    }
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_z = all_basis_[i].GetVzCoef()*all_basis_[i].GetDCTNorm();
    if (K_[2] < 1) {  // vz = 0.
      continue;
    }
    // Process z field. x: REDFT01, y: REDFT01 z: RODFT01.
    temp_[K_[0] + K_[1]*xRes_ + (K_[2] - 1)*SlabSize_] = coefficients[i] * C_z;
  }
  fftw_plan plan_z_3D;
  plan_z_3D = fftw_plan_r2r_3d(zRes_, yRes_, xRes_, temp_, temp_,
                               FFTW_RODFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_MEASURE);
  fftw_execute(plan_z_3D);
  fftw_destroy_plan(plan_z_3D);
  // Put the value back into field.
  for (uint i = 0; i < TotalSize_; i++) {
     (*field)[i][2] = temp_[i];
  }
#endif
}

void DirichletBasisSet3D::ForwardTransformtoFrequencyFast(
      const VECTOR3_FIELD_3D& field, Eigen::VectorXd* coefficients) {
   if (! basis_allocated_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  if (coefficients->size() != basis_dim_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  coefficients->setZero();
  ClearTemp();
  // vx...
  // Put the value into field.
  for (uint i = 0; i < TotalSize_; i++) {
     temp_[i] = field[i][0];
  }
   // x: RODFT10, y: REDFT10 z: REDFT10.
  block_DCT_V2F(xRes_, yRes_, zRes_, Nx_, Ny_, Nz_, temp_, tempY_, tempZ_,
               FFTW_RODFT10, FFTW_REDFT10, FFTW_REDFT10);
  const int zSlab = Nx_*zRes_;
  for (int i = 0; i < basis_dim_; i++) {
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_x = all_basis_[i].GetVxCoef();
    if (K_[0] < 1) {  // vx = 0., set coef as zero.
      continue;
    }
    (*coefficients)[i] += tempZ_[(K_[0] -1)*zRes_ + K_[1]*zSlab + K_[2]]*invTotalSize_*PI_CUBE*0.125*C_x;
  }
  
  // vy...
  ClearTemp();
  // Put the value into field.
  for (uint i = 0; i < TotalSize_; i++) {
     temp_[i] = field[i][1];
  }
  // x: REDFT10, y: RODFT10 z: REDFT10.
  block_DCT_V2F(xRes_, yRes_, zRes_, Nx_, Ny_, Nz_, temp_, tempY_, tempZ_,
               FFTW_REDFT10, FFTW_RODFT10, FFTW_REDFT10);
  
  for (int i = 0; i < basis_dim_; i++) {
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_y = all_basis_[i].GetVyCoef();
    if (K_[1] < 1) {  // vy = 0., set coef as zero.
      continue;
    }
    (*coefficients)[i] += tempZ_[K_[0]*zRes_ + (K_[1] - 1)*zSlab + K_[2]]*invTotalSize_*PI_CUBE*0.125*C_y;
  }
  // vz...
  ClearTemp();
  // Put the value into field.
  for (uint i = 0; i < TotalSize_; i++) {
     temp_[i] = field[i][2];
  }
  // x: REDFT10, y: REDFT10 z: RODFT10.
  block_DCT_V2F(xRes_, yRes_, zRes_, Nx_, Ny_, Nz_, temp_, tempY_, tempZ_,
                FFTW_REDFT10, FFTW_REDFT10, FFTW_RODFT10);
  
  for (int i = 0; i < basis_dim_; i++) {
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_z = all_basis_[i].GetVzCoef();
    if (K_[2] < 1) {  // vz = 0., set coef as zero.
      continue;
    }
    (*coefficients)[i] += tempZ_[K_[0]*zRes_ + K_[1]*zSlab + (K_[2] - 1)]*invTotalSize_*PI_CUBE*0.125*C_z;
  }
}

void DirichletBasisSet3D::InverseTramsformToVelocityFast(
      const Eigen::VectorXd& coefficients, VECTOR3_FIELD_3D* field) {

  if (! basis_allocated_) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: "  << "Basis not allocated." << std::endl; exit(0);}
  if (coefficients.size() != all_basis_.size()) {std::cout << "dirichlet_basis_set_3d.cpp " << __LINE__ << " FATAL: " << "The basis size and " << "coefficients size must be equal." << std::endl; exit(0);}
  field->clear();
  ClearTemp();
  const int yslabSize = zRes_*Nx_; 
  // vx..
  // Iterate over all the coefficients and put it into the field.
  for (int i = 0; i < basis_dim_; i++) {
    if (coefficients[i] == 0) {
      continue;
    }
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_x = all_basis_[i].GetVxCoef()*all_basis_[i].GetDCTNorm();
    if (K_[0] < 1) {  // vx = 0.
      continue;
    }
    // Process x field. x: RODFT01, y: REDFT01 z: REDFT01. 
    // tempX_[K_[0] -1 + K_[1]*xRes_ + K_[2]*xslabSize] = coefficients[i] * C_x;
    tempZ_[K_[2] + (K_[0] -1)*zRes_ + K_[1]*yslabSize] = coefficients[i] * C_x;
  }

  block_DCT_F2V(xRes_, yRes_, zRes_, Nx_, Ny_, Nz_, tempZ_, tempY_, temp_,
    FFTW_RODFT01, FFTW_REDFT01, FFTW_REDFT01);
  
  for (uint i = 0; i < TotalSize_; i++) {
     (*field)[i][0] = temp_[i];
  }
  
// vy...
  ClearTemp();
  for (int i = 0; i < basis_dim_; i++) {
    if (coefficients[i] == 0) {
      continue;
    }
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_y = all_basis_[i].GetVyCoef()*all_basis_[i].GetDCTNorm();
    if (K_[1] < 1) {  // vy = 0.
      continue;
    }
    // Process y field. x: REDFT01, y: RODFT01 z: REDFT01.
    tempZ_[K_[0]*zRes_ + (K_[1] - 1)*yslabSize + K_[2]] = coefficients[i] * C_y;
  }
  
  block_DCT_F2V(xRes_, yRes_, zRes_, Nx_, Ny_, Nz_, tempZ_, tempY_, temp_,
  FFTW_REDFT01, FFTW_RODFT01, FFTW_REDFT01);
  
  // Put the value back into field.
  for (uint i = 0; i < TotalSize_; i++) {
     (*field)[i][1] = temp_[i];
  }
  // vz...
  ClearTemp();
  for (int i = 0; i < basis_dim_; i++) {
    if (coefficients[i] == 0) {
      continue;
    }
    // Get wavenumber.
    VEC3I K_ = all_basis_[i].GetWaveNum();
    double C_z = all_basis_[i].GetVzCoef()*all_basis_[i].GetDCTNorm();
    if (K_[2] < 1) {  // vz = 0.
      continue;
    }
    // Process z field. x: REDFT01, y: REDFT01 z: RODFT01.
    tempZ_[K_[0]*zRes_ + K_[1]*yslabSize + (K_[2] - 1)] = coefficients[i] * C_z;
  }
  
  block_DCT_F2V(xRes_, yRes_, zRes_, Nx_, Ny_, Nz_, tempZ_, tempY_, temp_,
  FFTW_REDFT01, FFTW_REDFT01, FFTW_RODFT01);
  
  // Put the value back into field.
  for (uint i = 0; i < TotalSize_; i++) {
     (*field)[i][2] = temp_[i];
  }
}