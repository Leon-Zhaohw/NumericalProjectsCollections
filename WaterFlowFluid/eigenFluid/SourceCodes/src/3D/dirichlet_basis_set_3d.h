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

#ifndef DIRICHLET_BASIS_SET_3D_H
#define DIRICHLET_BASIS_SET_3D_H

#include <cstdlib>
#include "Eigen"
#include <fstream>
#include <vector>

#include "3D/laplacian_basis_set_3d.h"
#include "3D/dirichlet_basis_3d.h"
#include "3D/VECTOR3_FIELD_3D.h"

// #define TRANSFORM_FAST

class DirichletBasisSet3D : public LaplacianBasisSet3D {
public:
  DirichletBasisSet3D(const int des_basis_dim, const int xRes, 
                      const int yRes, const int zRes, const std::string
                      const_strategy):
                      LaplacianBasisSet3D(des_basis_dim, xRes, yRes, zRes,
                                          const_strategy){
                            
    temp_ = (Real*) malloc(sizeof(Real)*xRes_*yRes_*zRes_);
#ifdef TRANSFORM_FAST
    tempZ_ = (Real*) malloc(sizeof(Real)*zRes_*Ny_*Nx_);
    tempY_ = (Real*) malloc(sizeof(Real)*zRes_*yRes_*Nx_);
#endif
    SlabSize_ = xRes_*yRes_;
    TotalSize_ = xRes_*yRes_*zRes_;
    invTotalSize_ = 1.0 / TotalSize_;
    basis_allocated_ = false;
    ClearTemp();
  }
  ~DirichletBasisSet3D(){
    free(temp_);
#ifdef TRANSFORM_FAST
    free(tempZ_);
    free(tempY_);
#endif
  }
  void FillBasisNumerical(
      const Eigen::VectorXd& coefficients, VECTOR3_FIELD_3D* field) override;
  // Allocate all basis, initialize all wavenumbers and A_,B_,C_,
  // Default basis allocation strategy.
  int AllocateAllBasis() override;
  // initialize the basis vector using externally allocated basis.
  void AllocateAllBasis(const std::vector<DirichletBasis3D>& basis);
  void ComputeEigenValues(Eigen::VectorXd* eigenValue) override;
  
  void InverseTramsformToVelocity(
      const Eigen::VectorXd& coefficients, VECTOR3_FIELD_3D* field) override;
  // Input a vector field, transform and output the basis coefficients.
  void ForwardTransformtoFrequency(
      const VECTOR3_FIELD_3D& field, Eigen::VectorXd* coefficients) override;
  void FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) override;
  // P_x in doc.
  double ComputeXIntegral(const DirichletBasis3D& basis_i, const DirichletBasis3D& basis_j,
                          const DirichletBasis3D& basis_k);
  // P_y in doc.
  double ComputeYIntegral(const DirichletBasis3D& basis_i, const DirichletBasis3D& basis_j,
                          const DirichletBasis3D& basis_k);
  // P_z in doc.
  double ComputeZIntegral(const DirichletBasis3D& basis_i, const DirichletBasis3D& basis_j,
                          const DirichletBasis3D& basis_k);
  
  void WriteBasis(std::ofstream& out) override;
  int ReadBasis(std::ifstream& infile) override;
  
  Eigen::Vector3d GenerateCoefficients(const int k1, const int k2, const int k3);
  double ComputeBasisAt(const VEC3I& pos, const int basis_idx, const int mode) const override;
  void PrintDebugInfo(const Eigen::VectorXd& coefficients) override {};
  void InverseTramsformToVelocityFast(
      const Eigen::VectorXd& coefficients, VECTOR3_FIELD_3D* field);
  void ForwardTransformtoFrequencyFast(
      const VECTOR3_FIELD_3D& field, Eigen::VectorXd* coefficients);
  void ComputeBasisWeights(Eigen::VectorXd* weights, const Real weightMultiplierC) override;
  
protected:
  // The order stored in this vector must corresponds to the order of the
  // coefficients stored in LaplacianFluid3D.
  std::vector<DirichletBasis3D> all_basis_;
  void ClearTemp() {
    memset(temp_, 0x00, sizeof(Real)*xRes_*yRes_*zRes_);
#ifdef TRANSFORM_FAST
    memset(tempZ_, 0x00, sizeof(Real)*zRes_*Ny_*Nx_);
    memset(tempY_, 0x00, sizeof(Real)*zRes_*yRes_*Nx_);
#endif
  }
    
  const int Nx_ = 20;
  const int Ny_ = 20;
  const int Nz_ = 20;
  Real * tempZ_;
  Real * tempY_;
  
  int SlabSize_;
  uint TotalSize_;
  Real* temp_;
  double invTotalSize_;
  bool basis_allocated_;
};

#endif  // DIRICHLET_BASIS_SET_3D_H
