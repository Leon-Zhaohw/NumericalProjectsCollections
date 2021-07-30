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

#ifndef ONE_NEUMANN_BASIS_SET_3D_H
#define ONE_NEUMANN_BASIS_SET_3D_H

#include "Eigen"
#include <fstream>
#include <vector>

#include "3D/laplacian_basis_set_3d.h"
#include "3D/one_neumann_basis_3d.h"

class OneNeumannBasisSet3D: public LaplacianBasisSet3D {
public:
  OneNeumannBasisSet3D(const int basis_dim, const int xRes, 
                      const int yRes, const int zRes, const std::string
                      const_strategy):
                      LaplacianBasisSet3D(basis_dim, xRes, yRes, zRes, 
                                          const_strategy){
    // This one need twice memory along x ...
    temp_ = (Real*) malloc(sizeof(Real)*xRes_*yRes_*zRes_*2);
    memset(temp_, 0x00, sizeof(Real)*xRes_*yRes_*zRes_*2);
    SlabSize_ = xRes_*yRes_;
    TotalSize_ = xRes_*yRes_*zRes_;
    invTotalSize_ = 1.0 / TotalSize_;
    basis_allocated_ = false;
  }
  ~OneNeumannBasisSet3D(){
    free(temp_);
  }
  void FillBasisNumerical(
      const Eigen::VectorXd& coefficients, VECTOR3_FIELD_3D* field) override;
// Allocate all basis, initialize all wavenumbers and A_,B_,C_,
  // Default basis allocation strategy.
  int AllocateAllBasis() override;
  // initialize the basis vector using externally allocated basis.
  void AllocateAllBasis(const std::vector<OneNeumannBasis3D>& basis);
  
  void InverseTramsformToVelocity(
      const Eigen::VectorXd& coefficients, VECTOR3_FIELD_3D* field) override;
  // Input a vector field, transform and output the basis coefficients.
  void ForwardTransformtoFrequency(
      const VECTOR3_FIELD_3D& field, Eigen::VectorXd* coefficients) override;
  void FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) override;
  void ComputeEigenValues(Eigen::VectorXd* eigenValue) override;
  void WriteBasis(std::ofstream& out) override;
  int ReadBasis(std::ifstream& infile) override;
  Eigen::Vector3d GenerateCoefficients(const int k1, const int k2, const int k3);
  // P_x in doc.
  double ComputeXIntegral(const OneNeumannBasis3D& basis_i, const OneNeumannBasis3D& basis_j,
                          const OneNeumannBasis3D& basis_k);
  // P_y in doc.
  double ComputeYIntegral(const OneNeumannBasis3D& basis_i, const OneNeumannBasis3D& basis_j,
                          const OneNeumannBasis3D& basis_k);
  // P_z in doc.
  double ComputeZIntegral(const OneNeumannBasis3D& basis_i, const OneNeumannBasis3D& basis_j,
                          const OneNeumannBasis3D& basis_k);
  double ComputeBasisAt(const VEC3I& pos, const int basis_idx, const int mode) const override;
  void PrintDebugInfo(const Eigen::VectorXd& coefficients) override {};
  void ComputeBasisWeights(Eigen::VectorXd* weights, const Real weightMultiplierC) override;
protected:
  // The order stored in this vector must corresponds to the order of the
  // coefficients stored in LaplacianFluid3D.
  std::vector<OneNeumannBasis3D> all_basis_;
  void ClearTemp() {
    memset(temp_, 0x00, sizeof(Real)*xRes_*yRes_*zRes_*2);
  }
  int SlabSize_;
  uint TotalSize_;
  Real* temp_;
  double invTotalSize_;
  bool basis_allocated_;
};

#endif  // ONE_NEUMANN_BASIS_SET_3D_H
