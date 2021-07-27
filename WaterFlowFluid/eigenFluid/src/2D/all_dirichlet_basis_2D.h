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

#ifndef ALL_DIRICHLET_BASIS_2D_H
#define ALL_DIRICHLET_BASIS_2D_H

#include <cstring>
#include "Eigen"
#include <vector>

#include "2D/laplacian_basis_2D.h"
#include "2D/VFIELD2D.h"
#include "../setting.h"

class AllDirichletBasis2D : public LaplacianBasis2D {
  
public:
  AllDirichletBasis2D(int N, int root_num_basis): 
  LaplacianBasis2D(N, root_num_basis), num_basis_(root_num_basis*root_num_basis) {
    temp_ = (Real*) malloc(sizeof(Real)*N_*N_);
    std::memset(temp_, 0x00, sizeof(Real)*N_*N_);
  //  num_basis_square_ = num_basis_ * num_basis_;
  //  num_basis_cube_ = num_basis_ * num_basis_ * num_basis_;
  };
  ~AllDirichletBasis2D(){};
  void FillBasisNumerical(
      const Eigen::VectorXd& coefficients, VFIELD2D* field) override;
  double ComputeBasisAt(const int x_pos, const int y_pos,
                        const int band, const int mode) const override;
  void ProjectVelocityNumerical(const VFIELD2D& field,
      Eigen::VectorXd* coefficients) override;
  void InverseTramsformToVelocity(
      const Eigen::VectorXd& coefficients, VFIELD2D* field) override;
  // Input a vector field, transform and output the basis coefficients.
  void ForwardTransformtoFrequency(
      const VFIELD2D& field, Eigen::VectorXd* coefficients) override;
  // Fill the advection tensor for all Dirichlet case.
  void FillAdvectionTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) override;
  void ClearTemp() {
    std::memset(temp_, 0x00, sizeof(Real) * N_ * N_);
  }
  double compute_coefficient_terms(int i1, int j1, int i2, int j2, int mode);
  // void VerifyAntisymmetric(const std::vector<Adv_Tensor_Type>& C) override;
  void ComputeEigenValues(std::vector<double>* eigs, std::vector<double>* eigs_inv,
                          std::vector<double>* eigs_inv_root) override;
  void EnforceAntisymmetric(std::vector<Adv_Tensor_Type>* Adv_tensor) override;
 
  bool WithinRange(const int i1, const int i2);

  void FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) override;
  double GenerateVariationalTerms(const int j1, const int j2, const int k1, const int k2,
                                  const int index);

  double GenerateVelocityTensorTerms(const int i1, const int i2, const int j1,
                                     const int j2, const int index);
protected:
  const int num_basis_;
//  int64_t num_basis_square_;
//  int64_t num_basis_cube_;
  Real *temp_;
};

#endif  // ALL_DIRICHLET_BASIS_2D_H
