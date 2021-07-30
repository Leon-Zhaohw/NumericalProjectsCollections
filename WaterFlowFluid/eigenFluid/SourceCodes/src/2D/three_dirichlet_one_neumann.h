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

#ifndef THREE_DIRICHLET_ONE_NEUMANN_H
#define THREE_DIRICHLET_ONE_NEUMANN_H

#include "Eigen"
#include "Sparse"
#include <vector>
#include "../setting.h"
#include "2D/laplacian_basis_2D.h"
#include "2D/VFIELD2D.h"
class ThreeDirichletOneNeumann : public LaplacianBasis2D{

public:
  ThreeDirichletOneNeumann(int N, int root_num_basis): 
  LaplacianBasis2D(N, root_num_basis), num_basis_(root_num_basis * root_num_basis){
    temp_ = (Real*) malloc(sizeof(Real)*2*N_*N_);
    std::memset(temp_, 0x00, sizeof(Real)*2*N_*N_);
  };
  ~ThreeDirichletOneNeumann(){};
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
  void FillAdvectionTensor(
      std::vector<Adv_Tensor_Type> *Adv_tensor) override;
  double compute_coefficient_terms(int a, int b, int c, double i1j2, double i2j1,
                                   double inv_lambda_i, int mode);
  void ClearTemp() {
    std::memset(temp_, 0x00, sizeof(Real) * 2 * N_ * N_);
  }
  // void VerifyAntisymmetric(const std::vector<Adv_Tensor_Type>& C) override;
  void ComputeEigenValues(
  std::vector<double>* eigs, std::vector<double>* eigs_inv, 
  std::vector<double>* eigs_inv_root) override;
  void EnforceAntisymmetric(std::vector<Adv_Tensor_Type>* Adv_tensor) override;


  
//  void FillVelocityTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) override;
//  void FillThirdOrderTensor(CTF::Tensor<double> *C) override;
//  void FillFourthOrderTensor(CTF::Tensor<double> *Q) override;
  void FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) override;
  double ComputeVarVariationalTerms(const double j1k2, const double j2k1, const int i1,
                                  const int j1, const int k1, const int index);
protected:
  Real *temp_;
  int num_basis_;
};

#endif  // THREE_DIRICHLET_ONE_NEUMANN_H
