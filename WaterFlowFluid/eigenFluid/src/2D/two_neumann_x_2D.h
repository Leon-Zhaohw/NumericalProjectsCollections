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

#ifndef TWO_NEUMANN_X_2D_H
#define TWO_NEUMANN_X_2D_H

#include "Eigen"
// #include <glog/logging.h>
#include <math.h>
#include <vector>

#include "2D/laplacian_basis_2D.h"
#include "2D/two_neumann_basis_2D.h"
#include "2D/VFIELD2D.h"

/* New implementation of TwoNeumann basis along x direction.
 * The implementation in the two_newmann_x.h is old one.*/

class TwoNeumannX2D : public LaplacianBasis2D {
public:
  TwoNeumannX2D(int N, int root_basis_dim): 
  LaplacianBasis2D(N, root_basis_dim), basis_dim_(root_basis_dim*root_basis_dim),
  dxy_(VEC2(M_PI / static_cast<double>(N), M_PI / static_cast<double>(N))),
  root_basis_dim_(root_basis_dim)
  {
    temp_ = (Real*) malloc(sizeof(Real)*N_*N_);
    std::memset(temp_, 0x00, sizeof(Real)*N_*N_);
    inv_sqrt_two_ = 1.0 / sqrt(2.0);
    basis_allocated_ = false;
    if (AllocateAllBasis() != basis_dim_) {std::cout << "two_neumann_x_2D.h " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  };
  ~TwoNeumannX2D(){};
  void FillBasisNumerical(
      const Eigen::VectorXd& coefficients, VFIELD2D* field) override;
      
  double ComputeBasisAt(const int x_pos, const int y_pos, const int band, const int mode) const;
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
  
  // void VerifyAntisymmetric(const std::vector<Adv_Tensor_Type>& C) override;
  void ComputeEigenValues(std::vector<double>* eigs, std::vector<double>* eigs_inv,
                          std::vector<double>* eigs_inv_root) override;
                          
  void EnforceAntisymmetric(std::vector<Adv_Tensor_Type>* Adv_tensor) override;
  
  
  void FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) override;

  double ComputeZIntegral(const TwoNeumannBasis2D& basis_i, const TwoNeumannBasis2D& basis_j, const TwoNeumannBasis2D& basis_k);
  // Default basis allocation strategy.
  int AllocateAllBasis();
  void ComputeVorticity(const Eigen::VectorXd& coef, FIELD2D* vort) override;
  void ComputeBasisWeights(Eigen::VectorXd* weights) override;
  
protected:
  std::vector<TwoNeumannBasis2D> all_basis_;
  const int basis_dim_;
  const int root_basis_dim_;
  Real *temp_;
  double inv_sqrt_two_;
  bool basis_allocated_;
  const VEC2 dxy_;
};

#endif  // TWO_NEUMANN_X_2D_H
