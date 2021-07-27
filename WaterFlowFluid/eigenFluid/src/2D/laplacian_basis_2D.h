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

#ifndef LAPLACIAN_BASIS_2D_H
#define LAPLACIAN_BASIS_2D_H

#include <cstring>
#include "Eigen"
#include <stdlib.h>
#include <math.h>
#include <memory>
#include <unordered_map>
#include <vector>

#include "2D/basis_lookup_2D.h"
#include "setting.h"
#include "2D/VFIELD2D.h"

enum LaplacianBasisType {
  ALL_DIRICHLET = 0,
  THREE_DIRICHLET_ONE_NEUMANN = 1
};
// The base class for all kinds of laplacian basis with different 
// boundary conditions.
class LaplacianBasis2D {
public:
    LaplacianBasis2D(int N, int root_num_basis):N_(N), 
      root_num_basis_(root_num_basis) {
    basislookup_.reset(new BasisLookUp(root_num_basis));
    delta_ = 1.0 * M_PI / N_;
    inv_N_ = 1.0 / N_;
    inv_pi_2_ = 2.0 / M_PI;
    inv_pi_cube_8_ = 8.0 / (M_PI*M_PI*M_PI);
    inv_pi_square_ = 1.0 / (M_PI*M_PI); 
  };
  ~LaplacianBasis2D(){}
  std::unique_ptr<BasisLookUp> basislookup_;
  virtual void FillBasisNumerical(
      const Eigen::VectorXd& coefficients, VFIELD2D* field) = 0;
  // Compute the value of the basis at a particular position.
  virtual double ComputeBasisAt(const int x_pos, const int y_pos,
                               const int band, const int mode) const = 0;
  virtual void ProjectVelocityNumerical(const VFIELD2D& field,
      Eigen::VectorXd* coefficients) = 0;
  virtual void InverseTramsformToVelocity(
      const Eigen::VectorXd& coefficients, VFIELD2D* field) = 0;
  // Input a vector field, transform and output the basis coefficients.
  virtual void ForwardTransformtoFrequency(
      const VFIELD2D& field, Eigen::VectorXd* coefficients) = 0;
  // Fill the advection tensor.
  virtual void FillAdvectionTensor(
      std::vector<Adv_Tensor_Type> *Adv_tensor) = 0;
  // convert to vorticity field, using DCT.
  virtual void ComputeVorticity(const Eigen::VectorXd& coef, FIELD2D* vort) {
    std::cout << "not implemented" << std::endl;
  }
  
  int InverseLookUp(int i1, int i2) {
    return basislookup_.get()->InverseLookUp(i1,i2);
  }
  int LookUpBasis(int i, int mode) {
    return basislookup_.get()->LookUpBasis(i, mode);
  }
  int GetBasisNum() const {return root_num_basis_ * root_num_basis_;}
  int GetRootBasisNum() const {return root_num_basis_;}
  double GetDelta() const {return delta_;}
  // Verify the tensor is antisymmetric C_{ij}^k = - C_{ik}^j
  void VerifyAntisymmetric(const std::vector<Adv_Tensor_Type>& C);
  virtual void ComputeEigenValues( 
  std::vector<double>* eigs, std::vector<double>* eigs_inv, 
  std::vector<double>* eigs_inv_root) = 0;
  // This function enforce the advector tensor satisfy the antisymmetric relationship:
  // C_{i,j}^k = -C_{i,k}^j. As proved in the document. By this, we enforce the advection
  // to satisfy the energy conservation.
  virtual void EnforceAntisymmetric(std::vector<Adv_Tensor_Type>* Adv_tensor) = 0;

//  virtual void FillFourthOrderTensor(CTF::Tensor<double> *Q) = 0;
//  virtual void FillThirdOrderTensor(CTF::Tensor<double> *C) = 0;
//  virtual void FillVelocityTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) = 0;
  virtual void FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) = 0;
  virtual void ComputeBasisWeights(Eigen::VectorXd* weights) {
    weights->setOnes();
  }
  void MultiplyweightTensor(std::vector<Adv_Tensor_Type>& C, const Eigen::VectorXd& weights);
  void ElimiteDiagElements(std::vector<Adv_Tensor_Type>& C);
protected:
  inline bool check_in_range(int i) {
    if(i <= 0 || i > N_) {
      return false;
    } else 
      return true;
  }
  int N_;
  int root_num_basis_;
  double delta_;
  double inv_N_;
  double inv_pi_2_;
  double inv_pi_cube_8_;
  double inv_pi_square_;
  // std::unordered_map<std::pair<int,int>, int> indexmap_;
};

#endif
