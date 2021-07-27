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

#ifndef LAPLACIAN_BASIS_SET_3D_H
#define LAPLACIAN_BASIS_SET_3D_H

#include <fstream>
#include "Eigen"
#include <string>
#include <vector>

#include "3D/VECTOR3_FIELD_3D.h"
#include "Alg/VEC3.h"
#include "setting.h"
#include "util/util.h"

class LaplacianBasisSet3D {
  // des_basis_dim: desired basis_dimention.
public:
  LaplacianBasisSet3D(const int des_basis_dim, const int xRes, 
                      const int yRes, const int zRes, const std::string const_strategy):
                      des_basis_dim_(des_basis_dim),xRes_(xRes),
                      yRes_(yRes),zRes_(zRes), constant_init_strategy_(const_strategy),
                      dxyz_(Eigen::Vector3d(M_PI / static_cast<double>(xRes),
                                            M_PI / static_cast<double>(yRes),
                                            M_PI / static_cast<double>(zRes)))
  {// initialize the number of actually allocated basis to zero.
    basis_dim_ = 0;
    int_type_const_strategy_ = -1;
    if (constant_init_strategy_ == "principle_x") {
      int_type_const_strategy_ = 0;
    } else if (constant_init_strategy_ == "principle_y") {
      int_type_const_strategy_ = 1;
    } else if (constant_init_strategy_ == "principle_z") {
      int_type_const_strategy_ = 2;
    } else if (constant_init_strategy_ == "random") {
      int_type_const_strategy_ = 3;
    } else if (constant_init_strategy_ == "uniform") {
      int_type_const_strategy_ = 4;
    } else {
      std::cout << "laplacian_basis_set_3d.h " << __LINE__ << " FATAL: " <<  "Unknow constant_init_strategy: " << constant_init_strategy_ << std::endl; exit(0);
    }
  }
  ~LaplacianBasisSet3D(){}
  virtual void FillBasisNumerical(
      const Eigen::VectorXd& coefficients, VECTOR3_FIELD_3D* field) = 0;
  virtual int AllocateAllBasis() = 0;
  
  virtual void InverseTramsformToVelocity(
      const Eigen::VectorXd& coefficients, VECTOR3_FIELD_3D* field) = 0;
  // Input a vector field, transform and output the basis coefficients.
  virtual void ForwardTransformtoFrequency(
      const VECTOR3_FIELD_3D& field, Eigen::VectorXd* coefficients) = 0;
  // Verify the tensor is antisymmetric C_{ij}^k = - C_{ik}^j
  void VerifyAntisymmetric(const std::vector<Adv_Tensor_Type>& C);
  
  virtual void ComputeEigenValues(Eigen::VectorXd* eigenValue) = 0;
  
  virtual void FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) = 0;
  // For Dirichlet basis, the type is 0, one Neumann wall is 1, two neumann wall is 2.
  virtual void WriteBasis(std::ofstream& out) = 0;
  virtual int ReadBasis(std::ifstream& infile) = 0;
  int GetBasisDim() const {return basis_dim_;}
  virtual double ComputeBasisAt(const VEC3I& pos, const int basis_idx, const int mode) const = 0;
  virtual void PrintDebugInfo(const Eigen::VectorXd& coefficients) = 0;
  virtual void ComputeBasisWeights(Eigen::VectorXd* weights, const Real weightMultiplierC) = 0;
  void MultiplyweightTensor(std::vector<Adv_Tensor_Type>& C, const Eigen::VectorXd& weights);
protected:
  const int des_basis_dim_;
  int basis_dim_;
  const int xRes_;
  const int yRes_;
  const int zRes_;
  const Eigen::Vector3d dxyz_;
  const std::string constant_init_strategy_;
  
  int int_type_const_strategy_;
};

#endif  // LAPLACIAN_BASIS_SET_3D_H
