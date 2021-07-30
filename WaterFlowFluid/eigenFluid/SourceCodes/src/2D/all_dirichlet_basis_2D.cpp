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

#include <cmath>
#include "Eigen"
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "2D/all_dirichlet_basis_2D.h"
#include "setting.h"
#include "util/util.h"

// TODO: Modify the memory layout of VFIELD2D to do every thing inplace.
// Also Modify the memory structure to reduce the cost of shift the data
// in sine direction.
#ifndef USE_ROOT_LAMBDA
#define DEF_INV_LAMBDA(iddx, iddy) const double inv_lambda = inv_pi_2_ / (iddx*iddx + iddy*iddy);
#else 
#define DEF_INV_LAMBDA(iddx, iddy) const double inv_lambda = inv_pi_2_ / sqrt(iddx*iddx + iddy*iddy);
#endif

namespace {

double computeWeight(const int i1, const int i2) {
    //double weight = pow(i1*i1 + i2*i2, 0.01);
  //double x = (i1*i1 + i2*i2) / 3000.0;
  double weight = 1.0 - 0.9*(i1*i1 + i2*i2) / 12500.0;
  return weight;  
}

}  // namespace

void AllDirichletBasis2D::FillBasisNumerical(const Eigen::VectorXd& coefficients,
                                             VFIELD2D* field) {
  field->clear();
  // Iterate over all the coefficients.
  for (int i = 0; i < coefficients.size(); i++) {
    if (coefficients[i] == 0.0) {
      continue;
    }
    const double coef_val = coefficients[i];
    int i1 = basislookup_.get()->LookUpBasis(i, 0);
    int i2 = basislookup_.get()->LookUpBasis(i, 1);
    DEF_INV_LAMBDA(i1, i2)
    int index = 0;
    for (int j = 0; j < N_; j++) {
      for (int i = 0; i < N_; i++) {
        Real value_x = i2 * inv_lambda * sin((i + 0.5)*i1*delta_) * cos((j + 0.5)*i2*delta_);
        Real value_y = -i1 * inv_lambda * cos((i + 0.5)*i1*delta_) * sin((j + 0.5)*i2*delta_);
        (*field)[index][0] += value_x * coef_val;
        (*field)[index][1] += value_y * coef_val;
        index ++;
      }
    }
  }
}

double AllDirichletBasis2D::ComputeBasisAt(const int x_pos, const int y_pos,
                                           const int band, const int mode) const {
  if (x_pos < 0 || x_pos >= N_ || y_pos < 0 || y_pos >= N_) {
    std::cout << "all_dirichlet_basis_2D.cpp " << __LINE__ << " FATAL: " <<  "Index out of range" << std::endl; exit(0);
  }
  int i1 = basislookup_.get()->LookUpBasis(band, 0);
  int i2 = basislookup_.get()->LookUpBasis(band, 1);
  DEF_INV_LAMBDA(i1, i2)
  // X component.
  if (mode == 0){
    return i2 * inv_lambda * sin((x_pos + 0.5)*i1*delta_) * cos((y_pos + 0.5)*i2*delta_);
  } else {
    return -i1 * inv_lambda * cos((x_pos + 0.5)*i1*delta_) * sin((y_pos + 0.5)*i2*delta_);
  }
}

void AllDirichletBasis2D::ProjectVelocityNumerical(const VFIELD2D& field,
                                                   Eigen::VectorXd* coefficients) {
  if (coefficients->size() != num_basis_) {std::cout << "all_dirichlet_basis_2D.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  coefficients->setZero();
  // Iterate over all the coefficients.
  for (int k = 0; k < num_basis_; k++) {
    int i1 = LookUpBasis(k,0);
    int i2 = LookUpBasis(k,1);
    DEF_INV_LAMBDA(i1, i2)
    int index = 0;
    double contrib = 0.0;
    for (int j = 0; j < N_; j++) {
      for(int i = 0; i < N_; i++) {
        Real v_x = i2 * inv_lambda * sin ((i + 0.5) * i1 * delta_) * cos((j + 0.5) * i2 * delta_);
        Real v_y = - i1 * inv_lambda * cos((i + 0.5) * i1 * delta_) * sin((j + 0.5) * i2 * delta_);
        contrib += v_x * field[index][0] + v_y * field[index][1];
        index ++;
      }
    }
    (*coefficients)[k] = contrib * inv_N_ * inv_N_;
  }
}

void AllDirichletBasis2D::InverseTramsformToVelocity(
  const Eigen::VectorXd& coefficients, VFIELD2D* field) {
  field->clear();
  ClearTemp();
  // Iterate over all the coefficients and put it into the field.
  for (int i = 0; i < coefficients.size(); i ++) {
    if (coefficients[i] == 0) {
      continue;
    }
    int i1 = LookUpBasis(i, 0);
    int i2 = LookUpBasis(i, 1);
    DEF_INV_LAMBDA(i1, i2)
    // Process x field.
    temp_[i1 - 1 + i2 * N_] = coefficients[i] * 0.25 * inv_lambda * i2;
  }
  
  // Do the inverse transformation.
  fftw_plan plan_x_2D;
  // Inverse sine on x, inverse cosine on y.
  plan_x_2D = fftw_plan_r2r_2d(N_, N_, temp_, temp_, FFTW_REDFT01,
                                FFTW_RODFT01, FFTW_ESTIMATE);
  fftw_execute(plan_x_2D);  
  fftw_destroy_plan(plan_x_2D);
  // Put the value back into field.
  int index = 0;
  for (int j = 0; j < N_; j++) {
    for (int i = 0; i < N_; i++) {
      (*field)[index][0] = temp_[index];
      index ++;
    }
  }
  
  // Do same thing for y.
  ClearTemp();
  for (int i = 0; i < coefficients.size(); i++) {
    if (coefficients[i] == 0) {
      continue;
    }
    int i1 = LookUpBasis(i, 0);
    int i2 = LookUpBasis(i, 1);
    DEF_INV_LAMBDA(i1, i2)
    temp_[i1 + (i2 - 1)*N_] = -coefficients[i] * 0.25 * inv_lambda * i1;
  }
  // Do the transformation.
  fftw_plan plan_y_2D;
  // Inverse cosine transformation on y, inverse sine transformation on x.
  plan_y_2D = fftw_plan_r2r_2d(N_, N_, temp_, temp_,
                                FFTW_RODFT01, FFTW_REDFT01,FFTW_ESTIMATE);
  fftw_execute(plan_y_2D);  
  fftw_destroy_plan(plan_y_2D);
  // Put the value back.
  for (int j = 0; j < N_; j++) {
    for (int i = 0; i < N_; i++) {
      int index= i + j*N_;
      (*field)[index][1] = temp_[index];
    }
  }
}
// TODO: Compute basis dot product normalization term.
void AllDirichletBasis2D::ForwardTransformtoFrequency(
  const VFIELD2D& field, Eigen::VectorXd* coefficients) {
  if (coefficients->size() != num_basis_) {std::cout << "all_dirichlet_basis_2D.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  coefficients->setZero();
  ClearTemp();
  // Copy value to freq.
  int index = 0;
  for (int j = 0; j < N_; j++) {
    for (int i = 0; i < N_; i++) {
      temp_[index] = field[index][0];
      index++;
    }
  }
  // Do a transformation inplace on vx.
  fftw_plan plan_x_2D;

  // Sine transformation on X, cosine transformation on y.
  plan_x_2D = fftw_plan_r2r_2d(N_, N_, temp_, temp_, FFTW_REDFT10,
                                FFTW_RODFT10, FFTW_ESTIMATE);
  fftw_execute(plan_x_2D);
  fftw_destroy_plan(plan_x_2D);
  for (int i = 0; i < num_basis_; i++) {
    int i1 = LookUpBasis(i, 0);
    int i2 = LookUpBasis(i, 1);
    DEF_INV_LAMBDA(i1, i2)
    (*coefficients)[i] = temp_[i1 - 1 + i2 * N_] * inv_N_ * inv_N_ * 0.25 * inv_lambda * i2;
  }
  ClearTemp();
  // Copy Y value into Temp.
  index = 0;
  for (int j = 0; j < N_; j++) {
    for (int i = 0; i < N_; i++) {
      temp_[index] = field[index][1];
      index++;
    }
  }
  // Do a transformation inplace on vy.
  fftw_plan plan_y_2D;
  // Cosin transformation on X, sine transformation on y.
  plan_y_2D = fftw_plan_r2r_2d(N_, N_, temp_, temp_,FFTW_RODFT10,
                                FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(plan_y_2D);
  fftw_destroy_plan(plan_y_2D);
  for (int i = 0; i < num_basis_; i++) {
    int i1 = LookUpBasis(i, 0);
    int i2 = LookUpBasis(i, 1);
    DEF_INV_LAMBDA(i1, i2)
    (*coefficients)[i] -= temp_[i1 + (i2 - 1)*N_] * inv_N_ * inv_N_ * 0.25 * inv_lambda * i1;
  }
}
#undef DEF_INV_LAMBDA

// TODO: only compute part of the matrix since it's antisymmetric.
// Fill the tensor for all Dirichlet case.
void AllDirichletBasis2D::FillAdvectionTensor
    (std::vector<Adv_Tensor_Type> *Adv_tensor) {
  if(Adv_tensor != nullptr) { Adv_tensor->clear();};
  Adv_tensor->reserve(num_basis_);
  // Fill the tensor with zeros.
  for (int k = 0; k < num_basis_; k++) {
    Adv_Tensor_Type Ck(num_basis_, num_basis_);
    Ck.setZero();
    Adv_tensor->emplace_back(Ck);
  }
  for (int d1 = 0; d1 < num_basis_; d1++) {
    int i1 = LookUpBasis(d1,0);
    int i2 = LookUpBasis(d1,1);
#ifndef USE_ROOT_LAMBDA
    double inv_lambda_i = 1.0/(i1*i1 + i2*i2);
#else
    double inv_lambda_i = 1.0/sqrt(i1*i1 + i2*i2);
#endif
    for (int d2 = 0;d2 < num_basis_; d2++) {
      int j1 = LookUpBasis(d2,0);
      int j2 = LookUpBasis(d2,1);

      int k1 = d1;   // equals index of (i1, i2)
      int k2 = d2;   // equals index of (j1, j2)
      int antipairs[9][2];
      antipairs[0][0] = i1+j1; antipairs[0][1] = i2+j2;
      antipairs[1][0] = i1+j1; antipairs[1][1] = i2-j2;
      antipairs[2][0] = i1+j1; antipairs[2][1] = j2-i2;
      antipairs[3][0] = i1-j1; antipairs[3][1] = i2+j2;
      antipairs[4][0] = i1-j1; antipairs[4][1] = i2-j2;
      antipairs[5][0] = i1-j1; antipairs[5][1] = j2-i2;
      antipairs[6][0] = j1-i1; antipairs[6][1] = i2+j2;
      antipairs[7][0] = j1-i1; antipairs[7][1] = i2-j2;
      antipairs[8][0] = j1-i1; antipairs[8][1] = j2-i2;
      for (int c = 0; c < 9; c++) {
        int i = antipairs[c][0];
        int j = antipairs[c][1];

        int index = InverseLookUp(i,j);
        if (index != -1) {
          double coef = compute_coefficient_terms(i1,j1,i2,j2,c) * inv_lambda_i;
          //Adv_tensor.Ck[index].set(k1,k2,coef);
          IncrementMatrixEntry(&(*Adv_tensor)[index],k1,k2, coef);
        }
      }
    }
  }
    // Make makeCompressed.
  for (int k = 0; k < num_basis_; k++) {    
    (*Adv_tensor)[k].makeCompressed();
  }
}

double AllDirichletBasis2D::compute_coefficient_terms
    (int i1, int j1, int i2, int j2, int mode) {
      
  switch(mode) {
    case 0 : return 0.25 * (i1*j2 - i2*j1);  // mode 1
             break;
    case 1 : return 0.25 * (i1*j2 + i2*j1);  // mode 2
             break;
    case 2 : return 0.25 * -(i1*j2 + i2*j1);  // mode 3
             break;
    case 3 : return 0.25 * -(i1*j2 + i2*j1);  // mode 4
             break;
    case 4 : return 0.25 * (-i1*j2 + i2*j1);  // mode 5
             break;
    case 5 : return 0.25 * (i1*j2 - i2*j1);   // mode 6
             break;
    case 6 : return 0.25 * (i1*j2 + i2*j1);   // mode 7
             break;
    case 7 : return 0.25 * (i1*j2 - i2*j1);   // mode 8
             break;
    case 8 : return 0.25 * (-i1*j2 + i2*j1);  // mode 9
             break;
    default:
      std::cout << "all_dirichlet_basis_2D.cpp " << __LINE__ << " FATAL: " <<  "Unknow mode." << std::endl; exit(0);
      return 0.0;
  }
  return 0.0;
}

void AllDirichletBasis2D::EnforceAntisymmetric(std::vector<Adv_Tensor_Type>* Adv_tensor) {
    for (int k = 0; k < num_basis_; k++) {
      for (int i = 0; i < num_basis_; i ++) {
        for (int j = 0; j < num_basis_; j++) {
          double C_ijk = AccessMatrix((*Adv_tensor)[k], i, j);
          if (C_ijk!= 0) {
            double C_ikj = AccessMatrix((*Adv_tensor)[j], i, k);
            double sym = 0.5 * (C_ijk + C_ikj);
            if (sym > 1e-10) {
              std::cout <<  "Non antisymmetric entry found. value: " << sym << std::endl;
            }
            double anti_sym = 0.5 * (C_ijk - C_ikj);
            SetMatrixEntry(&(*Adv_tensor)[k], i, j, anti_sym);
            SetMatrixEntry(&(*Adv_tensor)[j], i, k, -anti_sym);
          }
        }
      }
    }
  return;
}

void AllDirichletBasis2D::ComputeEigenValues(
  std::vector<double>* eigs, std::vector<double>* eigs_inv, 
  std::vector<double>* eigs_inv_root) {
  for (int i = 0; i < num_basis_; i++) {
    int k1 = LookUpBasis(i,0);
    int k2 = LookUpBasis(i,1);
    (*eigs)[i] = (k1*k1 + k2*k2);
    (*eigs_inv)[i] = 1.0/(k1*k1+k2*k2);
    (*eigs_inv_root)[i] = 1.0/sqrt(k1*k1+k2*k2);;
  }
}

// Check wether a pair of integers within basis range.
bool AllDirichletBasis2D::WithinRange(const int i1, const int i2) {
  if (i1 <= 0 || i2 <= 0 || i1 > root_num_basis_ || i2 > root_num_basis_) {
    return false;
  } else {
    return true;
  }
}

double AllDirichletBasis2D::GenerateVariationalTerms(
    const int j1, const int j2, const int k1, const int k2, const int index) {
  double j1k2 = j1*k2;
  double j2k1 = j2*k1;
  
  switch (index) {
    case 0: // ++
      return j1k2 - j2k1;
      break;
    case 1: // +-
      return j1k2 + j2k1;
      break;
    case 2: // +- neg
      return -j1k2 - j2k1;
      break;
    case 3: // -+
      return -j1k2 - j2k1;
      break;
    case 4: // -+ neg
      return j1k2 + j2k1;
      break;
    case 5: // --
      return -j1k2 + j2k1;
      break;
    case 6: // -- neg
      return j1k2 - j2k1;
      break;
    case 7: // -- neg
      return j1k2 - j2k1;
      break;
    case 8: // -- pos
      return -j1k2 + j2k1;
      break;
    default:
      std::cout << "all_dirichlet_basis_2D.cpp " << __LINE__ << " FATAL: " <<  "Something is wrong" << std::endl; exit(0);
      break;
  }
}

void AllDirichletBasis2D::FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) {
  if(Adv_tensor != nullptr) { Adv_tensor->clear();};
  Adv_tensor->reserve(num_basis_);
  // Fill the tensor with zeros.
  for (int k = 0; k < num_basis_; k++) {
    Adv_Tensor_Type Ck(num_basis_, num_basis_);
    Ck.setZero();
    Adv_tensor->emplace_back(Ck);
  }
  const double pow_w = - 0.00;
  Real print_percentage = 0; 
  for (int k = 0; k < num_basis_; k++) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    
    int k1 = LookUpBasis(k, 0);
    int k2 = LookUpBasis(k, 1);
    int kk = k1*k1 + k2*k2;
    double weightk = computeWeight(k1, k2);
    const Real percentage = static_cast<Real>(k) / num_basis_ * 100.0;
    if (std::abs(percentage - print_percentage) > 5) {
      print_percentage = percentage;
      std::cout <<  "% "  << print_percentage << "Tensor computed." << std::endl;
    }
    
    for (int j = 0; j < num_basis_; j++) {
      int j1 = LookUpBasis(j, 0);
      int j2 = LookUpBasis(j, 1);
      int jj = j1*j1 + j2*j2;
      int antipairs[9][2];
      antipairs[0][0] = j1+k1; antipairs[0][1] = j2+k2;
      
      antipairs[1][0] = j1+k1; antipairs[1][1] = j2-k2;
      antipairs[2][0] = j1+k1; antipairs[2][1] = k2-j2;
      
      antipairs[3][0] = j1-k1; antipairs[3][1] = j2+k2;
      antipairs[4][0] = k1-j1; antipairs[4][1] = j2+k2;
      
      antipairs[5][0] = j1-k1; antipairs[5][1] = j2-k2;
      antipairs[6][0] = j1-k1; antipairs[6][1] = k2-j2;
      antipairs[7][0] = k1-j1; antipairs[7][1] = j2-k2;
      antipairs[8][0] = k1-j1; antipairs[8][1] = k2-j2;
      double weightj = computeWeight(j1, j2);
      for (int c = 0; c < 9; c++) {
        int i_ = antipairs[c][0];
        int j_ = antipairs[c][1];
        int ii = i_*i_ + j_*j_;
        int index = InverseLookUp(i_,j_);
        if (index != -1) {
          double inv_lambda_k = 1.0/sqrt(k1*k1 + k2*k2);
          double inv_lambda_j = 1.0/sqrt(j1*j1+j2*j2);
          double lambda_i = sqrt(i_*i_ + j_*j_);
          double weighti = computeWeight(i_, j_);
       //   const double invpi = 2.0 /( M_PI * M_PI * M_PI);
          double coef = - GenerateVariationalTerms(j1,j2,k1,k2,c) *
              lambda_i * inv_lambda_j * inv_lambda_k * 0.25 * inv_pi_2_ * weighti*weightj*weightk;
          // IncrementMatrixEntry(&(*Adv_tensor)[k],index, j, coef);
          //if (kk >= ii && kk >= jj) {
            tripletList.push_back(T(index, j, coef));
          //}
        }
      }
    }
    (*Adv_tensor)[k].setFromTriplets(tripletList.begin(), tripletList.end());
  }
  // Make makeCompressed.
  for (int k = 0; k < num_basis_; k++) {    
    (*Adv_tensor)[k].makeCompressed();
  }
}
