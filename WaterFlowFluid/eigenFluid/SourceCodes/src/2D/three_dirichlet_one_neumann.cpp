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
#include "Sparse"
#include <fftw3.h>
#include <omp.h>
// #include <glog/logging.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "2D/three_dirichlet_one_neumann.h"
#include "setting.h"
#include "util/util.h"

#ifndef USE_ROOT_LAMBDA
#define DEF_INV_LAMBDA(iddx, iddy) const double inv_lambda = inv_pi_2_ / (iddx*iddx + iddy*iddy);
#else 
#define DEF_INV_LAMBDA(iddx, iddy) const double inv_lambda = inv_pi_2_ / sqrt(iddx*iddx + iddy*iddy);
#endif

void ThreeDirichletOneNeumann::FillBasisNumerical(
  const Eigen::VectorXd& coefficients, VFIELD2D* field) {
  field->clear();
  // Iterate over all the coefficients.
  for (int i = 0; i < coefficients.size(); i++) {
    if (coefficients[i] == 0.0) {
      continue;
    }
    const double coef_val = coefficients[i];
    int a = basislookup_.get()->LookUpBasis(i, 0);
    Real i1 = a - 0.5;
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

double ThreeDirichletOneNeumann::ComputeBasisAt(const int x_pos, const int y_pos,
                                               const int band, const int mode) const {
  if (x_pos < 0 || x_pos >= N_ || y_pos < 0 || y_pos >= N_) {
    std::cout << "three_dirichlet_one_neumann.cpp " << __LINE__ << " FATAL: " <<  "Index out of range" << std::endl; exit(0);
  }
  int a = basislookup_.get()->LookUpBasis(band, 0);
  double i1 = a - 0.5;
  int i2 = basislookup_.get()->LookUpBasis(band, 1);
  DEF_INV_LAMBDA(i1, i2)
  // X component.
  if (mode == 0){
    return i2 * inv_lambda * sin((x_pos + 0.5)*i1*delta_) * cos((y_pos + 0.5)*i2*delta_);
  } else {
    return -i1 * inv_lambda * cos((x_pos + 0.5)*i1*delta_) * sin((y_pos + 0.5)*i2*delta_);
  }
}

void ThreeDirichletOneNeumann::ProjectVelocityNumerical(const VFIELD2D& field,
                                                        Eigen::VectorXd* coefficients) {
  if (coefficients->size() != num_basis_) {std::cout << "three_dirichlet_one_neumann.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  coefficients->setZero();
  // Iterate over all the coefficients.
  for (int k = 0; k < num_basis_; k++) {
    int a = LookUpBasis(k,0);
    double i1 = a - 0.5;
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

void ThreeDirichletOneNeumann::ForwardTransformtoFrequency(
    const VFIELD2D& field, Eigen::VectorXd* coefficients) {
  if (coefficients->size() != num_basis_) {std::cout << "three_dirichlet_one_neumann.cpp " << __LINE__ << " FATAL: "  << std::endl; exit(0);}
  coefficients->setZero();
  ClearTemp();
  // Copy value to freq.
  int field_idx = 0;
  for (int j = 0; j < N_; j++) {
    for (int i = 0; i < N_; i++) {
      int index = i + j*2*N_;
      temp_[index] = field[field_idx][0];
      field_idx ++;
    }
  }
  // Do a transformation inplace on vx.
  fftw_plan plan_x_2D;

  // Sine transformation on X, cosine transformation on y.
  plan_x_2D = fftw_plan_r2r_2d(N_, 2*N_, temp_, temp_, FFTW_REDFT10,
                                FFTW_RODFT10, FFTW_ESTIMATE);
  fftw_execute(plan_x_2D);
  fftw_destroy_plan(plan_x_2D);
  for (int i = 0; i < num_basis_; i++) {
    int a = LookUpBasis(i, 0);
    Real i1 = a - 0.5;
    int i2 = LookUpBasis(i, 1);
    DEF_INV_LAMBDA(i1, i2)
    (*coefficients)[i] = temp_[2*a - 2 + i2*2*N_] * inv_N_ * inv_N_ * 0.25 * inv_lambda * i2;
  }
  // Do same thing on vy.
  ClearTemp();
  field_idx = 0;
  for (int j = 0; j < N_; j++) {
    for (int i = 0; i < N_; i++) {
      int index = i + j*2*N_;
      temp_[index] = field[field_idx][1];
      field_idx ++;
    }
  }
  // Do a transformation inplace on vx.
  fftw_plan plan_y_2D;

  // Cosine transformation on X, sine transformation on y.
  plan_y_2D = fftw_plan_r2r_2d(N_, 2*N_, temp_, temp_,FFTW_RODFT10,
                                FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(plan_y_2D);
  fftw_destroy_plan(plan_y_2D);
  for (int i = 0; i < num_basis_; i++) {
    int a = LookUpBasis(i, 0);
    Real i1 = a - 0.5;
    int i2 = LookUpBasis(i, 1);
    DEF_INV_LAMBDA(i1, i2)
    (*coefficients)[i] -= temp_[2*a - 1 + (i2-1)*2*N_] * inv_N_ * inv_N_ * 0.25 * inv_lambda * i1;
  }
}

void ThreeDirichletOneNeumann::InverseTramsformToVelocity(
  const Eigen::VectorXd& coefficients, VFIELD2D* field) {
  field->clear();
  ClearTemp();
  // Iterate over all the coefficients and put it into the field.
  for (int i = 0; i < coefficients.size(); i ++) {
    if (coefficients[i] == 0) {
      continue;
    }
    int a = LookUpBasis(i, 0);
    Real i1 = a - 0.5;
    int i2 = LookUpBasis(i, 1);
    DEF_INV_LAMBDA(i1, i2)
    // Process x field.
    temp_[2*a - 2 +
    i2 *2 * N_] = coefficients[i] * 0.25 * inv_lambda * i2;
  }
  // Do the inverse transformation.
  fftw_plan plan_x_2D;
  // Inverse sine on x, inverse cosine on y.
  plan_x_2D = fftw_plan_r2r_2d(N_, 2*N_, temp_, temp_, FFTW_REDFT01,
                                FFTW_RODFT01, FFTW_ESTIMATE);
  fftw_execute(plan_x_2D);  
  fftw_destroy_plan(plan_x_2D);
  // Put the value back into field.
  int index = 0;
  for (int j = 0; j < N_; j++) {
    for (int i = 0; i < N_; i++) {
      int temp_ind = i + j * 2 * N_;
      (*field)[index][0] = temp_[temp_ind];
      index ++;
    }
  }
  
  // Do same thing for y.
  ClearTemp();
  for (int i = 0; i < coefficients.size(); i++) {
    if (coefficients[i] == 0) {
      continue;
    }
    int a = LookUpBasis(i, 0);
    Real i1 = a - 0.5;
    int i2 = LookUpBasis(i, 1);
    DEF_INV_LAMBDA(i1, i2)
    temp_[2 * a - 1 + (i2 - 1)*2*N_] = -coefficients[i] * 0.25 * inv_lambda * i1;
  }
  // Do the transformation.
  fftw_plan plan_y_2D;
  // Inverse cosine transformation on y, inverse sine transformation on x.
  plan_y_2D = fftw_plan_r2r_2d(N_, 2*N_, temp_, temp_,
                                FFTW_RODFT01, FFTW_REDFT01,FFTW_ESTIMATE);
  fftw_execute(plan_y_2D);  
  fftw_destroy_plan(plan_y_2D);
  // Put the value back.
  index = 0;
  for (int j = 0; j < N_; j++) {
    for (int i = 0; i < N_; i++) {
      int temp_idx = i + j * 2 * N_;
      (*field)[index][1] = temp_[temp_idx];
      index ++;
    }
  }
}

#undef DEF_INV_LAMBDA

void ThreeDirichletOneNeumann::FillAdvectionTensor(
    std::vector<Adv_Tensor_Type> *Adv_tensor) {
  if(Adv_tensor != nullptr) { Adv_tensor->clear();};
  Adv_tensor->reserve(num_basis_);
  // Fill the tensor with zeros.
  for (int k = 0; k < num_basis_; k++) {
    Adv_Tensor_Type Ck(num_basis_, num_basis_);
    Ck.setZero();
    Adv_tensor->emplace_back(Ck);
  }
  for (int d1 = 0;d1 < num_basis_;d1++) {
    // i1 = a - 0.5
    int a = LookUpBasis(d1,0);
    double i1 = a - 0.5;
    int i2 = LookUpBasis(d1,1);
    double i1_2 = i1*i1;
    // (a, i2) = d1;
#ifndef USE_ROOT_LAMBDA
    double inv_lambda_i = 1.0/(i1*i1 + i2*i2);
#else
    double inv_lambda_i = 1.0/sqrt(i1*i1 + i2*i2);
#endif
    for (int d2 = 0;d2 < num_basis_;d2++) {
      int b = LookUpBasis(d2,0);
      double j1 = b - 0.5;
      int j2 = LookUpBasis(d2,1);
      // (b, j2) = d2;
      // Iterate over every possible c.
      for(int c=1; c<=root_num_basis_; c++) {
        int k2 = i2 + j2;
        int index = InverseLookUp(c,k2);
        if (index != -1) {
          // Fill i2+j2 mode.
          double result = compute_coefficient_terms(
              a, b, c, i1*j2, i2*j1, inv_lambda_i, 0);
          IncrementMatrixEntry(&(*Adv_tensor)[index],d1, d2, result);
        }
        k2 = i2 - j2;
        index = InverseLookUp(c,k2);
        if (index != -1){
          // Fill i2-j2 mode.
          double result = compute_coefficient_terms(
              a, b, c, i1*j2, i2*j1, inv_lambda_i, 1);
          IncrementMatrixEntry(&(*Adv_tensor)[index],d1, d2, result);
        }
        k2 = j2 - i2;
        index = InverseLookUp(c,k2);
        if (index != -1) {
          // Fill j2-i2 mode.
          double result = compute_coefficient_terms(
              a, b, c, i1*j2, i2*j1, inv_lambda_i, 2);
          IncrementMatrixEntry(&(*Adv_tensor)[index],d1, d2, result);
        }
      }
    }
  }
    // Make makeCompressed.
  for (int k = 0; k < num_basis_; k++) {    
    (*Adv_tensor)[k].makeCompressed();
  }
}

double ThreeDirichletOneNeumann::compute_coefficient_terms(
    int a, int b, int c, double i1j2, double i2j1, double inv_lambda_i,
    int mode) {
  const double inv_pi = 1.0 / M_PI;
  double F1 = 1.0 / (a + b - c - 0.5);
  double F2 = 1.0 / (a + b + c - 1.5);
  double F3 = 1.0 / (a - b - c + 0.5);
  double F4 = 1.0 / (a - b + c - 0.5);

  if ((a + b - c) % 2 == 0) {
    F1 = - F1;
  }
  if ((a + b + c) % 2 != 0) {
    F2 = - F2;
  }
  if ((a - b - c) % 2 != 0) {
    F3 = - F3;
  }
  if ((a - b + c) % 2 == 0) {
    F4 = - F4;
  }
  double AA = i1j2 - i2j1;
  double BB = i1j2 + i2j1;
  double result = 0;
  // i2 + j2 mode,
  if (mode == 0) {
    result = AA*F1 - AA*F2 - BB*F3 + BB*F4;  
  }
  // i2 - j2 mode.
  if (mode == 1) {
    result = BB*F1 - BB*F2 - AA*F3 + AA*F4; 
  }
  // j2 - i2 mode.
  if (mode == 2) {
    result = -BB*F1 + BB*F2 + AA*F3 - AA*F4;
  }
  result = 0.25*inv_lambda_i*inv_pi*result;
  return result;
}

void ThreeDirichletOneNeumann::EnforceAntisymmetric(std::vector<Adv_Tensor_Type>* Adv_tensor) {
    for (int k = 0; k < num_basis_; k++) {
      for (int i = 0; i < num_basis_; i ++) {
        for (int j = 0; j < num_basis_; j++) {
          double C_ijk = AccessMatrix((*Adv_tensor)[k], i, j);
          if (C_ijk!= 0) {
            double C_ikj = AccessMatrix((*Adv_tensor)[j], i, k);
            double sym = 0.5 * (C_ijk + C_ikj);
            if (sym > 1e-10) {
//              std::cout <<  "Non antisymmetric entry found. value: " << sym << std::endl;
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

void ThreeDirichletOneNeumann::ComputeEigenValues(
    std::vector<double>* eigs, std::vector<double>* eigs_inv, 
    std::vector<double>* eigs_inv_root) {
  
  for (int i = 0; i < num_basis_; i++) {
    double k1 = LookUpBasis(i,0) - 0.5;
    double k2 = LookUpBasis(i,1);
    (*eigs)[i] = (k1*k1 + k2*k2);
    (*eigs_inv)[i] = 1.0/(k1*k1+k2*k2);
    (*eigs_inv_root)[i] = 1.0/sqrt(k1*k1+k2*k2);
  }
}

double ThreeDirichletOneNeumann::ComputeVarVariationalTerms(
    const double j1k2, const double j2k1, const int i1, const int j1, const int k1, const int index) {
  double F1 = 1.0 / (i1 - j1 - k1 + 0.5);
  double F2 = 1.0 / (i1 + j1 + k1 - 1.5);
  double F3 = 1.0 / (i1 - j1 + k1 - 0.5);
  double F4 = 1.0 / (i1 + j1 - k1 - 0.5);
  if ((i1 - j1 - k1) % 2 != 0) {
    F1 = -F1;
  }
  if ((i1 + j1 + k1) % 2 != 0) {
    F2 = -F2;
  }
  if ((i1 - j1 + k1) % 2 == 0) {
    F3 = -F3;
  }
  if ((i1 + j1 - k1) % 2 == 0) {
    F4 = -F4;
  }
  double A = 0,B = 0,C = 0,D = 0;
  switch (index) {
    case 0: // i2 = j2 + k2
      A = j1k2 - j2k1;
      C = -j1k2 - j2k1;
      return A*(F1-F2) + C*(F3-F4);
      break;
    case 1: // i2 = j2 - k2
      B = j1k2 + j2k1;
      D = -j1k2 + j2k1;
      return B*(F1-F2) + D*(F3-F4);
      break;
    case 2:
       B = j1k2 + j2k1;
       D = -j1k2 + j2k1;
       return -B*(F1-F2) + -D*(F3-F4);
      break;
    default:
      std::cout << "three_dirichlet_one_neumann.cpp " << __LINE__ << " FATAL: " <<  "Something is wrong." << std::endl; exit(0);
      return 0;
      break;
  }
}
#define FILL_TENSOR_INDX(idxx_)\
index = InverseLookUp(it,i2);\
if (index != -1) {\
  double lambda_i = sqrt(i1*i1 + i2*i2);\
  double coef = -ComputeVarVariationalTerms(j1*k2, j2*k1, it, jt, kt, idxx_) *\
      0.5*inv_pi_square_*inv_lambda_j*inv_lambda_k *lambda_i;\
  IncrementMatrixEntry(&(*Adv_tensor)[k],index, j, coef);\
}

void ThreeDirichletOneNeumann::FillVariationalTensor(std::vector<Adv_Tensor_Type> *Adv_tensor) {
 if(Adv_tensor != nullptr) { Adv_tensor->clear();};
  Adv_tensor->reserve(num_basis_);
  // Fill the tensor with zeros.
  for (int k = 0; k < num_basis_; k++) {
    Adv_Tensor_Type Ck(num_basis_, num_basis_);
    Ck.setZero();
    Adv_tensor->emplace_back(Ck);
  }
  Real print_percentage = 0; 
  for (int k = 0; k < num_basis_; k++) {
    int kt = LookUpBasis(k, 0);
    int k2 = LookUpBasis(k, 1);
    double k1 = kt - 0.5;
    double inv_lambda_k = 1.0 / sqrt(k1*k1 + k2*k2);
    const Real percentage = static_cast<Real>(k) / num_basis_ * 100.0;
    if (std::abs(percentage - print_percentage) > 5) {
      print_percentage = percentage;
      std::cout <<  "% "  << print_percentage << "Tensor computed." << std::endl;
    }
    
    for (int j = 0; j < num_basis_; j++) {
      int jt = LookUpBasis(j, 0);
      int j2 = LookUpBasis(j, 1);
      double j1 = jt - 0.5;
      double inv_lambda_j = 1.0 / sqrt(j1*j1 + j2*j2);
      // Dense in i1.
      int i2 = 0;
      int index = 0;
      for (int it = 1; it <= root_num_basis_; it++) {
        double i1 = it - 0.5;
        i2 = j2 + k2;
        FILL_TENSOR_INDX(0)
        i2 = j2 - k2;
        FILL_TENSOR_INDX(1)
        i2 = k2 - j2;
        FILL_TENSOR_INDX(2)
      }
    }
  }
  // Make makeCompressed.
  for (int k = 0; k < num_basis_; k++) {    
    (*Adv_tensor)[k].makeCompressed();
  }
}
#undef FILL_TENSOR_INDX
