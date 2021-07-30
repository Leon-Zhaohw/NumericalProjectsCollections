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
#include <fftw3.h>
//// #include <glog/logging.h>
#include <math.h>
#include <stdlib.h>

#include "2D/two_neumann_x_2D.h"
#include "2D/VFIELD2D.h"
// Test file for new implemented two Neumann 2D class.

namespace {

void testTransform() {
  const int root_basis_dim = 8;
  int basis_dim = root_basis_dim * root_basis_dim;
  const int N = 128;
  const double dx = M_PI / static_cast<double>(N);
  
  TwoNeumannX2D basisSet(N, root_basis_dim);
  // std::cout <<  basisSet.AllocateAllBasis() << std::endl;
  VFIELD2D numericalField(N, N);
  VFIELD2D transformField(N, N);
  Eigen::VectorXd randCoef(basis_dim);
  Eigen::VectorXd tranCoef(basis_dim);
  randCoef.setZero();
  tranCoef.setZero();
  for (int i = 0; i < basis_dim; i++) {
    randCoef[i] = rand() / static_cast<double>(RAND_MAX);
  }
  
  basisSet.FillBasisNumerical(randCoef, &numericalField);
  basisSet.InverseTramsformToVelocity(randCoef, &transformField);
  basisSet.ForwardTransformtoFrequency(transformField, &tranCoef);
  std::cout <<  "transformField: " << transformField.DotProduct(transformField)*dx*dx << std::endl;
  std::cout <<  "numericalField: " << numericalField.DotProduct(numericalField)*dx*dx << std::endl;
  std::cout <<  "random Coeff: " << randCoef.squaredNorm() << std::endl;
  std::cout <<  "transf Coeff: " << tranCoef.squaredNorm() << std::endl;
}

void testTensor() {
  const int root_basis_dim = 8;
  int basis_dim = root_basis_dim * root_basis_dim;
  const int N = 128;
  TwoNeumannX2D basisSet(N, root_basis_dim);
  // std::cout <<  basisSet.AllocateAllBasis() << std::endl;
  std::vector<Adv_Tensor_Type> Adv_tensor;
  basisSet.FillVariationalTensor(&Adv_tensor);
  basisSet.VerifyAntisymmetric(Adv_tensor);
  uint64_t nnzero;
  for (int i = 0; i < basis_dim; i++) {
    nnzero += Adv_tensor[i].nonZeros();
  }
  std::cout <<  "nnzero: " << nnzero << std::endl;
}

}  // namespace

int main(int argc, char ** argv) {
  //google::ParseCommandLineFlags(&argc, &argv, true);
  //google::InitGoogleLogging(argv[0]);
  
  // testTransform();
  testTensor();
  
  return 0;
}
