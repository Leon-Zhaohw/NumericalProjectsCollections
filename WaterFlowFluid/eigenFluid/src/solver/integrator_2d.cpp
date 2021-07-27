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
// Modified from Eigen source code.
 */

#include "solver/integrator_2d.h"
#include "Eigen"
#include "Sparse"

void Integrator2D::conjugate_gradient_NE(const Eigen::MatrixXd& mat, const Eigen::VectorXd& rhs,
                                         Eigen::VectorXd& x,
                        const Eigen::DiagonalPreconditioner<Eigen::VectorXd::RealScalar>& precond,
                        int& iters,
                        Eigen::VectorXd::RealScalar& tol_error)
{
  using std::sqrt;
  using std::abs;
  typedef typename Eigen::VectorXd::RealScalar RealScalar;
  typedef typename Eigen::VectorXd::Scalar Scalar;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> VectorType;
  
  RealScalar tol = tol_error;
  int maxIters = iters;
  
  int n = mat.cols();

  VectorType residual = rhs - mat.transpose()*(mat * x); //initial residual

  RealScalar rhsNorm2 = rhs.squaredNorm();
  if(rhsNorm2 == 0) 
  {
    x.setZero();
    iters = 0;
    tol_error = 0;
    return;
  }
  RealScalar threshold = tol*tol*rhsNorm2;
  RealScalar residualNorm2 = residual.squaredNorm();
  if (residualNorm2 < threshold)
  {
    iters = 0;
    tol_error = sqrt(residualNorm2 / rhsNorm2);
    return;
  }
  
  VectorType p(n);
  p = precond.solve(residual);      //initial search direction

  VectorType z(n), tmp(n), tmp1(n);
  RealScalar absNew = Eigen::numext::real(residual.dot(p));  // the square of the absolute value of r scaled by invM
  int i = 0;
  while(i < maxIters)
  {
    tmp1.noalias() = mat * p;
    tmp.noalias() = mat.transpose() * tmp1;              // the bottleneck of the algorithm

    Scalar alpha = absNew / p.dot(tmp);   // the amount we travel on dir
    x += alpha * p;                       // update solution
    residual -= alpha * tmp;              // update residue
    
    residualNorm2 = residual.squaredNorm();
    if(residualNorm2 < threshold)
      break;
    
    z = precond.solve(residual);          // approximately solve for "A z = residual"

    RealScalar absOld = absNew;
    absNew = Eigen::numext::real(residual.dot(z));     // update the absolute value of r
    RealScalar beta = absNew / absOld;            // calculate the Gram-Schmidt value used to create the new search direction
    p = z + beta * p;                             // update search direction
    i++;
  }
  tol_error = sqrt(residualNorm2 / rhsNorm2);
  iters = i;
}
