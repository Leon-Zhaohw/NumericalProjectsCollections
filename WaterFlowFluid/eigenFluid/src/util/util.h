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

#ifndef UTIL_H
#define UTIL_H

#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#include <GL/glu.h>
#endif

#include <math.h>
#include "Eigen"
#include <vector>
#include <fstream>

#include "setting.h"
#include "Alg/VEC3.h"

// #define PI_CUBE M_PI*M_PI*M_PI
// #define PI_SQUARE M_PI*M_PI
const double PI_CUBE = M_PI*M_PI*M_PI;
const double PI_SQUARE = M_PI*M_PI;
template<typename T> T inline clamp(T in, T _min, T _max) {
  if(in<_min)
    return _min;
  if(in>_max)
    return _max;
  return in;
}

template <typename T> T inline getsgn(T in) {
  if (in < 0) return -1;
    else return 1; 
}

// Generate a randomized normalized vector (norm is 1).
void GetNormalizedRandomVector(const int size, Eigen::VectorXd* vect);

inline int LookUpBasisFast(const int i, const int basis_num_root_, const int mode) {
  // Look up X.
 if (mode == 0) {
   return i % basis_num_root_ + 1;
 }
 // Look up Y.
 else if (mode == 1) {
   return i / basis_num_root_ + 1;
 }
};

inline double ComputeBasisAtFast(const int x_pos, const int y_pos,
                                 const int band, const int mode, const int root_num_basis_, const double  delta_) {
    int i1 = LookUpBasisFast(band, root_num_basis_, 0);
    int i2 = LookUpBasisFast(band, root_num_basis_, 1);
#ifndef USE_ROOT_LAMBDA
  const double inv_lambda = -1.0 / (i1*i1 + i2*i2);
#else
  const double inv_lambda = -1.0 / sqrt(i1*i1 + i2*i2);
#endif
  // X component.
  if (mode == 0){
    return i2 * inv_lambda * sin((x_pos + 0.5)*i1*delta_) * cos((y_pos + 0.5)*i2*delta_);
  } else {
    return -i1 * inv_lambda * cos((x_pos + 0.5)*i1*delta_) * sin((y_pos + 0.5)*i2*delta_);
  }
};

inline double AccessMatrix(const Adv_Tensor_Type& C, const int i , const int j) {
  double result = 0.;
#ifndef Adv_Tensor_Sparse
  result = C(i,j);
#else
  result = C.coeff(i,j);
#endif
  return result;
}

inline void SetMatrixEntry(Adv_Tensor_Type* C, const int i, const int j,const double& val) {
#ifndef Adv_Tensor_Sparse
  (*C)(i,j) = val;
#else
  (*C).coeffRef(i,j) = val;
#endif
}

inline void IncrementMatrixEntry(Adv_Tensor_Type* C, const int i, const int j,const double& val) {
#ifndef Adv_Tensor_Sparse
  (*C)(i,j) += val;
#else
  (*C).coeffRef(i,j) += val;
#endif
}

inline bool WithinRangeINT(const int a, const int min, const int max) {
  if ( a>= min && a <=max ) {
    return true;
  } else {
    return false;
  }
}

template<class Matrix>
void writeEigenDense_binary(std::ofstream& out, const Matrix& matrix){
    // std::ofstream out(filename,ios::out | ios::binary | ios::trunc);
    typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
    out.write((char*) (&rows), sizeof(typename Matrix::Index));
    out.write((char*) (&cols), sizeof(typename Matrix::Index));
    out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
    // out.close();
}

template<class Matrix>
void readEigenDense_binary(std::ifstream& in, Matrix& matrix){
    // std::ifstream in(filename,ios::in | std::ios::binary);
    typename Matrix::Index rows=0, cols=0;
    in.read((char*) (&rows),sizeof(typename Matrix::Index));
    in.read((char*) (&cols),sizeof(typename Matrix::Index));
    matrix.resize(rows, cols);
    in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
   // in.close();
}

inline float GenerateUnifomRnd() {
  return std::rand() / static_cast<float>(RAND_MAX);
}

void CompareTensor(const std::vector<Adv_Tensor_Type>& tensor1, 
                   const std::vector<Adv_Tensor_Type>& tensor2);
bool PointInBox(const VEC3F& point, const VEC3F& box_center, const VEC3F& length_);
bool BoxIntersect(const VEC3F& center1_, const VEC3F& length1_, const VEC3& center2_, const VEC3F& length2_);

//  Draw a bounding box.
void DrawCube(const float boxXLength, const float boxYLength, const float boxZLength);

#endif  // UTIL_H
