/* 
 * File: definition.hpp
 *
 * Copyright (C) 2019  Camille Schreck
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

#include <Eigen/Core>

using namespace Eigen;

//#define DOUBLE_PRECISION

#ifdef DOUBLE_PRECISION

#define FLOAT double
#define VEC2 Vector2d
#define VEC3 Vector3d
#define VEC4 Vector4d
#define VECX VectorXd
#define MAT2 Matrix2d
#define MAT3 Matrix3d
#define MAT4 Matrix4d
#define MATX MatrixXd
#define ANGLE_AXIS AngleAxisd
#define QUATERNION Quaterniond
#define COMPLEX std::complex<double>
#define VEC2C Vector2cd
#define VECXC VectorXcd
#define MATXC MatrixXcd

#else

#define FLOAT float
#define VEC2 Vector2f
#define VEC3 Vector3f
#define VEC4 Vector4f
#define VECX VectorXf
#define MAT2 Matrix2f
#define MAT3 Matrix3f
#define MAT4 Matrix4f
#define MATX MatrixXf
#define ANGLE_AXIS AngleAxisf
#define QUATERNION Quaternionf
#define COMPLEX std::complex<float>
#define VEC2C Vector2cf
#define VECXC VectorXcf
#define MATXC MatrixXcf

#endif

namespace definitions {
  extern COMPLEX i_;
  extern COMPLEX r_;
};


#endif
