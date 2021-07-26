///////////////////////////////////////////////////////////////////////////////////////////////////
// I. LICENSE CONDITIONS
//
// Copyright (c) 2018 by Disney-Pixar
//
// Permission is hereby granted to use this software solely for non-commercial applications
// and purposes including academic or industrial research, evaluation and not-for-profit media
// production. All other rights are retained by Pixar. For use for or in connection with
// commercial applications and purposes, including without limitation in or in connection with
// software products offered for sale or for-profit media production, please contact Pixar at
// tech-licensing@pixar.com.
//
// THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
// NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PIXAR OR ITS AFFILIATES BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H
#define SETTINGS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace CubeSim
{

using Scalar = double;

using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
using Matrix34 = Eigen::Matrix<Scalar, 3, 4>;
using Matrix38 = Eigen::Matrix<Scalar, 3, 8>;
using Matrix83 = Eigen::Matrix<Scalar, 8, 3>;
using Matrix912 = Eigen::Matrix<Scalar, 9, 12>;
using Matrix924 = Eigen::Matrix<Scalar, 9, 24>;
using Matrix9 = Eigen::Matrix<Scalar, 9, 9>;
using Matrix1212 = Eigen::Matrix<Scalar, 12, 12>;
using Matrix2424 = Eigen::Matrix<Scalar, 24, 24>;
using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
using Vector4 = Eigen::Matrix<Scalar, 4, 1>;
using Vector8 = Eigen::Matrix<Scalar, 8, 1>;
using Vector9 = Eigen::Matrix<Scalar, 9, 1>;
using Vector24 = Eigen::Matrix<Scalar, 24, 1>;

using Vector3i = Eigen::Matrix<int, 3, 1>;
using Vector4i = Eigen::Matrix<int, 4, 1>;

using Array3 = Eigen::Array<Scalar, 3, 1>;

using Matrix3X = Eigen::Matrix<Scalar, 3, Eigen::Dynamic>;
using Matrix8Xi = Eigen::Matrix<int, 8, Eigen::Dynamic>;

using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
using SparseMatrixX = Eigen::SparseMatrix<Scalar>;

constexpr auto GREEN_C = "\033[0;32m";
constexpr auto RED_C = "\033[0;31m";
constexpr auto NO_C = "\033[0m";

}

#endif
