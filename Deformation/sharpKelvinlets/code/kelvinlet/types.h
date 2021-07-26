///////////////////////////////////////////////////////////////////////////////////////////////////
// I. LICENSE CONDITIONS
//
// Copyright (c) 2019 by Disney-Pixar
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

#ifndef KELVINLET_TYPES_H
#define KELVINLET_TYPES_H

#include <Eigen/Dense>

#include <limits>

namespace Kelvinlet
{

using Scalar = double;

template <int r, int c> using MatrixType = Eigen::Matrix<Scalar, r, c>;
template <int d>        using VectorType = MatrixType<d, 1>;

using Vector2  = VectorType<2>;
using Vector3  = VectorType<3>;
using Matrix33 = MatrixType<3, 3>;

const Scalar ZeroRadius = 1.e-4;
const Scalar Infty = std::numeric_limits<Scalar>::max();

inline Matrix33 AssembleSkewSymMatrix(const Vector3& x)
{
    Matrix33 A = Matrix33::Zero();
    /*A(0,0) =  0.;*/  A(0,1) = -x[2];  A(0,2) =  x[1];
      A(1,0) =  x[2];/*A(1,1) =  0.;*/  A(1,2) = -x[0];
      A(2,0) = -x[1];  A(2,1) =  x[0];/*A(2,2) = 0.;*/
    return A;
}

} // namespace Kelvinlet

#endif // KELVINLET_TYPES_H
