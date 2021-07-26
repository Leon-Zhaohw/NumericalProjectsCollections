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

#include "kelvinlet/brushAffine.h"

namespace Kelvinlet 
{

void
BrushAffine::Copy(const BrushAffine& other)
{
    BrushBase::Copy(other);
    _force = other._force;
}

Vector3
BrushAffine::Eval(const Vector3& query) const
{
    Vector3 x = query - _point;
    Scalar  r = x.norm();
    if (r < ZeroRadius) return Vector3::Zero();

    Scalar re2 = r*r + _eps*_eps;
    Scalar re = std::sqrt(re2);
    Scalar re3 = re2 * re;

    Vector3 sum = Vector3::Zero();
    sum -= (_a / re3) * (1. + 1.5*_eps*_eps/re2) * (_force*x);
    sum += (_b / re3) * (_force + _force.transpose()) * x;
    sum += (_b / re3) * (_force.trace() - (3./re2)*x.dot(_force*x)) * x;
    return sum;
}

void
BrushAffine::Calibrate() 
{
    Scalar mu = GetMu();
    Scalar nu = GetNu();
    Scalar e3 = std::pow(_eps, 3);
    
    // Twist
    Vector3 q;
    q[0] = 0.5 * (_force(2,1) - _force(1,2));
    q[1] = 0.5 * (_force(0,2) - _force(2,0));
    q[2] = 0.5 * (_force(1,0) - _force(0,1));

    // Scale
    Scalar s = _force.trace() / 3.;

    // Pinch
    Matrix33 P = 0.5*(_force + _force.transpose()) - s*Matrix33::Identity();

    // See Suppl Material [de Goes and James 2017]
    q *= -0.4 * e3 / _a;
    if (nu == 0.5) s = 0.;
    else s *= 0.4 * e3 / (2.*_b-_a);
    P *= 2. * e3 / (4.*_b - 5*_a);

    // Reconstruct
    _force  = s * Matrix33::Identity();
    _force += AssembleSkewSymMatrix(q);
    _force += P;
}

void
BrushAffine::SetTwist(const Vector3& t)
{
    _force = AssembleSkewSymMatrix(t);
}

void 
BrushAffine::SetScale(const Vector3& s)
{
    _force.setZero();
    _force.diagonal() = s;
}

} // namespace Kelvinlet
