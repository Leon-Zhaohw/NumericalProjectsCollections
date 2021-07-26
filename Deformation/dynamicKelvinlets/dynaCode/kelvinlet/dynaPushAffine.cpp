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

#include "kelvinlet/dynaPushAffine.h"

namespace Kelvinlet 
{

void
DynaPushAffine::Copy(const DynaPushAffine& other)
{
    DynaBase::Copy(other);
    _force = other._force;
}

Vector3 
DynaPushAffine::_EvalDisp(
    const Vector3& query, 
    const Scalar   time) const 
{
    Scalar t = time - _time;
    if (t <= 0.) return Vector3::Zero();

    Vector3 x = query - _point;
    Scalar  r = x.norm();
    if (r < ZeroRadius) return Vector3::Zero();

    Cache values;
    _Compute(values, r, t);
    Vector3 Fx = _force * x;
    
    Vector3 result = Vector3::Zero();
    result += _EvalGradA(values) * (Fx / r);
    result += _EvalGradB(values) * x.dot(Fx) * (x / r);
    result += _EvalB(values) * (_force.transpose()*x + _force.trace()*x);
    return result;
}

void 
DynaPushAffine::_Compute(
    Cache& values, 
    const Scalar r, 
    const Scalar t) const
{
    values.r = r;

    values.Ua[1] = _EvalGradU(r, t, _alpha);
    values.Ua[2] = _EvalHessU(r, t, _alpha);
    
    values.Ub[1] = _EvalGradU(r, t, _beta);
    values.Ub[2] = _EvalHessU(r, t, _beta);
}

void
DynaPushAffine::Calibrate()
{
    Scalar mu = GetMu();
    Scalar nu = GetNu();
    Scalar a = 1./(4.*M_PI*mu);
    Scalar b = a/(4.*(1.-nu));
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
    q *= -0.4 * e3 / a;
    if (nu == 0.5) s = 0.;
    else s *= 0.4 * e3 / (2.*b-a);
    P *= 2. * e3 / (4.*b - 5*a);

    // Reconstruct
    _force  = s * Matrix33::Identity();
    _force += AssembleSkewSymMatrix(q);
    _force += P;
}

} // namespace Kelvinlet
