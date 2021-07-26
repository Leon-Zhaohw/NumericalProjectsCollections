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

#include "kelvinlet/dynaPulseAffine.h"

namespace Kelvinlet 
{

void
DynaPulseAffine::Copy(const DynaPulseAffine& other)
{
    DynaBase::Copy(other);
    _force = other._force;
}

Vector3 
DynaPulseAffine::_EvalDisp(
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
DynaPulseAffine::_Compute(
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
DynaPulseAffine::Calibrate()
{
    // Twist
    Vector3 q;
    q[0] = 0.5 * (_force(2,1) - _force(1,2));
    q[1] = 0.5 * (_force(0,2) - _force(2,0));
    q[2] = 0.5 * (_force(1,0) - _force(0,1));

    // Scale
    Scalar s = _force.trace() / 3.;

    // Pinch
    Matrix33 P = 0.5*(_force + _force.transpose()) - s*Matrix33::Identity();

    // At limit r->0, grad of pulse affine is:
    //   (GA/r)*F + B*(Ft + tr(F)*I)
    // = (GUa/r) * (  F + Ft + tr(F)*I)
    // + (GUb/r) * (4*F - Ft - tr(F)*I)
    //
    // Using Mathematica, we have:
    // \lim_{r->0} GUa/r 
    //    = (-21*t*e^4/(8*Pi)) * (e^2-2*a^2*t^2) / (e^2+a^2*t^2)^{11/2}
    //
    // argmax_t |_EvalDisp(_point,t)| = (e / _beta) / Sqrt[6]
    // At this time instance, we have: GUa/r ~ -1 / (10*a*e^4).
    //
    Scalar e4 = std::pow(_eps, 4);
    Scalar tFactor = - 10. * _beta * e4;

    // The actual calibration.
    q *= tFactor/5.;
    if (_alpha == Infty)
    {
        s = 0.;
        P *= tFactor/2.;
    }
    else
    {
        Scalar sFactor = - 10. * _alpha * e4;
        s *= sFactor/5.;
        P /= (2./tFactor + 3./sFactor);
    }

    // Reconstruct _force
    _force  = s * Matrix33::Identity();
    _force += AssembleSkewSymMatrix(q);
    _force += P;
}

} // namespace Kelvinlet
