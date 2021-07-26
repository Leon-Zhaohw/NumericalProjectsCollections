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

#include "kelvinlet/dynaPushGrab.h"

#include <cassert>

namespace Kelvinlet 
{

void
DynaPushGrab::Copy(const DynaPushGrab& other)
{
    DynaBase::Copy(other);
    _force = other._force;
}

Vector3 
DynaPushGrab::_EvalDisp(
    const Vector3& query, 
    const Scalar time) const 
{
    Scalar t = time - _time;
    if (t <= 0.) return Vector3::Zero();

    Vector3 x = query - _point;
    Scalar  r = x.norm();
    
    Cache values;
    _Compute(values, r, t);

    Scalar A = _EvalA(values);
    Scalar B = _EvalB(values);
    return A * _force + B * _force.dot(x) * x;
}

void 
DynaPushGrab::_Compute(
    Cache& values, 
    const Scalar r, 
    const Scalar t) const
{
    values.r = r;

    values.Ua[0] = _EvalU(r, t, _alpha);
    values.Ua[1] = _EvalGradU(r, t, _alpha);

    values.Ub[0] = _EvalU(r, t, _beta);
    values.Ub[1] = _EvalGradU(r, t, _beta);
}

void
DynaPushGrab::Calibrate() 
{
    // Eq. 7 in [de Goes and James 2017]
    Scalar mu = GetMu();
    Scalar nu = GetNu();
    Scalar a = 1./(4.*M_PI*mu);
    Scalar b = a/(4.*(1.-nu));
    _force *= 2.*_eps / (3.*a - 2.*b);
}

} // namespace Kelvinlet
