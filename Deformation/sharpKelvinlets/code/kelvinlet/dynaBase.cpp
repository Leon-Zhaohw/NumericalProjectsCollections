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

#include "kelvinlet/dynaBase.h"

#include <cassert>

namespace Kelvinlet 
{

void
DynaBase::SetMaterial(Scalar mu, Scalar nu)
{
    assert(mu > 0.);
    assert(nu <= 0.5);
    _beta = std::sqrt(mu);
    if (nu == 0.5) _alpha = Infty;
    else _alpha = _beta * std::sqrt(1. + 1./(1.-2.*nu));
}

Scalar 
DynaBase::GetNu() const
{
    if (_alpha == Infty) return 0.5;
    if (_alpha == _beta) return -Infty;
    assert(_alpha > _beta);
    Scalar c = std::pow(_alpha / _beta, 2);
    return 0.5*(1. + 1./(1.-c));
}

Scalar 
DynaBase::GetMu() const
{
    return _beta*_beta;
}

void
DynaBase::Copy(const DynaBase& other)
{
    _point = other._point;
    _alpha = other._alpha;
    _beta  = other._beta;
    _time  = other._time;
    _eps   = other._eps;
}

Scalar
DynaBase::_EvalA(const Cache& values) const
{
    // Eq. 9a in [de Goes and James 2018]
    Scalar r   = values.r;
    Scalar Ua  = values.Ua[0];
    Scalar Ub  = values.Ub[0];
    Scalar GUb = values.Ub[1];
    return Ua + 2.*Ub + r*GUb;
}

Scalar 
DynaBase::_EvalB(const Cache& values) const
{
    // Eq. 9b in [de Goes and James 2018]
    Scalar r = values.r;
    if (r < ZeroRadius) return 0.;
    
    Scalar GUa = values.Ua[1]; 
    Scalar GUb = values.Ub[1]; 
    return (GUa - GUb) / r;
}

Scalar
DynaBase::_EvalGradA(const Cache& values) const
{
    // See Appendix C in [de Goes and James 2018]
    Scalar r   = values.r;
    Scalar GUa = values.Ua[1]; 
    Scalar GUb = values.Ub[1]; 
    Scalar HUb = values.Ub[2]; 
    return GUa + 3.*GUb + r*HUb;
}

Scalar
DynaBase::_EvalGradB(const Cache& values) const
{
    // See Appendix C in [de Goes and James 2018]
    Scalar r   = values.r;
    Scalar HUa = values.Ua[2]; 
    Scalar HUb = values.Ub[2]; 
    Scalar B   = _EvalB(values);
    return (HUa - HUb - B) / r;
}

Vector3
DynaBase::EvalDispRK4(const Vector3& query, const Scalar time) const
{
    Vector3 v0 = _EvalDisp(query, time);
    Vector3 v1 = _EvalDisp(query + 0.5*v0, time);
    Vector3 v2 = _EvalDisp(query + 0.5*v1, time);
    Vector3 v3 = _EvalDisp(query + v2, time);
    return (v0 + 2.*v1 + 2.*v2 + v3) / 6.;
}

Vector3
DynaBase::EvalDisp(const Vector3& query, const Scalar time, const unsigned iters) const
{
    Vector3 p = query;
    for (unsigned i = 0; i < iters; ++i)
    {
        p += _EvalDisp(p, time) / iters;
    }
    return p - query;
}

} // namespace Kelvinlet
