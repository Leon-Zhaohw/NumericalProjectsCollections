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

#include "kelvinlet/dynaPushBase.h"

namespace Kelvinlet 
{

// Common factor used by most methods below.
static Scalar _EvalDen(const Scalar r, const Scalar a)
{
    return 16. * M_PI * std::pow(a,2) * std::pow(r,3);
}

Scalar
DynaPushBase::_EvalU(const Scalar r, const Scalar t, const Scalar a) const
{
    // Eq. 12a in [de Goes and James 2018]

    if (a == Infty) return 0.;

    if (t == Infty) return _EvalUStatic(r, a);

    if (r < ZeroRadius) return _EvalU0(t, a);

    return _EvalW(r, t, a) / _EvalDen(r, a);
}

Scalar 
DynaPushBase::_EvalGradU(const Scalar r, const Scalar t, const Scalar a) const 
{
    // See Appendix C in [de Goes and James 2018]

    if (a == Infty) return 0.;

    if (t == Infty) return _EvalGradUStatic(r, a);

    if (r < ZeroRadius) return 0.;

    Scalar f = _EvalW(r, t, a);
    Scalar g = _EvalGradW(r, t, a);
    return (g - 3.*f/r) / _EvalDen(r, a);
}

Scalar 
DynaPushBase::_EvalHessU(const Scalar r, const Scalar t, const Scalar a) const 
{
    // See Appendix C in [de Goes and James 2018]

    if (a == Infty) return 0.;

    if (t == Infty) return _EvalHessUStatic(r, a);

    if (r < ZeroRadius) return 0.;

    Scalar r2 = r*r;
    Scalar f = _EvalW(r, t, a);
    Scalar g = _EvalGradW(r, t, a);
    Scalar h = _EvalHessW(r, t, a);
    return (h - 6.*g/r + 12.*f/r2) / _EvalDen(r, a);
}

/////////////////////////////
// Statics and Derivatives //
/////////////////////////////

Scalar
DynaPushBase::_EvalUStatic(const Scalar r, const Scalar a) const
{
    // See Suppl Material in [de Goes and James 2018]
    Scalar re = std::sqrt(r*r + _eps*_eps);
    Scalar val = 8.*M_PI*a*a;
    return 1. / (val*re); 
}

Scalar
DynaPushBase::_EvalGradUStatic(const Scalar r, const Scalar a) const
{
    Scalar re2 = r*r + _eps*_eps;
    Scalar re = std::sqrt(re2);
    Scalar val = 8.*M_PI*a*a;
    return - r / (val*re2*re); 
}

Scalar
DynaPushBase::_EvalHessUStatic(const Scalar r, const Scalar a) const
{
    Scalar e2  = std::pow(_eps,2);
    Scalar r2  = std::pow(r,2);
    Scalar re2 = r2 + e2;
    Scalar re = std::sqrt(re2);
    Scalar val = 8.*M_PI*a*a;
    return (2.*r2 - e2) / (val*re2*re2*re); 
}

//////////////////
// Limit at R=0 //
//////////////////

Scalar
DynaPushBase::_EvalU0(const Scalar t, const Scalar a) const
{
    // Eq. 13 in [de Goes and James 2018]

    Scalar e2 = std::pow(_eps, 2);
    Scalar a2t2 = std::pow(a*t, 2);
    Scalar ate = std::sqrt(a2t2 + e2);

    Scalar val0 = 2. / _eps;
    Scalar val1 = 2.*e2*e2 / std::pow(ate, 5);
    Scalar val2 = 16.*M_PI*std::pow(a,2);
    return (val0 - val1) / val2;
}

///////////////////////
// W and Derivatives //
///////////////////////

Scalar 
DynaPushBase::_EvalW(const Scalar r, const Scalar t, const Scalar a) const
{
    // Eq. 12a in [de Goes and James 2018]

    Scalar r2 = std::pow(r, 2);
    Scalar e2 = std::pow(_eps, 2);
    Scalar re = std::sqrt(r2 + e2);

    Scalar Rp = r + a*t;
    Scalar Rn = r - a*t;
    Scalar Rpe = std::sqrt(Rp*Rp + e2);
    Scalar Rne = std::sqrt(Rn*Rn + e2);

    Vector3 val;
    val[0] = 2.*r2*r/re;
    val[1] = -Rp*Rn*(Rp/Rpe + Rn/Rne);
    val[2] = a*t*e2*(1./Rpe - 1./Rne);
    return val.sum();
}

Scalar 
DynaPushBase::_EvalGradW(const Scalar r, const Scalar t, const Scalar a) const
{
    // See Appendix C in [de Goes and James 2018]

    Scalar r2 = std::pow(r, 2);
    Scalar e2 = std::pow(_eps, 2);

    Scalar re2 = r2 + e2;
    Scalar re  = std::sqrt(re2);
    Scalar re3 = re * re2;

    Scalar Rp   = r + a*t;
    Scalar Rpe2 = Rp*Rp + e2;
    Scalar Rpe  = std::sqrt(Rpe2);
    Scalar Rpe3 = Rpe * Rpe2;

    Scalar Rn   = r - a*t;
    Scalar Rne2 = Rn*Rn + e2;
    Scalar Rne  = std::sqrt(Rne2);
    Scalar Rne3 = Rne * Rne2;

    Vector3 val;
    val[0] = 2.*r2*(3./re - r2/re3);
    val[1] = -2.*r*(Rn/Rne + Rp/Rpe);
    val[2] = -e2*r*(Rn/Rne3 + Rp/Rpe3);
    return val.sum();
}

Scalar
DynaPushBase::_EvalHessW(const Scalar r, const Scalar t, const Scalar a) const
{
    // See Appendix C in [de Goes and James 2018]

    Scalar r2 = std::pow(r, 2);
    Scalar e2 = std::pow(_eps, 2);

    Scalar re2 = r2 + e2;
    Scalar re  = std::sqrt(re2);
    Scalar re3 = re  * re2;
    Scalar re5 = re3 * re2;

    Scalar Rp   = r + a*t;
    Scalar Rpe2 = Rp*Rp + e2;
    Scalar Rpe  = std::sqrt(Rpe2);
    Scalar Rpe3 = Rpe  * Rpe2;
    Scalar Rpe5 = Rpe3 * Rpe2;

    Scalar Rn   = r - a*t;
    Scalar Rne2 = Rn*Rn + e2;
    Scalar Rne  = std::sqrt(Rne2);
    Scalar Rne3 = Rne  * Rne2;
    Scalar Rne5 = Rne3 * Rne2;

    Vector3 val;
    val[0] = 2.*(2.*r/re - Rn/Rne - Rp/Rpe);
    val[1] = e2*(2.*r/re3 - Rn/Rne3 - Rp/Rpe3);
    val[2] = 3.*e2*e2*r*(2./re5 - 1./Rne5 - 1./Rpe5);
    return val.sum();
}

} // namespace Kelvinlet 
