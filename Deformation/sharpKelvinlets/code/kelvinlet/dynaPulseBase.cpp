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

#include "kelvinlet/dynaPulseBase.h"

namespace Kelvinlet 
{

// Common factor used by most methods below.
static Scalar _EvalDen(const Scalar r, const Scalar a)
{
    return 16. * M_PI * a * std::pow(r,3);
}

Scalar 
DynaPulseBase::_EvalU(const Scalar r, const Scalar t, const Scalar a) const
{
    // Eq. 8a in [de Goes and James 2018]

    if (a == Infty) return 0.;

    if (t == Infty) return 0.;

    if (r < ZeroRadius) return _EvalU0(t, a);

    Scalar fP = _EvalW(r, r+a*t);
    Scalar fN = _EvalW(r, r-a*t);
    return (fP - fN) / _EvalDen(r, a);
}

Scalar 
DynaPulseBase::_EvalGradU(const Scalar r, const Scalar t, const Scalar a) const
{
    // See Appendix C in [de Goes and James 2018]

    if (a == Infty) return 0.;

    if (t == Infty) return 0.;

    if (r < ZeroRadius) return 0.;

    Scalar fP = _EvalW(r, r+a*t);
    Scalar fN = _EvalW(r, r-a*t);
    Scalar gP = _EvalGradW(r, r+a*t);
    Scalar gN = _EvalGradW(r, r-a*t);
    return ((gP - gN) - 3.*(fP - fN)/r) / _EvalDen(r, a);
}

Scalar
DynaPulseBase::_EvalHessU(const Scalar r, const Scalar t, const Scalar a) const
{
    // See Appendix C in [de Goes and James 2018]

    if (a == Infty) return 0.;

    if (t == Infty) return 0.;

    if (r < ZeroRadius) return 0.;

    Scalar r2 = r*r;
    Scalar fP = _EvalW(r, r+a*t);
    Scalar fN = _EvalW(r, r-a*t);
    Scalar gP = _EvalGradW(r, r+a*t);
    Scalar gN = _EvalGradW(r, r-a*t);
    Scalar hP = _EvalHessW(r, r+a*t);
    Scalar hN = _EvalHessW(r, r-a*t);
    return ((hP - hN) - 6.*(gP-gN)/r + 12.*(fP-fN)/r2) / _EvalDen(r, a);
}

Scalar
DynaPulseBase::_EvalU0(const Scalar t, const Scalar a) const
{
    // Eq. 10 in [de Goes and James 2018]

    Scalar a2t2 = std::pow(a*t, 2);
    Scalar e2   = std::pow(_eps, 2);
    Scalar num  = 5.*t*e2*e2;
    Scalar den  = 8.*M_PI*std::sqrt(std::pow(a2t2 + e2, 7));
    return num / den;
}

Scalar 
DynaPulseBase::_EvalW(const Scalar r, const Scalar R) const
{
    // Eq. 8b in [de Goes and James 2018]

    Scalar e2  = std::pow(_eps, 2);
    Scalar R2  = R*R;
    Scalar Re2 = R2 + e2;
    Scalar Re  = std::sqrt(Re2);

    Scalar num = e2 + 2.*R2 - r*R*(3. - R2/Re2);
    Scalar den = Re;
    return num / den;
}

Scalar
DynaPulseBase::_EvalGradW(const Scalar r, const Scalar R) const
{
    // See Appendix C in [de Goes and James 2018]

    Scalar e2  = std::pow(_eps, 2);
    Scalar R2  = R*R;
    Scalar Re2 = R2 + e2;
    Scalar Re  = std::sqrt(Re2);

    Scalar num = - 3.*e2*e2*r;
    Scalar den = Re2*Re2*Re;
    return num / den;
}

Scalar
DynaPulseBase::_EvalHessW(const Scalar r, const Scalar R) const
{
    // See Appendix C in [de Goes and James 2018]

    Scalar e2  = std::pow(_eps, 2);
    Scalar R2  = R*R;
    Scalar Re2 = R2 + e2;
    Scalar Re  = std::sqrt(Re2);

    Scalar num = -3.*e2*e2*(Re2 - 5.*r*R);
    Scalar den = Re2*Re2*Re2*Re;
    return num / den;
}

} // namespace Kelvinlet 
