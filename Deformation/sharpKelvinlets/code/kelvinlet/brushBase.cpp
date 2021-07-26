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

#include "kelvinlet/brushBase.h"

#include <cassert>

namespace Kelvinlet 
{

void
BrushBase::SetMaterial(const Scalar mu, const Scalar nu)
{
    assert(mu > 0.);
    assert(nu <= 0.5);
    _a = 0.25 / (M_PI * mu);
    _b = 0.25 * _a / (1.-nu);
}

Scalar 
BrushBase::GetNu() const
{
    return 1. - 0.25*(_a/_b);
}

Scalar 
BrushBase::GetMu() const
{
    return 0.25 / (M_PI * _a);
}

void
BrushBase::Copy(const BrushBase& other)
{
    _point = other._point;
    _eps = other._eps;
    _a = other._a;
    _b = other._b;
}

Vector3
BrushBase::EvalDispRK4(const Vector3& query) const
{
    Vector3 v0 = Eval(query);
    Vector3 v1 = Eval(query + 0.5*v0);
    Vector3 v2 = Eval(query + 0.5*v1);
    Vector3 v3 = Eval(query + v2);
    return (v0 + 2.*v1 + 2.*v2 + v3) / 6.;
}

Vector3
BrushBase::EvalDisp(const Vector3& query, const unsigned iters) const
{
    Vector3 p = query;
    for (unsigned i = 0; i < iters; ++i)
    {
        p += Eval(p) / iters;
    }
    return p - query;
}

} // namespace Kelvinlet
