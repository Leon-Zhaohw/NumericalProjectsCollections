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

#ifndef KELVINLET_BRUSH_GRAB_CUSP_H
#define KELVINLET_BRUSH_GRAB_CUSP_H

#include "kelvinlet/brushGrabBase.h"

namespace Kelvinlet 
{

class BrushGrabCusp : public BrushGrabBase
{
public:
    virtual Scalar EvalA(const Scalar r) const override
    {
        if (r < 1.e-8 * _eps)
        {
            return (3.*_a - 2.*_b) / _eps;
        }

        Scalar r2 = r * r;
        Scalar e2 = _eps * _eps;
        Scalar re2 = r2 + e2;
        Scalar re  = std::sqrt(re2);

        return (_a/e2) * (r2/re + 3.*re - 4.*r) + (2.*_b/e2)*(r - re);
    }

    virtual Scalar EvalB(const Scalar r) const override
    {
        if (r < 1.e-8 * _eps)
        {
            return 0.;
        }

        Scalar r2 = r * r;
        Scalar e2 = _eps * _eps;
        Scalar re2 = r2 + e2;
        Scalar re  = std::sqrt(re2);

        return (2.*_b/e2) * (1./r - 1./re);
    }
};

} // namespace Kelvinlet

#endif // KELVINLET_BRUSH_GRAB_CUSP_H
