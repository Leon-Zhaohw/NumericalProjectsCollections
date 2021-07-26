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

#ifndef KELVINLET_BRUSH_GRAB_TRISCALE_H
#define KELVINLET_BRUSH_GRAB_TRISCALE_H

#include "kelvinlet/brushGrab.h"

namespace Kelvinlet 
{

class BrushGrabTriScale : public BrushGrabBase
{
public:
    virtual Scalar EvalA(const Scalar r) const override
    {
        Vector3 w, A;
        Scalar x = 1.1;
        ComputeWeights(w, x);
        A[0] = BrushGrab::EvalRadial(r, _eps, _a, _b);
        A[1] = BrushGrab::EvalRadial(r, x*_eps, _a, _b);
        A[2] = BrushGrab::EvalRadial(r, x*x*_eps, _a, _b);
        return (w.array() * A.array()).sum();
    }

    virtual Scalar EvalB(const Scalar r) const override
    {
        Vector3 w, B;
        Scalar x = 1.1;
        ComputeWeights(w, x);
        B[0] = BrushGrab::EvalBulge(r, _eps, _b);
        B[1] = BrushGrab::EvalBulge(r, x*_eps, _b);
        B[2] = BrushGrab::EvalBulge(r, x*x*_eps, _b);
        return (w.array() * B.array()).sum();
    }

    void ComputeWeights(Vector3& w, const Scalar x) const
    {
        Scalar x2 = x*x;
        Scalar x4 = x2*x2;

        w[0] = 1.;
        w[1] = (1. - x4) / (x4 - x2);
        w[2] = (x2 - 1.) / (x4 - x2);
    }
};

} // namespace Kelvinlet

#endif // KELVINLET_BRUSH_GRAB_TRISCALE_H
