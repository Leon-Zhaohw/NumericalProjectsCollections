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

#ifndef KELVINLET_BRUSH_GRAB_CUSP_LAPLACIAN_H
#define KELVINLET_BRUSH_GRAB_CUSP_LAPLACIAN_H

#include "kelvinlet/brushGrabBase.h"

namespace Kelvinlet 
{

class BrushGrabCuspLaplacian : public BrushGrabBase
{
public:
    virtual Scalar EvalA(const Scalar r) const override
    {
        if (r < 1.e-8 * _eps)
        {
            Scalar e3 = _eps * _eps * _eps;
            return 10. * (3.*_a - 2.*_b) / e3;
        }

        Scalar e2 = _eps * _eps;
        Scalar e4 = e2 * e2;
        Scalar e6 = e4 * e2;

        Scalar r2 = r * r;
        Scalar r4 = r2 * r2;
        Scalar r6 = r4 * r2;

        Scalar re2 = r2 + e2;
        Scalar re  = std::sqrt(re2);
        Scalar re5 = re2 * re2 * re;

        Scalar A = 0.;
        A += (15.*_a - 10.*_b)*e6;
        A += (90.*_a - 88.*_b)*e4*r2;
        A += 120.*e2*r4*(_a-_b);
        A += 48.*r6*(_a-_b);
        A += 48.*r*re5*(_b-_a);
        A *= 2. / (e4 * re5);
        return A;
    }

    virtual Scalar EvalB(const Scalar r) const override
    {
        if (r < 1.e-8 * _eps)
        {
            return 0.;
        }

        Scalar e2 = _eps * _eps;
        Scalar e4 = e2 * e2;

        Scalar r2 = r * r;
        Scalar r4 = r2 * r2;

        Scalar re2 = r2 + e2;
        Scalar re  = std::sqrt(re2);
        Scalar re5 = re2 * re2 * re;

        Scalar B = 0.;
        B += 2.*e2*r2*(4.*re - 5.*r);
        B += 4.*r4*(re - r);
        B += e4*(4.*re - 7.*r);
        B *= -(12.*_b) / (e4 * r * re5);
        return B;
    }
};

} // namespace Kelvinlet

#endif // KELVINLET_BRUSH_GRAB_CUSP_LAPLACIAN_H
