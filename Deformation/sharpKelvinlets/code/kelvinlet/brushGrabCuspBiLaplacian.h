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

#ifndef KELVINLET_BRUSH_GRAB_CUSP_BILAPLACIAN_H
#define KELVINLET_BRUSH_GRAB_CUSP_BILAPLACIAN_H

#include "kelvinlet/brushGrabBase.h"

namespace Kelvinlet 
{

class BrushGrabCuspBiLaplacian : public BrushGrabBase
{
public:
    virtual Scalar EvalA(const Scalar r) const override
    {
        Scalar R = r / _eps;
        if (R >= 5.) return 0.;

        Scalar eps5 = std::pow(_eps, 5);

        Scalar sqR = R * R;
        Scalar R1 = std::sqrt(sqR + 1.);
        Scalar R1_3  = std::pow(R1,  3);
        Scalar R1_9  = std::pow(R1_3,3);
        Scalar R1_10 = R1_9 * R1;

        Scalar Iterm0 = 512.*sqR + 2304.;
        Iterm0 = Iterm0*sqR + 4032.;
        Iterm0 = Iterm0*sqR + 3360.;
        Iterm0 = Iterm0*sqR + 1260.;
        Iterm0 = Iterm0*sqR + 105.;
        Iterm0 = Iterm0*R1 - 512*R*R1_10;
        Iterm0 = Iterm0*_a;

        Scalar Iterm1 = 448. + 128.*sqR;
        Iterm1 = Iterm1*sqR + 560.;
        Iterm1 = Iterm1*sqR + 280.;
        Iterm1 = Iterm1*sqR + 35.;
        Iterm1 = 128.*R*R1_10 - Iterm1*R1_3;
        Iterm1 = Iterm1 * 2.*_b;

        return (Iterm0 + Iterm1) * 9. / (eps5*R1_10);
    }

    virtual Scalar EvalB(const Scalar r) const override
    {
        Scalar R = r / _eps;
        if (R >= 5.) return 0.;

        Scalar sqR = R * R;
        Scalar R1 = std::sqrt(sqR + 1.);
        Scalar R1_3  = std::pow(R1,  3);
        Scalar R1_9  = std::pow(R1_3,3);

        Scalar eps5 = std::pow(_eps, 5);
        Scalar eps7 = eps5 * _eps * _eps;

        Scalar Rterm = 576. + 128.*sqR;
        Rterm = Rterm*sqR + 1008.;
        Rterm = Rterm*sqR + 840.;
        Rterm = Rterm*sqR + 315.;
        Rterm = 128.*R1_9 - Rterm*R;
        return Rterm * (18.*_b) / (eps7*R*R1_9);
    }
};

} // namespace Kelvinlet

#endif // KELVINLET_BRUSH_GRAB_CUSP_BILAPLACIAN_H
