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

#ifndef KELVINLET_BRUSH_GRAB_LAPLACIAN_H
#define KELVINLET_BRUSH_GRAB_LAPLACIAN_H

#include "kelvinlet/brushGrabBase.h"

namespace Kelvinlet 
{

class BrushGrabLaplacian : public BrushGrabBase
{
public:
    virtual Scalar EvalA(const Scalar r) const override
    {
        Scalar eps2 = _eps * _eps;
        Scalar eps4 = eps2 * eps2;
        Scalar r2   = r * r;
        Scalar re2  = r2 + eps2;
        Scalar re   = std::sqrt(re2);
        Scalar re7  = std::pow(re, 7);
        return (_b*re2*(5.*eps2 + 2.*r2) - 7.5*_a*eps4) / re7;
    }

    virtual Scalar EvalB(const Scalar r) const override
    {
        Scalar eps2 = _eps * _eps;
        Scalar r2   = r * r;
        Scalar re2  = r2 + eps2;
        Scalar re   = std::sqrt(re2);
        Scalar re7  = std::pow(re, 7);
        return -3.*_b*(7.*eps2 + 2.*r2) / re7;
    }
};

} // namespace Kelvinlet

#endif // KELVINLET_BRUSH_GRAB_LAPLACIAN_H
