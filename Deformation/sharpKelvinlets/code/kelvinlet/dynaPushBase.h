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

#ifndef KELVINLET_DYNA_PUSH_BASE_H
#define KELVINLET_DYNA_PUSH_BASE_H

#include "kelvinlet/dynaBase.h"

///
/// \class DynaPusjBase
/// \brief Base class for a Push Dynamic Kelvinlet.
///
/// It corresponds to the response to a continuous time series of impulses,
/// ie, DynaPush(q,t) = \int_0^t DynaPulse(q,\tau) d\tau.
///
/// At time limit (ie, t = Infty), it reduces to 3D Regularized Kelvinlets
/// presented in [de Goes and James 2017].
/// 

namespace Kelvinlet 
{

class DynaPushBase : public DynaBase
{
protected:
    //////////////////////
    // Override Methods //
    //////////////////////

    virtual Scalar _EvalU(const Scalar r, const Scalar t, const Scalar a) const override;
    virtual Scalar _EvalGradU(const Scalar r, const Scalar t, const Scalar a) const override;
    virtual Scalar _EvalHessU(const Scalar r, const Scalar t, const Scalar a) const override;

    ////////////////////
    // Limit Response //
    ////////////////////

    Scalar _EvalU0(const Scalar t, const Scalar a) const;
    Scalar _EvalUStatic(const Scalar r, const Scalar a) const;
    Scalar _EvalGradUStatic(const Scalar r, const Scalar a) const;
    Scalar _EvalHessUStatic(const Scalar r, const Scalar a) const;

    ////////////////////
    // Push Potential //
    ////////////////////

    Scalar _EvalW(const Scalar r, const Scalar t, const Scalar a) const;
    Scalar _EvalGradW(const Scalar r, const Scalar t, const Scalar a) const;
    Scalar _EvalHessW(const Scalar r, const Scalar t, const Scalar a) const;
};

} // namespace Kelvinlet 

#endif // KELVINLET_DYNA_PUSH_BASE_H
