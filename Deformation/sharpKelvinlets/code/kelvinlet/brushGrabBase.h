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

#ifndef KELVINLET_BRUSH_GRAB_BASE_H
#define KELVINLET_BRUSH_GRAB_BASE_H

#include "kelvinlet/brushBase.h"

/// 
/// \class BrushGrabBase
/// \brief Abstract class for a grab-like Kelvinlet deformer.
///
/// Abstract Methods:
///    Scalar  EvalA(const Scalar r) const;
///    Scalar  EvalB(const Scalar r) const;
///

namespace Kelvinlet 
{

class BrushGrabBase : public BrushBase
{
public:
    using Force = Vector3;
    
    BrushGrabBase()
    : BrushBase()
    , _force(Vector3::Zero())
    { }

    BrushGrabBase(const BrushGrabBase& other)
    { 
        Copy(other); 
    }

    BrushGrabBase& operator=(const BrushGrabBase& other)
    {
        Copy(other);
        return *this;
    }

    inline const Vector3& GetForce() const { return _force; }
    inline void SetForce(const Vector3& val) { _force = val; }

    virtual void Calibrate() override;

    virtual void Copy(const BrushGrabBase& other);

    virtual Vector3 Eval(const Vector3& query) const override;

    //////////////////////
    // Abstract Methods //
    //////////////////////

    // Return radial amount.
    virtual Scalar EvalA(const Scalar r) const = 0;

    // Return bulding amount.
    virtual Scalar EvalB(const Scalar r) const = 0;

protected:
    /////////////
    // Members //
    /////////////

    Vector3 _force;
};

} // namespace Kelvinlet

#endif // KELVINLET_BRUSH_GRAB_BASE_H
