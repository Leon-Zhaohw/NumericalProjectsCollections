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

#ifndef KELVINLET_BRUSH_BASE_H
#define KELVINLET_BRUSH_BASE_H

#include "kelvinlet/types.h"

#include <memory>

/// 
/// \class BrushBase
/// \brief Abstract class for a Kelvinlet deformer.
///
/// Abstract Methods:
///    void    Calibrate();
///    Vector3 Eval(const Vector3& query) const;
///

namespace Kelvinlet 
{

class BrushBase
{
public:
    using Ptr = std::shared_ptr<BrushBase>;
    
    //////////
    // Ctor //
    //////////

    BrushBase()
    : _point(Vector3::Zero())
    , _eps(1.)
    , _a(0.)
    , _b(0.)
    { }

    BrushBase(const BrushBase& other)
    { 
        Copy(other); 
    }

    BrushBase& operator=(const BrushBase& other)
    {
        Copy(other);
        return *this;
    }

    virtual void Copy(const BrushBase& other);

    /////////////////////
    // Get/Set Methods //
    /////////////////////

    inline void SetPoint(const Vector3& val) { _point = val;  }
    inline const Vector3& GetPoint() const   { return _point; }

    inline void SetEps(const Scalar val) { _eps = val;  }
    inline Scalar GetEps() const         { return _eps; }

    /// Stiffness mu > 0 and Poisson ratio nu <= 0.5.
    void SetMaterial(const Scalar mu, const Scalar nu);

    Scalar GetNu() const;
    Scalar GetMu() const;

    //////////////////
    // Eval Methods //
    //////////////////

    /// Returns displacement at \p query advected by 4th order Runge-Kutta.
    Vector3 EvalDispRK4(const Vector3& query) const;

    /// Returns displacement at \p query advected by \p iters uniform steps.
    Vector3 EvalDisp(const Vector3& query, const unsigned iters = 1) const;

    //////////////////////
    // Abstract Methods //
    //////////////////////

    // Return displacement at \p query.
    virtual Vector3 Eval(const Vector3& query) const = 0;

    /// Adjust force so that its current value is the tip deformation.
    /// This is equivalenet to solve a point constraint at brush center.
    virtual void Calibrate() = 0;

protected:
    /////////////
    // Members //
    /////////////

    Vector3 _point; // brush center
    Scalar  _eps;   // brush scale

    Scalar  _a; // A amount
    Scalar  _b; // B amount
};

} // namespace Kelvinlet

#endif // KELVINLET_BRUSH_BASE_H
