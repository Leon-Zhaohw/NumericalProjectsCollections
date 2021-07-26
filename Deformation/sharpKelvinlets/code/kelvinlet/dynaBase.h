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

#ifndef KELVINLET_DYNA_BASE_H
#define KELVINLET_DYNA_BASE_H

#include "kelvinlet/types.h"

#include <memory>

/// 
/// \class DynaBase
/// \brief Abstract base class for Dynamic Kelvinlet.
///
/// It contains the common routines and formulae to evaluate the
/// time-dependent displacement in terms of the Pseudo-Potential
/// and its derivatives.
/// 
/// Abstract Methods:
///    void    Calibrate();
///    Vector3 _EvalDisp(const Vector3& query, const Scalar time) const;
///    void    _Compute(Cache& values, const Scalar r, const Scalar t) const;
///    Scalar  _EvalU(const Scalar r, const Scalar t, const Scalar a) const;
///    Scalar  _EvalGradU(const Scalar r, const Scalar t, const Scalar a) const;
///    Scalar  _EvalHessU(const Scalar r, const Scalar t, const Scalar a) const;
///

namespace Kelvinlet 
{

class DynaBase
{
public:
    using Ptr = std::shared_ptr<DynaBase>;
    
    /////////////////////
    // Internal Struct //
    /////////////////////
    
    struct Cache
    {   
        Cache()
        : Ua(Vector3::Zero())
        , Ub(Vector3::Zero())
        , r(0.)
        { }

        Cache(const Cache& other)
        : Ua(other.Ua)
        , Ub(other.Ub)
        , r(other.r)
        { }

        Cache& operator=(const Cache& other)
        {
            Ua = other.Ua;
            Ub = other.Ub;
            r  = other.r;
            return *this;
        }

        // Convention:
        // [0] - function
        // [1] - gradient
        // [2] - hessian 
        Vector3 Ua;
        Vector3 Ub;
        Scalar  r;
    };

    //////////
    // Ctor //
    //////////

    DynaBase()
    : _point(Vector3::Zero())
    , _time(0.)
    , _eps(1.)
    , _alpha(0.)
    , _beta(0.)
    { }

    DynaBase(const DynaBase& other)
    { 
        Copy(other); 
    }

    DynaBase& operator=(const DynaBase& other)
    {
        Copy(other);
        return *this;
    }

    virtual void Copy(const DynaBase& other);

    /////////////////////
    // Get/Set Methods //
    /////////////////////

    inline void SetPoint(const Vector3& val) { _point = val;  }
    inline const Vector3& GetPoint() const   { return _point; }

    inline void SetEps(const Scalar val) { _eps = val;  }
    inline Scalar GetEps() const         { return _eps; }

    inline void SetTime(const Scalar val) { _time = val;  }
    inline Scalar GetTime() const         { return _time; }

    /// Set wave speeds via the material parameters
    /// Stiffness mu > 0 and Poisson ratio nu <= 0.5.
    void SetMaterial(Scalar mu, Scalar nu);

    Scalar GetNu() const;
    Scalar GetMu() const;

    inline Scalar GetAlpha() const { return _alpha; }
    inline Scalar GetBeta()  const { return _beta;  }

    //////////////////
    // Eval Methods //
    //////////////////

    /// Returns displacement at \p query and \p time advected by 4th order Runge-Kutta.
    Vector3 EvalDispRK4(const Vector3& query, const Scalar time) const;

    /// Returns displacemetn at \p query and \p time advected by \p iters uniform steps.
    Vector3 EvalDisp(const Vector3& query, const Scalar time, const unsigned iters = 1) const;

protected:
    // Return radial (A) / bulging (B) values and their gradients.
    Scalar _EvalA(const Cache& values) const;
    Scalar _EvalB(const Cache& values) const;
    Scalar _EvalGradA(const Cache& values) const;
    Scalar _EvalGradB(const Cache& values) const;

    //////////////////////
    // Abstract Methods //
    //////////////////////

    /// Return disp at \p query and \p time.
    virtual Vector3 _EvalDisp(const Vector3& query, const Scalar time) const = 0;

    /// Cache \p values at radius \p r and time \p t.
    virtual void _Compute(Cache& values, const Scalar r, const Scalar t) const = 0;

    /// Eval Pseudo-Potential and Derivatives at radius \p r, time \p t, and speed \p a.
    virtual Scalar _EvalU(const Scalar r, const Scalar t, const Scalar a) const = 0;
    virtual Scalar _EvalGradU(const Scalar r, const Scalar t, const Scalar a) const = 0;
    virtual Scalar _EvalHessU(const Scalar r, const Scalar t, const Scalar a) const = 0;

public:
    //////////////////////
    // Abstract Methods //
    //////////////////////

    /// Adjust _force assuming its current value is the tip deformation.
    /// This is equivalenet to solve a point constraint at the brush center.
    virtual void Calibrate() = 0;

protected:
    /////////////
    // Members //
    /////////////

    Vector3 _point; // source location
    Scalar  _time;  // activation time
    Scalar  _eps;   // regularizer

    Scalar  _alpha; // P-wave speed
    Scalar  _beta;  // S-wave speed
};

} // namespace Kelvinlet

#endif // KELVINLET_DYNA_BASE_H
