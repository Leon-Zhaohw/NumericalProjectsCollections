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

#ifndef KELVINLET_BRUSH_AFFINE_H
#define KELVINLET_BRUSH_AFFINE_H

#include "kelvinlet/brushBase.h"

namespace Kelvinlet 
{

class BrushAffine : public BrushBase
{
public:
    using Force = Matrix33;

    BrushAffine()
    : BrushBase()
    , _force(Matrix33::Zero())
    { }

    BrushAffine(const BrushAffine& other)
    { 
        Copy(other); 
    }

    BrushAffine& operator=(const BrushAffine& other)
    {
        Copy(other);
        return *this;
    }

    inline const Matrix33& GetForce() const { return _force; }
    inline void SetForce(const Matrix33& val) { _force = val; }
    
    // Set force to be an infinitesimal rotation with angle*axis vector \p t.
    void SetTwist(const Vector3& t);

    // Set force to be a scaling with values \p s.
    void SetScale(const Vector3& s);

    virtual void Calibrate() override;

    virtual void Copy(const BrushAffine& other);

    virtual Vector3 Eval(const Vector3& query) const override;

protected:
    /////////////
    // Members //
    /////////////

    Matrix33 _force;
};

} // namespace Kelvinlet

#endif // KELVINLET_BRUSH_AFFINE_H
