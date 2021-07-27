
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_TEXTURES_FIELD_H
#define PBRT_TEXTURES_FIELD_H

// textures/checkerboard.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
#include "montecarlo.h"
#include "shape.h"
#include "parallel.h"
#include "progressreporter.h"
#include "FIELD_3D.h"

template <typename T> class Field3DTexture : public Texture<T> {
public:
    // Field3DTexture Public Methods
    Field3DTexture(TextureMapping3D *m, Reference<Texture<T> > c1,
                          Reference<Texture<T> > c2, string filename, float amplitude)
        : mapping(m), tex1(c1), tex2(c2), _field(NULL), _filename(filename), _amplitude(amplitude) {
      if (_filename.size() > 0)
        _field = new FIELD_3D(_filename.c_str());
    }
    ~Field3DTexture() {
        delete mapping;
        if (_field) delete _field;
    }
    T Evaluate(const DifferentialGeometry &dg) const {
        Vector dpdx, dpdy;
        Point p = mapping->Map(dg, &dpdx, &dpdy);

        VEC3F point(p.x, p.y, p.z);
        float value = (*_field).quarticLookup(point);

        return _amplitude * value;
    }

private:
    // Field3DTexture Private Data
    TextureMapping3D *mapping;
    Reference<Texture<T> > tex1, tex2;

    FIELD_3D* _field;
    string _filename;

    float _amplitude;
};


Texture<float> *CreateFieldFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
Texture<Spectrum> *CreateFieldSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_FIELD_H
