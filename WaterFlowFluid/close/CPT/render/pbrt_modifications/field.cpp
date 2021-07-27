
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


// textures/field.cpp*
#include "stdafx.h"
#include "textures/field.h"

#include <iostream>
using namespace std;

// CheckerboardTexture Method Definitions
Texture<float> *CreateFieldFloatTexture(const Transform &tex2world,
        const TextureParams &tp) {
    Reference<Texture<float> > tex1 = tp.GetFloatTexture("tex1", 1.f);
    Reference<Texture<float> > tex2 = tp.GetFloatTexture("tex2", 0.f);

    string filename = tp.FindString("filename", "");
    TextureMapping3D *map = new IdentityMapping3D(tex2world);

    float amplitude = tp.FindFloat("amplitude", 1.0);
    return new Field3DTexture<float>(map, tex1, tex2, filename, amplitude);
}



Texture<Spectrum> *CreateFieldSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp) {
    Reference<Texture<Spectrum> > tex1 = tp.GetSpectrumTexture("tex1", 1.f);
    Reference<Texture<Spectrum> > tex2 = tp.GetSpectrumTexture("tex2", 0.f);

    string filename = tp.FindString("filename", "");
    TextureMapping3D *map = new IdentityMapping3D(tex2world);
    float amplitude = tp.FindFloat("amplitude", 1.0);
    return new Field3DTexture<Spectrum>(map, tex1, tex2, filename, amplitude);
}


