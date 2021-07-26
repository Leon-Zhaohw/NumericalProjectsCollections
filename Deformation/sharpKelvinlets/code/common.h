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

#ifndef COMMON_H
#define COMMON_H

#include "kelvinlet/types.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using Scalar   = Kelvinlet::Scalar;
using Vector2  = Kelvinlet::Vector2;
using Vector3  = Kelvinlet::Vector3;
using Matrix33 = Kelvinlet::Matrix33;
using Points   = std::vector<Vector3>;

struct Mesh
{
    Points points;
    std::vector<int> faceSize;
    std::vector<int> faceIndices;
};

void SampleGrid(
    Mesh& mesh,
    const Vector2& bot,
    const Vector2& top,
    const size_t nx,
    const size_t ny)
{
    Vector2 diag = top - bot;
    Scalar dx = diag[0] / nx;
    Scalar dy = diag[1] / ny;

    for (size_t j = 0; j <= ny; ++j)
    for (size_t i = 0; i <= nx; ++i)
    {
        Vector2 p = bot + Vector2(i*dx, j*dy);
        Vector3 q(p[0], p[1], 0.);
        mesh.points.push_back(q);
    }

    for (size_t j = 1; j <= ny; ++j)
    for (size_t i = 1; i <= nx; ++i)
    {
        mesh.faceSize.push_back(4);
        mesh.faceIndices.push_back((j-1)*(nx+1) + i-1);
        mesh.faceIndices.push_back((j-1)*(nx+1) + i  );
        mesh.faceIndices.push_back(    j*(nx+1) + i  );
        mesh.faceIndices.push_back(    j*(nx+1) + i-1);
    }
}

void Write(
    const int suffix, 
    const Mesh& mesh,
    const char* prefix)
{
    std::stringstream filename;
    filename << prefix << "/output_" << (suffix+1) << ".obj";

    std::ofstream out(filename.str().c_str());
    for (const auto& p : mesh.points)
    {
        out << "v " << p.transpose() << "\n";
    }
    int counter = 0;
    for (int i = 0; i < mesh.faceSize.size(); ++i)
    {
        out << "f ";
        for (int j = 0; j < mesh.faceSize[i]; ++j)
        {
            out << mesh.faceIndices[counter++] + 1 << " ";
        }
        out << "\n";
    }
    out.close();

    std::cout << "Saved: " << filename.str() << "\n";
}

#endif // COMMON_H
