///////////////////////////////////////////////////////////////////////////////////////////////////
// I. LICENSE CONDITIONS
//
// Copyright (c) 2018 by Disney-Pixar
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

#include "kelvinlet/dynaBase.h"
#include "kelvinlet/dynaPulseAffine.h"
#include "kelvinlet/dynaPulseGrab.h"
#include "kelvinlet/dynaPushAffine.h"
#include "kelvinlet/dynaPushGrab.h"

#ifdef USE_TBB
    #include <tbb/parallel_for.h>
#endif

#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

//-------//
// Types //
//-------//

using Scalar   = Kelvinlet::Scalar;
using Vector2  = Kelvinlet::Vector2;
using Vector3  = Kelvinlet::Vector3;
using Matrix33 = Kelvinlet::Matrix33;
using Points   = std::vector<Vector3>;
using Deformer = std::vector<Kelvinlet::DynaBase::Ptr>;

struct Mesh
{
    Points points;
    std::vector<int> faceSize;
    std::vector<int> faceIndices;
};

//------------------//
// Static Functions //
//------------------//

static void _PrintHelp(const char* programName)
{
    std::cout 
    << "USAGE:\n"
    << programName << "\n"
    << "-a action [grab=0 | twist=1 | scale=2 | pinch=3] \n"
    << "-p [UsePush, otherwise Impulse] \n"
    << "-t time \n"
    << "-s step \n"
    << "-o prefix \n"
    << "-h \n"
    << "\n";
}

template<class DynaType>
static void _AddDeformer(
    Deformer& deformer,
    const typename DynaType::Force& force)
{
    //-----------------------//
    // Hard-coded Parameters //
    //-----------------------//
    Scalar nu  = .45; // <= 0.5
    Scalar mu  = 5.0; // > 0.
    Scalar eps = 1.5; // > 0.

    DynaType def;
    def.SetEps(eps);
    def.SetMaterial(mu, nu);
    def.SetForce(force);
    def.Calibrate();

    Kelvinlet::DynaBase::Ptr ptr;
    ptr.reset(new DynaType(def));
    deformer.push_back(ptr); 
}

static Vector3 _Deform(
    const Vector3& p, 
    const Scalar&  t,
    const Deformer& deformer)
{
    Vector3 u = Vector3::Zero();
    for (const auto& def : deformer)
    {
        u += def->EvalDispRK4(p, t);
    }
    return u;
}

static void _SampleGrid(
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

static void _Write(
    const int suffix, 
    const Mesh& mesh,
    const char* prefix)
{
    std::stringstream filename;
    filename << prefix << "/output_" << suffix+1 << ".obj";

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

//------//
// MAIN //
//------//

int main(int argc, char** argv)
{
    // Default Input //
    int   action = 0;
    float fTime  = 3.f;
    float dt     = 0.05f;
    bool  push   = false;
    char* prefix = (char*) ".";

    // Parse Cmd Line //
    int arg;
    while ((arg = getopt(argc,argv,"o:t:s:a:ph")) != -1) 
    {
        switch (arg) 
        {
            case 'o':
                prefix = optarg;
                break;
            case 't':
                fTime = atof(optarg);
                break;
            case 's':
                dt = atof(optarg);
                break;
            case 'a':
                action = atoi(optarg);
                break;
            case 'p':
                push = true;
                break;
            default:
                _PrintHelp(argv[0]);
                return 0;
        }
    }

    // Add deformer //
    Deformer deformer;
    if (action == 0) // grab
    {
        // Hard-coded force
        Vector3 f = 2.*Vector3::UnitY();
        if (push)
        {
            using DynaType = Kelvinlet::DynaPushGrab;
            _AddDeformer<DynaType>(deformer, f);
        }
        else
        {
            using DynaType = Kelvinlet::DynaPulseGrab;
            _AddDeformer<DynaType>(deformer, f);
        }
    }
    else // affine
    {
        // Hard-coded force
        Matrix33 F = Matrix33::Zero();
        if (action == 1) // twist
        {
            Vector3 axisAngle = 0.25 * M_PI * Vector3::UnitZ();
            F = Kelvinlet::AssembleSkewSymMatrix(axisAngle);
        }
        else if (action == 2) // uniform scale
        {
            F = Matrix33::Identity();
        }
        else if (action == 3) // pinch
        {
            F(0,0) =  1.;
            F(1,1) = -1.;
        }
        if (push) 
        {
            using DynaType = Kelvinlet::DynaPushAffine;
            _AddDeformer<DynaType>(deformer, F);
        }
        else
        {
            using DynaType = Kelvinlet::DynaPulseAffine;
            _AddDeformer<DynaType>(deformer, F);
        }
    }

    // Load data //
    Vector2 gDx = Vector2::Constant(5.);
    size_t gNum = 50;

    Mesh mesh;
    _SampleGrid(mesh, -gDx, gDx, gNum, gNum);

    int nPts   = mesh.points.size();
    int nSteps = std::ceil(fTime / dt);
    std::vector<Points> results(nSteps, mesh.points);

    // Deform //
#ifdef USE_TBB
    tbb::parallel_for(0, nSteps, [&](int i)
    {
        Scalar  time = i*dt;
        Points& points = results[i];
        tbb::parallel_for(0, nPts, [&](int j)
        {
            points[j] += _Deform(points[j], time, deformer);
        });
    });
#else
    for (int i = 0; i < nSteps; ++i)
    {
        Scalar  time = i*dt;
        Points& points = results[i];
        for(int j = 0; j < nPts; ++j)
        {
            points[j] += _Deform(points[j], time, deformer);
        }
    }
#endif

    // Output //
    for (int i = 0; i < results.size(); ++i)
    {
        mesh.points = results[i];
        _Write(i, mesh, prefix);
    }
    return 1;
}
