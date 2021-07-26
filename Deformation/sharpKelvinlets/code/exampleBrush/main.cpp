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

#include "common.h"

#include "kelvinlet/brushAffine.h"
#include "kelvinlet/brushGrab.h"
#include "kelvinlet/brushGrabBiLaplacian.h"
#include "kelvinlet/brushGrabBiScale.h"
#include "kelvinlet/brushGrabCusp.h"
#include "kelvinlet/brushGrabCuspBiLaplacian.h"
#include "kelvinlet/brushGrabCuspLaplacian.h"
#include "kelvinlet/brushGrabLaplacian.h"
#include "kelvinlet/brushGrabTriScale.h"

#ifdef USE_TBB
    #include <tbb/parallel_for.h>
#endif

#include <getopt.h>

using Deformer = std::vector<Kelvinlet::BrushBase::Ptr>;

static void _PrintHelp(const char* programName)
{
    std::cout 
    << "USAGE:\n"
    << programName << "\n"
    << "-a action \n"
    << "-o prefix \n"
    << "-h \n"
    << "\n"
    << "ACTION CAN BE\n"
    << "...0 - BrushGrab\n"
    << "...1 - BrushGrabBiScale\n"
    << "...2 - BrushGrabTriScale\n"
    << "...3 - BrushGrabLaplacian\n"
    << "...4 - BrushGrabBiLaplacian\n"
    << "...5 - BrushGrabCusp\n"
    << "...6 - BrushGrabCuspLaplacian\n"
    << "...7 - BrushGrabCuspBiLaplacian\n"
    << "...8 - BrushAffine (Twist)\n"
    << "...9 - BrushAffine (Scale)\n"
    << "..10 - BrushAffine (Pinch)\n";
}

template<class BrushType>
static void _AddDeformer(
    Deformer& deformer,
    const typename BrushType::Force& force)
{
    //-----------------------//
    // Hard-coded Parameters //
    //-----------------------//
    Scalar nu  = .45; // <= 0.5
    Scalar mu  = 5.0; // > 0.
    Scalar eps = 4.0; // > 0.

    BrushType def;
    def.SetEps(eps);
    def.SetMaterial(mu, nu);
    def.SetForce(force);
    def.Calibrate();

    Kelvinlet::BrushBase::Ptr ptr;
    ptr.reset(new BrushType(def));
    deformer.push_back(ptr); 
}

static Vector3 _Deform(
    const Vector3& p, 
    const Deformer& deformer)
{
    Vector3 u = Vector3::Zero();
    for (const auto& def : deformer)
    {
        u += def->Eval(p);
    }
    return u;
}

int main(int argc, char** argv)
{
    // Default Input //
    int   action = 0;
    char* prefix = (char*) ".";

    // Parse Cmd Line //
    int arg;
    while ((arg = getopt(argc,argv,"o:a:h")) != -1) 
    {
        switch (arg) 
        {
            case 'o':
                prefix = optarg;
                break;
            case 'a':
                action = atoi(optarg);
                break;
            default:
                _PrintHelp(argv[0]);
                return 0;
        }
    }

    // Add deformer //
    Deformer deformer;

    if (action == 0)
    {
        Vector3 f = Vector3::UnitY();
        using BrushType = Kelvinlet::BrushGrab;
        _AddDeformer<BrushType>(deformer, f);
    }
    else if (action == 1)
    {
        Vector3 f = Vector3::UnitY();
        using BrushType = Kelvinlet::BrushGrabBiScale;
        _AddDeformer<BrushType>(deformer, f);
    }
    else if (action == 2)
    {
        Vector3 f = Vector3::UnitY();
        using BrushType = Kelvinlet::BrushGrabTriScale;
        _AddDeformer<BrushType>(deformer, f);
    }
    else if (action == 3)
    {
        Vector3 f = Vector3::UnitY();
        using BrushType = Kelvinlet::BrushGrabLaplacian;
        _AddDeformer<BrushType>(deformer, f);
    }
    else if (action == 4)
    {
        Vector3 f = Vector3::UnitY();
        using BrushType = Kelvinlet::BrushGrabBiLaplacian;
        _AddDeformer<BrushType>(deformer, f);
    }
    else if (action == 5)
    {
        Vector3 f = Vector3::UnitY();
        using BrushType = Kelvinlet::BrushGrabCusp;
        _AddDeformer<BrushType>(deformer, f);
    }
    else if (action == 6)
    {
        Vector3 f = Vector3::UnitY();
        using BrushType = Kelvinlet::BrushGrabCuspLaplacian;
        _AddDeformer<BrushType>(deformer, f);
    }
    else if (action == 7)
    {
        Vector3 f = Vector3::UnitY();
        using BrushType = Kelvinlet::BrushGrabCuspBiLaplacian;
        _AddDeformer<BrushType>(deformer, f);
    }
    else if (action == 8) // twist
    {
        Vector3 axisAngle = 0.25 * M_PI * Vector3::UnitZ();
        Matrix33 F = Kelvinlet::AssembleSkewSymMatrix(axisAngle);
        using BrushType = Kelvinlet::BrushAffine;
        _AddDeformer<BrushType>(deformer, F);
    }
    else if (action == 9) // scale
    {
        Matrix33 F = Matrix33::Identity();
        using BrushType = Kelvinlet::BrushAffine;
        _AddDeformer<BrushType>(deformer, F);
    }
    else if (action == 10) // pinch
    {
        Matrix33 F = Matrix33::Zero();
        F(0,0) =  0.75;
        F(1,1) = -F(0,0);
        using BrushType = Kelvinlet::BrushAffine;
        _AddDeformer<BrushType>(deformer, F);
    }

    // Load data //
    Vector2 gDx = Vector2::Constant(5.);
    size_t gNum = 100;

    Mesh mesh;
    SampleGrid(mesh, -gDx, gDx, gNum, gNum);

    auto& points = mesh.points;
    int nPts = points.size();

    // Deform //
#ifdef USE_TBB
    tbb::parallel_for(0, nPts, [&](int j)
    {
        points[j] += _Deform(points[j], deformer);
    });
#else
    for(int j = 0; j < nPts; ++j)
    {
        points[j] += _Deform(points[j], deformer);
    }
#endif

    Write(action, mesh, prefix);
    return 1;
}
